# =============================================================================
# Proof-of-concept (B): canonical AnnData approach
#
# An alternative to poc_explore_dataset_2026-04-30.R (the four-SE-subclasses
# version). This POC normalises every input into the canonical AnnData
# layout (per anndata_spec_2026-04-17.qmd) at load time, and reads
# everything through ONE concrete class — `ExploreDatasetH5AD`.
#
# Architecture:
#
#   file.rds / .qs2 / .RData → load_se() → se_to_canonical_anndata(se) ─┐
#   file.h5ad                → anndataR::read_h5ad()                    ├→ ExploreDatasetH5AD
#                                                                       ─┘
#
# SE shape detection + extraction lives in se_to_canonical_anndata(),
# which dispatches to one of four per-shape converters. Each produces an
# in-memory AnnData object conforming to the canonical layout. The reader
# class doesn't know — or care — whether the AnnData came from disk or
# from an SE.
#
# Side benefit: the per-shape converters ARE the SE → .h5ad export
# pipeline that the AnnData migration plan needs anyway. Calling
# `anndataR::write_h5ad(adata, path)` on the converter output produces
# the canonical .h5ad file exactly as Phase 3 of the plan calls for.
#
# Spec: exploreDE/plans/anndata_migration/internal_api_spec_2026-04-30_thin.qmd
#       exploreDE/plans/anndata_migration/anndata_spec_2026-04-17.qmd
# Date: 2026-04-30
# =============================================================================

suppressPackageStartupMessages({
  library(R6)
  library(SummarizedExperiment)
  library(anndataR)
})

# =============================================================================
# 0. Helpers
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Pick the first column whose name matches a regex (case-insensitive).
poc_pick_col <- function(df, pattern) {
  hit <- grep(pattern, colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(hit)) df[[hit[1]]] else NULL
}

# Patterns lifted from R/helperFunctions.R::standardise_columns.
poc_de_patterns <- list(
  effect_log2fc = "^diff$|^log2[._]?(fold)?[._]?(change|ratio|fc)?$|^logfc$|^estimate$",
  score_pvalue  = "^p[._]?value$",
  score_fdr     = "^fdr$|^adj[._]?p[._]?value$|^padj$|^q[._]?value$|^p[._]?adjusted$"
)
poc_rowdata_patterns <- list(
  protein_Id  = "^protein[._]?id$",
  peptide_Id  = "^peptide[._]?id$|^peptide$",
  site        = "^site$",
  gene_name   = "^gene[._]?name$|^gene[._]?symbol$",
  description = "^description$|^desc$",
  nrPeptides  = "^nr[._]?peptides$|^n[._]?peptides$|^num[._]?peptides$",
  modelName   = "^model[._]?name$"
)
poc_standardise_columns <- function(df, patterns) {
  current <- colnames(df)
  for (canon in names(patterns)) {
    idx <- grep(patterns[[canon]], current, ignore.case = TRUE)
    if (length(idx) == 1L) current[idx] <- canon
  }
  colnames(df) <- current
  df
}

# Build a unique, human-readable feature_id for proteomics. Mirrors
# R/helperFunctions.R::build_feature_id.
poc_build_feature_id <- function(df) {
  if (!"protein_Id" %in% colnames(df)) {
    return(gsub("~[^~]+~", "_", rownames(df)))
  }
  pid <- df$protein_Id
  if (!anyDuplicated(pid)) return(pid)
  for (col in intersect(c("site", "peptide_Id"), colnames(df))) {
    combined <- paste(pid, gsub("~", "_", df[[col]]), sep = "_")
    if (!anyDuplicated(combined)) return(combined)
  }
  gsub("~[^~]+~", "_", rownames(df))
}

# Detect SE shape. Returns one of: rna_seq, count_qc, crispr_count_qc,
# proteomics. Identical logic to the subclass POC.
detect_se_shape <- function(se) {
  param  <- S4Vectors::metadata(se)$param %||% list()
  rd_nms <- names(SummarizedExperiment::rowData(se))

  has_real_comparison <- !is.null(param$comparison) &&
    is.character(param$comparison) &&
    length(param$comparison) == 1L &&
    !is.na(param$comparison) &&
    nzchar(param$comparison) &&
    grepl("--over--", param$comparison, fixed = TRUE)
  has_nested_contrasts <- any(grepl("^constrast_", rd_nms))

  if (has_nested_contrasts && !has_real_comparison) return("proteomics")
  if (has_real_comparison) return("rna_seq")

  rn <- rownames(se) %||% character()
  if (any(grepl("_(sgrna)?[0-9]+_[ACGTN]+$", rn, ignore.case = TRUE))) {
    return("crispr_count_qc")
  }
  "count_qc"
}

# Infer organism from refBuild prefix.
poc_infer_organism <- function(param) {
  sp <- param$species
  if (!is.null(sp) && !is.na(sp) && nzchar(sp)) return(sp)
  rb <- param$refBuild
  if (is.null(rb) || is.na(rb) || !nzchar(rb)) return(NA_character_)
  parts <- strsplit(rb, "/", fixed = TRUE)[[1]]
  if (!length(parts)) return(NA_character_)
  gsub("_", " ", parts[1], fixed = TRUE)
}


# =============================================================================
# 1. SE → canonical AnnData converter (the legacy tax payer)
# =============================================================================
# Top-level dispatcher + four per-shape implementations. Each builds an
# in-memory AnnData object that conforms to anndata_spec_2026-04-17.qmd.

se_to_canonical_anndata <- function(se) {
  shape <- detect_se_shape(se)
  switch(shape,
    rna_seq         = .se_to_canon_rnaseq(se),
    count_qc        = .se_to_canon_count_qc(se),
    crispr_count_qc = .se_to_canon_crispr_count_qc(se),
    proteomics      = .se_to_canon_proteomics(se),
    stop("Unknown SE shape: ", shape)
  )
}

# Build canonical obs (sample metadata). [Factor] suffix → factor_ prefix.
.canon_obs <- function(se) {
  cd <- data.frame(SummarizedExperiment::colData(se),
                   check.names = FALSE, stringsAsFactors = FALSE)
  factor_cols  <- grep("\\[Factor]", colnames(cd), value = TRUE)
  factor_names <- gsub(" \\[.*", "", factor_cols)
  factor_names <- gsub(" ",      "_", factor_names)
  obs_factors <- setNames(
    lapply(factor_cols, function(fc) factor(cd[[fc]])),
    paste0("factor_", factor_names)
  )
  obs_rn <- if (!is.null(colnames(se))) colnames(se) else cd$Name
  list(
    obs        = as.data.frame(obs_factors, row.names = obs_rn,
                               check.names = FALSE,
                               stringsAsFactors = FALSE),
    factor_names = factor_names
  )
}

# RNA-seq: top-level rowData has DE inline; counts in X; presentFlag in layers.
.se_to_canon_rnaseq <- function(se) {
  param   <- S4Vectors::metadata(se)$param %||% list()
  obs_blk <- .canon_obs(se)

  # var
  rd <- data.frame(SummarizedExperiment::rowData(se),
                   check.names = FALSE, stringsAsFactors = FALSE)
  rownames(rd) <- rownames(se)
  var <- data.frame(row.names = rownames(rd), check.names = FALSE)
  if ("gene_id"     %in% colnames(rd)) var$label_gene_id     <- rd$gene_id
  if ("gene_name"   %in% colnames(rd)) var$label_gene_symbol <- rd$gene_name
  if ("description" %in% colnames(rd)) var$description       <- gsub(" \\[.*", "", rd$description)
  if ("type"        %in% colnames(rd)) var$label_biotype     <- rd$type
  for (ont in c("GO BP","GO MF","GO CC")) {
    nm <- paste0("label_go_", tolower(sub("GO ", "", ont)))
    if (ont %in% colnames(rd)) var[[nm]] <- rd[[ont]]
  }

  # X = counts; layers = everything else (renamed to canonical vocabulary)
  ass_nms <- SummarizedExperiment::assayNames(se)
  X       <- t(SummarizedExperiment::assay(se, "counts"))
  layer_rename <- c(signal = "normalised", xNorm = "normalised",
                    presentFlag = "present_flag", isPresent = "present_flag")
  layers <- list()
  for (k in setdiff(ass_nms, "counts")) {
    canon_k <- layer_rename[k] %||% k
    layers[[canon_k]] <- t(SummarizedExperiment::assay(se, k))
  }

  # DE → varm + uns$de
  varm <- list(); uns_de <- list()
  comp <- param$comparison
  if (!is.null(comp) && grepl("--over--", comp, fixed = TRUE)) {
    parts    <- strsplit(comp, "--over--", fixed = TRUE)[[1]]
    test_grp <- parts[1]; ref_grp <- parts[2]
    cname    <- paste0(.poc_safe(test_grp), "_vs_", .poc_safe(ref_grp))
    de       <- data.frame(row.names = rownames(rd))
    de$effect_log2fc <- as.numeric(poc_pick_col(rd, poc_de_patterns$effect_log2fc) %||% NA_real_)
    de$score_pvalue  <- as.numeric(poc_pick_col(rd, poc_de_patterns$score_pvalue)  %||% NA_real_)
    de$score_fdr     <- as.numeric(poc_pick_col(rd, poc_de_patterns$score_fdr)     %||% NA_real_)
    varm[[paste0("de_", cname)]] <- as.matrix(de)
    uns_de[[cname]] <- list(
      test_group      = test_grp,
      reference_group = ref_grp,
      factor_name     = param$groupingName %||% NA_character_,
      method          = "DESeq2"
    )
  }

  # Enrichment: lift native clusterProfiler objects into uns$enrichment_native
  md         <- S4Vectors::metadata(se)
  uns_enrich <- list()
  if (length(varm)) {
    cname <- sub("^de_", "", names(varm)[1])
    if (!is.null(md$enrichResult))     uns_enrich[[cname]]$ora_combined <- md$enrichResult
    if (!is.null(md$enrichResultGSEA)) uns_enrich[[cname]]$gsea         <- md$enrichResultGSEA
  }

  uns <- list(
    schema_version    = "exploreDE/1.0",
    omics_type        = "rnaseq_bulk",
    hierarchy         = list(levels = param$featureLevel %||% "gene",
                             unit   = param$featureLevel %||% "gene"),
    design            = list(factors  = paste0("factor_", obs_blk$factor_names),
                             grouping = paste0("factor_", param$groupingName %||% obs_blk$factor_names[1])),
    de                = uns_de,
    enrichment_native = uns_enrich,
    reference_genome  = param$refBuild %||% NA_character_,
    organism          = poc_infer_organism(param),
    app_hints         = list(
      default_contrast = if (length(uns_de)) names(uns_de)[1] else NA_character_,
      default_assay    = if ("log2" %in% names(layers)) "log2" else NA_character_,
      log2_pseudocount = param$backgroundExpression %||% 10
    )
  )

  anndataR::AnnData(X = X, obs = obs_blk$obs, var = var,
                    layers = layers, varm = varm, uns = uns)
}

# CountQC: top-level rowData, no DE.
.se_to_canon_count_qc <- function(se) {
  param   <- S4Vectors::metadata(se)$param %||% list()
  obs_blk <- .canon_obs(se)

  rd <- data.frame(SummarizedExperiment::rowData(se),
                   check.names = FALSE, stringsAsFactors = FALSE)
  rownames(rd) <- rownames(se)
  var <- data.frame(row.names = rownames(rd), check.names = FALSE)
  if ("gene_id"     %in% colnames(rd)) var$label_gene_id     <- rd$gene_id
  if ("gene_name"   %in% colnames(rd)) var$label_gene_symbol <- rd$gene_name
  if ("description" %in% colnames(rd)) var$description       <- gsub(" \\[.*", "", rd$description)

  ass_nms <- SummarizedExperiment::assayNames(se)
  X       <- t(SummarizedExperiment::assay(se, "counts"))
  layer_rename <- c(signal = "normalised", xNorm = "normalised",
                    presentFlag = "present_flag", isPresent = "present_flag")
  layers <- list()
  for (k in setdiff(ass_nms, "counts")) {
    canon_k <- layer_rename[k] %||% k
    layers[[canon_k]] <- t(SummarizedExperiment::assay(se, k))
  }

  uns <- list(
    schema_version    = "exploreDE/1.0",
    omics_type        = "countqc",
    hierarchy         = list(levels = "gene", unit = "gene"),
    design            = list(factors  = paste0("factor_", obs_blk$factor_names),
                             grouping = paste0("factor_", param$groupingName %||% obs_blk$factor_names[1])),
    de                = list(),
    enrichment_native = list(),
    reference_genome  = param$refBuild %||% NA_character_,
    organism          = poc_infer_organism(param)
  )

  anndataR::AnnData(X = X, obs = obs_blk$obs, var = var,
                    layers = layers, uns = uns)
}

# CrisprCountQC: CountQC + sgRNA → gene mapping in label_target_gene.
.se_to_canon_crispr_count_qc <- function(se) {
  adata <- .se_to_canon_count_qc(se)
  rn    <- rownames(adata$var)
  adata$var$label_target_gene <- gsub("_(sgrna)?[0-9]+_[ACGTN]+$", "", rn,
                                      ignore.case = TRUE)
  adata$uns$omics_type <- "crispr_countqc"
  adata$uns$hierarchy  <- list(levels = c("gene", "guide"), unit = "guide")
  adata
}

# Proteomics: nested rowData$constrast_*; build feature_id; remap assays.
.se_to_canon_proteomics <- function(se) {
  param   <- S4Vectors::metadata(se)$param %||% list()

  # obs: factor detection by exclusion regex (proteomics convention)
  cd <- data.frame(SummarizedExperiment::colData(se),
                   check.names = FALSE, stringsAsFactors = FALSE)
  excluded     <- grep("^Sample$|^Name$|^names$|^channel$|^raw\\.file$",
                       colnames(cd), value = TRUE, ignore.case = TRUE)
  factor_cols  <- setdiff(colnames(cd), excluded)
  factor_names <- sub("^G_", "", factor_cols)
  factor_names[!nzchar(factor_names)] <- "Group"
  obs_factors <- setNames(
    lapply(factor_cols, function(fc) factor(cd[[fc]])),
    paste0("factor_", factor_names)
  )
  obs_rn <- if (!is.null(colnames(se))) colnames(se) else
            (cd$Name %||% cd$name %||% cd$Sample)
  obs <- as.data.frame(obs_factors, row.names = obs_rn,
                       check.names = FALSE, stringsAsFactors = FALSE)

  # Standardise the first contrast → reference for var + feature_id
  rd_nms    <- names(SummarizedExperiment::rowData(se))
  contrasts <- gsub("^constrast_", "",
                    grep("^constrast_", rd_nms, value = TRUE))
  if (!length(contrasts)) {
    stop("Proteomics SE has no constrast_* rowData columns")
  }

  raw_dfs <- setNames(
    lapply(contrasts, function(cn) {
      data.frame(SummarizedExperiment::rowData(se)[[paste0("constrast_", cn)]],
                 check.names = FALSE, stringsAsFactors = FALSE)
    }),
    contrasts
  )
  std_dfs <- lapply(raw_dfs, poc_standardise_columns,
                    patterns = poc_rowdata_patterns)

  # Reference DataFrame (UNSORTED — alignment with rownames(se) preserved).
  ref_df <- std_dfs[[1]]
  rownames(ref_df) <- rownames(se)[seq_len(nrow(ref_df))]
  ref_df$feature_id <- poc_build_feature_id(ref_df)

  # SE rowname → built feature_id (for assay remap)
  fid_lookup <- setNames(ref_df$feature_id, rownames(ref_df))

  # Hierarchy / omics-subtype detection
  hierarchy_levels <-
    if      ("site"       %in% colnames(ref_df)) c("protein", "site")
    else if ("peptide_Id" %in% colnames(ref_df)) c("protein", "peptide")
    else                                          "protein"

  # var (canonical labels)
  var <- data.frame(row.names = ref_df$feature_id, check.names = FALSE)
  if ("protein_Id"  %in% colnames(ref_df)) var$label_uniprot_id  <- ref_df$protein_Id
  if ("gene_name"   %in% colnames(ref_df)) var$label_gene_symbol <- ref_df$gene_name
  if ("description" %in% colnames(ref_df)) {
    var$description <- gsub(" \\[.*", "", as.character(ref_df$description))
  }
  if ("nrPeptides"  %in% colnames(ref_df)) var$label_nr_peptides <- ref_df$nrPeptides
  if ("site"        %in% colnames(ref_df)) var$hierarchy_site    <- ref_df$site
  if ("peptide_Id"  %in% colnames(ref_df)) var$hierarchy_peptide <- ref_df$peptide_Id
  if ("modelName"   %in% colnames(ref_df)) var$label_model_name  <- ref_df$modelName

  # X = primary assay (transformed/intensities); rownames remapped via fid_lookup
  ass_nms <- SummarizedExperiment::assayNames(se)
  primary_key <- intersect(c("transformedData", "transformed",
                             "intensities", "log2"), ass_nms)
  primary_key <- if (length(primary_key)) primary_key[1] else ass_nms[1]
  remap_rownames <- function(m) {
    new_rn <- unname(fid_lookup[rownames(m)])
    keep   <- !is.na(new_rn)
    m      <- m[keep, , drop = FALSE]
    rownames(m) <- new_rn[keep]
    m
  }
  X      <- t(remap_rownames(SummarizedExperiment::assay(se, primary_key)))
  layers <- list()
  for (k in setdiff(ass_nms, primary_key)) {
    canon_k <- if (k %in% c("rawData", "intensities_raw")) "intensities_raw" else k
    layers[[canon_k]] <- t(remap_rownames(SummarizedExperiment::assay(se, k)))
  }

  # varm + uns$de per contrast
  varm <- list(); uns_de <- list()
  for (i in seq_along(contrasts)) {
    cn  <- contrasts[i]
    df  <- std_dfs[[i]]
    df$feature_id <- ref_df$feature_id[seq_len(nrow(df))]

    de <- data.frame(row.names = df$feature_id)
    de$effect_log2fc <- as.numeric(poc_pick_col(df, poc_de_patterns$effect_log2fc) %||% NA_real_)
    de$score_pvalue  <- as.numeric(poc_pick_col(df, poc_de_patterns$score_pvalue)  %||% NA_real_)
    de$score_fdr     <- as.numeric(poc_pick_col(df, poc_de_patterns$score_fdr)     %||% NA_real_)
    if ("nrPeptides" %in% colnames(df)) de$label_n_peptides <- as.integer(df$nrPeptides)
    varm[[paste0("de_", cn)]] <- as.matrix(de)

    grp <- strsplit(cn, "_vs_|--over--")[[1]]
    uns_de[[cn]] <- list(
      test_group      = grp[1] %||% NA_character_,
      reference_group = grp[2] %||% NA_character_,
      factor_name     = factor_names[1] %||% NA_character_,
      method          = "prolfqua"
    )
  }

  uns <- list(
    schema_version    = "exploreDE/1.0",
    omics_type        = "proteomics",
    hierarchy         = list(levels = hierarchy_levels,
                             unit   = utils::tail(hierarchy_levels, 1)),
    design            = list(factors  = paste0("factor_", factor_names),
                             grouping = paste0("factor_", factor_names[1])),
    de                = uns_de,
    enrichment_native = list(),    # spec mandates exporter run clusterProfiler;
                                   # left empty in this POC since SUSHI doesn't
                                   # produce proteomics enrichment yet.
    reference_genome  = NA_character_,
    organism          = NA_character_,
    app_hints         = list(default_contrast = if (length(uns_de)) names(uns_de)[1] else NA_character_)
  )

  anndataR::AnnData(X = X, obs = obs, var = var,
                    layers = layers, varm = varm, uns = uns)
}

.poc_safe <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)


# =============================================================================
# 2. Abstract base class: ExploreDataset
# =============================================================================
# Identical to subclass POC's abstract base. Config-only methods implemented;
# backend-touching methods stub to stop("not implemented").

ExploreDataset <- R6::R6Class(
  classname = "ExploreDataset",

  public = list(

    initialize = function(...) {
      stop("ExploreDataset is abstract. ",
           "Use new_explore_dataset().")
    },

    # identity (config-only, derived from uns inside concrete subclasses)
    data_type        = function() stop("not implemented"),
    has_de           = function() length(self$contrasts()) > 0L,
    has_enrichment   = function() length(self$enrichment_kinds()) > 0L,
    organism         = function() stop("not implemented"),
    ref_build        = function() stop("not implemented"),
    feature_id_kind  = function() stop("not implemented"),
    backend_kind     = function() stop("not implemented"),
    backend          = function() stop("not implemented"),

    # samples
    samples          = function() stop("not implemented"),
    sample_ids       = function() rownames(self$samples()),
    n_samples        = function() nrow(self$samples()),
    factor_names     = function() stop("not implemented"),
    factor_levels    = function(factor) stop("not implemented"),
    primary_grouping = function() stop("not implemented"),
    sample_factor    = function(factor) self$samples()[[factor]],

    # features
    features         = function() stop("not implemented"),
    feature_ids      = function() rownames(self$features()),
    n_features       = function() nrow(self$features()),
    feature_display  = function(ids = NULL) {
      f <- self$features()
      if (is.null(ids)) ids <- rownames(f)
      setNames(f[ids, "display_name"], ids)
    },

    # assays
    assay_names        = function() stop("not implemented"),
    has_assay          = function(name) name %in% self$assay_names(),
    assay              = function(name, features = NULL, samples = NULL) {
      stop("not implemented")
    },
    raw_assay_name     = function() stop("not implemented"),
    primary_assay_name = function() stop("not implemented"),
    raw_assay          = function(...) self$assay(self$raw_assay_name(), ...),
    primary_assay      = function(...) self$assay(self$primary_assay_name(), ...),

    # DE
    contrasts        = function() stop("not implemented"),
    primary_contrast = function() stop("not implemented"),
    de               = function(contrast) stop("not implemented"),
    de_all           = function() {
      setNames(lapply(self$contrasts(), self$de), self$contrasts())
    },
    contrast_meta    = function(contrast) stop("not implemented"),

    # enrichment
    enrichment_kinds = function() stop("not implemented"),
    enrichment       = function(contrast, kind) stop("not implemented"),

    # source
    source           = function() stop("not implemented"),

    # print
    print = function(...) {
      cat("<", class(self)[1], "  (", self$backend_kind(), ")>\n", sep = "")
      cat("  data_type:        ", self$data_type(), "\n", sep = "")
      cat("  feature_id_kind:  ", self$feature_id_kind(), "\n", sep = "")
      cat("  organism:         ", self$organism(), "\n", sep = "")
      cat("  ref_build:        ", self$ref_build(), "\n", sep = "")
      cat("  samples × features: ",
          self$n_samples(), " × ", self$n_features(), "\n", sep = "")
      cat("  factors (", length(self$factor_names()), "):\n", sep = "")
      for (f in self$factor_names()) {
        lv <- self$factor_levels(f)
        cat("    ", f,
            if (identical(f, self$primary_grouping())) "  (primary)" else "",
            ": ",
            paste(utils::head(lv, 6), collapse = ", "),
            if (length(lv) > 6) paste0(", ... (+", length(lv) - 6, ")") else "",
            "\n", sep = "")
      }
      cat("  contrasts (",
          length(self$contrasts()), "): ",
          paste(head(self$contrasts(), 3), collapse = ", "),
          if (length(self$contrasts()) > 3) ", ..." else "",
          "\n", sep = "")
      cat("  primary_contrast: ",
          self$primary_contrast() %||% "(none)", "\n", sep = "")
      cat("  assays:           ",
          paste(self$assay_names(), collapse = ", "), "\n", sep = "")
      cat("  enrichment:       ",
          if (self$has_enrichment())
            paste(self$enrichment_kinds(), collapse = ", ")
          else "(none)",
          "\n", sep = "")
      invisible(self)
    }
  )
)


# =============================================================================
# 3. The single concrete reader: ExploreDatasetH5AD
# =============================================================================
# Reads an in-memory `anndataR::AnnData` object that conforms to the
# canonical layout (anndata_spec_2026-04-17.qmd). Doesn't know or care
# whether the AnnData was loaded from .h5ad or built by se_to_canonical_anndata.

ExploreDatasetH5AD <- R6::R6Class(
  classname = "ExploreDatasetH5AD",
  inherit   = ExploreDataset,

  public = list(

    initialize = function(adata) {
      stopifnot(inherits(adata, "AbstractAnnData"))
      private$.adata <- adata
      invisible(self)
    },

    backend      = function() private$.adata,
    backend_kind = function() "h5ad",

    data_type        = function() private$.adata$uns$omics_type %||% NA_character_,
    organism         = function() private$.adata$uns$organism %||% NA_character_,
    ref_build        = function() private$.adata$uns$reference_genome %||% NA_character_,
    feature_id_kind  = function() private$.adata$uns$hierarchy$unit %||% NA_character_,

    # samples: strip factor_ prefix from column names for display
    samples = function() {
      obs <- private$.adata$obs
      cn  <- colnames(obs)
      colnames(obs) <- sub("^factor_", "", cn)
      obs
    },

    factor_names = function() {
      sub("^factor_", "",
          grep("^factor_", colnames(private$.adata$obs), value = TRUE))
    },

    factor_levels = function(factor) {
      col <- private$.adata$obs[[paste0("factor_", factor)]]
      if (is.factor(col)) levels(col) else sort(unique(as.character(col)))
    },

    primary_grouping = function() {
      g <- private$.adata$uns$design$grouping %||% NA_character_
      sub("^factor_", "", g)
    },

    # features: var with feature_id (rownames) + display_name (label_gene_symbol →
    # label_uniprot_id → feature_id fallback chain)
    features = function() {
      var <- private$.adata$var
      var$feature_id <- rownames(var)

      display <- rep(NA_character_, nrow(var))
      for (col in c("label_gene_symbol", "label_uniprot_id")) {
        if (!col %in% colnames(var)) next
        blank <- is.na(display) | !nzchar(display)
        if (!any(blank)) break
        vals <- as.character(var[[col]])
        display[blank] <- vals[blank]
      }
      blank <- is.na(display) | !nzchar(display)
      display[blank] <- rownames(var)[blank]
      var$display_name <- display
      var
    },

    # assay_names: X-name-for-omics + layer keys
    assay_names = function() {
      x_name <- private$.x_name()
      c(x_name, names(private$.adata$layers %||% list()))
    },

    raw_assay_name = function() {
      switch(self$data_type(),
        rnaseq_bulk     = "counts",
        countqc         = "counts",
        crispr_countqc  = "counts",
        proteomics      = "intensities",
        private$.x_name()
      )
    },

    primary_assay_name = function() {
      hint <- private$.adata$uns$app_hints$default_assay
      if (!is.null(hint) && !is.na(hint) && nzchar(hint)) return(hint)
      self$raw_assay_name()
    },

    assay = function(name, features = NULL, samples = NULL) {
      x_name <- private$.x_name()
      m <- if (identical(name, x_name)) {
        t(private$.adata$X)
      } else if (name == "normalised_log2") {
        # Computed: log2(normalised + pseudocount)
        norm <- self$assay("normalised")
        log2(norm + (private$.adata$uns$app_hints$log2_pseudocount %||% 10))
      } else {
        layer <- private$.adata$layers[[name]]
        if (is.null(layer)) stop("Unknown assay: ", name)
        t(layer)
      }
      if (!is.null(features)) m <- m[features, , drop = FALSE]
      if (!is.null(samples))  m <- m[, samples,  drop = FALSE]
      m
    },

    # Contrasts: keys of varm matching ^de_
    contrasts = function() {
      sub("^de_", "",
          grep("^de_", names(private$.adata$varm %||% list()), value = TRUE))
    },

    primary_contrast = function() {
      hint <- private$.adata$uns$app_hints$default_contrast
      if (!is.null(hint) && !is.na(hint) && nzchar(hint)) return(hint)
      cs <- self$contrasts()
      if (length(cs)) cs[1] else NA_character_
    },

    de = function(contrast) {
      stopifnot(contrast %in% self$contrasts())
      m  <- private$.adata$varm[[paste0("de_", contrast)]]
      df <- as.data.frame(m, check.names = FALSE)
      df$feature_id <- rownames(df) %||% rownames(private$.adata$var)
      df
    },

    contrast_meta = function(contrast) {
      stopifnot(contrast %in% self$contrasts())
      private$.adata$uns$de[[contrast]]
    },

    # Enrichment: read uns$enrichment_native[<contrast>]; group by kind
    # ("ora" matches ^ora; "gsea" matches ^gsea).
    enrichment_kinds = function() {
      en <- private$.adata$uns$enrichment_native %||% list()
      if (!length(en)) return(character())
      methods <- unique(unlist(lapply(en, names)))
      kinds   <- character()
      if (any(grepl("^ora",  methods))) kinds <- c(kinds, "ora")
      if (any(grepl("^gsea", methods))) kinds <- c(kinds, "gsea")
      kinds
    },

    enrichment = function(contrast, kind) {
      stopifnot(kind %in% c("ora", "gsea"))
      en <- private$.adata$uns$enrichment_native[[contrast]]
      if (is.null(en)) return(list())
      pat <- if (kind == "ora") "^ora" else "^gsea"
      en[grepl(pat, names(en))]
    },

    source = function() {
      list(
        backend_kind   = "h5ad",
        schema_version = private$.adata$uns$schema_version %||% NA_character_,
        omics_type     = private$.adata$uns$omics_type     %||% NA_character_
      )
    }
  ),

  private = list(
    .adata = NULL,

    # Name of the canonical "X" assay for this omics_type. Determines what
    # layers are exposed under and whether eds$assay("counts") aliases to X.
    .x_name = function() {
      switch(private$.adata$uns$omics_type %||% "",
        rnaseq_bulk     = "counts",
        countqc         = "counts",
        crispr_countqc  = "counts",
        proteomics      = "intensities",
        metabolomics    = "intensities",
        lipidomics      = "intensities",
        "X"
      )
    }
  )
)


# =============================================================================
# 4. Factory + loader
# =============================================================================

new_explore_dataset <- function(input) {
  adata <- if (inherits(input, "AbstractAnnData")) {
    input
  } else if (inherits(input, "SummarizedExperiment")) {
    se_to_canonical_anndata(input)
  } else if (is.character(input) && length(input) == 1L) {
    if (grepl("\\.h5ad$", input, ignore.case = TRUE) && file.exists(input)) {
      anndataR::read_h5ad(input)
    } else {
      se <- poc_load(input)
      se_to_canonical_anndata(se)
    }
  } else {
    stop("Unrecognised input for new_explore_dataset()")
  }
  ExploreDatasetH5AD$new(adata)
}

poc_load <- function(input) {
  if (grepl("^https?://", input)) {
    message("[load] fetching URL: ", input)
    return(readRDS(url(input)))
  }
  if (grepl("\\.RData$", input, ignore.case = TRUE)) {
    dir_input <- gsub("/(result|counts)-.*-EzResult\\.RData$", "", input)
    return(poc_load(dir_input))
  }
  if (grepl("\\.rds$", input, ignore.case = TRUE) && file.exists(input)) {
    return(readRDS(input))
  }
  if (grepl("\\.qs2$", input, ignore.case = TRUE) && file.exists(input)) {
    return(qs2::qs_read(input, nthreads = 4))
  }
  if (dir.exists(input)) {
    for (fname in c("deResult.rds",        "deResult.qs2",
                    "rawData_final.rds",   "rawData_final.qs2",
                    "rawData.rds",         "rawData.qs2")) {
      fp <- file.path(input, fname)
      if (file.exists(fp)) {
        message("[load] reading ", fp)
        if (grepl("\\.qs2$", fname)) return(qs2::qs_read(fp, nthreads = 4))
        return(readRDS(fp))
      }
    }
    stop("No known SE artefact in directory: ", input)
  }
  stop("Cannot interpret input: ", input)
}


# =============================================================================
# 5. Demo runner — same shape as the subclass POC for head-to-head comparison
# =============================================================================

demo_paths <- list(
  rnaseq_teaching =
    "/srv/gstore/projects/p3009/o5638_DESeq2_diff--over--undiff_2024-11-20--12-29-40/diff--over--undiff/",
  mageck_count_qc =
    "/srv/gstore/projects/p36911/o39664_o37783_o39124_o38710_o41353_ExploreMageckCounts_2026-04-09--15-07-24/ExploreMageckCounts/counts-yuunbieawvpw-EzResult.RData",
  rnaseq_two_factors =
    "/srv/gstore/projects/p29779/o40890_DESeq2_shp65_w_HP--over--shp65_wo_HP_2026-04-16--13-38-37/shp65_w_HP--over--shp65_wo_HP",
  rnaseq_isoform =
    "/srv/gstore/projects/p34429/o35733_DESeq2_Palbo_day6--over--Palbo_day2_2024-08-22--11-21-44/Palbo_day6--over--Palbo_day2/",
  prot_protein =
    "https://fgcz-ms.uzh.ch/public/pStore/p23078/bfabric/Proteomics/exploreDE/2026/2026-03/2026-03-10/workunit_342063/3106962.rds",
  prot_phospho =
    "https://fgcz-ms.uzh.ch/public/pStore/p39819/bfabric/Proteomics/exploreDE/2025/2025-10/2025-10-22/workunit_335485/2997996.rds"
)

demo_one <- function(label, path) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("[", label, "]\n", sep = "")
  cat("path: ", path, "\n", sep = "")
  cat(strrep("-", 78), "\n", sep = "")

  eds <- tryCatch(
    new_explore_dataset(path),
    error = function(e) {
      cat("LOAD FAILED: ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  if (is.null(eds)) return(invisible(NULL))

  print(eds)

  cat("\nsamples() (head 3):\n")
  print(utils::head(eds$samples(), 3))

  cat("\nfeatures() (head 3, key cols):\n")
  f <- eds$features()
  show_cols <- intersect(c("feature_id", "display_name", "label_gene_id",
                           "label_gene_symbol", "label_uniprot_id",
                           "description"),
                         colnames(f))
  print(utils::head(f[, show_cols, drop = FALSE], 3))

  if (length(eds$assay_names()) >= 1) {
    a_name <- eds$raw_assay_name() %||% eds$assay_names()[1]
    cat("\nraw_assay '", a_name, "' dim: ",
        paste(dim(eds$assay(a_name)), collapse = " × "), "\n", sep = "")
  }

  if (eds$has_de()) {
    c1 <- eds$primary_contrast() %||% eds$contrasts()[1]
    cat("\nde('", c1, "') (head 3):\n", sep = "")
    de_df <- eds$de(c1)
    show <- intersect(c("feature_id", "effect_log2fc", "score_pvalue",
                        "score_fdr", "label_n_peptides"),
                      c(colnames(de_df), "feature_id"))
    print(utils::head(de_df[, show, drop = FALSE], 3))

    cm <- eds$contrast_meta(c1)
    cat("\ncontrast_meta('", c1, "'):\n", sep = "")
    cat("  test_group:      ", cm$test_group, "\n", sep = "")
    cat("  reference_group: ", cm$reference_group, "\n", sep = "")
    cat("  method:          ", cm$method, "\n", sep = "")
  }

  if (eds$has_enrichment()) {
    for (kind in eds$enrichment_kinds()) {
      e <- eds$enrichment(eds$primary_contrast() %||% "default", kind)
      cat("\nenrichment(", kind, "): list of ", length(e),
          " — names: ", paste(head(names(e), 5), collapse = ", "),
          "\n", sep = "")
      if (length(e)) {
        first <- e[[1]]
        # If the entry is itself a list (e.g. gsea = list(BP, MF, CC)), peek deeper.
        if (is.list(first) && !inherits(first, c("enrichResult","gseaResult"))
            && length(first)) {
          first <- first[[1]]
        }
        cat("  first leaf class: ",
            paste(class(first), collapse = "/"), "\n", sep = "")
      }
    }
  }
  invisible(eds)
}

run_all_demos <- function(which = names(demo_paths)) {
  results <- list()
  for (label in which) {
    results[[label]] <- demo_one(label, demo_paths[[label]])
  }
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("Summary: ", sum(!vapply(results, is.null, logical(1))),
      " / ", length(results), " demos loaded\n", sep = "")
  invisible(results)
}

# Uncomment to run:
results <- run_all_demos()

# Or run a single demo:
eds <- demo_one("rnaseq_teaching", demo_paths$rnaseq_teaching)

# Side benefit — the converter IS the .h5ad export path:
se    <- readRDS(demo_paths$rnaseq_teaching |> file.path("deResult.rds"))
adata <- se_to_canonical_anndata(se)
anndataR::write_h5ad(adata, "/tmp/teaching.h5ad")
