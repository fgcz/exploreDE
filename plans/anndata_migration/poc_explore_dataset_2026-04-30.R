# =============================================================================
# Proof-of-concept: ExploreDataset R6 API (Draft 2 / thin variant)
#
# Demonstrates the proposed storage-agnostic API end-to-end against 6 distinct
# real-world data shapes. Class hierarchy:
#
#   ExploreDataset (abstract)                 — config-only methods
#   └── ExploreDatasetSE (abstract)           — backend ref, shared SE plumbing
#       ├── ExploreDatasetSE_RNASeq           — top-level rowData + DE inline
#       ├── ExploreDatasetSE_CountQC          — top-level rowData, no DE
#       │   └── ExploreDatasetSE_CrisprCountQC — + sgrna_to_gene
#       └── ExploreDatasetSE_Proteomics       — nested rowData$constrast_*
#
# The factory `new_explore_dataset()` detects the SE shape and instantiates
# the appropriate subclass. Read methods (samples, features, de, assay,
# enrichment) are uniform across subclasses — the only branching is the
# detection switch in the factory.
#
# Each concrete subclass implements two private methods:
#   .build_config(se)  → returns the canonical config list
#   .prep()            → populates private$.features_df, .de_dfs, .fid_lookup
#                        from the SE
#
# Spec: exploreDE/plans/anndata_migration/internal_api_spec_2026-04-30_thin.qmd
# Date: 2026-04-30
# =============================================================================

suppressPackageStartupMessages({
  library(R6)
  library(SummarizedExperiment)
})

# =============================================================================
# 0. Helpers
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Recursive flatten of nested clusterProfiler results into a single named list.
poc_flatten_cp_results <- function(x, prefix = "") {
  if (inherits(x, c("enrichResult", "gseaResult"))) {
    nm <- if (nzchar(prefix)) prefix else "default"
    return(setNames(list(x), nm))
  }
  if (is.list(x) && length(x)) {
    out <- list()
    nms <- names(x) %||% as.character(seq_along(x))
    for (i in seq_along(x)) {
      sub_prefix <- if (nzchar(prefix)) paste0(prefix, "/", nms[i]) else nms[i]
      out <- c(out, poc_flatten_cp_results(x[[i]], sub_prefix))
    }
    return(out)
  }
  list()
}

# SUSHI param$species is rarely set; fall back to first slash-segment of
# refBuild ("Homo_sapiens/GENCODE/..." → "Homo sapiens").
poc_infer_organism <- function(param) {
  sp <- param$species
  if (!is.null(sp) && !is.na(sp) && nzchar(sp)) return(sp)
  rb <- param$refBuild
  if (is.null(rb) || is.na(rb) || !nzchar(rb)) return(NA_character_)
  parts <- strsplit(rb, "/", fixed = TRUE)[[1]]
  if (!length(parts)) return(NA_character_)
  gsub("_", " ", parts[1], fixed = TRUE)
}

# Standardise DE column names. Patterns lifted verbatim from
# R/helperFunctions.R::standardise_columns; canonical outputs in snake_case.
poc_de_column_patterns <- list(
  log2_fc = "^diff$|^log2[._]?(fold)?[._]?(change|ratio|fc)?$|^logfc$|^estimate$",
  p_value = "^p[._]?value$",
  p_adj   = "^fdr$|^adj[._]?p[._]?value$|^padj$|^q[._]?value$|^p[._]?adjusted$"
)

poc_standardise_de_columns <- function(df) {
  current <- colnames(df)
  for (canon in names(poc_de_column_patterns)) {
    pattern <- poc_de_column_patterns[[canon]]
    idx <- grep(pattern, current, ignore.case = TRUE)
    if (length(idx) >= 1L && !canon %in% colnames(df)) {
      df[[canon]] <- df[[idx[1]]]
    }
  }
  df
}

# Standardise proteomics rowData/contrast-DataFrame column names. Mirrors
# R/helperFunctions.R::standardise_columns but kept fresh for the POC.
poc_rowdata_column_patterns <- list(
  protein_Id  = "^protein[._]?id$",
  peptide_Id  = "^peptide[._]?id$|^peptide$",
  site        = "^site$",
  gene_name   = "^gene[._]?name$|^gene[._]?symbol$",
  description = "^description$|^desc$",
  fasta.id    = "^fasta[._]?id$",
  modelName   = "^model[._]?name$",
  nrPeptides  = "^nr[._]?peptides$|^n[._]?peptides$|^num[._]?peptides$"
)

poc_standardise_rowdata_columns <- function(df) {
  current <- colnames(df)
  for (canon in names(poc_rowdata_column_patterns)) {
    pattern <- poc_rowdata_column_patterns[[canon]]
    idx <- grep(pattern, current, ignore.case = TRUE)
    if (length(idx) == 1L) current[idx] <- canon
  }
  colnames(df) <- current
  df
}

# Build a unique, human-readable feature_id for proteomics features. Tries
# protein_Id alone, else protein_Id + site, else protein_Id + peptide_Id,
# else cleaned rownames. Mirrors R/helperFunctions.R::build_feature_id.
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

# Detect SE shape from its content. Returns one of:
#   "rna_seq", "count_qc", "crispr_count_qc", "proteomics".
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


# =============================================================================
# 1. Abstract base class: ExploreDataset
# =============================================================================
# Defines the contract. Methods that touch the storage backend are virtual
# (`stop("not implemented")`); methods that read only from `private$.config`
# are implemented here and inherited by every subclass.

ExploreDataset <- R6::R6Class(
  classname = "ExploreDataset",

  public = list(

    initialize = function(...) {
      stop("ExploreDataset is abstract. ",
           "Use new_explore_dataset() or one of the concrete subclasses.")
    },

    # ---- identity (config-only) -------------------------------------------
    data_type        = function() private$.config$data_type,
    has_de           = function() length(self$contrasts()) > 0L,
    has_enrichment   = function() {
      length(private$.config$enrichment_kinds %||% character()) > 0L
    },
    organism         = function() private$.config$organism %||% NA_character_,
    ref_build        = function() private$.config$ref_build %||% NA_character_,
    feature_id_kind  = function() private$.config$feature_id_kind,
    backend_kind     = function() stop("not implemented"),
    backend          = function() stop("not implemented"),

    # ---- samples (samples() virtual; rest derived) ------------------------
    samples          = function() stop("not implemented"),
    sample_ids       = function() rownames(self$samples()),
    n_samples        = function() nrow(self$samples()),
    factor_names     = function() private$.config$factors %||% character(),
    factor_levels    = function(factor) {
      private$.config$factor_levels[[factor]]
    },
    primary_grouping = function() private$.config$primary_grouping,
    sample_factor    = function(factor) self$samples()[[factor]],

    # ---- features (features() virtual; rest derived) ----------------------
    features         = function() stop("not implemented"),
    feature_ids      = function() rownames(self$features()),
    n_features       = function() nrow(self$features()),
    feature_display  = function(ids = NULL) {
      f <- self$features()
      if (is.null(ids)) ids <- rownames(f)
      setNames(f[ids, "display_name"], ids)
    },

    # ---- assays (assay_names + assay virtual; rest derived) ---------------
    assay_names        = function() stop("not implemented"),
    has_assay          = function(name) name %in% self$assay_names(),
    assay              = function(name, features = NULL, samples = NULL) {
      stop("not implemented")
    },
    raw_assay_name     = function() private$.config$raw_assay_name,
    primary_assay_name = function() private$.config$primary_assay_name,
    raw_assay          = function(...) self$assay(self$raw_assay_name(), ...),
    primary_assay      = function(...) self$assay(self$primary_assay_name(), ...),

    # ---- DE (de() virtual; rest config-only) ------------------------------
    contrasts        = function() {
      names(private$.config$contrast_meta) %||% character()
    },
    primary_contrast = function() private$.config$primary_contrast,
    de               = function(contrast) stop("not implemented"),
    de_all           = function() {
      setNames(lapply(self$contrasts(), self$de), self$contrasts())
    },
    contrast_meta    = function(contrast) {
      stopifnot(contrast %in% self$contrasts())
      private$.config$contrast_meta[[contrast]]
    },

    # ---- enrichment (all virtual except enrichment_kinds) -----------------
    enrichment_kinds = function() {
      private$.config$enrichment_kinds %||% character()
    },
    enrichment       = function(contrast, kind) stop("not implemented"),
    enrichment_input = function(contrast)        stop("not implemented"),

    # ---- domain extras (config-only) --------------------------------------
    sgrna_to_gene    = function() private$.config$sgrna_to_gene,
    proteomics_choices = function(kind = c("gene","protein")) {
      kind <- match.arg(kind)
      private$.config$proteomics[[kind]]
    },

    # ---- provenance --------------------------------------------------------
    source           = function() private$.config$source,

    # ---- lifecycle ---------------------------------------------------------
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
  ),

  private = list(
    .config = NULL
  )
)


# =============================================================================
# 2. Abstract intermediate: ExploreDatasetSE
# =============================================================================
# Holds the SE backend reference and the shared SE plumbing. Concrete
# subclasses override `.build_config(se)` and `.prep()` only — everything
# else is implemented uniformly.

ExploreDatasetSE <- R6::R6Class(
  classname = "ExploreDatasetSE",
  inherit   = ExploreDataset,

  public = list(

    initialize = function(se, config = NULL) {
      stopifnot(inherits(se, "SummarizedExperiment"))
      private$.se     <- se
      private$.config <- config %||% private$.build_config(se)
      private$.prep()
      invisible(self)
    },

    backend      = function() private$.se,
    backend_kind = function() "se",

    # samples(): uniform colData read + name cleanup. Both SUSHI [Factor]
    # suffix and proteomics G_ prefix are stripped here regardless of
    # data_type; whichever is present (if any) gets cleaned.
    samples = function() {
      cd <- data.frame(SummarizedExperiment::colData(private$.se), check.names = FALSE, stringsAsFactors = FALSE)
      colnames(cd) <- gsub(" \\[.*", "", colnames(cd))
      colnames(cd) <- gsub(" ", "_", colnames(cd))
      cn   <- colnames(cd)
      is_g <- grepl("^G_", cn)
      if (any(is_g)) {
        new <- sub("^G_", "", cn[is_g])
        new[!nzchar(new)] <- "Group"
        colnames(cd)[is_g] <- new
      }
      se_cn <- colnames(private$.se)
      if (!is.null(se_cn) && length(se_cn) == nrow(cd)) {
        rownames(cd) <- se_cn
      } else if ("Name" %in% colnames(cd)) {
        rownames(cd) <- as.character(cd$Name)
      } else if ("name" %in% colnames(cd)) {
        rownames(cd) <- as.character(cd$name)
      }
      cd
    },

    # features(), de(), assay_names(): uniform — read from .prep()'s output.
    features    = function() private$.features_df,

    de = function(contrast) {
      stopifnot(contrast %in% self$contrasts())
      private$.de_dfs[[contrast]] %||% data.frame()
    },

    assay_names = function() {
      backend_keys <- SummarizedExperiment::assayNames(private$.se)
      keys         <- names(private$.config$assay_map %||% list())
      keys[vapply(keys, function(k) {
        v <- private$.config$assay_map[[k]]
        is.null(v) || v %in% backend_keys
      }, logical(1))]
    },

    assay = function(name, features = NULL, samples = NULL) {
      v <- private$.config$assay_map[[name]]
      m <- if (is.null(v)) {
        # Computed assay
        switch(name,
          normalised_log2 = log2(
            self$assay("normalised") +
              (private$.config$background_expression %||% 10)
          ),
          stop("Unknown computed assay: ", name)
        )
      } else {
        SummarizedExperiment::assay(private$.se, v)
      }
      # Remap rownames via fid_lookup. Identity for genomics; built feature
      # IDs for proteomics. Drop rows whose source rowname has no mapping.
      lk <- private$.fid_lookup
      if (!is.null(lk)) {
        new_rn <- unname(lk[rownames(m)])
        keep <- !is.na(new_rn)
        m <- m[keep, , drop = FALSE]
        rownames(m) <- new_rn[keep]
      }
      if (!is.null(features)) m <- m[features, , drop = FALSE]
      if (!is.null(samples))  m <- m[, samples,  drop = FALSE]
      m
    },

    enrichment = function(contrast, kind) {
      stopifnot(kind %in% c("ora", "gsea"))
      loc <- private$.config$enrichment_locator[[kind]]$se
      if (is.null(loc)) return(list())
      obj <- S4Vectors::metadata(private$.se)[[loc]]
      if (is.null(obj)) return(list())
      poc_flatten_cp_results(obj)
    },

    enrichment_input = function(contrast) {
      ei <- S4Vectors::metadata(private$.se)$enrichInput
      if (is.null(ei)) {
        return(list(seq_anno = NULL,
                    ranked_l2fc = NULL,
                    selections = list()))
      }
      list(
        seq_anno    = ei$seqAnno,
        ranked_l2fc = ei$log2Ratio,
        selections  = ei$selections %||% list()
      )
    }
  ),

  private = list(
    .se          = NULL,
    .features_df = NULL,
    .de_dfs      = list(),
    .fid_lookup  = NULL,    # NULL → identity inside assay()

    # Subclasses must override these:
    .build_config = function(se) {
      stop("subclass of ExploreDatasetSE must implement .build_config(se)")
    },
    .prep = function() {
      stop("subclass of ExploreDatasetSE must implement .prep()")
    }
  )
)


# =============================================================================
# 3a. Concrete subclass: ExploreDatasetSE_RNASeq
# =============================================================================
# Genomics RNA-seq SE produced by SUSHI ezApps. Top-level rowData has
# feature annotations + DE columns inline. metadata(se)$param drives most
# config. Optional pre-computed enrichResult / enrichResultGSEA in metadata.

ExploreDatasetSE_RNASeq <- R6::R6Class(
  classname = "ExploreDatasetSE_RNASeq",
  inherit   = ExploreDatasetSE,

  private = list(

    .build_config = function(se) {
      param <- S4Vectors::metadata(se)$param %||% list()
      cd    <- data.frame(SummarizedExperiment::colData(se), check.names = FALSE, stringsAsFactors = FALSE)

      factor_cols  <- grep("\\[Factor]", colnames(cd), value = TRUE)
      factor_names <- gsub(" \\[.*", "", factor_cols)
      factor_names <- gsub(" ",      "_", factor_names)
      factor_levels <- setNames(
        lapply(factor_cols, function(fc) {
          v <- cd[[fc]]
          if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
        }),
        factor_names
      )

      comp <- param$comparison %||% "test--over--reference"
      parts <- strsplit(comp, "--over--", fixed = TRUE)[[1]]
      test_grp <- parts[1] %||% "test"
      ref_grp  <- parts[2] %||% "reference"
      contrast_name <- paste0(.poc_safe(test_grp), "_vs_", .poc_safe(ref_grp))

      cm <- list()
      cm[[contrast_name]] <- list(
        test_group      = test_grp,
        reference_group = ref_grp,
        factor_name     = param$groupingName,
        method          = "DESeq2",
        thresholds      = list(
          fdr    = param$fdrThreshORA   %||% 0.05,
          log2fc = param$log2RatioThreshGO %||% 0,
          go_p   = param$pValThreshGO   %||% 0.01
        )
      )

      enrichment_kinds <- character()
      md <- S4Vectors::metadata(se)
      if (isTRUE(param$runGO)) {
        if (!is.null(md$enrichResult))     enrichment_kinds <- c(enrichment_kinds, "ora")
        if (!is.null(md$enrichResultGSEA)) enrichment_kinds <- c(enrichment_kinds, "gsea")
      }

      norm_key <- intersect(c("signal", "xNorm", "normalizedCounts"),
                            SummarizedExperiment::assayNames(se))
      norm_key <- if (length(norm_key)) norm_key[1] else NA_character_

      list(
        schema_version        = "exploreDE/config/1.0-poc",
        data_type             = "rna_seq",
        organism              = poc_infer_organism(param),
        ref_build             = param$refBuild %||% NA_character_,
        feature_id_kind       = param$featureLevel %||% "gene",
        factors               = factor_names,
        factor_levels         = factor_levels,
        primary_grouping      = param$groupingName %||%
          (if (length(factor_names)) factor_names[[1]] else NA_character_),
        assay_map             = list(
          counts          = "counts",
          normalised      = norm_key,
          normalised_log2 = NULL                  # computed
        ),
        background_expression = param$backgroundExpression %||% 10,
        raw_assay_name        = "counts",
        primary_assay_name    = if (!is.na(norm_key)) "normalised_log2" else "counts",
        display_name_field    = "gene_name",
        description_field     = "description",
        primary_contrast      = contrast_name,
        contrast_meta         = cm,
        enrichment_kinds      = enrichment_kinds,
        enrichment_locator    = list(
          ora  = list(se = "enrichResult"),
          gsea = list(se = "enrichResultGSEA")
        ),
        sgrna_to_gene         = NULL,
        proteomics            = NULL,
        source = list(
          backend_kind   = "se",
          schema_version = "exploreDE/config/1.0-poc",
          built_at       = Sys.time()
        )
      )
    },

    .prep = function() {
      se  <- private$.se
      cfg <- private$.config

      rd <- data.frame(SummarizedExperiment::rowData(se), check.names = FALSE, stringsAsFactors = FALSE)
      rownames(rd)        <- rownames(se)
      rd$feature_id       <- rownames(rd)
      rd$display_name     <- if (cfg$display_name_field %in% colnames(rd)) {
        rd[[cfg$display_name_field]]
      } else {
        rownames(rd)
      }
      blank <- is.na(rd$display_name) | !nzchar(rd$display_name)
      rd$display_name[blank] <- rownames(rd)[blank]
      if ("description" %in% colnames(rd)) {
        rd$description <- gsub(" \\[.*", "", as.character(rd$description))
      }
      private$.features_df <- rd

      # RNA-seq has DE columns inline in rowData (one contrast).
      contrast <- cfg$primary_contrast
      private$.de_dfs <- setNames(
        list(poc_standardise_de_columns(rd)),
        contrast
      )

      # rownames(se) → rownames(se): identity lookup so assay() leaves
      # rownames untouched.
      private$.fid_lookup <-
        setNames(rownames(se), rownames(se))
    }
  )
)


# =============================================================================
# 3b. Concrete subclass: ExploreDatasetSE_CountQC
# =============================================================================
# Count QC SE — top-level rowData with feature annotations, no DE results.
# May include a "presentFlag" / "isPresent" assay for the QC denominator.

ExploreDatasetSE_CountQC <- R6::R6Class(
  classname = "ExploreDatasetSE_CountQC",
  inherit   = ExploreDatasetSE,

  private = list(

    .build_config = function(se) {
      param <- S4Vectors::metadata(se)$param %||% list()
      cd    <- data.frame(SummarizedExperiment::colData(se), check.names = FALSE, stringsAsFactors = FALSE)

      factor_cols  <- grep("\\[Factor]", colnames(cd), value = TRUE)
      factor_names <- gsub(" \\[.*", "", factor_cols)
      factor_names <- gsub(" ",      "_", factor_names)
      factor_levels <- setNames(
        lapply(factor_cols, function(fc) {
          v <- cd[[fc]]
          if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
        }),
        factor_names
      )

      ass_nms <- SummarizedExperiment::assayNames(se)
      norm_key <- intersect(c("signal", "xNorm"), ass_nms)
      norm_key <- if (length(norm_key)) norm_key[1] else NULL
      pf_key   <- intersect(c("presentFlag", "isPresent"), ass_nms)
      pf_key   <- if (length(pf_key)) pf_key[1] else NULL

      list(
        schema_version        = "exploreDE/config/1.0-poc",
        data_type             = "count_qc",
        organism              = poc_infer_organism(param),
        ref_build             = param$refBuild %||% NA_character_,
        feature_id_kind       = param$featureLevel %||% "gene",
        factors               = factor_names,
        factor_levels         = factor_levels,
        primary_grouping      = param$groupingName %||%
          (if (length(factor_names)) factor_names[[1]] else NA_character_),
        assay_map             = c(
          list(counts = "counts"),
          if (!is.null(norm_key)) list(normalised   = norm_key),
          if (!is.null(pf_key))   list(present_flag = pf_key)
        ),
        background_expression = param$backgroundExpression %||% 10,
        raw_assay_name        = "counts",
        primary_assay_name    = "counts",
        display_name_field    = "gene_name",
        description_field     = "description",
        primary_contrast      = NA_character_,
        contrast_meta         = list(),
        enrichment_kinds      = character(),
        enrichment_locator    = list(),
        sgrna_to_gene         = NULL,
        proteomics            = NULL,
        source = list(
          backend_kind   = "se",
          schema_version = "exploreDE/config/1.0-poc",
          built_at       = Sys.time()
        )
      )
    },

    .prep = function() {
      se  <- private$.se
      cfg <- private$.config

      rd <- data.frame(SummarizedExperiment::rowData(se), check.names = FALSE, stringsAsFactors = FALSE)
      rownames(rd)    <- rownames(se)
      rd$feature_id   <- rownames(rd)
      rd$display_name <- if (cfg$display_name_field %in% colnames(rd)) {
        rd[[cfg$display_name_field]]
      } else {
        rownames(rd)
      }
      blank <- is.na(rd$display_name) | !nzchar(rd$display_name)
      rd$display_name[blank] <- rownames(rd)[blank]
      private$.features_df <- rd

      private$.de_dfs     <- list()                      # no DE
      private$.fid_lookup <- setNames(rownames(se), rownames(se))
    }
  )
)


# =============================================================================
# 3c. Concrete subclass: ExploreDatasetSE_CrisprCountQC
# =============================================================================
# Identical to CountQC plus a sgRNA → gene mapping built from the feature
# IDs (assumed to follow GENE_<n>_<SEQUENCE> convention).

ExploreDatasetSE_CrisprCountQC <- R6::R6Class(
  classname = "ExploreDatasetSE_CrisprCountQC",
  inherit   = ExploreDatasetSE_CountQC,

  private = list(

    .build_config = function(se) {
      cfg <- super$.build_config(se)
      cfg$data_type       <- "crispr_count_qc"
      cfg$feature_id_kind <- "guide"
      rn <- rownames(se) %||% character()
      cfg$sgrna_to_gene <- setNames(
        gsub("_(sgrna)?[0-9]+_[ACGTN]+$", "", rn, ignore.case = TRUE),
        rn
      )
      cfg
    }
    # .prep() inherited from CountQC unchanged
  )
)


# =============================================================================
# 3d. Concrete subclass: ExploreDatasetSE_Proteomics
# =============================================================================
# Proteomics SE produced by the FGCZ proteomics pipeline. Top-level
# rowData(se) is sparse — feature metadata and DE stats live nested inside
# rowData(se)$constrast_<name> per-contrast DataFrames. The first contrast
# is the source of truth for feature metadata; feature_id is BUILT via
# build_feature_id (protein_Id [+ site / peptide_Id]). Assay rownames need
# remapping via .fid_lookup since the SE rownames don't match the built
# feature_ids.

ExploreDatasetSE_Proteomics <- R6::R6Class(
  classname = "ExploreDatasetSE_Proteomics",
  inherit   = ExploreDatasetSE,

  private = list(

    .build_config = function(se) {
      cd <- data.frame(SummarizedExperiment::colData(se), check.names = FALSE, stringsAsFactors = FALSE)

      # Factors: exclusion regex (everything that isn't a sample id or
      # technical column).
      excluded     <- grep("^Sample$|^Name$|^names$|^channel$|^raw\\.file$",
                           colnames(cd), value = TRUE, ignore.case = TRUE)
      factor_cols  <- setdiff(colnames(cd), excluded)
      # G_ prefix: drop, blank → "Group"
      factor_names <- sub("^G_", "", factor_cols)
      factor_names[!nzchar(factor_names)] <- "Group"
      factor_levels <- setNames(
        lapply(factor_cols, function(fc) {
          v <- cd[[fc]]
          if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
        }),
        factor_names
      )

      # Contrasts from rowData ^constrast_ columns.
      rd_nms    <- names(SummarizedExperiment::rowData(se))
      contrasts <- gsub("^constrast_", "",
                        grep("^constrast_", rd_nms, value = TRUE))

      # Inspect first contrast to determine flavour (protein/phospho/peptide).
      feature_id_kind <- "protein"
      if (length(contrasts)) {
        df1 <- poc_standardise_rowdata_columns(
          data.frame(SummarizedExperiment::rowData(se)[[paste0("constrast_", contrasts[1])]], check.names = FALSE, stringsAsFactors = FALSE)
        )
        if      ("site"       %in% colnames(df1)) feature_id_kind <- "phospho"
        else if ("peptide_Id" %in% colnames(df1)) feature_id_kind <- "peptide"
      }

      cm <- setNames(
        lapply(contrasts, function(cn) {
          grp <- strsplit(cn, "_vs_|--over--")[[1]]
          list(
            test_group      = grp[1] %||% NA_character_,
            reference_group = grp[2] %||% NA_character_,
            factor_name     = factor_names[1] %||% NA_character_,
            method          = "prolfqua",
            thresholds      = list(fdr = 0.05, log2fc = 1)
          )
        }),
        contrasts
      )

      ass_nms     <- SummarizedExperiment::assayNames(se)
      primary_key <- intersect(c("transformedData", "transformed",
                                 "intensities", "log2"), ass_nms)
      primary_key <- if (length(primary_key)) primary_key[1] else ass_nms[1]
      raw_key     <- intersect(c("intensities_raw", "raw", "intensities"),
                               ass_nms)
      raw_key     <- if (length(raw_key)) raw_key[1] else primary_key

      list(
        schema_version        = "exploreDE/config/1.0-poc",
        data_type             = "proteomics",
        organism              = NA_character_,
        ref_build             = NA_character_,
        feature_id_kind       = feature_id_kind,
        factors               = factor_names,
        factor_levels         = factor_levels,
        primary_grouping      = factor_names[1] %||% NA_character_,
        assay_map             = list(
          intensities     = primary_key,
          intensities_raw = if (raw_key != primary_key) raw_key else NULL
        ),
        background_expression = 1,
        raw_assay_name        = "intensities",
        primary_assay_name    = "intensities",
        display_name_field    = NA_character_,    # handled in .prep()
        description_field     = "description",
        primary_contrast      = if (length(contrasts)) contrasts[1] else NA_character_,
        contrast_meta         = cm,
        enrichment_kinds      = character(),
        enrichment_locator    = list(),
        sgrna_to_gene         = NULL,
        proteomics            = NULL,             # filled by .prep()
        source = list(
          backend_kind   = "se",
          schema_version = "exploreDE/config/1.0-poc",
          built_at       = Sys.time()
        )
      )
    },

    .prep = function() {
      se        <- private$.se
      contrasts <- self$contrasts()

      if (!length(contrasts)) {
        private$.features_df <- data.frame(feature_id   = character(),
                                           display_name = character())
        private$.de_dfs      <- list()
        private$.fid_lookup  <- setNames(rownames(se), rownames(se))
        return(invisible(NULL))
      }

      # Extract and standardise each contrast DataFrame. Build feature_id
      # from the first contrast (all share the same row order at extract
      # time).
      raw_dfs <- setNames(
        lapply(contrasts, function(cn) {
          data.frame(SummarizedExperiment::rowData(se)[[paste0("constrast_", cn)]], check.names = FALSE, stringsAsFactors = FALSE)
        }),
        contrasts
      )
      std_dfs <- lapply(raw_dfs, poc_standardise_rowdata_columns)

      # Reference DataFrame (first contrast, NOT yet sorted by p-value)
      # so feature_id ↔ rownames(se) alignment is preserved for the
      # fid_lookup.
      ref_df <- std_dfs[[1]]
      rownames(ref_df) <- rownames(se)[seq_len(nrow(ref_df))]
      ref_df$feature_id <- poc_build_feature_id(ref_df)

      # fid_lookup: SE rowname → built feature_id
      private$.fid_lookup <- setNames(ref_df$feature_id, rownames(ref_df))

      # features_df: ref_df with display_name (gene_name → protein_Id →
      # feature_id fallback chain).
      display_priority <- c("gene_name", "protein_Id")
      display_priority <- intersect(display_priority, colnames(ref_df))
      display <- rep(NA_character_, nrow(ref_df))
      for (col in display_priority) {
        blank <- is.na(display) | !nzchar(display)
        if (!any(blank)) break
        vals <- as.character(ref_df[[col]])
        display[blank] <- vals[blank]
      }
      blank <- is.na(display) | !nzchar(display)
      display[blank] <- ref_df$feature_id[blank]
      ref_df$display_name <- display

      # Tidy description if present.
      if ("description" %in% colnames(ref_df)) {
        ref_df$description <- gsub(" \\[.*", "",
                                   as.character(ref_df$description))
      }
      rownames(ref_df) <- ref_df$feature_id
      private$.features_df <- ref_df

      # de_dfs: standardised columns + feature_id (built, not raw rownames),
      # sorted by p-value where possible.
      de_dfs <- setNames(
        lapply(seq_along(contrasts), function(i) {
          d <- std_dfs[[i]]
          # Use the same feature_ids as ref_df (per-row alignment to SE).
          d$feature_id <- ref_df$feature_id[seq_len(nrow(d))]
          d <- poc_standardise_de_columns(d)
          if ("p_value" %in% colnames(d)) {
            d <- d[order(d$p_value, na.last = TRUE), , drop = FALSE]
          }
          rownames(d) <- d$feature_id
          d
        }),
        contrasts
      )
      private$.de_dfs <- de_dfs

      # gene_lookup / protein_lookup (proteomics_choices accessor).
      sa1 <- de_dfs[[1]]
      gene_lookup    <- NULL
      gene_choices   <- NULL
      protein_lookup <- NULL
      protein_choices <- NULL
      if ("gene_name" %in% colnames(sa1) && !identical(sa1$feature_id, sa1$gene_name)) {
        gene_lookup  <- split(sa1$feature_id, sa1$gene_name)
        gene_choices <- unique(sa1$gene_name)
      }
      if ("protein_Id" %in% colnames(sa1)) {
        protein_lookup  <- split(sa1$feature_id, sa1$protein_Id)
        protein_choices <- unique(sa1$protein_Id)
      }
      private$.config$proteomics <- list(
        gene    = list(choices = gene_choices,    lookup = gene_lookup),
        protein = list(choices = protein_choices, lookup = protein_lookup)
      )
    }
  )
)


# =============================================================================
# 4. Helpers used in subclasses
# =============================================================================

.poc_safe <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)


# =============================================================================
# 5. Factory: new_explore_dataset()
# =============================================================================

new_explore_dataset <- function(input, config = NULL) {
  se <- if (inherits(input, "SummarizedExperiment")) {
    input
  } else if (is.character(input) && length(input) == 1L) {
    poc_load(input)
  } else {
    stop("Unrecognised input for new_explore_dataset()")
  }

  cls <- switch(detect_se_shape(se),
    rna_seq         = ExploreDatasetSE_RNASeq,
    count_qc        = ExploreDatasetSE_CountQC,
    crispr_count_qc = ExploreDatasetSE_CrisprCountQC,
    proteomics      = ExploreDatasetSE_Proteomics,
    stop("Unknown SE shape")
  )
  cls$new(se, config = config)
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
        if (grepl("\\.qs2$", fname)) {
          return(qs2::qs_read(fp, nthreads = 4))
        }
        return(readRDS(fp))
      }
    }
    stop("No known SE artefact in directory: ", input)
  }
  stop("Cannot interpret input: ", input)
}


# =============================================================================
# 6. Demo runner
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
  # prot_peptide =
  #   "https://fgcz-ms.uzh.ch/public/pStore/p41085/bfabric/Proteomics/exploreDE/2026/2026-03/2026-03-10/workunit_342065/3106968.rds"
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
  show_cols <- intersect(c("feature_id", "display_name", "gene_id",
                           "gene_name", "protein_Id", "description"),
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
    show <- intersect(c("feature_id", "display_name", "gene_name",
                        "protein_Id", "log2_fc", "p_value", "p_adj"),
                      c(colnames(de_df), "feature_id"))
    if (!"feature_id" %in% colnames(de_df))
      de_df$feature_id <- rownames(de_df)
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
        cat("  first element class: ",
            paste(class(e[[1]]), collapse = "/"), "\n", sep = "")
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
