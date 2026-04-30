#!/usr/bin/env Rscript
# ============================================================================
#  build_example_anndata.R
#
#  READING THIS FILE IS PART OF THE DOCUMENTATION.
#
#  Purpose
#  -------
#  Convert the exploreDE demo SUSHI outputs into AnnData (.h5ad) files that
#  conform to the exploreDE/1.0 schema defined in
#    exploreDE/plans/anndata_migration/anndata_spec_2026-04-17.qmd
#
#  This script is written as a verbose, narrated walk-through.  Each section
#  explains:
#    (a) WHAT is being read from the source SummarizedExperiment,
#    (b) WHERE it lands in the AnnData object, and
#    (c) WHY the mapping was chosen (with a reference back to the spec).
#
#  When run, it prints the same narrative to the console so you can follow
#  the conversion live:
#
#      module load Dev/R/4.5.0
#      Rscript exploreDE/plans/anndata_migration/build_example_anndata.R [out_dir]
#
#  If you only want to READ the conversion, open this file top-to-bottom —
#  no code execution needed.
#
#  Source inputs (the demo suite)
#  ------------------------------
#      1. Teaching RNA-seq (DESeq2, 1 factor, gene-level)      — p3009
#      2. CRISPR screen CountQC (Mageck .RData, no DE)         — p36911
#      3. CRISPR screen CountQC (Mageck .RData, no DE)         — p29779
#      4. RNA-seq DESeq2, two factors                          — p29779
#      5. RNA-seq DESeq2, transcript/isoform level             — p34429
#      6. Proteomics (prolfqua, protein-level, many contrasts) — p23078
#      7. Proteomics (prolfqua, phospho sites)                 — p39819
#      8. Proteomics (prolfqua, peptide-level)                 — p41085
#
#  Outputs
#  -------
#      <out_dir>/teaching_rna_p3009.h5ad
#      <out_dir>/mageck_p36911.h5ad
#      …
# ============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
})

# --------------------------------------------------------------------------
#  Pretty console helpers — used everywhere below to narrate what is
#  happening at each step.  These are not part of the schema; they exist so
#  a human watching the script run can see the mapping.
# --------------------------------------------------------------------------

.rule     <- function(char = "-", n = 74) message(strrep(char, n))
.h1       <- function(txt) { .rule("="); message(toupper(txt)); .rule("=") }
.h2       <- function(txt) { .rule("-"); message(txt);          .rule("-") }
.note     <- function(...) message("  . ", ...)
.mapline  <- function(src, dst, why = "") {
  # Format: "SE: foo.bar            ->  ADATA: obs['bar']   (why …)"
  message(sprintf("  %-42s ->  %-32s %s",
                  paste0("SE: ",    src),
                  paste0("ADATA: ", dst),
                  if (nchar(why)) paste0("(", why, ")") else ""))
}
.warn     <- function(...) message("  ! ", ...)

# --------------------------------------------------------------------------
#  0. Setup — arguments, output directory, optional dependencies
# --------------------------------------------------------------------------

args     <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(args) >= 1) args[[1]] else "/scratch/peterleary/anndata_examples"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SCHEMA_VERSION <- "exploreDE/1.0"   # see anndata_spec_2026-04-17.qmd §uns$schema_version

.h1("Setup")
.note("schema version : ", SCHEMA_VERSION)
.note("output dir     : ", out_dir)

# anndataR is the R-native AnnData implementation (scverse/anndataR).
# We auto-install via pak if missing — matches the project convention
# (CLAUDE.md: "R packages: Use pak::pkg_install(); faster than install.packages").
if (!requireNamespace("anndataR", quietly = TRUE)) {
  .warn("anndataR not installed — installing via pak from scverse/anndataR")
  if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak", repos = "https://cloud.r-project.org")
  pak::pkg_install("scverse/anndataR")
}
suppressPackageStartupMessages(library(anndataR))

has_cp  <- requireNamespace("clusterProfiler", quietly = TRUE) &&
           requireNamespace("org.Hs.eg.db",     quietly = TRUE)
has_scu <- requireNamespace("scuttle",          quietly = TRUE)
has_jsn <- requireNamespace("jsonlite",         quietly = TRUE)

.note("clusterProfiler + org.Hs.eg.db : ", if (has_cp)  "available" else "missing (enrichment skipped)")
.note("scuttle                         : ", if (has_scu) "available" else "missing (feature_id fallback)")
.note("jsonlite                        : ", if (has_jsn) "available" else "missing (no enrichment_native)")

# --------------------------------------------------------------------------
#  1. The demo suite — one entry per .h5ad we will produce
# --------------------------------------------------------------------------

DATASETS <- list(
  teaching_rna_p3009 =
    "/srv/gstore/projects/p3009/o5638_DESeq2_diff--over--undiff_2024-11-20--12-29-40/diff--over--undiff/",
  mageck_p36911 =
    "/srv/gstore/projects/p36911/o39664_o37783_o39124_o38710_o41353_ExploreMageckCounts_2026-04-09--15-07-24/ExploreMageckCounts/counts-yuunbieawvpw-EzResult.RData",
  mageck_p29779 =
    "/srv/gstore/projects/p29779/o41056_ExploreMageckCounts_2026-04-10--09-29-57/ExploreMageckCounts/counts-fzvwfvxsycwz-EzResult.RData",
  rna_two_factor_p29779 =
    "/srv/gstore/projects/p29779/o40890_DESeq2_shp65_w_HP--over--shp65_wo_HP_2026-04-16--13-38-37/shp65_w_HP--over--shp65_wo_HP",
  rna_isoform_p34429 =
    "/srv/gstore/projects/p34429/o35733_DESeq2_Palbo_day6--over--Palbo_day2_2024-08-22--11-21-44/Palbo_day6--over--Palbo_day2/",
  proteomics_protein_p23078 =
    "https://fgcz-ms.uzh.ch/public/pStore/p23078/bfabric/Proteomics/exploreDE/2026/2026-03/2026-03-10/workunit_342063/3106962.rds",
  proteomics_phospho_p39819 =
    "https://fgcz-ms.uzh.ch/public/pStore/p39819/bfabric/Proteomics/exploreDE/2025/2025-10/2025-10-22/workunit_335485/2997996.rds",
  proteomics_peptide_p41085 =
    "https://fgcz-ms.uzh.ch/public/pStore/p41085/bfabric/Proteomics/exploreDE/2026/2026-03/2026-03-10/workunit_342065/3106968.rds"
)

# --------------------------------------------------------------------------
#  2. Source loader
#
#  SUSHI ships its data in three different envelopes.  The loader hides
#  that heterogeneity — downstream code only ever sees a SummarizedExperiment.
#
#     Genomics DESeq2    : directory  containing  deResult.rds
#     CRISPR  CountQC    : .RData file with ≥1 SE in the top environment
#     Proteomics         : .rds over HTTPS (fgcz-ms.uzh.ch pStore)
#
#  Detaching from SUSHI means the PUBLIC loader will only ever accept a
#  local .h5ad.  For THIS script we're going the other way (SUSHI → h5ad),
#  so we still read the SUSHI envelopes.
# --------------------------------------------------------------------------

slurp_source <- function(path) {
  if (grepl("\\.RData$", path, ignore.case = TRUE)) {
    env <- new.env(); load(path, envir = env)
    # Mageck EzResult files contain several objects (param, rawData, result,
    # sceset, se) — we want the SE.
    if ("se" %in% ls(env)) return(env$se)
    for (n in ls(env)) if (inherits(env[[n]], "SummarizedExperiment")) return(env[[n]])
    stop("No SummarizedExperiment found in ", path)
  }
  if (grepl("^https?://", path))                    return(readRDS(url(path)))
  if (grepl("\\.rds$", path, ignore.case = TRUE))   return(readRDS(path))

  # Directory — look for the DESeq2 convention deResult.rds
  cand <- file.path(path, "deResult.rds")
  if (file.exists(cand)) return(readRDS(cand))
  stop("Unrecognised source: ", path)
}

# Omics-type detection — we consult three signals, in order of reliability:
#   1. metadata(se)$param$appName     (set by SUSHI)
#   2. substrings of the path hint    (SUSHI folder/URL conventions)
#   3. structural clues in rowData    (constrast_ columns → proteomics)
detect_omics_type <- function(se, path_hint) {
  md <- metadata(se)
  if (!is.null(md$param$appName) &&
      grepl("Mageck|CountQC|Crispr", md$param$appName, ignore.case = TRUE)) return("crispr_countqc")
  if (grepl("ExploreMageckCounts|CountQC|Crispr", path_hint, ignore.case = TRUE)) return("crispr_countqc")
  if (grepl("Proteomics|prolfqua|pStore|/bfabric/", path_hint, ignore.case = TRUE)) return("proteomics")
  if (any(grepl("^constrast_|^contrast_", colnames(rowData(se))))) return("proteomics")
  "rnaseq_bulk"
}

# --------------------------------------------------------------------------
#  3. Shared helpers
# --------------------------------------------------------------------------

# SUSHI writes descriptions like "Tumor protein p53 [Source:HGNC;Acc:11998]".
# The " [Source:…]" tail pollutes search results — the current app strips
# it on load, and so does the new schema.
strip_source_tail <- function(x) ifelse(is.na(x), x, sub(" \\[.*", "", x))

# Contrast / factor names must survive as h5ad keys.  The schema reserves
# "[a-z0-9_]+" only (see §open-questions — to be enforced in the validator).
sanitise_name <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  sub("^_+|_+$", "", x)
}

# `%||%`: tiny null-coalesce operator.  Defined before first use.
`%||%` <- function(a, b) if (!is.null(a) && length(a)) a else b

# --------------------------------------------------------------------------
#  FACTOR STANDARDISATION
#
#  SUSHI                         -->  exploreDE/1.0
#  --------------------------------   -----------------------------------
#  colData column  "Condition [Factor]"    obs column  "factor_condition"
#  colData column  "Time point [Factor]"   obs column  "factor_time_point"
#
#  Why: the "[Factor]" suffix is SUSHI-specific and requires string parsing
#  in every downstream callsite.  Moving to a "factor_" PREFIX lets the
#  loader detect factors with one grep, works natively in Python, and
#  survives h5ad round-trips without quoting.
# --------------------------------------------------------------------------

standardise_factor_cols <- function(col_data) {
  df <- as.data.frame(col_data, check.names = FALSE, stringsAsFactors = FALSE)
  factor_orig <- grep("\\[Factor\\]", colnames(df), value = TRUE)
  new_names <- vapply(factor_orig, function(c)
    paste0("factor_", sanitise_name(sub(" \\[Factor\\]", "", c))), character(1))
  for (i in seq_along(factor_orig)) {
    old <- factor_orig[i]; new <- new_names[i]
    colnames(df)[colnames(df) == old] <- new
    if (!is.factor(df[[new]])) df[[new]] <- factor(df[[new]])
    .mapline(paste0("colData$", old), paste0("obs$", new),
             "rename [Factor]-suffix to factor_ prefix")
  }
  list(df = df, factors = unname(new_names))
}

# --------------------------------------------------------------------------
#  PROTEOMICS COLUMN STANDARDISATION
#
#  Proteomics DE tables arrive with wildly varied column names depending on
#  the pipeline (`diff` vs `log2FC` vs `logFC`; `p.value` vs `pValue`; etc.)
#  The app's own `standardise_columns()` fuzzy-matches to a canonical schema
#  BEFORE any other processing.  We reproduce that here so the mapping to
#  the spec's `effect_*` / `score_*` is deterministic.
# --------------------------------------------------------------------------

standardise_prot_columns <- function(df) {
  patterns <- c(
    protein_Id  = "^protein[._]?id$",
    peptide_Id  = "^peptide[._]?id$|^peptide$",
    site        = "^site$",
    gene_name   = "^gene[._]?name$|^gene[._]?symbol$",
    description = "^description$|^desc$",
    fasta.id    = "^fasta[._]?id$",
    modelName   = "^model[._]?name$",
    nrPeptides  = "^nr[._]?peptides$|^n[._]?peptides$|^num[._]?peptides$",
    log2Ratio   = "^diff$|^log2[._]?(fold)?[._]?(change|ratio|fc)?$|^logfc$",
    pValue      = "^p[._]?value$",
    fdr         = "^fdr$|^adj[._]?p[._]?value$|^padj$|^q[._]?value$"
  )
  cur <- colnames(df)
  for (canonical in names(patterns)) {
    idx <- grep(patterns[[canonical]], cur, ignore.case = TRUE)
    if (length(idx) == 1) cur[idx] <- canonical
  }
  colnames(df) <- cur
  df
}

# --------------------------------------------------------------------------
#  FEATURE_ID CONSTRUCTION (proteomics)
#
#  Preserves the 2026-03-12 refactor decision (see feature_id_refactor.qmd):
#
#    * start from protein_Id
#    * append the next-finer hierarchy column (site → peptide → peptidoform
#      → ion) until the resulting vector is unique
#
#  Examples (from real data):
#    Protein-level (p23078):  feature_id = protein_Id                   = "A0A024RBG1"
#    Phospho       (p39819):  feature_id = protein_Id_site              = "A0A2R8Y4L2_S199"
#    Peptide       (p23078):  feature_id = protein_Id_peptide_Id        = "A0A024RBG1_LLGIFEQNQDR"
# --------------------------------------------------------------------------

build_feature_id_prot <- function(df) {
  if (!"protein_Id" %in% colnames(df)) return(make.unique(rownames(df)))
  candidates <- intersect(c("site", "peptide_Id", "peptidoform", "ion"), colnames(df))
  fid <- as.character(df$protein_Id)
  if (!anyDuplicated(fid)) return(fid)
  for (col in candidates) {
    fid <- paste(df$protein_Id, df[[col]], sep = "_")
    if (!anyDuplicated(fid)) return(fid)
  }
  make.unique(fid)
}

# --------------------------------------------------------------------------
#  DE TABLE CANONICALISATION
#
#  The spec mandates three columns on every varm["de_<contrast>"] table:
#
#      effect_log2fc  (float)   -- volcano x-axis
#      score_pvalue   (float)   -- volcano y-axis
#      score_fdr      (float)   -- significance threshold
#
#  SUSHI convention                -->  exploreDE/1.0
#  ------------------------------       -----------------------
#  rowData$log2Ratio                    varm[de_<c>]$effect_log2fc
#  rowData$pValue                       varm[de_<c>]$score_pvalue
#  rowData$fdr                          varm[de_<c>]$score_fdr
#  rowData$nrPeptides                   varm[de_<c>]$label_nr_peptides
#  rowData$modelName                    varm[de_<c>]$label_model_name
#  rowData$usedInTest                   varm[de_<c>]$label_used_in_test
#
#  Per-assay group means (<assay>_<group>_Mean) the current app computes on
#  load are deliberately DROPPED here — the spec says those are derived on
#  the fly from layers + uns$de[<c>] when needed, not baked into varm.
# --------------------------------------------------------------------------

canonicalise_de <- function(df) {
  out <- data.frame(row.names = rownames(df))
  if ("log2Ratio"  %in% colnames(df)) out$effect_log2fc     <- as.numeric(df$log2Ratio)
  if ("pValue"     %in% colnames(df)) out$score_pvalue      <- as.numeric(df$pValue)
  if ("fdr"        %in% colnames(df)) out$score_fdr         <- as.numeric(df$fdr)
  if ("nrPeptides" %in% colnames(df)) out$label_nr_peptides <- as.numeric(df$nrPeptides)
  for (n in c("modelName", "fasta.id", "usedInTest"))
    if (n %in% colnames(df)) out[[paste0("label_", sanitise_name(n))]] <- df[[n]]
  out
}

# --------------------------------------------------------------------------
#  ENRICHMENT
#
#  Spec decision (see §`uns["enrichment"]`):
#     - First-class requirement ACROSS ALL OMICS, not genomics only.
#     - Reference implementation is clusterProfiler (enrichResult /
#       gseaResult).
#     - Store both NATIVE (base64-serialised R object → uns$enrichment_native)
#       AND TIDY (data.frame → uns$enrichment).
#       R loader prefers native; Python readers read tidy.
#     - If the upstream pipeline already ran clusterProfiler (genomics does,
#       proteomics does not), we lift the pre-computed objects.
#     - Otherwise the exporter runs clusterProfiler from scratch here.
# --------------------------------------------------------------------------

species_from_refbuild <- function(refBuild) {
  if (is.null(refBuild)) return(NA_character_)
  rb <- tolower(as.character(refBuild))
  if (grepl("homo_sapiens|grch|hg19|hg38", rb)) return("Hs")
  if (grepl("mus_musculus|grcm|mm10|mm39",   rb)) return("Mm")
  if (grepl("rattus|rn6|rn7",                rb)) return("Rn")
  if (grepl("drosophila|dm6",                rb)) return("Dm")
  NA_character_
}

orgdb_for <- function(species_code) {
  db <- switch(species_code %||% "", Hs = "org.Hs.eg.db", Mm = "org.Mm.eg.db",
                                     Rn = "org.Rn.eg.db", Dm = "org.Dm.eg.db",
                                     NA_character_)
  if (is.na(db) || !requireNamespace(db, quietly = TRUE)) return(NULL)
  db
}

tidy_enrich <- function(obj, method) {
  if (is.null(obj)) return(NULL)
  df <- tryCatch(as.data.frame(obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  # clusterProfiler columns --> spec tidy columns
  ren <- c(ID = "id", Description = "description", GeneRatio = "gene_ratio",
           BgRatio = "bg_ratio", pvalue = "pvalue", p.adjust = "p_adjust",
           qvalue = "q_value", Count = "n_genes_mapped",
           setSize = "n_genes_in_set", enrichmentScore = "enrichment_score",
           NES = "nes", geneID = "feature_ids", core_enrichment = "feature_ids",
           ONTOLOGY = "ontology")
  for (old in names(ren))
    if (old %in% colnames(df)) colnames(df)[colnames(df) == old] <- ren[[old]]
  if ("feature_ids" %in% colnames(df))
    df$feature_ids <- strsplit(as.character(df$feature_ids), "/", fixed = TRUE)
  df$method   <- method
  df$category <- df$ontology %||% NA_character_
  df
}

run_cp_for_contrast <- function(varm_de, feature_ids, entrez_map, species_code,
                                fdr_thresh = 0.1, l2fc_thresh = 1) {
  if (!has_cp || is.null(entrez_map) || is.null(orgdb_for(species_code))) {
    .warn("clusterProfiler or OrgDb unavailable — skipping enrichment")
    return(list(native = list(), tidy = list()))
  }
  orgdb <- orgdb_for(species_code)
  df <- varm_de
  df$entrez <- entrez_map[match(feature_ids, names(entrez_map))]
  df <- df[!is.na(df$entrez) & !is.na(df$effect_log2fc) & !is.na(df$score_fdr), , drop = FALSE]
  if (!nrow(df)) return(list(native = list(), tidy = list()))
  up   <- df$entrez[df$effect_log2fc >  l2fc_thresh & df$score_fdr < fdr_thresh]
  down <- df$entrez[df$effect_log2fc < -l2fc_thresh & df$score_fdr < fdr_thresh]
  safely <- function(expr) tryCatch(expr, error = function(e) { .warn("cP: ", conditionMessage(e)); NULL })

  .note(sprintf("clusterProfiler: %d up, %d down, %d ranked genes (species=%s)",
                length(up), length(down), nrow(df), species_code))

  ora_up <- if (length(up)   >= 5) safely(clusterProfiler::enrichGO(gene = up,   OrgDb = orgdb, ont = "BP", readable = TRUE)) else NULL
  ora_dn <- if (length(down) >= 5) safely(clusterProfiler::enrichGO(gene = down, OrgDb = orgdb, ont = "BP", readable = TRUE)) else NULL
  ranked <- setNames(df$effect_log2fc, df$entrez)
  ranked <- sort(ranked[!duplicated(names(ranked))], decreasing = TRUE)
  gsea   <- if (length(ranked) >= 50)
              safely(clusterProfiler::gseGO(ranked, OrgDb = orgdb, ont = "BP", verbose = FALSE)) else NULL

  list(
    native = Filter(Negate(is.null), list(ora_up = ora_up, ora_down = ora_dn, gsea = gsea)),
    tidy   = Filter(Negate(is.null), list(
      ora_up   = tidy_enrich(ora_up, "ora_up"),
      ora_down = tidy_enrich(ora_dn, "ora_down"),
      gsea     = tidy_enrich(gsea,   "gsea")))
  )
}

# --------------------------------------------------------------------------
#  4. RNA-seq converter   (DESeq2 / edgeR / limma bulk)
#
#  SE → AnnData mapping (full):
#
#    assay(counts)             -->  layers$counts    (+ X)     (features × samples → samples × features  via t())
#    assay(signal) / xNorm     -->  layers$normalised          (axis flipped)
#    assay(presentFlag)        -->  layers$present_flag        (optional)
#    derived log2(norm + bg)   -->  layers$log2                (bg from param$backgroundExpression)
#
#    colData$Condition[Factor] -->  obs$factor_condition       (rename)
#    rownames(colData)         -->  obs_names                  (same)
#    param$featureLevel        -->  obs$feature_level          (copied per sample)
#    param$refBuild            -->  obs$organism (head), uns$reference_genome
#
#    rowData$gene_id           -->  var$label_gene_id
#    rowData$gene_name         -->  var$label_gene_symbol
#    rowData$transcript_id     -->  var$label_transcript_id
#    rowData$description       -->  var$description            (Source tail stripped)
#    rowData$type              -->  var$label_biotype
#    rowData$GO BP/MF/CC       -->  var$label_go_bp/mf/cc      (col names normalised)
#    rowData$entrez_id         -->  var$label_entrez_id
#    scuttle::uniquify(gene_id, gene_name)  -->  var_names (= feature_id)
#
#    rowData$log2Ratio/pValue/fdr  -->  varm$de_<c>$effect_log2fc / score_pvalue / score_fdr
#    (single contrast per SE for genomics — name from param$comparison)
#
#    param$sampleGroup / refGroup  -->  uns$de[<c>]$test_group / reference_group
#    param$comparison              -->  contrast KEY of uns$de AND varm key
#    metadata$method               -->  uns$de[<c>]$model
#
#    metadata$enrichResult     -->  uns$enrichment_native$<c>$ora_all
#                                   uns$enrichment$<c>$ora_all               (tidy mirror)
#    metadata$enrichResultGSEA -->  uns$enrichment_native$<c>$gsea
#                                   uns$enrichment$<c>$gsea                  (tidy mirror)
#
#    metadata$enrichInput      -->  DROPPED (derived on demand from varm + var)
#    param (everything else)   -->  DROPPED (SUSHI-specific; not needed)
# --------------------------------------------------------------------------

convert_rnaseq <- function(se, path_hint, omics_type = "rnaseq_bulk") {
  md    <- metadata(se)
  param <- md$param %||% list()

  # -- OBS (samples × obs-columns) ----------------------------------------
  .h2("obs (sample metadata)")
  fac <- standardise_factor_cols(colData(se))
  obs <- fac$df
  obs$sample_label <- rownames(obs)
  .mapline("rownames(colData)", "obs_names", "sample identifiers")
  .mapline("(new)", "obs$sample_label", "= obs_names by default")
  if (!is.null(param$featureLevel)) {
    obs$feature_level <- param$featureLevel
    .mapline("param$featureLevel", "obs$feature_level",
             paste0("value = '", param$featureLevel, "'"))
  }
  if (!is.null(param$refBuild)) {
    obs$organism <- sub("/.*", "", param$refBuild)
    .mapline("param$refBuild (head)", "obs$organism",
             paste0("value = '", obs$organism[[1]], "'"))
  }

  # -- VAR (features × var-columns) ---------------------------------------
  .h2("var (feature annotation)")
  rd <- as.data.frame(rowData(se), check.names = FALSE, stringsAsFactors = FALSE)
  feature_level <- param$featureLevel %||% "gene"

  # feature_id: preserves the 2026-03-12 refactor decision.  var_names in
  # the resulting AnnData is this column.
  fid <- if (has_scu && all(c("gene_id", "gene_name") %in% colnames(rd))) {
    if (feature_level == "isoform" && "transcript_id" %in% colnames(rd)) {
      .mapline("rowData$transcript_id + gene_name",
               "var_names (= feature_id)", "uniquifyFeatureNames (isoform)")
      scuttle::uniquifyFeatureNames(rd$transcript_id, rd$gene_name)
    } else {
      .mapline("rowData$gene_id + gene_name",
               "var_names (= feature_id)", "uniquifyFeatureNames (gene)")
      scuttle::uniquifyFeatureNames(rd$gene_id, rd$gene_name)
    }
  } else {
    .warn("scuttle unavailable — falling back to rownames for feature_id")
    make.unique(rownames(rd))
  }

  var <- data.frame(row.names = fid,
                    description = strip_source_tail(rd$description %||% rep(NA, nrow(rd))),
                    stringsAsFactors = FALSE)
  .mapline("rowData$description", "var$description", "[Source:...] tail stripped")

  lab_map <- c(gene_id = "label_gene_id", gene_name = "label_gene_symbol",
               transcript_id = "label_transcript_id",
               type = "label_biotype", biotypes = "label_biotypes",
               entrez_id = "label_entrez_id",
               `GO BP` = "label_go_bp", `GO MF` = "label_go_mf", `GO CC` = "label_go_cc",
               `GO.BP` = "label_go_bp", `GO.MF` = "label_go_mf", `GO.CC` = "label_go_cc",
               gc = "label_gc", featWidth = "label_feat_width")
  for (old in names(lab_map)) {
    if (old %in% colnames(rd)) {
      var[[lab_map[[old]]]] <- rd[[old]]
      .mapline(paste0("rowData$", old), paste0("var$", lab_map[[old]]), "identifier / annotation")
    }
  }

  hierarchy <- if (feature_level == "isoform")
    list(levels = c("gene", "transcript"), unit = "transcript") else
    list(levels = "gene", unit = "gene")
  .note("uns$hierarchy = { levels: [", paste(hierarchy$levels, collapse = ", "),
        "], unit: '", hierarchy$unit, "' }")

  # -- LAYERS (axis flip!) ------------------------------------------------
  .h2("X and layers (abundance matrices)")
  .note("AXIS FLIP: SE is features×samples; AnnData is samples×features. We t()")
  assay_map <- c(counts = "counts", xNorm = "normalised", signal = "normalised",
                 presentFlag = "present_flag", isPresent = "present_flag")
  layers <- list()
  for (a in assayNames(se)) {
    key <- assay_map[a] %||% a
    m <- as.matrix(assay(se, a)); rownames(m) <- fid
    layers[[key]] <- t(m)
    .mapline(paste0("assay(se, '", a, "')"),
             paste0("layers$", key),
             sprintf("%dx%d  ->  %dx%d  (t())",
                     nrow(m), ncol(m), nrow(layers[[key]]), ncol(layers[[key]])))
  }
  if ("normalised" %in% names(layers) && !"log2" %in% names(layers)) {
    bg <- param$backgroundExpression %||% 10
    layers[["log2"]] <- log2(layers[["normalised"]] + bg)
    .mapline("derived log2(normalised + bg)", "layers$log2",
             sprintf("bg = param$backgroundExpression = %s", bg))
  }
  X <- layers[["counts"]] %||% layers[[1]]
  .note("X = layers$counts (raw counts — canonical for rnaseq_bulk per spec)")

  # -- VARM + uns$de (per-contrast DE statistics) -------------------------
  .h2("varm / uns$de (differential expression)")
  contrast_name <- sanitise_name(param$comparison %||% "de")
  varm   <- list()
  uns_de <- list()
  de_cols <- intersect(c("log2Ratio", "pValue", "fdr", "usedInTest", "modelName"), colnames(rd))
  if (length(de_cols)) {
    de_df <- rd[, de_cols, drop = FALSE]; rownames(de_df) <- fid
    varm[[paste0("de_", contrast_name)]] <- canonicalise_de(de_df)
    uns_de[[contrast_name]] <- list(
      assay_used       = "counts",
      factor_used      = fac$factors,
      test_group       = as.character(param$sampleGroup %||% NA),
      reference_group  = as.character(param$refGroup    %||% NA),
      contrast_formula = param$comparison %||% NA_character_,
      model            = md$method %||% "DESeq2",
      n_test           = sum(obs[[fac$factors[[1]]]] == (param$sampleGroup %||% ""), na.rm = TRUE),
      n_reference      = sum(obs[[fac$factors[[1]]]] == (param$refGroup    %||% ""), na.rm = TRUE)
    )
    .note(sprintf("contrast '%s': test=%s (n=%d), reference=%s (n=%d)",
                  contrast_name,
                  uns_de[[contrast_name]]$test_group,      uns_de[[contrast_name]]$n_test,
                  uns_de[[contrast_name]]$reference_group, uns_de[[contrast_name]]$n_reference))
    .mapline(paste(de_cols, collapse = " / "),
             paste0("varm$de_", contrast_name, "$effect_log2fc / score_pvalue / score_fdr / label_*"),
             "canonicalise column names")
    .mapline("param$sampleGroup / refGroup",
             paste0("uns$de$", contrast_name, "$test_group / reference_group"),
             "replaces SUSHI-specific param")
  } else {
    .warn("no DE columns (log2Ratio/pValue/fdr) found in rowData — varm stays empty")
  }

  # -- ENRICHMENT ---------------------------------------------------------
  .h2("uns$enrichment + uns$enrichment_native")
  enr_native <- list(); enr_tidy <- list()
  if (!is.null(md$enrichResult) || !is.null(md$enrichResultGSEA)) {
    .note("pre-computed clusterProfiler objects present in metadata — lifting")
    native <- Filter(Negate(is.null),
                     list(ora_all = md$enrichResult, gsea = md$enrichResultGSEA))
    tidy <- Filter(Negate(is.null),
                   list(ora_all = tidy_enrich(md$enrichResult,     "ora_all"),
                        gsea    = tidy_enrich(md$enrichResultGSEA, "gsea")))
    enr_native[[contrast_name]] <- native
    enr_tidy[[contrast_name]]   <- tidy
    .mapline("metadata$enrichResult",     paste0("uns$enrichment_native$", contrast_name, "$ora_all"), "native R object (base64 in h5ad)")
    .mapline("metadata$enrichResultGSEA", paste0("uns$enrichment_native$", contrast_name, "$gsea"),    "native R object (base64 in h5ad)")
    .mapline("(tidy mirror)",              paste0("uns$enrichment$", contrast_name, "$ora_all / gsea"),"portable data.frame")
  } else if (length(varm)) {
    species <- species_from_refbuild(param$refBuild)
    entrez  <- if ("label_entrez_id" %in% colnames(var))
                 setNames(as.character(var$label_entrez_id), rownames(var)) else NULL
    if (!is.null(entrez)) {
      .note("no pre-computed enrichment — running clusterProfiler from varm")
      out <- run_cp_for_contrast(varm[[paste0("de_", contrast_name)]],
                                 rownames(var), entrez, species)
      if (length(out$native)) enr_native[[contrast_name]] <- out$native
      if (length(out$tidy))   enr_tidy[[contrast_name]]   <- out$tidy
    }
  }

  # -- assembly -----------------------------------------------------------
  list(X = X, layers = layers, obs = obs, var = var, varm = varm,
       uns = list(
         schema_version    = SCHEMA_VERSION,
         omics_type        = omics_type,
         hierarchy         = hierarchy,
         design            = list(factors = fac$factors,
                                  grouping = param$grouping %||% fac$factors[[1]] %||% NA),
         de                = uns_de,
         enrichment        = enr_tidy,
         enrichment_native = enr_native,
         reference_genome  = as.character(param$refBuild %||% NA),
         provenance        = list(
           created_utc       = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
           source_format     = "SummarizedExperiment",
           source_path       = path_hint,
           converter_version = "build_example_anndata.R/0.1")
       ))
}

# --------------------------------------------------------------------------
#  5. CRISPR CountQC converter  (Mageck .RData)
#
#  Special features vs generic RNA-seq:
#
#    * No DE.  CountQC is a quality-control view over sgRNA counts — there
#      are no log2Ratio / pValue / fdr columns in rowData.  varm stays empty;
#      uns$de stays empty.
#    * Features are sgRNAs, not genes.  rownames look like
#      "CXCR6_1_CCCGCCCCTGTGTAACGCGT" = "<target_gene>_<guide_num>_<sequence>".
#      We parse these into explicit var columns.
#    * assay(presentFlag) is critical for the present-gene filter.
#      In the current app it drives server-summary.R library-size calcs.
#
#    rowData rownames       -->  parsed into  var$label_target_gene,
#                                              var$label_sgrna_number,
#                                              var$label_sgrna_sequence
#    assay(counts)          -->  layers$counts (= X)
#    assay(signal)          -->  layers$normalised
#    assay(presentFlag)     -->  layers$present_flag
#    colData$Condition[F]   -->  obs$factor_condition       (same rule as RNA)
#    uns$hierarchy          -->  { levels: ["gene","sgrna"], unit: "sgrna" }
# --------------------------------------------------------------------------

convert_crispr_countqc <- function(se, path_hint) {
  md    <- metadata(se)
  param <- md$param %||% list()

  .h2("obs (sample metadata)")
  fac <- standardise_factor_cols(colData(se))
  obs <- fac$df; obs$sample_label <- rownames(obs); obs$feature_level <- "sgrna"
  .mapline("(implicit)", "obs$feature_level", "value = 'sgrna'")

  .h2("var (feature annotation)")
  rd <- as.data.frame(rowData(se), check.names = FALSE, stringsAsFactors = FALSE)
  parts <- strsplit(rownames(rd), "_", fixed = TRUE)
  target_gene <- vapply(parts, `[`, character(1), 1)
  guide_num   <- vapply(parts, function(p) if (length(p) >= 2) p[[2]] else NA_character_, character(1))
  sg_seq      <- vapply(parts, function(p) if (length(p) >= 3) paste(p[-(1:2)], collapse = "_") else NA_character_, character(1))
  fid <- make.unique(rownames(rd))
  .mapline("rowData rowname (parsed)",  "var$label_target_gene",     "gene_symbol — first underscore-separated token")
  .mapline("rowData rowname (parsed)",  "var$label_sgrna_number",    "second token (guide index)")
  .mapline("rowData rowname (parsed)",  "var$label_sgrna_sequence",  "tokens 3..N")
  .mapline("rowData rownames",          "var_names (= feature_id)",  "unique sgRNA id")

  var <- data.frame(row.names = fid,
                    description = strip_source_tail(rd$description %||% rep(NA, nrow(rd))),
                    label_target_gene    = target_gene,
                    label_sgrna_number   = guide_num,
                    label_sgrna_sequence = sg_seq,
                    stringsAsFactors = FALSE)
  lab_map <- c(gene_name = "label_gene_symbol", entrez_id = "label_entrez_id",
               type = "label_biotype", gc = "label_gc",
               isValid = "label_is_valid", isPresentProbe = "label_is_present")
  for (old in names(lab_map))
    if (old %in% colnames(rd)) {
      var[[lab_map[[old]]]] <- rd[[old]]
      .mapline(paste0("rowData$", old), paste0("var$", lab_map[[old]]), "annotation")
    }

  .h2("X and layers")
  assay_map <- c(counts = "counts", signal = "normalised", xNorm = "normalised",
                 presentFlag = "present_flag", isPresent = "present_flag")
  layers <- list()
  for (a in assayNames(se)) {
    key <- assay_map[a] %||% a
    m <- as.matrix(assay(se, a)); rownames(m) <- fid
    layers[[key]] <- t(m)
    .mapline(paste0("assay(se, '", a, "')"),
             paste0("layers$", key),
             sprintf("%dx%d flipped", nrow(m), ncol(m)))
  }
  X <- layers[["counts"]] %||% layers[[1]]
  .note("X = layers$counts (raw sgRNA counts)")
  .note("uns$de = {} and varm = {} — CountQC ships no DE contrasts")

  list(X = X, layers = layers, obs = obs, var = var, varm = list(),
       uns = list(
         schema_version    = SCHEMA_VERSION,
         omics_type        = "crispr_countqc",
         hierarchy         = list(levels = c("gene", "sgrna"), unit = "sgrna"),
         design            = list(factors = fac$factors,
                                  grouping = param$grouping %||% fac$factors[[1]] %||% NA),
         de                = list(),
         enrichment        = list(),
         enrichment_native = list(),
         reference_genome  = as.character(param$refBuild %||% NA),
         provenance = list(
           created_utc       = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
           source_format     = "SummarizedExperiment (Mageck CountQC)",
           source_path       = path_hint,
           converter_version = "build_example_anndata.R/0.1")
       ))
}

# --------------------------------------------------------------------------
#  6. Proteomics converter   (prolfqua)
#
#  Key differences from genomics:
#
#   (a) MANY contrasts per object.  Genomics has ONE (param$comparison);
#       proteomics has N, stored in SE rowData as NESTED DataFrames under
#       columns named "constrast_<name>" — note the SUSHI typo, preserved
#       here so existing data loads.  Each contrast becomes one
#       varm["de_<name>"] + one uns$de[<name>] entry.
#
#   (b) No [Factor] markup.  Proteomics colData has plain-named columns
#       (e.g. "Group", "Treatment") and the app treats everything that
#       isn't Sample/Name/channel/raw.file as a factor.  We do the same.
#
#   (c) Hierarchy columns move OUT of the DE tables and INTO shared var
#       columns (hierarchy_protein / hierarchy_peptide / hierarchy_site) —
#       annotation belongs to the feature, not to the contrast.
#
#   (d) No pre-computed enrichment.  We run clusterProfiler from scratch
#       via UniProt → Entrez mapping (bitr).
#
#  SE → AnnData mapping:
#
#    colData$Sample/Name/sampleName  -->  obs_names
#    colData$<everything else>       -->  obs$factor_<sanitised>
#
#    rowData$constrast_<c> (DataFrame) with cols:
#        protein_Id / peptide_Id / site / description / gene_name / ...
#                                    --> (shared) --> var$description,
#                                                     var$hierarchy_*,
#                                                     var$label_gene_symbol,
#                                                     var$label_uniprot_id,
#                                                     var$label_nr_peptides
#        log2Ratio / pValue / fdr / nrPeptides / modelName
#                                    --> varm$de_<c>$effect_log2fc / score_pvalue
#                                                 / score_fdr / label_nr_peptides
#
#    assay(transformedData)     -->  X                   (analysis-ready)
#    assay(rawData)             -->  layers$intensities_raw
#    assay(<others>)            -->  layers$<sanitised>
# --------------------------------------------------------------------------

convert_proteomics <- function(se, path_hint) {
  md  <- metadata(se)

  .h2("obs (sample metadata)")
  obs_raw <- as.data.frame(colData(se), check.names = FALSE, stringsAsFactors = FALSE)
  sample_name_col <- grep("^Sample$|^Name$|^sampleName$", colnames(obs_raw),
                          value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(sample_name_col)) {
    rownames(obs_raw) <- obs_raw[[sample_name_col]]
    .mapline(paste0("colData$", sample_name_col), "obs_names",
             "proteomics uses Name / Sample as the sample identifier")
  }
  keep <- grep("sample|name|names|channel|raw\\.file",
               colnames(obs_raw), invert = TRUE, ignore.case = TRUE, value = TRUE)
  obs <- data.frame(row.names = rownames(obs_raw),
                    sample_label = rownames(obs_raw),
                    stringsAsFactors = FALSE)
  factor_cols <- character(0)
  for (col in keep) {
    new <- paste0("factor_", sanitise_name(col))
    obs[[new]] <- factor(obs_raw[[col]])
    factor_cols <- c(factor_cols, new)
    .mapline(paste0("colData$", col), paste0("obs$", new),
             "proteomics factor — treated as categorical")
  }

  # rowData → VAR + per-contrast VARM ------------------------------------
  .h2("var + varm (contrasts are nested in rowData)")
  rd <- rowData(se)
  contrast_cols <- grep("^constrast_|^contrast_", colnames(rd), value = TRUE)
  if (!length(contrast_cols)) stop("No constrast_* columns in proteomics rowData")
  .note(length(contrast_cols), " contrasts found: ", paste(contrast_cols, collapse = ", "))

  first <- as.data.frame(rd[[contrast_cols[1]]], check.names = FALSE, stringsAsFactors = FALSE)
  first <- standardise_prot_columns(first)
  rownames(first) <- if (is.null(rownames(rd))) make.unique(rownames(first)) else rownames(rd)
  fid   <- build_feature_id_prot(first)
  .mapline("rowData[[first contrast]]$protein_Id (+ site/peptide if needed)",
           "var_names (= feature_id)",
           "build_feature_id_prot() — preserves 2026-03-12 refactor")

  var <- data.frame(row.names = fid,
                    description = strip_source_tail(first$description %||% rep(NA, nrow(first))),
                    stringsAsFactors = FALSE)
  if ("protein_Id" %in% colnames(first)) {
    var$hierarchy_protein <- first$protein_Id
    .mapline("rowData[[c]]$protein_Id", "var$hierarchy_protein", "shared across contrasts")
  }
  if ("peptide_Id" %in% colnames(first)) {
    var$hierarchy_peptide <- first$peptide_Id
    .mapline("rowData[[c]]$peptide_Id", "var$hierarchy_peptide", "sub-feature lineage")
  }
  if ("site" %in% colnames(first)) {
    var$hierarchy_site <- first$site
    .mapline("rowData[[c]]$site", "var$hierarchy_site", "PTM site lineage")
  }
  if ("gene_name"  %in% colnames(first)) { var$label_gene_symbol <- first$gene_name;
    .mapline("rowData[[c]]$gene_name",  "var$label_gene_symbol", "identifier") }
  if ("protein_Id" %in% colnames(first)) { var$label_uniprot_id  <- first$protein_Id;
    .mapline("rowData[[c]]$protein_Id", "var$label_uniprot_id",  "identifier") }
  if ("nrPeptides" %in% colnames(first)) { var$label_nr_peptides <- first$nrPeptides;
    .mapline("rowData[[c]]$nrPeptides", "var$label_nr_peptides", "annotation") }
  if ("fasta.id"   %in% colnames(first)) { var$label_fasta_id    <- first$fasta.id;
    .mapline("rowData[[c]]$fasta.id",   "var$label_fasta_id",    "annotation") }

  hier_levels <- intersect(c("protein", "peptide", "site"),
                           sub("^hierarchy_", "", grep("^hierarchy_", colnames(var), value = TRUE)))
  if (!length(hier_levels)) hier_levels <- "protein"
  hierarchy <- list(levels = hier_levels, unit = tail(hier_levels, 1))
  .note("uns$hierarchy = { levels: [", paste(hierarchy$levels, collapse = ", "),
        "], unit: '", hierarchy$unit, "' }")

  # LAYERS ---------------------------------------------------------------
  .h2("X and layers")
  .note("AXIS FLIP: SE features×samples -> AnnData samples×features (t())")
  layers <- list()
  an <- assayNames(se)
  primary <- intersect(c("transformedData", "lfq", "LFQ"), an)[1]
  if (is.na(primary)) primary <- an[1]
  .note(sprintf("primary (X): assay '%s' (analysis-ready intensities)", primary))
  for (a in an) {
    key <- switch(a,
                  rawData         = "intensities_raw",
                  transformedData = "intensities_transformed",
                  lfq = , LFQ     = "lfq",
                  sanitise_name(a))
    m <- as.matrix(assay(se, a)); rownames(m) <- fid
    layers[[key]] <- t(m)
    .mapline(paste0("assay(se, '", a, "')"), paste0("layers$", key),
             sprintf("%dx%d flipped", nrow(m), ncol(m)))
  }
  X <- t(as.matrix(assay(se, primary))); colnames(X) <- fid
  .mapline(paste0("assay(se, '", primary, "')"), "X", "primary analysis matrix")

  # VARM — one per contrast, plus uns$de entries ------------------------
  .h2("varm / uns$de (one entry per contrast)")
  varm   <- list()
  uns_de <- list()
  for (cc in contrast_cols) {
    df <- as.data.frame(rd[[cc]], check.names = FALSE, stringsAsFactors = FALSE)
    df <- standardise_prot_columns(df); rownames(df) <- fid
    de_df <- canonicalise_de(df)
    cname <- sanitise_name(sub("^constrast_|^contrast_", "", cc))
    varm[[paste0("de_", cname)]] <- de_df
    parts <- strsplit(cname, "_vs_", fixed = TRUE)[[1]]
    uns_de[[cname]] <- list(
      assay_used       = primary,
      factor_used      = factor_cols,
      test_group       = if (length(parts) == 2) parts[[1]] else NA_character_,
      reference_group  = if (length(parts) == 2) parts[[2]] else NA_character_,
      contrast_formula = cname,
      model            = md$method %||% "prolfqua"
    )
    .mapline(paste0("rowData$", cc, " (nested DataFrame)"),
             paste0("varm$de_", cname, " + uns$de$", cname),
             "constrast_<c> -> de_<c> (SUSHI typo preserved on input)")
  }

  # ENRICHMENT — proteomics ships none; generate via clusterProfiler ---
  .h2("uns$enrichment + uns$enrichment_native (generated)")
  enr_native <- list(); enr_tidy <- list()
  species <- "Hs"   # demo assumption; real exporter infers from organism metadata
  if (has_cp && !is.null(var$label_uniprot_id)) {
    orgdb <- orgdb_for(species)
    if (!is.null(orgdb)) {
      uni <- sub(";.*", "", as.character(var$label_uniprot_id))
      map <- tryCatch(
        clusterProfiler::bitr(unique(uni), fromType = "UNIPROT",
                              toType = "ENTREZID", OrgDb = orgdb),
        error = function(e) { .warn("bitr failed: ", conditionMessage(e)); NULL })
      if (!is.null(map)) {
        entrez <- setNames(map$ENTREZID, map$UNIPROT)[uni]
        entrez <- setNames(entrez, rownames(var))
        cname  <- sub("^de_", "", names(varm)[1])
        .note("running clusterProfiler on first contrast '", cname, "' as demo")
        out <- run_cp_for_contrast(varm[[paste0("de_", cname)]],
                                   rownames(var), entrez, species)
        if (length(out$native)) enr_native[[cname]] <- out$native
        if (length(out$tidy))   enr_tidy[[cname]]   <- out$tidy
      }
    }
  }

  list(X = X, layers = layers, obs = obs, var = var, varm = varm,
       uns = list(
         schema_version    = SCHEMA_VERSION,
         omics_type        = "proteomics",
         hierarchy         = hierarchy,
         design            = list(factors = factor_cols,
                                  grouping = factor_cols[[1]] %||% NA),
         de                = uns_de,
         enrichment        = enr_tidy,
         enrichment_native = enr_native,
         provenance        = list(
           created_utc       = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
           source_format     = "SummarizedExperiment (prolfqua)",
           source_path       = path_hint,
           converter_version = "build_example_anndata.R/0.1")
       ))
}

# --------------------------------------------------------------------------
#  7. AnnData assembly + h5ad write
#
#  enrichResult / gseaResult are native R S4 objects — h5ad has no schema
#  for those.  We serialise to raw bytes and base64-encode so they travel
#  as strings inside uns.  The R loader decodes them on read; Python
#  readers ignore the "_native" slot and rely on the tidy mirror.
# --------------------------------------------------------------------------

encode_native <- function(obj) {
  raw <- serialize(obj, NULL)
  jsonlite::base64_enc(raw)
}

flatten_native_enrichment <- function(nested) {
  if (!length(nested) || !has_jsn) return(list())
  lapply(nested, function(cl) lapply(cl, function(obj) if (is.null(obj)) NULL else encode_native(obj)))
}

write_anndata <- function(parts, out_path) {
  parts$uns$enrichment_native <- flatten_native_enrichment(parts$uns$enrichment_native)

  # anndataR accepts: X (matrix), obs / var (data.frames with row.names),
  # layers / obsm / varm / uns as named lists.  obs_names / var_names come
  # from the row.names of obs / var respectively.
  adata <- anndataR::AnnData(
    X      = parts$X,
    obs    = parts$obs,
    var    = parts$var,
    layers = parts$layers,
    varm   = parts$varm,
    uns    = parts$uns
  )
  anndataR::write_h5ad(adata, out_path)
  invisible(out_path)
}

summarise_adata <- function(parts, out_path) {
  .h2("summary")
  .note(sprintf("obs     : %d samples x %d columns", nrow(parts$obs), ncol(parts$obs)))
  .note(sprintf("var     : %d features x %d columns", nrow(parts$var), ncol(parts$var)))
  .note(sprintf("layers  : %s", paste(names(parts$layers), collapse = ", ")))
  .note(sprintf("varm    : %s", if (length(parts$varm)) paste(names(parts$varm), collapse = ", ") else "<none>"))
  .note(sprintf("uns$de  : %s", if (length(parts$uns$de)) paste(names(parts$uns$de), collapse = ", ") else "<none>"))
  .note(sprintf("enrich  : %s", if (length(parts$uns$enrichment)) paste(names(parts$uns$enrichment), collapse = ", ") else "<none>"))
  .note(sprintf("wrote   : %s", out_path))
}

# --------------------------------------------------------------------------
#  8. Main dispatch
# --------------------------------------------------------------------------

convert_one <- function(name, path) {
  .h1(paste("converting:", name))
  .note("source: ", path)
  se <- tryCatch(slurp_source(path), error = function(e) { .warn("load failed: ", conditionMessage(e)); NULL })
  if (is.null(se)) return(invisible(NULL))
  omics <- detect_omics_type(se, path)
  .note(sprintf("detected omics_type = '%s'  (SE dim = %s)",
                omics, paste(dim(se), collapse = " x ")))
  parts <- switch(omics,
    rnaseq_bulk    = convert_rnaseq(se, path, "rnaseq_bulk"),
    crispr_countqc = convert_crispr_countqc(se, path),
    proteomics     = convert_proteomics(se, path),
    stop("Unknown omics type: ", omics))
  out_path <- file.path(out_dir, paste0(name, ".h5ad"))
  tryCatch({
    write_anndata(parts, out_path)
    summarise_adata(parts, out_path)
  }, error = function(e) {
    .warn("h5ad write failed: ", conditionMessage(e))
    .warn("(the parts object is still returned in-memory for inspection)")
  })
  invisible(parts)
}

if (!interactive()) {
  for (n in names(DATASETS))
    tryCatch(convert_one(n, DATASETS[[n]]),
             error = function(e) .warn("FAILED [", n, "]: ", conditionMessage(e)))
  .h1("done")
  .note("outputs in: ", out_dir)
}
