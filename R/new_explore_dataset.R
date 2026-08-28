#' Detect the shape of a `SummarizedExperiment` for subclass dispatch
#'
#' Classifies a SUSHI / FGCZ-proteomics SE into one of four shapes:
#' `"rna_seq"`, `"count_qc"`, `"crispr_count_qc"`, or `"proteomics"`.
#' Decision is made from `metadata(se)$param$comparison`,
#' `rowData(se)` column names, and `rownames(se)` patterns. Internal
#' helper for [new_explore_dataset()]; exported so tests and producer-
#' side migrations can verify dispatch logic.
#'
#' @param se A `SummarizedExperiment`.
#' @return Character scalar: one of `"rna_seq"`, `"count_qc"`,
#'   `"crispr_count_qc"`, `"proteomics"`.
#' @export
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
  # Guide feature IDs, two known conventions, both ending in the protospacer:
  #   GENE_<n>_<SEQ>                       (Mageck: GENE_1_ACGT… / GENE_sgrna1_ACGT…)
  #   GENE_<strand>_<locus>_<SEQ>          (e.g. Ptplad2_-_88439147.23-P1P2_GATCC…)
  if (any(grepl("_(sgrna)?[0-9]+_[ACGTN]+$", rn, ignore.case = TRUE)) ||
      any(grepl("_[+-]_[^_]+_[ACGTN]+$", rn, ignore.case = TRUE))) {
    return("crispr_count_qc")
  }
  "count_qc"
}

#' Construct an `ExploreDataset` instance from a path, URL, or SE
#'
#' Factory function: loads the input (or accepts a live
#' `SummarizedExperiment`), classifies its shape via
#' [detect_se_shape()], and returns the appropriate concrete subclass.
#'
#' Accepted inputs:
#' - A live `SummarizedExperiment` object.
#' - A character path to an `.rds` / `.qs2` file.
#' - A character path to a `.RData` file (Mageck EzResult convention).
#' - A character path to a directory containing one of the standard
#'   SUSHI artefacts (`deResult.rds`, `rawData_final.rds`, etc.).
#' - A character path to a parquet directory written by
#'   [export_h5ad_to_parquet()] (the exploreDElight bridge format).
#' - An `http(s)://` URL pointing at an `.rds` file.
#'
#' @param input See above.
#' @param config Optional pre-built config list; advanced use only.
#' @return An object inheriting from [ExploreDataset].
#' @export
new_explore_dataset <- function(input, config = NULL) {
  # Canonical AnnData input (a live anndataR object or a `.h5ad` path) routes
  # straight to the H5AD backend — the "bring your own data" path.
  if (inherits(input, "AbstractAnnData") ||
      (is.character(input) && length(input) == 1L &&
       grepl("\\.h5ad$", input, ignore.case = TRUE))) {
    return(ExploreDatasetH5AD$new(input, config = config))
  }

  # A parquet directory written by export_h5ad_to_parquet() (or by the Rust
  # bridge in exploreDElight). Checked BEFORE the SE branch, because that one
  # also accepts a directory — eds_load_se() would go looking for deResult.rds
  # and report a confusing "no SUSHI artefact" rather than the real shape.
  if (is.character(input) && length(input) == 1L && dir.exists(input) &&
      file.exists(file.path(input, "manifest.json")) &&
      file.exists(file.path(input, "obs.parquet"))) {
    return(ExploreDatasetParquet$new(input, config = config))
  }

  se <- if (inherits(input, "SummarizedExperiment")) {
    input
  } else if (is.character(input) && length(input) == 1L) {
    eds_load_se(input)
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
