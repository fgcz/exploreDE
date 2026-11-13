# Feature ids and table columns: standardised DE columns, proteomics feature ids, guides to genes.

#' standardise_columns
#' Fuzzy-match incoming column names to canonical names.
#' Upstream data may call columns anything; this maps them to a known schema
#' so downstream code can rely on consistent names. Patterns come from the one
#' column dictionary (`inst/dictionaries/column_patterns.json`); DE concepts
#' are renamed to the legacy names the SE-backed app code reads
#' (`log2Ratio` / `pValue` / `fdr`). A column is only renamed when the
#' deciding pattern matches exactly one column.
#' @param df A data.frame (typically from rowData contrast)
#' @return The same data.frame with standardised column names
#' @noRd
standardise_columns <- function(df) {
  # dictionary concept -> legacy column name it is renamed to
  rowdata_concepts <- c(
    "protein_Id",
    "peptide_Id",
    "site",
    "gene_name",
    "description",
    "fasta.id",
    "modelName",
    "nrPeptides"
  )
  de_concepts <- c(effect_log2fc = "log2Ratio", score_pvalue = "pValue", score_fdr = "fdr")

  orig <- colnames(df)
  current <- orig
  rename <- function(current, patterns, to) {
    src <- eds_match_col(orig, patterns, unique = TRUE)
    if (!is.null(src)) {
      current[match(src, orig)] <- to
    }
    current
  }
  for (concept in rowdata_concepts) {
    current <- rename(current, eds_rowdata_column_patterns[[concept]], concept)
  }
  for (concept in names(de_concepts)) {
    current <- rename(current, eds_de_column_patterns[[concept]], de_concepts[[concept]])
  }
  colnames(df) <- current
  df
}

#' build_feature_id
#' Build a unique, human-readable feature identifier for proteomics data.
#' Uses canonical column names (after standardisation).
#' @param df A standardised data.frame with canonical column names
#' @return Character vector of unique feature IDs
#' @noRd
build_feature_id <- function(df) {
  if (!"protein_Id" %in% colnames(df)) {
    stop(
      "build_feature_id(): no protein_Id column; columns are: ",
      paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }
  # A protein_Id can have gaps (SAINT): use the cleaned row name there, as the
  # all-rows fallback below does, so no feature gets an NA id.
  clean_rownames <- gsub("~[^~]+~", "_", rownames(df))
  pid <- as.character(df$protein_Id)
  missing_pid <- is.na(pid) | !nzchar(pid)
  pid[missing_pid] <- clean_rownames[missing_pid]

  # If protein_Id is already unique, use it directly
  if (!anyDuplicated(pid)) {
    return(pid)
  }

  # Try known sub-feature columns (canonical names) in priority order
  # Always prepend protein_Id so users see which protein a sub-feature belongs to
  for (col in intersect(c("site", "peptide_Id"), colnames(df))) {
    combined <- paste(pid, gsub("~", "_", df[[col]]), sep = "_")
    if (!anyDuplicated(combined)) return(combined)
  }

  # Fallback: cleaned raw rownames
  clean_rownames
}

#' collapse_guides_to_genes
#' Collapse sgRNA-level count matrix to gene-level using row sums.
#' Filters out negative control genes (matching ^negCtrl).
#' @param counts A numeric matrix with sgRNA rows and sample columns
#' @param sgRNA_to_gene Named character vector: names are sgRNA feature IDs, values are gene names
#' @return A numeric matrix with gene-level rows
#' @noRd
collapse_guides_to_genes <- function(counts, sgRNA_to_gene, exclude = NULL) {
  common <- intersect(rownames(counts), names(sgRNA_to_gene))
  counts <- counts[common, , drop = FALSE]
  # A logical assay (present flag) sums to the number of guides present per gene
  if (is.logical(counts)) {
    storage.mode(counts) <- "integer"
  }
  collapsed <- rowsum(counts, group = sgRNA_to_gene[common], reorder = FALSE)
  if (length(exclude) > 0) {
    collapsed <- collapsed[!rownames(collapsed) %in% exclude, , drop = FALSE]
  }
  collapsed
}

#' Collapse an assay to gene level for download, respecting its scale
#'
#' Linear assays are summed per gene. Log2-scale assays are summed on the
#' linear scale and logged again (`log2(sum(2^x))`), because a sum of logs is
#' not a meaningful gene-level value.
#' @param m Numeric matrix, sgRNA rows and sample columns.
#' @param sgRNA_to_gene Named character vector mapping sgRNA ids to genes.
#' @param is_log2 Whether `m` is on the log2 scale.
#' @return A numeric matrix with gene rows.
#' @noRd
collapse_assay_to_genes <- function(m, sgRNA_to_gene, is_log2 = FALSE) {
  m <- as.matrix(m)
  if (!isTRUE(is_log2)) {
    return(collapse_guides_to_genes(m, sgRNA_to_gene))
  }
  log2(collapse_guides_to_genes(2^m, sgRNA_to_gene))
}
