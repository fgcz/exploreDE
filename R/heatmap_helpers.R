# Heatmap helpers (server-heatmap2.R): count heatmap annotation and rows, table export, sample ratios.

#' heatmap_export_table
#' The table a heatmap's "Download ... Counts" button writes: the plotted
#' values per feature and, when asked for, the contrast's DE statistics.
#' @param counts The plotted matrix, features x samples, feature ids as rownames.
#' @param de The DE table (with `feature_id`), or `NULL` when there is none.
#' @param include_de Whether to add `log2Ratio`, `pValue` and `fdr`
#'   ("Include DE statistics").
#' @return A data frame: `feature_id`, the sample columns, then the statistics
#'   the DE table has.
#' @noRd
heatmap_export_table <- function(counts, de, include_de) {
  out <- tibble::rownames_to_column(as.data.frame(counts), "feature_id")
  stats <- if (isTRUE(include_de) && !is.null(de)) intersect(c("log2Ratio", "pValue", "fdr"), colnames(de))
  if (length(stats) == 0) {
    return(out)
  }
  dplyr::left_join(out, de[, c("feature_id", stats), drop = FALSE], by = "feature_id")
}

#' heatmap_term_matrix
#' The rows of a GO term's or KEGG pathway's genes in a heatmap matrix: the
#' genes the matrix has, in their order, with an absolute row median of at
#' least `min_abs_median` (no filter when `NULL`). Always a matrix, also with
#' one row or none.
#' @param counts Heatmap matrix, features x columns.
#' @param genes Feature ids of the term or pathway.
#' @param min_abs_median Minimum absolute row median, or `NULL`.
#' @return The matching rows of `counts`.
#' @noRd
heatmap_term_matrix <- function(counts, genes, min_abs_median = NULL) {
  m <- counts[intersect(genes, rownames(counts)), , drop = FALSE]
  if (!is.null(min_abs_median)) {
    m <- m[which(abs(matrixStats::rowMedians(m)) >= min_abs_median), , drop = FALSE]
  }
  m
}

#' heatmap_top_annotation
#' The sample annotation above the count heatmaps: one track per factor, each
#' level in the colour given. With every colour given, ComplexHeatmap draws no
#' random colours.
#' @param samples Sample table, one row per heatmap column.
#' @param factors Columns of `samples` to show.
#' @param colours Named list by factor of colours named by level.
#' @param legend_param `annotation_legend_param`, a list by factor.
#' @param gp Cell border `gpar()`.
#' @param show_legend Whether ComplexHeatmap draws the legends.
#' @return A `HeatmapAnnotation`.
#' @noRd
heatmap_top_annotation <- function(samples, factors, colours, legend_param, gp, show_legend = TRUE) {
  do.call(
    ComplexHeatmap::HeatmapAnnotation,
    c(
      lapply(stats::setNames(factors, factors), function(f) samples[[f]]),
      list(
        col = colours[factors],
        show_legend = show_legend,
        annotation_legend_param = legend_param,
        gp = gp
      )
    )
  )
}

#' select_genes_by_sd
#' Pick the top-N most variable features based on row-wise standard deviation.
#' Useful when ranking CountQC features without DE statistics.
#' @param counts Numeric matrix (features x samples)
#' @param top_n Number of features to keep (default 2000)
#' @return Character vector of rownames for the top features
#' @noRd
select_genes_by_sd <- function(counts, top_n = 2000) {
  gene_sds <- apply(counts, 1, stats::sd, na.rm = TRUE)
  top_genes <- names(utils::head(sort(gene_sds, decreasing = TRUE), n = top_n))
  return(top_genes)
}

#' sample_log2_ratio
#' The log2 ratio of two sample sets: per feature, the mean of the numerator
#' samples minus the mean of the denominator samples, on a log2-scale matrix.
#' Missing values are left out of the means.
#' @param mat Numeric matrix (features x samples) on the log2 scale
#' @param numerator,denominator Column names of `mat`
#' @return Named numeric vector, one value per row of `mat`
#' @noRd
sample_log2_ratio <- function(mat, numerator, denominator) {
  unknown <- setdiff(c(numerator, denominator), colnames(mat))
  if (length(unknown)) {
    stop("not samples of the count matrix: ", paste(unknown, collapse = ", "))
  }
  num <- rowMeans(mat[, numerator, drop = FALSE], na.rm = TRUE)
  den <- rowMeans(mat[, denominator, drop = FALSE], na.rm = TRUE)
  stats::setNames(num - den, rownames(mat))
}

#' ratio_heatmap_matrix
#' The features x ratios matrix of the Ratios heatmap: among the features with
#' a finite value in every ratio, the `top_n` with the largest variance across
#' the ratios or the largest absolute ratio. A single ratio is always ranked by
#' absolute value. The matrix has no missing values, so it can be clustered.
#' @param ratios Named list of named numeric vectors, one per ratio
#' @param top_n Number of features to keep
#' @param rank_by `"variance"` or `"absolute"`
#' @return Numeric matrix (features x ratios); 0 rows when nothing is left
#' @noRd
ratio_heatmap_matrix <- function(ratios, top_n, rank_by = c("variance", "absolute")) {
  rank_by <- match.arg(rank_by)
  ratios <- Filter(Negate(is.null), ratios)
  if (!length(ratios)) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  features <- Reduce(intersect, lapply(ratios, names))
  m <- do.call(cbind, lapply(ratios, function(r) unname(r[features])))
  m <- matrix(as.numeric(m), nrow = length(features), ncol = length(ratios))
  dimnames(m) <- list(if (length(features)) features else NULL, names(ratios))
  m <- m[rowSums(!is.finite(m)) == 0, , drop = FALSE]
  if (!nrow(m)) {
    return(m)
  }
  score <- if (rank_by == "variance" && ncol(m) > 1) {
    apply(m, 1, stats::var)
  } else {
    apply(abs(m), 1, max)
  }
  keep <- utils::head(order(score, decreasing = TRUE), top_n)
  m[keep, , drop = FALSE]
}

#' contrast_log2_ratio
#' One contrast's log2 fold change as a named vector, offered on the Ratios
#' heatmap next to the sample ratios.
#' @param eds An `ExploreDataset`
#' @param contrast A contrast name, one of `eds$contrasts()`
#' @return Named numeric vector (by feature id), or `NULL` when the contrast's
#'   DE table has no `log2Ratio`
#' @noRd
contrast_log2_ratio <- function(eds, contrast) {
  de <- eds$de(contrast)
  if (!all(c("feature_id", "log2Ratio") %in% colnames(de))) {
    return(NULL)
  }
  stats::setNames(de$log2Ratio, de$feature_id)
}
