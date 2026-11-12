#' classify_volcano_status
#' The significance class of each Volcano point, from its y value (`-log10`
#' p-value, or a score) and log2 ratio: significant means `log10p >=
#' y_threshold`, a fold change means `|log2ratio| >= lfc`. A missing value is
#' not significant. With `lfc = 0` a zero ratio counts as down.
#' @param log10p Plotted y values.
#' @param log2ratio Log2 ratios.
#' @param y_threshold Significance threshold on the y scale.
#' @param lfc Log2 fold-change threshold.
#' @return Character vector: `FoldChangeSignificantUp/Down`, `FoldChangeUp/Down`,
#'   `SignificantUp/Down` or `NotSignificant`.
#' @noRd
classify_volcano_status <- function(log10p, log2ratio, y_threshold, lfc) {
  dplyr::case_when(
    log10p >= y_threshold & log2ratio <= -lfc ~ "FoldChangeSignificantDown",
    log10p >= y_threshold & log2ratio >= lfc ~ "FoldChangeSignificantUp",
    log10p < y_threshold & log2ratio <= -lfc ~ "FoldChangeDown",
    log10p < y_threshold & log2ratio >= lfc ~ "FoldChangeUp",
    log10p >= y_threshold & log2ratio < 0 ~ "SignificantDown",
    log10p >= y_threshold & log2ratio > 0 ~ "SignificantUp",
    .default = "NotSignificant"
  )
}
