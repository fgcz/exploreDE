# QC scatter plots of the Count QC and DE QC tabs: drawing, colour-by classes, signal alignment.

#' make_qc_scatter
#' Build a QC scatter plot in the style of ezRun's ezXYScatter.2:
#' log10 axes, three typed colour layers, an optional 1:1 solid guide line,
#' and optional N-fold dashed guide lines. Used by the DE QC and Count QC
#' tabs for between-group and intra-group scatter plots.
#' @param x Numeric vector (linear scale, NOT log2)
#' @param y Numeric vector (linear scale, same length as x)
#' @param types_vec Character vector of category labels, same length as x/y
#' @param types_colours Named character vector mapping level -> colour
#' @param xlab,ylab Axis labels
#' @param title Plot title (NULL for none)
#' @param point_size,alpha Point aesthetics
#' @param base_size Font size passed to ggprism::theme_prism
#' @param draw_order Optional character vector of level names defining z-order
#'   (last level = drawn on top). Defaults to names(types_colours).
#' @param guide_colour Colour for the 1:1 and fold-change guide lines.
#'   Default "blue".
#' @param guide_fold Fold-change threshold for the dashed guide lines.
#'   Default 2 (gives ±log10(2) dashed lines).
#' @param show_identity Whether to draw the solid 1:1 identity line.
#'   Default TRUE.
#' @param show_fold_bands Whether to draw the dashed fold-change guide lines.
#'   Default TRUE.
#' @return ggplot object
#' @noRd
make_qc_scatter <- function(
  x,
  y,
  types_vec,
  types_colours,
  xlab,
  ylab,
  title = NULL,
  point_size = 1,
  alpha = 0.6,
  base_size = 12,
  draw_order = NULL,
  guide_colour = "blue",
  guide_fold = 2,
  show_identity = TRUE,
  show_fold_bands = TRUE,
  log_scale = TRUE
) {
  # log_scale = TRUE: x/y are LINEAR values (e.g. RNA-seq group-mean counts)
  #   plotted on log10 axes; a `guide_fold`-fold band is log10(guide_fold) apart.
  # log_scale = FALSE: x/y are ALREADY in log space (e.g. proteomics log2
  #   intensities, centred near 0, can be negative) — plot on linear axes and
  #   space the fold band by log2(guide_fold). Applying log10 to already-log
  #   data drops negatives and mislocates the identity/fold guides.
  if (is.null(draw_order)) {
    draw_order <- names(types_colours)
  }
  df <- tibble::tibble(x = x, y = y, types = types_vec)
  df$types <- factor(df$types, levels = draw_order)
  df <- df[order(df$types), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = .data$types)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::scale_color_manual(values = types_colours)
  if (isTRUE(log_scale)) {
    p <- p + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
  }

  if (isTRUE(show_identity)) {
    p <- p + ggplot2::geom_abline(intercept = 0, colour = guide_colour)
  }
  if (isTRUE(show_fold_bands) && is.numeric(guide_fold) && guide_fold > 1) {
    off <- if (isTRUE(log_scale)) log10(guide_fold) else log2(guide_fold)
    p <- p +
      ggplot2::geom_abline(intercept = off, colour = guide_colour, linetype = "dashed") +
      ggplot2::geom_abline(intercept = -off, colour = guide_colour, linetype = "dashed")
  }

  p +
    ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1)
      )
    ) +
    ggplot2::labs(x = xlab, y = ylab, title = title)
}

#' qc_types_significance
#' Colour-by classes "Standard" of the DE QC and Count QC scatters: Present /
#' Absent by `present_mask`, and Significants (p below `p_thresh`, |log2
#' ratio| at least `lfc_thresh`) when `rd` has `pValue` and `log2Ratio`.
#' @param rd Feature table aligned with the scatter points.
#' @param present_mask Logical per feature (RNA-seq: usedInTest; Count QC:
#'   presentFlag mean).
#' @param p_thresh,lfc_thresh Significance thresholds.
#' @param cols List with `Absent`, `Present`, `Significants` colours.
#' @return `list(types, colours, draw_order)` for make_qc_scatter().
#' @noRd
qc_types_significance <- function(
  rd,
  present_mask,
  p_thresh = 0.05,
  lfc_thresh = 1,
  cols = list(Absent = "grey70", Present = "black", Significants = "#A6CEE3")
) {
  types <- ifelse(present_mask, "Present", "Absent")
  if (all(c("log2Ratio", "pValue") %in% colnames(rd))) {
    is_sig <- !is.na(rd$pValue) & rd$pValue < p_thresh & !is.na(rd$log2Ratio) & abs(rd$log2Ratio) >= lfc_thresh
    types[is_sig] <- "Significants"
  }
  list(
    types = types,
    colours = c(Absent = cols$Absent, Present = cols$Present, Significants = cols$Significants),
    draw_order = c("Absent", "Present", "Significants")
  )
}

#' qc_types_binned
#' Colour-by classes of the QC scatters binning a feature value (GC content,
#' gene width) below `lo`, between, and above `hi`; a missing value counts as
#' between. Labels read "<prefix>< lo<suffix>", "lo-hi<suffix>",
#' "<prefix>> hi<suffix>", with `lo` and `hi` formatted by `fmt`.
#' @param values Value per feature.
#' @param lo,hi Bin limits.
#' @param cols List with `low`, `mid`, `high` colours.
#' @param fmt `sprintf()` format of the limits.
#' @param prefix,suffix Label prefix ("GC ") and suffix (" nt").
#' @return `list(types, colours, draw_order)` for make_qc_scatter().
#' @noRd
qc_types_binned <- function(
  values,
  lo,
  hi,
  cols = list(low = "#1F78B4", mid = "black", high = "#FF7F00"),
  fmt = "%.2f",
  prefix = "",
  suffix = ""
) {
  labels <- c(
    low = paste0(prefix, "< ", sprintf(fmt, lo), suffix),
    mid = paste0(sprintf(fmt, lo), "-", sprintf(fmt, hi), suffix),
    high = paste0(prefix, "> ", sprintf(fmt, hi), suffix)
  )
  types <- rep(labels[["mid"]], length(values))
  types[!is.na(values) & values < lo] <- labels[["low"]]
  types[!is.na(values) & values > hi] <- labels[["high"]]
  list(
    types = types,
    colours = stats::setNames(c(cols$low, cols$mid, cols$high), labels),
    draw_order = unname(labels[c("mid", "low", "high")])
  )
}

#' transform_feature_id
#' Mirror the character substitutions the RNA-seq backend applies to
#' isoform-level feature ids (replacing ":", "/", "(", ")",
#' "-", and spaces with "."). Gene-level feature_ids are unaffected. Use
#' this before matching seqAnno$feature_id against signal rownames.
#' @param x Character vector of feature_ids
#' @return Character vector with init-compatible substitutions applied
#' @noRd
transform_feature_id <- function(x) {
  gsub(" ", ".", gsub("[:/()-]", ".", x), fixed = FALSE)
}

#' align_signal_to_seqanno
#' Drop-through alignment of a count/signal matrix to a seqAnno data frame
#' by feature_id. Falls back to the isoform-style transformation when the
#' direct match leaves NAs. Returns the matrix reordered to the seqAnno rows
#' that successfully matched, along with the subset of seqAnno kept.
#' @param signal Count/signal matrix with rownames to align
#' @param rd A data frame with a `feature_id` column
#' @return list(signal = aligned matrix, rd = subset of rd that matched)
#' @noRd
align_signal_to_seqanno <- function(signal, rd) {
  align <- match(rd$feature_id, rownames(signal))
  if (any(is.na(align))) {
    alt <- transform_feature_id(rd$feature_id)
    align_alt <- match(alt, rownames(signal))
    if (sum(!is.na(align_alt)) > sum(!is.na(align))) align <- align_alt
  }
  keep <- !is.na(align)
  list(
    signal = signal[align[keep], , drop = FALSE],
    rd = rd[keep, , drop = FALSE]
  )
}
