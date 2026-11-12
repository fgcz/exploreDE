# Ordination (server-ordination.R): PCA / PCoA / GLM-PCA, biplot data, the preprocessing cascade steps
# and the plots. plot_ordination_samples() and plot_ordination_biplot() read their settings from
# `input` (or a list with the same names).

#' calculate_pca
#' Run base R prcomp on a sample-by-feature matrix and return scores / loadings / variance.
#' @param counts Matrix with samples as rows, features as columns
#' @param center Logical
#' @param scale Logical
#' @return List with scores, loadings, variance, method, component_names, removed_features
#' @noRd
calculate_pca <- function(counts, center = TRUE, scale = FALSE) {
  # Remove features (columns) with missing values; proteomics has NAs in every
  # sample, so dropping samples would leave none.
  complete_features <- colSums(is.na(counts)) == 0
  counts_clean <- counts[, complete_features, drop = FALSE]

  if (nrow(counts_clean) < 3 || ncol(counts_clean) < 3) {
    return(NULL)
  }

  # Run PCA
  pca_result <- tryCatch(
    {
      stats::prcomp(counts_clean, center = center, scale. = scale)
    },
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(pca_result)) {
    return(NULL)
  }

  # Calculate variance explained
  variance <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

  list(
    scores = tibble::as_tibble(pca_result$x, rownames = "Name"),
    loadings = tibble::as_tibble(pca_result$rotation, rownames = "gene"),
    variance = variance,
    method = "PCA",
    component_names = paste0("PC", seq_along(variance)),
    removed_features = sum(!complete_features)
  )
}

#' calculate_pcoa
#' Run PCoA (via cmdscale) or NMDS (via vegan::metaMDS) on a distance matrix.
#' @param counts Matrix with samples as rows, features as columns
#' @param method Distance method passed to vegan::vegdist
#' @param nmds Logical, run NMDS instead of PCoA
#' @return List with scores, loadings (NULL), variance / stress, method, component_names
#' @noRd
calculate_pcoa <- function(counts, method = "euclidean", nmds = FALSE) {
  # Remove features (columns) with missing values, as calculate_pca() does.
  complete_features <- colSums(is.na(counts)) == 0
  counts_clean <- counts[, complete_features, drop = FALSE]

  if (nrow(counts_clean) < 3 || ncol(counts_clean) < 3) {
    return(NULL)
  }

  # Calculate distance matrix
  dist_mat <- tryCatch(
    {
      vegan::vegdist(counts_clean, method = method)
    },
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(dist_mat)) {
    return(NULL)
  }

  if (nmds) {
    # NMDS
    nmds_result <- tryCatch(
      {
        vegan::metaMDS(dist_mat, k = 2, trymax = 100, trace = 0)
      },
      error = function(e) {
        return(NULL)
      }
    )

    if (is.null(nmds_result)) {
      return(NULL)
    }

    scores_mat <- vegan::scores(nmds_result)

    # Name the columns before converting to tibble
    component_names <- paste0("NMDS", seq_len(ncol(scores_mat)))
    colnames(scores_mat) <- component_names

    list(
      scores = tibble::as_tibble(scores_mat, rownames = "Name"),
      loadings = NULL,
      variance = NULL,
      stress = nmds_result$stress,
      method = "NMDS",
      component_names = component_names,
      removed_features = sum(!complete_features)
    )
  } else {
    # Classical PCoA
    pcoa_result <- tryCatch(
      {
        stats::cmdscale(dist_mat, k = min(10, nrow(counts_clean) - 1), eig = TRUE)
      },
      error = function(e) {
        return(NULL)
      }
    )

    if (is.null(pcoa_result)) {
      return(NULL)
    }

    # Name the columns before converting to tibble
    n_dims <- ncol(pcoa_result$points)
    component_names <- paste0("PCo", seq_len(n_dims))
    colnames(pcoa_result$points) <- component_names

    # Calculate variance explained - only for the dimensions we kept
    # cmdscale returns all eigenvalues, but we only use the first n_dims
    variance_explained <- (pcoa_result$eig[seq_len(n_dims)] / sum(abs(pcoa_result$eig))) * 100

    list(
      scores = tibble::as_tibble(pcoa_result$points, rownames = "Name"),
      loadings = NULL,
      variance = variance_explained,
      stress = NULL,
      method = "PCoA",
      component_names = component_names,
      removed_features = sum(!complete_features)
    )
  }
}

# nolint start: cyclocomp_linter. Input checks before glmpca.
#' calculate_glmpca
#' Run GLM-PCA on raw counts using the glmpca package.
#' @param counts Matrix with features as rows, samples as columns (RAW counts)
#' @param L Number of latent dimensions
#' @param family One of "poi", "nb", "binom"
#' @param X Optional covariate matrix
#' @return List with scores, loadings (NULL), variance (NULL), method, component_names
#' @noRd
calculate_glmpca <- function(counts, L = 10, family = "nb", X = NULL) {
  if (nrow(counts) < 3 || ncol(counts) < 3) {
    return(NULL)
  }

  # Validate family parameter
  if (is.null(family) || !is.character(family) || length(family) != 1) {
    warning("GLM-PCA: Invalid family parameter (must be a single character string)")
    return(NULL)
  }

  valid_families <- c("poi", "nb", "binom")
  if (!family %in% valid_families) {
    warning("GLM-PCA: family must be one of: poi, nb, binom. Got: ", family)
    return(NULL)
  }

  # Remove genes with all zeros (required by glmpca)
  row_sums <- rowSums(counts)
  non_zero_rows <- row_sums > 0
  if (sum(non_zero_rows) < 3) {
    warning("GLM-PCA: Too few non-zero genes remaining")
    return(NULL)
  }
  counts <- counts[non_zero_rows, ]

  if (!requireNamespace("glmpca", quietly = TRUE)) {
    warning("GLM-PCA needs the glmpca package, which is not installed")
    return(NULL)
  }

  # GLM-PCA needs raw/normalized counts (not log-transformed)
  glmpca_result <- tryCatch(
    {
      glmpca::glmpca(counts, L = L, fam = family, X = X)
    },
    error = function(e) {
      warning("GLM-PCA failed: ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(glmpca_result)) {
    return(NULL)
  }

  # Name the columns before converting to tibble
  component_names <- paste0("dim", seq_len(L))
  colnames(glmpca_result$factors) <- component_names

  list(
    scores = tibble::as_tibble(glmpca_result$factors, rownames = "Name"),
    loadings = NULL, # GLM-PCA loadings are on a different scale, biplot not supported
    variance = NULL, # GLM-PCA doesn't have simple variance explained like PCA
    method = "GLM-PCA",
    component_names = component_names,
    removed_samples = 0 # GLM-PCA handles missing values internally
  )
}
# nolint end

#' prepare_biplot_data
#' Pick which features to overlay as arrows on a biplot based on selection method.
#' @param ordination_result Result from calculate_* function
#' @param top_n Number of top features
#' @param custom_genes Manual selection
#' @param pc_x PC index or name for X
#' @param pc_y PC index or name for Y
#' @param selection_method One of "overall", "per_pc", "manual"
#' @return List with sample_data, gene_data, pc_x_name, pc_y_name
#' @noRd
prepare_biplot_data <- function(
  ordination_result,
  top_n = 10,
  custom_genes = NULL,
  pc_x = 1,
  pc_y = 2,
  selection_method = "overall"
) {
  if (is.null(ordination_result$loadings)) {
    # Return NULL gracefully instead of stopping - biplot not available for distance-based methods
    return(NULL)
  }

  loadings <- ordination_result$loadings

  # Get component names (pc_x and pc_y can be either indices or component name strings)
  pc_x_name <- if (is.numeric(pc_x)) ordination_result$component_names[pc_x] else pc_x
  pc_y_name <- if (is.numeric(pc_y)) ordination_result$component_names[pc_y] else pc_y

  # Calculate loading magnitude
  loadings_subset <- loadings |>
    dplyr::mutate(magnitude = sqrt(.data[[pc_x_name]]^2 + .data[[pc_y_name]]^2))

  # Select genes based on selection method
  if (selection_method == "manual") {
    # Manual selection: use custom_genes
    if (!is.null(custom_genes) && length(custom_genes) > 0) {
      genes_to_plot <- loadings_subset |>
        dplyr::filter(.data$gene %in% custom_genes)
    } else {
      genes_to_plot <- loadings_subset[0, ] # Empty tibble if no genes specified
    }
  } else if (selection_method == "per_pc") {
    # Per-PC selection: top N from each PC separately
    top_x <- loadings_subset |>
      dplyr::arrange(dplyr::desc(abs(.data[[pc_x_name]]))) |>
      utils::head(top_n) |>
      dplyr::mutate(pc_source = pc_x_name)

    top_y <- loadings_subset |>
      dplyr::arrange(dplyr::desc(abs(.data[[pc_y_name]]))) |>
      utils::head(top_n) |>
      dplyr::mutate(pc_source = pc_y_name)

    # Combine and deduplicate - keep first occurrence (prioritize PC X if gene is in both)
    genes_to_plot <- dplyr::bind_rows(top_x, top_y) |>
      dplyr::distinct(.data$gene, .keep_all = TRUE)
  } else {
    # Overall selection: top N by combined magnitude (default behavior)
    genes_to_plot <- loadings_subset |>
      dplyr::arrange(dplyr::desc(.data$magnitude)) |>
      utils::head(top_n)
  }

  list(
    sample_data = ordination_result$scores,
    gene_data = genes_to_plot,
    pc_x_name = pc_x_name,
    pc_y_name = pc_y_name
  )
}

#' step_group_subset
#' Cascade step: subset counts and dataset to the specified groups.
#' @noRd
step_group_subset <- function(counts, dataset, factor_col, groups) {
  dataset <- dataset[dataset[[factor_col]] %in% groups, ]
  counts <- counts[, rownames(dataset), drop = FALSE]
  list(counts = counts, dataset = dataset)
}

#' step_peptide_filter
#' Cascade step: keep features with >= min_peptides entries (proteomics only).
#' @noRd
step_peptide_filter <- function(counts, dataset, de_df, min_peptides) {
  if (!is.null(de_df) && !is.null(min_peptides) && "nrPeptides" %in% colnames(de_df)) {
    keep <- de_df$feature_id[de_df$nrPeptides >= min_peptides]
    counts <- counts[rownames(counts) %in% keep, , drop = FALSE]
  }
  list(counts = counts, dataset = dataset)
}

#' step_top_n_sd
#' Cascade step: keep top-N features by row-wise SD.
#' @noRd
step_top_n_sd <- function(counts, dataset, top_n) {
  top_n <- min(top_n, nrow(counts))
  if (nrow(counts) > top_n) {
    sds <- apply(counts, 1, stats::sd, na.rm = FALSE)
    top_features <- names(utils::head(sort(sds, decreasing = TRUE), n = top_n))
    counts <- counts[top_features, , drop = FALSE]
  }
  list(counts = counts, dataset = dataset)
}

#' ordination_axis_labels
#' Axis titles of an ordination plot: the explained variance per axis (PCA,
#' PCoA), or the stress on the x axis (NMDS), else the axis names.
#' @param ordination_result List with `variance` (percent per axis) or `stress`.
#' @param pc_x,pc_y Axis names, e.g. "PC1"; their number indexes `variance`.
#' @return `list(x, y)`.
#' @noRd
ordination_axis_labels <- function(ordination_result, pc_x, pc_y) {
  variance <- ordination_result$variance
  if (!is.null(variance)) {
    pct <- function(pc) paste0(pc, ": ", round(variance[as.numeric(gsub("\\D", "", pc))]), "% variance")
    return(list(x = pct(pc_x), y = pct(pc_y)))
  }
  if (!is.null(ordination_result$stress)) {
    return(list(x = paste0(pc_x, " (stress: ", round(ordination_result$stress, 3), ")"), y = pc_y))
  }
  list(x = pc_x, y = pc_y)
}

#' ordination_theme
#' Theme settings of the ordination plots, added after `theme_bw()`.
#' @param bold Whether axis texts and titles are bold.
#' @param show_lines Whether grid lines are drawn.
#' @param show_axes Whether axes, ticks, axis texts and the panel border are drawn.
#' @param rich Whether legend labels are rich text (rich_label_theme()).
#' @return A list of ggplot2 theme components.
#' @noRd
ordination_theme <- function(bold, show_lines, show_axes, rich = FALSE) {
  face <- if (isTRUE(bold)) "bold" else "plain"
  blank <- ggplot2::element_blank()
  c(
    list(ggplot2::theme(
      axis.text = ggplot2::element_text(face = face),
      axis.title = ggplot2::element_text(face = face)
    )),
    if (rich) list(rich_label_theme(legend = TRUE)),
    if (!show_lines) {
      list(ggplot2::theme(
        panel.grid.major.x = blank,
        panel.grid.minor.x = blank,
        panel.grid.major.y = blank,
        panel.grid.minor.y = blank
      ))
    },
    if (!show_axes) {
      list(ggplot2::theme(
        axis.line = blank,
        axis.text.x = blank,
        axis.text.y = blank,
        axis.ticks = blank,
        axis.title.x = blank,
        axis.title.y = blank,
        panel.border = blank
      ))
    }
  )
}

#' plot_ordination_samples
#' Render the sample-only ordination scatter plot (no gene arrows).
#' @noRd
plot_ordination_samples <- function(ordination_result, pc_scores, input, coloursPCA, levelLabels = NULL) {
  fillLabels <- if (is.null(levelLabels$fill)) ggplot2::waiver() else levelLabels$fill
  shapeLabels <- if (is.null(levelLabels$shape)) ggplot2::waiver() else levelLabels$shape
  pc_x <- input$pcaX
  pc_y <- input$pcaY

  # Base plot
  p <- ggplot2::ggplot(pc_scores, ggplot2::aes(x = .data[[pc_x]], y = .data[[pc_y]], fill = .data[[input$pcaFactor1]]))

  if (input$pcaAddEllipses) {
    p <- p +
      ggplot2::stat_ellipse(
        ggplot2::aes(colour = .data[[input$pcaFactor1]], fill = .data[[input$pcaFactor1]]),
        geom = "polygon",
        alpha = input$pcaEllipsesAlpha,
        show.legend = FALSE
      )
    p <- p + ggplot2::scale_colour_manual(breaks = names(coloursPCA), values = as.character(coloursPCA))
  }

  if (input$pcaFactor2 == "None") {
    p <- p +
      ggplot2::geom_point(size = input$dotSizePCA, alpha = input$alphaPCA, shape = 21, stroke = input$dotBorderPCA)
  } else {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(shape = .data[[input$pcaFactor2]]),
        size = input$dotSizePCA,
        alpha = input$alphaPCA,
        stroke = input$dotBorderPCA
      )
    p <- p + ggplot2::scale_shape_manual(values = c(21, 22, 23, 24, 25, 8, 3, 4), labels = shapeLabels)
  }

  if (input$pcaShowNames) {
    p <- p +
      ggrepel::geom_label_repel(
        ggplot2::aes(label = .data$Name, fill = .data[[input$pcaFactor1]]),
        force = 2,
        nudge_y = as.numeric(input$geneLabelNudgePCAY) / 10,
        nudge_x = as.numeric(input$geneLabelNudgePCAX) / 10,
        size = (input$textSizePCA / 3),
        max.overlaps = input$pcaLabelMaxOverlap,
        fontface = "bold",
        show.legend = FALSE,
        na.rm = TRUE
      )
  }

  axis_labels <- ordination_axis_labels(ordination_result, pc_x, pc_y)
  p <- p + ggplot2::labs(x = axis_labels$x, y = axis_labels$y, fill = input$pcaFactor1)

  if (input$pcaAxesProp && !is.null(ordination_result$variance)) {
    pc_x_num <- as.numeric(gsub("\\D", "", pc_x))
    pc_y_num <- as.numeric(gsub("\\D", "", pc_y))
    p <- p + ggplot2::coord_fixed(ratio = ordination_result$variance[pc_y_num] / ordination_result$variance[pc_x_num])
  }

  p <- p +
    ggplot2::scale_fill_manual(breaks = names(coloursPCA), values = as.character(coloursPCA), labels = fillLabels)
  p +
    ggplot2::theme_bw(base_size = input$textSizePCA) +
    ordination_theme(input$boldPCA, input$showLinesPCA, input$showAxesPCA, isTRUE(levelLabels$rich))
}

#' plot_ordination_biplot
#' Render the ordination biplot (samples + feature arrows).
#' @noRd
plot_ordination_biplot <- function(
  ordination_result,
  pc_scores,
  input,
  coloursPCA,
  biplot_data,
  arrow_scale = 1,
  arrow_colour = "red",
  selection_method = "overall",
  levelLabels = NULL
) {
  fillLabels <- if (is.null(levelLabels$fill)) ggplot2::waiver() else levelLabels$fill
  pc_x <- biplot_data$pc_x_name
  pc_y <- biplot_data$pc_y_name

  # Scale arrows
  sample_range <- max(abs(c(biplot_data$sample_data[[pc_x]], biplot_data$sample_data[[pc_y]])), na.rm = TRUE)
  gene_range <- max(abs(c(biplot_data$gene_data[[pc_x]], biplot_data$gene_data[[pc_y]])), na.rm = TRUE)
  scale_factor <- (sample_range / gene_range) * arrow_scale

  biplot_data$gene_data <- biplot_data$gene_data |>
    dplyr::mutate(pc_x_scaled = .data[[pc_x]] * scale_factor, pc_y_scaled = .data[[pc_y]] * scale_factor)

  # Add arrow color column for per_pc mode
  if (selection_method == "per_pc" && "pc_source" %in% names(biplot_data$gene_data)) {
    biplot_data$gene_data <- biplot_data$gene_data |>
      dplyr::mutate(arrow_color = ifelse(.data$pc_source == pc_x, arrow_colour$x, arrow_colour$y))
  } else {
    # Single color for all arrows
    biplot_data$gene_data$arrow_color <- if (is.list(arrow_colour)) arrow_colour$x else arrow_colour
  }

  p <- ggplot2::ggplot()
  p <- p +
    ggplot2::geom_point(
      data = pc_scores,
      ggplot2::aes(x = .data[[pc_x]], y = .data[[pc_y]], fill = .data[[input$pcaFactor1]]),
      size = input$dotSizePCA,
      alpha = input$alphaPCA,
      shape = 21,
      stroke = input$dotBorderPCA
    )
  p <- p +
    ggplot2::geom_segment(
      data = biplot_data$gene_data,
      ggplot2::aes(x = 0, y = 0, xend = .data$pc_x_scaled, yend = .data$pc_y_scaled, colour = .data$arrow_color),
      arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "closed"),
      linewidth = 0.5,
      alpha = 0.7
    )
  # Arrow labels: optional display label (e.g. gene symbol), else the feature_id
  if (!"label" %in% names(biplot_data$gene_data)) {
    biplot_data$gene_data$label <- biplot_data$gene_data$gene
  }
  p <- p +
    ggrepel::geom_text_repel(
      data = biplot_data$gene_data,
      ggplot2::aes(x = .data$pc_x_scaled, y = .data$pc_y_scaled, label = .data$label, colour = .data$arrow_color),
      size = input$textSizePCA / 4,
      fontface = "bold",
      max.overlaps = 20
    )
  p <- p + ggplot2::scale_colour_identity() # Use the actual color values

  axis_labels <- ordination_axis_labels(ordination_result, pc_x, pc_y)
  p <- p + ggplot2::labs(x = axis_labels$x, y = axis_labels$y, fill = input$pcaFactor1)
  p <- p +
    ggplot2::scale_fill_manual(breaks = names(coloursPCA), values = as.character(coloursPCA), labels = fillLabels)
  p +
    ggplot2::theme_bw(base_size = input$textSizePCA) +
    ordination_theme(input$boldPCA, input$showLinesPCA, input$showAxesPCA, isTRUE(levelLabels$rich))
}
