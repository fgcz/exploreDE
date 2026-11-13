# Replicability tab plots (server-replicability.R): correlation heatmap and dendrogram.

#' group_colours
#' The Colours-tab colours of a factor's levels (`GroupColour<factor><level>`
#' inputs) as "#RRGGBBFF", named by level; grey for a level without a picker.
#' @param input Shiny `input`, or a list with the same names.
#' @param factor Factor name.
#' @param levels Its levels.
#' @return Named character vector.
#' @noRd
group_colours <- function(input, factor, levels) {
  vapply(
    levels,
    function(l) {
      col <- input[[paste0("GroupColour", factor, l)]]
      if (is.null(col)) "#808080FF" else paste0(col2hex(col), "FF")
    },
    character(1)
  )
}

#' make_rep_heatmap
#' Render the replicability correlation heatmap.
#' @param data List containing `corr` (correlation matrix) and `dataset` (sample metadata)
#' @param input Shiny `input` (or a list with the same names): the Replicability
#'   settings and the Colours-tab colours and labels.
#' @param title_suffix String appended to heatmap title
#' @noRd
make_rep_heatmap <- function(data, input, title_suffix = "") {
  shiny::req(data)
  size <- input$replicabilityTextSize
  names_gp <- grid::gpar(fontsize = size)
  anno <- list(top = NULL, left = NULL)
  if (!is.null(input$replicabilityFactors)) {
    anno <- rep_annotations(data, input)
    names_gp <- grid::gpar(fontsize = size, col = anno$sample_colours)
  }
  hm <- ComplexHeatmap::Heatmap(
    data$corr,
    name = "Correlation",
    col = rep_heatmap_colours(input),
    top_annotation = anno$top,
    left_annotation = anno$left,
    show_row_names = input$replicabilityShowNames,
    show_column_names = input$replicabilityShowNames,
    row_names_gp = names_gp,
    column_names_gp = names_gp,
    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = size, fontface = "bold")),
    rect_gp = grid::gpar(col = input$replicabilityBorderColor2, lwd = input$replicabilityBorderSize2 / 10),
    column_title = paste0("Correlation Heatmap: ", title_suffix),
    column_title_gp = grid::gpar(fontsize = size + 2, fontface = "bold")
  )
  draw_heatmap(
    hm,
    merge_legend = TRUE,
    padding = grid::unit(c(10, 10, 10, 10), "mm"),
    legend_dx = input$replicabilityLegendDX,
    legend_dy = input$replicabilityLegendDY
  )
}

# Annotation bars of the selected factors (top and left) and the sample label
# colours, which follow the first selected factor.
rep_annotations <- function(data, input) {
  factors <- input$replicabilityFactors
  dataset <- data$dataset[match(rownames(data$corr), data$dataset$Name), ]
  colour_list <- stats::setNames(
    lapply(factors, function(f) group_colours(input, f, levels(as.factor(dataset[[f]])))),
    factors
  )
  anno_df <- dataset[, factors, drop = FALSE]
  list(
    top = rep_top_annotation(anno_df, colour_list, input),
    left = ComplexHeatmap::rowAnnotation(
      df = anno_df,
      col = colour_list,
      show_legend = FALSE,
      annotation_name_gp = grid::gpar(fontsize = input$replicabilityTextSize, fontface = "bold"),
      gp = grid::gpar(col = input$replicabilityBorderColor1, lwd = input$replicabilityBorderSize1 / 10)
    ),
    sample_colours = stats::setNames(colour_list[[1]][as.character(dataset[[factors[1]]])], dataset$Name)
  )
}

# The top annotation. Its legend shows the Colours-tab display labels: rich
# ones via legends built at draw time by draw_heatmap() (see R/level_labels.R).
rep_top_annotation <- function(anno_df, colour_list, input) {
  factors <- names(colour_list)
  size <- input$replicabilityTextSize
  rich <- any(vapply(
    factors,
    function(f) is_rich_labels(level_labels(input, f, names(colour_list[[f]]))),
    logical(1)
  ))
  legend_param <- stats::setNames(
    lapply(factors, function(f) {
      levs <- names(colour_list[[f]])
      list(at = levs, labels = unname(level_labels(input, f, levs, mode = "plain")))
    }),
    factors
  )
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = anno_df,
    col = colour_list,
    show_legend = !rich,
    annotation_legend_param = legend_param,
    annotation_name_gp = grid::gpar(fontsize = size, fontface = "bold"),
    gp = grid::gpar(col = input$replicabilityBorderColor1, lwd = input$replicabilityBorderSize1 / 10)
  )
  if (rich) {
    attr(ha, "rich_legends") <- lapply(factors, function(f) {
      heatmap_legend_spec(
        title = f,
        labels = level_labels(input, f, names(colour_list[[f]]), force_rich = TRUE),
        fill = colour_list[[f]],
        title_gp = grid::gpar(fontsize = size, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = size)
      )
    })
  }
  ha
}

# Correlation colour scale: a Brewer palette, else the three picked colours.
rep_heatmap_colours <- function(input) {
  at <- c(input$replicabilityAtLow, input$replicabilityAtMid, input$replicabilityAtHigh)
  if (is.null(input$replicabilityPalette) || input$replicabilityPalette == "None") {
    return(circlize::colorRamp2(
      at,
      c(input$replicabilityColorLow, input$replicabilityColorMid, input$replicabilityColorHigh)
    ))
  }
  selected <- RColorBrewer::brewer.pal(9, input$replicabilityPalette)[c(1, 5, 9)]
  circlize::colorRamp2(at, if (input$replicabilityPaletteRev) selected else rev(selected))
}

#' make_rep_dendrogram
#' Render the replicability dendrogram with coloured labels / coloured bars.
#' @param data List containing `hclust` (hclust object) and `dataset` (sample metadata)
#' @param input Shiny `input` (or a list with the same names), as for
#'   make_rep_heatmap().
#' @param title_suffix String appended to title
#' @noRd
make_rep_dendrogram <- function(data, input, title_suffix = "") {
  shiny::req(data)
  dend <- stats::as.dendrogram(data$hclust)

  if (!is.null(input$replicabilityFactors)) {
    dataset <- data$dataset
    factor_col <- input$replicabilityFactors[1]
    colors <- group_colours(input, factor_col, levels(as.factor(dataset[[factor_col]])))

    sample_colors <- colors[as.character(dataset[[factor_col]])]
    names(sample_colors) <- dataset$Name
    dendextend::labels_colors(dend) <- sample_colors[labels(dend)]

    colors_mat <- as.data.frame(dataset[, input$replicabilityFactors, drop = FALSE])
    for (f in input$replicabilityFactors) {
      colors_mat[[f]] <- group_colours(input, f, as.character(colors_mat[[f]]))
    }
    colors_mat <- colors_mat[match(labels(dend), dataset$Name), , drop = FALSE]
  }

  old_par <- graphics::par(mar = c(20, 8, 6, 4))
  on.exit(graphics::par(old_par), add = TRUE)
  plot(
    dend,
    main = paste0(
      "Cluster Dendrogram (",
      input$replicabilityDistMethod,
      "/",
      input$replicabilityClustMethod,
      "): ",
      title_suffix
    ),
    xlab = "",
    sub = "",
    cex = input$replicabilityTextSize / 10
  )

  if (!is.null(input$replicabilityFactors)) {
    dendextend::colored_bars(colors = colors_mat, dend = dend, sort_by_labels_order = FALSE)
  }
}
