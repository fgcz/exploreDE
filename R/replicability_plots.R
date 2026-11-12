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

  # Prepare annotations
  ha_top <- NULL
  ha_left <- NULL
  row_names_gp <- grid::gpar(fontsize = input$replicabilityTextSize)
  column_names_gp <- grid::gpar(fontsize = input$replicabilityTextSize)

  if (!is.null(input$replicabilityFactors)) {
    dataset <- data$dataset
    sample_names <- rownames(data$corr)
    dataset <- dataset[match(sample_names, dataset$Name), ]

    # Get colors
    colourList <- lapply(input$replicabilityFactors, function(f) {
      group_colours(input, f, levels(as.factor(dataset[[f]])))
    }) |>
      stats::setNames(input$replicabilityFactors)

    # Color labels based on the first selected factor
    factor_col <- input$replicabilityFactors[1]
    colors <- group_colours(input, factor_col, levels(as.factor(dataset[[factor_col]])))

    sample_colors <- colors[as.character(dataset[[factor_col]])]
    names(sample_colors) <- dataset$Name

    row_names_gp <- grid::gpar(fontsize = input$replicabilityTextSize, col = sample_colors)
    column_names_gp <- grid::gpar(fontsize = input$replicabilityTextSize, col = sample_colors)

    anno_df <- dataset[, input$replicabilityFactors, drop = FALSE]

    # Legend shows the Colours-tab display labels: rich ones via legends
    # built at draw time by draw_heatmap() (see R/level_labels.R)
    repRich <- any(vapply(
      input$replicabilityFactors,
      function(f) {
        is_rich_labels(level_labels(input, f, names(colourList[[f]])))
      },
      logical(1)
    ))
    annoLegendParam <- lapply(input$replicabilityFactors, function(f) {
      levs <- names(colourList[[f]])
      list(at = levs, labels = unname(level_labels(input, f, levs, mode = "plain")))
    }) |>
      stats::setNames(input$replicabilityFactors)

    ha_top <- ComplexHeatmap::HeatmapAnnotation(
      df = anno_df,
      col = colourList,
      show_legend = !repRich,
      annotation_legend_param = annoLegendParam,
      annotation_name_gp = grid::gpar(fontsize = input$replicabilityTextSize, fontface = "bold"),
      gp = grid::gpar(col = input$replicabilityBorderColor1, lwd = input$replicabilityBorderSize1 / 10)
    )
    if (repRich) {
      attr(ha_top, "rich_legends") <- lapply(input$replicabilityFactors, function(f) {
        heatmap_legend_spec(
          title = f,
          labels = level_labels(input, f, names(colourList[[f]]), force_rich = TRUE),
          fill = colourList[[f]],
          title_gp = grid::gpar(fontsize = input$replicabilityTextSize, fontface = "bold"),
          labels_gp = grid::gpar(fontsize = input$replicabilityTextSize)
        )
      })
    }

    ha_left <- ComplexHeatmap::rowAnnotation(
      df = anno_df,
      col = colourList,
      show_legend = FALSE,
      annotation_name_gp = grid::gpar(fontsize = input$replicabilityTextSize, fontface = "bold"),
      gp = grid::gpar(col = input$replicabilityBorderColor1, lwd = input$replicabilityBorderSize1 / 10)
    )
  }

  # Get color breaks
  at <- c(input$replicabilityAtLow, input$replicabilityAtMid, input$replicabilityAtHigh)

  # Determine color scheme - use palette if selected, otherwise use manual colors
  if (!is.null(input$replicabilityPalette) && input$replicabilityPalette != "None") {
    palette_colours <- RColorBrewer::brewer.pal(9, input$replicabilityPalette)
    selected_colours <- palette_colours[c(1, 5, 9)]
    if (input$replicabilityPaletteRev) {
      heatmap_colours <- circlize::colorRamp2(at, selected_colours)
    } else {
      heatmap_colours <- circlize::colorRamp2(at, rev(selected_colours))
    }
  } else {
    heatmap_colours <- circlize::colorRamp2(
      at,
      c(input$replicabilityColorLow, input$replicabilityColorMid, input$replicabilityColorHigh)
    )
  }

  hm <- ComplexHeatmap::Heatmap(
    data$corr,
    name = "Correlation",
    col = heatmap_colours,
    top_annotation = ha_top,
    left_annotation = ha_left,
    show_row_names = input$replicabilityShowNames,
    show_column_names = input$replicabilityShowNames,
    row_names_gp = row_names_gp,
    column_names_gp = column_names_gp,
    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = input$replicabilityTextSize, fontface = "bold")),
    rect_gp = grid::gpar(col = input$replicabilityBorderColor2, lwd = input$replicabilityBorderSize2 / 10),
    column_title = paste0("Correlation Heatmap: ", title_suffix),
    column_title_gp = grid::gpar(fontsize = input$replicabilityTextSize + 2, fontface = "bold")
  )

  draw_heatmap(
    hm,
    merge_legend = TRUE,
    padding = grid::unit(c(10, 10, 10, 10), "mm"),
    legend_dx = input$replicabilityLegendDX,
    legend_dy = input$replicabilityLegendDY
  )
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
