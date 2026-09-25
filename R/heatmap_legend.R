# =============================================================================
# Heatmap legend offset: nudge a ComplexHeatmap's legends around by hand.
#
# ComplexHeatmap has no argument for moving the legends, so draw_heatmap()
# draws the heatmap off-screen (grid.grabExpr), shifts the legend grobs by
# (legend_dx, legend_dy) mm, and draws the result. With both offsets at 0 it is
# exactly ComplexHeatmap::draw(). Every tab with a ComplexHeatmap exposes the
# offsets through legend_offset_ui().
# =============================================================================

#' Legend offset sliders for a ComplexHeatmap tab
#'
#' Two sliders, `<prefix>LegendDX` and `<prefix>LegendDY` (mm), both 0 by
#' default. Pass their values to [draw_heatmap()].
#' @param prefix Input-id prefix for the tab (e.g. `"heatmap"`).
#' @return A shiny tag list.
#' @keywords internal
#' @export
legend_offset_ui <- function(prefix) {
  shiny::splitLayout(
    shiny::sliderInput(paste0(prefix, "LegendDX"), "Legend: move right (mm)",
                       min = -100, max = 100, value = 0, step = 1, ticks = FALSE, width = "85%"),
    shiny::sliderInput(paste0(prefix, "LegendDY"), "Legend: move up (mm)",
                       min = -100, max = 100, value = 0, step = 1, ticks = FALSE, width = "85%")
  )
}

# Shift every legend grob in a grabbed heatmap gTree by (dx, dy) mm. Legends
# are recognised by their viewport path (".../heatmap_legend/..." or
# ".../annotation_legend/..."), which holds whatever the legend layout
# (single legend, packed legends, right or bottom side).
.shift_heatmap_legends <- function(g, dx, dy) {
  move <- grid::viewport(x = grid::unit(0.5, "npc") + grid::unit(dx, "mm"),
                         y = grid::unit(0.5, "npc") + grid::unit(dy, "mm"))
  for (i in seq_along(g$children)) {
    ch <- g$children[[i]]
    if (inherits(ch$vp, "vpPath") &&
        grepl("(heatmap|annotation)_legend", as.character(ch$vp))) {
      g$children[[i]] <- grid::gTree(children = grid::gList(grid::editGrob(ch, vp = NULL)),
                                     vp = grid::vpStack(ch$vp, move), name = ch$name)
    }
  }
  g
}

# Redrawing a grabbed gTree navigates to viewports by name, and grid matches
# those names as regular expressions: a heatmap called "transformedData
# [Centred]" names its body viewport "transformedData [Centred]_heatmap_body_1_1",
# which then fails to match itself ("Viewport ... was not found"). Replace
# regex metacharacters in every viewport name, consistently in the viewport
# tree and in the grobs' viewport paths. Only viewport names change; visible
# titles and labels are separate text grobs and are untouched.
.vp_safe <- function(x) gsub("[][(){}*+?^$|\\\\]", "_", x)

.sanitise_vp <- function(vp) {
  if (is.null(vp)) return(vp)
  if (inherits(vp, "vpPath")) {
    return(do.call(grid::vpPath, as.list(.vp_safe(strsplit(as.character(vp), "::", fixed = TRUE)[[1]]))))
  }
  if (inherits(vp, "vpTree")) {
    vp$parent   <- .sanitise_vp(vp$parent)
    vp$children <- .sanitise_vp(vp$children)
    return(vp)
  }
  if (inherits(vp, c("vpList", "vpStack"))) {
    cls <- class(vp)
    vp  <- lapply(vp, .sanitise_vp)
    class(vp) <- cls
    return(vp)
  }
  if (inherits(vp, "viewport")) vp$name <- .vp_safe(vp$name)
  vp
}

.sanitise_grob_vps <- function(g) {
  g$vp <- .sanitise_vp(g$vp)
  if (inherits(g, "gTree")) {
    g$childrenvp <- .sanitise_vp(g$childrenvp)
    for (i in seq_along(g$children)) g$children[[i]] <- .sanitise_grob_vps(g$children[[i]])
  }
  g
}

#' Draw a ComplexHeatmap, optionally with its legends moved
#'
#' Drop-in replacement for `ComplexHeatmap::draw(ht, ...)`. With non-zero
#' offsets the legends are moved right/up by `legend_dx`/`legend_dy` mm
#' (negative = left/down). Moving a right-hand legend further right, or a
#' bottom legend further down, adds the same amount of padding on that side so
#' the legend stays on the canvas.
#' @param ht A `Heatmap` / `HeatmapList`.
#' @param ... Passed to `ComplexHeatmap::draw()` (e.g. `merge_legend`, `padding`).
#' @param legend_dx,legend_dy Legend offsets in mm (`NULL` = 0).
#' @keywords internal
#' @export
draw_heatmap <- function(ht, ..., legend_dx = 0, legend_dy = 0) {
  dx <- as.numeric(legend_dx %||% 0)
  dy <- as.numeric(legend_dy %||% 0)
  args <- list(...)
  # Start the page ourselves: building the rich legends below measures text
  # on the device, which makes grid open a page -- so draw() must then NOT
  # start another, or PDFs get a blank first page.
  grid::grid.newpage()
  # Rich level-label legends/column labels (R/level_labels.R), built here on
  # the real device, before any off-screen grab below.
  prep <- prepare_rich_heatmap(ht, args)
  ht <- prep$ht
  args <- prep$args
  if (!dx && !dy) {
    return(invisible(do.call(ComplexHeatmap::draw, c(list(ht), args, list(newpage = FALSE)))))
  }

  # padding is c(bottom, left, top, right)
  side <- args$heatmap_legend_side %||% "right"
  pad  <- args$padding %||% grid::unit(c(2, 2, 2, 2), "mm")
  if (side == "right"  && dx > 0) pad <- pad + grid::unit(c(0, 0, 0, dx), "mm")
  if (side == "bottom" && dy < 0) pad <- pad + grid::unit(c(-dy, 0, 0, 0), "mm")
  args$padding <- pad

  sz <- grDevices::dev.size("in")
  g <- grid::grid.grabExpr(do.call(ComplexHeatmap::draw, c(list(ht), args)),
                           width = sz[1], height = sz[2])
  grid::grid.draw(.shift_heatmap_legends(.sanitise_grob_vps(g), dx, dy))
  invisible(NULL)
}
