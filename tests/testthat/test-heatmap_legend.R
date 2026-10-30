library(ComplexHeatmap)

.legend_ht <- function() {
  set.seed(1)
  m <- matrix(rnorm(100), 10)
  Heatmap(m, name = "z", top_annotation = HeatmapAnnotation(grp = rep(c("a", "b"), 5)))
}

.legend_children <- function(g) {
  Filter(function(ch) grepl("(heatmap|annotation)_legend", as.character(ch$vp %||% "")), g$children)
}

test_that(".shift_heatmap_legends() wraps only the legend grobs, adding the offset viewport", {
  g <- grid::grid.grabExpr(draw(.legend_ht(), merge_legend = TRUE), width = 7, height = 7)
  n_before <- length(g$children)
  expect_gt(length(.legend_children(g)), 0)

  shifted <- .shift_heatmap_legends(g, dx = 10, dy = -5)
  expect_length(shifted$children, n_before)
  moved <- Filter(function(ch) inherits(ch$vp, "vpStack"), shifted$children)
  expect_length(moved, length(.legend_children(g)))
})

test_that("draw_heatmap() with zero / NULL offsets is a plain draw()", {
  f <- withr::local_tempfile(fileext = ".png")
  grDevices::png(f, 500, 500)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s4_class(draw_heatmap(.legend_ht(), legend_dx = NULL, legend_dy = 0), "HeatmapList")
})

test_that("draw_heatmap() with offsets draws without error on right and bottom legends", {
  f <- withr::local_tempfile(fileext = ".png")
  grDevices::png(f, 500, 500)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(draw_heatmap(.legend_ht(), merge_legend = TRUE, legend_dx = 20, legend_dy = 10))
  expect_no_error(draw_heatmap(
    .legend_ht(),
    merge_legend = TRUE,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    padding = grid::unit(c(2, 2, 2, 2), "cm"),
    legend_dx = -20,
    legend_dy = -15
  ))
})

test_that("draw_heatmap() copes with regex characters in heatmap / annotation names", {
  # Regression: "transformedData [Centred]" made redrawing the grabbed heatmap
  # fail with "Viewport '..._heatmap_body_1_1' was not found".
  set.seed(1)
  m <- matrix(rnorm(100), 10)
  ha <- HeatmapAnnotation("Condition [Factor]" = rep(c("a", "b"), 5))
  f <- withr::local_tempfile(fileext = ".png")
  grDevices::png(f, 500, 500)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (nm in c("transformedData \n[Centred]", "p (a|b)*")) {
    ht <- Heatmap(m, name = nm, top_annotation = ha, row_split = rep(1:2, 5))
    expect_no_error(draw_heatmap(ht, merge_legend = TRUE, legend_dx = 10, legend_dy = 5))
  }
})
