# ============================================================================
# make_rep_heatmap / make_rep_dendrogram — smoke tests
#
# Visualisation-heavy: they read the Replicability settings from `input` (a
# plain list here) and call shiny::req() / col2hex() / dendextend internals.
# The tests assert they run to completion. Visual output is verified by eye /
# shinytest2 end-to-end, not here.
# ============================================================================

# replicability_plots.R calls functions from these packages without namespace
# qualifiers (it was written for a Shiny app context where all are attached).
# Attach them explicitly here since devtools::load_all() doesn't attach Imports.
skip_if_not_installed("dendextend")
library(grid) # gpar()
library(dendextend) # labels_colors<-
library(ComplexHeatmap) # HeatmapAnnotation, Heatmap, draw
library(circlize) # colorRamp2
library(RColorBrewer) # brewer.pal

make_rep_input <- function(overrides = list()) {
  base <- list(
    replicabilityTextSize = 12,
    replicabilityFactors = NULL, # start simple: no factor annotations
    replicabilityPalette = "None",
    replicabilityPaletteRev = FALSE,
    replicabilityAtLow = -1,
    replicabilityAtMid = 0,
    replicabilityAtHigh = 1,
    replicabilityColorLow = "#1F78B4",
    replicabilityColorMid = "white",
    replicabilityColorHigh = "#E31A1C",
    replicabilityBorderColor1 = "black",
    replicabilityBorderColor2 = "black",
    replicabilityBorderSize1 = 2,
    replicabilityBorderSize2 = 2,
    replicabilityShowNames = TRUE,
    replicabilityDistMethod = "euclidean",
    replicabilityClustMethod = "complete"
  )
  modifyList(base, overrides)
}

make_rep_data <- function(n = 5) {
  set.seed(1)
  m <- matrix(rnorm(n * n), nrow = n)
  # Force symmetry (correlation-like)
  corr <- (m + t(m)) / 2
  diag(corr) <- 1
  rownames(corr) <- colnames(corr) <- paste0("s", seq_len(n))

  list(
    corr = corr,
    hclust = hclust(as.dist(1 - abs(corr))),
    dataset = data.frame(
      Name = paste0("s", seq_len(n)),
      cond = rep(c("A", "B"), length.out = n),
      row.names = paste0("s", seq_len(n))
    )
  )
}


# --- make_rep_heatmap -------------------------------------------------------
test_that("make_rep_heatmap runs without error when replicabilityFactors is NULL", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input()
  rep_fn <- make_rep_heatmap

  # draw() writes to a graphics device — capture it
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(rep_fn(data, input, title_suffix = "Smoke"), NA) # NA = no error
})

test_that("make_rep_heatmap runs with replicabilityFactors and GroupColour inputs", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input(list(
    replicabilityFactors = "cond",
    GroupColourcondA = "red",
    GroupColourcondB = "blue"
  ))
  rep_fn <- make_rep_heatmap

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(rep_fn(data, input, title_suffix = "Smoke-with-factor"), NA)
})

test_that("make_rep_heatmap runs with Colours-tab level labels", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input(list(
    replicabilityFactors = "cond",
    GroupColourcondA = "red",
    GroupColourcondB = "blue",
    # No Greek: this helper draws plain labels straight onto the base pdf() device below
    GroupLabelcondA = "<i>Il4ra</i><sup>fl/fl</sup>"
  ))
  rep_fn <- make_rep_heatmap

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(rep_fn(data, input, title_suffix = "Smoke-with-labels"), NA)
})

test_that("make_rep_heatmap runs with a BrewerPal-derived palette", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input(list(replicabilityPalette = "RdBu"))
  rep_fn <- make_rep_heatmap

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(rep_fn(data, input, title_suffix = "Smoke-with-palette"), NA)
})


# --- make_rep_dendrogram ----------------------------------------------------
test_that("make_rep_dendrogram runs without error when replicabilityFactors is NULL", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input()
  rep_fn <- make_rep_dendrogram

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(rep_fn(data, input, title_suffix = "Smoke"), NA)
})

test_that("make_rep_dendrogram runs with replicabilityFactors and GroupColour inputs", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input(list(
    replicabilityFactors = "cond",
    GroupColourcondA = "red",
    GroupColourcondB = "blue"
  ))
  rep_fn <- make_rep_dendrogram

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(rep_fn(data, input, title_suffix = "Smoke-with-factor"), NA)
})


# --- group_colours ----------------------------------------------------------
test_that("group_colours takes each level's colour from the Colours tab, grey without a picker", {
  input <- list(GroupColourcondA = "red", GroupColourcondB = "#0000FF")
  expect_identical(
    group_colours(input, "cond", c("A", "B", "C")),
    c(A = "#FF0000FF", B = "#0000FFFF", C = "#808080FF")
  )
})
