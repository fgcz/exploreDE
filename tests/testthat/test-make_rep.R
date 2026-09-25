# ============================================================================
# make_rep_heatmap / make_rep_dendrogram — smoke tests
#
# These are visualisation-heavy and tightly Shiny-coupled (reference
# `input$...` from enclosing scope and call shiny::req() / col2hex() /
# dendextend internals). We inject a complete mock `input` via with_mock_env()
# and assert the functions run to completion without erroring. Visual output
# is verified by eye / shinytest2 end-to-end, not here.
# ============================================================================

# helperFunctions.R calls functions from these packages without namespace
# qualifiers (it was written for a Shiny app context where all are attached).
# Attach them explicitly here since devtools::load_all() doesn't attach Imports.
library(grid)           # gpar()
library(dendextend)     # labels_colors<-
library(ComplexHeatmap) # HeatmapAnnotation, Heatmap, draw
library(circlize)       # colorRamp2
library(RColorBrewer)   # brewer.pal

make_rep_input <- function(overrides = list()) {
  base <- list(
    replicabilityTextSize      = 12,
    replicabilityFactors       = NULL,  # start simple: no factor annotations
    replicabilityPalette       = "None",
    replicabilityPaletteRev    = FALSE,
    replicabilityAtLow         = -1,
    replicabilityAtMid         = 0,
    replicabilityAtHigh        = 1,
    replicabilityColorLow      = "#1F78B4",
    replicabilityColorMid      = "white",
    replicabilityColorHigh     = "#E31A1C",
    replicabilityBorderColor1  = "black",
    replicabilityBorderColor2  = "black",
    replicabilityBorderSize1   = 2,
    replicabilityBorderSize2   = 2,
    replicabilityShowNames     = TRUE,
    replicabilityDistMethod    = "euclidean",
    replicabilityClustMethod   = "complete"
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
    corr    = corr,
    hclust  = hclust(as.dist(1 - abs(corr))),
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
  fn <- with_mock_env(make_rep_heatmap, input = input)

  # draw() writes to a graphics device — capture it
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(fn(data, title_suffix = "Smoke"), NA)  # NA = no error
})

test_that("make_rep_heatmap runs with replicabilityFactors and GroupColour inputs", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input(list(
    replicabilityFactors = "cond",
    GroupColourcondA     = "red",
    GroupColourcondB     = "blue"
  ))
  fn <- with_mock_env(make_rep_heatmap, input = input)

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(fn(data, title_suffix = "Smoke-with-factor"), NA)
})

test_that("make_rep_heatmap runs with Colours-tab level labels", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input(list(
    replicabilityFactors = "cond",
    GroupColourcondA     = "red",
    GroupColourcondB     = "blue",
    GroupLabelcondA      = "<i>Il4ra</i><sup>fl/fl</sup> &Delta;"
  ))
  fn <- with_mock_env(make_rep_heatmap, input = input)

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(fn(data, title_suffix = "Smoke-with-labels"), NA)
})

test_that("make_rep_heatmap runs with a BrewerPal-derived palette", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input(list(replicabilityPalette = "RdBu"))
  fn <- with_mock_env(make_rep_heatmap, input = input)

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(fn(data, title_suffix = "Smoke-with-palette"), NA)
})


# --- make_rep_dendrogram ----------------------------------------------------
test_that("make_rep_dendrogram runs without error when replicabilityFactors is NULL", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input()
  fn <- with_mock_env(make_rep_dendrogram, input = input)

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(fn(data, title_suffix = "Smoke"), NA)
})

test_that("make_rep_dendrogram runs with replicabilityFactors and GroupColour inputs", {
  skip_on_cran()
  data <- make_rep_data()
  input <- make_rep_input(list(
    replicabilityFactors = "cond",
    GroupColourcondA     = "red",
    GroupColourcondB     = "blue"
  ))
  fn <- with_mock_env(make_rep_dendrogram, input = input)

  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  expect_error(fn(data, title_suffix = "Smoke-with-factor"), NA)
})
