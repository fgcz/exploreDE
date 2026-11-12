# -----------------------------------------------------------------------------
# Multi DE Bubble (helper-app.R): the pathway table renders for ORA and GSEA,
# GSEA with its enrichment-score columns.
# -----------------------------------------------------------------------------

# Bubble plots compare the pathway results of at least two datasets: the demo
# (with ORA and GSEA) opened twice, as exploreDE(c(...)) does.
test_that("the Multi DE Bubble table shows ORA and GSEA results", {
  demo <- system.file("examples", "genomics", "deResult.rds", package = "exploreDE")
  copy <- file.path(withr::local_tempdir(), "copy.rds")
  file.copy(demo, copy)
  app <- app_start("multi-bubble", options = list(exploreDE.file = c(demo, copy)))
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  expect_true("multiBubbleTab" %in% app_menu_tabs(app))
  # The pathway lists are built on the Upset/Venn tab
  app_open_tab(app, "multiUpsetPlotTab")
  app_open_tab(app, "multiBubbleTab")
  headers <- function() {
    unlist(app$get_js("Array.from(document.querySelectorAll('#BubbleTable thead th')).map(th => th.innerText)"))
  }
  rows <- function() app$get_js("document.querySelectorAll('#BubbleTable tbody tr').length")

  app$set_inputs(dataToBubble = "ORA")
  app$set_inputs(bubblePathwayNumber = 5)
  app$wait_for_idle(duration = 1500)
  expect_gt(rows(), 0)
  expect_true(all(c("pvalue", "p.adjust") %in% headers()))

  app$set_inputs(dataToBubble = "GSEA")
  app$wait_for_idle(duration = 1500)
  expect_gt(rows(), 0)
  expect_true(all(c("enrichmentScore", "NES") %in% headers()))
  expect_identical(app_output_errors(app), character(0))
})
