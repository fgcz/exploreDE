# -----------------------------------------------------------------------------
# Count QC tab (helper-app.R), on the demo RNA-seq result without its
# comparison, which exploreDE opens as a CountQC result: plot downloads follow
# the figure size set at the time of each download.
# -----------------------------------------------------------------------------

test_that("Count QC downloads take the figure size set at download time", {
  skip_if_not_installed("png")
  se <- readRDS(system.file("examples", "genomics", "deResult.rds", package = "exploreDE"))
  md <- S4Vectors::metadata(se)
  md$param$comparison <- NULL
  S4Vectors::metadata(se) <- md
  f <- file.path(withr::local_tempdir(), "countqc.rds")
  saveRDS(se, f)

  app <- app_start("countqc-download", options = list(exploreDE.file = f))
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  app_open_tab(app, "countQcTab")
  expect_identical(app_output_errors(app), character(0))
  # A download button gets its link only while shown: open its Download tab
  app$run_js("$('#shiny-tab-countQcTab a[data-value=\"Download\"]').click()")
  app$wait_for_js("!!$('#dlPlot_countQcTotalReads').attr('href')")
  app$set_inputs(countQcDlFormat = "PNG", countQcDlDpi = "72")
  app$wait_for_idle()
  png_width <- function() dim(png::readPNG(app$get_download("dlPlot_countQcTotalReads")))[2]

  expect_equal(png_width(), 800 / 96 * 72, tolerance = 1)
  app$set_inputs(countQcFigWidth = 400)
  app$wait_for_idle()
  expect_equal(png_width(), 400 / 96 * 72, tolerance = 1)
})
