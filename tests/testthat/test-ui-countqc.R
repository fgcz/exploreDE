# -----------------------------------------------------------------------------
# Count QC tab (helper-app.R), on an uploaded CRISPR guide-count file, which
# has Count QC and no DE tabs: plot downloads follow the figure size set at
# the time of each download.
# -----------------------------------------------------------------------------

test_that("Count QC downloads take the figure size set at download time", {
  skip_if_not_installed("png")
  skip_if_not_installed("anndataR")
  p <- make_toy_pieces()
  f <- file.path(withr::local_tempdir(), "guides.h5ad")
  build_explore_h5ad(p$counts, out = f, omics = "crispr", coldata = p$coldata)

  app <- app_session("countqc-download")
  app$upload_file(uploadH5ad = f, timeout_ = 180e3)
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  expect_false("volcanoTab" %in% app_menu_tabs(app))
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
