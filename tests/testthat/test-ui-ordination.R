# -----------------------------------------------------------------------------
# Ordination (helper-app.R): "Log2 counts?" logs raw proteomics intensities
# before batch correction, as it does RNA-seq counts.
# -----------------------------------------------------------------------------

test_that("Ordination logs raw intensities before batch correction", {
  skip_if_not_installed("anndataR")
  pr <- make_toy_prot("protein")
  set.seed(7)
  x <- matrix(rnorm(length(pr$intensities), 20, 1), nrow(pr$intensities), dimnames = dimnames(pr$intensities))
  coldata <- data.frame(row.names = colnames(x), group = rep(c("A", "B"), each = 3), batch = rep(c("x", "y"), 3))
  path <- file.path(withr::local_tempdir(), "prot.h5ad")
  # log2(raw + 1) is the transformed layer
  build_explore_h5ad(
    list(raw = 2^x - 1, transformed = x),
    out = path,
    omics = "proteomics",
    coldata = coldata,
    rowdata = pr$rowdata
  )

  app <- app_session("ordination-log2")
  app$upload_file(uploadH5ad = path, timeout_ = 180e3)
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  app_open_tab(app, "ordinationTab")
  app$set_inputs(pcaBatch = "batch", pcaBatchDesign = "group")
  app$wait_for_idle()
  app$run_js("$('#shiny-tab-ordinationTab a[data-value=\"Download\"]').click()")
  app$wait_for_js("!!$('#dlTable_ordinationCoords').attr('href')")
  coords <- function() openxlsx::read.xlsx(app$get_download("dlTable_ordinationCoords"))

  expect_true(app$get_value(input = "pcaLog2"))
  transformed <- coords()
  app$set_inputs(pcaCounts = "intensities_raw")
  app$wait_for_idle()
  expect_equal(coords(), transformed)
  expect_identical(app_output_errors(app), character(0))
})
