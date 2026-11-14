# -----------------------------------------------------------------------------
# Several uploaded .h5ad files compared in Multi DE (helper-app.R): the
# landing page takes more than one file, Multi DE's Data Selection tab adds
# more and refuses a file without contrasts, and the Multi DE tabs render.
# exploreDE(c(...)) does the same with launch files (test-ui-launch.R).
# -----------------------------------------------------------------------------

test_that("the contrasts of several uploaded files are compared in Multi DE", {
  skip_if_not_installed("anndataR")
  dir <- withr::local_tempdir()
  design1 <- toy_h5ad(dir, "design1", seed = 1)
  design2 <- toy_h5ad(dir, "design2", seed = 2)
  design3 <- toy_h5ad(dir, "design3", seed = 3)
  counts_only <- file.path(dir, "counts_only.h5ad")
  p <- make_toy_pieces()
  build_explore_h5ad(p$counts, out = counts_only, omics = "rnaseq", coldata = p$coldata)
  own <- c("control_vs_treated", "treated_vs_control")

  app <- app_session("multi-upload")
  app$upload_file(uploadH5ad = c(design1, design2), timeout_ = 180e3)
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  expect_true("multiDataSelectionTab" %in% app_menu_tabs(app))

  app_open_tab(app, "multiDataSelectionTab")
  expect_setequal(multi_datasets(app), c(own, paste0("design2_", own)))

  app$upload_file(multiUploadH5ad = design3, timeout_ = 180e3)
  app$wait_for_idle()
  expect_setequal(multi_datasets(app), c(own, paste0("design2_", own), paste0("design3_", own)))

  # A file without contrasts adds nothing
  app$upload_file(multiUploadH5ad = counts_only, timeout_ = 180e3)
  app$wait_for_idle()
  expect_setequal(multi_datasets(app), c(own, paste0("design2_", own), paste0("design3_", own)))
  expect_identical(app_output_errors(app), character(0))

  multi_tabs <- setdiff(grep("^multi", app_menu_tabs(app), value = TRUE), "multiDataSelectionTab")
  expect_true(length(multi_tabs) > 0)
  for (tab in multi_tabs) {
    app_open_tab(app, tab)
    expect_identical(app_output_errors(app), character(0), label = paste("R errors on", tab))
  }
  expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})

# The same contrast loaded twice gives identical tile columns; ComplexHeatmap
# cannot cluster those ("invalid 'yscale' in viewport").
test_that("DE Tile Plots draw the same contrast loaded twice", {
  skip_if_not_installed("anndataR")
  dir <- withr::local_tempdir()
  dds <- make_toy_dds(n_genes = 80, seed = 1)
  one <- file.path(dir, "one.h5ad")
  build_explore_h5ad(
    dds,
    out = one,
    omics = "rnaseq",
    de = list(treated_vs_control = DESeq2::results(dds, contrast = c("condition", "treated", "control"))),
    organism = "Homo sapiens"
  )
  copy <- file.path(dir, "copy.h5ad")
  file.copy(one, copy)

  app <- app_session("multi-tile-identical")
  app$upload_file(uploadH5ad = c(one, copy), timeout_ = 180e3)
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  app_open_tab(app, "multiDataSelectionTab")
  expect_setequal(multi_datasets(app), c("treated_vs_control", "copy_treated_vs_control"))

  app_open_tab(app, "multiTilePlotTab")
  expect_identical(app_output_errors(app), character(0))
  expect_false(is.null(app$get_value(output = "TopNTileMap")$src))
})

# An .rds is opened with readRDS(), which can run code: uploads take .h5ad only.
test_that("an uploaded .rds is refused, on the landing page and in Multi DE", {
  skip_if_not_installed("anndataR")
  rds <- file.path(withr::local_tempdir(), "deResult.rds")
  file.copy(system.file("examples", "genomics", "deResult.rds", package = "exploreDE"), rds)
  app <- app_session("upload-rds")
  app$upload_file(uploadH5ad = rds, timeout_ = 180e3)
  app$wait_for_idle()
  expect_identical(app_menu_tabs(app), "landingTab")

  design1 <- toy_h5ad(withr::local_tempdir(), "design1", seed = 1)
  app$upload_file(uploadH5ad = design1, timeout_ = 180e3)
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  app_open_tab(app, "multiDataSelectionTab")
  before <- multi_datasets(app)
  app$upload_file(multiUploadH5ad = rds, timeout_ = 180e3)
  app$wait_for_idle()
  expect_setequal(multi_datasets(app), before)
})
