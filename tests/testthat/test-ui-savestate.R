# -----------------------------------------------------------------------------
# Save / Load Session (helper-app.R): a setting saved to the .json comes back
# when the file is loaded again, and the tab says what was restored.
# -----------------------------------------------------------------------------

test_that("a saved session restores its settings on the same dataset", {
  app <- app_session("save-restore")
  app_load_demo(app)
  app_open_tab(app, "heatmapTab")
  app$set_inputs(lfcHeatmap = 1.5)
  app$wait_for_idle()

  app_open_tab(app, "saveStateTab")
  app$wait_for_js("!!$('#download_state_json').attr('href')")
  saved <- app$get_download("download_state_json")
  expect_identical(jsonlite::fromJSON(saved)$inputs$lfcHeatmap, 1.5)

  app_open_tab(app, "heatmapTab")
  app$set_inputs(lfcHeatmap = 0)
  app$wait_for_idle()
  app_open_tab(app, "saveStateTab")
  app$upload_file(load_state_file = saved)
  app$wait_for_idle()
  expect_match(app$get_value(output = "load_state_status"), "Restored")
  app_open_tab(app, "heatmapTab")
  app$wait_for_idle()
  expect_identical(app$get_value(input = "lfcHeatmap"), 1.5)
  expect_identical(app_output_errors(app), character(0))
})

test_that("a session of a proteomics dataset is saved", {
  f <- system.file("examples", "proteomics", "SummarizedExperiment.rds", package = "exploreDE")
  app <- app_session("save-proteomics", options = list(exploreDE.file = f))
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  app_open_tab(app, "saveStateTab")
  app$wait_for_js("!!$('#download_state_json').attr('href')")
  saved <- jsonlite::fromJSON(app$get_download("download_state_json"))
  expect_identical(saved$meta$fingerprint$omics, "proteomics")
  expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})

test_that("restoring a session brings back the loaded dataset's own contrast in Multi DE", {
  skip_if_not_installed("anndataR")
  dir <- withr::local_tempdir()
  dds <- make_toy_dds(n_genes = 80, seed = 1)
  main <- file.path(dir, "main.h5ad")
  build_explore_h5ad(
    dds,
    out = main,
    omics = "rnaseq",
    de = list(treated_vs_control = DESeq2::results(dds, contrast = c("condition", "treated", "control"))),
    organism = "Homo sapiens"
  )
  other <- file.path(dir, "B.h5ad")
  file.copy(main, other)

  # The saving session is closed (local()) before the restoring one opens, as
  # when a user saves, closes the app and comes back later.
  saved <- local({
    app <- app_session("save-own-contrast")
    app$upload_file(uploadH5ad = main, timeout_ = 180e3)
    app$wait_for_value(output = "settingsTable", timeout = 180e3)
    app$wait_for_idle()
    app_open_tab(app, "multiDataSelectionTab")
    app$upload_file(multiUploadH5ad = other, timeout_ = 180e3)
    app$wait_for_idle()
    expect_setequal(multi_datasets(app), c("treated_vs_control", "B_treated_vs_control"))
    app_open_tab(app, "saveStateTab")
    app$wait_for_js("!!$('#download_state_json').attr('href')")
    app$get_download("download_state_json")
  })

  fresh <- app_session("restore-own-contrast")
  fresh$upload_file(uploadH5ad = main, timeout_ = 180e3)
  fresh$wait_for_value(output = "settingsTable", timeout = 180e3)
  fresh$wait_for_idle()
  app_open_tab(fresh, "saveStateTab")
  fresh$upload_file(load_state_file = saved)
  fresh$wait_for_idle()
  status <- fresh$get_value(output = "load_state_status")
  expect_match(status, "1 of 2 datasets restored", fixed = TRUE)
  expect_match(status, "add B again", fixed = TRUE)
  app_open_tab(fresh, "multiDataSelectionTab")
  expect_identical(multi_datasets(fresh), "treated_vs_control")
})
