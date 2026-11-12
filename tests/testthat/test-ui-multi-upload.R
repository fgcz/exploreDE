# -----------------------------------------------------------------------------
# Several uploaded .h5ad files compared in Multi DE (helper-app.R): the
# landing page takes more than one file, Multi DE's Data Selection tab adds
# more and refuses a file without contrasts, and the Multi DE tabs render.
# exploreDE(c(...)) does the same with any file the app opens, an .rds too.
# -----------------------------------------------------------------------------

# A two-contrast RNA-seq .h5ad named <name>.h5ad in `dir`.
toy_h5ad <- function(dir, name, seed) {
  dds <- make_toy_dds(n_genes = 80, seed = seed)
  out <- file.path(dir, paste0(name, ".h5ad"))
  build_explore_h5ad(
    dds,
    out = out,
    omics = "rnaseq",
    de = list(
      control_vs_treated = DESeq2::results(dds, contrast = c("condition", "control", "treated")),
      treated_vs_control = DESeq2::results(dds, contrast = c("condition", "treated", "control"))
    ),
    organism = "Homo sapiens"
  )
  out
}

# The Multi DE datasets, as listed on the Data Selection tab.
multi_datasets <- function(app) {
  unlist(app$get_value(input = "multiContrastOrder")) %||% character(0)
}

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

  app <- app_start("multi-upload")
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

  app <- app_start("multi-tile-identical")
  app$upload_file(uploadH5ad = c(one, copy), timeout_ = 180e3)
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  app_open_tab(app, "multiDataSelectionTab")
  expect_setequal(multi_datasets(app), c("treated_vs_control", "copy_treated_vs_control"))

  app_open_tab(app, "multiTilePlotTab")
  expect_identical(app_output_errors(app), character(0))
  expect_false(is.null(app$get_value(output = "TopNTileMap")$src))
})

# exploreDE(c(...)) sets the option exploreDE.file to its paths; set here
# directly in the app's process.
test_that("exploreDE() with several files opens the first and compares the others' contrasts", {
  skip_if_not_installed("anndataR")
  design1 <- toy_h5ad(withr::local_tempdir(), "design1", seed = 1)
  examples <- system.file("examples", package = "exploreDE")
  files <- c(
    design1,
    file.path(examples, "genomics", "deResult.rds"),
    file.path(examples, "proteomics", "SummarizedExperiment.rds")
  )

  app <- app_start("multi-files-option", options = list(exploreDE.file = files))
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  expect_true("multiDataSelectionTab" %in% app_menu_tabs(app))

  # The SummarizedExperiment .rds joins; the proteomics one is refused, since
  # the opened file holds RNA-seq data.
  app_open_tab(app, "multiDataSelectionTab")
  expect_setequal(multi_datasets(app), c("control_vs_treated", "treated_vs_control", "deResult_diff_vs_undiff"))
  expect_identical(app_output_errors(app), character(0))
  expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})

# An .rds is opened with readRDS(), which can run code: uploads take .h5ad only.
test_that("an uploaded .rds is refused, on the landing page and in Multi DE", {
  skip_if_not_installed("anndataR")
  rds <- file.path(withr::local_tempdir(), "deResult.rds")
  file.copy(system.file("examples", "genomics", "deResult.rds", package = "exploreDE"), rds)
  app <- app_start("upload-rds")
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

test_that("exploreDE() with several files opens the first that loads and compares the rest", {
  skip_if_not_installed("anndataR")
  dir <- withr::local_tempdir()
  bad <- file.path(dir, "bad.h5ad")
  writeLines("not an h5ad file", bad)
  design1 <- toy_h5ad(dir, "design1", seed = 1)
  design2 <- toy_h5ad(dir, "design2", seed = 2)

  app <- app_start("multi-files-first-bad", options = list(exploreDE.file = c(bad, design1, design2)))
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  expect_match(app$get_js("document.querySelector('.modal-body').innerText"), "bad.h5ad", fixed = TRUE)
  app_open_tab(app, "multiDataSelectionTab")
  own <- c("control_vs_treated", "treated_vs_control")
  expect_setequal(multi_datasets(app), c(own, paste0("design2_", own)))
})
