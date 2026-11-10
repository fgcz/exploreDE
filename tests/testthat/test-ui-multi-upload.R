# -----------------------------------------------------------------------------
# Several uploaded .h5ad files compared in Multi DE (helper-app.R): the
# landing page takes more than one file, Multi DE's Data Selection tab adds
# more and refuses a file without contrasts, and the Multi DE tabs render.
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

notifications <- function(app) {
  app$get_js("Array.from(document.querySelectorAll('.shiny-notification')).map(n => n.innerText).join(' | ')")
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
  expect_match(notifications(app), "Added 2 contrasts from 1 file.", fixed = TRUE)

  app$upload_file(multiUploadH5ad = counts_only, timeout_ = 180e3)
  app$wait_for_idle()
  expect_length(multi_datasets(app), 6L)
  expect_match(notifications(app), "counts_only.h5ad has no DE contrasts", fixed = TRUE)
  expect_identical(app_output_errors(app), character(0))

  multi_tabs <- setdiff(grep("^multi", app_menu_tabs(app), value = TRUE), "multiDataSelectionTab")
  expect_true(length(multi_tabs) > 0)
  for (tab in multi_tabs) {
    app_open_tab(app, tab)
    expect_identical(app_output_errors(app), character(0), label = paste("R errors on", tab))
  }
  expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})
