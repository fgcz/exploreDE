# -----------------------------------------------------------------------------
# The GSEA tab of an uploaded .h5ad with GSEA results on one contrast only
# (helper-app.R): switching to the other contrast shows the "no results"
# message, and switching back shows the table again.
# -----------------------------------------------------------------------------

gsea_pane_text <- function(app) {
  app$get_js("document.getElementById('shiny-tab-gseaTab').innerText")
}

test_that("the GSEA table comes back after a contrast without results", {
  skip_if_not_installed("anndataR")
  dds <- make_toy_dds(n_genes = 80)
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(
    dds,
    out = h5ad,
    omics = "rnaseq",
    de = list(
      control_vs_treated = DESeq2::results(dds, contrast = c("condition", "control", "treated")),
      treated_vs_control = DESeq2::results(dds, contrast = c("condition", "treated", "control"))
    ),
    pathways = list(control_vs_treated = list(gsea = list(BP = .mini_cp("gseaResult", 3)))),
    organism = "Homo sapiens"
  )

  app <- app_session("gsea")
  app$upload_file(uploadH5ad = h5ad)
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  app_open_tab(app, "gseaTab")
  app$set_inputs(contrastSelected_gsea = "control_vs_treated")
  app$wait_for_idle()
  expect_match(gsea_pane_text(app), "term 1", fixed = TRUE)

  app$set_inputs(contrastSelected_gsea = "treated_vs_control")
  app$wait_for_idle()
  app$run_js("$('.modal').modal('hide')")
  expect_match(gsea_pane_text(app), "No GSEA results for this contrast", fixed = TRUE)

  app$set_inputs(contrastSelected_gsea = "control_vs_treated")
  app$wait_for_idle()
  expect_match(gsea_pane_text(app), "term 1", fixed = TRUE)
  expect_identical(app_output_errors(app), character(0))
})
