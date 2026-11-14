# -----------------------------------------------------------------------------
# A prolfquapp 2.10 DEA result (fixtures/prolfquapp_2.10_AnnData.h5ad) in the
# app (helper-app.R): every tab renders without an R error, and Ordination
# says why it draws nothing when the peptide filter leaves too few features.
# All 40 proteins of the fixture have one peptide; the filter's minimum
# defaults to 2.
# -----------------------------------------------------------------------------

# The validation message an output shows, or "" when it shows none.
output_message <- function(app, id) {
  app$get_js(sprintf("$('#%s.shiny-output-error-validation').text()", id))
}

test_that("a prolfquapp result opens in every tab; Ordination and Heatmap name the peptide filter", {
  skip_if_not_installed("anndataR")
  app <- app_session("prolfquapp")
  app$upload_file(uploadH5ad = test_path("fixtures", "prolfquapp_2.10_AnnData.h5ad"), timeout_ = 180e3)
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()

  for (tab in app_menu_tabs(app)) {
    app_open_tab(app, tab)
    expect_identical(app_output_errors(app), character(0), label = paste("R errors on", tab))
  }

  app_open_tab(app, "ordinationTab")
  expect_match(output_message(app, "pcaStatic"), "Minimum number of peptides")

  app$set_inputs(nrPeptidePCA = 1)
  app$wait_for_idle()
  expect_identical(output_message(app, "pcaStatic"), "")
  expect_false(is.null(app$get_value(output = "pcaStatic")$src))

  app_open_tab(app, "heatmapTab")
  expect_match(
    app$get_js("document.getElementById('heatmapDesignBoth Directions').innerText"),
    "'Minimum number of peptides' \\(2\\) keeps 0 of 40 features"
  )
  expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})
