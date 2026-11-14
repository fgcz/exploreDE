# -----------------------------------------------------------------------------
# One session on a prolfquapp 2.10 DEA result
# (fixtures/prolfquapp_2.10_AnnData.h5ad, helper-app.R): every tab renders
# without an R error, Ordination and Heatmap say why they draw nothing when
# the peptide filter leaves too few features, and the DE heatmap table keeps
# full ids. All 40 proteins of the fixture have one peptide; the filter's
# minimum defaults to 2. Their ids (Protein_1-NoChange1, Protein_1-NoChange2,
# ...) share their first 15 characters.
# -----------------------------------------------------------------------------

# The validation message an output shows, or "" when it shows none.
output_message <- function(app, id) {
  app$get_js(sprintf("$('#%s.shiny-output-error-validation').text()", id))
}

skip_if_not_installed("anndataR")
fixture <- test_path("fixtures", "prolfquapp_2.10_AnnData.h5ad")
app <- app_session("prolfquapp", .local_envir = environment())
app$upload_file(uploadH5ad = fixture, timeout_ = 180e3)
app$wait_for_value(output = "settingsTable", timeout = 180e3)
app$wait_for_idle()

test_that("a prolfquapp result opens in every tab; Ordination and Heatmap name the peptide filter", {
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

test_that("the DE heatmap table has the plotted features under their full ids, and none when empty", {
  ids <- new_explore_dataset(fixture)$feature_ids()
  app_open_tab(app, "heatmapTab")
  app$run_js("$('#shiny-tab-heatmapTab a[data-value=\"Download\"]').click()")
  app$set_inputs(nrPeptideHeatmap = 1, includeDEStatsHeatmap = TRUE)
  app$wait_for_idle(duration = 1500)
  app$wait_for_js("!!$('#dlTable_heatmapDE').attr('href')")
  table <- function() openxlsx::read.xlsx(app$get_download("dlTable_heatmapDE"))

  shown <- table()
  expect_gt(nrow(shown), 0)
  expect_true(all(shown$feature_id %in% ids))
  expect_false(anyNA(shown$log2Ratio))

  app$set_inputs(nrPeptideHeatmap = 2)
  app$wait_for_idle(duration = 1500)
  expect_false(any(nzchar(table()$feature_id)))
})
