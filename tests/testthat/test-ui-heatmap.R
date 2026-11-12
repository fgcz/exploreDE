# -----------------------------------------------------------------------------
# Heatmaps (helper-app.R): what a redraw recomputes, and the table downloads.
# -----------------------------------------------------------------------------

# Shown notifications whose text contains `text`.
notifications <- function(app, text) {
  app$get_js(sprintf(
    "Array.from(document.querySelectorAll('.shiny-notification')).filter(n => n.innerText.includes('%s')).length",
    text
  ))
}

test_that("a font size change redraws the heatmap without preparing the counts again", {
  app <- app_start("heatmap-prepared")
  app_load_demo(app)
  app_open_tab(app, "heatmapTab")
  # Batch correction on raw counts is refused with a message each time the
  # counts are prepared
  app$set_inputs(heatmapCounts = "counts", heatmapLog2 = FALSE, heatmapBatch = "Condition")
  app$wait_for_idle(duration = 1500)
  expect_gt(notifications(app, "requires log-scale"), 0)

  app$run_js("document.querySelectorAll('.shiny-notification').forEach(n => n.remove())")
  app$set_inputs(textSizeHeatmap = 13)
  app$wait_for_idle(duration = 1500)
  expect_equal(notifications(app, "requires log-scale"), 0)
  expect_identical(app_output_errors(app), character(0))
})

# The prolfquapp fixture's ids (Protein_1-NoChange1, Protein_1-NoChange2, ...)
# share their first 15 characters; all 40 proteins have one peptide, so the
# peptide filter's default of 2 keeps none.
test_that("the DE heatmap table has the plotted features under their full ids, and none when empty", {
  skip_if_not_installed("anndataR")
  fixture <- test_path("fixtures", "prolfquapp_2.10_AnnData.h5ad")
  ids <- new_explore_dataset(fixture)$feature_ids()
  app <- app_start("heatmap-table")
  app$upload_file(uploadH5ad = fixture, timeout_ = 180e3)
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
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
