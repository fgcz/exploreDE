# -----------------------------------------------------------------------------
# Heatmap table downloads (helper-app.R): the DE statistics are in the file
# only when "Include DE statistics" is ticked.
# -----------------------------------------------------------------------------

test_that("the DE heatmap download has DE statistics only when asked for", {
  app <- app_start("heatmap-export")
  app_load_demo(app)
  app_open_tab(app, "heatmapTab")
  # A download button gets its link only while shown: open Settings > Download
  app$run_js("$('#shiny-tab-heatmapTab a[data-value=\"Download\"]').click()")
  app$wait_for_js("!!$('#dlTable_heatmapDE').attr('href')")

  statistics <- c("log2Ratio", "pValue", "fdr")
  without <- openxlsx::read.xlsx(app$get_download("dlTable_heatmapDE"))
  expect_true("feature_id" %in% colnames(without))
  expect_gt(nrow(without), 0)
  expect_length(intersect(statistics, colnames(without)), 0)

  app$set_inputs(includeDEStatsHeatmap = TRUE)
  app$wait_for_idle()
  with <- openxlsx::read.xlsx(app$get_download("dlTable_heatmapDE"))
  expect_true(all(statistics %in% colnames(with)))
  expect_identical(with$feature_id, without$feature_id)
})
