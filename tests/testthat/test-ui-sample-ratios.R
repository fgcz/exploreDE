# -----------------------------------------------------------------------------
# Heatmaps > Settings > Sample ratios and the Ratios heatmap, on the demo
# dataset (helper-app.R): a ratio is defined, plotted next to the contrast's
# log2 fold change, downloaded, and removed again.
# -----------------------------------------------------------------------------

ratio_list_text <- function(app) {
  app$get_js("document.getElementById('heatmapRatioListUI').innerText")
}

add_ratio <- function(app, name, numerator, denominator) {
  app$set_inputs(
    heatmapRatioName = name,
    heatmapRatioNumerator = numerator,
    heatmapRatioDenominator = denominator
  )
  app$click("heatmapRatioAdd")
  app$wait_for_idle()
}

# The Ratios heatmap's matrix and figure are debounced by a second each.
ratio_heatmap_src <- function(app) {
  app$wait_for_idle(duration = 2500)
  app$get_value(output = "heatmapRatios")$src %||% ""
}

test_that("a sample ratio is defined, plotted next to a contrast and removed", {
  app <- app_session("sample-ratios")
  app_load_demo(app)
  app_open_tab(app, "heatmapTab")
  app$click(selector = "a[data-value='Sample ratios']")
  app$wait_for_idle()
  expect_match(ratio_list_text(app), "No sample ratios defined yet.", fixed = TRUE)

  add_ratio(app, "diff_45", c("S4_diff", "S5_diff"), c("S1_undiff", "S2_undiff"))
  expect_match(ratio_list_text(app), "mean(S4_diff, S5_diff) - mean(S1_undiff, S2_undiff)", fixed = TRUE)
  expect_identical(app$get_value(input = "heatmapRatioPlot"), "diff_45")

  # A name already used, by a ratio or a contrast, is refused
  add_ratio(app, "diff_45", "S6_diff", "S3_undiff")
  add_ratio(app, "diff_vs_undiff", "S6_diff", "S3_undiff")
  # feature_id names the id column of the table download
  add_ratio(app, "feature_id", "S6_diff", "S3_undiff")
  expect_identical(app$get_value(input = "heatmapRatioPlot"), "diff_45")
  expect_no_match(ratio_list_text(app), "S6_diff", fixed = TRUE)

  app$set_inputs(sig = "Ratios")
  expect_match(ratio_heatmap_src(app), "^data:image/png")
  expect_identical(app_output_errors(app), character(0))

  app$set_inputs(heatmapRatioPlot = c("diff_45", "diff_vs_undiff"))
  expect_match(ratio_heatmap_src(app), "^data:image/png")
  # A download link gets its address once its panel is shown
  app$click(selector = "#shiny-tab-heatmapTab a[data-value='Download']")
  app$wait_for_idle()
  values <- openxlsx::read.xlsx(app$get_download("dlTable_heatmapRatios"))
  expect_identical(colnames(values), c("feature_id", "diff_45", "diff_vs_undiff"))
  expect_identical(nrow(values), 50L)
  expect_false(anyNA(values))
  # Fewer than 2 features is not a heatmap: the smallest is 2
  app$set_inputs(heatmapRatioTopN = 1)
  ratio_heatmap_src(app)
  expect_identical(nrow(openxlsx::read.xlsx(app$get_download("dlTable_heatmapRatios"))), 2L)

  app$click(selector = "a[data-value='Sample ratios']")
  app$set_inputs(heatmapRatioRemove = "diff_45")
  app$click("heatmapRatioRemoveBtn")
  app$wait_for_idle()
  expect_match(ratio_list_text(app), "No sample ratios defined yet.", fixed = TRUE)
  expect_identical(app_output_errors(app), character(0))
})

test_that("a ratio on raw counts is the difference of mean log2(count + 1)", {
  app <- app_session("sample-ratio-values")
  app_load_demo(app)
  app_open_tab(app, "heatmapTab")
  app$set_inputs(heatmapCounts = "counts")
  app$wait_for_idle()
  app$click(selector = "a[data-value='Sample ratios']")
  app$wait_for_idle()
  num <- c("S4_diff", "S5_diff")
  den <- c("S1_undiff", "S2_undiff")
  add_ratio(app, "r_counts", num, den)
  app$set_inputs(heatmapRatioPlot = "r_counts", sig = "Ratios")
  ratio_heatmap_src(app)
  app$click(selector = "#shiny-tab-heatmapTab a[data-value='Download']")
  app$wait_for_idle()
  values <- openxlsx::read.xlsx(app$get_download("dlTable_heatmapRatios"))

  eds <- new_explore_dataset(readRDS(system.file("examples", "genomics", "deResult.rds", package = "exploreDE")))
  counts <- log2(eds$assay("counts") + 1)
  expected <- rowMeans(counts[, num]) - rowMeans(counts[, den])
  expect_gt(nrow(values), 0)
  expect_equal(values$r_counts, unname(expected[values$feature_id]), tolerance = 1e-10)
})
