# -----------------------------------------------------------------------------
# One session on the demo dataset in the Heatmaps tab (helper-app.R): the DE
# heatmap's table download, what a redraw recomputes, and sample ratios. The
# tests run in this order on the same session.
# -----------------------------------------------------------------------------

# Shown notifications whose text contains `text`.
notifications <- function(app, text) {
  app$get_js(sprintf(
    "Array.from(document.querySelectorAll('.shiny-notification')).filter(n => n.innerText.includes('%s')).length",
    text
  ))
}

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

app <- app_session("demo-heatmaps", .local_envir = environment())
app_load_demo(app)

test_that("the DE heatmap download has DE statistics only when asked for", {
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

test_that("a font size change redraws the heatmap without preparing the counts again", {
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
  # Back to the defaults, as a user turns batch correction off again
  app$set_inputs(heatmapBatch = "None", heatmapCounts = "normalised_log2", heatmapLog2 = TRUE)
  app$wait_for_idle(duration = 1500)
})

# The prolfquapp fixture's ids (Protein_1-NoChange1, Protein_1-NoChange2, ...)
# share their first 15 characters; all 40 proteins have one peptide, so the
# peptide filter's default of 2 keeps none.

test_that("a sample ratio is defined, plotted next to a contrast and removed", {
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
