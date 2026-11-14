# -----------------------------------------------------------------------------
# One session on the demo dataset (helper-app.R), as a user explores a result:
# every tab renders, then the Volcano thresholds, the DE QC colour modes, the
# feature-pattern box, Replicability's Log2 option, the Methods download, the
# colour template and saving and restoring the session. The tests run in this
# order on the same session; each leaves the state a user would.
# -----------------------------------------------------------------------------

# Type a pattern; the box is debounced, so wait for longer than that.
set_pattern <- function(app, pattern, regex = FALSE) {
  app$set_inputs(fb_pattern_regex_boxplot = regex, fb_pattern_boxplot = pattern)
  app$wait_for_idle(duration = 1500)
}

bucket <- function(app) {
  unlist(app$get_value(input = "fb_keep_boxplot")) %||% character(0)
}

pattern_info <- function(app) {
  app$get_value(output = "fb_pattern_info_boxplot")
}

app <- app_session("demo-explore", .local_envir = environment())

test_that("every tab of the demo dataset renders without an error", {
  expect_identical(app_menu_tabs(app), "landingTab")

  app_load_demo(app)
  tabs <- app_menu_tabs(app)
  expect_true(all(c("summaryTab", "heatmapTab", "boxplotTab", "volcanoTab", "oraTab", "gseaTab") %in% tabs))
  for (tab in tabs) {
    app_open_tab(app, tab)
    expect_identical(app$get_value(input = "tabs"), tab)
    expect_identical(app_output_errors(app), character(0), label = paste("R errors on", tab))
  }
  expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})

test_that("the Volcano tab names the contrast and classes every point in its table", {
  app_open_tab(app, "volcanoTab")
  app$run_js("$('#volcanoOverviewTable').closest('.box').find('[data-widget=\"collapse\"]').click()")
  app$wait_for_idle(duration = 1500)
  expect_identical(app$get_text("#volcanoDesign"), "diff_vs_undiff")

  statuses <- function() {
    unlist(app$get_js(paste0(
      "Array.from(new Set($('#volcanoOverviewTable').find('table').DataTable()",
      ".column(-1).data().toArray()))"
    )))
  }
  expect_true("FoldChangeSignificantUp" %in% statuses())
  # A fold-change threshold no ratio reaches leaves no fold-change classes
  app$set_inputs(lfcVolcano = 100)
  app$wait_for_idle(duration = 1500)
  expect_false(any(grepl("^FoldChange", statuses())))
  expect_identical(app_output_errors(app), character(0))
})

test_that("the DE QC scatter draws in each colour-by mode", {
  app_open_tab(app, "deQcTab")
  modes <- unlist(app$get_js(
    "Array.from(document.querySelectorAll('input[name=\"deQcScatterColourBy\"]')).map(e => e.value)"
  ))
  expect_true("Standard" %in% modes)
  for (mode in modes) {
    app$set_inputs(deQcScatterColourBy = mode)
    app$wait_for_idle(duration = 1500)
    expect_identical(app_output_errors(app), character(0), label = paste("R errors in mode", mode))
  }
})

test_that("a pattern fills the feature bucket with its matches", {
  app_open_tab(app, "boxplotTab")

  set_pattern(app, "hspa1*")
  expect_identical(pattern_info(app), "7 names match.")
  expect_setequal(bucket(app), c("HSPA12A", "HSPA14", "HSPA1A", "HSPA1B", "HSPA1L", "HSPA13", "HSPA12B"))
  expect_identical(app_output_errors(app), character(0))

  # A new pattern replaces the old matches
  set_pattern(app, "CXCL?")
  expect_setequal(bucket(app), c("CXCL6", "CXCL5", "CXCL1", "CXCL3", "CXCL2", "CXCL8", "CXCL9"))

  set_pattern(app, "^CXCL1[0-9]$", regex = TRUE)
  expect_identical(pattern_info(app), "7 names match.")
  expect_setequal(bucket(app), c("CXCL16", "CXCL17", "CXCL10", "CXCL14", "CXCL13", "CXCL11", "CXCL12"))

  set_pattern(app, "CXCL1(", regex = TRUE)
  expect_identical(pattern_info(app), "Not a valid pattern.")
  expect_identical(bucket(app), character(0))

  set_pattern(app, "hspa1a")
  expect_identical(bucket(app), "HSPA1A")
  set_pattern(app, "")
  expect_identical(pattern_info(app), "")
  expect_identical(bucket(app), character(0))
  expect_identical(app_output_errors(app), character(0))
})

test_that("Replicability offers Log2 for every count method and logs linear ones only", {
  app_open_tab(app, "replicabilityTab")
  app$wait_for_js("!!$('#replicabilityDownloadCorr').attr('href')")
  corr <- function() openxlsx::read.xlsx(app$get_download("replicabilityDownloadCorr"))

  app$set_inputs(replicabilityCounts = "counts")
  app$wait_for_idle()
  expect_true(app$get_js("$('#replicabilityCountsLog').is(':visible')"))

  app$set_inputs(replicabilityCounts = "normalised_log2")
  app$wait_for_idle()
  expect_true(app$get_js("$('#replicabilityCountsLog').is(':visible')"))
  as_is <- corr()
  app$set_inputs(replicabilityCountsLog = TRUE)
  app$wait_for_idle()
  expect_equal(corr(), as_is)

  app$set_inputs(replicabilityCounts = "counts")
  app$wait_for_idle()
  logged <- corr()
  app$set_inputs(replicabilityCountsLog = FALSE)
  app$wait_for_idle()
  expect_false(isTRUE(all.equal(corr(), logged)))
  expect_identical(app_output_errors(app), character(0))
})

test_that("the Methods download of the demo is a Word document with its sections", {
  skip_if_not_installed("officer")
  app_open_tab(app, "summaryTab")
  # A download button gets its link only while shown: open the Methods tab
  app$run_js("$('#shiny-tab-summaryTab a[data-value=\"Methods\"]').click()")
  app$wait_for_js("!!$('#dlMethods_summary').attr('href')")

  text <- officer::docx_summary(officer::read_docx(app$get_download("dlMethods_summary")))$text
  expect_true(all(
    c("Methods", "Read Mapping and Quantification", "Differential Expression Analysis", "References") %in% text
  ))
  expect_true(any(grepl("DESeq2", text)))
})

test_that("the colour template round-trips through download and import", {
  app_open_tab(app, "colourTab")
  app$wait_for_js("!!$('#dlColourTemplate').attr('href')")

  template <- openxlsx::read.xlsx(app$get_download("dlColourTemplate"))
  expect_identical(colnames(template), c("Level", "Colour", "Label"))
  expect_gt(nrow(template), 0)
  expect_true(all(grepl("^GroupColour", template$Level)))

  template$Colour[1] <- "#123456"
  f <- file.path(withr::local_tempdir(), "colours.xlsx")
  openxlsx::write.xlsx(template, f)
  app$upload_file(importColourFile = f)
  app$wait_for_idle()
  expect_identical(tolower(app$get_value(input = template$Level[1])), "#123456")
  expect_identical(app_output_errors(app), character(0))
})

test_that("a saved session restores its settings on the same dataset", {
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
