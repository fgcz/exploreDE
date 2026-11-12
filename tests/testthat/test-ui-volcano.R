# -----------------------------------------------------------------------------
# Volcano (helper-app.R): the contrast label and the Results Summary table,
# whose classes follow the thresholds.
# -----------------------------------------------------------------------------

test_that("the Volcano tab names the contrast and classes every point in its table", {
  app <- app_start("volcano")
  app_load_demo(app)
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
