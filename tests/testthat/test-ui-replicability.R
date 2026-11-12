# -----------------------------------------------------------------------------
# Replicability (helper-app.R): "Log2 Counts?" is offered for every count
# method and does not log a count method that is already log2.
# -----------------------------------------------------------------------------

test_that("Replicability offers Log2 for every count method and logs linear ones only", {
  app <- app_start("replicability-log2")
  app_load_demo(app)
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
