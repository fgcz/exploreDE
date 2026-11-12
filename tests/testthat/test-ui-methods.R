# -----------------------------------------------------------------------------
# Results Summary > Downloads > Methods (helper-app.R): the demo RNA-seq result
# gives a Word document with the analysis sections and references.
# -----------------------------------------------------------------------------

test_that("the Methods download of the demo is a Word document with its sections", {
  skip_if_not_installed("officer")
  app <- app_start("methods-download")
  app_load_demo(app)
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
