# -----------------------------------------------------------------------------
# DE QC (helper-app.R): the between-group scatter draws in every colour-by
# mode the dataset offers.
# -----------------------------------------------------------------------------

test_that("the DE QC scatter draws in each colour-by mode", {
  app <- app_session("deqc-modes")
  app_load_demo(app)
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
