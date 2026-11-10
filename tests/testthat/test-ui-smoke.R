# -----------------------------------------------------------------------------
# The app in headless Chrome (helper-app.R): the demo dataset loads, and every
# tab of its menu opens and renders without an R error.
# -----------------------------------------------------------------------------

test_that("every tab of the demo dataset renders without an error", {
  app <- app_start("smoke")
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
