# -----------------------------------------------------------------------------
# UI tests (test-ui-*.R): the app in headless Chrome, driven by shinytest2.
# app_start() runs inst/app: the source tree under devtools::test(), the
# installed app under R CMD check. shinytest2's AppDriver skips on CRAN.
# -----------------------------------------------------------------------------

# Start the app; it is stopped when the calling test ends.
app_start <- function(name) {
  skip_if_not_installed("shinytest2")
  app <- shinytest2::AppDriver$new(
    system.file("app", package = "exploreDE"),
    name = name,
    load_timeout = 180e3,
    timeout = 60e3
  )
  withr::defer(app$stop(), envir = parent.frame())
  app
}

# Press the landing page's demo button and wait until the Results Summary shows.
app_load_demo <- function(app) {
  app$click("loadDemo")
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  invisible(app)
}

# The tabName of every entry of the sidebar menu.
app_menu_tabs <- function(app) {
  unlist(app$get_js(
    "Array.from(document.querySelectorAll('.sidebar-menu a[data-value]')).map(a => a.dataset.value)"
  ))
}

app_open_tab <- function(app, tab) {
  app$set_inputs(tabs = tab)
  app$wait_for_idle()
  invisible(app)
}

# The outputs showing an R error, as "<id>: <message>". A validate() or req()
# message is not an error. An htmlwidget (DT table, plotly) shows its error in
# an element without id right after the widget, so that one is named by the
# widget's id.
app_output_errors <- function(app) {
  unlist(app$get_js(paste0(
    "Array.from(document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)'))",
    ".map(e => (e.id || e.previousElementSibling.id) + ': ' + e.innerText)"
  ))) %||%
    character(0)
}
