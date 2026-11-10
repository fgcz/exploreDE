# -----------------------------------------------------------------------------
# UI tests (test-ui-*.R): the app in headless Chrome, driven by shinytest2.
# app_start() runs inst/app: the source tree under devtools::test(), the
# installed app under R CMD check. shinytest2's AppDriver skips on CRAN.
# R CMD check on CI skips them too: on the GitHub runners a dataset load
# sometimes never finishes, for reasons not yet known. The UI test workflow
# (.github/workflows/ui-tests.yaml) runs them with EXPLOREDE_UI_TESTS=true,
# and with EXPLOREDE_UI_ARTIFACTS set, which keeps the app's log, a
# screenshot and the shown tab of every test there.
# -----------------------------------------------------------------------------

# Start the app; it is stopped when the calling test ends. `options` are set
# in the app's R process (exploreDE.file, as exploreDE(file) sets it).
app_start <- function(name, options = list()) {
  if (!identical(Sys.getenv("EXPLOREDE_UI_TESTS"), "true")) {
    skip_on_ci()
  }
  skip_if_not_installed("shinytest2")
  # Kept for inspection: every message between app and browser, including
  # when each output starts and finishes recalculating.
  if (nzchar(Sys.getenv("EXPLOREDE_UI_ARTIFACTS"))) {
    options$shiny.trace <- TRUE
  }
  app <- shinytest2::AppDriver$new(
    system.file("app", package = "exploreDE"),
    name = name,
    load_timeout = 180e3,
    timeout = 60e3,
    options = options
  )
  withr::defer(
    {
      app_save_state(app, name)
      app$stop()
    },
    envir = parent.frame()
  )
  app
}

# With EXPLOREDE_UI_ARTIFACTS set, write the app's log (R and browser), the
# shown tab and a screenshot to that directory, as <name>-log.txt,
# <name>-state.json and <name>.png.
app_save_state <- function(app, name) {
  dir <- Sys.getenv("EXPLOREDE_UI_ARTIFACTS")
  if (!nzchar(dir)) {
    return(invisible())
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(utils::capture.output(print(app$get_logs())), file.path(dir, paste0(name, "-log.txt")))
  writeLines(
    app$get_js(paste0(
      "(() => { const v = Shiny.shinyapp.$inputValues;",
      " const get = k => v[Object.keys(v).find(n => n.split(':')[0] === k)];",
      " return JSON.stringify({tabs: get('tabs'),",
      " shown: $('.tab-pane.active[id^=shiny-tab-]').attr('id'),",
      " settingsTable_hidden: get('.clientdata_output_settingsTable_hidden'),",
      " waiter_overlays: document.querySelectorAll('.waiter-overlay').length,",
      " busy: $('html').hasClass('shiny-busy'),",
      " connected: Shiny.shinyapp.isConnected()}); })()"
    )),
    file.path(dir, paste0(name, "-state.json"))
  )
  app$get_screenshot(file.path(dir, paste0(name, ".png")))
  # The pid names the app's R profile, Rprof-<pid>.out, when ui-tests.yaml
  # turns profiling on.
  pid <- app$.__enclos_env__$private$shiny_process$get_pid()
  writeLines(c(paste("pid", pid), app_proc_snapshot(pid)), file.path(dir, paste0(name, "-proc.txt")))
  invisible()
}

# On Linux, two snapshots 5 s apart of the app's R process from /proc: its
# resident and peak memory, per thread its state, user and system CPU ticks
# and kernel wait function, and the child processes. Tells a blocked process (and what it waits on) from a
# busy one. Empty elsewhere.
app_proc_snapshot <- function(pid) {
  if (!dir.exists(file.path("/proc", pid))) {
    return(character(0))
  }
  script <- sprintf(
    paste(
      "for i in 1 2; do date +%%T; grep -E '^(VmRSS|VmHWM)' /proc/%1$d/status;",
      "for t in /proc/%1$d/task/*; do",
      "echo \"$(basename $t) $(cat $t/comm) state/utime/stime=$(cut -d' ' -f3,14,15 $t/stat) wchan=$(cat $t/wchan)\";",
      "done;",
      "ps -o pid,stat,pcpu,etime,wchan:32,args --ppid %1$d;",
      "if [ $i = 1 ]; then sleep 5; fi;",
      "done"
    ),
    as.integer(pid)
  )
  system2("sh", c("-c", shQuote(script)), stdout = TRUE, stderr = TRUE)
}

# Press the landing page's demo button and wait until the Results Summary shows.
app_load_demo <- function(app) {
  app$click("loadDemo")
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  invisible(app)
}

# The tabName of every entry of the sidebar menu. The menu is an output, drawn
# after the app has started, so wait for it.
app_menu_tabs <- function(app) {
  app$wait_for_js("document.querySelectorAll('.sidebar-menu a[data-value]').length > 0")
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
