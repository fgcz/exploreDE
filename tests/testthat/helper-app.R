# -----------------------------------------------------------------------------
# UI tests (test-ui-*.R): the app in headless Chrome, driven by shinytest2.
# The app runs in a few background R processes, one per set of launch options
# (at most three), each started on first use and stopped after the last test
# file. Every test opens its own browser tab on one of them, which is a new
# Shiny session: the R process, its packages and its compiled code are
# shared, the session state is not. inst/app is the source tree under
# devtools::test() and the installed app under R CMD check.
# R CMD check on CI skips the UI tests: on the GitHub runners a dataset load
# sometimes never finishes, for reasons not yet known. The UI test workflow
# (.github/workflows/ui-tests.yaml) runs them with EXPLOREDE_UI_TESTS=true,
# and with EXPLOREDE_UI_ARTIFACTS set, which keeps there each app process's
# R log, profile and process snapshot, and the browser log, a screenshot and
# the shown tab of every test.
# -----------------------------------------------------------------------------

.ui_apps <- new.env(parent = emptyenv())

# A new session (browser tab) on the app process started with `options`
# (exploreDE.file, as exploreDE(file) sets it); closed when `.local_envir`
# (the calling function) ends.
app_session <- function(name, options = list(), .local_envir = parent.frame()) {
  if (!identical(Sys.getenv("EXPLOREDE_UI_TESTS"), "true")) {
    skip_on_ci()
  }
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("callr")
  app <- shinytest2::AppDriver$new(app_process_url(options), name = name, load_timeout = 180e3, timeout = 60e3)
  withr::defer(
    {
      app_save_state(app, name)
      app$stop()
    },
    envir = .local_envir
  )
  app
}

# The address of the app process started with `options`: started on first
# use, again if it has died, and stopped after the last test file.
app_process_url <- function(options) {
  key <- paste(deparse(options), collapse = "")
  running <- .ui_apps[[key]]
  if (!is.null(running) && running$process$is_alive()) {
    return(running$url)
  }
  # Test mode lets shinytest2 read outputs and inputs, as it sets when it
  # starts an app itself.
  options$shiny.testmode <- TRUE
  # Kept for inspection: every message between app and browser, including
  # when each output starts and finishes recalculating.
  if (nzchar(Sys.getenv("EXPLOREDE_UI_ARTIFACTS"))) {
    options$shiny.trace <- TRUE
  }
  log <- app_log_file()
  app_dir <- system.file("app", package = "exploreDE")
  process <- callr::r_bg(
    function(app_dir, options) {
      options(options)
      # The app directory's .Rprofile, as when shinytest2 starts an app there
      # (ui-tests.yaml puts one there that profiles the app)
      profile <- file.path(app_dir, ".Rprofile")
      if (file.exists(profile)) {
        sys.source(profile, envir = globalenv())
      }
      shiny::runApp(app_dir, launch.browser = FALSE)
    },
    args = list(app_dir = app_dir, options = options),
    stdout = log,
    stderr = "2>&1"
  )
  url <- app_wait_for_url(process, log)
  .ui_apps[[key]] <- list(process = process, url = url, log = log)
  withr::defer(app_process_stop(key), envir = teardown_env())
  url
}

# Where an app process writes its output: the artifacts directory, or a
# temporary file.
app_log_file <- function() {
  dir <- Sys.getenv("EXPLOREDE_UI_ARTIFACTS")
  if (!nzchar(dir)) {
    return(tempfile("app-", fileext = ".log"))
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("app-%d-R.log", length(list.files(dir, "^app-[0-9]+-R\\.log$")) + 1L))
}

# The address shiny::runApp() reports once it listens.
app_wait_for_url <- function(process, log) {
  deadline <- Sys.time() + 180
  repeat {
    lines <- if (file.exists(log)) readLines(log, warn = FALSE) else character(0)
    hit <- regmatches(lines, regexpr("http://127\\.0\\.0\\.1:[0-9]+", lines))
    if (length(hit)) {
      return(hit[[1]])
    }
    if (!process$is_alive() || Sys.time() > deadline) {
      stop("The app did not start:\n", paste(lines, collapse = "\n"), call. = FALSE)
    }
    Sys.sleep(0.5)
  }
}

# Stop an app process: interrupted, so that runApp() returns and R writes its
# profile, then killed if it is still running. With EXPLOREDE_UI_ARTIFACTS
# set, a snapshot of the process is kept first.
app_process_stop <- function(key) {
  running <- .ui_apps[[key]]
  if (is.null(running) || !running$process$is_alive()) {
    return(invisible())
  }
  if (nzchar(Sys.getenv("EXPLOREDE_UI_ARTIFACTS"))) {
    pid <- running$process$get_pid()
    writeLines(c(paste("pid", pid), app_proc_snapshot(pid)), sub("-R\\.log$", "-proc.txt", running$log))
  }
  running$process$interrupt()
  running$process$wait(5000)
  running$process$kill()
  invisible()
}

# With EXPLOREDE_UI_ARTIFACTS set, write the browser's log, the shown tab and
# a screenshot of a session to that directory, as <name>-log.txt,
# <name>-state.json and <name>.png. The app's R log is the process's
# app-<n>-R.log.
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

# A two-contrast RNA-seq .h5ad named <name>.h5ad in `dir`.
toy_h5ad <- function(dir, name, seed) {
  dds <- make_toy_dds(n_genes = 80, seed = seed)
  out <- file.path(dir, paste0(name, ".h5ad"))
  build_explore_h5ad(
    dds,
    out = out,
    omics = "rnaseq",
    de = list(
      control_vs_treated = DESeq2::results(dds, contrast = c("condition", "control", "treated")),
      treated_vs_control = DESeq2::results(dds, contrast = c("condition", "treated", "control"))
    ),
    organism = "Homo sapiens"
  )
  out
}

# The Multi DE datasets, as listed on the Data Selection tab.
multi_datasets <- function(app) {
  unlist(app$get_value(input = "multiContrastOrder")) %||% character(0)
}
