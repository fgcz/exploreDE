# ---------------------------------------------------------------------------
# Landing screen (BYO data, M2)
#
# Two ways to populate the app when it was opened without a ?data= URL:
#   - upload a canonical .h5ad (built with build_explore_h5ad()), or
#   - load the demo dataset.
# load_h5ad_upload() / load_demo() live in server-initInputData.R and publish
# into the shared edsState reactive; once set, the deferred module-source
# observer (app.R) sources the tab servers and the tabs light up.
#
# The full-screen spinner is shown here the instant a load starts and is hidden
# by that deferred observer once the modules are sourced — a single continuous
# loading screen from click to ready. On a FAILED load (nothing published), we
# hide it here, because the deferred observer won't fire.
# ---------------------------------------------------------------------------

.landing_spinner <- function(msg) {
  waiter::waiter_show(
    html = tagList(waiter::spin_fading_circles(), tags$h4(msg, style = "color:#fff; margin-top:1em;")),
    color = "#527085"
  )
}

# Demo dataset button.
observeEvent(input$loadDemo, {
  .landing_spinner("Loading demo dataset…")
  load_demo()
  # If the load failed (no dataset published), the deferred observer never fires,
  # so unstick the spinner here.
  if (is.null(edsState())) waiter::waiter_hide()
})

# Upload canonical .h5ad files. Shiny stores them at temp datapaths and cleans
# them up itself — nothing is written to the server beyond those ephemeral
# copies. The first file that loads opens in the tabs; the contrasts of the
# others join it in Multi DE.
observeEvent(input$uploadH5ad, {
  f <- h5ad_uploads_only(input$uploadH5ad)
  req(nrow(f) > 0)
  .landing_spinner(if (nrow(f) > 1) "Loading your files…" else "Loading your data…")
  failed <- character()
  opened <- NULL
  for (i in seq_len(nrow(f))) {
    res <- load_h5ad_upload(f$datapath[i])
    if (isTRUE(res$ok)) {
      opened <- i
      break
    }
    failed <- c(failed, stats::setNames(res$msg, f$name[i]))
  }
  if (!is.null(opened)) {
    # Success: edsState is set → the deferred observer sources modules and hides
    # the spinner. Just surface a status message.
    msg <- "Loaded successfully — your data is ready. Explore the tabs on the left."
    rest <- f[-seq_len(opened), , drop = FALSE]
    if (nrow(f) > 1) {
      added <- list(contrasts = 0L, files = 0L)
      if (nrow(rest)) {
        added <- tryCatch(multi_add_uploads(rest), error = function(e) {
          showNotification(paste("Could not add the other files to Multi DE:", conditionMessage(e)), type = "error")
          added
        })
      }
      msg <- sprintf(
        "%s is open in the tabs; %d contrast%s from %d other file%s joined it in Multi DE Analysis.",
        f$name[opened],
        added$contrasts,
        if (added$contrasts == 1) "" else "s",
        added$files,
        if (added$files == 1) "" else "s"
      )
    }
    if (length(failed)) {
      msg <- paste0(
        msg,
        "\n\nCould not load ",
        paste(names(failed), collapse = ", "),
        ":\n",
        paste(failed, collapse = "\n")
      )
    }
    output$upload_h5ad_status <- renderText(msg)
    showNotification("Data loaded successfully.", type = "message", duration = 6)
  } else {
    waiter::waiter_hide()
    output$upload_h5ad_status <- renderText(paste0(
      if (length(failed) > 1) "Could not load these files.\n\n" else "Could not load this file.\n\n",
      paste(paste0(names(failed), ": ", failed), collapse = "\n\n"),
      "\n\nMake sure the files were built with build_explore_h5ad()."
    ))
    showNotification("Upload failed — see the message on the landing screen.", type = "error", duration = 8)
  }
})
