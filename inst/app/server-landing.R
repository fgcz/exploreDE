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
    html  = tagList(waiter::spin_fading_circles(),
                    tags$h4(msg, style = "color:#fff; margin-top:1em;")),
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

# Upload a canonical .h5ad. Shiny stores the file at a temp datapath and cleans
# it up itself — nothing is written to the server beyond that ephemeral copy.
observeEvent(input$uploadH5ad, {
  f <- input$uploadH5ad
  req(f)
  .landing_spinner("Loading your data…")
  res <- load_h5ad_upload(f$datapath)
  if (isTRUE(res$ok)) {
    # Success: edsState is set → the deferred observer sources modules and hides
    # the spinner. Just surface a status message.
    output$upload_h5ad_status <- renderText(
      "Loaded successfully — your data is ready. Explore the tabs on the left."
    )
    showNotification("Data loaded successfully.", type = "message", duration = 6)
  } else {
    waiter::waiter_hide()
    output$upload_h5ad_status <- renderText(paste0(
      "Could not load this file.\n\n", res$msg,
      "\n\nMake sure it was built with build_explore_h5ad()."
    ))
    showNotification("Upload failed — see the message on the landing screen.",
                     type = "error", duration = 8)
  }
})
