# ═══════════════════════════════════════════════════════════════════════════
# SAVE / RESTORE SESSION STATE  (v1 — single-dataset; multi-DE excluded)
# ---------------------------------------------------------------------------
# Serialise the app's widget state + a few reactiveValues that live outside
# input$ to a single .rds the user downloads, and restore it on re-upload.
# See plans/archive/save_state_plan.md for the full design. Restore uses three paths:
#   (1) a sendInputMessage retry burst for standard widgets,
#   (2) the restore_seed store (consulted by restore_default() in every rebuilt
#       renderUI / child observer) for widgets that can't be message-set,
#   (3) direct reactiveValues injection for the feature bucket + GO/KEGG genes,
#       and DT::selectRows-by-Description for the ORA/GSEA pathway selections.
# The dataset itself is NOT saved — it comes from the ?data= URL; a bundle is
# only valid against the same dataset.
# ═══════════════════════════════════════════════════════════════════════════

# The pure helpers live in R/utils_savestate.R (exported, unit-tested):
#   save_state_filter_inputs(x), save_state_fingerprint(eds, dataDir),
#   save_state_fp_match(a, b), restore_default_value(...).
# current_fingerprint() wraps save_state_fingerprint() with this session's eds.
current_fingerprint <- function() {
  idr <- tryCatch(inputDataReactive(), error = function(e) NULL)
  save_state_fingerprint(
    tryCatch(idr$.eds,             error = function(e) NULL),
    tryCatch(attr(idr, "dataDir"), error = function(e) NULL))
}

# ── SAVE ────────────────────────────────────────────────────────────────────
output$download_state <- downloadHandler(
  filename = function() {
    dt <- tryCatch(inputDataReactive()$.eds$data_type(), error = function(e) "data")
    sprintf("exploreDE_state_%s_%s.rds", dt, format(Sys.time(), "%Y%m%d_%H%M%S"))
  },
  content = function(file) {
    eds <- tryCatch(inputDataReactive()$.eds, error = function(e) NULL)

    # ORA / GSEA selected pathway rows. The enrichResult behind these tables is
    # only materialised after a selection fires the plot observer, so it isn't
    # available to map Descriptions on restore — but for the same dataset the
    # table renders in identical row order, so the saved row indices re-apply
    # exactly (re-selection is keyed on oraType/oraDirection/gseaType, which are
    # restored too). Store indices; restore re-selects them on tab open.
    ora_sel  <- tryCatch(input$selectedTable_ORA_rows_selected,  error = function(e) NULL)
    gsea_sel <- tryCatch(input$selectedTable_GSEA_rows_selected, error = function(e) NULL)

    bundle <- list(
      meta = list(
        app         = "exploreDE",
        version     = 1L,
        saved_at    = Sys.time(),
        dataUrl     = tryCatch(isolate(session$clientData$url_search), error = function(e) NA_character_),
        fingerprint = current_fingerprint()   # exact-dataset identity (see restore)
      ),
      inputs = save_state_filter_inputs(reactiveValuesToList(input, all.names = TRUE)),
      rv = list(
        featureBucket = tryCatch(featureBucket$features, error = function(e) NULL),
        goHeatmap = tryCatch(list(
          genes = goHeatmapReactive$currentGenes,
          name  = goHeatmapReactive$currentGoName,
          term  = goHeatmapReactive$currentGoTerm), error = function(e) NULL),
        keggHeatmap = tryCatch(list(
          genes = keggHeatmapReactive$currentGenes,
          name  = keggHeatmapReactive$currentPathwayName), error = function(e) NULL),
        goTileKegg = tryCatch(list(
          genes  = goTileKeggReactive$currentGenes,
          keggId = goTileKeggReactive$currentKeggId,
          name   = goTileKeggReactive$currentPathwayName), error = function(e) NULL),
        oraSelected  = ora_sel,
        gseaSelected = gsea_sel,
        gseaType     = tryCatch(isolate(input$gseaType), error = function(e) NULL)
      )
    )
    saveRDS(bundle, file)
  }
)

output$save_state_status <- renderText({
  n <- tryCatch(length(save_state_filter_inputs(reactiveValuesToList(input, all.names = TRUE))),
                error = function(e) 0L)
  if (isTRUE(tryCatch(!is.null(inputDataReactive()), error = function(e) FALSE))) {
    sprintf("Ready to save — %d settings will be captured.", n)
  } else {
    "Waiting for data to load…"
  }
})

# ── RESTORE machinery ───────────────────────────────────────────────────────
# Retry burst: sendInputMessage each queued value until input$<id> matches it.
# Widgets on not-yet-opened tabs mount when the tab opens, so the burst is
# re-armed on every input$tabs change (below). While the burst is live the
# `restoring` flag is TRUE, which suppresses the M2 pure-clear observers.
pending_inputs   <- reactiveVal(list())
restore_ticks    <- reactiveVal(0L)
restore_deadline <- reactiveVal(NULL)   # bounds the whole restore window

observe({
  pend  <- pending_inputs()
  ticks <- restore_ticks()
  dl    <- isolate(restore_deadline())
  expired <- !is.null(dl) && Sys.time() > dl
  # Restore over (queue drained) or window elapsed — leave restore mode; drop any
  # stragglers that never applied so normal use isn't impaired long-term.
  if (length(pend) == 0 || expired) {
    if (isolate(restoring())) restoring(FALSE)
    if (expired && length(pend) > 0) pending_inputs(list())
    return()
  }
  if (ticks >= 40L) { if (isolate(restoring())) restoring(FALSE); return() }  # burst spent
  if (!isolate(restoring())) restoring(TRUE)
  invalidateLater(250, session)
  still <- list()
  for (id in names(pend)) {
    v <- pend[[id]]
    if (is.null(v)) next
    cur <- isolate(input[[id]])
    if (!is.null(cur) && identical(cur, v)) next  # already applied
    tryCatch(session$sendInputMessage(id, list(value = v)), error = function(e) NULL)
    still[[id]] <- v
  }
  if (length(still) != length(pend)) pending_inputs(still)
  restore_ticks(ticks + 1L)
})

# Opening a tab mounts its widgets — re-arm the burst so their saved values land,
# but only within the restore window (after that, stragglers are abandoned).
observeEvent(input$tabs, {
  dl <- restore_deadline()
  if (length(pending_inputs()) > 0 && !is.null(dl) && Sys.time() <= dl) restore_ticks(0L)
}, ignoreInit = TRUE)

# ORA / GSEA DT re-selection: re-apply the saved row indices once the user opens
# the tab (the table renders then). shinyjs::delay lets the DataTable mount before
# selectRows fires; each is consumed once so later manual clicks aren't disturbed.
save_state_pending_ora  <- reactiveVal(NULL)
save_state_pending_gsea <- reactiveVal(NULL)

observeEvent(input$tabs, {
  if (identical(input$tabs, "oraTab") && length(save_state_pending_ora()) > 0) {
    idx <- save_state_pending_ora(); save_state_pending_ora(NULL)
    shinyjs::delay(400, tryCatch(
      DT::selectRows(DT::dataTableProxy("selectedTable_ORA"), idx), error = function(e) NULL))
  }
  if (identical(input$tabs, "gseaTab") && length(save_state_pending_gsea()) > 0) {
    idx <- save_state_pending_gsea(); save_state_pending_gsea(NULL)
    shinyjs::delay(400, tryCatch(
      DT::selectRows(DT::dataTableProxy("selectedTable_GSEA"), idx), error = function(e) NULL))
  }
}, ignoreInit = TRUE)

# ── RESTORE handler ─────────────────────────────────────────────────────────
observeEvent(input$load_state_file, {
  f <- input$load_state_file
  req(f)
  bundle <- tryCatch(readRDS(f$datapath),
                     error = function(e) list(.error = conditionMessage(e)))
  if (!is.null(bundle$.error) || !is.list(bundle) ||
      !identical(bundle$meta$app, "exploreDE")) {
    msg <- paste0("✗ Not a valid exploreDE session file",
                  if (!is.null(bundle$.error)) paste0(" (", bundle$.error, ")") else "", ".")
    output$load_state_status <- renderText(msg)
    showNotification(msg, type = "error", duration = 8)
    return()
  }

  # Exact-dataset check — a session may ONLY be restored on the dataset it was
  # saved from (see save_state_fingerprint). Refuse cleanly on any mismatch rather
  # than half-applying stale factors/contrasts/genes into a foreign dataset.
  cur_fp   <- current_fingerprint()
  saved_fp <- bundle$meta$fingerprint
  if (!isTRUE(save_state_fp_match(saved_fp, cur_fp))) {
    sdt <- if (!is.null(saved_fp) && !is.null(saved_fp$dataType) && !is.na(saved_fp$dataType)) saved_fp$dataType else "unknown"
    cdt <- if (!is.null(cur_fp)   && !is.null(cur_fp$dataType)   && !is.na(cur_fp$dataType))   cur_fp$dataType   else "unknown"
    msg <- sprintf(paste0(
      "✗ This session was saved from a different dataset (%s) and can't be restored ",
      "on the one loaded here (%s). Open the app with the original dataset link, then restore."),
      sdt, cdt)
    output$load_state_status <- renderText(msg)
    showNotification(msg, type = "error", duration = 12)
    return()
  }

  inputs_bundle <- if (is.list(bundle$inputs)) bundle$inputs else list()

  # Enter restore mode (drives the M2 guards; the burst observe keeps it live).
  restoring(TRUE)
  restore_deadline(Sys.time() + 120)  # restore window: re-arm on tab opens for 2 min

  # (1) Seed store — consulted by restore_default() in every rebuilt widget.
  for (id in names(inputs_bundle)) restore_seed[[id]] <- inputs_bundle[[id]]

  # (2) reactiveValues injection — feature bucket + gated gene sets. The GO/KEGG
  #     heatmap renders READ these rv (they don't re-fetch), so plots rebuild.
  #     exists() guards: goTileKeggReactive only exists for RNASeq datasets, so a
  #     bundle carrying it must not crash a restore on another data type.
  if (!is.null(bundle$rv$featureBucket) && exists("featureBucket")) {
    featureBucket$features <- bundle$rv$featureBucket
    featureBucket$version  <- isolate(featureBucket$version) + 1L
  }
  if (!is.null(bundle$rv$goHeatmap) && exists("goHeatmapReactive")) {
    goHeatmapReactive$currentGenes  <- bundle$rv$goHeatmap$genes
    goHeatmapReactive$currentGoName <- bundle$rv$goHeatmap$name
    goHeatmapReactive$currentGoTerm <- bundle$rv$goHeatmap$term
  }
  if (!is.null(bundle$rv$keggHeatmap) && exists("keggHeatmapReactive")) {
    keggHeatmapReactive$currentGenes       <- bundle$rv$keggHeatmap$genes
    keggHeatmapReactive$currentPathwayName <- bundle$rv$keggHeatmap$name
  }
  if (!is.null(bundle$rv$goTileKegg) && exists("goTileKeggReactive")) {
    goTileKeggReactive$currentGenes       <- bundle$rv$goTileKegg$genes
    goTileKeggReactive$currentKeggId      <- bundle$rv$goTileKegg$keggId
    goTileKeggReactive$currentPathwayName <- bundle$rv$goTileKegg$name
  }

  # (3) Retry burst for standard widgets (active tab handled separately, last).
  queue <- inputs_bundle
  queue[["tabs"]] <- NULL
  pending_inputs(queue)
  restore_ticks(0L)

  # (4) Queue ORA/GSEA pathway re-selection (applied when their tables render).
  save_state_pending_ora(bundle$rv$oraSelected)
  save_state_pending_gsea(bundle$rv$gseaSelected)

  # (5) Active tab last — only if valid for this data type and not a multi tab.
  saved_tab <- inputs_bundle[["tabs"]]
  if (!is.null(saved_tab) && !grepl("[Mm]ulti", saved_tab)) {
    session$onFlushed(function() {
      tryCatch(updateTabItems(session, "tabs", saved_tab), error = function(e) NULL)
    }, once = TRUE)
  }

  saved_at <- tryCatch(format(bundle$meta$saved_at, "%Y-%m-%d %H:%M:%S"),
                       error = function(e) "unknown")
  msg <- sprintf("✓ Restored from '%s' (saved %s). %d settings queued.",
                 f$name, saved_at, length(inputs_bundle))
  output$load_state_status <- renderText(msg)
  showNotification(msg, type = "message", duration = 8)
})
