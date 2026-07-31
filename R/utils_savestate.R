# =============================================================================
# Save/restore session-state — pure helpers (unit-testable; no Shiny session)
# Used by inst/app/server-saveState.R and the restore_default() closure in app.R.
# Kept in R/ (exported) so they are callable bare from the app and covered by the
# testthat suite, mirroring new_explore_dataset() etc.
# =============================================================================

#' Filter a Shiny input snapshot down to the restorable widgets
#'
#' Drops file uploads (stale datapath), action buttons, DT/plotly UI echoes, the
#' save-state tab's own widgets, download buttons, and every `multi*` input
#' (out of scope for save-state v1).
#'
#' @param x Named list, typically `reactiveValuesToList(input, all.names = TRUE)`.
#' @return `x` with the non-restorable entries removed.
#' @export
save_state_filter_inputs <- function(x) {
  ids <- names(x)
  skip_pat <- paste0(
    "^(load_state_file|importColourFile|multiImportFile|multiImportColourFile",
    "|fb_reset_|fb_add_|resetGeneBucket|resetDEGTableSelectRows|multi)",
    "|(_btn|_cell_edit|_cell_clicked|_rows_selected|_row_last_clicked",
    "|_state|_search|_hover|_click|_relayout|_brush|_selected)$")
  keep <- !grepl(skip_pat, ids)
  nv_exact  <- c("download_state", "sidebarMenu", "load_settings_file")
  nv_prefix <- grepl("^(dl|download)", ids)
  nv_suffix <- grepl("(Go|Button|Brush)$", ids)
  keep <- keep & !(ids %in% nv_exact) & !nv_prefix & !nv_suffix
  x[keep]
}

#' Build a dataset-identity fingerprint for save-state
#'
#' Identity = resolved `dataDir` plus a content fingerprint built entirely from
#' base-class `ExploreDataset` accessors (uniform across RNASeq / CountQC /
#' CrisprCountQC / proteomics; `factor_names`/`contrasts` default to empty where
#' a concept doesn't apply). Every field is `tryCatch`-wrapped so no data form can
#' make it throw, and collections are sorted so ordering can't cause a false
#' mismatch.
#'
#' @param eds An `ExploreDataset` (or `NULL`).
#' @param dataDir Resolved dataset path/URL, or `NULL`.
#' @return A named list of identity fields, or `NULL` if `eds` is `NULL`.
#' @export
save_state_fingerprint <- function(eds, dataDir = NULL) {
  if (is.null(eds)) return(NULL)
  list(
    dataDir   = if (is.null(dataDir)) NA_character_ else as.character(dataDir),
    dataType  = tryCatch(as.character(eds$data_type()),         error = function(e) NA_character_),
    omics     = tryCatch(as.character(eds$omics_class()),       error = function(e) NA_character_),
    nSamples  = tryCatch(as.integer(eds$n_samples()),           error = function(e) NA_integer_),
    nFeatures = tryCatch(as.integer(length(eds$feature_ids())), error = function(e) NA_integer_),
    factors   = tryCatch(sort(as.character(eds$factor_names())), error = function(e) character()),
    contrasts = tryCatch(sort(as.character(eds$contrasts())),    error = function(e) character())
  )
}

#' Do two fingerprints denote the exact same dataset?
#'
#' `dataDir` must match when both sides have it; the content fields must always
#' match (this also covers the NA-dataDir teaching/local case, and catches a
#' dataset regenerated in place).
#'
#' @param a,b Fingerprints from [save_state_fingerprint()] (or `NULL`).
#' @return `TRUE` only if the two denote the same dataset.
#' @export
save_state_fp_match <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  da <- a$dataDir; db <- b$dataDir
  if (!is.null(da) && !is.null(db) && !is.na(da) && !is.na(db) && !identical(da, db)) return(FALSE)
  identical(a$dataType,  b$dataType)  &&
  identical(a$omics,     b$omics)     &&
  identical(a$nSamples,  b$nSamples)  &&
  identical(a$nFeatures, b$nFeatures) &&
  identical(a$factors,   b$factors)   &&
  identical(a$contrasts, b$contrasts)
}

#' Resolve a widget default during a session-state restore
#'
#' During a restore (`restoring()` TRUE) returns the saved value for `id` from
#' `restore_seed` (optionally intersected with `valid` = the widget's current
#' choices); otherwise the computed default. A no-op outside restore, so live
#' behaviour is unchanged. The `restore_default()` closure in app.R delegates here
#' so the logic is unit-testable.
#'
#' @param id Input id (character).
#' @param computed_default Value to use when not restoring / nothing saved.
#' @param valid Optional vector of currently-valid choices to intersect with.
#' @param restoring A `reactiveVal` returning `TRUE` during a restore.
#' @param restore_seed A `reactiveValues` holding the saved input snapshot.
#' @return The saved value while restoring (if valid), else `computed_default`.
#' @export
restore_default_value <- function(id, computed_default, valid = NULL, restoring, restore_seed) {
  if (!isTRUE(shiny::isolate(restoring()))) return(computed_default)
  v <- shiny::isolate(restore_seed[[id]])
  if (is.null(v)) return(computed_default)
  if (!is.null(valid)) {
    v <- v[v %in% valid]
    if (!length(v)) return(computed_default)
  }
  v
}
