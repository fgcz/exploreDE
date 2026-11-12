# =============================================================================
# Save/restore session-state — pure helpers (unit-testable; no Shiny session)
# Used by inst/app/server-saveState.R and the restore_default() closure in app.R.
# Kept in R/ (exported) so they are callable bare from the app and covered by the
# testthat suite, mirroring new_explore_dataset() etc.
# =============================================================================

#' Filter a Shiny input snapshot down to the restorable widgets
#'
#' Drops file uploads (stale datapath), action buttons, DT/plotly UI echoes, the
#' save-state tab's own widgets, download buttons, every `multi*` input
#' (out of scope for save-state v1), and the per-tab `contrastSelected_<tab>`
#' mirrors (only the canonical `contrastSelected` is saved; restoring it syncs
#' the mirrors), and the sample-ratio widgets (ratios are not saved).
#'
#' @param x Named list, typically `reactiveValuesToList(input, all.names = TRUE)`.
#' @return `x` with the non-restorable entries removed.
#' @noRd
save_state_filter_inputs <- function(x) {
  ids <- names(x)
  skip_pat <- paste0(
    "^(load_state_file|uploadH5ad|importColourFile|multiImportFile|multiImportColourFile",
    "|fb_reset_|fb_add_|resetGeneBucket|resetDEGTableSelectRows|multi",
    "|contrastSelected_)",
    "|(_btn|_cell_edit|_cell_clicked|_rows_selected|_row_last_clicked",
    "|_state|_search|_hover|_click|_relayout|_brush|_selected)$"
  )
  keep <- !grepl(skip_pat, ids)
  # Sample ratios are not saved, so neither are the widgets that name them
  ratio_widgets <- paste0("heatmapRatio", c("Plot", "Name", "Numerator", "Denominator", "Remove"))
  nv_exact <- c("download_state", "sidebarMenu", "load_settings_file", ratio_widgets)
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
#' a concept doesn't apply). Collections are sorted so ordering can't cause a
#' false mismatch.
#'
#' @param eds An `ExploreDataset` (or `NULL`).
#' @param dataDir Resolved dataset path/URL, or `NULL`.
#' @return A named list of identity fields, or `NULL` if `eds` is `NULL`.
#' @noRd
save_state_fingerprint <- function(eds, dataDir = NULL) {
  if (is.null(eds)) {
    return(NULL)
  }
  list(
    dataDir = if (is.null(dataDir)) NA_character_ else as.character(dataDir),
    dataType = as.character(eds$data_type()),
    omics = as.character(eds$omics_class()),
    nSamples = as.integer(eds$n_samples()),
    nFeatures = length(eds$feature_ids()),
    factors = sort(as.character(eds$factor_names())),
    contrasts = sort(as.character(eds$contrasts()))
  )
}

# nolint start: cyclocomp_linter. Field-by-field comparison.
#' Do two fingerprints denote the exact same dataset?
#'
#' `dataDir` must match when both sides have it; the content fields must always
#' match (this also covers the NA-dataDir teaching/local case, and catches a
#' dataset regenerated in place).
#'
#' @param a,b Fingerprints from `save_state_fingerprint()` (or `NULL`).
#' @return `TRUE` only if the two denote the same dataset.
#' @noRd
save_state_fp_match <- function(a, b) {
  if (is.null(a) || is.null(b)) {
    return(FALSE)
  }
  da <- a$dataDir
  db <- b$dataDir
  if (!is.null(da) && !is.null(db) && !is.na(da) && !is.na(db) && !identical(da, db)) {
    return(FALSE)
  }
  identical(a$dataType, b$dataType) &&
    identical(a$omics, b$omics) &&
    identical(a$nSamples, b$nSamples) &&
    identical(a$nFeatures, b$nFeatures) &&
    identical(a$factors, b$factors) &&
    identical(a$contrasts, b$contrasts)
}
# nolint end

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
#' @noRd
restore_default_value <- function(id, computed_default, valid = NULL, restoring, restore_seed) {
  if (!isTRUE(shiny::isolate(restoring()))) {
    return(computed_default)
  }
  v <- shiny::isolate(restore_seed[[id]])
  if (is.null(v)) {
    return(computed_default)
  }
  if (!is.null(valid)) {
    v <- v[v %in% valid]
    if (!length(v)) return(computed_default)
  }
  v
}

# ---- Multi DE --------------------------------------------------------------

#' Pick the restorable Multi DE widgets from a Shiny input snapshot
#'
#' Keeps the `multi*` inputs `save_state_filter_inputs()` leaves out, minus the
#' ones a restore must not replay: File uploads, the project browser's
#' checkboxes, the order bucket (the dataset order is restored directly),
#' table-selection and download echoes, and action buttons.
#'
#' @param x Named list, typically `reactiveValuesToList(input, all.names = TRUE)`.
#' @return The restorable `multi*` entries of `x`.
#' @noRd
save_state_multi_inputs <- function(x) {
  ids <- names(x)
  keep <- grepl("^multi", ids) &
    !ids %in%
      c("multiImportFile", "multiImportColourFile", "multiUploadH5ad", "multiDatasetsCheckbox", "multiContrastOrder") &
    !grepl("^multiDl", ids) &
    !grepl(
      paste0(
        "(_rows_selected|_rows_current|_rows_all|_row_last_clicked|_state|_search|_cell_clicked|_cell_edit|",
        "_cells_selected)$"
      ),
      ids
    )
  x <- x[keep]
  Filter(function(v) !inherits(v, "shinyActionButtonValue"), x)
}

# nolint start: cyclocomp_linter. One branch per record kind.
#' Describe the loaded Multi DE datasets for a saved session
#'
#' One entry per dataset, in display order: Its display name and where it came
#' from, a gStore link (`source = "link"`, root-relative), a contrast of the
#' loaded dataset (`source = "contrast"`) or a contrast of an uploaded file
#' (`source = "upload"`, with the file's label). Datasets without a recorded
#' source are left out.
#'
#' @param results `multiDeReactive$results`: A named list of Multi DE objects,
#'   each carrying `$source` (`list(kind, link)`, `list(kind, contrast)` or
#'   `list(kind, file, contrast)`).
#' @return A list of `list(name, source, link)` / `list(name, source, contrast)`
#'   / `list(name, source, file, contrast)`.
#' @noRd
save_state_multi_datasets <- function(results) {
  out <- list()
  for (nm in names(results)) {
    src <- results[[nm]][["source"]]
    kind <- src[["kind"]]
    if (identical(kind, "link") && is.character(src[["link"]])) {
      out[[length(out) + 1L]] <- list(name = nm, source = "link", link = src[["link"]])
    } else if (identical(kind, "contrast") && is.character(src[["contrast"]])) {
      out[[length(out) + 1L]] <- list(name = nm, source = "contrast", contrast = src[["contrast"]])
    } else if (identical(kind, "upload") && is.character(src[["file"]]) && is.character(src[["contrast"]])) {
      out[[length(out) + 1L]] <- list(name = nm, source = "upload", file = src[["file"]], contrast = src[["contrast"]])
    }
  }
  out
}
# nolint end

# ---- JSON session files (format exploreDE-session/2) -----------------------
# See inst/spec/exploreDE-session-2.md. JSON can only hold text, numbers,
# true/false, null, arrays, and objects, so reading a session file can never
# run code. (The .rds format it replaced could; it was retired in 3.9.0.)
# Every value is written as an array (no unboxing), so single values and
# multi-selects read back the same way.

SAVE_STATE_JSON_FORMAT <- "exploreDE-session/2"

#' Write a saved-session bundle as JSON
#'
#' @param bundle The session bundle: `list(meta, inputs, rv, multi)`.
#' @return The JSON text (character(1)).
#' @noRd
save_state_to_json <- function(bundle) {
  drop_buttons <- function(x) Filter(function(v) !inherits(v, "shinyActionButtonValue"), x)
  bundle$inputs <- drop_buttons(bundle$inputs %||% list())
  if (!is.null(bundle$multi$inputs)) {
    bundle$multi$inputs <- drop_buttons(bundle$multi$inputs)
  }
  bundle$meta$format <- SAVE_STATE_JSON_FORMAT
  as.character(jsonlite::toJSON(
    bundle,
    auto_unbox = FALSE,
    null = "null",
    na = "null",
    digits = NA,
    POSIXt = "ISO8601",
    pretty = TRUE
  ))
}

# nolint start: cyclocomp_linter. One guard per session-file field.
#' Read a JSON saved-session file
#'
#' Refuses anything larger than `max_bytes`, anything that isn't an exploreDE
#' session in format `exploreDE-session/2`, and anything nested deeper than a
#' session ever is. Returns the bundle in the shape the restore handler uses
#' (`inst/app/server-saveState.R`).
#'
#' @param path Path to the uploaded `.json` file.
#' @param max_bytes Size limit; 5 MB by default.
#' @return The bundle: `list(meta, inputs, rv, multi)`. Stops with a short
#'   message on any problem.
#' @noRd
save_state_from_json <- function(path, max_bytes = 5 * 1024^2) {
  size <- file.size(path)
  if (is.na(size) || size > max_bytes) {
    stop(sprintf("the file is larger than %d MB", max_bytes %/% 1024^2))
  }
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  x <- tryCatch(
    jsonlite::fromJSON(txt, simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE),
    error = function(e) stop("not valid JSON: ", conditionMessage(e))
  )
  depth <- function(v) if (is.list(v) && length(v)) 1L + max(vapply(v, depth, integer(1))) else 0L
  if (!is.list(x) || is.null(names(x)) || depth(x) > 6L) {
    stop("not an exploreDE session file")
  }
  meta <- x[["meta"]]
  if (!identical(meta[["app"]], "exploreDE")) {
    stop("not an exploreDE session file")
  }
  if (!identical(meta[["format"]], SAVE_STATE_JSON_FORMAT)) {
    stop(sprintf("unsupported session format (expected %s)", SAVE_STATE_JSON_FORMAT))
  }

  if (!is.null(meta[["fingerprint"]])) {
    meta$fingerprint <- .session_fingerprint_from_json(meta[["fingerprint"]])
  }
  meta$saved_at <- as.POSIXct(.json_one(meta[["saved_at"]]), format = "%Y-%m-%dT%H:%M:%S")
  list(
    meta = meta,
    inputs = .json_named_list(x[["inputs"]]),
    rv = .session_rv_from_json(x[["rv"]] %||% list()),
    multi = if (!is.null(x[["multi"]])) .session_multi_from_json(x[["multi"]])
  )
}
# nolint end

# ---- JSON read-back helpers ---------------------------------------------------

# [] reads back as list(); the restore expects character(0) there.
.json_empty_to_chr <- function(v) {
  if (is.list(v) && !length(v) && is.null(names(v))) {
    return(character())
  }
  if (is.list(v)) {
    v[] <- lapply(v, .json_empty_to_chr)
  }
  v
}
.json_chr <- function(v) as.character(unlist(v))
.json_one <- function(v, as = as.character) if (length(v)) as(unlist(v))[1] else as(NA)
# A named list of widget values, or list() for anything else
.json_named_list <- function(v) {
  if (is.list(v) && length(v) && !is.null(names(v))) .json_empty_to_chr(v) else list()
}

.session_fingerprint_from_json <- function(fp) {
  list(
    dataDir = .json_one(fp[["dataDir"]]),
    dataType = .json_one(fp[["dataType"]]),
    omics = .json_one(fp[["omics"]]),
    nSamples = .json_one(fp[["nSamples"]], as.integer),
    nFeatures = .json_one(fp[["nFeatures"]], as.integer),
    factors = .json_chr(fp[["factors"]]),
    contrasts = .json_chr(fp[["contrasts"]])
  )
}

.session_rv_from_json <- function(rv) {
  if (!is.null(rv[["featureBucket"]])) {
    rv$featureBucket <- .json_chr(rv[["featureBucket"]])
  }
  for (k in c("oraSelected", "gseaSelected")) {
    if (!is.null(rv[[k]])) rv[[k]] <- as.integer(unlist(rv[[k]]))
  }
  rv
}

# Multi DE datasets: each a link, a contrast of the loaded dataset, or an
# uploaded file's contrast; an entry missing its references is dropped.
.session_multi_from_json <- function(multi) {
  ds <- lapply(multi[["datasets"]] %||% list(), function(d) {
    kind <- .json_one(d[["source"]])
    name <- .json_one(d[["name"]])
    if (is.na(name) || !kind %in% c("link", "contrast", "upload")) {
      return(NULL)
    }
    refs <- switch(kind, link = "link", contrast = "contrast", upload = c("file", "contrast"))
    out <- list(name = name, source = kind)
    for (r in refs) {
      v <- .json_one(d[[r]])
      if (is.na(v)) {
        return(NULL)
      }
      out[[r]] <- v
    }
    out
  })
  list(datasets = Filter(Negate(is.null), ds), inputs = .json_named_list(multi[["inputs"]]))
}
