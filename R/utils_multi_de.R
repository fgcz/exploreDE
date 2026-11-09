# =============================================================================
# Multi-DE dataset links — helpers (unit-testable; no Shiny session)
# Used by inst/app/server-multi-data-selection.R. A link is a gStore-relative
# result directory ("p1000/DESeq2_.../SDY--over--Control/"), an exploreDE URL
# carrying one (".../?data=p1000/..."), or an absolute path under one of the
# data roots (the project browser's checkbox values).
# =============================================================================

#' Resolve a multi-DE dataset link to a result directory under a data root
#'
#' The check is lexical, so it never follows the link anywhere first: A link
#' may not contain a `..` component, and an absolute path must lie under one of
#' `roots`. The first component below the root must be a project (`p` followed
#' by digits), and the directory must hold `deResult.rds` or `deResult.qs2`.
#'
#' @param link Character(1): A gStore-relative path, an exploreDE `?data=` URL,
#'   or an absolute path under one of `roots`.
#' @param roots Character: The data roots, e.g. `"/srv/gstore/projects"`.
#' @return `list(path, root, project, link)` (`link` is root-relative, the
#'   form a saved session stores), or `NULL` when the link is rejected or holds
#'   no DE result.
#' @noRd
multi_resolve_dataset <- function(link, roots) {
  if (!is.character(link) || length(link) != 1L || is.na(link)) {
    return(NULL)
  }
  link <- gsub("\\s+", "", link)
  if (grepl("^https?://", link)) {
    link <- sub(".*\\?data=", "", link)
  }
  if (!nzchar(link) || any(strsplit(link, "/", fixed = TRUE)[[1]] == "..")) {
    return(NULL)
  }
  for (root in sub("/+$", "", roots)) {
    if (startsWith(link, "/")) {
      if (!startsWith(link, paste0(root, "/"))) {
        next
      }
      rel <- substring(link, nchar(root) + 2L)
    } else {
      rel <- link
    }
    parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts) & parts != "."]
    if (!length(parts) || !grepl("^p[0-9]+$", parts[[1]])) {
      next
    }
    path <- paste(c(root, parts), collapse = "/")
    if (
      file.exists(file.path(path, "deResult.rds")) ||
        file.exists(file.path(path, "deResult.qs2"))
    ) {
      return(list(path = path, root = root, project = parts[[1]], link = paste(parts, collapse = "/")))
    }
  }
  NULL
}

#' Does an LDAP lookup allow access to a project?
#'
#' The same rule the `?data=` loader applies: Employees (role `R_2`) see every
#' project, users (`R_3`) only the projects they belong to, anyone else none.
#'
#' @param lookup `list(roles = character, projects = character)`, as returned
#'   by the app's `ldap_lookup()`.
#' @param project Character(1), e.g. `"p1000"`.
#' @return `TRUE` or `FALSE`.
#' @noRd
ldap_allows_project <- function(lookup, project) {
  roles <- lookup$roles %||% character()
  if ("R_2" %in% roles) {
    return(TRUE)
  }
  if ("R_3" %in% roles) {
    return(isTRUE(project %in% (lookup$projects %||% character())))
  }
  FALSE
}

#' Load one SUSHI DE result directory for the multi-DE tabs
#'
#' Reads `deResult.rds` (or `deResult.qs2`) and builds its Multi-DE record with
#' `legacy_multi_contrast_from_eds()`, the converter every dataset in the
#' Multi-DE tabs goes through. Only call it on a directory returned by
#' `multi_resolve_dataset()`.
#'
#' @param path A result directory holding `deResult.rds` or `deResult.qs2`.
#' @return The Multi-DE record, with `ora` / `gsea` when GO was run, or `NULL`
#'   when the directory holds no DE result.
#' @noRd
multi_load_de_result <- function(path) {
  if (file.exists(file.path(path, "deResult.rds"))) {
    se <- readRDS(file.path(path, "deResult.rds"))
  } else if (file.exists(file.path(path, "deResult.qs2"))) {
    se <- qs2::qs_read(file.path(path, "deResult.qs2"), nthreads = 4)
  } else {
    return(NULL)
  }
  eds <- new_explore_dataset(se)
  legacy_multi_contrast_from_eds(eds, eds$primary_contrast())
}

#' Label of a dataset added from an uploaded file
#'
#' The file name without its extension, reduced to letters, digits, `_` and
#' `-` (the label becomes part of Multi DE input ids) and made unique against
#' the labels already in use.
#'
#' @param filename The uploaded file's name.
#' @param existing Labels already in use.
#' @return Character(1).
#' @noRd
multi_upload_label <- function(filename, existing = character()) {
  label <- tools::file_path_sans_ext(basename(filename))
  label <- gsub("[^A-Za-z0-9_-]+", "_", label)
  label <- gsub("^_+|_+$", "", label)
  if (!nzchar(label)) {
    label <- "upload"
  }
  make.unique(c(existing, label), sep = "_")[length(existing) + 1L]
}

#' Multi DE records of every contrast of a dataset
#'
#' One record per contrast, built with `legacy_multi_contrast_from_eds()`.
#' With a `label` the records are named `<label>_<contrast>`, so contrasts of
#' several datasets can sit side by side; without one they carry the bare
#' contrast names. The source recorded for saved sessions is the upload
#' (`kind = "upload"`, with the file's label) when `uploaded` is `TRUE`, else
#' the loaded dataset's own contrast (`kind = "contrast"`). A contrast whose
#' record cannot be built is left out and named, with the error message, in
#' the `failed` attribute.
#'
#' @param eds An `ExploreDataset`.
#' @param label The dataset's label, or `NULL`.
#' @param uploaded Whether the dataset came from a file added for Multi DE.
#' @return A named list of Multi DE records (empty when the dataset has no
#'   contrasts), with attribute `failed`: a character vector of error messages
#'   named by contrast, empty when every record was built.
#' @noRd
multi_records_from_dataset <- function(eds, label = NULL, uploaded = !is.null(label)) {
  out <- list()
  failed <- character()
  for (con in eds$contrasts()) {
    obj <- tryCatch(legacy_multi_contrast_from_eds(eds, con), error = function(e) e)
    if (inherits(obj, "error")) {
      failed[[con]] <- conditionMessage(obj)
      next
    }
    obj$source <- if (uploaded) {
      list(kind = "upload", file = label, contrast = con)
    } else {
      list(kind = "contrast", contrast = con)
    }
    out[[if (is.null(label)) con else paste0(label, "_", con)]] <- obj
  }
  structure(out, failed = failed)
}
