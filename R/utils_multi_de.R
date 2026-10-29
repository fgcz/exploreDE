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
#' @export
multi_resolve_dataset <- function(link, roots) {
  if (!is.character(link) || length(link) != 1L || is.na(link)) return(NULL)
  link <- gsub("\\s+", "", link)
  if (grepl("^https?://", link)) link <- sub(".*\\?data=", "", link)
  if (!nzchar(link) || any(strsplit(link, "/", fixed = TRUE)[[1]] == "..")) return(NULL)
  for (root in sub("/+$", "", roots)) {
    if (startsWith(link, "/")) {
      if (!startsWith(link, paste0(root, "/"))) next
      rel <- substring(link, nchar(root) + 2L)
    } else {
      rel <- link
    }
    parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts) & parts != "."]
    if (!length(parts) || !grepl("^p[0-9]+$", parts[[1]])) next
    path <- paste(c(root, parts), collapse = "/")
    if (file.exists(file.path(path, "deResult.rds")) ||
        file.exists(file.path(path, "deResult.qs2"))) {
      return(list(path = path, root = root, project = parts[[1]],
                  link = paste(parts, collapse = "/")))
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
#' @export
ldap_allows_project <- function(lookup, project) {
  roles <- lookup$roles %||% character()
  if ("R_2" %in% roles) return(TRUE)
  if ("R_3" %in% roles) return(isTRUE(project %in% (lookup$projects %||% character())))
  FALSE
}

#' Load one SUSHI DE result directory for the multi-DE tabs
#'
#' Reads `deResult.rds` (or `deResult.qs2`), converts it with
#' [convert_genomics_se()], makes sure a `log2Ratio` column exists, and adds a
#' `Label` column ("ID\\nDescription") to its ORA and GSEA results. Only call it
#' on a directory returned by [multi_resolve_dataset()].
#'
#' @param path A result directory holding `deResult.rds` or `deResult.qs2`.
#' @return The converted list, with `ora` / `gsea` when GO was run, or `NULL`
#'   when the directory holds no DE result.
#' @export
multi_load_de_result <- function(path) {
  if (file.exists(file.path(path, "deResult.rds"))) {
    se <- readRDS(file.path(path, "deResult.rds"))
  } else if (file.exists(file.path(path, "deResult.qs2"))) {
    se <- qs2::qs_read(file.path(path, "deResult.qs2"), nthreads = 4)
  } else {
    return(NULL)
  }
  formatted <- convert_genomics_se(se, path)

  # Compatibility fix: Ensure log2Ratio exists (for the multiDEG tabs)
  if (!"log2Ratio" %in% colnames(formatted$seqAnno)) {
    lfcCol <- grep("^(log2?[._ ]?(fold)?[._ ]?(change|ratio|fc)|logfc)$",
                   colnames(formatted$seqAnno), value = TRUE, ignore.case = TRUE)
    if (length(lfcCol) > 0) {
      formatted$seqAnno$log2Ratio <- formatted$seqAnno[[lfcCol[1]]]
    }
  }

  # "ID\nDescription" labels for the ORA / GSEA plots (as multiDEG did)
  add_labels <- function(er) {
    er@result$Label <- ""
    for (j in seq_along(er@result$Description)) {
      desc <- er@result$Description[j]
      if (nchar(desc) > 35) {
        desc <- substr(desc, 1, 45)
      }
      er@result$Label[j] <- paste0(er@result$ID[j], "\n", desc)
    }
    er
  }
  if (!is.null(formatted$param$runGO) && formatted$param$runGO) {
    ora <- S4Vectors::metadata(se)$enrichResult
    if (!is.null(ora)) {
      for (ont in names(ora)) {
        for (sig in names(ora[[ont]])) {
          if (nrow(ora[[ont]][[sig]]) >= 2) ora[[ont]][[sig]] <- add_labels(ora[[ont]][[sig]])
        }
      }
      formatted$ora <- ora
    }
    gsea <- S4Vectors::metadata(se)$enrichResultGSEA
    if (!is.null(gsea)) {
      for (ont in names(gsea)) {
        if (nrow(gsea[[ont]]) >= 2) gsea[[ont]] <- add_labels(gsea[[ont]])
      }
      formatted$gsea <- gsea
    }
  }
  formatted
}
