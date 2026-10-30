#' Runtime validator for a live `ExploreDataset` instance
#'
#' Checks subclass-agnostic invariants by walking the accessor API:
#'
#' 1. The runtime class is a concrete `ExploreDataset` subclass — not
#'    the abstract base or intermediate.
#' 2. Cross-references resolve: `factor_names()` are columns of
#'    `samples()`; every `contrasts()` entry has a working `de()`;
#'    `primary_grouping()` is a member of `factor_names()`;
#'    `primary_contrast()` is in `contrasts()` (or `NA`).
#' 3. Sample / feature alignment: assay rownames / colnames agree
#'    with `sample_ids()` / `feature_ids()`.
#'
#' Cost is proportional to the number of contrasts (one `de(c)`
#' round-trip each) plus a small constant of other accessor calls.
#'
#' @param eds An `ExploreDataset` instance.
#' @param strict If `TRUE` (default), `stop()` on the first failure.
#'   If `FALSE`, returns `list(ok, errors)`.
#' @return Invisibly, `TRUE` on success (strict mode); or
#'   `list(ok = logical(1), errors = character)` (non-strict mode).
#' @export
validate_explore_dataset <- function(eds, strict = TRUE) {
  errors <- character()

  add <- function(msg) errors[[length(errors) + 1L]] <<- msg

  # 1. Concrete subclass
  if (!inherits(eds, "ExploreDataset")) {
    add("Object does not inherit from 'ExploreDataset'.")
  }
  abstract <- c("ExploreDataset", "ExploreDatasetSE")
  rt_class <- class(eds)[1]
  if (rt_class %in% abstract) {
    add(sprintf("Runtime class '%s' is abstract; expected a concrete subclass.", rt_class))
  }

  # 2. Cross-references
  fnames <- tryCatch(eds$factor_names(), error = function(e) {
    add(paste0("factor_names() failed: ", conditionMessage(e)))
    character()
  })
  samples_df <- tryCatch(eds$samples(), error = function(e) {
    add(paste0("samples() failed: ", conditionMessage(e)))
    data.frame()
  })
  missing_factors <- setdiff(fnames, colnames(samples_df))
  if (length(missing_factors)) {
    add(sprintf("factor_names not in samples() columns: %s", paste(missing_factors, collapse = ", ")))
  }

  pg <- tryCatch(eds$primary_grouping(), error = function(e) {
    add(paste0("primary_grouping() failed: ", conditionMessage(e)))
    NA_character_
  })
  if (!is.na(pg) && !pg %in% fnames) {
    add(sprintf("primary_grouping '%s' is not a member of factor_names.", pg))
  }

  cs <- tryCatch(eds$contrasts(), error = function(e) {
    add(paste0("contrasts() failed: ", conditionMessage(e)))
    character()
  })
  pc <- tryCatch(eds$primary_contrast(), error = function(e) {
    add(paste0("primary_contrast() failed: ", conditionMessage(e)))
    NA_character_
  })
  if (!is.na(pc) && !pc %in% cs) {
    add(sprintf("primary_contrast '%s' is not a member of contrasts().", pc))
  }
  for (c in cs) {
    de <- tryCatch(eds$de(c), error = function(e) {
      add(sprintf("de('%s') failed: %s", c, conditionMessage(e)))
      NULL
    })
    if (is.null(de)) {
      next
    }
    if (!is.data.frame(de)) {
      add(sprintf("de('%s') is not a data.frame.", c))
    }
  }

  # 3. Sample / feature alignment (against raw assay)
  raw_name <- tryCatch(eds$raw_assay_name(), error = function(e) {
    add(paste0("raw_assay_name() failed: ", conditionMessage(e)))
    NA_character_
  })
  if (!is.na(raw_name) && raw_name %in% (tryCatch(eds$assay_names(), error = function(e) character()))) {
    a <- tryCatch(eds$raw_assay(), error = function(e) {
      add(paste0("raw_assay() failed: ", conditionMessage(e)))
      NULL
    })
    if (!is.null(a)) {
      sids <- eds$sample_ids()
      fids <- eds$feature_ids()
      if (!is.null(colnames(a)) && !all(colnames(a) %in% sids)) {
        add("raw_assay() colnames are not a subset of sample_ids().")
      }
      if (!is.null(rownames(a)) && !all(rownames(a) %in% fids)) {
        add("raw_assay() rownames are not a subset of feature_ids().")
      }
    }
  }

  ok <- length(errors) == 0L
  if (strict && !ok) {
    stop(
      "validate_explore_dataset(): ",
      length(errors),
      " failure(s)\n  ",
      paste(errors, collapse = "\n  "),
      call. = FALSE
    )
  }
  if (strict) {
    return(invisible(TRUE))
  }
  list(ok = ok, errors = errors)
}
