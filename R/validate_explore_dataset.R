#' Runtime validator for a live `ExploreDataset` instance
#'
#' Checks subclass-agnostic invariants by walking the accessor API:
#'
#' 1. The runtime class is a concrete `ExploreDataset` subclass — not
#'    the abstract base or intermediate.
#' 2. Cross-references resolve: `factor_names()` are columns of
#'    `samples()`; every `contrasts()` entry has a working `de()` with
#'    the `effect_log2fc` and `score_fdr` columns;
#'    `primary_grouping()` is a member of `factor_names()`;
#'    `primary_contrast()` is in `contrasts()` (or `NA`).
#' 3. Assays: `raw_assay_name()` and `primary_assay_name()` are in
#'    `assay_names()`, every assay has a log2 fact (`assay_is_log2()`), and the raw assay's rownames / colnames agree
#'    with `feature_ids()` / `sample_ids()`.
#'
#' Cost is proportional to the number of contrasts (one `de(c)`
#' round-trip each) plus a small constant of other accessor calls.
#'
#' @param eds An `ExploreDataset` instance.
#' @param strict If `TRUE` (default), `stop()` listing every failure.
#'   If `FALSE`, returns `list(ok, errors)` and never throws: a failing
#'   accessor is recorded as an error.
#' @return Invisibly, `TRUE` on success (strict mode); or
#'   `list(ok = logical(1), errors = character)` (non-strict mode).
#' @examples
#' se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
#' eds <- new_explore_dataset(se)
#' validate_explore_dataset(eds)
#' validate_explore_dataset(eds, strict = FALSE)
#' @export
validate_explore_dataset <- function(eds, strict = TRUE) {
  checked <- .vx_factors(eds)
  errors <- c(.vx_class(eds), checked$errors, .vx_contrasts(eds), .vx_assays(eds, checked$samples))
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

# Calls one accessor: its value, or `default` and the failure message.
.vx_call <- function(f, what, default) {
  tryCatch(
    list(value = f(), error = character()),
    error = function(e) list(value = default, error = paste0(what, " failed: ", conditionMessage(e)))
  )
}

# 1. A concrete subclass
.vx_class <- function(eds) {
  errors <- character()
  if (!inherits(eds, "ExploreDataset")) {
    errors <- c(errors, "Object does not inherit from 'ExploreDataset'.")
  }
  rt_class <- class(eds)[1]
  if (rt_class %in% c("ExploreDataset", "ExploreDatasetSE")) {
    errors <- c(errors, sprintf("Runtime class '%s' is abstract; expected a concrete subclass.", rt_class))
  }
  errors
}

# 2. Factors and the primary grouping are columns of samples(); returns the
# messages and samples() for the assay checks.
.vx_factors <- function(eds) {
  fnames <- .vx_call(function() eds$factor_names(), "factor_names()", character())
  samples <- .vx_call(function() eds$samples(), "samples()", data.frame())
  errors <- c(fnames$error, samples$error)
  missing_factors <- setdiff(fnames$value, colnames(samples$value))
  if (length(missing_factors)) {
    errors <- c(errors, sprintf("factor_names not in samples() columns: %s", paste(missing_factors, collapse = ", ")))
  }
  pg <- .vx_call(function() eds$primary_grouping(), "primary_grouping()", NA_character_)
  errors <- c(errors, pg$error)
  if (!is.na(pg$value) && !pg$value %in% fnames$value) {
    errors <- c(errors, sprintf("primary_grouping '%s' is not a member of factor_names.", pg$value))
  }
  list(errors = errors, samples = samples$value)
}

# 2. Every contrast has a DE table with the two DE statistics the spec
# requires (section 5); the app's volcano, tables and enrichment inputs fail
# without them.
.vx_contrasts <- function(eds) {
  cs <- .vx_call(function() eds$contrasts(), "contrasts()", character())
  pc <- .vx_call(function() eds$primary_contrast(), "primary_contrast()", NA_character_)
  errors <- c(cs$error, pc$error)
  if (!is.na(pc$value) && !pc$value %in% cs$value) {
    errors <- c(errors, sprintf("primary_contrast '%s' is not a member of contrasts().", pc$value))
  }
  for (cn in cs$value) {
    de <- .vx_call(function() eds$de(cn), sprintf("de('%s')", cn), NULL)
    errors <- c(errors, de$error)
    if (is.null(de$value)) {
      next
    }
    if (!is.data.frame(de$value)) {
      errors <- c(errors, sprintf("de('%s') is not a data.frame.", cn))
      next
    }
    missing_stats <- setdiff(c("effect_log2fc", "score_fdr"), colnames(de$value))
    if (length(missing_stats)) {
      errors <- c(errors, sprintf("de('%s') lacks column(s): %s", cn, paste(missing_stats, collapse = ", ")))
    }
  }
  errors
}

# 3. The raw and primary assays are offered, and every assay states its scale.
.vx_assays <- function(eds, samples_df) {
  an <- .vx_call(function() eds$assay_names(), "assay_names()", character())
  errors <- an$error
  for (acc in c("raw_assay_name", "primary_assay_name")) {
    nm <- .vx_call(function() eds[[acc]](), paste0(acc, "()"), NA_character_)
    errors <- c(errors, nm$error)
    if (!is.na(nm$value) && !nm$value %in% an$value) {
      errors <- c(errors, sprintf("%s() '%s' is not in assay_names().", acc, nm$value))
    }
  }
  for (a in an$value) {
    log2 <- .vx_call(function() eds$assay_is_log2(a), sprintf("assay_is_log2('%s')", a), FALSE)
    errors <- c(errors, log2$error)
    if (!isTRUE(log2$value) && !isFALSE(log2$value)) {
      errors <- c(errors, sprintf("assay_is_log2('%s') is not TRUE or FALSE.", a))
    }
  }
  c(errors, .vx_raw_assay(eds, an$value, samples_df))
}

# 3. The raw assay is aligned to the samples and features.
.vx_raw_assay <- function(eds, assay_names, samples_df) {
  raw_name <- .vx_call(function() eds$raw_assay_name(), "raw_assay_name()", NA_character_)$value
  if (is.na(raw_name) || !raw_name %in% assay_names) {
    return(character())
  }
  a <- .vx_call(function() eds$raw_assay(), "raw_assay()", NULL)
  fids <- .vx_call(function() eds$feature_ids(), "feature_ids()", NULL)
  errors <- c(a$error, fids$error)
  if (!is.null(a$value)) {
    if (!is.null(colnames(a$value)) && !all(colnames(a$value) %in% rownames(samples_df))) {
      errors <- c(errors, "raw_assay() colnames are not a subset of sample_ids().")
    }
    if (!is.null(fids$value) && !is.null(rownames(a$value)) && !all(rownames(a$value) %in% fids$value)) {
      errors <- c(errors, "raw_assay() rownames are not a subset of feature_ids().")
    }
  }
  errors
}
