# Log2 and batch correction of an assay before plotting (prepare_expression()), and user messages.

# nolint start: cyclocomp_linter. One guard per unusable design.
#' batch_correction_design
#' Build a design matrix for limma::removeBatchEffect that preserves biological
#' signal. Without a design matrix, removeBatchEffect assumes a one-group
#' experiment and can remove biology that correlates with batch.
#'
#' Never throws — returns a tagged result so callers can notify the user and
#' decide whether to skip batch correction or proceed without a design matrix.
#'
#' @param dataset Sample metadata data.frame
#' @param designColumn Column name to build the design from. "None" / empty
#'   means "no design" (batch correction still runs, without preserving any
#'   biology).
#' @param batchColumn Column name of the batch factor — required for the
#'   collinearity check.
#' @return list(design = <matrix or NULL>, skip = <logical>, msg = <chr or NULL>)
#'   - skip = TRUE  → caller must abort batch correction entirely
#'   - skip = FALSE, design = NULL → run removeBatchEffect without a design
#'   - skip = FALSE, design = <matrix> → run with the design
#' @noRd
batch_correction_design <- function(dataset, designColumn, batchColumn) {
  null_design <- list(design = NULL, skip = FALSE, msg = NULL)
  skipped <- function(msg) list(design = NULL, skip = TRUE, msg = msg)
  warned <- function(msg) list(design = NULL, skip = FALSE, msg = msg)

  # model.matrix() drops rows with NA, so a design or batch with missing values
  # no longer lines up with the samples: skip instead.
  has_batch <- !is.null(batchColumn) &&
    nzchar(batchColumn) &&
    batchColumn != "None" &&
    batchColumn %in% colnames(dataset)
  if (has_batch && anyNA(dataset[[batchColumn]])) {
    return(skipped(sprintf(
      "Batch column '%s' has missing values \u2014 batch correction skipped.",
      batchColumn
    )))
  }
  if (is.null(designColumn) || !nzchar(designColumn) || designColumn == "None") {
    return(null_design)
  }
  if (!designColumn %in% colnames(dataset)) {
    return(skipped(sprintf(
      "Design column '%s' missing from sample metadata \u2014 batch correction skipped.",
      designColumn
    )))
  }
  if (anyNA(dataset[[designColumn]])) {
    return(skipped(sprintf(
      "Design column '%s' has missing values \u2014 batch correction skipped.",
      designColumn
    )))
  }

  grp <- factor(dataset[[designColumn]])
  if (nlevels(grp) < 2) {
    return(warned(sprintf(
      "Design factor '%s' has only one level \u2014 batch correction applied without a design matrix.",
      designColumn
    )))
  }

  design <- stats::model.matrix(~grp)

  if (has_batch) {
    batch <- factor(dataset[[batchColumn]])
    if (nlevels(batch) >= 2) {
      combined <- cbind(stats::model.matrix(~batch), design[, -1, drop = FALSE])
      if (qr(combined)$rank < ncol(combined)) {
        return(warned(sprintf(
          "Batch ('%s') and design ('%s') are confounded \u2014 batch correction applied without a design matrix.",
          batchColumn,
          designColumn
        )))
      }
    }
  }

  list(design = design, skip = FALSE, msg = NULL)
}
# nolint end

# nolint start: cyclocomp_linter. Guards before limma.
#' apply_batch_correction
#' Graceful-failure wrapper around limma::removeBatchEffect. Never throws:
#' on any validation or limma error, returns the original counts unchanged
#' and surfaces a Shiny notification.
#'
#' @param counts Numeric matrix (samples in columns)
#' @param dataset Sample metadata (nrow = ncol(counts))
#' @param batchCol Column name of batch factor, or "None" to skip entirely
#' @param designCol Column name of design factor, or "None"
#' @return The (possibly batch-corrected) counts matrix
#' @noRd
apply_batch_correction <- function(counts, dataset, batchCol, designCol) {
  if (is.null(batchCol) || !nzchar(batchCol) || batchCol == "None") {
    return(counts)
  }

  if (
    !batchCol %in% colnames(dataset) ||
      nlevels(factor(dataset[[batchCol]])) < 2
  ) {
    notify_user(sprintf("Batch factor '%s' has <2 levels \u2014 batch correction skipped.", batchCol))
    return(counts)
  }

  info <- batch_correction_design(dataset, designCol, batchCol)
  if (!is.null(info$msg)) {
    notify_user(info$msg, type = if (isTRUE(info$skip)) "error" else "warning")
  }
  if (isTRUE(info$skip)) {
    return(counts)
  }

  args <- list(x = counts, batch = dataset[[batchCol]])
  if (!is.null(info$design)) {
    args$design <- info$design
  }

  tryCatch(
    do.call(limma::removeBatchEffect, args),
    error = function(e) {
      notify_user(
        sprintf("Batch correction failed: %s \u2014 showing uncorrected counts.", conditionMessage(e)),
        type = "error"
      )
      counts
    }
  )
}
# nolint end

# A message to the user: a notification in a Shiny session, a warning outside.
notify_user <- function(msg, type = "warning") {
  if (!is.null(shiny::getDefaultReactiveDomain())) {
    shiny::showNotification(msg, type = type, duration = 8)
  } else {
    warning(msg, call. = FALSE)
  }
}

# The log2 and batch steps a plotting tab applies to an assay, in this order.
# `is_log2`: the assay is on the log2 scale (eds$assay_is_log2()); it is not
# logged again when `take_log2` asks for log2. Batch correction (`batch` a
# sample column, or "None") needs log2 data: on linear data the user is told
# and the matrix kept. Returns the matrix, whether it is on the log2 scale, and
# the batch column applied ("None" when none was, or the correction could not
# run).
prepare_expression <- function(counts, samples, is_log2, take_log2 = FALSE, batch = "None", design = "None") {
  if (take_log2 && !is_log2) {
    counts <- log2(counts + 1)
    is_log2 <- TRUE
  }
  applied <- "None"
  if (!batch %in% c("", "None")) {
    if (!is_log2) {
      notify_user(
        "Batch correction requires log-scale data. Please select a logged count type or enable the log2 transform."
      )
    } else {
      corrected <- apply_batch_correction(counts, samples, batch, design)
      if (!identical(corrected, counts)) {
        counts <- corrected
        applied <- batch
      }
    }
  }
  list(counts = counts, is_log2 = is_log2, batch = applied)
}
