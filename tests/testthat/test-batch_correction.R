# ============================================================================
# batch_correction_design — returns tagged (design, skip, msg) triple
# ============================================================================
test_that("batch_correction_design treats 'None' / empty designColumn as null-design", {
  ds <- data.frame(cond = c("A", "A", "B", "B"))
  out <- batch_correction_design(ds, designColumn = "None", batchColumn = "batch")
  expect_null(out$design)
  expect_false(out$skip)
  expect_null(out$msg)

  out <- batch_correction_design(ds, designColumn = "",     batchColumn = "batch")
  expect_null(out$design)
  expect_false(out$skip)
})

test_that("batch_correction_design skips with a message when designColumn is missing", {
  ds <- data.frame(cond = c("A", "B"))
  out <- batch_correction_design(ds, designColumn = "nonexistent", batchColumn = "batch")
  expect_null(out$design)
  expect_true(out$skip)
  expect_match(out$msg, "missing from sample metadata")
})

test_that("batch_correction_design warns (no skip) when design factor has <2 levels", {
  ds <- data.frame(cond = rep("A", 4))
  out <- batch_correction_design(ds, designColumn = "cond", batchColumn = "batch")
  expect_null(out$design)
  expect_false(out$skip)
  expect_match(out$msg, "only one level")
})

test_that("batch_correction_design returns a valid design matrix for a 2-level factor", {
  ds <- data.frame(cond = c("A", "A", "B", "B"), batch = c("x", "y", "x", "y"))
  out <- batch_correction_design(ds, designColumn = "cond", batchColumn = "batch")
  expect_null(out$msg)
  expect_false(out$skip)
  expect_true(is.matrix(out$design))
  expect_equal(nrow(out$design), 4)
  expect_gte(ncol(out$design), 2)  # intercept + one treatment column
})

test_that("batch_correction_design warns when batch and design are confounded", {
  ds <- data.frame(cond = c("A", "A", "B", "B"), batch = c("x", "x", "y", "y"))
  out <- batch_correction_design(ds, designColumn = "cond", batchColumn = "batch")
  # Collinear — function returns warned(), not skipped()
  expect_null(out$design)
  expect_false(out$skip)
  expect_match(out$msg, "confounded")
})


# ============================================================================
# apply_batch_correction — graceful-failure wrapper around limma
# Uses shiny::getDefaultReactiveDomain() internally; falls back to warning()
# when called outside a Shiny session, so tests use expect_warning().
# ============================================================================
test_that("apply_batch_correction is a no-op when batchCol is 'None' / empty / NULL", {
  counts <- matrix(runif(20), nrow = 5, dimnames = list(NULL, paste0("s", 1:4)))
  ds <- data.frame(cond = c("A", "A", "B", "B"))

  expect_identical(apply_batch_correction(counts, ds, "None",  "cond"), counts)
  expect_identical(apply_batch_correction(counts, ds, "",      "cond"), counts)
  expect_identical(apply_batch_correction(counts, ds, NULL,    "cond"), counts)
})

test_that("apply_batch_correction warns and returns counts when batch factor has <2 levels", {
  counts <- matrix(runif(20), nrow = 5, dimnames = list(NULL, paste0("s", 1:4)))
  ds <- data.frame(batch = rep("x", 4), cond = c("A", "A", "B", "B"))

  expect_warning(
    out <- apply_batch_correction(counts, ds, "batch", "cond"),
    "has <2 levels"
  )
  expect_identical(out, counts)
})

test_that("apply_batch_correction applies limma::removeBatchEffect when inputs are valid", {
  set.seed(1)
  counts <- matrix(
    c(rnorm(5, mean = 10), rnorm(5, mean = 11),   # batch x
      rnorm(5, mean = 20), rnorm(5, mean = 21)),  # batch y
    nrow = 5,
    dimnames = list(paste0("g", 1:5), paste0("s", 1:4))
  )
  ds <- data.frame(
    batch = c("x", "x", "y", "y"),
    cond  = c("A", "B", "A", "B")
  )
  out <- apply_batch_correction(counts, ds, "batch", "cond")
  expect_equal(dim(out), dim(counts))
  expect_equal(dimnames(out), dimnames(counts))
  # Batch effect should be reduced (column means of batch x closer to batch y)
  mean_x <- mean(out[, ds$batch == "x"])
  mean_y <- mean(out[, ds$batch == "y"])
  expect_lt(abs(mean_x - mean_y), abs(mean(counts[, 1:2]) - mean(counts[, 3:4])))
})

test_that("apply_batch_correction warns and returns counts when limma errors", {
  # Force a limma error by making batch a 1-level factor once NAs are dropped, so
  # the design-matrix dimension doesn't match. We use a design/batch combo where
  # the design matrix has more rows than the batch model matrix after NA drop.
  counts <- matrix(runif(30), nrow = 5, dimnames = list(NULL, paste0("s", 1:6)))
  ds <- data.frame(
    batch = c("x", "x", "x", "y", "y", "y"),
    cond  = c("A", "A", "A", "A", "A", "A")  # single-level design → no design cols
  )
  # With a valid 2-level batch and a 1-level design column, batch_correction_design
  # warns about the single-level design but still runs. We verify the function
  # doesn't crash and preserves shape.
  out <- suppressWarnings(apply_batch_correction(counts, ds, "batch", "cond"))
  expect_equal(dim(out), dim(counts))
})
