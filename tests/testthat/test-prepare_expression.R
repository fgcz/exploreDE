# prepare_expression(): the log2 and batch steps every plotting tab applies
# to an assay. Outside a Shiny session its messages are warnings.

test_that("prepare_expression logs linear data, never log2 data", {
  m <- matrix(c(1, 3, 7, 15), 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  s <- data.frame(row.names = colnames(m))
  expect_identical(
    prepare_expression(m, s, is_log2 = FALSE, take_log2 = TRUE),
    list(counts = log2(m + 1), is_log2 = TRUE, batch = "None")
  )
  expect_identical(prepare_expression(m, s, is_log2 = TRUE, take_log2 = TRUE)$counts, m)
  expect_identical(prepare_expression(m, s, is_log2 = FALSE, take_log2 = FALSE)$counts, m)
})

test_that("prepare_expression batch-corrects log2 data and refuses linear data", {
  set.seed(3)
  m <- matrix(rnorm(20, 10), 5, dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  s <- data.frame(row.names = colnames(m), batch = c("x", "x", "y", "y"), cond = c("A", "B", "A", "B"))
  corrected <- apply_batch_correction(m, s, "batch", "cond")
  expect_false(identical(corrected, m))

  out <- prepare_expression(m, s, is_log2 = TRUE, batch = "batch", design = "cond")
  expect_identical(out$counts, corrected)
  expect_identical(out$batch, "batch")

  # Linear data is logged first when asked for, then corrected
  expect_equal(
    prepare_expression(2^m - 1, s, is_log2 = FALSE, take_log2 = TRUE, batch = "batch", design = "cond")$counts,
    corrected
  )

  expect_warning(
    lin <- prepare_expression(2^m, s, is_log2 = FALSE, batch = "batch", design = "cond"),
    "requires log-scale"
  )
  expect_identical(lin$counts, 2^m)
  expect_identical(lin$batch, "None")

  # No batch column: nothing to say
  expect_no_warning(none <- prepare_expression(2^m, s, is_log2 = FALSE, batch = "None"))
  expect_identical(none$counts, 2^m)
})

test_that("prepare_expression reports a batch correction that could not run as not applied", {
  m <- matrix(runif(20), 5, dimnames = list(NULL, paste0("s", 1:4)))
  s <- data.frame(row.names = colnames(m), batch = c("x", NA, "y", "y"), cond = c("A", "A", "B", "B"))
  expect_warning(out <- prepare_expression(m, s, is_log2 = TRUE, batch = "batch", design = "cond"), "missing values")
  expect_identical(out$counts, m)
  expect_identical(out$batch, "None")
})
