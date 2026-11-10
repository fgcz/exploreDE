# -----------------------------------------------------------------------------
# sample_log2_ratio(), ratio_heatmap_matrix(), contrast_log2_ratio()
# (R/helperFunctions.R): the Heatmap tab's "Sample ratios" setting and its
# "Ratios" heatmap.
# -----------------------------------------------------------------------------

ratio_mat <- function() {
  matrix(
    c(1, 2, 3, 4, 5, 5, 1, 1, NA, 2, 2, NA, NA, NA, 1, 1),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3", "g4"), c("A1", "A2", "B1", "B2"))
  )
}

test_that("sample_log2_ratio is the difference of the sample-set means, per feature", {
  r <- sample_log2_ratio(ratio_mat(), c("A1", "A2"), c("B1", "B2"))
  expect_identical(names(r), c("g1", "g2", "g3", "g4"))
  expect_equal(r[["g1"]], 1.5 - 3.5)
  expect_equal(r[["g2"]], 5 - 1)
  # a missing value is left out of the mean
  expect_equal(r[["g3"]], 2 - 2)
  # no numerator value at all
  expect_true(is.nan(r[["g4"]]))
  # one sample on each side
  expect_equal(sample_log2_ratio(ratio_mat(), "A1", "B1")[["g1"]], 1 - 3)
})

test_that("sample_log2_ratio refuses a sample that is not a column", {
  expect_error(sample_log2_ratio(ratio_mat(), "A1", "C9"), "C9")
})

test_that("ratio_heatmap_matrix aligns the ratios on their shared features and ranks them", {
  ratios <- list(
    r1 = c(g1 = 1, g2 = 0, g3 = 3, g4 = NA),
    r2 = c(g1 = -1, g2 = 0, g3 = 3, g5 = 10)
  )
  m <- ratio_heatmap_matrix(ratios, top_n = 10, rank_by = "variance")
  expect_identical(colnames(m), c("r1", "r2"))
  # g4 is in r1 only and g5 in r2 only
  expect_setequal(rownames(m), c("g1", "g2", "g3"))
  expect_identical(rownames(m)[1], "g1")
  expect_equal(m["g1", ], c(r1 = 1, r2 = -1))

  m2 <- ratio_heatmap_matrix(ratios, top_n = 2, rank_by = "absolute")
  expect_identical(rownames(m2), c("g3", "g1"))
  expect_identical(dim(m2), c(2L, 2L))
})

test_that("ratio_heatmap_matrix keeps only features finite in every ratio and copes with nothing left", {
  ratios <- list(r1 = c(g1 = NA, g2 = 2, g3 = 1), r2 = c(g1 = NaN, g2 = 1, g3 = Inf))
  m <- ratio_heatmap_matrix(ratios, top_n = 5)
  expect_identical(rownames(m), "g2")
  expect_false(anyNA(m))

  empty <- ratio_heatmap_matrix(list(r1 = c(g1 = 1), r2 = c(g2 = 1)), top_n = 5)
  expect_identical(nrow(empty), 0L)
  expect_identical(colnames(empty), c("r1", "r2"))
  expect_identical(nrow(ratio_heatmap_matrix(list(), top_n = 5)), 0L)
  expect_identical(nrow(ratio_heatmap_matrix(list(r1 = NULL), top_n = 5)), 0L)
})

test_that("ratio_heatmap_matrix ranks a single ratio by its absolute value", {
  m <- ratio_heatmap_matrix(list(r1 = c(g1 = 0.5, g2 = -4, g3 = 2)), top_n = 2, rank_by = "variance")
  expect_identical(rownames(m), c("g2", "g3"))
})

test_that("contrast_log2_ratio is the contrast's log2 fold change, named by feature id", {
  eds <- new_explore_dataset(fx_rnaseq_gene()$se)
  con <- eds$primary_contrast()
  r <- contrast_log2_ratio(eds, con)
  de <- eds$de(con)
  expect_identical(names(r), de$feature_id)
  expect_equal(unname(r), de$log2Ratio)
})

test_that("contrast_log2_ratio is NULL when the DE table has no log2Ratio", {
  eds <- list(de = function(contrast) data.frame(feature_id = c("g1", "g2"), pValue = c(0.1, 0.2)))
  expect_null(contrast_log2_ratio(eds, "A_vs_B"))
})
