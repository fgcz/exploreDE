test_that("eds_finalise_de adds the legacy aliases from the canonical columns", {
  d <- data.frame(effect_log2fc = c(1, -1), score_pvalue = c(0.2, 0.01), score_fdr = c(0.3, 0.02))
  out <- eds_finalise_de(d)
  expect_equal(out$log2Ratio, c(-1, 1))
  expect_equal(out$pValue, c(0.01, 0.2))
  expect_equal(out$fdr, c(0.02, 0.3))
})

test_that("eds_finalise_de keeps legacy columns that are already there", {
  d <- data.frame(log2Ratio = 5, effect_log2fc = 1, pValue = 0.5, fdr = 0.5)
  expect_equal(eds_finalise_de(d)$log2Ratio, 5)
})

test_that("eds_finalise_de reads used_in_test / usedInTest flags as logical", {
  d <- data.frame(pValue = c(0.1, 0.2, 0.3), used_in_test = c(1, 0, 1))
  expect_identical(eds_finalise_de(d)$usedInTest, c(TRUE, FALSE, TRUE))
  d <- data.frame(pValue = c(0.1, 0.2, 0.3), usedInTest = c("True", "False", "TRUE"))
  expect_identical(eds_finalise_de(d)$usedInTest, c(TRUE, FALSE, TRUE))
})

test_that("eds_finalise_de derives usedInTest from the p-value, else the fdr, else TRUE", {
  d <- data.frame(pValue = c(0.1, NA), fdr = c(0.2, 0.3))
  expect_identical(eds_finalise_de(d)$usedInTest, c(TRUE, FALSE))
  d <- data.frame(pValue = c(NA_real_, NA_real_), fdr = c(0.2, NA))
  expect_identical(eds_finalise_de(d)$usedInTest, c(TRUE, FALSE))
  d <- data.frame(log2Ratio = c(1, 2))
  expect_identical(eds_finalise_de(d)$usedInTest, c(TRUE, TRUE))
})

test_that("eds_finalise_de orders by p-value, and by fdr when no p-value is usable", {
  d <- data.frame(feature_id = c("a", "b", "c"), pValue = c(0.5, 0.1, NA), fdr = c(0.1, 0.2, 0.3))
  expect_equal(eds_finalise_de(d)$feature_id, c("b", "a", "c"))
  d$pValue <- NA_real_
  expect_equal(eds_finalise_de(d)$feature_id, c("a", "b", "c"))
  d$fdr <- c(0.3, 0.2, 0.1)
  expect_equal(eds_finalise_de(d)$feature_id, c("c", "b", "a"))
})

test_that("eds_finalise_de handles a table without rows", {
  d <- data.frame(pValue = numeric(), fdr = numeric())
  out <- eds_finalise_de(d)
  expect_equal(nrow(out), 0L)
  expect_true("usedInTest" %in% colnames(out))
})
