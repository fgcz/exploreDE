counts <- matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
de <- data.frame(
  feature_id = c("g2", "g1", "g3"),
  log2Ratio = c(1, -1, 0),
  pValue = c(0.01, 0.2, 1),
  fdr = c(0.05, 0.4, 1)
)

test_that("heatmap_export_table gives the plotted values without DE statistics unless asked", {
  out <- heatmap_export_table(counts, de, FALSE)
  expect_identical(colnames(out), c("feature_id", "s1", "s2"))
  expect_identical(out$feature_id, c("g1", "g2"))
  expect_identical(colnames(heatmap_export_table(counts, de, NULL)), c("feature_id", "s1", "s2"))
})

test_that("heatmap_export_table adds the DE statistics by feature id", {
  out <- heatmap_export_table(counts, de, TRUE)
  expect_identical(colnames(out), c("feature_id", "s1", "s2", "log2Ratio", "pValue", "fdr"))
  expect_identical(out$feature_id, c("g1", "g2"))
  expect_identical(out$log2Ratio, c(-1, 1))
})

test_that("heatmap_export_table adds only the statistics the DE table has, and none without one", {
  expect_identical(
    colnames(heatmap_export_table(counts, de[, c("feature_id", "fdr")], TRUE)),
    c("feature_id", "s1", "s2", "fdr")
  )
  expect_identical(colnames(heatmap_export_table(counts, NULL, TRUE)), c("feature_id", "s1", "s2"))
})
