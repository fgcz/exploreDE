# ============================================================================
# Feature bucket: bucket_add() / bucket_remove() return the bucket after the
# change. server-featureBucket.R stores it (fb_add / fb_remove / fb_clear).
# ============================================================================

test_that("bucket_add adds only known features", {
  expect_equal(bucket_add(character(0), c("A", "B", "not_a_gene", ""), known = c("A", "B", "C")), c("A", "B"))
})

test_that("bucket_add does not add a feature twice", {
  expect_equal(bucket_add("A", c("A", "B"), known = c("A", "B", "C")), c("A", "B"))
  expect_equal(bucket_add(c("A", "B"), "A", known = c("A", "B")), c("A", "B"))
})

test_that("bucket_add keeps the order the new features come in", {
  expect_equal(
    bucket_add(character(0), c("gene3", "gene1", "gene4"), known = paste0("gene", 1:5)),
    c("gene3", "gene1", "gene4")
  )
})

test_that("bucket_remove drops features; a feature not in the bucket changes nothing", {
  expect_equal(bucket_remove(c("A", "B", "C"), "B"), c("A", "C"))
  expect_equal(bucket_remove("A", "Z"), "A")
})
