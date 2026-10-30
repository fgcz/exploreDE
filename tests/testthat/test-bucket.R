# ============================================================================
# Feature-bucket API tests — bucket_add / bucket_remove / bucket_clear
#
# These functions reference `featureBucket`, `seqAnnoReactive` and
# `inputDataReactive` via lexical scope. We rebind the function's environment
# to one holding our mock reactives (via with_mock_env, defined in
# helper-mock_reactives.R).
#
# All reactiveValues reads need a reactive context, so every call is wrapped
# in shiny::isolate().
# ============================================================================

make_bucket_mocks <- function(feature_ids = c("A", "B", "C"), sgRNA_to_gene = NULL, initial = character(0)) {
  list(
    featureBucket = shiny::reactiveValues(features = initial, version = 0L),
    seqAnnoReactive = shiny::reactiveValues(sa = data.frame(feature_id = feature_ids)),
    inputDataReactive = function() list(sgRNA_to_gene = sgRNA_to_gene)
  )
}

# --- bucket_add -------------------------------------------------------------
test_that("bucket_add adds only valid features from the seqAnno feature_id set", {
  m <- make_bucket_mocks()
  fn <- with_mock_env(bucket_add, !!!m)

  shiny::isolate(fn(c("A", "B", "not_a_gene", "")))
  expect_equal(shiny::isolate(m$featureBucket$features), c("A", "B"))
  expect_equal(shiny::isolate(m$featureBucket$version), 1L)
})

test_that("bucket_add does not re-add features already in the bucket (idempotent)", {
  m <- make_bucket_mocks(initial = c("A"))
  fn <- with_mock_env(bucket_add, !!!m)

  shiny::isolate(fn(c("A", "B")))
  expect_equal(shiny::isolate(m$featureBucket$features), c("A", "B"))
  # Version bumped once (for the real addition of "B", not for "A")
  expect_equal(shiny::isolate(m$featureBucket$version), 1L)

  # Add "A" again → no change, no version bump
  shiny::isolate(fn("A"))
  expect_equal(shiny::isolate(m$featureBucket$features), c("A", "B"))
  expect_equal(shiny::isolate(m$featureBucket$version), 1L)
})

test_that("bucket_add preserves insertion order of new features", {
  m <- make_bucket_mocks(feature_ids = paste0("gene", 1:5))
  fn <- with_mock_env(bucket_add, !!!m)

  shiny::isolate(fn(c("gene3", "gene1", "gene4")))
  expect_equal(shiny::isolate(m$featureBucket$features), c("gene3", "gene1", "gene4"))
})

test_that("bucket_add accepts gene-level names when sgRNA_to_gene is present (CrisprCountQC)", {
  m <- make_bucket_mocks(
    feature_ids = c("sg1", "sg2"),
    sgRNA_to_gene = c(sg1 = "GENE_A", sg2 = "GENE_B")
  )
  fn <- with_mock_env(bucket_add, !!!m)

  shiny::isolate(fn(c("GENE_A", "sg2")))
  expect_setequal(shiny::isolate(m$featureBucket$features), c("GENE_A", "sg2"))
})


# --- bucket_remove ----------------------------------------------------------
test_that("bucket_remove drops features from the bucket", {
  m <- make_bucket_mocks(initial = c("A", "B", "C"))
  fn <- with_mock_env(bucket_remove, !!!m)

  shiny::isolate(fn("B"))
  expect_equal(shiny::isolate(m$featureBucket$features), c("A", "C"))
  expect_equal(shiny::isolate(m$featureBucket$version), 1L)
})

test_that("bucket_remove is a no-op when the feature is not in the bucket", {
  m <- make_bucket_mocks(initial = c("A"))
  fn <- with_mock_env(bucket_remove, !!!m)

  shiny::isolate(fn("Z"))
  expect_equal(shiny::isolate(m$featureBucket$features), "A")
  expect_equal(shiny::isolate(m$featureBucket$version), 0L) # no bump
})


# --- bucket_clear -----------------------------------------------------------
test_that("bucket_clear empties the bucket and bumps the version", {
  m <- make_bucket_mocks(initial = c("A", "B", "C"))
  fn <- with_mock_env(bucket_clear, !!!m)

  shiny::isolate(fn())
  expect_equal(shiny::isolate(m$featureBucket$features), character(0))
  expect_equal(shiny::isolate(m$featureBucket$version), 1L)
})
