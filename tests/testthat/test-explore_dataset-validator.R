# -----------------------------------------------------------------------------
# validate_explore_dataset() — covered indirectly by the per-subclass tests
# via expect_silent(validate_explore_dataset(eds)). This file adds direct
# coverage of failure paths and the non-strict return shape.
# -----------------------------------------------------------------------------

test_that("validator: non-strict returns list(ok, errors)", {
  fx <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)
  out <- validate_explore_dataset(eds, strict = FALSE)
  expect_type(out, "list")
  expect_named(out, c("ok", "errors"))
  expect_true(out$ok)
  expect_identical(out$errors, character(0))
})

test_that("validator: rejects non-ExploreDataset input", {
  expect_error(
    validate_explore_dataset(list(a = 1), strict = TRUE),
    "Object does not inherit from 'ExploreDataset'"
  )
})

test_that("validator: catches mismatched primary_grouping", {
  fx <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)

  # Patch the private config to introduce a deliberate mismatch.
  eds$.__enclos_env__$private$.config$primary_grouping <- "no_such_factor"
  out <- validate_explore_dataset(eds, strict = FALSE)
  expect_false(out$ok)
  expect_true(any(grepl("primary_grouping", out$errors)))
})

test_that("validator: a DE table without the spec's required statistics fails", {
  testthat::skip_if_not_installed("anndataR")
  dds <- make_toy_dds()
  res <- DESeq2::results(dds, contrast = c("condition", "treated", "control"))
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(dds, out = out, omics = "rnaseq", de = list(treated_vs_control = res))
  ad <- anndataR::read_h5ad(out, as = "InMemoryAnnData")
  v <- ad$varm
  v$de_treated_vs_control <- v$de_treated_vs_control[, "score_pvalue", drop = FALSE]
  ad$varm <- v
  out2 <- validate_explore_dataset(new_explore_dataset(ad), strict = FALSE)
  expect_false(out2$ok)
  expect_true(any(grepl("lacks column\\(s\\): effect_log2fc, score_fdr", out2$errors)))
})

test_that("validator: non-strict mode records failing accessors instead of throwing", {
  testthat::skip_if_not_installed("anndataR")
  p <- make_toy_pieces()
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata)
  Broken <- R6::R6Class(
    "Broken",
    inherit = ExploreDatasetH5AD,
    public = list(
      samples = function() stop("samples boom"),
      assay_names = function() stop("assay_names boom")
    )
  )
  res <- validate_explore_dataset(Broken$new(out), strict = FALSE)
  expect_false(res$ok)
  expect_true(any(grepl("samples boom", res$errors)))
  expect_true(any(grepl("assay_names boom", res$errors)))
})
