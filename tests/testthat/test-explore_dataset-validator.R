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
