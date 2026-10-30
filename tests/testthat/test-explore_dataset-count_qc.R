# -----------------------------------------------------------------------------
# ExploreDatasetSE_CountQC + ExploreDatasetSE_CrisprCountQC tests
# -----------------------------------------------------------------------------

test_that("CountQC: subclass + no DE", {
  fx <- fx_countqc()
  eds <- new_explore_dataset(fx$se)

  expect_s3_class(eds, "ExploreDatasetSE_CountQC")
  expect_identical(eds$omics_class(), "genomics")
  expect_identical(eds$data_type(), "CountQC")
  expect_false(eds$has_de())
  expect_identical(length(eds$contrasts()), 0L)
  expect_silent(validate_explore_dataset(eds))
})

test_that("CountQC: assay accessor returns counts", {
  fx <- fx_countqc()
  eds <- new_explore_dataset(fx$se)
  expect_true("counts" %in% eds$assay_names())
  m <- eds$assay("counts")
  expect_identical(ncol(m), eds$n_samples())
  expect_identical(nrow(m), eds$n_features())
})

test_that("CrisprCountQC: subclass + sgrna_to_gene populated", {
  fx <- fx_crispr_countqc()
  eds <- new_explore_dataset(fx$se)

  expect_s3_class(eds, "ExploreDatasetSE_CrisprCountQC")
  expect_s3_class(eds, "ExploreDatasetSE_CountQC") # inherits
  expect_identical(eds$omics_class(), "genomics")
  expect_identical(eds$data_type(), "CrisprCountQC")
  expect_identical(eds$feature_id_kind(), "guide")

  sg <- eds$sgrna_to_gene()
  expect_false(is.null(sg))
  expect_identical(length(sg), eds$n_features())
  expect_identical(names(sg), eds$feature_ids())
})

test_that("CrisprCountQC: validator passes", {
  fx <- fx_crispr_countqc()
  eds <- new_explore_dataset(fx$se)
  expect_silent(validate_explore_dataset(eds))
})
