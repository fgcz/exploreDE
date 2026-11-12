# -----------------------------------------------------------------------------
# Factory + SE-shape dispatch tests
#
# Verifies detect_se_shape() classifies every demo SE correctly and that
# new_explore_dataset() instantiates the right concrete subclass for each.
# -----------------------------------------------------------------------------

test_that("detect_se_shape classifies RNA-seq", {
  fx <- fx_rnaseq_gene()
  expect_identical(detect_se_shape(fx$se), "rna_seq")
})

test_that("detect_se_shape classifies isoform RNA-seq", {
  fx <- fx_rnaseq_isoform()
  expect_identical(detect_se_shape(fx$se), "rna_seq")
})

test_that("detect_se_shape classifies CountQC", {
  fx <- fx_countqc()
  expect_identical(detect_se_shape(fx$se), "count_qc")
})

test_that("detect_se_shape classifies CRISPR CountQC", {
  fx <- fx_crispr_countqc()
  expect_identical(detect_se_shape(fx$se), "crispr_count_qc")
})

test_that("detect_se_shape classifies proteomics", {
  fx <- fx_proteomics_protein()
  expect_identical(detect_se_shape(fx$se), "proteomics")
})

test_that("new_explore_dataset() picks ExploreDatasetSE_RNASeq", {
  fx <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)
  expect_s3_class(eds, "ExploreDatasetSE_RNASeq")
  expect_s3_class(eds, "ExploreDatasetSE")
  expect_s3_class(eds, "ExploreDataset")
})

test_that("new_explore_dataset() picks ExploreDatasetSE_CountQC", {
  fx <- fx_countqc()
  eds <- new_explore_dataset(fx$se)
  expect_s3_class(eds, "ExploreDatasetSE_CountQC")
})

test_that("new_explore_dataset() picks ExploreDatasetSE_CrisprCountQC", {
  fx <- fx_crispr_countqc()
  eds <- new_explore_dataset(fx$se)
  expect_s3_class(eds, "ExploreDatasetSE_CrisprCountQC")
  expect_s3_class(eds, "ExploreDatasetSE_CountQC") # inherits
})

test_that("new_explore_dataset() picks ExploreDatasetSE_Proteomics", {
  fx <- fx_proteomics_protein()
  eds <- new_explore_dataset(fx$se)
  expect_s3_class(eds, "ExploreDatasetSE_Proteomics")
})

test_that("ExploreDataset abstract base rejects direct instantiation", {
  expect_error(ExploreDataset$new(), "abstract")
})

test_that("ExploreDatasetSE intermediate rejects shape-blind use", {
  fx <- fx_rnaseq_gene()
  expect_error(ExploreDatasetSE$new(fx$se), "must implement .build_config")
})

test_that("eds_load_se refuses an .RData path that is not a SUSHI EzResult", {
  expect_error(eds_load_se(file.path(tempdir(), "analysis.RData")), "not a SUSHI EzResult")
})

test_that("read_se_any reads <stem>.rds or <stem>.qs2 from a directory, NULL when neither exists", {
  dir <- withr::local_tempdir()
  expect_null(read_se_any(dir, "deResult"))
  saveRDS(list(kind = "rds"), file.path(dir, "deResult.rds"))
  expect_identical(read_se_any(dir, "deResult")$kind, "rds")
  skip_if_not_installed("qs2")
  qs2::qs_save(list(kind = "qs2"), file.path(dir, "rawData.qs2"))
  expect_identical(read_se_any(dir, "rawData")$kind, "qs2")
})
