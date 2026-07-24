# -----------------------------------------------------------------------------
# ExploreDatasetSE_RNASeq tests
# -----------------------------------------------------------------------------

test_that("RNA-seq gene-level: accessors return populated values", {
  fx  <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)

  expect_identical(eds$omics_class(), "genomics")
  expect_identical(eds$data_type(), "RNASeq")
  expect_identical(eds$feature_id_kind(), "gene")
  expect_identical(eds$backend_kind(), "se")

  expect_gt(eds$n_samples(), 0L)
  expect_gt(eds$n_features(), 0L)
  expect_gt(length(eds$factor_names()), 0L)

  s <- eds$samples()
  expect_s3_class(s, "data.frame")
  expect_identical(nrow(s), eds$n_samples())

  f <- eds$features()
  expect_s3_class(f, "data.frame")
  expect_true("feature_id"   %in% colnames(f))
  expect_true("display_name" %in% colnames(f))
})

test_that("RNA-seq gene-level: DE accessor", {
  fx  <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)

  expect_gt(length(eds$contrasts()), 0L)
  expect_true(eds$has_de())

  c1 <- eds$primary_contrast()
  expect_true(c1 %in% eds$contrasts())

  de <- eds$de(c1)
  expect_s3_class(de, "data.frame")
  expect_gt(nrow(de), 0L)
  expect_true("effect_log2fc" %in% colnames(de))
  expect_true("score_pvalue"  %in% colnames(de))
  expect_true("score_fdr"     %in% colnames(de))

  cm <- eds$contrast_meta(c1)
  expect_identical(cm$method, "DESeq2")
  expect_false(is.null(cm$test_group))
  expect_false(is.null(cm$reference_group))
})

test_that("RNA-seq gene-level: assay accessor", {
  fx  <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)

  an <- eds$assay_names()
  expect_true("counts" %in% an)

  counts <- eds$assay("counts")
  expect_true(is.matrix(counts) || methods::is(counts, "Matrix"))
  expect_identical(ncol(counts), eds$n_samples())
  expect_identical(nrow(counts), eds$n_features())
})

test_that("RNA-seq gene-level: validator passes", {
  fx  <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)
  expect_silent(validate_explore_dataset(eds))
})

test_that("RNA-seq isoform-level: feature_id_kind = isoform", {
  fx  <- fx_rnaseq_isoform()
  eds <- new_explore_dataset(fx$se)
  expect_identical(eds$omics_class(), "genomics")
  expect_identical(eds$data_type(), "RNASeq")
  expect_identical(eds$feature_id_kind(), "isoform")
  expect_silent(validate_explore_dataset(eds))
})
