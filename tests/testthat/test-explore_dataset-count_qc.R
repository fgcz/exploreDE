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

test_that("CountQC isoform-level: feature_id keyed on transcript_id, no duplicates", {
  # Isoform CountQC (e.g. Kallisto) has several rows per gene. Keying on
  # gene_id gave duplicate feature_ids and new_explore_dataset() errored with
  # "duplicate 'row.names' are not allowed" (p38719 o38990, 2026-11-04).
  fx <- fx_countqc()
  se <- fx$se[1:3, ]
  se <- se[rep(seq_len(nrow(se)), each = 2), ] # 2 transcripts per gene
  tx <- paste0("ENST", sprintf("%011d", seq_len(nrow(se))))
  SummarizedExperiment::rowData(se)$transcript_id <- tx
  rownames(se) <- tx
  md <- S4Vectors::metadata(se)
  md$param$featureLevel <- "isoform"
  S4Vectors::metadata(se) <- md

  eds <- new_explore_dataset(se)
  expect_s3_class(eds, "ExploreDatasetSE_CountQC")
  expect_identical(eds$feature_id_kind(), "isoform")
  fid <- eds$feature_ids()
  expect_identical(length(fid), nrow(se))
  expect_false(anyDuplicated(fid) > 0)
  expect_true(all(grepl("_ENST", fid))) # every gene here has 2 transcripts
  expect_identical(sort(rownames(eds$assay("counts"))), sort(fid))
  expect_silent(validate_explore_dataset(eds))
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
