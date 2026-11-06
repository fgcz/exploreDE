# legacy_multi_contrast_from_eds(): the one converter from an ExploreDataset to
# the per-contrast record the Multi DE modules read. Every backend must give
# the fields those modules use.

expect_multi_record <- function(rec, data_type) {
  expect_identical(rec$dataType, data_type)
  expect_true(all(c("feature_id", "log2Ratio", "pValue", "fdr", "usedInTest") %in% colnames(rec$seqAnno)))
  expect_gt(length(rec$countList), 0L)
  expect_true(all(c("dataset", "factors", "groupingName", "genes", "param", "design") %in% names(rec)))
  expect_match(rec$param$comparison, "--over--", fixed = TRUE)
  expect_true(all(rec$genes %in% rownames(rec$countList[[1]])))
}

build_toy_prot_h5ad <- function(env = parent.frame()) {
  p <- make_toy_prot("protein")
  rd <- p$rowdata
  rd$protein_Id <- sprintf("P%03d", seq_len(nrow(rd)))
  de <- data.frame(
    row.names = rownames(p$intensities),
    diff = seq(-2, 2, length.out = nrow(rd)),
    p.value = seq(0.001, 0.5, length.out = nrow(rd)),
    FDR = seq(0.01, 0.9, length.out = nrow(rd))
  )
  out <- withr::local_tempfile(fileext = ".h5ad", .local_envir = env)
  build_explore_h5ad(
    list(intensities = p$intensities),
    out = out,
    omics = "proteomics",
    coldata = p$coldata,
    rowdata = rd,
    de = list(B_vs_A = de),
    de_columns = list(effect_log2fc = "diff", score_pvalue = "p.value", score_fdr = "FDR")
  )
  out
}

test_that("legacy_multi_contrast_from_eds builds a record for a canonical proteomics .h5ad", {
  skip_if_not_installed("anndataR")
  eds <- new_explore_dataset(build_toy_prot_h5ad())
  rec <- legacy_multi_contrast_from_eds(eds, "B_vs_A")
  expect_multi_record(rec, "proteomics")
  expect_true("intensities" %in% names(rec$countList))
})

test_that("legacy_multi_contrast_from_eds builds a record for a proteomics Parquet directory", {
  skip_if_not_installed("anndataR")
  if (!requireNamespace("nanoparquet", quietly = TRUE) && !requireNamespace("arrow", quietly = TRUE)) {
    skip("neither nanoparquet nor arrow is installed")
  }
  pq <- withr::local_tempdir()
  export_h5ad_to_parquet(build_toy_prot_h5ad(), pq)
  eds <- new_explore_dataset(pq)
  expect_identical(eds$backend_kind(), "parquet")
  expect_multi_record(legacy_multi_contrast_from_eds(eds, "B_vs_A"), "proteomics")
})

test_that("legacy_multi_contrast_from_eds builds a record for a prolfquapp result (SummarizedExperiment)", {
  skip_if_not_installed("anndataR")
  eds <- new_explore_dataset(test_path("fixtures", "prolfquapp_2.10_AnnData.h5ad"))
  expect_identical(eds$backend_kind(), "se")
  rec <- legacy_multi_contrast_from_eds(eds, eds$contrasts()[1])
  expect_multi_record(rec, "proteomics")
  expect_identical(rec$param$comparison, "G_2--over--G_1")
})
