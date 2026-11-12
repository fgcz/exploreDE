# Which assays are on the log2 scale (eds$assay_is_log2()), per backend.

log2_assays <- function(eds) {
  log2 <- vapply(eds$assay_names(), eds$assay_is_log2, logical(1))
  names(log2)[log2]
}

test_that("RNA-seq and Count QC: log2 counts and VST are log2, the other assays are not", {
  rna <- new_explore_dataset(system.file("examples", "genomics", "deResult.rds", package = "exploreDE"))
  expect_setequal(log2_assays(rna), c("normalised_log2", "vst"))

  m <- matrix(rpois(60, 50), 10, dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
  cd <- data.frame("Condition [Factor]" = rep(c("A", "B"), each = 3), check.names = FALSE, row.names = colnames(m))
  qc <- new_explore_dataset(SummarizedExperiment::SummarizedExperiment(assays = list(counts = m), colData = cd))
  expect_s3_class(qc, "ExploreDatasetSE_CountQC")
  expect_setequal(log2_assays(qc), "vst")

  skip_if_not_installed("anndataR")
  p <- make_toy_pieces()
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata)
  h5 <- new_explore_dataset(out)
  expect_setequal(log2_assays(h5), c("normalised_log2", "vst"))
  expect_false(h5$assay_is_log2("counts"))
})

test_that("an assay the dataset does not have is an error", {
  rna <- new_explore_dataset(system.file("examples", "genomics", "deResult.rds", package = "exploreDE"))
  expect_error(rna$assay_is_log2("intensities"))
})

test_that("proteomics SE: the transformed intensities are log2; a first assay taken in their place is not", {
  se <- readRDS(system.file("examples", "proteomics", "SummarizedExperiment.rds", package = "exploreDE"))
  expect_true(new_explore_dataset(se)$assay_is_log2("intensities"))
  # Without transformedData the first assay (rawData) becomes "intensities"
  raw_only <- se
  SummarizedExperiment::assays(raw_only) <- SummarizedExperiment::assays(se)["rawData"]
  eds <- new_explore_dataset(raw_only)
  expect_s3_class(eds, "ExploreDatasetSE_Proteomics")
  expect_false(eds$assay_is_log2("intensities"))
})

test_that("proteomics .h5ad: transformed intensities are log2, raw ones are not, also in a raw-only file", {
  skip_if_not_installed("anndataR")
  pr <- make_toy_prot("protein")
  both <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(
    list(raw = 2^pr$intensities, transformed = pr$intensities),
    out = both,
    omics = "proteomics",
    coldata = pr$coldata,
    rowdata = pr$rowdata
  )
  eds <- new_explore_dataset(both)
  expect_true(eds$assay_is_log2("intensities"))
  expect_false(eds$assay_is_log2("intensities_raw"))

  # Raw intensities only: no layer is log2, whatever it is offered as
  raw <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(
    list(raw = 2^pr$intensities),
    out = raw,
    omics = "proteomics",
    coldata = pr$coldata,
    rowdata = pr$rowdata
  )
  eds <- new_explore_dataset(raw)
  expect_true("intensities" %in% eds$assay_names())
  expect_identical(log2_assays(eds), character(0))
})
