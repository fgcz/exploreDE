# Assays computed from raw counts: VST (eds_vst_matrix()), the normalised
# counts of an .h5ad built from raw counts alone, and the normalised layer
# build_explore_h5ad() takes from a DESeqDataSet.

toy_counts <- function(n, mean, seed = 1) {
  set.seed(seed)
  matrix(
    stats::rnbinom(n * 6, mu = mean, size = 10),
    n,
    6,
    dimnames = list(paste0("g", seq_len(n)), paste0("s", 1:6))
  )
}

test_that("eds_vst_matrix uses vst() when 1000 features have a mean normalised count above 5", {
  big <- toy_counts(1200, 50)
  expect_equal(eds_vst_matrix(big), DESeq2::vst(big, blind = FALSE))
  # Too few expressed features for vst(): the full transformation
  low <- toy_counts(1200, 2)
  expect_error(DESeq2::vst(low, blind = FALSE), "nsub")
  expect_equal(eds_vst_matrix(low), DESeq2::varianceStabilizingTransformation(low, blind = FALSE))
  small <- toy_counts(200, 50)
  expect_equal(eds_vst_matrix(small), DESeq2::varianceStabilizingTransformation(small, blind = FALSE))
})

test_that("an .h5ad of raw counts is normalised by size factors, or by library size when every feature has a zero", {
  skip_if_not_installed("anndataR")
  normalised <- function(counts, coldata) {
    out <- withr::local_tempfile(fileext = ".h5ad")
    build_explore_h5ad(counts, out = out, omics = "rnaseq", coldata = coldata)
    new_explore_dataset(out)$assay("normalised")
  }
  p <- make_toy_pieces()
  cts <- p$counts
  expect_equal(normalised(cts, p$coldata), sweep(cts, 2, DESeq2::estimateSizeFactorsForMatrix(cts), "/"))

  zeros <- cts
  zeros[cbind(seq_len(nrow(zeros)), rep_len(seq_len(ncol(zeros)), nrow(zeros)))] <- 0L
  expect_equal(normalised(zeros, p$coldata), sweep(zeros, 2, colSums(zeros) / mean(colSums(zeros)), "/"))
})

test_that("build_explore_h5ad() stores a DESeqDataSet's normalised counts when it has size factors", {
  skip_if_not_installed("anndataR")
  layers <- function(dds) {
    out <- withr::local_tempfile(fileext = ".h5ad")
    build_explore_h5ad(dds, out = out, omics = "rnaseq")
    names(anndataR::read_h5ad(out, as = "InMemoryAnnData")$layers)
  }
  dds <- make_toy_dds()
  expect_true("normalised" %in% layers(dds))
  # colData without the sizeFactor column DESeq() added
  cd <- SummarizedExperiment::colData(dds)[, "condition", drop = FALSE]
  bare <- DESeq2::DESeqDataSetFromMatrix(DESeq2::counts(dds), cd, ~condition)
  expect_false("normalised" %in% layers(bare))
})
