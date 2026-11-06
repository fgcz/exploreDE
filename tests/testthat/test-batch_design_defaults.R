# batch_design_defaults(): the batch / design dropdowns of the Boxplot,
# Correlation, Heatmap and Ordination tabs (one rule for all four).

test_that("two-factor RNA-seq: batch is the second factor, design locked to the grouping", {
  testthat::skip_if_not_installed("anndataR")
  p <- make_toy_pieces()
  cd <- p$coldata
  cd$batch <- rep(c("b1", "b2"), length.out = nrow(cd))
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = cd, factors = c("condition", "batch"))
  bd <- batch_design_defaults(new_explore_dataset(out))
  expect_identical(bd$batch_choices, c("None", "condition", "batch"))
  expect_identical(bd$batch_default, "batch")
  expect_identical(bd$design_choices, c("None", "condition"))
  expect_identical(bd$design_default, "condition")
})

test_that("two-factor phospho data follows the proteomics rule (not the free-design one)", {
  testthat::skip_if_not_installed("anndataR")
  pr <- make_toy_prot("phospho")
  cd <- pr$coldata
  cd$batch <- rep(c("b1", "b2"), length.out = nrow(cd))
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(
    pr$intensities,
    out = out,
    omics = "proteomics",
    coldata = cd,
    rowdata = pr$rowdata,
    factors = c("group", "batch")
  )
  eds <- new_explore_dataset(out)
  expect_identical(eds$data_type(), "phosphoproteomics")
  bd <- batch_design_defaults(eds)
  expect_identical(bd$batch_default, "batch")
  expect_identical(bd$design_choices, c("None", "group"))
  expect_identical(bd$design_default, "group")
})

test_that("count QC data: batch None by default, any factor can be the design", {
  testthat::skip_if_not_installed("anndataR")
  p <- make_toy_pieces()
  cd <- p$coldata
  cd$batch <- rep(c("b1", "b2"), length.out = nrow(cd))
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(p$counts, out = out, omics = "crispr", coldata = cd, factors = c("condition", "batch"))
  eds <- new_explore_dataset(out)
  expect_identical(eds$data_type(), "CrisprCountQC")
  bd <- batch_design_defaults(eds)
  expect_identical(bd$batch_default, "batch")
  expect_identical(bd$design_choices, c("None", "condition", "batch"))
  expect_identical(bd$design_default, "None")
})

test_that("assay_is_log2: proteomics raw intensities are linear, like genomics counts", {
  testthat::skip_if_not_installed("anndataR")
  pr <- make_toy_prot("protein")
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(
    list(raw = 2^pr$intensities, transformed = pr$intensities),
    out = out,
    omics = "proteomics",
    coldata = pr$coldata,
    rowdata = pr$rowdata
  )
  prot <- new_explore_dataset(out)
  expect_true(assay_is_log2(prot, "intensities"))
  expect_false(assay_is_log2(prot, "intensities_raw"))

  p <- make_toy_pieces()
  out2 <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(p$counts, out = out2, omics = "rnaseq", coldata = p$coldata)
  rna <- new_explore_dataset(out2)
  expect_true(assay_is_log2(rna, "vst"))
  expect_true(assay_is_log2(rna, "normalised_log2"))
  expect_false(assay_is_log2(rna, "counts"))
})
