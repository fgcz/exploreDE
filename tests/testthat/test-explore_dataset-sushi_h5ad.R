# -----------------------------------------------------------------------------
# A SUSHI twoGroups result written as an .h5ad, as ezRun's DESeq2 app does
# (writeAnnData), opens as the SUSHI SummarizedExperiment does: same factor
# and contrast names, sample columns, DE table and pathway terms.
# -----------------------------------------------------------------------------

skip_if_not_installed("anndataR")

# The ezRun call (R/app-deseq2TwoGroups.R), on the demo SUSHI result.
sushi_h5ad <- function(se, env = parent.frame()) {
  param <- S4Vectors::metadata(se)$param
  md <- S4Vectors::metadata(se)
  out <- withr::local_tempfile(fileext = ".h5ad", .local_envir = env)
  build_explore_h5ad(
    x = list(
      raw = SummarizedExperiment::assays(se)$counts,
      normalised = SummarizedExperiment::assays(se)$xNorm
    ),
    out = out,
    omics = "rnaseq",
    feature_level = param$featureLevel,
    ref_build = param$refBuild,
    coldata = SummarizedExperiment::colData(se),
    rowdata = SummarizedExperiment::rowData(se),
    factors = paste(param$groupingName, "[Factor]"),
    de = stats::setNames(list(SummarizedExperiment::rowData(se)), param$comparison),
    pathways = stats::setNames(list(list(ora = md$enrichResult, gsea = md$enrichResultGSEA)), param$comparison),
    enrichment_fdr = list(ora = param$fdrThreshORA, gsea = param$fdrThreshGSEA)
  )
  out
}

se <- readRDS(system.file("examples", "genomics", "deResult.rds", package = "exploreDE"))
from_se <- new_explore_dataset(se)
from_h5ad <- new_explore_dataset(sushi_h5ad(se))

test_that("factor and contrast names are the SUSHI route's", {
  expect_identical(from_h5ad$factor_names(), from_se$factor_names())
  expect_identical(from_h5ad$primary_grouping(), from_se$primary_grouping())
  expect_identical(from_h5ad$contrasts(), from_se$contrasts())
  cm <- from_h5ad$contrast_meta(from_h5ad$primary_contrast())
  expect_identical(
    c(cm$test_group, cm$reference_group),
    unlist(from_se$contrast_meta(from_se$primary_contrast())[c("test_group", "reference_group")], use.names = FALSE)
  )
})

test_that("sample columns lose their SUSHI tags", {
  cols <- colnames(from_h5ad$samples())
  expect_false(any(grepl("\\[", cols)))
  expect_true(all(setdiff(colnames(from_se$samples()), "sf") %in% cols))
})

test_that("GO annotation columns are stored once", {
  cols <- colnames(from_h5ad$features())
  expect_false(anyDuplicated(cols) > 0)
  expect_true(all(c("GO BP", "GO MF", "GO CC") %in% cols))
})

test_that("the DE table has the Type column", {
  c_se <- from_se$primary_contrast()
  c_h5 <- from_h5ad$primary_contrast()
  d_se <- from_se$de(c_se)
  d_h5 <- from_h5ad$de(c_h5)
  expect_true("type" %in% colnames(d_h5))
  expect_identical(d_h5$type[match(d_se$feature_id, d_h5$feature_id)], as.character(d_se$type))
})

test_that("every pathway term is kept, judged at the SUSHI thresholds", {
  for (kind in from_se$enrichment_kinds()) {
    expect_identical(
      nrow(from_h5ad$enrichment_table(from_h5ad$primary_contrast(), kind)),
      nrow(from_se$enrichment_table(from_se$primary_contrast(), kind)),
      info = kind
    )
    expect_identical(.enr_fdr_thresh(from_h5ad, kind), .enr_fdr_thresh(from_se, kind), info = kind)
  }
})
