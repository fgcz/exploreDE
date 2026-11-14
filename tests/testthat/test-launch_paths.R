# -----------------------------------------------------------------------------
# exploreDE(file) with a live object: .launch_paths() writes it to a temporary
# file that the app opens like a path.
# -----------------------------------------------------------------------------

test_that("a SummarizedExperiment becomes an .rds the app opens as that dataset", {
  se <- readRDS(system.file("examples", "proteomics", "SummarizedExperiment.rds", package = "exploreDE"))
  path <- .launch_paths(se)
  expect_match(path, "\\.rds$")
  expect_s3_class(new_explore_dataset(path), "ExploreDatasetSE_Proteomics")
})

test_that("an AnnData object becomes an .h5ad the app opens", {
  skip_if_not_installed("anndataR")
  p <- make_toy_pieces()
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(p$counts, out = h5ad, omics = "rnaseq", coldata = p$coldata)
  adata <- anndataR::read_h5ad(h5ad, as = "InMemoryAnnData")
  path <- .launch_paths(adata)
  expect_match(path, "\\.h5ad$")
  expect_s3_class(new_explore_dataset(path), "ExploreDatasetH5AD")
})

test_that("paths pass through, normalised, and must exist", {
  h5ad <- withr::local_tempfile(fileext = ".h5ad")
  writeLines("x", h5ad)
  expect_identical(.launch_paths(h5ad), normalizePath(h5ad))
  expect_error(.launch_paths(file.path(tempdir(), "does-not-exist.h5ad")))
})
