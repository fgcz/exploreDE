# -----------------------------------------------------------------------------
# prolfquapp DEA AnnData (.h5ad) -> ExploreDatasetSE_Proteomics
# Reads committed fixtures (tests/testthat/fixtures, made by
# trim_prolfquapp_fixtures.R), so the tests run without prolfquapp installed.
# -----------------------------------------------------------------------------

fixture_210 <- function() testthat::test_path("fixtures", "prolfquapp_2.10_AnnData.h5ad")
fixture_211 <- function() testthat::test_path("fixtures", "prolfquapp_2.11_AnnData.h5ad")

# The SE rebuilt from the 2.10 AnnData, before the legacy-layout step.
se_210 <- function() {
  testthat::skip_if_not_installed("anndataR")
  se_from_prolfquapp_anndata(.h5ad_read(fixture_210()))
}

test_that("a prolfquapp 2.10 h5ad loads as proteomics with annotated contrasts", {
  testthat::skip_if_not_installed("anndataR")
  eds <- new_explore_dataset(fixture_210())

  expect_s3_class(eds, "ExploreDatasetSE_Proteomics")
  expect_identical(eds$omics_class(), "proteomics")
  cons <- grep("^constrast_", names(SummarizedExperiment::rowData(se_210())), value = TRUE)
  expect_setequal(eds$contrasts(), sub("^constrast_", "", cons))
  expect_true(validate_explore_dataset(eds, strict = FALSE)$ok)

  de <- eds$de(eds$primary_contrast())
  expect_true(all(c("description", "log2Ratio", "pValue", "fdr") %in% colnames(de)))
  expect_true(any(!is.na(de$description)))
})

test_that("the rebuilt contrast frame puts the annotation in front of the stats", {
  se_in <- se_210()
  rd_in <- SummarizedExperiment::rowData(se_in)
  rd_out <- SummarizedExperiment::rowData(.prolfquapp_legacy_layout(se_in))
  con <- grep("^constrast_", names(rd_in), value = TRUE)[1]

  annotation <- as.data.frame(rd_in$annotation)
  stats <- as.data.frame(rd_in[[con]])
  rebuilt <- as.data.frame(rd_out[[con]])
  expect_identical(colnames(rebuilt), c(colnames(annotation), setdiff(colnames(stats), colnames(annotation))))
  expect_equal(rebuilt$description, annotation$description)
  expect_equal(rebuilt$diff, stats$diff)
})

test_that("the modelled factor leads colData, so it is the default grouping", {
  se_in <- se_210()
  se <- .prolfquapp_legacy_layout(se_in)
  declared <- names(S4Vectors::metadata(se_in)$analysis_configuration$factors)
  expect_identical(colnames(SummarizedExperiment::colData(se))[seq_along(declared)], declared)
})

test_that("a 2.10.0-2.10.4 SE, contrasts without feature keys, is joined on feature ids", {
  se <- se_210()
  keys <- unlist(S4Vectors::metadata(se)$feature_keys)
  rd <- SummarizedExperiment::rowData(se)
  for (name in grep("^constrast_", names(rd), value = TRUE)) {
    frame <- rd[[name]]
    SummarizedExperiment::rowData(se)[[name]] <- frame[, setdiff(colnames(frame), keys), drop = FALSE]
  }
  out <- as.data.frame(SummarizedExperiment::rowData(.prolfquapp_legacy_layout(se))[[name]])
  expect_true(all(keys %in% colnames(out)))
  expect_equal(out$description, as.data.frame(rd$annotation)$description)
  expect_s3_class(new_explore_dataset(se), "ExploreDatasetSE_Proteomics")
})

test_that("a prolfquapp >= 2.11 h5ad is not read yet (spec section 10)", {
  testthat::skip_if_not_installed("anndataR")
  # Pins the documented refusal; reading this layout will change this test.
  expect_error(new_explore_dataset(fixture_211()), "varm_columns names 0")
})
