# -----------------------------------------------------------------------------
# prolfquapp DEA AnnData (.h5ad) -> ExploreDatasetSE_Proteomics
# Needs prolfquapp to write a DEA AnnData; skipped where it is not installed.
# -----------------------------------------------------------------------------

prolfquapp_h5ad <- function() {
  testthat::skip_if_not_installed("prolfquapp")
  dea <- prolfquapp::example_deanalyse(Nprot = 30)
  se <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config, name = "")$make_SummarizedExperiment()
  path <- file.path(withr::local_tempdir(.local_envir = parent.frame()), "AnnData.h5ad")
  prolfquapp:::write_summarized_experiment_h5ad(se, path)
  list(se = se, path = path)
}

test_that("a prolfquapp h5ad loads as proteomics with annotated contrasts", {
  fx <- prolfquapp_h5ad()
  eds <- new_explore_dataset(fx$path)

  expect_s3_class(eds, "ExploreDatasetSE_Proteomics")
  expect_identical(eds$omics_class(), "proteomics")
  cons <- grep("^constrast_", names(SummarizedExperiment::rowData(fx$se)), value = TRUE)
  expect_setequal(eds$contrasts(), sub("^constrast_", "", cons))
  expect_true(validate_explore_dataset(eds, strict = FALSE)$ok)

  de <- eds$de(eds$primary_contrast())
  expect_true(all(c("description", "log2Ratio", "pValue", "fdr") %in% colnames(de)))
  expect_true(any(!is.na(de$description)))
})

test_that("the rebuilt contrast frame puts the annotation in front of the stats", {
  fx <- prolfquapp_h5ad()
  se <- .prolfquapp_legacy_layout(se_from_prolfquapp_anndata(.h5ad_read(fx$path)))
  rd_in <- SummarizedExperiment::rowData(fx$se)
  rd_out <- SummarizedExperiment::rowData(se)
  con <- grep("^constrast_", names(rd_in), value = TRUE)[1]

  annotation <- as.data.frame(rd_in$annotation)
  stats <- as.data.frame(rd_in[[con]])
  rebuilt <- as.data.frame(rd_out[[con]])
  expect_identical(colnames(rebuilt), c(colnames(annotation), setdiff(colnames(stats), colnames(annotation))))
  expect_equal(rebuilt$description, annotation$description)
  expect_equal(rebuilt$diff, stats$diff)
})

test_that("the modelled factor leads colData, so it is the default grouping", {
  fx <- prolfquapp_h5ad()
  se <- .prolfquapp_legacy_layout(se_from_prolfquapp_anndata(.h5ad_read(fx$path)))
  declared <- names(S4Vectors::metadata(fx$se)$analysis_configuration$factors)
  expect_identical(colnames(SummarizedExperiment::colData(se))[seq_along(declared)], declared)
})

test_that("a 2.10.0-2.10.4 SE, contrasts without feature keys, is joined on feature ids", {
  fx <- prolfquapp_h5ad()
  se <- fx$se
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
