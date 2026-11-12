# RNA-seq and Count QC SummarizedExperiments read factors, features and
# samples the same way (shared genomics helpers).

test_that("Count QC reads factors, features and samples as RNA-seq does", {
  de_result <- readRDS(system.file("examples", "genomics", "deResult.rds", package = "exploreDE"))
  rna <- ExploreDatasetSE_RNASeq$new(de_result)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = SummarizedExperiment::assay(de_result, "counts")),
    colData = SummarizedExperiment::colData(de_result),
    rowData = SummarizedExperiment::rowData(de_result)[, c("gene_id", "gene_name", "description")]
  )
  qc <- ExploreDatasetSE_CountQC$new(se)

  cols <- c("feature_id", "display_name", "description")
  f_rna <- rna$features()[, cols]
  expect_identical(qc$features()[rownames(f_rna), cols], f_rna)
  expect_false(any(grepl("[Source:", qc$features()$description, fixed = TRUE)))
  expect_identical(qc$factor_names(), rna$factor_names())
  for (f in qc$factor_names()) {
    expect_identical(qc$factor_levels(f), rna$factor_levels(f))
  }
  expect_identical(qc$samples(), rna$samples())
  expect_identical(rownames(qc$assay("counts")), rownames(rna$assay("counts")))
})
