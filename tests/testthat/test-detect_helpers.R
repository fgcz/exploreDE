# --- detectAligner on real fixture ------------------------------------------
test_that("detectAligner identifies the aligner from the real RNASeq fixture's input_dataset.tsv", {
  fx <- fx_rnaseq_gene_gstore()
  aligner <- detectAligner(fx$dataDir)
  expect_true(aligner %in% c("STAR", "Kallisto", "FeatureCounts"), info = paste("got:", aligner))
})

test_that("detectAligner returns 'Unknown' when input_dataset.tsv is missing", {
  td <- tempfile("detectAligner_empty_")
  dir.create(file.path(td, "sub"), recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE))
  # Passing a path whose ../input_dataset.tsv does not exist
  expect_equal(detectAligner(file.path(td, "sub")), "Unknown")
})

test_that("detectAligner picks up Kallisto / FeatureCounts from Count File column", {
  td <- tempfile("detectAligner_")
  dir.create(file.path(td, "sub"), recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE))
  # Kallisto
  readr::write_tsv(
    data.frame(Name = "s1", `Count [File]` = "Kallisto/s1.h5", check.names = FALSE),
    file.path(td, "input_dataset.tsv")
  )
  expect_equal(detectAligner(file.path(td, "sub")), "Kallisto")

  # FeatureCounts
  readr::write_tsv(
    data.frame(Name = "s1", `Count [File]` = "FeatureCounts/s1.counts", check.names = FALSE),
    file.path(td, "input_dataset.tsv")
  )
  expect_equal(detectAligner(file.path(td, "sub")), "FeatureCounts")

  # STAR — falls back to BAM column presence
  readr::write_tsv(
    data.frame(Name = "s1", `BAM [File]` = "bams/s1.bam", check.names = FALSE),
    file.path(td, "input_dataset.tsv")
  )
  expect_equal(detectAligner(file.path(td, "sub")), "STAR")
})


# --- detectDEMethod ---------------------------------------------------------
test_that("detectDEMethod picks DESeq2 / edgeR from the directory path", {
  expect_equal(detectDEMethod("/srv/gstore/foo/DESeq2_run/bar", list()), "DESeq2")
  expect_equal(detectDEMethod("/srv/gstore/foo/EdgeR_run/bar", list()), "edgeR")
  expect_equal(detectDEMethod("/srv/gstore/foo/edgeR_run/bar", list()), "edgeR")
})

test_that("detectDEMethod uses normMethod hint when the path is uninformative", {
  expect_equal(detectDEMethod("/foo/bar", list(normMethod = "TMM")), "edgeR")
  expect_equal(detectDEMethod("/foo/bar", list(normMethod = "RLE")), "edgeR")
})

test_that("detectDEMethod uses testMethod hint when path / normMethod are uninformative", {
  expect_equal(detectDEMethod("/foo/bar", list(testMethod = "glm")), "edgeR")
  expect_equal(detectDEMethod("/foo/bar", list(testMethod = "exacttest")), "edgeR")
  expect_equal(detectDEMethod("/foo/bar", list(testMethod = "Wald using some model")), "DESeq2")
})

test_that("detectDEMethod defaults to DESeq2 when no signals present", {
  expect_equal(detectDEMethod("/foo/bar", list()), "DESeq2")
})

test_that("count_result_summary gives the settings table of the shipped SUSHI example", {
  se <- readRDS(system.file("examples", "genomics", "deResult.rds", package = "exploreDE"))
  ours <- count_result_summary(S4Vectors::metadata(se)$param, se)
  # count_result_summary() ports ezRun::makeCountResultSummary(). ezRun is not
  # on CRAN or Bioconductor, so it cannot be a test dependency: the expected
  # table is pinned from exploreDE 3.10.0 instead of compared live.
  expected <- data.frame(
    settings = c(
      "NGS two group analysis",
      "Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30",
      "gene",
      "est_counts",
      "deseq2",
      "diff--over--undiff",
      "DESeq2_MedianRatio",
      "3.322",
      "10"
    ),
    row.names = c(
      "Analysis:",
      "Reference:",
      "Feature level:",
      "Data Column Used:",
      "Method:",
      "Comparison:",
      "Normalization:",
      "Log2 signal threshold:",
      "Linear signal threshold:"
    )
  )
  expect_identical(ours, expected)
})
