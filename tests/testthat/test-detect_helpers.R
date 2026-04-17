# --- detectAligner on real fixture ------------------------------------------
test_that("detectAligner identifies the aligner from the real RNASeq fixture's input_dataset.tsv", {
  fx <- fx_rnaseq_gene()
  aligner <- detectAligner(fx$dataDir)
  expect_true(aligner %in% c("STAR", "Kallisto", "FeatureCounts"),
    info = paste("got:", aligner))
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
  expect_equal(detectDEMethod("/srv/gstore/foo/EdgeR_run/bar",  list()), "edgeR")
  expect_equal(detectDEMethod("/srv/gstore/foo/edgeR_run/bar",  list()), "edgeR")
})

test_that("detectDEMethod uses normMethod hint when the path is uninformative", {
  expect_equal(detectDEMethod("/foo/bar", list(normMethod = "TMM")), "edgeR")
  expect_equal(detectDEMethod("/foo/bar", list(normMethod = "RLE")), "edgeR")
})

test_that("detectDEMethod uses testMethod hint when path / normMethod are uninformative", {
  expect_equal(detectDEMethod("/foo/bar", list(testMethod = "glm")),       "edgeR")
  expect_equal(detectDEMethod("/foo/bar", list(testMethod = "exacttest")), "edgeR")
  expect_equal(detectDEMethod("/foo/bar", list(testMethod = "Wald using some model")), "DESeq2")
})

test_that("detectDEMethod defaults to DESeq2 when no signals present", {
  expect_equal(detectDEMethod("/foo/bar", list()), "DESeq2")
})
