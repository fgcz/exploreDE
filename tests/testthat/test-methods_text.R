# ============================================================================
# buildCitationsList
# ============================================================================
test_that("buildCitationsList always includes SUSHI / FastQC / MultiQC / fastp / R / exploreDE", {
  cits <- buildCitationsList(aligner = "STAR", de_method = "DESeq2")
  expect_true(all(c("SUSHI", "FastQC", "MultiQC", "fastp", "R", "exploreDE") %in% names(cits)))
})

test_that("buildCitationsList adds aligner-specific citations", {
  expect_true("STAR" %in% names(buildCitationsList("STAR", "DESeq2")))
  expect_true("Kallisto" %in% names(buildCitationsList("Kallisto", "DESeq2")))
  # FeatureCounts pulls in both Rsubread and STAR
  fc <- buildCitationsList("FeatureCounts", "DESeq2")
  expect_true(all(c("Rsubread", "STAR") %in% names(fc)))
})

test_that("buildCitationsList adds DESeq2 / edgeR citations per DE method", {
  expect_true("DESeq2" %in% names(buildCitationsList("STAR", "DESeq2")))
  expect_true("edgeR" %in% names(buildCitationsList("STAR", "edgeR")))
})

test_that("buildCitationsList adds clusterProfiler only when run_go = TRUE", {
  with_go <- buildCitationsList("STAR", "DESeq2", run_go = TRUE)
  without_go <- buildCitationsList("STAR", "DESeq2", run_go = FALSE)
  expect_true("clusterProfiler" %in% names(with_go))
  expect_false("clusterProfiler" %in% names(without_go))
})

test_that("buildCitationsList adds GENCODE citation only when ref_build references GENCODE", {
  with_g <- buildCitationsList("STAR", "DESeq2", ref_build = "Homo_sapiens/GENCODE/GRCh38")
  without_g <- buildCitationsList("STAR", "DESeq2", ref_build = "Homo_sapiens/Ensembl/GRCh38")
  expect_true("GENCODE" %in% names(with_g))
  expect_false("GENCODE" %in% names(without_g))
})


# ============================================================================
# getFastpVersion
# ============================================================================
test_that("getFastpVersion returns NULL when no *_preprocessing.log exists", {
  td <- tempfile("fastp_empty_")
  dir.create(file.path(td, "sub"), recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE))
  expect_null(getFastpVersion(file.path(td, "sub")))
})

test_that("getFastpVersion extracts the version from a preprocessing log", {
  td <- tempfile("fastp_ok_")
  dir.create(file.path(td, "sub"), recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE))
  writeLines(
    c("Running fastp v0.23.4 ...", "some other output"),
    file.path(td, "sample_preprocessing.log")
  )
  expect_equal(getFastpVersion(file.path(td, "sub")), "0.23.4")
})


# ============================================================================
# generateMethodsText — end-to-end on a synthetic param bundle
# ============================================================================
test_that("generateMethodsText produces the expected sections for RNASeq/STAR/DESeq2", {
  param <- list(
    refBuild = "Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30",
    featureLevel = "gene",
    transcriptTypes = "protein-coding",
    runGO = TRUE,
    fdrThreshORA = 0.05,
    fdrThreshGSEA = 0.05
  )
  pkg_versions <- c(DESeq2 = "1.42.0", edgeR = "version unavailable", clusterProfiler = "4.10.0")

  sections <- generateMethodsText(
    param = param,
    aligner = "STAR",
    de_method = "DESeq2",
    fastp_version = "0.23.4",
    pkg_versions = pkg_versions,
    se = NULL,
    raw_reads_text = ""
  )
  expect_true(all(c("trimming", "mapping", "de_analysis", "go_analysis", "data_processing") %in% names(sections)))
  expect_match(sections$trimming, "fastp.*0\\.23\\.4")
  expect_match(sections$mapping, "STAR")
  expect_match(sections$mapping, "GENCODE human genome build GRCh38")
  expect_match(sections$de_analysis, "DESeq2")
  expect_match(sections$go_analysis, "clusterProfiler")
  expect_match(sections$data_processing, "SUSHI")
})

test_that("generateMethodsText omits go_analysis when runGO is FALSE", {
  param <- list(refBuild = "GRCh38", runGO = FALSE)
  sections <- generateMethodsText(
    param = param,
    aligner = "STAR",
    de_method = "DESeq2",
    fastp_version = NULL,
    pkg_versions = c(DESeq2 = "1.0", edgeR = "1.0", clusterProfiler = "1.0"),
    se = NULL
  )
  expect_false("go_analysis" %in% names(sections))
})

test_that("generateMethodsText picks the Kallisto / edgeR wording when configured that way", {
  param <- list(refBuild = "Mus_musculus/GENCODE/GRCm39", runGO = FALSE)
  sections <- generateMethodsText(
    param = param,
    aligner = "Kallisto",
    de_method = "edgeR",
    fastp_version = NULL,
    pkg_versions = c(DESeq2 = "version unavailable", edgeR = "4.0.0"),
    se = NULL
  )
  expect_match(sections$mapping, "Kallisto")
  expect_match(sections$de_analysis, "edgeR")
  expect_match(sections$de_analysis, "TMM")
})
