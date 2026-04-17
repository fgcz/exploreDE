# --- generateParams ---------------------------------------------------------
test_that("generateParams returns a list with all parameter defaults", {
  p <- generateParams()
  expect_type(p, "list")
  expect_true(all(c("sampleGroup", "refGroup", "groupingName", "comparison",
                    "featureLevel", "runGO", "fdrThreshORA", "fdrThreshGSEA",
                    "pValThreshGO", "log2RatioThreshGO", "backgroundExpression",
                    "refBuild") %in% names(p)))
  expect_equal(p$sampleGroup, "diff")
  expect_equal(p$refGroup, "undiff")
  expect_equal(p$featureLevel, "gene")
  expect_true(p$runGO)
})

test_that("generateParams applies overrides passed as named arguments", {
  p <- generateParams(sampleGroup = "tumour", refGroup = "normal",
                      featureLevel = "isoform", runGO = FALSE)
  expect_equal(p$sampleGroup, "tumour")
  expect_equal(p$refGroup, "normal")
  expect_equal(p$featureLevel, "isoform")
  expect_false(p$runGO)
})


# --- transform_feature_id ---------------------------------------------------
test_that("transform_feature_id replaces special characters with '.'", {
  expect_equal(transform_feature_id("ENST00000001:001"),  "ENST00000001.001")
  expect_equal(transform_feature_id("gene(X)-var"),      "gene.X..var")
  expect_equal(transform_feature_id("foo bar"),           "foo.bar")
  expect_equal(transform_feature_id("a/b/c"),             "a.b.c")
})

test_that("transform_feature_id leaves plain gene-symbol-like strings unchanged", {
  expect_equal(transform_feature_id("BRCA1"),       "BRCA1")
  expect_equal(transform_feature_id("ENSG0000001"), "ENSG0000001")
})

test_that("transform_feature_id is vectorised", {
  expect_equal(
    transform_feature_id(c("A:1", "B/2", "C")),
    c("A.1", "B.2", "C")
  )
})


# --- align_signal_to_seqanno ------------------------------------------------
test_that("align_signal_to_seqanno reorders a signal matrix to match seqAnno rows", {
  signal <- matrix(1:12, nrow = 3,
                   dimnames = list(c("G1", "G2", "G3"), paste0("s", 1:4)))
  rd <- data.frame(feature_id = c("G3", "G1"))

  out <- align_signal_to_seqanno(signal, rd)
  expect_equal(rownames(out$signal), c("G3", "G1"))
  expect_equal(out$rd$feature_id, c("G3", "G1"))
})

test_that("align_signal_to_seqanno drops rd rows that don't match any signal row", {
  signal <- matrix(1:8, nrow = 2, dimnames = list(c("G1", "G2"), paste0("s", 1:4)))
  rd <- data.frame(feature_id = c("G1", "ORPHAN", "G2"))

  out <- align_signal_to_seqanno(signal, rd)
  expect_setequal(rownames(out$signal), c("G1", "G2"))
  expect_setequal(out$rd$feature_id, c("G1", "G2"))
})

test_that("align_signal_to_seqanno falls back to isoform-style substitutions when needed", {
  # Signal rownames are already in the transformed form (e.g. exported by the app)
  # "ENST00000001:001" → "ENST00000001.001"; "gene(X)" → "gene.X."
  signal <- matrix(1:4, nrow = 2,
                   dimnames = list(c("ENST00000001.001", "gene.X."), paste0("s", 1:2)))
  rd <- data.frame(feature_id = c("ENST00000001:001", "gene(X)"))

  out <- align_signal_to_seqanno(signal, rd)
  # Fallback matches both
  expect_equal(nrow(out$signal), 2)
  expect_equal(nrow(out$rd), 2)
})
