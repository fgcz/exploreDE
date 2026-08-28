# -----------------------------------------------------------------------------
# Conformance test: inst/dictionaries/column_patterns.json -- the AUTOMATIC
# half of exploreDE's dictionary (plans/anndata_manifest_design.md). When no
# explicit de_columns/rowdata_columns mapping is given, THIS is what figures
# out which column is which for an arbitrary incoming data.frame. Lives as
# data, validated against inst/schema/exploreDE_column_patterns.schema.json,
# so widening coverage to a new pipeline's naming convention is a JSON edit,
# not an R code change.
# -----------------------------------------------------------------------------

patterns_path <- system.file("dictionaries", "column_patterns.json", package = "exploreDE")
schema_path   <- system.file("schema", "exploreDE_column_patterns.schema.json", package = "exploreDE")

.patterns_validator <- function() {
  testthat::skip_if_not_installed("jsonvalidate")
  if (!nzchar(schema_path)) testthat::skip("exploreDE_column_patterns.schema.json not found")
  jsonvalidate::json_validator(schema_path, engine = "ajv")
}

test_that("column_patterns.json conforms to its own schema", {
  if (!nzchar(patterns_path)) testthat::skip("column_patterns.json not found")
  v <- .patterns_validator()
  expect_true(v(patterns_path))
})

test_that("eds_load_column_patterns() returns non-empty de/rowdata groups", {
  loaded <- eds_load_column_patterns()
  expect_true(length(loaded$de) > 0)
  expect_true(length(loaded$rowdata) > 0)
  expect_true(all(c("effect_log2fc", "score_pvalue", "score_fdr") %in% names(loaded$de)))
})

test_that("eds_de_column_patterns / eds_rowdata_column_patterns are populated from the JSON file at load time", {
  expect_identical(eds_de_column_patterns, eds_load_column_patterns()$de)
  expect_identical(eds_rowdata_column_patterns, eds_load_column_patterns()$rowdata)
})

# ---- concrete before/after wins: names the OLD hardcoded regex missed ------

test_that("limma-style 'adj.P.Val' is now recognised as score_fdr (previously unmatched: needed literal 'value')", {
  d <- data.frame(logFC = 1, `adj.P.Val` = 0.03, check.names = FALSE)
  out <- eds_standardise_de_columns(d)
  expect_identical(out$score_fdr, 0.03)
})

test_that("'qval' (no separator, short form) is now recognised as score_fdr (previously required literal 'value')", {
  d <- data.frame(logFC = 1, qval = 0.02)
  out <- eds_standardise_de_columns(d)
  expect_identical(out$score_fdr, 0.02)
})

test_that("standard DESeq2 names still resolve unchanged", {
  d <- data.frame(log2FoldChange = 1.5, pvalue = 0.01, padj = 0.04)
  out <- eds_standardise_de_columns(d)
  expect_identical(out$effect_log2fc, 1.5)
  expect_identical(out$score_pvalue, 0.01)
  expect_identical(out$score_fdr, 0.04)
})

test_that("SAINTexpress's log2_EFCs/SaintScore/BFDR now auto-resolve with ZERO explicit mapping", {
  # The flagship motivating case (plans/anndata_manifest_design.md) -- these
  # three previously needed a hand-written de_columns= override because the
  # old patterns were fully-anchored (^fdr$ requires the name to be EXACTLY
  # "fdr"). BFDR/log2_EFCs/SaintScore all carry the concept as a SUBSTRING,
  # which the broadened patterns now catch.
  d <- data.frame(log2_EFCs = 1.5, SaintScore = 0.9, BFDR = 0.03)
  out <- eds_standardise_de_columns(d)
  expect_identical(out$effect_log2fc, 1.5)
  expect_identical(out$score_fdr, 0.03)
  expect_identical(out$score, 0.9)
  expect_false("score_pvalue" %in% colnames(out))   # SAINT genuinely has no p-value
})

test_that("a truly opaque column name (no textual signal at all) is left alone, not guessed", {
  d <- data.frame(xq7_metric = 1, zzflag = 0.9, weird_col = 0.05)
  out <- eds_standardise_de_columns(d)
  expect_false("effect_log2fc" %in% colnames(out))
  expect_false("score_pvalue" %in% colnames(out))
  expect_false("score_fdr" %in% colnames(out))
  expect_false("score" %in% colnames(out))
})

# ---- the schema itself is actually strict ----------------------------------

test_that("schema rejects a pattern group with an empty regex array", {
  v <- .patterns_validator()
  expect_false(v('{"de":{"score_fdr":[]},"rowdata":{"protein_Id":["^protein_id$"]}}'))
})

test_that("schema rejects an unknown top-level key", {
  v <- .patterns_validator()
  expect_false(v('{"de":{"score_fdr":["^fdr$"]},"rowdata":{"protein_Id":["^protein_id$"]},"assay":{}}'))
})

test_that("schema requires both de and rowdata groups", {
  v <- .patterns_validator()
  expect_false(v('{"de":{"score_fdr":["^fdr$"]}}'))
})

# ---- laziness (exploreDElight / wasm) ---------------------------------------

test_that("sourcing utils_explore_dataset.R does not require the dictionary", {
  # Regression: the pattern objects used to be assigned eagerly, so merely
  # source()ing this file demanded the JSON on disk -- and when it was absent
  # the whole file failed to load, taking eds_flatten_cp_results() and
  # .enr_tidy_from_native() with it. exploreDElight sources these files into
  # webR, where there is no installed package and no inst/ tree, so enrichment
  # broke over a dictionary Light never consults.
  #
  # Simulated by shadowing system.file() inside the evaluation environment, so
  # the loader cannot find the JSON by any route.
  src <- testthat::test_path("..", "..", "R", "utils_explore_dataset.R")
  skip_if_not(file.exists(src))
  # Resolve BEFORE with_dir(): arguments are lazy, so a relative path would be
  # normalised against the temp dir and fail to resolve.
  src <- normalizePath(src)

  e <- new.env(parent = globalenv())
  e$system.file <- function(...) ""
  tmp <- withr::local_tempdir()

  # Both the source AND the force happen in the temp dir. Not tidiness: the
  # promise is forced later than it is created, and the loader's last-resort
  # fallback is getwd()-relative -- so forcing back in the package root would
  # find the real inst/dictionaries and the test would pass vacuously. (That
  # this matters at all is itself evidence the load really is deferred.)
  withr::with_dir(tmp, {
    sys.source(src, envir = e, keep.source = FALSE)

    # Guard: if sourcing silently failed, force(NULL) would not error and the
    # dictionary assertion below would pass for the wrong reason.
    expect_true(exists("eds_de_column_patterns", envir = e, inherits = FALSE))

    # The helpers that used to go down with the file are present ...
    expect_true(is.function(e$eds_flatten_cp_results))
    expect_true(is.function(e$.enr_tidy_from_native))
    expect_identical(e$eds_safe("a b/c"), "a_b_c")

    # ... and the dictionary still fails LOUDLY when something actually asks,
    # rather than degrading to an empty pattern set that would silently stop
    # recognising column names.
    expect_error(force(e$eds_de_column_patterns), "column_patterns")
  })
})
