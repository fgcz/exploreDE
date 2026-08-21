# -----------------------------------------------------------------------------
# Conformance test: uns['exploreDE']['mapping'] -- the declared "dictionary"
# WEW proposed (see plans/anndata_manifest_design.md) -- as written by
# build_explore_h5ad(), validates against inst/schema/exploreDE_mapping.schema.json.
#
# Scope: this schema covers ONLY uns$exploreDE$mapping, not the whole h5ad
# contract. It is the dictionary itself -- for THIS dataset, which raw column
# fed each canonical DE target (log2fc/p-value/fdr/extra score), plus the
# feature id/symbol/description column when auto-detection needed an explicit
# override -- kept in the file rather than discarded once used at build time.
# validate_explore_dataset() operates on the live accessor object and never
# inspects this block, so this is the only place the on-disk mapping shape is
# checked against a published, external-facing contract (design doc §4).
# -----------------------------------------------------------------------------

schema_path <- system.file("schema", "exploreDE_mapping.schema.json", package = "exploreDE")

.mapping_validator <- function() {
  testthat::skip_if_not_installed("jsonvalidate")
  if (!nzchar(schema_path)) testthat::skip("exploreDE_mapping.schema.json not found")
  jsonvalidate::json_validator(schema_path, engine = "ajv")
}

# Real pipeline output (an R list, e.g. eds$backend()$uns$exploreDE$mapping) -> JSON -> validate.
validate_mapping <- function(mapping_list) {
  v <- .mapping_validator()
  v(jsonlite::toJSON(mapping_list, auto_unbox = TRUE))
}

# Hand-crafted fixtures as raw JSON text, for testing the schema's own
# strictness. (Sidesteps jsonlite's empty-list {} vs [] ambiguity, which never
# matters for the real writer -- build_explore_h5ad() never emits an empty
# semantics object -- but would make a hand-authored "empty object" case
# ambiguous.)
validate_mapping_json <- function(json_text) {
  v <- .mapping_validator()
  v(json_text)
}

# ---- real pipeline output conforms ------------------------------------------

test_that("default-only DE mapping (generic, non-canonical column names) conforms", {
  testthat::skip_if_not_installed("anndataR")
  p <- make_toy_pieces()                                 # de_df has LFC / P / Q
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata,
                     rowdata = p$rowdata, de = list(treated_vs_control = p$de_df),
                     de_columns = list(effect_log2fc = "LFC", score_pvalue = "P",
                                       score_fdr = "Q", score = "P",
                                       score_kind = "p_value_low_sig",
                                       fdr_kind = "fdr_low_sig"))

  eds <- new_explore_dataset(out)
  mp <- eds$backend()$uns$exploreDE$mapping
  expect_identical(mp$de$default$score_label, "P")
  expect_identical(mp$de$default$effect_log2fc_col, "LFC")
  expect_identical(mp$de$default$score_pvalue_col, "P")
  expect_identical(mp$de$default$score_fdr_col, "Q")
  expect_identical(mp$de$default$score_col, "P")
  expect_true(validate_mapping(mp))
})

test_that("per-contrast override mapping (no shared default) conforms", {
  testthat::skip_if_not_installed("anndataR")
  n <- 12L; ns <- 8L
  m <- matrix(seq_len(n * ns) + 0.5, nrow = n, ncol = ns)
  rownames(m) <- sprintf("F%02d", seq_len(n)); colnames(m) <- sprintf("S%02d", seq_len(ns))
  coldata <- data.frame(row.names = colnames(m),
                        group = rep(c("A", "B", "C", "ctrl"), each = ns / 4L))
  mk_de <- function(seed) data.frame(row.names = rownames(m),
                                     lfc = seq(-2, 2, length.out = n) + seed,
                                     fdrcol = seq(0, 0.9, length.out = n),
                                     extra = seq(0, 1, length.out = n))
  de <- list(A_vs_ctrl = mk_de(0), B_vs_ctrl = mk_de(0.2))
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(m, out = out, omics = "rnaseq", coldata = coldata, de = de,
                     de_columns = list(
                       A_vs_ctrl = list(effect_log2fc = "lfc", score_fdr = "fdrcol",
                                        fdr_kind = "fdr_low_sig", score = "extra",
                                        score_kind = "p_value_low_sig"),
                       B_vs_ctrl = list(effect_log2fc = "lfc", score_fdr = "fdrcol",
                                        fdr_kind = "fdr_low_sig", score = "extra",
                                        score_kind = "probability_high_sig")))

  eds <- new_explore_dataset(out)
  mp <- eds$backend()$uns$exploreDE$mapping
  expect_null(mp$de$default)
  expect_length(mp$de$contrasts, 2L)
  expect_identical(mp$de$contrasts$a_vs_ctrl$score_kind, "p_value_low_sig")
  expect_identical(mp$de$contrasts$b_vs_ctrl$score_kind, "probability_high_sig")
  expect_true(validate_mapping(mp))
})

test_that("no declared mapping -> uns$exploreDE is absent (canonical files stay lean)", {
  testthat::skip_if_not_installed("anndataR")
  p <- make_toy_pieces()
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata,
                     rowdata = p$rowdata, de = list(treated_vs_control = p$de_df))

  eds <- new_explore_dataset(out)
  expect_null(eds$backend()$uns$exploreDE)
})

test_that("explicit rowdata_columns override -> mapping$feature conforms, and actually drives the feature id (not just a label)", {
  testthat::skip_if_not_installed("anndataR")
  n <- 6L; ns <- 4L
  m <- matrix(seq_len(n * ns) + 0.5, nrow = n, ncol = ns)
  rownames(m) <- sprintf("F%02d", seq_len(n)); colnames(m) <- sprintf("S%02d", seq_len(ns))
  coldata <- data.frame(row.names = colnames(m), group = rep(c("A", "B"), each = ns / 2L))
  # Column names a foreign proteomics pipeline might use instead of the
  # built-in protein_Id/gene_name detection.
  rowdata <- data.frame(row.names = rownames(m),
                        UniProt_Accession = sprintf("P%05d", seq_len(n)),
                        Display_Name = sprintf("GENE%d", seq_len(n)),
                        stringsAsFactors = FALSE)
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(m, out = out, omics = "proteomics", coldata = coldata, rowdata = rowdata,
                     rowdata_columns = list(id_col = "UniProt_Accession", symbol_col = "Display_Name"))

  eds <- new_explore_dataset(out)
  mp <- eds$backend()$uns$exploreDE$mapping
  expect_identical(mp$feature$id_col, "UniProt_Accession")
  expect_identical(mp$feature$symbol_col, "Display_Name")
  expect_null(mp$de)
  expect_true(validate_mapping(mp))
  # Not just persisted as a label — it drove the actual protein-level feature_id.
  expect_setequal(eds$feature_ids(), sprintf("P%05d", seq_len(n)))
})

# ---- the schema itself is actually strict, not a rubber stamp --------------

test_that("schema accepts the SAINT-shaped default mapping directly", {
  ok <- validate_mapping_json(
    '{"de":{"default":{"score_label":"SaintScore","score_kind":"probability_high_sig","fdr_kind":"fdr_low_sig"}}}')
  expect_true(ok)
})

test_that("schema accepts a feature mapping alongside de", {
  ok <- validate_mapping_json(
    '{"de":{"default":{"score_label":"x"}},"feature":{"id_col":"UniProt_Accession"}}')
  expect_true(ok)
})

test_that("schema rejects a bad score_kind enum value", {
  ok <- validate_mapping_json('{"de":{"default":{"score_kind":"not_a_real_kind"}}}')
  expect_false(ok)
})

test_that("schema rejects an undeclared top-level key (e.g. a not-yet-built 'assay' block)", {
  ok <- validate_mapping_json(
    '{"de":{"default":{"score_label":"x"}},"assay":{"counts_layer":"raw"}}')
  expect_false(ok)
})

test_that("schema rejects an empty semantics object", {
  ok <- validate_mapping_json('{"de":{"default":{}}}')
  expect_false(ok)
})
