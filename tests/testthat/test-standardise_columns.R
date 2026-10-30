test_that("standardise_columns maps canonical names for common proteomics columns", {
  df <- data.frame(
    "Protein.Id" = c("P1", "P2"),
    "Peptide" = c("AAA", "BBB"),
    "site" = c("S1", "S2"),
    "gene_symbol" = c("G1", "G2"),
    "description" = c("d1", "d2"),
    "nrPeptides" = c(3, 4),
    "diff" = c(1.1, -0.8),
    "p.value" = c(0.01, 0.5),
    "FDR" = c(0.02, 0.6),
    check.names = FALSE
  )
  out <- standardise_columns(df)
  expect_setequal(
    colnames(out),
    c("protein_Id", "peptide_Id", "site", "gene_name", "description", "nrPeptides", "log2Ratio", "pValue", "fdr")
  )
})

test_that("standardise_columns preserves unknown columns unchanged", {
  df <- data.frame(
    Protein_ID = "P1",
    mystery_col = "x",
    check.names = FALSE
  )
  out <- standardise_columns(df)
  expect_true("mystery_col" %in% colnames(out))
  expect_true("protein_Id" %in% colnames(out))
})

test_that("standardise_columns does not rename when pattern matches multiple columns", {
  # If two columns both match a canonical pattern, neither should be renamed
  # (the function only acts when exactly one column matches).
  df <- data.frame(
    p_value = 0.01,
    pValue = 0.02,
    check.names = FALSE
  )
  out <- standardise_columns(df)
  # Both original names preserved; neither collapsed
  expect_true(all(c("p_value", "pValue") %in% colnames(out)))
})

test_that("standardise_columns is case-insensitive", {
  df <- data.frame(
    PROTEIN_ID = "P1",
    GENE_NAME = "G1",
    PADJ = 0.01,
    check.names = FALSE
  )
  out <- standardise_columns(df)
  expect_true("protein_Id" %in% colnames(out))
  expect_true("gene_name" %in% colnames(out))
  expect_true("fdr" %in% colnames(out))
})

test_that("standardise_columns accepts log2FC variants", {
  for (variant in c("log2FC", "log2_fold_change", "logFC", "log2.ratio")) {
    df <- data.frame(x = variant, check.names = FALSE)
    colnames(df) <- variant
    out <- standardise_columns(df)
    expect_true(
      "log2Ratio" %in% colnames(out),
      info = paste("variant that failed to map:", variant)
    )
  }
})
