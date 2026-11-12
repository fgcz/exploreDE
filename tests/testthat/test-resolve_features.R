# ============================================================================
# resolve_features — expands gene symbols / protein IDs to feature_ids with the
# lookup tables of eds$proteomics_choices().
# ============================================================================

known <- c("P1_pep1", "P1_pep2", "P2_pep1")

test_that("resolve_features with feature_id mode returns input unchanged", {
  expect_equal(resolve_features(c("P1_pep1", "P2_pep1"), "feature_id", list(), known), c("P1_pep1", "P2_pep1"))
})

test_that("resolve_features with gene_name mode expands gene symbols to feature_ids", {
  lookups <- list(gene = list(GENE_A = c("P1_pep1", "P1_pep2"), GENE_B = "P2_pep1"))
  expect_setequal(
    resolve_features(c("GENE_A", "GENE_B"), "gene_name", lookups, known),
    c("P1_pep1", "P1_pep2", "P2_pep1")
  )
})

test_that("resolve_features gene_name mode passes through direct feature_id hits", {
  lookups <- list(gene = list(GENE_A = c("P1_pep1", "P1_pep2")))
  # A gene symbol and a feature id
  expect_setequal(
    resolve_features(c("GENE_A", "P2_pep1"), "gene_name", lookups, known),
    c("P1_pep1", "P1_pep2", "P2_pep1")
  )
})

test_that("resolve_features with protein_Id mode expands via the protein lookup", {
  lookups <- list(protein = list(P1 = c("P1_pep1", "P1_pep2"), P2 = "P2_pep1"))
  expect_setequal(resolve_features("P1", "protein_Id", lookups, known), c("P1_pep1", "P1_pep2"))
})

test_that("resolve_features returns input unchanged when the relevant lookup is absent", {
  expect_equal(resolve_features("GENE_A", "gene_name", list(gene = NULL), known), "GENE_A")
})
