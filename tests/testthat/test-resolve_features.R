# ============================================================================
# resolve_features — expands gene symbols / protein IDs to feature_ids based
# on lookup tables attached to the inputDataReactive() return value.
# ============================================================================

make_resolve_mocks <- function(feature_ids = c("P1_pep1", "P1_pep2", "P2_pep1"),
                               gene_lookup = NULL, protein_lookup = NULL) {
  list(
    seqAnnoReactive   = shiny::reactiveValues(sa = data.frame(feature_id = feature_ids)),
    inputDataReactive = function() list(gene_lookup = gene_lookup,
                                         protein_lookup = protein_lookup)
  )
}

test_that("resolve_features with feature_id mode returns input unchanged", {
  m <- make_resolve_mocks()
  fn <- with_mock_env(resolve_features, !!!m)

  out <- shiny::isolate(fn(c("P1_pep1", "P2_pep1")))
  expect_equal(out, c("P1_pep1", "P2_pep1"))
})

test_that("resolve_features with gene_name mode expands gene symbols to feature_ids", {
  gene_lookup <- list(GENE_A = c("P1_pep1", "P1_pep2"),
                      GENE_B = c("P2_pep1"))
  m <- make_resolve_mocks(gene_lookup = gene_lookup)
  fn <- with_mock_env(resolve_features, !!!m)

  out <- shiny::isolate(fn(c("GENE_A", "GENE_B"), search_mode = "gene_name"))
  expect_setequal(out, c("P1_pep1", "P1_pep2", "P2_pep1"))
})

test_that("resolve_features gene_name mode passes through direct feature_id hits", {
  gene_lookup <- list(GENE_A = c("P1_pep1", "P1_pep2"))
  m <- make_resolve_mocks(gene_lookup = gene_lookup)
  fn <- with_mock_env(resolve_features, !!!m)

  # User entered a mix: gene symbol + a direct feature_id
  out <- shiny::isolate(fn(c("GENE_A", "P2_pep1"), search_mode = "gene_name"))
  expect_setequal(out, c("P1_pep1", "P1_pep2", "P2_pep1"))
})

test_that("resolve_features with protein_Id mode expands via protein_lookup", {
  protein_lookup <- list(P1 = c("P1_pep1", "P1_pep2"),
                         P2 = c("P2_pep1"))
  m <- make_resolve_mocks(protein_lookup = protein_lookup)
  fn <- with_mock_env(resolve_features, !!!m)

  out <- shiny::isolate(fn("P1", search_mode = "protein_Id"))
  expect_setequal(out, c("P1_pep1", "P1_pep2"))
})

test_that("resolve_features returns input unchanged when the relevant lookup is absent", {
  # gene_name mode but no gene_lookup → pass-through
  m <- make_resolve_mocks(gene_lookup = NULL)
  fn <- with_mock_env(resolve_features, !!!m)
  out <- shiny::isolate(fn("GENE_A", search_mode = "gene_name"))
  expect_equal(out, "GENE_A")
})
