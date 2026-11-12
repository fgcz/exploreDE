test_that("kegg_species_label maps the three supported organisms", {
  expect_identical(kegg_species_label("Homo_sapiens/GENCODE38/Genes"), "Human")
  expect_identical(kegg_species_label("Mus_musculus/GENCODE_M29/Genes"), "Mouse")
  expect_identical(kegg_species_label("Rattus_norvegicus/Ensembl/Genes"), "Rat")
})

test_that("kegg_species_label matches on the leading Genus_species segment alone", {
  expect_identical(kegg_species_label("Homo_sapiens"), "Human")
})

test_that("kegg_species_label returns NULL for unsupported organisms", {
  expect_null(kegg_species_label("Danio_rerio/GRCz11/Genes"))
})

test_that("kegg_species_label returns NULL for NA/NULL/empty input", {
  expect_null(kegg_species_label(NA_character_))
  expect_null(kegg_species_label(NULL))
  expect_null(kegg_species_label(""))
})

test_that("kegg_org_code gives the KEGG code of human, mouse and rat builds only", {
  expect_identical(kegg_org_code("Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30"), "hsa")
  expect_identical(kegg_org_code("Mus_musculus/GENCODE_M29/Genes"), "mmu")
  expect_identical(kegg_org_code("Rattus_norvegicus/Ensembl/Genes"), "rno")
  # Matched on the species, not on substrings such as "mm" in Mmul_10
  expect_null(kegg_org_code("Macaca_mulatta/Ensembl/Mmul_10"))
  expect_null(kegg_org_code(NA_character_))
  expect_null(kegg_org_code(NULL))
})

test_that("kegg_fetch_pathway returns the pathway name and its genes' Entrez ids", {
  skip_if_not_installed("KEGGREST")
  local_mocked_bindings(
    keggGet = function(id) list(list(NAME = "Cell cycle - Homo sapiens (human)")),
    keggLink = function(target, source) c("path:hsa04110" = "hsa:1017", "path:hsa04110" = "hsa:1019"),
    .package = "KEGGREST"
  )
  res <- kegg_fetch_pathway("hsa04110", "hsa")
  expect_identical(res$name, "Cell cycle")
  expect_identical(unname(res$entrez), c("1017", "1019"))
})

test_that("kegg_fetch_pathway passes the server's errors on", {
  skip_if_not_installed("KEGGREST")
  local_mocked_bindings(keggGet = function(id) stop("HTTP 404"), .package = "KEGGREST")
  expect_error(kegg_fetch_pathway("hsa99999", "hsa"), "HTTP 404")
})

test_that("kegg_map_genes matches the entrez_id column when the table has one", {
  de <- data.frame(feature_id = c("ENSG1", "ENSG2", "ENSG3"), entrez_id = c("1017", "999", "1019"))
  expect_identical(kegg_map_genes(de, c("1017", "1019", "5"), "hsa"), c("ENSG1", "ENSG3"))
})
