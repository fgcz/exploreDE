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
