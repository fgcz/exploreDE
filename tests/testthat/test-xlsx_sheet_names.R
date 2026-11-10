# xlsx_sheet_names() (R/downloadHelpers.R): the ORA / GSEA results files name a
# sheet per ontology and direction, and any stored ontology name can reach it.

test_that("xlsx_sheet_names keeps valid names", {
  expect_identical(xlsx_sheet_names(c("BP_upGenes", "MF_bothGenes")), c("BP_upGenes", "MF_bothGenes"))
})

test_that("xlsx_sheet_names replaces forbidden characters and truncates to 31", {
  expect_identical(xlsx_sheet_names("MSigDB:C2/[x]?*\\"), "MSigDB_C2__x____")
  long <- xlsx_sheet_names("WikiPathways_2024_Human_bothGenes")
  expect_identical(nchar(long), 31L)
  expect_identical(long, substr("WikiPathways_2024_Human_bothGenes", 1, 31))
  expect_identical(xlsx_sheet_names(""), "Sheet")
})

test_that("xlsx_sheet_names makes names unique, ignoring case", {
  x <- c("WikiPathways_2024_Human_upGenes_A", "WikiPathways_2024_Human_upGenes_B", "bp", "BP")
  out <- xlsx_sheet_names(x)
  expect_false(anyDuplicated(tolower(out)) > 0)
  expect_true(all(nchar(out) <= 31))
  expect_identical(out[3:4], c("bp", "BP_1"))
})
