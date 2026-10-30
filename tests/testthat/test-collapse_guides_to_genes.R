test_that("collapse_guides_to_genes sums sgRNA counts to gene-level", {
  counts <- matrix(
    c(
      10,
      20,
      30,
      40, # sgRNA1 → GENE_A
      50,
      60,
      70,
      80, # sgRNA2 → GENE_A
      1,
      2,
      3,
      4
    ), # sgRNA3 → GENE_B
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("sgRNA1", "sgRNA2", "sgRNA3"), paste0("s", 1:4))
  )
  map <- c(sgRNA1 = "GENE_A", sgRNA2 = "GENE_A", sgRNA3 = "GENE_B")

  out <- collapse_guides_to_genes(counts, map)
  expect_equal(nrow(out), 2)
  expect_setequal(rownames(out), c("GENE_A", "GENE_B"))
  expect_equal(unname(out["GENE_A", ]), c(60, 80, 100, 120))
  expect_equal(unname(out["GENE_B", ]), c(1, 2, 3, 4))
})

test_that("collapse_guides_to_genes respects the exclude argument", {
  counts <- matrix(1:12, nrow = 3, dimnames = list(c("sgRNA1", "sgRNA2", "sgRNA3"), paste0("s", 1:4)))
  map <- c(sgRNA1 = "GENE_A", sgRNA2 = "GENE_B", sgRNA3 = "negCtrl_1")

  out <- collapse_guides_to_genes(counts, map, exclude = "negCtrl_1")
  expect_false("negCtrl_1" %in% rownames(out))
  expect_setequal(rownames(out), c("GENE_A", "GENE_B"))
})

test_that("collapse_guides_to_genes intersects rownames with the sgRNA_to_gene map", {
  counts <- matrix(1:12, nrow = 3, dimnames = list(c("sgRNA1", "sgRNA2", "orphan"), paste0("s", 1:4)))
  map <- c(sgRNA1 = "GENE_A", sgRNA2 = "GENE_B") # "orphan" not mapped

  out <- collapse_guides_to_genes(counts, map)
  expect_setequal(rownames(out), c("GENE_A", "GENE_B"))
  expect_equal(nrow(out), 2)
})
