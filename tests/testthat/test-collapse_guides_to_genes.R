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

test_that("collapse_assay_to_genes sums linear assays per gene", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  map <- c(g1 = "A", g2 = "A", g3 = "B")
  out <- collapse_assay_to_genes(m, map, is_log2 = FALSE)
  expect_equal(unname(out["A", ]), c(3, 9))
  expect_equal(unname(out["B", ]), c(3, 6))
})

test_that("collapse_assay_to_genes sums log2 assays on the linear scale", {
  lin <- matrix(c(1, 3, 8, 4, 4, 16), nrow = 3, dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  map <- c(g1 = "A", g2 = "A", g3 = "B")
  out <- collapse_assay_to_genes(log2(lin), map, is_log2 = TRUE)
  expect_equal(unname(out["A", ]), log2(c(4, 8)))
  expect_equal(unname(out["B", ]), log2(c(8, 16)))
})

test_that("collapse_guides_to_genes counts present guides for a logical assay", {
  pf <- matrix(c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE), nrow = 3, dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  out <- collapse_guides_to_genes(pf, c(g1 = "A", g2 = "A", g3 = "B"))
  expect_equal(unname(out["A", ]), c(1L, 1L))
  expect_equal(unname(out["B", ]), c(1L, 0L))
})

test_that("CrisprCountQC takes the gene from rowData when the guide name defeats the regex", {
  rn <- c(
    "Chd9_+_90828456.23-ENSMUST00000109614.2_ENSMUST00000048665.6_GAACACCCTGGGAGCATCCG",
    "Chd9_-_90828401.23-P1P2_GAGGGTGGGAGGTGACGTGA",
    "Actb_1_ACGTACGTACGTACGTACGT"
  )
  counts <- matrix(1:6, nrow = 3, dimnames = list(rn, c("s1", "s2")))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(Gene = c("Chd9", "Chd9", NA))
  )
  map <- ExploreDatasetSE_CrisprCountQC$new(se)$sgrna_to_gene()
  expect_equal(unname(map[rn]), c("Chd9", "Chd9", "Actb")) # NA in rowData falls back to the regex
})
