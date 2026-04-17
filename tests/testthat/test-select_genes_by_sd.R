test_that("select_genes_by_sd returns top N rownames by descending SD", {
  set.seed(42)
  m <- matrix(rnorm(50), nrow = 10)
  # Force gene3 to have the largest SD and gene7 the second
  m[3, ] <- c(-100, 100, -100, 100, -100)
  m[7, ] <- c(-50, 50, -50, 50, -50)
  rownames(m) <- paste0("gene", 1:10)

  top2 <- select_genes_by_sd(m, top_n = 2)
  expect_equal(top2, c("gene3", "gene7"))
})

test_that("select_genes_by_sd returns all rownames when top_n >= nrow", {
  m <- matrix(1:20, nrow = 5, dimnames = list(letters[1:5], NULL))
  out <- select_genes_by_sd(m, top_n = 10)
  expect_setequal(out, letters[1:5])
})

test_that("select_genes_by_sd ignores NA values when computing SD", {
  m <- matrix(c(1, 2, 3, 4,  NA, 10, NA, 20,  5, 5, 5, 5),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("g1", "g2", "g3"), paste0("s", 1:4)))
  # SDs (with na.rm): g1 ≈ 1.29, g2 ≈ 7.07, g3 = 0
  top1 <- select_genes_by_sd(m, top_n = 1)
  expect_equal(top1, "g2")
})
