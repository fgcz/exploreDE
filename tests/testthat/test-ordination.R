# Small deterministic count matrix reused across several ordination tests.
# Samples as ROWS (that's the convention calculate_pca / calculate_pcoa expect).
make_sample_counts <- function(n_samples = 20, n_features = 50, seed = 1) {
  set.seed(seed)
  m <- matrix(
    rpois(n_samples * n_features, lambda = 50),
    nrow = n_samples,
    ncol = n_features,
    dimnames = list(paste0("s", seq_len(n_samples)), paste0("g", seq_len(n_features)))
  )
  m
}


# ============================================================================
# calculate_pca
# ============================================================================
test_that("calculate_pca returns a well-shaped ordination result", {
  m <- make_sample_counts()
  res <- calculate_pca(m)
  expect_type(res, "list")
  expect_equal(res$method, "PCA")
  expect_true("Name" %in% colnames(res$scores))
  expect_true("gene" %in% colnames(res$loadings))
  expect_equal(length(res$variance), nrow(m))
  expect_equal(res$component_names, paste0("PC", 1:nrow(m)))
  expect_equal(round(sum(res$variance), 2), 100)
})

test_that("calculate_pca drops samples with NAs and reports removed_samples", {
  m <- make_sample_counts()
  m[3, 1] <- NA
  res <- calculate_pca(m)
  expect_equal(res$removed_samples, 1)
  expect_false("s3" %in% res$scores$Name)
})

test_that("calculate_pca returns NULL when matrix is too small", {
  expect_null(calculate_pca(matrix(1:4, nrow = 2)))
})


# ============================================================================
# calculate_pcoa (classical PCoA + NMDS paths)
# ============================================================================
test_that("calculate_pcoa returns PCoA result with scores and variance", {
  m <- make_sample_counts()
  res <- calculate_pcoa(m, method = "euclidean", nmds = FALSE)
  expect_equal(res$method, "PCoA")
  expect_true(all(grepl("^PCo", res$component_names)))
  expect_true("Name" %in% colnames(res$scores))
  expect_null(res$loadings) # no loadings for distance-based methods
})

test_that("calculate_pcoa NMDS mode returns stress instead of variance", {
  m <- make_sample_counts(n_samples = 15)
  res <- calculate_pcoa(m, method = "euclidean", nmds = TRUE)
  expect_equal(res$method, "NMDS")
  expect_null(res$variance)
  expect_true(is.numeric(res$stress))
})


# ============================================================================
# calculate_glmpca
# ============================================================================
test_that("calculate_glmpca runs on a small raw count matrix", {
  # glmpca expects genes x samples, raw counts
  set.seed(2)
  m <- matrix(rpois(20 * 10, lambda = 20), nrow = 20, ncol = 10, dimnames = list(paste0("g", 1:20), paste0("s", 1:10)))
  res <- calculate_glmpca(m, L = 2, family = "nb")
  expect_false(is.null(res))
  expect_equal(res$method, "GLM-PCA")
  expect_equal(res$component_names, c("dim1", "dim2"))
  expect_true("Name" %in% colnames(res$scores))
})

test_that("calculate_glmpca returns NULL for invalid family parameter", {
  m <- matrix(rpois(100, lambda = 10), nrow = 10)
  expect_warning(res <- calculate_glmpca(m, L = 2, family = "bogus"))
  expect_null(res)
})

test_that("calculate_glmpca returns NULL when too few non-zero genes remain", {
  m <- matrix(0, nrow = 10, ncol = 5)
  m[1, ] <- 1
  expect_warning(res <- calculate_glmpca(m, L = 2, family = "nb"))
  expect_null(res)
})


# ============================================================================
# prepare_biplot_data
# ============================================================================
test_that("prepare_biplot_data returns NULL for distance-based results without loadings", {
  m <- make_sample_counts()
  res <- calculate_pcoa(m, nmds = FALSE) # loadings = NULL
  expect_null(prepare_biplot_data(res))
})

test_that("prepare_biplot_data 'overall' picks top N by loading magnitude", {
  m <- make_sample_counts()
  pca <- calculate_pca(m)
  bp <- prepare_biplot_data(pca, top_n = 5, selection_method = "overall")
  expect_equal(nrow(bp$gene_data), 5)
  expect_equal(bp$pc_x_name, "PC1")
  expect_equal(bp$pc_y_name, "PC2")
})

test_that("prepare_biplot_data 'per_pc' selects top N from each axis and deduplicates", {
  m <- make_sample_counts()
  pca <- calculate_pca(m)
  bp <- prepare_biplot_data(pca, top_n = 3, selection_method = "per_pc")
  expect_lte(nrow(bp$gene_data), 6) # at most 2 × 3
  expect_true("pc_source" %in% colnames(bp$gene_data))
})

test_that("prepare_biplot_data 'manual' selects only the named genes", {
  m <- make_sample_counts()
  pca <- calculate_pca(m)
  picks <- c("g1", "g5", "g42")
  bp <- prepare_biplot_data(pca, custom_genes = picks, selection_method = "manual")
  expect_setequal(bp$gene_data$gene, picks)
})


# ============================================================================
# step_* cascade helpers
# ============================================================================
test_that("step_group_subset subsets counts and dataset to selected groups", {
  counts <- matrix(1:20, nrow = 5, dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  # Multi-column data.frame — single-column df would drop to a vector on row subset
  ds <- data.frame(
    Name = paste0("s", 1:4),
    cond = c("A", "A", "B", "B"),
    row.names = paste0("s", 1:4)
  )

  out <- step_group_subset(counts, ds, "cond", "A")
  expect_equal(ncol(out$counts), 2)
  expect_equal(rownames(out$dataset), c("s1", "s2"))
})

test_that("step_peptide_filter keeps only features with nrPeptides >= threshold", {
  counts <- matrix(1:12, nrow = 3, dimnames = list(c("f1", "f2", "f3"), paste0("s", 1:4)))
  sa <- list(data.frame(feature_id = c("f1", "f2", "f3"), nrPeptides = c(3, 1, 5)))
  ds <- data.frame(row.names = paste0("s", 1:4), cond = c("A", "A", "B", "B"))

  out <- step_peptide_filter(counts, ds, sa, min_peptides = 3)
  expect_setequal(rownames(out$counts), c("f1", "f3"))
})

test_that("step_peptide_filter is a no-op when seqAnnoList or min_peptides is NULL", {
  counts <- matrix(1:12, nrow = 3, dimnames = list(c("f1", "f2", "f3"), paste0("s", 1:4)))
  ds <- data.frame(row.names = paste0("s", 1:4), cond = rep("A", 4))
  expect_equal(step_peptide_filter(counts, ds, NULL, 3)$counts, counts)
  expect_equal(step_peptide_filter(counts, ds, list(), NULL)$counts, counts)
})

test_that("step_log2 log-transforms raw / normalised / TPM / FPKM counts only", {
  counts <- matrix(1:4, nrow = 2)
  ds <- data.frame(cond = c("A", "B"))
  expect_equal(step_log2(counts, ds, "Raw")$counts, log2(counts + 1))
  expect_equal(step_log2(counts, ds, "Normalised")$counts, log2(counts + 1))
  expect_equal(step_log2(counts, ds, "TPM")$counts, log2(counts + 1))
  # Already-log counts: no-op
  expect_equal(step_log2(counts, ds, "Normalised + Log2")$counts, counts)
  expect_equal(step_log2(counts, ds, "VST")$counts, counts)
})

test_that("step_top_n_sd keeps the top N most variable features", {
  set.seed(1)
  m <- matrix(rnorm(50), nrow = 10)
  rownames(m) <- paste0("g", 1:10)
  # Force g4 to have largest SD
  m[4, ] <- c(-100, 100, -100, 100, -100)
  ds <- data.frame(cond = rep("A", 5))

  out <- step_top_n_sd(m, ds, top_n = 3)
  expect_equal(nrow(out$counts), 3)
  expect_true("g4" %in% rownames(out$counts))
})

test_that("step_top_n_sd is a no-op when top_n >= nrow", {
  m <- matrix(1:20, nrow = 5, dimnames = list(letters[1:5], NULL))
  ds <- data.frame(cond = rep("A", 4))
  expect_equal(step_top_n_sd(m, ds, top_n = 100)$counts, m)
})

test_that("step_batch_correction applies when batch_column points to a valid factor", {
  set.seed(3)
  counts <- matrix(
    c(rnorm(5, 10), rnorm(5, 11), rnorm(5, 20), rnorm(5, 21)),
    nrow = 5,
    dimnames = list(paste0("g", 1:5), paste0("s", 1:4))
  )
  ds <- data.frame(row.names = paste0("s", 1:4), batch = c("x", "x", "y", "y"), cond = c("A", "B", "A", "B"))
  out <- step_batch_correction(counts, ds, "batch", "cond")
  expect_equal(dim(out$counts), dim(counts))
  expect_equal(dimnames(out$counts), dimnames(counts))
})

test_that("step_batch_correction is a no-op when batch_column is 'None' or NULL", {
  counts <- matrix(runif(20), nrow = 5, dimnames = list(NULL, paste0("s", 1:4)))
  ds <- data.frame(row.names = paste0("s", 1:4), cond = c("A", "A", "B", "B"))
  expect_identical(step_batch_correction(counts, ds, "None", "cond")$counts, counts)
  expect_identical(step_batch_correction(counts, ds, NULL, "cond")$counts, counts)
})
