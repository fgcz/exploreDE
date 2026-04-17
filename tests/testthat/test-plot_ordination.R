# ============================================================================
# plot_ordination_samples / plot_ordination_biplot
#
# These take `input` as a parameter, so we just build a mock list with every
# reactive-input field the function reads. Smoke-level coverage: the function
# returns a ggplot object without erroring.
# ============================================================================

make_mock_pca_input <- function(overrides = list()) {
  base <- list(
    pcaX = "PC1", pcaY = "PC2",
    pcaFactor1 = "cond", pcaFactor2 = "None",
    pcaAddEllipses = FALSE, pcaEllipsesAlpha = 0.2,
    dotSizePCA = 3, alphaPCA = 0.8, dotBorderPCA = 0.3,
    pcaShowNames = FALSE,
    geneLabelNudgePCAY = 0, geneLabelNudgePCAX = 0,
    textSizePCA = 12, pcaLabelMaxOverlap = 10,
    pcaAxesProp = FALSE,
    boldPCA = "FALSE",
    showLinesPCA = TRUE, showAxesPCA = TRUE
  )
  modifyList(base, overrides)
}


test_that("plot_ordination_samples returns a ggplot for a PCA result", {
  m <- matrix(rnorm(200), nrow = 10,
              dimnames = list(paste0("s", 1:10), paste0("g", 1:20)))
  pca <- calculate_pca(m)

  scores <- dplyr::left_join(
    pca$scores,
    tibble::tibble(Name = paste0("s", 1:10), cond = rep(c("A", "B"), each = 5)),
    by = "Name"
  )
  coloursPCA <- c(A = "red", B = "blue")
  input <- make_mock_pca_input()

  p <- plot_ordination_samples(pca, scores, input, coloursPCA)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$x, "PC1")
  expect_match(p$labels$y, "PC2")
})

test_that("plot_ordination_samples honours pcaShowNames and pcaAddEllipses toggles", {
  m <- matrix(rnorm(200), nrow = 10,
              dimnames = list(paste0("s", 1:10), paste0("g", 1:20)))
  pca <- calculate_pca(m)
  scores <- dplyr::left_join(
    pca$scores,
    tibble::tibble(Name = paste0("s", 1:10), cond = rep(c("A", "B"), each = 5)),
    by = "Name"
  )
  coloursPCA <- c(A = "red", B = "blue")

  input <- make_mock_pca_input(list(pcaShowNames = TRUE, pcaAddEllipses = TRUE))
  p <- plot_ordination_samples(pca, scores, input, coloursPCA)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ordination_samples renders a non-PCA result (NMDS-style)", {
  m <- matrix(rnorm(200), nrow = 10,
              dimnames = list(paste0("s", 1:10), paste0("g", 1:20)))
  nmds <- calculate_pcoa(m, nmds = TRUE)
  scores <- dplyr::left_join(
    nmds$scores,
    tibble::tibble(Name = paste0("s", 1:10), cond = rep(c("A", "B"), each = 5)),
    by = "Name"
  )
  coloursPCA <- c(A = "red", B = "blue")
  input <- make_mock_pca_input(list(pcaX = "NMDS1", pcaY = "NMDS2"))

  p <- plot_ordination_samples(nmds, scores, input, coloursPCA)
  expect_s3_class(p, "ggplot")
  # NMDS axes show stress rather than variance
  expect_match(p$labels$x, "stress")
})


test_that("plot_ordination_biplot renders a biplot with gene arrows", {
  m <- matrix(rnorm(200), nrow = 10,
              dimnames = list(paste0("s", 1:10), paste0("g", 1:20)))
  pca <- calculate_pca(m)
  scores <- dplyr::left_join(
    pca$scores,
    tibble::tibble(Name = paste0("s", 1:10), cond = rep(c("A", "B"), each = 5)),
    by = "Name"
  )
  coloursPCA <- c(A = "red", B = "blue")
  bp <- prepare_biplot_data(pca, top_n = 5)

  input <- make_mock_pca_input()
  p <- plot_ordination_biplot(pca, scores, input, coloursPCA, bp,
                              arrow_scale = 1, arrow_colour = "red")
  expect_s3_class(p, "ggplot")
  # Gene-arrow layers: point + segment + text_repel = 3
  expect_gte(length(p$layers), 3)
})
