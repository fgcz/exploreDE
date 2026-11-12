# Heatmap helpers: the top annotation of the count heatmaps.

test_that("the top annotation gives each level its colour, drawing no random colours", {
  samples <- data.frame(
    Condition = factor(c("a", "a", "b")),
    Batch = factor(c("x", "y", "y"))
  )
  colours <- list(
    Condition = c(a = "#FF0000FF", b = "#0000FFFF"),
    Batch = c(x = "#00FF00FF", y = "#000000FF")
  )
  set.seed(1)
  seed <- .Random.seed
  ha <- heatmap_top_annotation(
    samples,
    c("Condition", "Batch"),
    colours,
    legend_param = list(),
    gp = grid::gpar(col = "white")
  )
  expect_identical(.Random.seed, seed)
  expect_s4_class(ha, "HeatmapAnnotation")
  for (f in names(colours)) {
    expect_identical(ha@anno_list[[f]]@color_mapping@colors, colours[[f]])
  }
})

test_that("heatmap_term_matrix stays a matrix for one row or none", {
  m <- matrix(
    c(1, 2, 3, -5, -6, -7, 0.1, 0.2, 0.1),
    3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("s1", "s2", "s3"))
  )
  expect_identical(heatmap_term_matrix(m, c("g3", "gX", "g1")), m[c("g3", "g1"), ])
  expect_identical(heatmap_term_matrix(m, "g2"), m["g2", , drop = FALSE])
  # |median| >= 5 leaves g2 only
  expect_identical(heatmap_term_matrix(m, c("g1", "g2", "g3"), min_abs_median = 5), m["g2", , drop = FALSE])
  expect_identical(dim(heatmap_term_matrix(m, c("g1", "g3"), min_abs_median = 10)), c(0L, 3L))
  expect_identical(dim(heatmap_term_matrix(m, character(0))), c(0L, 3L))
})
