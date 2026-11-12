# ============================================================================
# qc_types_significance / qc_types_binned — the colour-by classes of the QC
# scatters (Standard; GC content and Gene width)
# ============================================================================
test_that("qc_types_significance labels features Present / Absent / Significants", {
  rd <- data.frame(
    feature_id = paste0("g", 1:5),
    log2Ratio = c(2.0, -3.0, 0.1, 1.5, NA),
    pValue = c(0.01, 0.01, 0.5, 0.1, NA)
  )
  res <- qc_types_significance(rd, c(TRUE, TRUE, TRUE, FALSE, TRUE), p_thresh = 0.05, lfc_thresh = 1)
  # g1, g2: p < 0.05 and |lfc| >= 1; g3: not significant; g4: absent; g5: NA p
  expect_equal(res$types, c("Significants", "Significants", "Present", "Absent", "Present"))
  expect_setequal(names(res$colours), c("Absent", "Present", "Significants"))
  expect_equal(res$draw_order, c("Absent", "Present", "Significants"))
})

test_that("qc_types_significance without DE columns marks presence only", {
  rd <- data.frame(feature_id = paste0("g", 1:3))
  expect_equal(qc_types_significance(rd, c(TRUE, TRUE, FALSE))$types, c("Present", "Present", "Absent"))
})

test_that("qc_types_binned bins GC content and gene width, NA in the middle bin", {
  gc <- qc_types_binned(c(0.30, 0.50, 0.70, NA), 0.4, 0.6, fmt = "%.2f", prefix = "GC ")
  expect_equal(gc$types, c("GC < 0.40", "0.40-0.60", "GC > 0.60", "0.40-0.60"))
  expect_equal(gc$draw_order, c("0.40-0.60", "GC < 0.40", "GC > 0.60"))
  width <- qc_types_binned(c(300, 2000, 8000), 500, 5000, fmt = "%.0f", suffix = " nt")
  expect_equal(width$types, c("< 500 nt", "500-5000 nt", "> 5000 nt"))
  expect_equal(names(width$colours), c("< 500 nt", "500-5000 nt", "> 5000 nt"))
})


# ============================================================================
# make_qc_scatter — wraps ggplot construction of the QC scatter
# ============================================================================
test_that("make_qc_scatter produces a ggplot object with the expected layers", {
  x <- c(10, 100, 1000, 10000)
  y <- c(11, 90, 1100, 9000)
  types <- c("Absent", "Present", "Present", "Significants")
  colours <- c(Absent = "grey", Present = "black", Significants = "red")

  p <- make_qc_scatter(x, y, types, colours, xlab = "Group A", ylab = "Group B")
  expect_s3_class(p, "ggplot")

  # Drawn layers: geom_point + geom_abline identity + 2x fold-band abline = 4
  expect_gte(length(p$layers), 4)
  expect_equal(p$labels$x, "Group A")
  expect_equal(p$labels$y, "Group B")
})

test_that("make_qc_scatter omits identity / fold-band lines when disabled", {
  p <- make_qc_scatter(
    c(1, 2),
    c(1, 2),
    c("A", "A"),
    c(A = "black"),
    xlab = "x",
    ylab = "y",
    show_identity = FALSE,
    show_fold_bands = FALSE
  )
  # Only geom_point — no ablines
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false(any(grepl("GeomAbline", layer_classes)))
})

test_that("make_qc_scatter respects an explicit draw_order for z-ordering", {
  x <- 1:3
  y <- 1:3
  types <- c("A", "B", "C")
  colours <- c(A = "red", B = "blue", C = "green")

  p <- make_qc_scatter(x, y, types, colours, xlab = "x", ylab = "y", draw_order = c("C", "A", "B"))
  expect_s3_class(p, "ggplot")
  # The drawing order factor reflects draw_order
  expect_equal(levels(p$data$types), c("C", "A", "B"))
})
