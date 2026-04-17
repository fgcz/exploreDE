# ============================================================================
# build_qc_scatter_types — classifies features into colour-by categories
# ============================================================================
test_that("build_qc_scatter_types Standard mode labels features Present / Absent / Significants", {
  rd <- data.frame(
    feature_id = paste0("g", 1:5),
    log2Ratio  = c( 2.0, -3.0, 0.1,  1.5, NA),
    pValue     = c(0.01, 0.01, 0.5, 0.1, NA)
  )
  present_mask <- c(TRUE, TRUE, TRUE, FALSE, TRUE)

  res <- build_qc_scatter_types(rd, mode = "Standard", present_mask = present_mask,
                                p_thresh = 0.05, lfc_thresh = 1)
  # g1: present + sig (p<0.05 & |lfc|>=1) → Significants
  # g2: present + sig → Significants
  # g3: present + not-sig → Present
  # g4: absent → Absent
  # g5: NA p → Present (no override)
  expect_equal(res$types, c("Significants", "Significants", "Present", "Absent", "Present"))
  expect_setequal(names(res$colours), c("Absent", "Present", "Significants"))
  expect_equal(res$draw_order, c("Absent", "Present", "Significants"))
})

test_that("build_qc_scatter_types GC content mode bins by gc_lo / gc_hi", {
  rd <- data.frame(feature_id = paste0("g", 1:4),
                   gc = c(0.30, 0.50, 0.70, NA))
  res <- build_qc_scatter_types(rd, mode = "GC content", present_mask = rep(TRUE, 4),
                                gc_lo = 0.4, gc_hi = 0.6)
  # g1: low, g2: mid, g3: high, g4: mid (NA → falls into default label)
  expect_equal(res$types[1], "GC < 0.40")
  expect_equal(res$types[2], "0.40-0.60")
  expect_equal(res$types[3], "GC > 0.60")
  expect_equal(res$types[4], "0.40-0.60")
})

test_that("build_qc_scatter_types Gene width mode bins by w_lo / w_hi", {
  rd <- data.frame(feature_id = paste0("g", 1:4),
                   featWidth = c(300, 2000, 8000, NA))
  res <- build_qc_scatter_types(rd, mode = "Gene width", present_mask = rep(TRUE, 4),
                                w_lo = 500, w_hi = 5000)
  expect_equal(res$types[1], "< 500 nt")
  expect_equal(res$types[2], "500-5000 nt")
  expect_equal(res$types[3], "> 5000 nt")
})

test_that("build_qc_scatter_types default-falls-through to 'All' for unknown mode", {
  rd <- data.frame(feature_id = paste0("g", 1:3))
  res <- build_qc_scatter_types(rd, mode = "mystery", present_mask = rep(TRUE, 3))
  expect_equal(res$types, rep("All", 3))
  expect_equal(res$colours, c(All = "black"))
})

test_that("build_qc_scatter_types handles Standard mode without log2Ratio/pValue columns", {
  rd <- data.frame(feature_id = paste0("g", 1:3))  # no DE columns
  res <- build_qc_scatter_types(rd, mode = "Standard", present_mask = c(TRUE, TRUE, FALSE))
  expect_equal(res$types, c("Present", "Present", "Absent"))
})


# ============================================================================
# make_qc_scatter — wraps ggplot construction of the QC scatter
# ============================================================================
test_that("make_qc_scatter produces a ggplot object with the expected layers", {
  x <- c(10, 100, 1000, 10000)
  y <- c(11,  90, 1100,  9000)
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
  p <- make_qc_scatter(c(1, 2), c(1, 2), c("A", "A"), c(A = "black"),
                       xlab = "x", ylab = "y",
                       show_identity = FALSE, show_fold_bands = FALSE)
  # Only geom_point — no ablines
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false(any(grepl("GeomAbline", layer_classes)))
})

test_that("make_qc_scatter respects an explicit draw_order for z-ordering", {
  x <- 1:3
  y <- 1:3
  types <- c("A", "B", "C")
  colours <- c(A = "red", B = "blue", C = "green")

  p <- make_qc_scatter(x, y, types, colours, xlab = "x", ylab = "y",
                       draw_order = c("C", "A", "B"))
  expect_s3_class(p, "ggplot")
  # The drawing order factor reflects draw_order
  expect_equal(levels(p$data$types), c("C", "A", "B"))
})
