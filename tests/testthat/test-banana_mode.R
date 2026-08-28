# Banana mode pure helpers (R/banana_mode.R). The marker PNG is embedded, so
# no fixtures / options are needed. Rendering (banana_grid, banana_layer output)
# is device work and covered by manual render checks, not here.

test_that("banana_png_path decodes the embedded PNG to a real file", {
  p <- banana_png_path()
  expect_true(file.exists(p))
  expect_gt(file.size(p), 1000)
  expect_identical(readBin(p, "raw", 8),
                   as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a)))  # PNG magic
  expect_identical(banana_png_path(), p)  # cached
})

test_that("banana_tier maps p-values to 0-4 on the shared cutpoints", {
  expect_identical(banana_tier(1e-6),  4L)
  expect_identical(banana_tier(1e-4),  4L)
  expect_identical(banana_tier(5e-4),  3L)
  expect_identical(banana_tier(5e-3),  2L)
  expect_identical(banana_tier(0.03),  1L)
  expect_identical(banana_tier(0.2),   0L)
  expect_identical(banana_tier(NA_real_), 0L)
  expect_identical(banana_tier(numeric(0)), 0L)
})

test_that("banana_on is just the toggle state (asset always embedded)", {
  expect_true(banana_on(TRUE))
  expect_false(banana_on(FALSE))
  expect_false(banana_on(NA))
  expect_false(banana_on(NULL))
})

test_that("banana_md builds n <img> tags, empty when n < 1", {
  expect_equal(banana_md(0), "")
  expect_equal(banana_md(-1), "")
  expect_equal(banana_md(NA_integer_), "")
  three <- banana_md(3, px = 20)
  expect_equal(lengths(regmatches(three, gregexpr("<img ", three))), 3L)
  expect_match(three, "height='20'")
  expect_match(three, "\\.png'")
})

test_that("banana_layer returns NULL when nothing is significant or cols missing", {
  d_ns <- data.frame(Feature = c("A", "B"), fdr = c(0.5, 0.9))
  expect_null(banana_layer(d_ns, x = 1.5, y = 10, fdr_col = "fdr"))
  expect_null(banana_layer(data.frame(x = 1), x = 1, y = 1, fdr_col = "fdr"))

  d_sig <- data.frame(Feature = c("A", "B"), fdr = c(1e-5, 0.9))
  lyr <- banana_layer(d_sig, x = 1.5, y = c(10, 12), fdr_col = "fdr")
  expect_s3_class(lyr, "LayerInstance")   # one row (A) survives
  expect_null(banana_layer(d_sig, x = 1.5, y = c(NA_real_, NA_real_), fdr_col = "fdr"))
})
