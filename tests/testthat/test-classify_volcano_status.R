# classify_volcano_status(): the significance class of each Volcano point.

test_that("classify_volcano_status classes points by significance and fold change", {
  y <- c(3, 3, 3, 3, 3, 3, 0.5, 0.5, 0.5, 0.5)
  r <- c(2, 1, 0.5, 0, -0.5, -1, 2, -2, 0.5, 1)
  expect_identical(
    classify_volcano_status(y, r, y_threshold = 1.3, lfc = 1),
    c(
      "FoldChangeSignificantUp", # beyond lfc, significant
      "FoldChangeSignificantUp", # at lfc
      "SignificantUp",
      "NotSignificant", # significant, no change
      "SignificantDown",
      "FoldChangeSignificantDown",
      "FoldChangeUp",
      "FoldChangeDown",
      "NotSignificant",
      "FoldChangeUp"
    )
  )
})

test_that("classify_volcano_status: with lfc 0 a zero ratio counts as down; NA is not significant", {
  expect_identical(
    classify_volcano_status(c(3, 0.5), c(0, 0), y_threshold = 1.3, lfc = 0),
    c("FoldChangeSignificantDown", "FoldChangeDown")
  )
  expect_identical(
    classify_volcano_status(c(NA, 3, NA), c(2, NA, NA), y_threshold = 1.3, lfc = 1),
    rep("NotSignificant", 3)
  )
})
