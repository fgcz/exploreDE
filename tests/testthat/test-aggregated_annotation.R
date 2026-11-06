design <- data.frame(
  Name = paste0("s", 1:6),
  Condition = rep(c("WT", "KO"), each = 3),
  Batch = rep(c("b1", "b2"), each = 3),
  Day = c("d1", "d2", "d1", "d1", "d1", "d1")
)

test_that("aggregated_annotation gives one value per level for a factor constant within levels", {
  ann <- aggregated_annotation(design, "Condition", "Batch", c("WT", "KO"))
  expect_false(ann$varies)
  expect_identical(ann$values, c(WT = "b1", KO = "b2"))
})

test_that("aggregated_annotation flags a factor that varies within a level", {
  ann <- aggregated_annotation(design, "Condition", "Day", c("WT", "KO"))
  expect_true(ann$varies)
})
