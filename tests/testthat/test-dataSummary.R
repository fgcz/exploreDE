test_that("dataSummary computes mean / sd / se per group", {
  df <- data.frame(
    group = rep(c("A", "B"), each = 4),
    value = c(1, 2, 3, 4,  10, 12, 14, 16)
  )
  out <- dataSummary(df, varname = "value", groupnames = "group")

  # Column layout: group, value (the mean, renamed), sd, se
  expect_true(all(c("group", "value", "sd", "se") %in% colnames(out)))
  expect_equal(out$value[out$group == "A"], 2.5)
  expect_equal(out$value[out$group == "B"], 13)
  expect_equal(out$sd[out$group == "A"],    sd(c(1, 2, 3, 4)))
  expect_equal(out$sd[out$group == "B"],    sd(c(10, 12, 14, 16)))
  expect_equal(out$se[out$group == "A"],    sd(c(1, 2, 3, 4))  / sqrt(4))
  expect_equal(out$se[out$group == "B"],    sd(c(10, 12, 14, 16)) / sqrt(4))
})

test_that("dataSummary handles NAs in the measured column", {
  df <- data.frame(
    group = c("A", "A", "A", "B", "B"),
    value = c(1, NA, 3, 5, 7)
  )
  out <- dataSummary(df, varname = "value", groupnames = "group")
  expect_equal(out$value[out$group == "A"], mean(c(1, 3)))
  expect_equal(out$value[out$group == "B"], 6)
})

test_that("dataSummary supports multi-column grouping", {
  df <- data.frame(
    group1 = rep(c("A", "B"), each = 4),
    group2 = rep(c("X", "Y"), times = 4),
    value  = 1:8
  )
  out <- dataSummary(df, varname = "value", groupnames = c("group1", "group2"))
  expect_equal(nrow(out), 4)
  expect_true(all(c("group1", "group2") %in% colnames(out)))
})
