test_that("col2hex returns #RRGGBB for colour names and hex strings, dropping alpha", {
  expect_identical(col2hex(c("red", "grey87", "#1B9E77", "#FF000080")), c("#FF0000", "#DEDEDE", "#1B9E77", "#FF0000"))
  expect_error(col2hex("notacolour"))
})
