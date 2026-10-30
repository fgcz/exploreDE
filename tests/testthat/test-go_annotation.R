# go_annotation_columns() / features_in_go_term(): GO tabs must work with any
# subset of GO BP / MF / CC (non-reference organisms often have one ontology).

go_table <- function() {
  data.frame(
    feature_id = c("g1", "g2", "g3", "g4"),
    `GO BP` = c("GO:0006955; GO:0008150", "GO:0008150", NA, ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

test_that("go_annotation_columns finds any subset, in BP/MF/CC order", {
  expect_identical(go_annotation_columns(go_table()), "GO BP")
  d <- go_table()
  d$`GO CC` <- NA_character_
  d$`GO MF` <- "GO:0003674"
  expect_identical(go_annotation_columns(d), c("GO BP", "GO MF", "GO CC"))
  expect_identical(go_annotation_columns(data.frame(feature_id = "g1")), character(0))
})

test_that("features_in_go_term works with only GO BP present", {
  expect_identical(features_in_go_term(go_table(), "GO:0006955"), c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(features_in_go_term(go_table(), "GO:0008150"), c(TRUE, TRUE, FALSE, FALSE))
})

test_that("features_in_go_term looks across every GO column present", {
  d <- go_table()
  d$`GO MF` <- c(NA, NA, "GO:0003674", NA)
  expect_identical(features_in_go_term(d, "GO:0003674"), c(FALSE, FALSE, TRUE, FALSE))
})

test_that("features_in_go_term is all FALSE without GO columns, and matches literally", {
  expect_identical(features_in_go_term(data.frame(feature_id = c("a", "b")), "GO:0006955"), c(FALSE, FALSE))
  expect_false(any(features_in_go_term(go_table(), "GO:000695."))) # no regex
})
