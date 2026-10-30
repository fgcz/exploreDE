# The column dictionary when exploreDE is source()'d from a checkout, as app.R
# does at start-up: The checkout's own copy must be found, from any working
# directory, even though the patterns load lazily after source() returns.

test_that("a source()d utils_explore_dataset.R records its own path", {
  f <- normalizePath(test_path("..", "..", "R", "utils_explore_dataset.R"))
  skip_if_not(file.exists(f))
  env <- new.env()
  withr::with_dir(tempdir(), source(f, local = env))
  expect_identical(env$.eds_source_file, f)
})

test_that("the lazily loaded dictionary comes from the checkout, not the working directory", {
  f <- normalizePath(test_path("..", "..", "R", "utils_explore_dataset.R"))
  skip_if_not(file.exists(f))
  env <- new.env()
  withr::with_dir(tempdir(), source(f, local = env))
  pats <- withr::with_dir(tempdir(), env$eds_load_column_patterns())
  json <- jsonlite::fromJSON(file.path(dirname(dirname(f)), "inst", "dictionaries", "column_patterns.json"),
                             simplifyVector = FALSE)
  expect_identical(names(pats$de), names(json$de))
  expect_identical(pats$de$score_fdr, as.character(unlist(json$de$score_fdr)))
})
