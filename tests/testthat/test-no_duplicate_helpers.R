# Every helper has one definition, in R/. The app sources R/ first and then server-*.R into the
# server scope, so a function redefined there silently replaces the tested R/ version (issue #10).

app_function_names <- function(files) {
  is_fun_def <- function(e) {
    is.call(e) &&
      is.symbol(e[[1]]) &&
      as.character(e[[1]]) %in% c("<-", "=") &&
      is.symbol(e[[2]]) &&
      is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))
  }
  unlist(lapply(files, function(f) {
    defs <- Filter(is_fun_def, as.list(parse(f, keep.source = FALSE)))
    vapply(defs, function(e) as.character(e[[2]]), character(1))
  }))
}

test_that("app files do not redefine package functions", {
  app_dir <- system.file("app", package = "exploreDE")
  files <- list.files(app_dir, pattern = "^(server-.*|app)\\.R$", full.names = TRUE)
  expect_gt(length(files), 10)
  # Still defined in both places: they read server-scope reactives (next step of issue #10)
  pending <- c(
    "bucket_add",
    "bucket_remove",
    "bucket_clear",
    "resolve_features",
    "make_rep_heatmap",
    "make_rep_dendrogram"
  )
  redefined <- intersect(app_function_names(files), ls(asNamespace("exploreDE"), all.names = TRUE))
  expect_setequal(setdiff(redefined, pending), character(0))
})
