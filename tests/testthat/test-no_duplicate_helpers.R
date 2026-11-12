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
  redefined <- intersect(app_function_names(files), ls(asNamespace("exploreDE"), all.names = TRUE))
  expect_setequal(redefined, character(0))
})

# app.R sources every server module into ONE environment, so a top-level name
# assigned in two modules is one variable: the module sourced last replaces the
# other's value for both (designReactive had five definitions with two
# different bodies). `eds` is the exception: every module sets it to the same
# inputDataReactive()$.eds.
test_that("no two server modules assign the same top-level name", {
  app_dir <- system.file("app", package = "exploreDE")
  app_src <- readLines(file.path(app_dir, "app.R"))
  modules <- unique(regmatches(app_src, regexpr("server-[A-Za-z0-9-]+\\.R", app_src)))
  expect_gt(length(modules), 10)
  is_assign <- function(e) {
    is.call(e) && is.symbol(e[[1]]) && as.character(e[[1]]) %in% c("<-", "=") && is.symbol(e[[2]])
  }
  owners <- list()
  for (m in modules) {
    exprs <- Filter(is_assign, as.list(parse(file.path(app_dir, m), keep.source = FALSE)))
    for (nm in unique(vapply(exprs, function(e) as.character(e[[2]]), character(1)))) {
      owners[[nm]] <- c(owners[[nm]], m)
    }
  }
  owners$eds <- NULL
  shared <- Filter(function(x) length(x) > 1, owners)
  expect_identical(
    vapply(shared, paste, character(1), collapse = ", "),
    setNames(character(0), character(0))
  )
})
