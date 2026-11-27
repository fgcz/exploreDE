# Runs the UI tests (ui-tests.yaml) and reports every failure and error as a
# GitHub annotation, so the failing test, its line and its message can be read
# from the check run without opening the job log.
results <- testthat::test_dir(
  "tests/testthat",
  package = "exploreDE",
  load_package = "installed",
  filter = Sys.getenv("UI_TEST_FILTER", "^ui-"),
  stop_on_failure = FALSE
)

# A workflow command's message must have its %, CR and LF escaped.
escape <- function(x) {
  x <- gsub("%", "%25", x, fixed = TRUE)
  x <- gsub("\r", "%0D", x, fixed = TRUE)
  gsub("\n", "%0A", x, fixed = TRUE)
}

n_bad <- 0L
for (test in results) {
  for (e in test$results) {
    if (!inherits(e, c("expectation_failure", "expectation_error"))) {
      next
    }
    n_bad <- n_bad + 1L
    line <- if (!is.null(e$srcref)) e$srcref[1] else 1L
    msg <- substr(paste0(test$test, ": ", conditionMessage(e)), 1L, 4000L)
    cat(sprintf("::error file=tests/testthat/%s,line=%d::%s\n", test$file, line, escape(msg)))
  }
}
if (n_bad > 0L) {
  cat(n_bad, "UI test failure(s)\n")
  quit(status = 1L)
}
