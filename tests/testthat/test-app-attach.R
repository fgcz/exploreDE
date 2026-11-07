# -----------------------------------------------------------------------------
# The app attaches packages by name (packagesToLoad in app.R) and calls their
# functions unqualified. Each such call must reach a package in that list, not
# one a listed package happens to attach (tidyverse attached dplyr, tidyr and
# tibble that way), so dropping a package cannot silently break a call.
# -----------------------------------------------------------------------------

app_packages <- function() {
  for (e in parse(system.file("app", "app.R", package = "exploreDE"))) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("packagesToLoad"))) {
      return(eval(e[[3]], baseenv()))
    }
  }
}

# Unqualified function calls in inst/app, minus names the app or the package defines.
app_unqualified_calls <- function() {
  calls <- character()
  defined <- character()
  for (f in list.files(system.file("app", package = "exploreDE"), "\\.R$", full.names = TRUE)) {
    pd <- utils::getParseData(parse(f, keep.source = TRUE))
    pd <- pd[pd$terminal, ]
    pd <- pd[order(pd$line1, pd$col1), ]
    prev <- c("", utils::head(pd$token, -1))
    nxt <- c(pd$token[-1], "")
    is_call <- pd$token %in% c("SYMBOL_FUNCTION_CALL", "SPECIAL") & !prev %in% c("NS_GET", "NS_GET_INT", "'$'", "'@'")
    calls <- c(calls, pd$text[is_call])
    is_def <- (pd$token == "SYMBOL" & nxt %in% c("LEFT_ASSIGN", "EQ_ASSIGN")) | pd$token == "SYMBOL_FORMALS"
    defined <- c(defined, pd$text[is_def])
  }
  setdiff(unique(calls), c(defined, ls(asNamespace("exploreDE"), all.names = TRUE), telemetry_module))
}

# Defined by /srv/shiny-server/fgcz_apps/telemetry_module.R, which app.R sources when it exists
telemetry_module <- c("use_telemetry_js", "setup_telemetry_global", "start_telemetry_server")

exports_of <- function(pkg) suppressPackageStartupMessages(getNamespaceExports(pkg))

default_packages <- c("stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")

test_that("every unqualified call in inst/app reaches a package app.R attaches by name", {
  exported <- unique(unlist(lapply(c(app_packages(), default_packages), exports_of)))
  expect_setequal(setdiff(app_unqualified_calls(), exported), character(0))
})
