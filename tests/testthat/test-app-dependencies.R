# -----------------------------------------------------------------------------
# Every package the Shiny app uses is declared in DESCRIPTION. R CMD check does
# not look inside inst/app, so an undeclared package only shows up when a user
# on a fresh install reaches the code (pathview broke app startup that way).
# -----------------------------------------------------------------------------

test_that("packages used by inst/app are declared in DESCRIPTION", {
  app_dir <- system.file("app", package = "exploreDE")
  desc <- system.file("DESCRIPTION", package = "exploreDE")
  used <- character()
  for (f in list.files(app_dir, "\\.R$", full.names = TRUE)) {
    exprs <- parse(f, keep.source = TRUE)
    pd <- utils::getParseData(exprs)
    used <- c(used, pd$text[pd$token == "SYMBOL_PACKAGE"])
    # app.R attaches packages by name: packagesToLoad <- c("shiny", ...)
    for (e in exprs) {
      if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("packagesToLoad"))) {
        used <- c(used, eval(e[[3]], baseenv()))
      }
    }
  }
  fields <- read.dcf(desc, fields = c("Depends", "Imports", "Suggests"))
  declared <- trimws(sub("\\(.*", "", unlist(strsplit(paste(fields[!is.na(fields)], collapse = ","), ","))))
  base <- rownames(utils::installed.packages(priority = "base"))
  expect_setequal(setdiff(unique(used), c(declared, base, "exploreDE")), character(0))
})
