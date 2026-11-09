#' Start the exploreDE app
#'
#' Opens the app locally. With no `file`, it starts on the landing page, where
#' you can upload one or more `.h5ad` files built with [build_explore_h5ad()]
#' or load the demo dataset. With `file`, it opens straight into that dataset;
#' with several files, the first opens in the tabs and the contrasts of the
#' others join it in the Multi DE tabs.
#'
#' @param file Optional path to an exploreDE `.h5ad`, a SummarizedExperiment
#'   `.rds`/`.qs2`, or a SUSHI result directory; or several paths, the first of
#'   which opens in the tabs, the others going to Multi DE.
#' @param ... Passed to [shiny::runApp()] (e.g. `port`, `launch.browser`).
#' @return Called to run the app; returns the value passed to [shiny::stopApp()]
#'   when the app stops (usually `NULL`), as [shiny::runApp()] does.
#' @export
#' @examples
#' if (interactive()) {
#'   exploreDE()
#'   exploreDE("my_data.h5ad")
#'   exploreDE(c("design1.h5ad", "design2.h5ad"))
#' }
exploreDE <- function(file = NULL, ...) {
  if (!is.null(file)) {
    file <- normalizePath(file, mustWork = TRUE)
  }
  old <- options(exploreDE.file = file)
  on.exit(options(old), add = TRUE)
  shiny::runApp(appDir = system.file("app", package = "exploreDE"), ...)
}
