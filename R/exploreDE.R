#' Start the exploreDE app
#'
#' Opens the app locally. With no `file`, it starts on the landing page, where
#' you can upload an `.h5ad` built with [build_explore_h5ad()] or load the demo
#' dataset. With `file`, it opens straight into that dataset.
#'
#' @param file Optional path to an exploreDE `.h5ad`, a SummarizedExperiment
#'   `.rds`/`.qs2`, or a SUSHI result directory.
#' @param ... Passed to [shiny::runApp()] (e.g. `port`, `launch.browser`).
#' @export
#' @examples
#' \dontrun{
#' exploreDE()
#' exploreDE("my_data.h5ad")
#' }
exploreDE <- function(file = NULL, ...) {
  if (!is.null(file)) file <- normalizePath(file, mustWork = TRUE)
  old <- options(exploreDE.file = file)
  on.exit(options(old), add = TRUE)
  shiny::runApp(appDir = system.file("app", package = "exploreDE"), ...)
}
