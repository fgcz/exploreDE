#' Start the exploreDE app
#'
#' Opens the app locally. With no `file`, it starts on the landing page, where
#' you can upload one or more `.h5ad` files built with [build_explore_h5ad()]
#' or load the demo dataset. With `file`, it opens straight into that dataset;
#' with several files, the first opens in the tabs and the contrasts of the
#' others join it in the Multi DE tabs.
#'
#' A `SummarizedExperiment` opens as [new_explore_dataset()] classifies it: a
#' SUSHI RNA-seq or prolfqua proteomics result with all its tabs, any other as
#' a Count QC dataset. To see the DE tabs for another `SummarizedExperiment`,
#' build an `.h5ad` from it with [build_explore_h5ad()] and open that.
#'
#' @param file Optional dataset: a `SummarizedExperiment` or AnnData object; a
#'   path to an exploreDE `.h5ad`, a SummarizedExperiment `.rds`/`.qs2`, or a
#'   SUSHI result directory; or several paths, the first of which opens in the
#'   tabs, the others going to Multi DE. An object is written to a temporary
#'   file, which the app opens.
#' @param ... Passed to [shiny::runApp()] (e.g. `port`, `launch.browser`).
#' @return Called to run the app; returns the value passed to [shiny::stopApp()]
#'   when the app stops (usually `NULL`), as [shiny::runApp()] does.
#' @export
#' @examples
#' if (interactive()) {
#'   exploreDE()
#'   exploreDE("my_data.h5ad")
#'   exploreDE(c("design1.h5ad", "design2.h5ad"))
#'   se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
#'   exploreDE(se)
#' }
exploreDE <- function(file = NULL, ...) {
  if (!is.null(file)) {
    file <- .launch_paths(file)
  }
  old <- options(exploreDE.file = file)
  on.exit(options(old), add = TRUE)
  shiny::runApp(appDir = system.file("app", package = "exploreDE"), ...)
}

# The paths exploreDE(file) hands to the app: a SummarizedExperiment is saved
# as a temporary .rds, an AnnData object written as a temporary .h5ad.
.launch_paths <- function(file) {
  if (inherits(file, "SummarizedExperiment")) {
    path <- tempfile(fileext = ".rds")
    saveRDS(file, path)
    file <- path
  } else if (inherits(file, "AbstractAnnData")) {
    path <- tempfile(fileext = ".h5ad")
    anndataR::write_h5ad(file, path, mode = "w")
    file <- path
  }
  normalizePath(file, mustWork = TRUE)
}
