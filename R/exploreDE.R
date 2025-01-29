#' start exploreDE application
#' @export
#' @examples
#' \dontrun{
#' myDir <- system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE")
#' exploreDE(myDir)
#' }
exploreDE <- function(fileSE = NA) {
  folder <- system.file("app", package = "exploreDE")
  message(folder)
  fileSE <- tools::file_path_as_absolute(fileSE)
  shiny::runApp(appDir = folder)
}
