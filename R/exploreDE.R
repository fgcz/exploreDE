#' start exploreDE application
#' @export
#' @examples
#' \dontrun{
#' myDir <- "/srv/GT/analysis/peter/rs_connect_apps2/exploreDE/SummarizedExperiment.rds"
#' exploreDE(myDir)
#' }
exploreDE <- function(fileSE = NA) {
  folder <- system.file("app", package = "exploreDE")
  message(folder)
  fileSE <- tools::file_path_as_absolute(fileSE)
  shiny::runApp(appDir = folder)
}
