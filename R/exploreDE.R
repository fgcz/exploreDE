#' start exploreDE application
#' @export
#' @examples
#' \dontrun{
#' myDir <- "/srv/GT/analysis/peter/rs_connect_apps2/exploreDE/SummarizedExperiment.rds"
#' exploreDE(myDir)
#' }
exploreDE <- function(dir = NA) {
  folder <- system.file("app", package = "exploreDE")
  shiny::runApp(appDir = folder)
  }