#' start exploreDE application
#' @export
#' @examples
#' \dontrun{
#' myDir <- system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE")
#' exploreDE(myDir)
#' }
exploreDE <- function(fileSE = NA, demo = NA) {
  folder <- system.file("app", package = "exploreDE")
  if (!is.na(demo) & is.na(fileSE)) {
    if (demo == "genomics") {
      fileSE <- system.file("examples/genomics", package = "exploreDE")
    } else if (demo == "proteomics") {
      fileSE <- system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE")
    }
  }
  if (is.na(demo) & is.na(fileSE)) {
    fileSE <- system.file("examples/genomics", package = "exploreDE")
  }
  fileSE <<- tools::file_path_as_absolute(fileSE)
  shiny::runApp(appDir = folder)
}
