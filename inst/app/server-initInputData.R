# Start from here:
queryList = parseQueryString(session$clientData$url_search) 
if (is.list(queryList)){
  dataUrl <- queryList$data
} else {
  dataUrl <- NULL
}

if (!is.null(dataUrl)) {
  # Added this so we can read from both genomics and proteomics servers
  if (!grepl("Proteomics|prolfqua", dataUrl, ignore.case = TRUE)) {
    urlDataRoot = c("/srv/gstore/projects", "/srv/GT/analysis/course_sushi/public/gstore/projects")
    dataDir <- file.path(urlDataRoot, dataUrl)
    dataDir <- dataDir[file.exists(dataDir)][1]
    if (!file.exists(dataDir)){
      showModal(modalDialog(
        title = "Something went wrong",
        "It looks like either the dataset you're looking for doesn't exist, or has not finished being processed in SUSHI yet."
      ))
      stopApp(returnValue = invisible())
    }
  } else {
    dataDir <- paste0("https://fgcz-ms.uzh.ch/public/pStore/", dataUrl)
  }
} else if (is.null(dataUrl) & !exists("fileSE")) {
  dataDir <- "/srv/gstore/projects/p3009/o5638_DESeq2_diff--over--undiff_2024-11-20--12-29-40/diff--over--undiff/"
  showNotification("Since you did not specify a dataset in the URL, you are seeing a demo dataset.", type = "message", duration = NULL, closeButton = TRUE)
}

# 2025-01-29: Read in local proteomics file if specified 
if(exists("fileSE")) {
  dataDir <- fileSE
}

if(!exists("dataDir")) {
  showModal(modalDialog(
    title = "Something went wrong",
    "It looks like either the dataset you're looking for doesn't exist, or has not finished being processed in SUSHI yet."
  ))
  stop()
}

# Import proteomics data from pStore ----
is_url <- function(dataDir) {
  return(grepl("^https?://", dataDir))  # Checks if it starts with http:// or https://
}
# Import proteomics data from pStore
if (grepl("rds", dataDir) & grepl("Proteomics|prolfqua", dataDir) || grepl("rds", dataDir) & exists("fileSE")) {
  if (is_url(dataDir)) {
    se <- readRDS(url(dataDir))  # If it's a URL
  } else {
    se <- readRDS(dataDir)       # If it's a local file path
  }
} else if (grepl("zip", dataDir)) {
  myTempFile <- tempfile(fileext=".rds", tmpdir = ".") #
  res = system(
    paste("ssh fgcz-ms.uzh.ch unzip -c -qq /srv/www/htdocs/p30385/bfabric/Proteomics/DEA_FragPipe-DiaNN/2023/2023-09/2023-09-22/workunit_293901/2360252.zip '\\*rds' > ", myTempFile))
  if (res != 0){
    message("download failed")
  }
  se <- readRDS(myTempFile)
}
if (grepl("rds|zip", dataDir) & grepl("Proteomics|prolfqua", dataDir) || grepl("rds|zip", dataDir) & exists("fileSE")) {
  contrasts <- names(rowData(se))[grep("^constrast_", names(rowData(se)))] %>% gsub("constrast_", "", .)
  output$proteomicsContrastSelectorUI <- renderUI({
    selectInput(inputId = "contrastSelected", label = "Select contrast to view", choices = contrasts, selected = contrasts[1], multiple = F, selectize = T)
  })
}

# Import RNA seq data from SUSHI ----
if (grepl("gstore", dataDir) | exists("fileSE")) {
  if (grepl("EzResult.RData", dataDir)) {
    dataDir <- gsub("\\/result-.*.-EzResult.RData", "", dataDir)
  }
  if (file.exists(file.path(dataDir, "deResult.rds"))) {
    se <- readRDS(file.path(dataDir, "deResult.rds"))
  } else {
    showModal(modalDialog(
      title = "The file does not exist", 
      "Either the analysis has not yet finished running, you have made a mistake in the URL, or you have not pointed to any dataset. Please try again! If the issue persists, email peter.leary@uzh.ch",
      easyClose = TRUE,
      footer = NULL
    ))
  }
}

# Generate inputDataReactive ----
inputDataReactive <- reactive({
  waiter <- waiter::Waiter$new(fadeout = TRUE, color = "#86a5bf")
  waiter$show()
  on.exit(waiter$hide())
  
  # Load proteomics data ----
  if (grepl("rds|zip", dataDir)) {
    return(convert_proteomics_se(se))
  }
  
  # Load RNA Seq data ----
  if (file.exists(file.path(dataDir, "deResult.rds"))) {
    return(convert_genomics_se(se, dataDir))
  }
})
inputDataReactive()$dataType

figuresDataReactive <- reactiveValues(
  "volcanoStatic" = NULL,
  "heatmapBoth Directions" = NULL, "heatmapUp-Regulated" = NULL, "heatmapDown-Regulated" = NULL, "heatmapCustom" = NULL, "heatmapGO" = NULL,
  "pcaStatic" = NULL, "pcaPaired" = NULL,
  "boxplotStatic" = NULL, "barplotStatic" = NULL
)

output$test <- renderText({
  inputDataReactive()$dataType
})