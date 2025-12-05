tabItem(
  tabName = "multiDETableTab",
  fluidRow(
    column(
      width = 12,
      box(
        title = "Integrated DE Table",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        helpText("Combined differential expression results from all loaded contrasts. Search and filter across all datasets."),
        DT::dataTableOutput("multiDETable"),
        downloadButton(outputId = "dlMultiDETable", label = "Download Table")
      )
    )
  )
)
