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
        br(),
        div(
          style = "display: inline-block; vertical-align: top;",
          textInput(inputId = "filename_multiTable", label = NULL, value = "Multi_DE_Table", placeholder = "Filename")
        ),
        div(
          style = "display: inline-block; vertical-align: top; margin-top: 0px;",
          downloadButton(outputId = "dlTable_multiDE", label = "Download Table")
        )
      )
    )
  )
)
