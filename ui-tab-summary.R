use_waiter()
tabItem(
  tabName = "summaryTab",
  fluidPage(
    fluidRow(
      column(
        width = 6,
        box(
          title = paste("Settings"),
          width = 12,
          solidHeader = TRUE,
          status = "primary",
          verbatimTextOutput("test"),
          h5("Main settings"),
          uiOutput(outputId = "proteomicsContrastSelectorUI", inline = T),
          tableOutput("settingsTable"),
          br(), 
          tableOutput("summaryTable")
        )
      ),
      column(
        width = 6,
        box(
          title = paste("Download Files"),
          width = 12,
          solidHeader = TRUE,
          status = "primary",
          downloadButton("downloadInputs", "Download all input options"), br(),br(),
          downloadButton("downloadData", "Download app data"), br(),
          helpText("App data will download as a qs file, which can be imported with `qs::qread()`"), br(),
          selectizeInput(
            inputId = "downloadCount", 
            "Select Count Table", 
            choices = c()),
          downloadButton("downloadTable", "Download selected count file")
        ),
        uiOutput(outputId = "linkToMultiDEG"),
        uiOutput(outputId = "referencesText"),
      ),
      column(
        width = 12,
        box(
          title = paste("Input dataset"),
          width = 6,
          solidHeader = TRUE,
          status = "primary",
          tableOutput("inputTable")
        ),
        box(
          title = c("Group Colours"),
          width = 6, 
          solidHeader = TRUE,
          status = "primary",
          uiOutput("colourPaletteUI"),
          uiOutput(outputId = "colourPickerUI", inline = T)
        )
      )
    )
  )
)