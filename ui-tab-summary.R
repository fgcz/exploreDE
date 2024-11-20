use_waiter()
tabItem(
  tabName = "summaryTab",
  fluidPage(
    fluidRow(
      column(
        width = 6,
        box(
          title = paste("Settings"),
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          verbatimTextOutput("test"),
          h5("Main settings"),
          uiOutput(outputId = "proteomicsContrastSelectorUI", inline = T),
          tableOutput("settingsTable"),
          br(), 
          tableOutput("summaryTable"),
        ),
        uiOutput(outputId = "linkToMultiDEG")
      ),
      column(
        width = 6,
        box(
          title = paste("Download Files"),
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          h4("Download App Data"),
          helpText("App data will download as a qs file, which can be imported into R with `qs::qread()`"),
          downloadButton("downloadData", "Download app data"),
          hr(style = "border-top: 1px solid #000000;"), h4("Download App Settings"),
          helpText("Download all the app settings as either an Excel sheet, or as a qs file of the settings as a list."),
          downloadButton("downloadInputsE", "Download settings (Excel)"), downloadButton("downloadInputsR", "Download settings (qs)"),
          hr(style = "border-top: 1px solid #000000;"), h4("Download Count Files"),
          helpText("Select the type of count data you wish to download for all samples in the dataset."),
          selectizeInput(
            inputId = "downloadCount", 
            "Select Count Table", 
            choices = c(), 
            width = "33%"
            ),
          downloadButton("downloadTable", "Download selected count file")
        ),
        uiOutput(outputId = "referencesText"),
      ),
      column(
        width = 12,
        box(
          title = paste("Input dataset"),
          width = 6,
          solidHeader = TRUE,
          status = "primary",
          helpText("The list of samples that are included in the entire dataset, and their factor columns."),
          tableOutput("inputTable")
        ),
        box(
          title = c("Group Colours"),
          width = 6,
          solidHeader = TRUE,
          status = "primary",
          downloadButton(outputId = "dlColourTemplate", label = paste("Download Colour Template")),
          fileInput(inputId = "importColourFile", accept = ".xlsx", width = "85%", label = "Upload Colour Template"),
          uiOutput("colourPaletteUI"),
          uiOutput(outputId = "colourPickerUI", inline = T)
        )
      )
    )
  )
)