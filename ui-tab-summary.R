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
          h5("New! (November 2024)"),
          helpText(
            "Now you can more easily set custom colours without having to choose them each time you load the app. 
            Download the colour template Excel file and change the colours for each factor level in the second column. 
            You can use hex codes (e.g., #750912, #59a7cf) or ggplot colour names (e.g., dodgerblue, plum4, etc.).
            Then, upload that file and watch the colour palette and your figures automagically update to the new colours!
            Don't change the values in the first column, or it won't work!
            Hopefully this is a much easier way to set colours once."
          ),
          splitLayout(
            downloadButton(outputId = "dlColourTemplate", label = paste("Download Colour Template")),
            fileInput(inputId = "importColourFile", accept = ".xlsx", width = "85%", label = NULL, placeholder = "Upload Colour Template")
          ),
          uiOutput("colourPaletteUI"),
          uiOutput(outputId = "colourPickerUI", inline = T)
        )
      )
    )
  )
)