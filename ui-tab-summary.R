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
          # h5("Threshold"),
          tableOutput("summaryTable")
        ) # end of box
      ),
      column(
        width = 6,
        box(
          title = paste("Download Count Files"),
          width = 12,
          solidHeader = TRUE,
          status = "primary",
          selectizeInput(
            inputId = "downloadCount", 
            "Select Count Table", 
            choices = c()),
          downloadButton("downloadTable", "Download Selected Count File")
        ), # end of box
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
        ), # end of box
        box(
          title = c("Group Colours"),
          width = 6, 
          solidHeader = TRUE,
          status = "primary",
          selectInput(
            inputId = "colourPalette", 
            label = "Colour Palette",
            choices = c("House colours", "Paired", "Set1", "Set2", "Set3", "Dark2"), 
            selected = "House colours", 
            multiple = TRUE
            ),
          uiOutput(outputId = "colourPickerUI", inline = T)
        )
      ) # end of column
    )
  ) # end of fluid page
) # end of tabItem: parameters