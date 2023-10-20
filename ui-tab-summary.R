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
          # Colour picker for each of the groups in Condition:
          lapply(1:50, function(i) {
            colourpicker::colourInput(
              inputId = paste0("GroupColour", i),
              label = "",
              value =  rep(c(
                "indianred", "steelblue", "chartreuse4", "gray30", 
                "goldenrod3", "indianred4", "royalblue4", "mediumorchid3",
                "turquoise4", "darkolivegreen", "thistle4", "darkorange3", 
                "hotpink2", "burlywood3", "cadetblue4", "chocolate4", "firebrick"
              ), times = 5)[i],
              palette = "square",
              closeOnClick = TRUE,
              returnName = TRUE
            )
          })
        )
      ) # end of column
    )
  ) # end of fluid page
) # end of tabItem: parameters