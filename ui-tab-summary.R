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
        width = 4,
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
        box(
          title = "Compare Multiple DE Tests!",
          width = 12,
          solidHeader = TRUE,
          status = "primary",
          tags$p("Do you want to compare the results in this DE test with other DE tests you have run? Then, good news, because there is an app for just that! Head to https://fgcz-shiny.uzh.ch/multiDEG to find out more!"),
          h4("Try the MultiDEG app:", a("MultiDEG", href="https://fgcz-shiny.uzh.ch/multiDEG/", target = "_blank"))
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
          tableOutput("inputTable")
        ), # end of box
        box(
          title = c("Group Colours"),
          width = 4, 
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