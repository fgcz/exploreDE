tabItem(
  tabName = "parametersTab",
  fluidPage(
    fluidRow(
      column(
        width = 3,
        box(
          title = paste("Colours:"),
          width = NULL,
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
        ) # end of box
      ) # end of column
    )
  ) # end of fluid page
) # end of tabItem: parameters