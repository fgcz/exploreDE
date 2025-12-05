tabItem(
    tabName = "multiDataSelectionTab",
    fluidRow(
        column(
            width = 12,
            uiOutput("multiDataSelectionBoxUI"),
            box(
                title = "Selected Datasets",
                width = 6,
                solidHeader = TRUE,
                status = "success",
                DT::dataTableOutput("multiSelectedDatasetsTable"),
                actionButton(
                    inputId = "multiRemoveSelected",
                    label = "Remove Selected",
                    icon = icon("trash"),
                    class = "btn-danger"
                ),
                actionButton(
                    inputId = "multiClearAll",
                    label = "Clear All",
                    icon = icon("trash-alt"),
                    class = "btn-danger"
                )
            ),
            box(
              title = "Choose colours",
              width = 6,
              solidHeader = TRUE,
              collapsible = TRUE,
              status = "success",
              uiOutput("multiColourPaletteUI"),
              uiOutput(outputId = "multiColourPickerUI", inline = T)
            )
        )
    )
)
