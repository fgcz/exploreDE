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
          uiOutput("workunitsLinkUI"),
          h5("Main settings"),
          uiOutput(outputId = "contrastSelectorUI_summary", inline = T),
          tableOutput("settingsTable"),
          br(),
          tableOutput("summaryTable"),
        ),
        # uiOutput(outputId = "linkToMultiDEG"),
        box(
          title = paste("Input dataset"),
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          helpText("The list of samples that are included in the entire dataset, and their factor columns."),
          tableOutput("inputTable")
        )
      ),
      column(
        width = 6,
        box(
          title = paste("Download Files"),
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          tabsetPanel(
            tabPanel(
              title = "Counts",
              br(),
              uiOutput("countsDownloadUI")
            ),
            tabPanel(
              title = "Methods",
              br(),
              uiOutput("methodsDownloadUI")
            )
          )
        ),
        uiOutput(outputId = "referencesText"),
        box(
          title = "Updates, October 2026",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          tags$b("Yer jokin'! Not another one?!"),
          br(),
          tags$p("ExploreDE 3.9.0! This app gets more updates than the UK does prime ministers.", tags$sub(tags$em("[Spain enters the chat]"), style = "font-size: 60%;")),
          tags$ul(
            tags$li(
              tags$b("Finally!"), " You can relabel groups!",
              tags$ul(
                tags$li("Under the 'Colours / Label' tab, you can rename the group labels to whatever you want!"),
                tags$li("It supports italics, super/subscript, and Greek letters..."),
                tags$li("... mostly! They might not all show in the Shiny app, but the PDFs should render fine.")
              )
            ),
            tags$li(
              "We've also spent some time tidying some plots, fixing some little niggles (like the feature bucket), and adding some extra little sliders."
            )
          ),
          tags$p(tags$em("As always, feedback on the app (good and bad), feature requests, and bug reports are greatly appreciated!"))
        )
      )
    )
  )
)