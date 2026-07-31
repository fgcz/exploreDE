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
          uiOutput(outputId = "proteomicsContrastSelectorUI", inline = T),
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
            # tabPanel(
            #   title = "Settings",
            #   br(),
            #   downloadButton("dlInputsQS_summary", "Download settings (qs)")
            # ),
            tabPanel(
              title = "Counts",
              br(),
              uiOutput("countsDownloadUI")
            ),
            tabPanel(
              title = "Methods",
              br(),
              uiOutput("methodsDownloadUI")
            ),
            tabPanel(
              title = "App Data",
              br(),
              h4("App data"),
              helpText("Deprecated, please use the 'Save/Load' tab."),
              downloadButton("dlData_summary", "Download app data"),
              helpText("Download all the app settings as either an Excel sheet, or as a qs file of the settings as a list."),
              downloadButton("dlInputsExcel_summary", "Download settings (Excel)"),
            )
          )
        ),
        uiOutput(outputId = "referencesText"),
        box(
          title = "Updates, August 2026",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          tags$b("More new features? Surely you can't be serious?!"),
          br(),
          tags$p("I am serious, and don't call me Shirley. ExploreDE v3.7.0 has crash landed! Otto himself got this down."),
          tags$ul(
            tags$li(
              "The main new feature is a ", tags$b("save-state function"), ". Now you can save all the settings you have carefully tweaked as a single .rds file, and then upload that again next time.",
              tags$ul(
                tags$li("It's very much a beta, and possibly quite buggy — please report any issues (and send links and .rds files too, please). One .rds file works for one exploreDE URL only — no crossing the streams here."),
                tags$li("It lives in its own tab under the now-labelled 'Session' section, called 'Save/Load Session (beta)'. Try it out!"),
                tags$li("It also doesn't save anything for the multiDE section yet; that will hopefully come shortly.")
              )
            ),
            tags$li(
              "We've also rejigged the ", tags$b("main landing page"), " of the app (i.e. just ",
              tags$a(href = "https://fgcz-shiny.uzh.ch/exploreDE", "https://fgcz-shiny.uzh.ch/exploreDE"),
              ") as we're starting to build out support for custom datasets.",
              tags$ul(
                tags$li("The spec isn't finalised yet, so documentation will follow. But very soon you'll be able to easily make an exploreDE object of your own data, wherever it's from. Exciting!")
              )
            ),
            tags$li(
              "Lastly, some KEGG/GO pathway plots have been expanded.",
              tags$ul(
                tags$li("Specifically to heatmaps and multiDE, you can now select pathways independent of statistics. Might be useful to target specific pathways.")
              )
            )
          ),
          tags$p(tags$em("As always, feedback on the app (good and bad), feature requests, and bug reports are greatly appreciated!")),
          tags$p(tags$em("P.S., if you got the film reference, then the app will run better for you."))
        )
      )
    )
  )
)