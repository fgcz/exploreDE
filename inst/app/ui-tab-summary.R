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
          h4("Download App Data"),
          helpText("App data will download as a qs file, which can be imported into R with `qs::qread()`"),
          downloadButton("dlData_summary", "Download app data"),
          hr(style = "border-top: 1px solid #000000;"), h4("Download App Settings"),
          helpText("Download all the app settings as either an Excel sheet, or as a qs file of the settings as a list."),
          downloadButton("dlInputsExcel_summary", "Download settings (Excel)"), downloadButton("dlInputsQS_summary", "Download settings (qs)"),
          uiOutput("methodsDownloadUI"),
          uiOutput("countsDownloadUI")
        ),
        uiOutput(outputId = "referencesText"),
        box(
          title = "Updates, December 2025",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          tags$b("Exciting new updates to your favourite shiny app, to help you have even more fun with your data into 2026!"),
          br(), br(),
          tags$p("Introducing explodeDE v 3.5.0! What's new?"),
          tags$p("- App now supports CountQC data and more proteomics data types! As such, the tabs have been rejigged."),
          tags$p("- MultiDE has been subsumed! You can now select other DE results from your project directly in exploreDE."),
          tags$p("- PCA has been reborn as ordination, and now supports MDS, NMDS, and my favourite new feature, biplot mode!"),
          tags$p("- A new 'Replicability' tab has been added to show sample-to-sample correlations and replicate consistency."),
          tags$p("- In the heatmaps tab, you can now group the samples by a factor column, and plot the mean or median per group."),
          helpText("It's hidden right at the bottom of the settings box..."),
          tags$p("- In the volcano tab, if you select a feature to highlight, but then drop it into the exclude bucket, it'll remove the feature from the plot entirely."),
          tags$p("- In the ORA and GSEA tabs, you can now filter which genes are highlighted/displayed in the figures. The box is hidden at the bottom of the page."),
          tags$p("- Now, in places where you can download the feature counts (summary, boxplot, heatmap), you can choose to append the DE results (RNA seq only for now.)"),
          tags$p("Most of these new features are here because of user requests, so please keep the feedback coming!")
        )
      )
    )
  )
)