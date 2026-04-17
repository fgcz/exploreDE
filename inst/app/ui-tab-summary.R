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
          title = "Updates, April 2026",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          tags$b("Your favourite shiny app continues to get new features, how fab!"),
          br(), br(),
          tags$p("Say hello to explodeDE v 3.6.0! What treats have we added for you this time?"),
          tags$p("- The app now supports Mageck counts from Crispr screens! So now you can boxplot your favourite guides, and hope they expanded."),
          tags$p("- Batch correction has been modified to include the ability to pass the design. In DE tests, this will be defaulted to the column tested. If you applied a batch correction, this will be passed as default. It'll help make batch correction "),
          tags$p("- Proteomics: You can now select features using the gene symbol. So, if you have peptide level data and select geneA, you'll see all the peptides for that gene."),
          tags$p("- Genomics: A whole bunch of QC plots have been brough over the CountQC and DE static reports, with two new tabs. Play around!"),
          tags$p("As always, feedback on the app (good and bad), feature requests, and bug reports are greatly appreciated!")
        )
      )
    )
  )
)