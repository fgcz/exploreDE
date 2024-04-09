tabItem(
  tabName = "degTableTab",
  fluidPage(
    fluidRow(
      column(
        width = 12,
        box(
          title = "Results Summary",
          width = 8,
          solidHeader = TRUE,
          status = "primary",
          collapsible = TRUE,
          collapsed = FALSE,
          DT::dataTableOutput("resultsTable")
        ),
        box(
          title = "DE Calculator",
          width = 4,
          solidHeader = TRUE,
          status = "primary",
          collapsible = TRUE,
          collapsed = FALSE,
          # h4("Summary"),
          tags$p(
            "Use this little box to calculate how many features are up- or down-regulated based on given thresholds."
          ),
          br(),
          selectInput(inputId = "featureCalcPType", label = "P-value type", choices = c("FDR", "Raw"), selected = "FDR"),
          selectInput(inputId = "featureCalcPValue", label = "P-value threshold", choices = c(0.1, 0.05, 0.01, 0.001, 1e-04, 1e-05), selected = 0.05),
          sliderInput(inputId = "featureCalcLFC", label = "Log2 FC threshold (+/-)", min = 0, max = 10, value = 0.5, step = 0.25),
          br(),
          tableOutput(outputId = "featureCalcTable")
        )
      ), 
      column(
        width = 12,
        box(
          title = "Results Table",
          width = 9,
          solidHeader = TRUE,
          status = "primary",
          textOutput("DEGDesign"), br(),
          downloadButton(
            outputId = "dlDEGTable",
            label = "Download Full DE Results"),
          downloadButton(
            outputId = "dlDEGTable2",
            label = "Download Selected DE Results"),
          br(), br(), 
          tags$p(
            "You can click on any feature in this table, and they will automatically be added to the selected features buckets in volcano, heatmaps, and boxplot tabs. So you can sort for features based on name, type, TPM, log2 fold change, whatever you like, and visualise those features instantly!"),
          DT::dataTableOutput("degTable")
        ),
        box(
          title = "Feature Bucket",
          width = 3,
          solidHeader = TRUE,
          status = "primary",
          uiOutput("geneBucketDEG"),
        )
      ) 
    ) 
  ) 
) 