tabItem(
  tabName = "multiDataSelectionTab",
  fluidRow(
    column(
      width = 12,
      box(
        title = "Multi-Contrast Analysis",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        tags$b(
          "Easily compare differential expression results across multiple contrasts or analysis methods!"
        ),
        tags$p(
          "Load datasets from your current project to visualise overlapping genes, pathway enrichments, and expression patterns across different comparisons. "
        ),
        tags$p(
          "If you load multiple DE results from the same FeatureCounts dataset (RNA-seq) or multiple proteomics datasets, you can create integrated boxplots to compare expression levels directly. "
        ),
        tags$p(
          "Select contrasts below to begin your integrated analysis."
        )
      ),
      shinycssloaders::withSpinner(
        uiOutput("multiDataSelectionBoxUI"),
        type = 4,
        color = "#00a65a"
      ),
      box(
        title = "Selected Datasets",
        width = 6,
        solidHeader = TRUE,
        status = "success",
        DT::dataTableOutput("multiSelectedDatasetsTable"),
        br(),
        uiOutput("multiContrastOrderBucket"),
        actionButton(
          inputId = "multiApplyOrder",
          label = "Apply Order",
          icon = icon("sort"),
          class = "btn-info"
        ),
        br(),
        br(),
        uiOutput("multiRenameInputsUI"),
        actionButton(
          inputId = "multiUpdateNames",
          label = "Update Names",
          icon = icon("check"),
          class = "btn-success"
        ),
        br(),
        br(),
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
        helpText("Upload an Excel file with 'Contrast' and 'Colour' columns to set custom colours."),
        splitLayout(
          downloadButton(outputId = "multiDlColourTemplate", label = "Template"),
          fileInput(inputId = "multiImportColourFile", accept = ".xlsx", label = NULL, placeholder = "Upload Colours")
        ),
        uiOutput(outputId = "multiColourPickerUI", inline = T)
      )
    )
  )
)
