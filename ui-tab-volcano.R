tabItem(
  tabName = "volcanoTab",
  fluidRow(
    column(
      width = 4,
      box(
        title = "Info",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        collapsed = TRUE,
        tags$p(
          "Volcano plots show the p-value and log fold change results
          from the differential expression test. The p-value is first
          converted onto log10 scale (to maximise the distance between
          very many, very low p-values), then negated so that the 
          lowest p-values are now at the top of the axis."
        ),
        tags$p(
          "Specific features can be highlighted either by selecting them from the
           drop down or pasting the feature name the box. Highlighted features will be
           marked as red dots irrespective of p-value and log fold change 
           values"
        ),
        tags$p(
          "The Results Summary tab shows the grouping of each feature
          based on the p-value and log fold change parameters you have
          selected for the plot."),
        tags$p(
          "The download button will download a PDF version of the plot 
          exactly as it is displayed with the settings you have chosen."
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        h4("Volcano plot settings"),
        tabsetPanel(
          tabPanel(
            title = "Main settings",
            numericInput(inputId = "lfcVolcano", label = "Log2 Fold Change Threshold:", value = 0.5, min = 0, max = 20, step = 0.25),
            selectInput(inputId = "pTypeVolcano", label = "P-Value", choices = c("FDR", "Raw"), selected = "FDR"),
            selectInput(inputId = "pThresholdVolcano", label = "P-value Threshold", choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 1e-6), selected = 0.05),
            numericInput(inputId = "xLimVolcano", label = "x axis (log2 FC) limits (both +/-)", value = 5, min = 1, max = 20, step = 0.5),
            numericInput(inputId = "yLimVolcano", label = "y axis (-log 10 p-value) limits", value = 50, min = 10, max = 250, step = 10),
            hr(style = "border-top: 1px solid #000000;"), h4("Volcano plot annotation"),
            h5("Use this checkbox to turn on/off all annotations:"),
            checkboxInput(inputId = "volcanoShowGenes", label = "Annotate volcano?", value = TRUE),
            h5("Label all up and/or down regulated features (based on current settings)?"),
            checkboxInput(inputId = "volcanoLabelAllUp", label = "Label all Up", value = FALSE),
            checkboxInput(inputId = "volcanoLabelAllDown", label = "Label all Down", value = FALSE),
            helpText("You can both label all up/down features and select specific features."),
            tags$b("Annotate volcano plot with specific features?"),
            selectizeInput(inputId = "volcanoGenes", label = "Select features to annotate", choices = "", selected = "", multiple = TRUE),
            textAreaInput(inputId = "volcanoGenesText", label = "Or, paste list of features:", placeholder = "feature1 feature2 feature3 feature4", cols = 1),
            h5(
              "Features will be added to this bucket as you select them from the DE table 
              tab and from the inputs in this tab. You can use this bucket to quickly 
              re-order and exclude these features as you need to by dragging and dropping 
              them in order or into the exclude bucket."),
            uiOutput("geneBucket2"),
            hr(style = "border-top: 1px solid #000000;"), h4("Volcano colours"),
            checkboxInput(inputId = "volcanoAnnotationHighlightColour", label = "Use special highlight colour?", value = F),
            uiOutput("volcanoColourPicker")
          ),
          tabPanel(
            title = "Figure settings",
            splitLayout(
              checkboxInput(inputId = "showBorderVolcano", label = "Show cell border?", value = TRUE),
              checkboxInput(inputId = "showAxesVolcano", label = "Show axes lines?", value = TRUE),
              checkboxInput(inputId = "boldVolcano", label = "Use bold font?", value = TRUE),
              checkboxInput(inputId = "showLinesVolcano", label = "Show grid lines?", value = TRUE)
            ),
            
            sliderInput(inputId = "dotSizeVolcano", label = "Dot size", min = 1, max = 10, value = 3, step = 0.5, width = "33%", ticks = TRUE),
            sliderInput(inputId = "alphaVolcano", label = "Point alpha", min = 0.1, max = 1, value = 0.9, step = 0.1, width = "33%", ticks = TRUE),
            numericInput(inputId = "volcanoLabelMaxOverlap", label = "Number of max overlapping labels", min = 1, max = 1e4, value = 10, step = 1, width = "33%"),
            helpText("Increasing the number of overlapping labels will label more genes, but can take a *very* long time to generate"),
            sliderInput(inputId = "geneLabelNudgeVolcanoX", label = "Nudge Gene Labels X", min = -10, max = 10, value = 0, step = 1, width = "33%", ticks = TRUE),
            sliderInput(inputId = "geneLabelNudgeVolcanoY", label = "Nudge Gene Labels Y", min = -10, max = 10, value = 0, step = 1, width = "33%", ticks = TRUE),
            sliderInput(inputId = "geneLabelSizeVolcano", label = "Gene Label Size", min = 4, max = 30, value = 12, step = 0.5, width = "33%", ticks = TRUE),
            sliderInput(inputId = "textSizeVolcano", label = "Figure Font Size", min = 4, max = 30, value = 12, step = 0.5, width = "33%", ticks = TRUE),
            numericInput(inputId = "figWidthVolcano", label = "Figure Width", min = 100, max = 2000, value = 800, step = 10),
            numericInput(inputId = "figHeightVolcano", label = "Figure Height", min = 100, max = 2000, value = 600, step = 10)
          ),
        )
      )
    ),
    column(
      width = 8,
      box(
        title = "Results Summary",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        collapsed = TRUE,
        # tableOutput(outputId = "volcanoOverview"),
        # br(),
        DT::dataTableOutput("volcanoOverviewTable")
      ),
      box(
        title = "Volcano Plot",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        textOutput("volcanoDesign"), br(),
        downloadButton(outputId = "dlVolcanoPlotButton", HTML("Download Volcano Plot (PDF)")),
        downloadButton(outputId = "dlVolcanoDFButton", HTML("Download Volcano Results Table (Excel)")),
        br(), br(),
        plotOutput("volcanoStatic", inline = TRUE, brush = "volcanoBrush"),
        br(), br(),
        DT::dataTableOutput("volcanoBrushTable")
      ),
      box(
        title = "MA Plot", 
        width = NULL,
        solidHeader = TRUE, 
        status = "primary", 
        downloadButton(outputId = "dlMAPlotButton", HTML("Download MA Plot (PDF)")),
        br(), br(),
        plotOutput("MAPlot", inline = TRUE, brush = "MABrush"),
        DT::dataTableOutput("MABrushTable")
      )
    ) 
  ) 
)