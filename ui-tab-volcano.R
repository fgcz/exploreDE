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
        tabsetPanel(
          tabPanel(
            title = "Main settings",
            numericInput(inputId = "lfcVolcano", label = "Log2 Fold Change Threshold:", value = 0.5, min = 0, max = 20, step = 0.25),
            selectInput(inputId = "pTypeVolcano", label = "P-Value", choices = c("FDR", "Raw"), selected = "FDR"),
            selectInput(inputId = "pThresholdVolcano", label = "P-value Threshold", choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 1e-6), selected = 0.05),
            uiOutput(outputId = "nrPeptidesVolcanoUI"),
            numericInput(inputId = "xLimVolcano", label = "x axis (log2 FC) limits (both +/-)", value = 5, min = 1, max = 20, step = 0.5),
            numericInput(inputId = "yLimVolcano", label = "y axis (-log 10 p-value) limits", value = 50, min = 10, max = 250, step = 10),
            uiOutput(outputId = "showImputedVolcanoUI"),
            hr(style = "border-top: 1px solid #000000;"), h4("Volcano plot annotation"),
            uiOutput(outputId = "highlightImputedVolcanoUI"),
            h5("Use this checkbox to turn on/off all annotations:"),
            checkboxInput(inputId = "volcanoShowGenes", label = "Annotate volcano?", value = TRUE),
            h5("Label all up and/or down regulated features (based on current settings)?"),
            checkboxInput(inputId = "volcanoLabelAllUp", label = "Label all Up", value = FALSE),
            checkboxInput(inputId = "volcanoLabelAllDown", label = "Label all Down", value = FALSE),
            helpText("You can both label all up/down features and select specific features."),
            tags$b("Annotate volcano plot with specific features?"),
            selectizeInput(inputId = "volcanoGenes", label = "Select features to annotate", choices = "", selected = "", multiple = TRUE),
            textAreaInput(inputId = "volcanoGenesText", label = "Or, paste list of features:", placeholder = "feature1 feature2 feature3 feature4", cols = 1),
            hr(style = "border-top: 1px solid #000000;"), h4("Feature Bucket"),
            h5(
              "Features will be added to this bucket as you select them from the DE table tab and from the inputs in various tabs. You can use this bucket to quickly 
              include/exclude these features from highlighting as you need to by dragging and dropping them in order or into the exclude bucket."),
            uiOutput("geneBucket2"),
            hr(style = "border-top: 1px solid #000000;"), h4("Volcano colours"),
            checkboxInput(inputId = "volcanoAnnotationHighlightColour", label = "Use special highlight colour?", value = F),
            uiOutput("volcanoColourPicker")
          ),
          tabPanel(
            title = "Figure settings",
            h4("Plot settings"),
            splitLayout(
              sliderInput(inputId = "figWidthVolcano", label = "Width", min = 100, max = 2000, value = 800, step = 10, width = "85%"),
              sliderInput(inputId = "figHeightVolcano", label = "Height", min = 100, max = 2000, value = 600, step = 10, width = "85%"),
              sliderInput(inputId = "textSizeVolcano", label = "Font Size", min = 4, max = 30, value = 12, step = 0.5, width = "85%", ticks = TRUE)
            ),
            splitLayout(
              sliderInput(inputId = "dotSizeVolcano", label = "Dot size", min = 1, max = 10, value = 3, step = 0.5, width = "85%", ticks = TRUE),
              sliderInput(inputId = "alphaVolcano", label = "Dot alpha", min = 0.1, max = 1, value = 0.9, step = 0.1, width = "85%", ticks = TRUE),
              sliderInput(inputId = "volcanoPointBorder", label = "Dot border", min = 0, max = 1, value = 0.2, step = 0.1, width = "85%"),
            ),
            br(),
            splitLayout(
              checkboxInput(inputId = "showAxesVolcano", label = "Axes lines?", value = TRUE),
              checkboxInput(inputId = "boldVolcano", label = "Bold?", value = TRUE),
              checkboxInput(inputId = "showLinesVolcano", label = "Grid lines?", value = TRUE)
            ),
            hr(style = "border-top: 1px solid #000000;"), h4("Annotation settings"),
            numericInput(inputId = "volcanoLabelMaxOverlap", label = "Number of max overlapping labels", min = 1, max = 1e4, value = 10, step = 1, width = "33%"),
            helpText("Increasing the number of overlapping labels will label more features, but can take a *very* long time to generate"),
            splitLayout(
              sliderInput(inputId = "geneLabelNudgeVolcanoX", label = "Nudge Labels X", min = -10, max = 10, value = 0, step = 1, width = "85%", ticks = TRUE),
              sliderInput(inputId = "geneLabelNudgeVolcanoY", label = "Nudge Labels Y", min = -10, max = 10, value = 0, step = 1, width = "85%", ticks = TRUE),
              sliderInput(inputId = "geneLabelSizeVolcano", label = "Label Size", min = 4, max = 30, value = 12, step = 0.5, width = "85%", ticks = TRUE)
            )
          ),
          tabPanel(
            title = "Download settings",
            selectInput(inputId = "downloadFormatVolcano", label = "Select format", choices = c("PDF", "SVG", "PNG"), selected = "PDF"),
            selectInput(inputId = "dpiVolcano", label = "PNG DPI", choices = c(72, 150, 300, 600, 1000), selected = 600),
            textInput(inputId = "filnameVolcano", label = "Enter filename", value = "Volcano")
          )
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
        DT::dataTableOutput("volcanoOverviewTable")
      ),
      box(
        title = "Volcano Plot",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        h5("New feature alert! (February 2024)"),
        helpText("If you click and drag on the volcano plot to select features, they will automatically be added to your feature bucket in other tabs too, e.g., boxplots, heatmaps."),
        textOutput("volcanoDesign"),
        downloadButton(outputId = "dlVolcanoPlotButton", HTML("Download Volcano Plot")),
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