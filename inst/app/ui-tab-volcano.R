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
          selected for the plot."
        ),
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
        uiOutput(outputId = "contrastSelectorUI_volcano"),
        tabsetPanel(
          tabPanel(
            title = "Main settings",
            numericInput(
              inputId = "lfcVolcano",
              label = "Log2 Fold Change Threshold:",
              value = 0.5,
              min = 0,
              max = 20,
              step = 0.25
            ),
            selectInput(inputId = "pTypeVolcano", label = "P-Value", choices = c("FDR", "Raw"), selected = "FDR"),
            selectInput(
              inputId = "pThresholdVolcano",
              label = "P-value Threshold",
              choices = c(1.0, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 1e-6),
              selected = 0.05
            ),
            uiOutput(outputId = "nrPeptidesVolcanoUI"),
            numericInput(
              inputId = "xLimVolcano",
              label = "x axis (log2 FC) limits (both +/-)",
              value = 5,
              min = 1,
              max = 20,
              step = 0.5
            ),
            numericInput(
              inputId = "yLimVolcano",
              label = "y axis (-log 10 p-value) limits",
              value = 50,
              min = 10,
              max = 250,
              step = 10
            ),
            uiOutput(outputId = "showImputedVolcanoUI"),
            hr(style = "border-top: 1px solid #000000;"),
            h4("Volcano plot annotation"),
            uiOutput(outputId = "highlightImputedVolcanoUI"),
            h5("Use this checkbox to turn on/off all annotations:"),
            checkboxInput(inputId = "volcanoShowGenes", label = "Annotate volcano?", value = TRUE),
            uiOutput("fb_label_mode_volcano"),
            h5("Label all up and/or down regulated features (based on current settings)?"),
            checkboxInput(inputId = "volcanoLabelAllUp", label = "Label all Up", value = FALSE),
            checkboxInput(inputId = "volcanoLabelAllDown", label = "Label all Down", value = FALSE),
            helpText("You can both label all up/down features and select specific features."),
            tags$b("Annotate volcano plot with specific features?"),
            uiOutput("fb_search_mode_volcano"),
            selectizeInput(
              inputId = "fb_select_volcano",
              label = "Select features to annotate",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list(
                placeholder = 'Select features',
                plugins = list('remove_button', 'drag_drop', 'restore_on_backspace', 'clear_button')
              )
            ),
            textAreaInput(
              inputId = "fb_text_volcano",
              label = "Or, paste list of features:",
              placeholder = "feature1 feature2 feature3 feature4",
              cols = 1
            ),
            hr(style = "border-top: 1px solid #000000;"),
            h4("Feature Bucket"),
            h5(
              "Features will be added to this bucket as you select them from the DE table tab and from the inputs in various tabs. You can use this bucket to quickly 
              include/exclude these features from highlighting as you need to by dragging and dropping them in order or into the exclude bucket."
            ),
            uiOutput("fb_bucket_volcano"),
            actionButton(inputId = "fb_add_brushed", label = "Add brushed genes to bucket", icon = icon("arrow-right")),
            actionButton(inputId = "fb_reset_volcano", label = "Empty the bucket?", icon = icon("bucket")),
            hr(style = "border-top: 1px solid #000000;"),
            h4("Volcano colours"),
            checkboxInput(
              inputId = "volcanoAnnotationHighlightColour",
              label = "Use special highlight colour?",
              value = F
            ),
            uiOutput("volcanoColourPicker")
          ),
          tabPanel(
            title = "Figure settings",
            h4("Plot settings"),
            splitLayout(
              sliderInput(
                inputId = "figWidthVolcano",
                label = "Width",
                min = 100,
                max = 2000,
                value = 700,
                step = 10,
                width = "85%",
                ticks = FALSE
              ),
              sliderInput(
                inputId = "figHeightVolcano",
                label = "Height",
                min = 100,
                max = 2000,
                value = 480,
                step = 10,
                width = "85%",
                ticks = FALSE
              ),
              sliderInput(
                inputId = "textSizeVolcano",
                label = "Font Size",
                min = 4,
                max = 30,
                value = 14,
                step = 0.5,
                width = "85%",
                ticks = FALSE
              )
            ),
            checkboxInput(inputId = "volcanoHaloStyle", label = "Halo-style dot border (new)?", value = TRUE),
            splitLayout(
              sliderInput(
                inputId = "dotSizeVolcano",
                label = "Dot size",
                min = 1,
                max = 10,
                value = 1.5,
                step = 0.5,
                width = "85%",
                ticks = FALSE
              ),
              sliderInput(
                inputId = "alphaVolcano",
                label = "Dot alpha",
                min = 0.1,
                max = 1,
                value = 1.0,
                step = 0.1,
                width = "85%",
                ticks = FALSE
              ),
              tagList(
                conditionalPanel(
                  condition = "input.volcanoHaloStyle == true",
                  sliderInput(
                    inputId = "volcanoHaloMult",
                    label = "Halo thickness (× dot size)",
                    min = 0,
                    max = 3,
                    value = 1,
                    step = 0.1,
                    width = "85%",
                    ticks = FALSE
                  )
                ),
                conditionalPanel(
                  condition = "input.volcanoHaloStyle == false",
                  sliderInput(
                    inputId = "volcanoPointBorder",
                    label = "Dot border",
                    min = 0,
                    max = 1,
                    value = 0.2,
                    step = 0.1,
                    width = "85%",
                    ticks = FALSE
                  )
                )
              )
            ),
            helpText(
              "Halo style draws a black dot underneath each point, larger by this multiplier, so only the rim shows through — keeps the plot clean even where points pile up densely. Works best with Dot alpha near 1."
            ),
            br(),
            splitLayout(
              checkboxInput(inputId = "showAxesVolcano", label = "Axes lines?", value = TRUE),
              checkboxInput(inputId = "borderVolcano", label = "Plot Border?", value = TRUE),
              checkboxInput(inputId = "showLinesVolcano", label = "Grid lines?", value = FALSE)
            ),
            hr(style = "border-top: 1px solid #000000;"),
            h4("Quadrant counts"),
            checkboxInput(
              inputId = "volcanoShowCounts",
              label = "Show n = for the top-left/top-right quadrants?",
              value = TRUE
            ),
            splitLayout(
              sliderInput(
                inputId = "volcanoCountsX",
                label = "Position X (log2 FC)",
                min = 0,
                max = 7.5,
                value = 3,
                step = 0.05,
                width = "85%",
                ticks = FALSE
              ),
              sliderInput(
                inputId = "volcanoCountsY",
                label = "Position Y (y axis)",
                min = 0,
                max = 75,
                value = 52.5,
                step = 0.5,
                width = "85%",
                ticks = FALSE
              ),
              sliderInput(
                inputId = "volcanoCountsSize",
                label = "Font Size",
                min = 4,
                max = 30,
                value = 14,
                step = 0.5,
                width = "85%",
                ticks = FALSE
              )
            ),
            helpText(
              "Counts are the features passing both the log2 fold change and p-value thresholds (down on the left, up on the right). Positions are in the plot's own axis units and can go past the axis limits (default: above the dots); the plot extends to fit them while the axis ticks stay within the limits. X is mirrored for the left count."
            ),
            hr(style = "border-top: 1px solid #000000;"),
            h4("Annotation settings"),
            numericInput(
              inputId = "volcanoLabelMaxOverlap",
              label = "Number of max overlapping labels",
              min = 1,
              max = 1e4,
              value = 10,
              step = 1,
              width = "33%"
            ),
            helpText(
              "Increasing the number of overlapping labels will label more features, but can take a *very* long time to generate"
            ),
            splitLayout(
              sliderInput(
                inputId = "geneLabelNudgeVolcanoX",
                label = "Nudge Labels X",
                min = -10,
                max = 10,
                value = 0,
                step = 1,
                width = "85%",
                ticks = FALSE
              ),
              sliderInput(
                inputId = "geneLabelNudgeVolcanoY",
                label = "Nudge Labels Y",
                min = -10,
                max = 10,
                value = 0,
                step = 1,
                width = "85%",
                ticks = FALSE
              ),
              sliderInput(
                inputId = "geneLabelSizeVolcano",
                label = "Label Size",
                min = 4,
                max = 30,
                value = 12,
                step = 0.5,
                width = "85%",
                ticks = FALSE
              )
            )
          ),
          tabPanel(
            title = "Download",
            h4("Download Settings"),
            splitLayout(
              selectInput(
                inputId = "downloadFormatVolcano",
                label = "Select format",
                choices = c("PDF", "SVG", "PNG"),
                selected = "PDF"
              ),
              selectInput(
                inputId = "dpiVolcano",
                label = "PNG DPI",
                choices = c(72, 150, 300, 600, 1000),
                selected = 600
              )
            ),
            textInput(inputId = "filename_volcano", label = "Enter filename", value = "Volcano"),
            downloadButton(outputId = "dlPlot_volcano", HTML("Download Volcano Plot")),
            downloadButton(outputId = "dlPlot_MA", HTML("Download MA Plot")),
            downloadButton(outputId = "dlTable_volcano", HTML("Download Volcano Results Table (Excel)"))
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
        tabsetPanel(
          tabPanel(
            title = "Static",
            helpText(
              "Click and drag on the volcano plot to select features. Brushed genes are annotated on the plot immediately. Use the 'Add brushed genes to bucket' button to push them to the shared feature bucket for boxplots, heatmaps, etc."
            ),
            textOutput("volcanoDesign"),
            plotOutput("volcanoStatic", inline = TRUE, brush = "volcanoBrush"),
            br(),
            br(),
            DT::dataTableOutput("volcanoBrushTable")
          ),
          tabPanel(
            title = "Interactive",
            plotlyOutput("volcanoPlotly", inline = TRUE, height = "auto")
          )
        )
      ),
      box(
        title = "MA Plot",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        plotOutput("MAPlot", inline = TRUE, brush = "MABrush"),
        DT::dataTableOutput("MABrushTable")
      )
    )
  )
)
