tabItem(
  tabName = "deQcTab",
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
          "The DE QC tab provides diagnostic plots for the differential
          expression analysis. These help assess whether the test was
          well calibrated and whether results are confounded by technical
          features such as GC content or gene length."
        ),
        tags$p(
          tags$b("Between-group scatter:"),
          " mean expression of the reference group (x) against the sample
          group (y) on a log10 scale. Points are coloured by expression
          status and significance. Blue lines mark 1:1 and the 2-fold
          fold-change thresholds."
        ),
        tags$p(
          tags$b("P-value distribution:"),
          " a well-calibrated test produces a uniform distribution of raw
          p-values with a spike near zero. The dashed line marks the
          expected count per bin under the null. Modes elsewhere suggest
          model misspecification."
        ),
        tags$p(
          tags$b("GC / length bias:"),
          " scatter of log2 fold change against GC content (or log10 gene
          width), with an optional LOESS trend. If the trend deviates from
          y = 0, the DE results may be confounded by that technical feature."
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          id = "deQcSettingsTabs",
          tabPanel(
            title = "Between-group",
            # Proteomics-only: pick the contrast + count type + group pair
            uiOutput("deQcProteomicsSelectorsUI"),
            radioButtons(
              inputId = "deQcScatterColourBy",
              label = "Colour by:",
              choices = c("Standard", "GC content", "Gene width"),
              selected = "Standard", inline = TRUE
            ),
            helpText("GC and Gene-width modes require 'gc' and 'featWidth' in
                     seqAnno (RNASeq only)."),
            conditionalPanel(
              condition = "input.deQcScatterColourBy == 'Standard'",
              numericInput(
                inputId = "deQcScatterPval",
                label = "P-value threshold for 'Significants'",
                value = 0.05, min = 0, max = 1, step = 0.01
              ),
              numericInput(
                inputId = "deQcScatterLFC",
                label = "Log2 fold-change threshold for 'Significants'",
                value = 1, min = 0, max = 10, step = 0.25
              ),
              splitLayout(
                colourpicker::colourInput("deQcColAbsent", "Absent", value = "grey70"),
                colourpicker::colourInput("deQcColPresent", "Present", value = "black"),
                colourpicker::colourInput("deQcColSig", "Significants", value = "#A6CEE3"),
                cellArgs = list(style = "overflow: visible;")
              )
            ),
            conditionalPanel(
              condition = "input.deQcScatterColourBy == 'GC content'",
              splitLayout(
                numericInput("deQcScatterGCLow",  "GC low",  value = 0.4,
                             min = 0, max = 1, step = 0.05, width = "85%"),
                numericInput("deQcScatterGCHigh", "GC high", value = 0.6,
                             min = 0, max = 1, step = 0.05, width = "85%"),
                cellArgs = list(style = "overflow: visible;")
              ),
              splitLayout(
                colourpicker::colourInput("deQcColGCLow",  "GC low bin",  value = "#1F78B4"),
                colourpicker::colourInput("deQcColGCMid",  "GC mid bin",  value = "black"),
                colourpicker::colourInput("deQcColGCHigh", "GC high bin", value = "#FF7F00"),
                cellArgs = list(style = "overflow: visible;")
              )
            ),
            conditionalPanel(
              condition = "input.deQcScatterColourBy == 'Gene width'",
              splitLayout(
                numericInput("deQcScatterWidthLow",  "Width low (nt)",  value = 500,
                             min = 0, step = 100, width = "85%"),
                numericInput("deQcScatterWidthHigh", "Width high (nt)", value = 5000,
                             min = 0, step = 100, width = "85%"),
                cellArgs = list(style = "overflow: visible;")
              ),
              splitLayout(
                colourpicker::colourInput("deQcColWidthLow",  "Width low bin",  value = "#B15928"),
                colourpicker::colourInput("deQcColWidthMid",  "Width mid bin",  value = "black"),
                colourpicker::colourInput("deQcColWidthHigh", "Width high bin", value = "#33A02C"),
                cellArgs = list(style = "overflow: visible;")
              )
            ),
            hr(),
            splitLayout(
              sliderInput("deQcScatterPointSize", "Point size",
                          min = 0.1, max = 5, value = 1, step = 0.1,
                          width = "85%", ticks = FALSE),
              sliderInput("deQcScatterAlpha", "Alpha",
                          min = 0.1, max = 1, value = 0.6, step = 0.05,
                          width = "85%", ticks = FALSE),
              cellArgs = list(style = "overflow: visible;")
            ),
            hr(),
            h5("Guide lines"),
            helpText("Solid line: 1:1 identity. Dashed lines: N-fold change bands."),
            splitLayout(
              colourpicker::colourInput("deQcColGuide", "Guide colour",
                                        value = "blue", width = "85%"),
              numericInput("deQcGuideFold", "Fold-change threshold",
                           value = 2, min = 1, max = 20, step = 0.5, width = "85%"),
              cellArgs = list(style = "overflow: visible;")
            ),
            splitLayout(
              checkboxInput("deQcShowIdentity", "Show 1:1 line", value = TRUE),
              checkboxInput("deQcShowFoldBands", "Show fold-change bands", value = TRUE)
            )
          ),
          tabPanel(
            title = "P-value hist",
            sliderInput(
              inputId = "deQcPhistBins",
              label = "Number of bins", min = 50, max = 500, value = 200, step = 10,
              width = "85%", ticks = FALSE
            ),
            radioButtons(
              inputId = "deQcPhistType",
              label = "P-value type:",
              choices = c("Raw" = "pValue", "FDR" = "fdr"),
              selected = "pValue"
            ),
            splitLayout(
              colourpicker::colourInput("deQcColPhistUsed", "Used in test", value = "#1965B0"),
              colourpicker::colourInput("deQcColPhistNotUsed", "Not expressed", value = "#FF7B00"),
              cellArgs = list(style = "overflow: visible;")
            )
          ),
          tabPanel(
            title = "Bias",
            numericInput(
              inputId = "deQcBiasPval",
              label = "P-value threshold for 'Significants'",
              value = 0.05, min = 0, max = 1, step = 0.01
            ),
            numericInput(
              inputId = "deQcBiasLFC",
              label = "Log2 fold-change threshold for 'Significants'",
              value = 1, min = 0, max = 10, step = 0.25
            ),
            checkboxInput(
              inputId = "deQcBiasLoess",
              label = "Show LOESS trend",
              value = TRUE
            ),
            splitLayout(
              colourpicker::colourInput("deQcColBiasPresent", "Present", value = "black"),
              colourpicker::colourInput("deQcColBiasSig", "Significants", value = "#A6CEE3"),
              cellArgs = list(style = "overflow: visible;")
            ),
            splitLayout(
              sliderInput("deQcBiasPointSize", "Point size",
                          min = 0.1, max = 3, value = 0.8, step = 0.1,
                          width = "85%", ticks = FALSE),
              sliderInput("deQcBiasAlpha", "Alpha",
                          min = 0.05, max = 1, value = 0.3, step = 0.05,
                          width = "85%", ticks = FALSE),
              cellArgs = list(style = "overflow: visible;")
            )
          ),
          tabPanel(
            title = "Highlight",
            helpText(
              "Highlight specific features on the between-group scatter and
              bias plots. Type or select feature IDs / gene symbols, or paste
              a list."
            ),
            selectizeInput(
              inputId = "deQcHighlightSelect",
              label = "Select features to highlight",
              choices = NULL, selected = NULL, multiple = TRUE,
              options = list(
                placeholder = "Select features",
                plugins = list("remove_button", "drag_drop", "restore_on_backspace", "clear_button")
              )
            ),
            textAreaInput(
              inputId = "deQcHighlightText",
              label = "Or, paste list of features:",
              placeholder = "feature1 feature2 feature3",
              cols = 1
            ),
            splitLayout(
              colourpicker::colourInput("deQcColHighlight", "Highlight colour",
                                        value = "goldenrod3", width = "85%"),
              sliderInput("deQcHighlightSize", "Highlight point size",
                          min = 1, max = 8, value = 3, step = 0.5,
                          width = "85%", ticks = FALSE),
              cellArgs = list(style = "overflow: visible;")
            ),
            checkboxInput("deQcHighlightLabel", "Label highlighted features", value = TRUE),
            sliderInput(
              inputId = "deQcHighlightLabelSize",
              label = "Label size", min = 2, max = 14, value = 4, step = 0.5,
              width = "85%", ticks = FALSE
            )
          ),
          tabPanel(
            title = "Figure",
            splitLayout(
              sliderInput("deQcTextSize", "Font size",
                          min = 4, max = 30, value = 12, step = 0.5,
                          width = "85%", ticks = FALSE),
              sliderInput("deQcFigWidth", "Figure width (px)",
                          min = 300, max = 2000, value = 700, step = 10,
                          width = "85%", ticks = FALSE),
              sliderInput("deQcFigHeight", "Figure height (px)",
                          min = 300, max = 2000, value = 500, step = 10,
                          width = "85%", ticks = FALSE),
              cellArgs = list(style = "overflow: visible;")
            )
          ),
          tabPanel(
            title = "Download",
            selectInput(
              inputId = "deQcDlFormat",
              label = "Format",
              choices = c("PDF", "SVG", "PNG"),
              selected = "PDF"
            ),
            selectInput(
              inputId = "deQcDlDpi",
              label = "PNG DPI",
              choices = c(72, 150, 300, 600),
              selected = 300
            ),
            textInput(
              inputId = "deQcDlName",
              label = "Filename prefix",
              value = "DE_QC"
            ),
            downloadButton("dlPlot_betweenGroupScatter", "Download between-group scatter"),
            br(), br(),
            downloadButton("dlPlot_pvalHistogram", "Download p-value histogram"),
            br(), br(),
            downloadButton("dlPlot_gcBias", "Download GC bias"),
            br(), br(),
            downloadButton("dlPlot_lengthBias", "Download length bias")
          )
        )
      )
    ),
    column(
      width = 8,
      box(
        title = "Results",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          id = "deQcPlotsTabs",
          tabPanel(
            title = "Between-group Scatter",
            textOutput("deQcScatterDesign"),
            helpText("Click and drag to select features; the table below lists
                     the selected points."),
            withSpinner(plotOutput("betweenGroupScatter", inline = TRUE,
                                    brush = brushOpts(id = "deQcScatterBrush")),
                        color = "#6E8D9E"),
            br(),
            DT::dataTableOutput("deQcScatterBrushTable")
          ),
          tabPanel(
            title = "P-value Distribution",
            withSpinner(plotOutput("pvalHistogram", inline = TRUE), color = "#6E8D9E")
          ),
          tabPanel(
            title = "GC Bias",
            helpText("Click and drag to select features; table below lists the selection."),
            withSpinner(plotOutput("gcBiasPlot", inline = TRUE,
                                    brush = brushOpts(id = "deQcGCBrush")),
                        color = "#6E8D9E"),
            br(),
            DT::dataTableOutput("deQcGCBrushTable")
          ),
          tabPanel(
            title = "Length Bias",
            helpText("Click and drag to select features; table below lists the selection."),
            withSpinner(plotOutput("lengthBiasPlot", inline = TRUE,
                                    brush = brushOpts(id = "deQcLengthBrush")),
                        color = "#6E8D9E"),
            br(),
            DT::dataTableOutput("deQcLengthBrushTable")
          )
        )
      )
    )
  )
)
