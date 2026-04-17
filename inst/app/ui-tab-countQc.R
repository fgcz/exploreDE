tabItem(
  tabName = "countQcTab",
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
          "The Count QC tab provides diagnostic plots on the counts
          themselves (not the DE results). These help spot outlier
          samples, dominant / contaminating genes, and bias by GC
          content or gene length."
        ),
        tags$p(
          tags$b("Total reads / features above threshold:"),
          " per-sample bar charts. The percentage label on the second
          chart shows how many features pass the present-flag threshold."
        ),
        tags$p(
          tags$b("Dominant genes:"),
          " boxplot of per-gene expression fraction, with the top-N
          most over-represented genes labelled per sample. Also
          highlights known gene families (mitochondrial, ribosomal,
          haemoglobin, etc.) for spotting common contamination."
        ),
        tags$p(
          tags$b("Intra-group scatter:"),
          " for a selected group level, each replicate (y) is plotted
          against the group geometric mean (x) on log10 axes. Points
          can be coloured by Standard (Absent/Present/Significants),
          GC bins, or feature-width bins. Points drifting from the 1:1
          line reveal sample-specific or technical bias."
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          id = "countQcSettingsTabs",
          tabPanel(
            title = "Bars",
            selectInput(
              inputId = "countQcBarColourBy",
              label = "Colour bars by factor",
              choices = NULL, selected = NULL
            ),
            helpText("Colours are taken from the Colours tab (GroupColour settings).
                     Drag levels between the Include / Exclude lists to filter and
                     reorder groups; samples follow the group order."),
            uiOutput("countQcBarBucketUI"),
            hr(),
            selectInput(
              inputId = "countQcBarFacetBy",
              label = "Facet bars by factor (optional)",
              choices = "(none)", selected = "(none)"
            ),
            helpText("If set, bar charts split into panels by this factor.
                     Drag the facet levels below to pick which panels appear
                     and their order."),
            uiOutput("countQcBarFacetBucketUI")
          ),
          tabPanel(
            title = "Dominant",
            numericInput(
              inputId = "countQcDomTopN",
              label = "Top-N genes per sample to label",
              value = 3, min = 0, max = 20, step = 1
            ),
            checkboxGroupInput(
              inputId = "countQcDomPresets",
              label = "Gene-family presets",
              choices = c(
                "Mitochondrial"    = "Mitochondrial",
                "Ribosomal (cyto)" = "Ribosomal (cyto)",
                "Ribosomal (mito)" = "Ribosomal (mito)",
                "Haemoglobin"      = "Haemoglobin",
                "Immunoglobulin"   = "Immunoglobulin",
                "TCR"              = "TCR",
                "Heat shock"       = "Heat shock",
                "Histone"          = "Histone"
              ),
              selected = NULL
            ),
            textInput(
              inputId = "countQcDomPattern",
              label = "Custom regex (OR'd with presets)",
              placeholder = "e.g. ^MT- or ^HB[ABZ]"
            ),
            sliderInput(
              inputId = "countQcDomFamilyLabelsPerSample",
              label = "Family labels per sample (top N by fraction)",
              min = 0, max = 20, value = 3, step = 1,
              width = "85%", ticks = FALSE
            ),
            helpText("All family-matching points are drawn as triangles;
                     labels are capped to the top-N per sample to keep the
                     plot readable."),
            radioButtons(
              inputId = "countQcDomColourBy",
              label = "Colour top-N points by",
              choices = c("Rank" = "Rank", "Gene family" = "Family"),
              selected = "Rank", inline = TRUE
            ),
            helpText("Rank: interpolates between the three rank colours so all
                     top-N points get a unique shade. Gene family: each top-N
                     point is coloured by which preset/custom family it belongs
                     to ('Other' = unmatched)."),
            splitLayout(
              colourpicker::colourInput("countQcColDomRank1", "Rank 1 / Low", value = "red"),
              colourpicker::colourInput("countQcColDomRank2", "Rank mid", value = "orange"),
              colourpicker::colourInput("countQcColDomRank3", "Rank N / High", value = "gold"),
              colourpicker::colourInput("countQcColDomFamily", "Family hit / Other", value = "#1965B0"),
              cellArgs = list(style = "overflow: visible;")
            )
          ),
          tabPanel(
            title = "Intra-group",
            uiOutput("countQcIntraFactorUI"),
            uiOutput("countQcIntraGroupUI"),
            selectInput(
              inputId = "countQcIntraCountType",
              label = "Counts to use",
              choices = NULL, selected = NULL
            ),
            radioButtons(
              inputId = "countQcScatterColourBy",
              label = "Colour by:",
              choices = c("Standard", "GC content", "Gene width"),
              selected = "Standard", inline = TRUE
            ),
            conditionalPanel(
              condition = "input.countQcScatterColourBy == 'Standard'",
              numericInput("countQcScatterPval", "P-value threshold",
                           value = 0.05, min = 0, max = 1, step = 0.01),
              numericInput("countQcScatterLFC", "Log2 FC threshold",
                           value = 1, min = 0, max = 10, step = 0.25),
              splitLayout(
                colourpicker::colourInput("countQcColAbsent",  "Absent",       value = "grey70"),
                colourpicker::colourInput("countQcColPresent", "Present",      value = "black"),
                colourpicker::colourInput("countQcColSig",     "Significants", value = "#A6CEE3"),
                cellArgs = list(style = "overflow: visible;")
              ),
              helpText("CountQC data has no DE results, so 'Significants' is unused there.")
            ),
            conditionalPanel(
              condition = "input.countQcScatterColourBy == 'GC content'",
              splitLayout(
                numericInput("countQcScatterGCLow",  "GC low",  value = 0.4, min = 0, max = 1, step = 0.05),
                numericInput("countQcScatterGCHigh", "GC high", value = 0.6, min = 0, max = 1, step = 0.05),
                cellArgs = list(style = "overflow: visible;")
              ),
              splitLayout(
                colourpicker::colourInput("countQcColGCLow",  "GC low bin",  value = "#1F78B4"),
                colourpicker::colourInput("countQcColGCMid",  "GC mid bin",  value = "black"),
                colourpicker::colourInput("countQcColGCHigh", "GC high bin", value = "#FF7F00"),
                cellArgs = list(style = "overflow: visible;")
              )
            ),
            conditionalPanel(
              condition = "input.countQcScatterColourBy == 'Gene width'",
              splitLayout(
                numericInput("countQcScatterWidthLow",  "Width low (nt)",  value = 500,  min = 0, step = 100),
                numericInput("countQcScatterWidthHigh", "Width high (nt)", value = 5000, min = 0, step = 100),
                cellArgs = list(style = "overflow: visible;")
              ),
              splitLayout(
                colourpicker::colourInput("countQcColWidthLow",  "Width low bin",  value = "#B15928"),
                colourpicker::colourInput("countQcColWidthMid",  "Width mid bin",  value = "black"),
                colourpicker::colourInput("countQcColWidthHigh", "Width high bin", value = "#33A02C"),
                cellArgs = list(style = "overflow: visible;")
              )
            ),
            splitLayout(
              sliderInput("countQcScatterPointSize", "Point size",
                          min = 0.1, max = 3, value = 0.6, step = 0.1,
                          width = "85%", ticks = FALSE),
              sliderInput("countQcScatterAlpha", "Alpha",
                          min = 0.05, max = 1, value = 0.4, step = 0.05,
                          width = "85%", ticks = FALSE),
              cellArgs = list(style = "overflow: visible;")
            ),
            numericInput("countQcIntraNcol", "Max columns in grid",
                         value = 4, min = 1, max = 8, step = 1, width = "85%"),
            hr(),
            h5("Guide lines"),
            helpText("Solid line: 1:1 identity (sample = group mean). Dashed lines: N-fold change bands."),
            splitLayout(
              colourpicker::colourInput("countQcColGuide", "Guide colour",
                                        value = "blue", width = "85%"),
              numericInput("countQcGuideFold", "Fold-change threshold",
                           value = 2, min = 1, max = 20, step = 0.5, width = "85%"),
              cellArgs = list(style = "overflow: visible;")
            ),
            splitLayout(
              checkboxInput("countQcShowIdentity", "Show 1:1 line", value = TRUE),
              checkboxInput("countQcShowFoldBands", "Show fold-change bands", value = TRUE)
            )
          ),
          tabPanel(
            title = "Highlight",
            helpText("Highlight specific features across the intra-group scatter."),
            selectizeInput(
              inputId = "countQcHighlightSelect",
              label = "Select features",
              choices = NULL, selected = NULL, multiple = TRUE,
              options = list(
                placeholder = "Select features",
                plugins = list("remove_button", "drag_drop", "restore_on_backspace", "clear_button")
              )
            ),
            textAreaInput(
              inputId = "countQcHighlightText",
              label = "Or, paste list of features:",
              placeholder = "feature1 feature2 feature3",
              cols = 1
            ),
            splitLayout(
              colourpicker::colourInput("countQcColHighlight", "Highlight colour",
                                        value = "goldenrod3", width = "85%"),
              sliderInput("countQcHighlightSize", "Highlight point size",
                          min = 1, max = 8, value = 3, step = 0.5,
                          width = "85%", ticks = FALSE),
              cellArgs = list(style = "overflow: visible;")
            ),
            checkboxInput("countQcHighlightLabel", "Label highlighted features", value = TRUE),
            sliderInput("countQcHighlightLabelSize", "Label size",
                        min = 2, max = 14, value = 3.5, step = 0.5,
                        width = "85%", ticks = FALSE)
          ),
          tabPanel(
            title = "Figure",
            splitLayout(
              sliderInput("countQcTextSize", "Font size",
                          min = 4, max = 30, value = 12, step = 0.5,
                          width = "85%", ticks = FALSE),
              sliderInput("countQcFigWidth", "Figure width (px)",
                          min = 300, max = 2400, value = 800, step = 10,
                          width = "85%", ticks = FALSE),
              sliderInput("countQcFigHeight", "Figure height (px)",
                          min = 300, max = 2400, value = 500, step = 10,
                          width = "85%", ticks = FALSE),
              cellArgs = list(style = "overflow: visible;")
            )
          ),
          tabPanel(
            title = "Download",
            selectInput("countQcDlFormat", "Format",
                        choices = c("PDF", "SVG", "PNG"), selected = "PDF"),
            selectInput("countQcDlDpi", "PNG DPI",
                        choices = c(72, 150, 300, 600), selected = 300),
            textInput("countQcDlName", "Filename prefix", value = "Count_QC"),
            downloadButton("dlPlot_countQcTotalReads", "Download total-reads bar"),
            br(), br(),
            downloadButton("dlPlot_countQcFeaturesAbove", "Download features-above-threshold bar"),
            br(), br(),
            downloadButton("dlPlot_countQcDominant", "Download dominant-genes boxplot"),
            br(), br(),
            downloadButton("dlPlot_countQcIntra", "Download intra-group scatter")
          )
        )
      )
    ),
    column(
      width = 8,
      box(
        title = "Results",
        width = NULL, solidHeader = TRUE, status = "primary",
        tabsetPanel(
          id = "countQcPlotsTabs",
          tabPanel(
            title = "Total Reads",
            withSpinner(plotOutput("countQcTotalReads", inline = TRUE), color = "#6E8D9E")
          ),
          tabPanel(
            title = "Features Above Threshold",
            withSpinner(plotOutput("countQcFeaturesAbove", inline = TRUE), color = "#6E8D9E")
          ),
          tabPanel(
            title = "Dominant Genes",
            withSpinner(plotOutput("countQcDominant", inline = TRUE), color = "#6E8D9E")
          ),
          tabPanel(
            title = "Intra-group Scatter",
            withSpinner(plotOutput("countQcIntraScatter", inline = TRUE), color = "#6E8D9E")
          )
        )
      )
    )
  )
)
