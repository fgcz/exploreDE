tabItem(
  tabName = "heatmapTab",
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
          "These heatmaps show the expression counts of the most variable or
          differentially expressed features across the sample groups selected."
        ),
        tags$p(
          tags$b("For RNA-Seq and proteomics data:"), " Features are selected based on
          differential expression p-values. The plots are split into 'Both' (features
          with significant changes in either direction), 'Up' (positive log fold changes),
          and 'Down' (negative log fold changes). The 'Custom' tab allows you to select
          specific features of interest."
        ),
        tags$p(
          tags$b("For CountQC data:"), " Features are selected based on standard deviation
          (variability) across all samples, showing the most variable features regardless
          of their association with any particular experimental condition."
        ),
        tags$p(
          "Values on the heatmap can be on a continuous or diverging
          scale. If continuous is selected, the values shown are the 
          normalised and log2-transformed counts for each feature. 
          If diverging is selected, then for each feature, the mean 
          value of that feature across all samples is subtracted from the
          log2-transformed count, meaning some values become negative.
          These diverging log2 counts are not the same as the log2 fold
          change values from the differential expression test.
          If a diverging scale is selected, it is possible to change 
          the limits of this scale, to limit the impact of few features
          having very extreme values. Viewing feature counts on this 
          diverging scale helps to exaggerate differences."
        ),
        tags$p(
          "Increasing the number of features displayed adds additional 
          features in order of lowest p-value from the differential
          expression test."
        ),
        tags$p(
          "Individual sample groups can be added/removed from the 
          plot and reordered by dragging and dropping groups in the bucket
          list."
        ),
        tags$p(
          "By default, features and samples are clustered via Ward's hierarchical 
          clustering. Either clustering can be turned off via the tick 
          boxes, in which case, features will be displayed in order of 
          p-value along the y-axis, and samples in order as they are in 
          the bucket list."
        ),
        tags$p(
          "Sample names from the x-axis can be removed using the tick box."
        ),
        tags$p(
          "Changes to group colours are not yet automatically reflected in the
           heatmaps (or boxplots and correlation plots). You will have to change
           a setting on the heatmaps and change it back again for the new group
           colours to be updated for now.")
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          tabPanel(
            title = "Feature settings",
            h4("Feature selection"),
            conditionalPanel(
              condition = "output.test == 'RNASeq' || output.test == 'proteomics'",
              helpText("Plot the feature expression per sample of DE features, based on thresholds. If you don't see a plot, start by turning the p-value to raw."),
              numericInput(inputId = "heatmapGeneNumber", label = "Number of Features on Heatmap", min = 5, max = 2000, value = 50, step = 1),
              numericInput(inputId = "lfcHeatmap", label = "Log2 Fold Change Threshold:", value = 0, min = 0, max = 10, step = 0.25),
              selectInput(inputId = "pTypeHeatmap", choices = c("FDR", "Raw"), label = "P-Value:", selected = "Raw"),
              selectInput(inputId = "pThresholdHeatmap", choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), label = "P-Value Threshold:", selected = 0.05),
              selectInput(inputId = "heatmapFeatureDirection", label = "Feature DE direction", choices = c("Both", "Up-regulated", "Down-regulated"), selected = "Both", multiple = FALSE),
              uiOutput(outputId = "nrPeptidesHeatmapUI")
            ),
            conditionalPanel(
              condition = "output.test == 'CountQC'",
              helpText("Features are selected by standard deviation (variability across samples). More variable features will be shown."),
              numericInput(inputId = "heatmapGeneNumberCountQC", label = "Number of Features on Heatmap", min = 5, max = 5000, value = 2000, step = 10),
              helpText("Plotting > 1000 genes will take a few moments...")
            ),
            
            hr(style = "border-top: 1px solid #000000;"), h4("Count settings"),
            selectInput(inputId = "heatmapCounts", label = "Select count method to plot", choices = "", selected = ""),
            selectInput(inputId = "heatmapBatch", label = "Apply a batch-correction to the counts?", choices = "", multiple = FALSE),
            checkboxInput(inputId = "heatmapLog2", label = "Log2 transform counts?", value = TRUE),
            helpText("Only applies to counts not already in log space."),
            radioButtons(inputId = "heatmapZScore", label = "Centre or use Z-Scores?", choices = c("None", "Centred", "Z-Score"), selected = "Centred"),
            helpText("Calculates per-gene Z-Scores. Works on any count method."),

            hr(style = "border-top: 1px solid #000000;"), h4("Sample aggregation"),
            checkboxInput(inputId = "heatmapShowMeans", label = "Show per-condition aggregation (instead of individual samples)", value = FALSE),
            conditionalPanel(
              condition = "input.heatmapShowMeans == true",
              radioButtons(inputId = "heatmapAggMethod", label = "Aggregation method:", choices = c("Mean", "Median"), selected = "Mean", inline = TRUE),
              helpText("When enabled, samples are grouped by the first selected factor and their aggregated expression is displayed."),
              helpText("Additionally, not all secondary factors will be addable to the heatmap annotation (since there can be multiple per merged group).")
            ),
          ),
          
          tabPanel(
            title = "Heatmap settings",
            helpText("Select the factors to add as annotations. Click and drag to change their order. Sample order is determined by sorting in the factor orders. The first factor is considered the 'main' one and you can change order/remove groups using the bucket below."),
            uiOutput(outputId = "heatmapColumnSelectUI", inline = TRUE),
            uiOutput(outputId = "clusterWarningProteomics", inline = TRUE),
            splitLayout(
              checkboxInput(inputId = "clusterColsHeatmap", label = "Cluster samples?", value = TRUE),
              checkboxInput(inputId = "showClusterColDend", label = "Show sample dendrogram?", value = TRUE)
            ), 
            splitLayout(
              checkboxInput(inputId = "clusterRowsHeatmap", label = "Cluster features?", value = TRUE),
              checkboxInput(inputId = "showClusterRowDend", label = "Show feature dendrogram?", value = TRUE)
            ),
            tags$b("Turn off/on names of rows/columns"),
            splitLayout(
              checkboxInput(inputId = "colnamesHeatmap", label = "Show sample names?", value = TRUE),
              checkboxInput(inputId = "geneNamesHeatmap", label = "Show feature names?", value = TRUE)
            ),
            
            hr(style = "border-top: 1px solid #000000;"), h4("Groups to plot"),
            helpText("Drag and drop groups in order to be plotted"),
            lapply(1:10, function(i) {
              uiOutput(paste0("heatmapBucket", i))
            })
          ),
          
          tabPanel(
            title = "Figure settings",
            h4("Heatmap colours"),
            splitLayout(
              colourpicker::colourInput(inputId = "heatmapColourBlue", label = "Low", value = "dodgerblue4", width = "85%"),
              colourpicker::colourInput(inputId = "heatmapColourWhite", label = "Mid", value = "white", width = "85%"),
              colourpicker::colourInput(inputId = "heatmapColourRed", label = "High", value = "firebrick4", width = "85%"),
              cellArgs = list(style = "overflow: visible;")
            ),
            splitLayout(
              numericInput(inputId = "heatmapAtLow", label = "Low break", min = -10, max = 10, value = -1, step = 0.5, width = "85%"),
              numericInput(inputId = "heatmapAtMid", label = "Mid break", min = -10, max = 10, value = 0, step = 0.5, width = "85%"),
              numericInput(inputId = "heatmapAtHigh", label = "High break", min = -10, max = 10, value = 1, step = 0.5, width = "85%"),
              cellArgs = list(style = "overflow: visible;")
            ),
            selectInput(inputId = "heatmapLimitCPalette", label = "Or, select palette", choices = c("None", rownames(brewer.pal.info)[which(brewer.pal.info$category %in% c("div", "seq"))]), selected = ""),
            checkboxInput(inputId = "heatmapLimitCPaletteRev", "Reverse palette", value = FALSE),
            splitLayout(
              colourpicker::colourInput(inputId = "cellBorderColourHeatmap1", label = "Anno. border colour", value = "black", width = "85%", palette = "limited"),
              sliderInput(inputId = "cellBorderSizeHeatmap1", label = "Anno. border size", min = 0, max = 20, value = 5, step = 1, ticks = F, width = "85%"),
            ),
            splitLayout(
              colourpicker::colourInput(inputId = "cellBorderColourHeatmap2", label = "Cell border colour", value = "black", width = "85%", palette = "limited"),
              sliderInput(inputId = "cellBorderSizeHeatmap2", label = "Cell border size", min = 0, max = 20, value = 2, step = 1, ticks = F, width = "85%"),
            ),
            splitLayout(
              sliderInput(inputId = "textSizeHeatmap", label = "Figure Font Size", min = 4, max = 30, value = 12, step = 0.5, width = "85%", ticks = FALSE),
              sliderInput(inputId = "figWidthHeatmap", label = "Figure Width", min = 100, max = 2000, value = 800, step = 10, width = "85%", ticks = FALSE),
              sliderInput(inputId = "figHeightHeatmap", label = "Figure Height", min = 100, max = 2000, value = 800, step = 10, width = "85%", ticks = FALSE)
            )
          ),
          
          tabPanel(
            title = "Download",
            h4("Download Count Tables"),
            downloadButton(outputId = "dlTable_heatmapDE", label = "Download DE Features Counts"),
            helpText("Downloads counts from the DE Features heatmap (whichever view you have selected: Both, Up, or Down)."),
            checkboxInput(inputId = "includeDEStatsHeatmap", label = "Include DE statistics in count downloads?", value = FALSE),
            helpText("Adds log2FC, p-value, and FDR columns to downloaded count files."),
            br(),
            downloadButton(outputId = "dlTable_heatmapCustom", label = "Download Custom Heatmap Counts"),
            helpText("Downloads counts from the Custom heatmap."),
            br(),
            downloadButton(outputId = "dlTable_heatmapGO", label = "Download GO Heatmap Counts"),
            helpText("Downloads counts from the GO term heatmap."),
            hr(style = "border-top: 1px solid #000000;"),
            h4("Download Heatmap Plots"),
            helpText("Whichever heatmap you are currently viewing will be the one downloaded."),
            downloadButton(outputId = "dlPlot_heatmapDE", label = "Download DE Features Heatmap"),
            br(), br(),
            downloadButton(outputId = "dlPlot_heatmapCustom", label = "Download Custom Heatmap"),
            br(), br(),
            downloadButton(outputId = "dlPlot_heatmapGO", label = "Download GO Heatmap"),
            hr(style = "border-top: 1px solid #000000;"),
            h4("Download Settings"),
            splitLayout(
              selectInput(inputId = "heatmapDownloadFormat", label = "Format", choices = c("PDF", "SVG", "PNG"), selected = "PDF"),
              selectInput(inputId = "heatmapDPI", label = "DPI", choices = c(72, 150, 300, 600, 1000), selected = 600)
            ),
            textInput(inputId = "heatmapFilename", label = "Enter filename", value = "heatmap", placeholder = "up_regulated_heatmap"),
            selectInput(inputId = "heatmapFont", label = "Font", choices = c("serif", "sans", "mono"), selected = "sans")
          )
        )
      )
    ), 
    column(
      width = 8,
      box(
        title = "Heatmaps",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          id = "sig",
          tabPanel(
            title = uiOutput("heatmapTabTitle", inline = TRUE),
            textOutput("heatmapDesignBoth Directions"),
            br(),
            withSpinner(
              plotOutput("heatmapDEFeatures", inline = TRUE),
              color = "#6E8D9E"
            ),
            textOutput("hGNDE")
          ),
          
          tabPanel(
            title = "Custom", 
            column(
              width = 3,
              helpText("Plot the expression of any feature detected, irrespective of the DE test. Requires at least three features."),
              uiOutput("fb_search_mode_heatmap"),
              selectizeInput(
                inputId = "fb_select_heatmap", label = "Features to heatmap:", multiple = TRUE, choices = NULL, selected = NULL,
                options = list(placeholder = 'Select features', plugins = list('remove_button', 'drag_drop', 'restore_on_backspace', 'clear_button'))
              ),
              textAreaInput(inputId = "fb_text_heatmap", label = "Or, paste list of features:", placeholder = "feature1 feature2 feature3 feature4 ", cols = 1),
              textAreaInput(inputId = "heatmapCustomTitle", label = "Enter a title for the heatmap", value = "Custom Heatmap"),
              hr(style = "border-top: 1px solid #000000;"), h4("Feature Bucket"),
              uiOutput("fb_bucket_heatmap"),
              actionButton(inputId = "fb_reset_heatmap", label = "Empty the bucket?", icon = icon("bucket")),
            ),
            column(
              width = 9,
              withSpinner(
                plotOutput(paste0("heatmapCustom"), inline = TRUE),
                color = "#6E8D9E"
              )
            )
          ),
          tabPanel(
            title = "GO", 
            column(
              width = 3,
              helpText("You can also plot the expression of features annotated to a given GO term (again, irrespective of p-value from the DE test)."),
              selectizeInput(
                inputId = "heatmapGoInput",
                label = "Or, select the GO Term you wish to plot",
                choices = NULL,
                selected = NULL
              ),
              uiOutput("heatmapGOUI2")
            ),
            column(
              width = 9,
              withSpinner(
                plotOutput(paste0("heatmapGO"), inline = TRUE),
                color = "#6E8D9E"
              )
            )
          )
        )
      )
    )
  )
)