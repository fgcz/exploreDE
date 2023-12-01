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
          "These heatmaps show the expression counts of the most 
          differentially expressed features (i.e., lowest p-value) 
          across the sample groups selected. The plots are split
          into four main sections. The first, 'Both', shows the most 
          differentially expressed features in either direction, i.e.,
          both positive and negative log fold changes. The second
          and third show features with either positive ('Up') or negative
          ('Down') log fold changes only. The last, 'Custom', will show
          features selected from the custom heatmap input. When selecting 
          features for the custom heatmap, the features in the list are in 
          order of p-value, and can be searched for by typing in 
          the feature name. They will automatically be added to the plot. 
          feature names can also be pasted in the text box."
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
            title = "Main settings",
          uiOutput(outputId = "heatmapProteomicsColumnSelectUI", inline = TRUE),
          helpText("Plot the feature expression per sample of DE features, based on thresholds."),
          numericInput(inputId = "heatmapGeneNumber", label = "Number of Features on Heatmap", min = 5, max = 5000, value = 50, step = 1),
          numericInput(inputId = "lfcHeatmap", label = "Log2 Fold Change Threshold:", value = 0.5, min = 0, max = 10, step = 0.25),
          selectInput(inputId = "pTypeHeatmap", choices = c("FDR", "Raw"), label = "P-Value:", selected = "FDR"),
          selectInput(inputId = "pThresholdHeatmap", choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), label = "P-Value Threshold:", selected = 0.05),
          hr(style = "border-top: 1px solid #000000;"), h4("Count settings"),
          selectInput(inputId = "heatmapCounts", label = "Select count method to plot", choices = "", selected = ""),
          selectInput(inputId = "heatmapBatch", label = "Apply a batch-correction to the counts?", choices = "", multiple = TRUE),
          hr(style = "border-top: 1px solid #000000;"), h4("Scale settings"),
          radioButtons(inputId = "heatmapScale", label = "Scale", choices = c("Continuous", "Diverging"), selected = "Diverging"),
          sliderInput(inputId = "heatmapLimitD", label = "Scale limit for diverging heatmap (+/-)", value = 4, min = 1, max = 20, step = 0.5),
          selectInput(inputId = "heatmapLimitCPalette", label = "Continuous heatmap palette", choices = rownames(brewer.pal.info)[which(brewer.pal.info$category %in% c("div", "seq"))], selected = "YlGnBu"),
          checkboxInput(inputId = "heatmapLimitCPaletteRev", label = "Click to reverse palette", value = FALSE),
          numericInput(inputId = "heatmapLimitCHigh", label = "Scale limit for continuous heatmap", value = 1e3, min = 1, max = 1e5, step = 1),
          # h5("Select high, mid, and low scale limits for continuous heatmap"),
          # splitLayout(
          #   sliderInput(inputId = "heatmapLimitCMid", label = "Mid", value = 1e3, min = 1, max = 1e5, step = 1),
          #   sliderInput(inputId = "heatmapLimitCLow", label = "Low", value = 1e2, min = 1, max = 1e5, step = 1)
          # ),
          hr(style = "border-top: 1px solid #000000;"), h4("Groups to plot"),
          uiOutput("heatmapBucket")
          ),
          tabPanel(
            title = "Custom heatmap",
            helpText("This allows you to plot the expression of any feature detected, irrespective of the DE test."),
            selectInput(inputId = "heatmapGenes", label = "Select features for custom heatmap:", multiple = TRUE, choices = "", selected = ""),
            textAreaInput(inputId = "heatmapGenesText", label = "Or, paste list of features:", placeholder = "feature1 feature2 feature3 feature4 ", cols = 1),
            h5(
              "Features will be added to this bucket as you select them from the DE table 
              tab and from the inputs in this tab. You can use this bucket to quickly 
              re-order and exclude these features as you need to by dragging and dropping 
              them in order or into the exclude bucket."),
            uiOutput("geneBucket3")
          ),
          tabPanel(
            title = "GO heatmap",
            helpText("You can also plot the expression of features annotated to a given GO term (again, irrespective of p-value from the DE test)."),
            uiOutput("heatmapGOUI1"),
            uiOutput("heatmapGOUI2")
          ),
          tabPanel(
            title = "Figure settings",
            tags$b("Turn off/on clustering of rows/columns"),
            uiOutput(outputId = "clusterWarningProteomics", inline = TRUE),
            checkboxInput(inputId = "clusterColsHeatmap", label = "Cluster samples?", value = TRUE),
            checkboxInput(inputId = "showClusterColDend", label = "Show sample dendrogram?", value = TRUE),
            checkboxInput(inputId = "clusterRowsHeatmap", label = "Cluster features?", value = TRUE),
            checkboxInput(inputId = "showClusterRowDend", label = "Show feature dendrogram?", value = TRUE),
            tags$b("Turn off/on names of rows/columns"),
            checkboxInput(inputId = "colnamesHeatmap", label = "Show sample names?", value = TRUE),
            checkboxInput(inputId = "geneNamesHeatmap", label = "Show feature names?", value = TRUE),
            uiOutput("heatmapFactors"),
            checkboxGroupInput(inputId = "heatmapFactor2", label = "Select which additional factors should be displayed", choices = "", selected = ""),
            colourpicker::colourInput(inputId = "heatmapColourRed", label = "Heatmap Red (+)", value = "#801717"),
            colourpicker::colourInput(inputId = "heatmapColourWhite", label = "Heatmap White (0)", value = "#FFFFFF"),
            colourpicker::colourInput(inputId = "heatmapColourBlue", label = "Heatmap Blue (-)", value = "#113D69"),
            numericInput(inputId = "textSizeHeatmap", label = "Figure Font Size", min = 4, max = 30, value = 12, step = 0.5),
            numericInput(inputId = "figWidthHeatmap", label = "Figure Width", min = 100, max = 2000, value = 800, step = 10),
            numericInput(inputId = "figHeightHeatmap", label = "Figure Height", min = 100, max = 2000, value = 800, step = 10)
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
        textOutput("heatmapDesign"), br(),
        # Generate tabs, output, and download buttons for
        # 3 x heatmaps: Both, up, and down regulated features:
        do.call(tabsetPanel, c(
          id = "sig",
          lapply(c("Both Directions", "Up-Regulated", "Down-Regulated", "Custom", "GO"), function(sig) {
              tabPanel(
                paste0(sig, " Features"),
                downloadButton(outputId = paste0("dlHeatmap", sig, "Button"), label = paste0("Download ", sig, " Heatmap (PDF)")),
                downloadButton(outputId = paste0("dlHeatmap", sig, "ButtonPNG"), label = paste0("Download ", sig, " Heatmap (PNG)")),
                downloadButton(outputId = paste0("dlHeatmap", sig, "DFButton"), label = paste0("Download ", sig, " Heatmap Results (Excel)")),
                br(),
                plotOutput(paste0("heatmap", sig), inline = TRUE),
                br(),
                textOutput(paste0(sig, "HeatmapZeroMeanMessage"))
              )
            }
          )
        ))
      )
    )
  )
)