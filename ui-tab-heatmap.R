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
          uiOutput(outputId = "nrPeptidesHeatmapUI"),
          hr(style = "border-top: 1px solid #000000;"), h4("Count settings"),
          selectInput(inputId = "heatmapCounts", label = "Select count method to plot", choices = "", selected = ""),
          selectInput(inputId = "heatmapBatch", label = "Apply a batch-correction to the counts?", choices = "", multiple = TRUE),
          hr(style = "border-top: 1px solid #000000;"), h4("Scale settings"),
          radioButtons(inputId = "heatmapScale", label = "Scale", choices = c("Continuous", "Diverging"), selected = "Diverging"),
          sliderInput(inputId = "heatmapLimitD", label = "Scale limit for diverging heatmap (+/-)", value = 4, min = 1, max = 20, step = 0.5),
          selectInput(inputId = "heatmapLimitCPalette", label = "Continuous heatmap palette", choices = rownames(brewer.pal.info)[which(brewer.pal.info$category %in% c("div", "seq"))], selected = "YlGnBu"),
          checkboxInput(inputId = "heatmapLimitCPaletteRev", label = "Click to reverse palette", value = FALSE),
          numericInput(inputId = "heatmapLimitCHigh", label = "Scale limit for continuous heatmap", value = 1e3, min = 1, max = 1e5, step = 1),
          hr(style = "border-top: 1px solid #000000;"), h4("Groups to plot"),
          uiOutput("heatmapBucket")
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
            # uiOutput("heatmapFactors"),
            checkboxGroupInput(inputId = "heatmapFactors", label = "Select which factors should be displayed", choices = "", selected = ""),
            colourpicker::colourInput(inputId = "heatmapColourRed", label = "Heatmap Red (+)", value = "darkred"),
            colourpicker::colourInput(inputId = "heatmapColourWhite", label = "Heatmap White (0)", value = "#FFFFFF"),
            colourpicker::colourInput(inputId = "heatmapColourBlue", label = "Heatmap Blue (-)", value = "deepskyblue4"),
            sliderInput(inputId = "textSizeHeatmap", label = "Figure Font Size", min = 4, max = 30, value = 12, step = 0.5, width = "85%"),
            sliderInput(inputId = "figWidthHeatmap", label = "Figure Width", min = 100, max = 2000, value = 800, step = 10, width = "85%"),
            sliderInput(inputId = "figHeightHeatmap", label = "Figure Height", min = 100, max = 2000, value = 800, step = 10, width = "85%")
          ),
          tabPanel(
            title = "Download Settings",
            helpText("Which ever heatmap you are currently viewing (e.g., both, up, or down) will be the one downloaded."),
            selectInput(inputId = "heatmapDownloadFormat", label = "Select format", choices = c("PDF", "SVG", "PNG"), selected = "PDF"),
            # numericInput(inputId = "heatmapDownloadWidth", label = "Figure Width", min = 100, max = 2000, value = 800, step = 10),
            # numericInput(inputId = "heatmapDownloadHeight", label = "Figure Height", min = 100, max = 2000, value = 800, step = 10),
            selectInput(inputId = "heatmapDPI", label = "PNG DPI", choices = c(72, 150, 300, 600, 1000), selected = 600),
            textInput(inputId = "heatmapFilename", label = "Enter filename", value = "heatmap", placeholder = "up_regulated_heatmap"),
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
            title = "Both Directions",
            textOutput("heatmapDesignBoth Directions"),
            downloadButton(outputId = paste0("dlHeatmapButtonBoth Directions"), label = "Download Heatmap"),
            downloadButton(outputId = paste0("dlHeatmapDFButtonBoth Directions"), label = "Download Counts (Excel)"), 
            withSpinner(
              plotOutput(paste0("heatmapBoth Directions"), inline = TRUE)
            )
          ),
          tabPanel(
            title = "Up-Regulated",
            textOutput("heatmapDesignUp-Regulated"), br(),
            downloadButton(outputId = paste0("dlHeatmapButtonUp-Regulated"), label = "Download Heatmap"),
            downloadButton(outputId = paste0("dlHeatmapDFButtonUp-Regulated"), label = "Download Counts (Excel)"), br(),
            withSpinner(
              plotOutput(paste0("heatmapUp-Regulated"), inline = TRUE)
            )
          ),
          tabPanel(
            title = "Down-Regulated",
            textOutput("heatmapDesignDown-Regulated"), br(),
            downloadButton(outputId = paste0("dlHeatmapButtonDown-Regulated"), label = "Download Heatmap"),
            downloadButton(outputId = paste0("dlHeatmapDFButtonDown-Regulated"), label = "Download Counts (Excel)"), br(),
            withSpinner(
              plotOutput(paste0("heatmapDown-Regulated"), inline = TRUE)
            )
          ),
          tabPanel(
            title = "Custom", 
            column(
              width = 3,
              helpText("Plot the expression of any feature detected, irrespective of the DE test. Requires at least three features."),
              selectInput(inputId = "heatmapGenes", label = "Select features for custom heatmap:", multiple = TRUE, choices = "", selected = ""),
              textAreaInput(inputId = "heatmapGenesText", label = "Or, paste list of features:", placeholder = "feature1 feature2 feature3 feature4 ", cols = 1),
              textAreaInput(inputId = "heatmapCustomTitle", label = "Enter a title for the heatmap", value = "Custom Heatmap"),
              hr(style = "border-top: 1px solid #000000;"), h4("Feature Bucket"),
              uiOutput("geneBucket3"),
              actionButton(inputId = "resetGeneBucketHeatmap", label = "Empty the bucket?", icon = icon("bucket")),
            ),
            column(
              width = 9,
              downloadButton(outputId = paste0("dlHeatmapButtonCustom"), label = "Download Heatmap"),
              downloadButton(outputId = paste0("dlHeatmapDFButtonCustom"), label = "Download Counts (Excel)"), br(),
              withSpinner(
                plotOutput(paste0("heatmapCustom"), inline = TRUE)
              )
            )
          ),
          tabPanel(
            title = "GO", 
            column(
              width = 3,
              helpText("You can also plot the expression of features annotated to a given GO term (again, irrespective of p-value from the DE test)."),
            uiOutput("heatmapGOUI1"),
            uiOutput("heatmapGOUI2")
            ),
            column(
              width = 9,
              downloadButton(outputId = paste0("dlHeatmapButtonGO"), label = "Download Heatmap"),
              downloadButton(outputId = paste0("dlHeatmapDFButtonGO"), label = "Download Counts (Excel)"), br(),
              withSpinner(
                plotOutput(paste0("heatmapGO"), inline = TRUE)
              )
            )
          )
        )
      )
    )
  )
)