tabItem(
    tabName = "replicabilityTab",
    fluidRow(
        column(
            width = 4,
            box(
                title = "Settings",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                selectInput(inputId = "replicabilityCounts", label = "Select count method", choices = NULL, selected = NULL),
                conditionalPanel(
                    condition = "input.replicabilityCounts == 'TPM' || input.replicabilityCounts == 'FPKM' || input.replicabilityCounts == 'Normalised' || input.replicabilityCounts == 'Raw'",
                    checkboxInput(inputId = "replicabilityCountsLog", label = "Log2 Counts?", value = FALSE),
                    helpText("Only applies to counts that are not already on in log space.")
                ),
                uiOutput("replicabilityFactorSelectUI"),
                sliderInput(inputId = "replicabilityTopGenes", label = "Number of top variable genes", min = 100, max = 5000, value = 1000, step = 100),
                selectInput(inputId = "replicabilityCorrMethod", label = "Correlation Method", choices = c("pearson", "spearman", "kendall"), selected = "pearson"),
                selectInput(inputId = "replicabilityDistMethod", label = "Distance Method", choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), selected = "euclidean"),
                selectInput(inputId = "replicabilityClustMethod", label = "Clustering Method", choices = c("ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"), selected = "ward.D2"),
                hr(),
                h4("Heatmap colours"),
                splitLayout(
                    colourpicker::colourInput(inputId = "replicabilityColorLow", label = "Low", value = "dodgerblue4", width = "85%"),
                    colourpicker::colourInput(inputId = "replicabilityColorMid", label = "Mid", value = "white", width = "85%"),
                    colourpicker::colourInput(inputId = "replicabilityColorHigh", label = "High", value = "firebrick", width = "85%"),
                    cellArgs = list(style = "overflow: visible;")
                ),
                splitLayout(
                    numericInput(inputId = "replicabilityAtLow", label = "Low break", min = -1, max = 1, value = -1, step = 0.1, width = "85%"),
                    numericInput(inputId = "replicabilityAtMid", label = "Mid break", min = -1, max = 1, value = 0, step = 0.1, width = "85%"),
                    numericInput(inputId = "replicabilityAtHigh", label = "High break", min = -1, max = 1, value = 1, step = 0.1, width = "85%"),
                    cellArgs = list(style = "overflow: visible;")
                ),
                selectInput(inputId = "replicabilityPalette", label = "Or, select palette", choices = c("None", rownames(brewer.pal.info)[which(brewer.pal.info$category %in% c("div", "seq"))]), selected = "None"),
                checkboxInput(inputId = "replicabilityPaletteRev", "Reverse palette", value = FALSE),
                splitLayout(
                    colourpicker::colourInput(inputId = "replicabilityBorderColor1", label = "Anno. border colour", value = "black", width = "85%", palette = "limited"),
                    sliderInput(inputId = "replicabilityBorderSize1", label = "Anno. border size", min = 0, max = 20, value = 5, step = 1, ticks = FALSE, width = "85%")
                ),
                splitLayout(
                    colourpicker::colourInput(inputId = "replicabilityBorderColor2", label = "Main border colour", value = "white", width = "85%", palette = "limited"),
                    sliderInput(inputId = "replicabilityBorderSize2", label = "Main border size", min = 0, max = 20, value = 5, step = 1, ticks = FALSE, width = "85%")
                ),
                checkboxInput(inputId = "replicabilityShowNames", label = "Show sample names", value = TRUE),
                splitLayout(
                    sliderInput(inputId = "replicabilityTextSize", label = "Text Size", min = 4, max = 20, value = 11, step = 1, width = "85%"),
                    sliderInput(inputId = "figWidthReplicability", label = "Figure Width", min = 400, max = 1200, value = 800, step = 10, width = "85%"),
                    sliderInput(inputId = "figHeightReplicability", label = "Figure Height", min = 400, max = 1200, value = 600, step = 10, width = "85%"),
                    cellArgs = list(style = "overflow: visible;")
                )
            )
        ),
        column(
            width = 8,
            tabsetPanel(
                id = "replicabilityTabs",
                tabPanel(
                    "All Genes",
                    br(),
                    box(
                        title = "Sample Correlation Heatmap",
                        width = NULL,
                        solidHeader = TRUE,
                        status = "primary",
                        plotOutput("replicabilityHeatmap", height = "auto")
                    ),
                    box(
                        title = "Sample Clustering Dendrogram",
                        width = NULL,
                        solidHeader = TRUE,
                        status = "primary",
                        plotOutput("replicabilityDendrogram", height = "auto", inline = TRUE)
                    )
                ),
                tabPanel(
                    "Top 100 Abundant",
                    br(),
                    box(
                        title = "Sample Correlation Heatmap (Top 100)",
                        width = NULL,
                        solidHeader = TRUE,
                        status = "primary",
                        plotOutput("replicabilityHeatmapTop100", height = "auto")
                    ),
                    box(
                        title = "Sample Clustering Dendrogram (Top 100)",
                        width = NULL,
                        solidHeader = TRUE,
                        status = "primary",
                        plotOutput("replicabilityDendrogramTop100", height = "auto", inline = TRUE)
                    )
                )
            ),
            box(
                title = "Download",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                downloadButton("replicabilityDownloadHeatmap", "Download Heatmap (PDF)"),
                downloadButton("replicabilityDownloadDendrogram", "Download Dendrogram (PDF)"),
                downloadButton("replicabilityDownloadCorr", "Download Correlation Matrix (Excel)")
            )
        )
    )
)
