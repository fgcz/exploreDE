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
                h4("Heatmap Appearance"),
                fluidRow(
                    column(4, colourpicker::colourInput("replicabilityColorLow", "Low", "blue")),
                    column(4, colourpicker::colourInput("replicabilityColorMid", "Mid", "white")),
                    column(4, colourpicker::colourInput("replicabilityColorHigh", "High", "red"))
                ),
                sliderInput(inputId = "replicabilityBorderSize", label = "Cell Border Size", min = 0, max = 5, value = 0.5, step = 0.1),
                colourpicker::colourInput("replicabilityBorderColor", "Cell Border Color", "white"),
                checkboxInput(inputId = "replicabilityShowNames", label = "Show sample names", value = TRUE),
                sliderInput(inputId = "replicabilityTextSize", label = "Text Size", min = 4, max = 20, value = 10, step = 1),
                sliderInput(inputId = "figWidthReplicability", label = "Figure Width", min = 400, max = 1200, value = 800, step = 10),
                sliderInput(inputId = "figHeightReplicability", label = "Figure Height", min = 400, max = 1200, value = 600, step = 10)
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
