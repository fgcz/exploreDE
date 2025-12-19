tabItem(
    tabName = "multiHeatmapTab",
    fluidRow(
        column(
            width = 3,
            box(
                title = "Info",
                width = NULL,
                solidHeader = TRUE,
                status = "success",
                collapsible = TRUE,
                collapsed = TRUE,
                tags$p(
                    "This heatmap shows pathway enrichment results across multiple contrasts. Each row represents a pathway, and each column represents a contrast. The heatmap can be coloured by p-value, FDR, or NES (for GSEA only)."
                ),
                h4("Heatmap Interpretation"),
                tags$ul(
                    tags$li("Rows: Pathways from enrichment analysis"),
                    tags$li("Columns: Contrasts loaded in the app"),
                    tags$li("ORA Mode: Colour represents -log10(p-value) or -log10(FDR)"),
                    tags$li("GSEA Mode: Colour can represent -log10(p-value), -log10(FDR), or NES")
                )
            ),
            box(
                title = "Settings",
                width = NULL,
                solidHeader = TRUE,
                status = "success",
                tabsetPanel(
                    tabPanel(
                        title = "Input Options",
                        selectInput(
                            inputId = "dataToHeatmap",
                            label = "Select between ORA or GSEA Results",
                            choices = c("ORA", "GSEA"),
                            selected = "ORA"
                        ),
                        selectInput(
                            inputId = "subdataToHeatmap",
                            label = "Select the direction/ontology to plot",
                            choices = c(),
                            selected = NULL
                        ),
                        hr(style = "border-top: 1px solid #000000;"),
                        h4("Pathway Selection"),
                        sliderInput(
                            inputId = "heatmapPathwayNumber",
                            label = "Number of top pathways from each contrast",
                            min = 5,
                            max = 100,
                            value = 10,
                            step = 5,
                            ticks = FALSE
                        ),
                        selectInput(
                            inputId = "heatmapColourBy",
                            label = "Colour heatmap by:",
                            choices = c("p-value" = "pvalue", "FDR" = "p.adjust"),
                            selected = "p.adjust"
                        ),
                        uiOutput("heatmapNESOption"),
                        hr(style = "border-top: 1px solid #000000;"),
                        h4("Display Options"),
                        checkboxInput(
                            inputId = "colnamesHeatmap",
                            label = "Show contrast names?",
                            value = TRUE
                        ),
                        checkboxInput(
                            inputId = "rownamesHeatmap",
                            label = "Show pathway names?",
                            value = TRUE
                        ),
                        checkboxInput(
                            inputId = "clusterRowsHeatmap",
                            label = "Cluster rows (pathways)?",
                            value = TRUE
                        ),
                        checkboxInput(
                            inputId = "clusterColsHeatmap",
                            label = "Cluster columns (contrasts)?",
                            value = FALSE
                        )
                    ),
                    tabPanel(
                        title = "Figure Settings",
                        h4("Colour Scheme"),
                        colourpicker::colourInput(
                            inputId = "heatmapColourLow",
                            label = "Low value colour",
                            value = "#FFFFFF"
                        ),
                        colourpicker::colourInput(
                            inputId = "heatmapColourHigh",
                            label = "High value colour",
                            value = "#DC050C"
                        ),
                        hr(style = "border-top: 1px solid #000000;"),
                        h4("Size & Fonts"),
                        sliderInput(
                            inputId = "textSizeHeatmap",
                            label = "Base font size",
                            min = 6,
                            max = 24,
                            value = 12,
                            step = 1,
                            ticks = FALSE
                        ),
                        sliderInput(
                            inputId = "figWidthHeatmap",
                            label = "Figure width (px)",
                            min = 800,
                            max = 4000,
                            value = 1000,
                            step = 100,
                            ticks = FALSE
                        ),
                        sliderInput(
                            inputId = "figHeightHeatmap",
                            label = "Figure height (px)",
                            min = 1000,
                            max = 4000,
                            value = 800,
                            step = 100,
                            ticks = FALSE
                        )
                    )
                )
            )
        ),
        column(
            width = 9,
            box(
                title = "Pathway Enrichment Heatmap",
                width = NULL,
                solidHeader = TRUE,
                status = "success",
                h4("Download Settings"),
                splitLayout(
                  selectInput(inputId = "downloadFormatMultiHeatmap", label = "Format", choices = c("PDF", "SVG", "PNG"), selected = "PDF"),
                  selectInput(inputId = "dpiMultiHeatmap", label = "DPI", choices = c(72, 150, 300, 600, 1000), selected = 600)
                ),
                textInput(inputId = "filename_multiHeatmap", label = "Filename", value = "MultiHeatmap"),
                downloadButton(
                    outputId = "dlPlot_multiHeatmap",
                    label = "Download Heatmap"
                ),
                downloadButton(
                    outputId = "dlTable_multiHeatmap",
                    label = "Download Data (Excel)"
                ),
                br(), br(),
                withSpinner(plotOutput(outputId = "pathwayHeatmapPlot", inline = TRUE), color = "#00a65a")
            ),
            box(
                title = "Pathway Data Table",
                width = NULL,
                solidHeader = TRUE,
                status = "success",
                DT::dataTableOutput(outputId = "pathwayHeatmapTable")
            )
        )
    )
)
