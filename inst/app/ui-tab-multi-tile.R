tabItem(
    tabName = "multiTilePlotTab",
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
                h4("Intro"),
                tags$p("This tab is split into three main sections: The top differentially expressed genes; a plot of selected genes, and; a plot by GO pathway genes."),
                tags$p(
                    "Each column is one of the contrasts you uploaded, and each row is a gene.
          Each cell is coloured by the log2 fold change of that gene for that contrast.
          If a gene is DE in that contrast, it will have stars reflecting the p-value.
          If a gene is not DE in that contrast, the cell will have no stars.
          The colour of the rowname denotes which contrast that gene is from. If a gene name is black, it came from multiple contrasts (e.g., it was in the top 10 of more than one contrast).
          "
                ),
                h4("Top DEGs tab"),
                tags$p(
                    "The top DEG plot has three main settings to determine the top n DEGs. The first is to sort each contrast by p-value, and take the top n.
          The second is to sort each contrast by absolute log2FC (so both positive and negative in a single, decending order) and then take the top n.
          The third option is to sort each contrast by both positive and negative LFC and take the top n/2 from each. Here, if you select an odd number,
          it will round down (e.g., if you select n = 25, it will take top 12 + and top 12 -).
          In all instances, it will first filter each contrast by the p-value and LFC thresholds in the settings box."
                ),
                h4("Custom Genes tab"),
                tags$p(
                    "From the input boxes in the settings pane, you can select any gene to create a custom tile plot.
          The plot needs at least 3 genes before it will be shown. You can
          start typing in gene names in the top drop down input and click the ones you want.
          Or, in the input text box underneath, you can paste in gene names from other locations. Gene names can be either
          space-separated, new line-separated, or '/' separated (e.g., if you are pasting them in from
          the pathway results tab)."
                ),
                h4("GO Pathway tab"),
                tags$p(
                    "Here, you can select a GO pathway from your results to plot the all genes for that pathway.
          Some pathways have thousands of genes, which will make your heatmap a bit hard to interpret,
          so we can filter the pathway gene lists in two ways. These two filters can be run
          either separately or additionally. The first is to select genes that are DE in at least
          one of your contrasts, by ticking the 'Yes, limit my pathway' checkbox. The second
          is to take the top n genes for that pathway from each contrast, in order of p-value. You can enter any number above 1.
          If you both check the box and enter a number, the filtering will be done in that order:
          First only genes that are DE in at least one contrast, and then take the top n."
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
                        h4("Top DEGs Options"),
                        selectInput(inputId = "tileSortBy", label = "Sort DEG table by:", choices = c("P-value", "Absolute Log2FC", "Top n +/- Log2FC"), selected = "P-value"),
                        numericInput(inputId = "tileDEGNumber", label = "Number Top DEGs from each contrast", min = 5, max = 5000, value = 10, step = 1),
                        sliderInput(inputId = "lfcTile", label = "Log2 Fold Change Threshold:", value = 0.5, min = 0, max = 10, step = 0.25, ticks = F, width = "85%"),
                        selectInput(inputId = "pTypeTile", choices = c("FDR", "Raw"), label = "P-Value:", selected = "Raw"),
                        selectInput(inputId = "pThresholdTile", choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), label = "P-Value Threshold:", selected = 0.05)
                    ),
                    tabPanel(
                        title = "Figure settings",
                        h4("Tile Plot Labels/Clustering"),
                        tags$b("Turn off/on clustering of rows/columns"),
                        checkboxInput(inputId = "clusterRowsTile", label = "Cluster genes?", value = TRUE),
                        checkboxInput(inputId = "clusterColsTile", label = "Cluster contrasts?", value = TRUE),
                        tags$b("Turn off/on names of rows/columns"),
                        checkboxInput(inputId = "colnamesTile", label = "Show contrast names?", value = FALSE),
                        # checkboxInput(inputId = "geneNamesTile", label = "Show gene names?", value = TRUE),

                        hr(style = "border-top: 1px solid #000000;"), h4("Tile Plot Colours"),
                        splitLayout(
                            colourpicker::colourInput(inputId = "tileColourBlue", label = "Low", value = "dodgerblue4", width = "85%"),
                            colourpicker::colourInput(inputId = "tileColourWhite", label = "Mid", value = "white", width = "85%"),
                            colourpicker::colourInput(inputId = "tileColourRed", label = "High", value = "firebrick4", width = "85%"),
                            cellArgs = list(style = "overflow: visible;")
                        ),
                        splitLayout(
                            numericInput(inputId = "tileAtLow", label = "Low break", min = -10, max = 10, value = -4, step = 0.5, width = "85%"),
                            numericInput(inputId = "tileAtMid", label = "Mid break", min = -10, max = 10, value = 0, step = 0.5, width = "85%"),
                            numericInput(inputId = "tileAtHigh", label = "High break", min = -10, max = 10, value = 4, step = 0.5, width = "85%"),
                            cellArgs = list(style = "overflow: visible;")
                        ),
                        selectInput(inputId = "tileLimitCPalette", label = "Or, select palette", choices = c("None", rownames(brewer.pal.info)[which(brewer.pal.info$category %in% c("div", "seq"))]), selected = ""),
                        checkboxInput(inputId = "tileLimitCPaletteRev", "Reverse palette", value = FALSE),
                        splitLayout(
                            colourpicker::colourInput(inputId = "cellBorderColourTile1", label = "Anno. border colour", value = "black", width = "85%", palette = "limited"),
                            sliderInput(inputId = "cellBorderSizeTile1", label = "Anno. border size", min = 0, max = 20, value = 5, step = 1, ticks = F, width = "85%"),
                        ),
                        splitLayout(
                            colourpicker::colourInput(inputId = "cellBorderColourTile2", label = "Cell border colour", value = "black", width = "85%", palette = "limited"),
                            sliderInput(inputId = "cellBorderSizeTile2", label = "Cell border size", min = 0, max = 20, value = 2, step = 1, ticks = F, width = "85%"),
                        ),
                        hr(style = "border-top: 1px solid #000000;"), h4("Tile Plot Size"),
                        splitLayout(
                            sliderInput(inputId = "textSizeTile", label = "Figure Font Size", min = 4, max = 30, value = 12, step = 0.5, width = "85%", ticks = FALSE),
                            sliderInput(inputId = "figWidthTile", label = "Figure Width", min = 100, max = 2000, value = 800, step = 10, width = "85%", ticks = FALSE),
                            sliderInput(inputId = "figHeightTile", label = "Figure Height", min = 100, max = 2000, value = 800, step = 10, width = "85%", ticks = FALSE)
                        ),
                        hr(style = "border-top: 1px solid #000000;"), h4("Download Settings"),
                        splitLayout(
                             selectInput(inputId = "downloadFormatMultiTile", label = "Format", choices = c("PDF", "SVG", "PNG"), selected = "PDF", width = "85%"),
                             selectInput(inputId = "dpiMultiTile", label = "DPI", choices = c(72, 150, 300, 600, 1000), selected = 600, width = "85%")
                        )
                    )
                )
            )
        ),
        column(
            width = 9,
            box(
                title = "Tile plots",
                width = NULL,
                solidHeader = TRUE,
                status = "success",
                tabsetPanel(
                    id = "multiTileTabs",
                    tabPanel(
                        title = "Top DEGs",
                        h3("Tile Plot"),
                        textInput(inputId = "filename_tileTop", label = NULL, value = "Top_DEGs_Tile", placeholder = "Filename"),
                        downloadButton(outputId = "dlPlot_tileTop", label = "Download This Plot"),
                        downloadButton(outputId = "dlTable_tileTopRes", label = "Download Results Matrix (Excel)"),
                        br(), br(),
                        plotOutput(outputId = "TopNTileMap", inline = TRUE),
                        br(), br(),
                        h3("Table"),
                        downloadButton(outputId = "dlTable_tileTop", label = "Download This Table (Excel)"),
                        br(), br(),
                        DT::dataTableOutput(outputId = "TopNTileTopTable")
                    ),
                    tabPanel(
                        title = "Custom Genes",
                        column(
                            width = 3,
                            helpText("Plot the DE results of any feature detected, irrespective of ranking."),
                            selectizeInput(
                                inputId = "tileGenes", label = "Features for custom tile:", multiple = TRUE, choices = NULL, selected = NULL,
                                options = list(placeholder = "Select features", plugins = list("remove_button", "drag_drop", "restore_on_backspace", "clear_button"))
                            ),
                            textAreaInput(inputId = "tileGenesText", label = "Or, paste list of features:", placeholder = "feature1 feature2 feature3 feature4 ", cols = 1),
                            textAreaInput(inputId = "tileCustomTitle", label = "Enter a title for the heatmap", value = "Custom Heatmap"),
                            hr(style = "border-top: 1px solid #000000;")
                        ),
                        column(
                            width = 9,
                            h3("Tile Plot"),
                            textInput(inputId = "filename_tileCustom", label = NULL, value = "Custom_Tile", placeholder = "Filename"),
                            downloadButton(outputId = "dlPlot_tileCustom", label = "Download This Plot"),
                            br(), br(),
                            plotOutput(outputId = "CustomTileMap", inline = TRUE),
                            br(), br(),
                            h3("Table"),
                            downloadButton(outputId = "dlTable_tileCustom", label = "Download This Table (Excel)"),
                            br(), br(),
                            DT::dataTableOutput(outputId = "TopNTileCustomTable")
                        )
                    ),
                    tabPanel(
                        title = "GO Pathway Genes",
                        column(
                            width = 3,
                            helpText("You can also plot the DE results of features annotated to a given GO term."),
                            selectInput(
                                inputId = "tileGenesPathway",
                                label = "Select Pathway for custom tile plot:",
                                multiple = FALSE,
                                choices = "",
                                selected = ""
                            ),
                            h5("Selected GO Pathway has too many genes? Want to only show genes that are DE in at least one contrast?"),
                            checkboxInput(inputId = "GOTooMany", label = "Yes, limit my pathway", value = FALSE),
                            h5("Still too many genes? Want to limit to the first n by p-value? If so, enter a number here:"),
                            numericInput(inputId = "GOTooManyN", label = "Top genes by p-value", value = "", min = 5, max = 1000)
                        ),
                        column(
                            width = 9,
                            h3("Tile Plot"),
                            textInput(inputId = "filename_tileGO", label = NULL, value = "GO_Tile", placeholder = "Filename"),
                            downloadButton(outputId = "dlPlot_tileGO", label = "Download This Plot"),
                            br(), br(),
                            plotOutput(outputId = "GOTileMap", inline = TRUE),
                            br(), br(),
                            h3("Table"),
                            downloadButton(outputId = "dlTable_tileGO", label = "Download This Table (Excel)"),
                            br(), br(),
                            DT::dataTableOutput(outputId = "TopNTileGOTable")
                        )
                    )
                )
            )
        )
    )
)
