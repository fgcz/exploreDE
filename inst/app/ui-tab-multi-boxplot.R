tabItem(
    tabName = "multiBoxplotTab",
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
                    "This tab allows you to visualize expression patterns across multiple contrasts on the same boxplot.
          This only works when all loaded contrasts come from the same experiment (same samples, same count matrix).
          The app will automatically detect if your contrasts are compatible."
                ),
                h4("Requirements"),
                tags$ul(
                    tags$li("All contrasts must share the same samples"),
                    tags$li("All contrasts must have the same underlying count data"),
                    tags$li("At least 2 contrasts must be loaded")
                ),
                h4("P-value Display"),
                tags$p(
                    "Each contrast's statistical comparison will be displayed with stacked brackets above the boxplot.
          Brackets are coloured according to the contrast colours you selected in the Import tab.
          Stars indicate significance: **** (p≤0.0001), *** (p≤0.001), ** (p≤0.01), * (p≤0.05), ns (not significant)."
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
                        h4("Gene Selection"),
                        selectizeInput(
                            inputId = "multiBoxplotGenes",
                            label = "Features for boxplots:",
                            multiple = TRUE,
                            choices = NULL,
                            selected = NULL,
                            options = list(placeholder = "Select features", plugins = list("remove_button", "drag_drop", "restore_on_backspace", "clear_button"))
                        ),
                        textAreaInput(
                            inputId = "multiBoxplotGenesText",
                            label = "Or, paste list of features:",
                            placeholder = "feature1 feature2 feature3",
                            rows = 3
                        ),
                        uiOutput("multiGeneBucket"),
                        hr(style = "border-top: 1px solid #000000;"),
                        h4("Count Data Options"),
                        selectInput(
                            inputId = "multiBoxplotCounts",
                            label = "Which count assay to plot:",
                            choices = c(
                                "Normalised + Log2" = "Normalised + Log2",
                                "Raw" = "Raw",
                                "Normalised" = "Normalised",
                                "TPM" = "TPM",
                                "FPKM" = "FPKM",
                                "VST" = "VST"
                            ),
                            selected = "Normalised + Log2"
                        ),
                        checkboxInput(
                            inputId = "multiBoxplotCountsLog",
                            label = "Apply Log2 transformation?",
                            value = FALSE
                        ),
                        checkboxInput(
                            inputId = "multiBoxplotZScore",
                            label = "Plot Z-scores (gene-wise)?",
                            value = FALSE
                        ),
                        hr(style = "border-top: 1px solid #000000;"),
                        h4("Sample Grouping"),
                        helpText("The grouping factor is automatically determined from the DE analysis."),
                        uiOutput("multiBoxplotBucket"),
                        hr(style = "border-top: 1px solid #000000;"),
                        h4("P-value Options"),
                        checkboxInput(
                            inputId = "multiBoxplotShowP",
                            label = "Show p-values on plot?",
                            value = TRUE
                        ),
                        selectInput(
                            inputId = "multiBoxplotPType",
                            label = "P-value type:",
                            choices = c("FDR" = "fdr", "Raw p-value" = "pValue"),
                            selected = "fdr"
                        ),
                        sliderInput(
                            inputId = "multiBoxplotMaxGenes",
                            label = "Max genes to plot",
                            value = 20,
                            min = 1,
                            max = 50,
                            step = 1,
                            ticks = FALSE
                        )
                    ),
                    tabPanel(
                        title = "Figure Settings",
                        h4("Plot Style"),
                        splitLayout(
                            checkboxInput(
                                inputId = "multiBoxplotShowBox",
                                label = "Show boxplot",
                                value = FALSE
                            ),
                            checkboxInput(
                                inputId = "multiBoxplotShowViolin",
                                label = "Show violin plot",
                                value = FALSE
                            ),
                            checkboxInput(
                                inputId = "multiBoxplotShowPoint",
                                label = "Show individual points",
                                value = TRUE
                            )
                        ),
                        splitLayout(
                            checkboxInput(
                                inputId = "multiBoxplotShowMeanBar",
                                label = "Show mean bar",
                                value = TRUE
                            ),
                            checkboxInput(
                                inputId = "multiBoxplotMeanBarFront",
                                label = "Bring mean bar to front",
                                value = FALSE
                            )
                        ),
                        splitLayout(
                            sliderInput(
                                inputId = "multiBoxplotPointSize",
                                label = "Point size",
                                min = 0.5, max = 5, value = 4, step = 0.25, ticks = FALSE, width = "85%"
                            ),
                            sliderInput(
                                inputId = "multiBoxplotPointAlpha",
                                label = "Point alpha",
                                min = 0, max = 1, value = 0.9, step = 0.05, ticks = FALSE, width = "85%"
                            ),
                            sliderInput(
                                inputId = "multiBoxplotBoxAlpha",
                                label = "Box alpha",
                                min = 0, max = 1, value = 0.2, step = 0.05, ticks = FALSE, width = "85%"
                            )
                        ),
                        sliderInput(
                            inputId = "multiBoxplotMeanLine",
                            label = "Mean bar width",
                            min = 0, max = 1, value = 0.5, step = 0.1, ticks = FALSE, width = "80%"
                        ),
                        hr(style = "border-top: 1px solid #000000;"),
                        h4("P-value Bracket Settings"),
                        checkboxInput(
                            inputId = "multiBoxplotPBlack",
                            label = "Use black brackets",
                            value = TRUE
                        ),
                        splitLayout(
                            sliderInput(
                                inputId = "multiBoxplotPLabelSize",
                                label = "Text size",
                                min = 2, max = 10, value = 5, step = 0.5, ticks = FALSE, width = "85%"
                            ),
                            sliderInput(
                                inputId = "multiBoxplotPBracketSize",
                                label = "Line width",
                                min = 0.1, max = 2, value = 0.6, step = 0.1, ticks = FALSE, width = "85%"
                            ),
                            sliderInput(
                                inputId = "multiBoxplotPTipSize",
                                label = "Tip length",
                                min = 0, max = 0.1, value = 0.02, step = 0.005, ticks = FALSE, width = "85%"
                            )
                        ),
                        checkboxInput(
                            inputId = "multiBoxplotPStack",
                            label = "Stack brackets",
                            value = TRUE
                        ),
                        splitLayout(
                            sliderInput(
                                inputId = "multiBoxplotPOffset",
                                label = "Stack spacing",
                                min = 0, max = 20, value = 8, step = 1, ticks = FALSE, width = "85%"
                            ),
                            sliderInput(
                                inputId = "multiBoxplotPYOffset",
                                label = "Y position (%)",
                                min = 0, max = 50, value = 5, step = 1, ticks = FALSE, width = "85%"
                            ),
                            sliderInput(
                                inputId = "multiBoxplotPWidth",
                                label = "Width (%)",
                                min = 50, max = 100, value = 100, step = 5, ticks = FALSE, width = "85%"
                            )
                        ),
                        helpText("Reduce bracket width to prevent overlapping when brackets are at the same height (non-stacked)."),
                        hr(style = "border-top: 1px solid #000000;"),
                        h4("Layout"),
                        splitLayout(
                            sliderInput(
                                inputId = "multiBoxplotNCol",
                                label = "Columns",
                                value = 3,
                                min = 1,
                                max = 10,
                                step = 1,
                                ticks = FALSE,
                                width = "85%"
                            ),
                            sliderInput(
                                inputId = "multiTextSizeBoxplot",
                                label = "Font size",
                                min = 6, max = 24, value = 12, step = 1, ticks = FALSE, width = "85%"
                            ),
                            sliderInput(
                                inputId = "multiBoxplotConditionAngle",
                                label = "X-axis angle",
                                min = 0, max = 90, value = 45, step = 45, ticks = FALSE, width = "85%"
                            )
                        ),
                        splitLayout(
                            sliderInput(
                                inputId = "multiFigWidthBoxplot",
                                label = "Width (px)",
                                value = 850,
                                min = 400,
                                max = 2400,
                                step = 100,
                                ticks = FALSE,
                                width = "85%"
                            ),
                            sliderInput(
                                inputId = "multiFigHeightBoxplot",
                                label = "Height (px)",
                                value = 600,
                                min = 400,
                                max = 2400,
                                step = 100,
                                ticks = FALSE,
                                width = "85%"
                            )
                        )
                    ),
                    tabPanel(
                        title = "Condition Colours",
                        h4("Group Colours"),
                        helpText("Set colours for each group in your grouping factor. These colours are used for boxplot fills and points."),
                        uiOutput("multiBoxplotGroupColours")
                    )
                )
            )
        ),
        column(
            width = 9,
            box(
                title = "Data Compatibility Check",
                width = NULL,
                solidHeader = TRUE,
                status = "warning",
                collapsible = TRUE,
                collapsed = FALSE,
                uiOutput("boxplotCompatibility")
            ),
            box(
                title = "Multi-Contrast Boxplots",
                width = NULL,
                solidHeader = TRUE,
                status = "success",
                downloadButton(
                    outputId = "dlBoxplot",
                    label = "Download Plot (PDF)"
                ),
                downloadButton(
                    outputId = "dlBoxplotTable",
                    label = "Download Data (Excel)"
                ),
                br(), br(),
                withSpinner(plotOutput(outputId = "boxplotPlot", inline = TRUE), color = "#b5d17b"),
                br(), br(),
                h3("Plot Data Table"),
                DT::dataTableOutput(outputId = "boxplotTable")
            )
        )
    )
)
