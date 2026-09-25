tabItem(
    tabName = "multiBubbleTab",
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
                    "Similar to the DE tile plot, contrasts are columns, and pathways are rows. The top n pathways, based on p-value, have been plotted. In ORA mode, the colour of the dot is the contrast, and in GSEA mode, the colour is the
          normalised enrichment score. The size of the dot is the number of genes in that pathway (DE in ORA mode, total set size in GSEA mode). The rowname colour denotes which contrast the pathway came from. If it is black, it came from
          multiple contrasts."
                )
            ),
            box(
                title = "Settings",
                width = NULL,
                solidHeader = TRUE,
                status = "success",
                selectInput(inputId = "dataToBubble", label = "Select between plotting ORA or GSEA Results", choices = c("ORA", "GSEA"), selected = "ORA", width = "50%"),
                selectInput(inputId = "subdataToBubble", label = "Select the direction/ontology to plot", choices = c(), width = "50%"),
                hr(style = "border-top: 1px solid #000000;"), h4("Input Settings"),
                # The number of genes displayed on the bubble map
                numericInput(
                    inputId = "bubblePathwayNumber",
                    label = "Number Top Pathways from each contrast",
                    min = 5, max = 5000, value = 5, step = 1
                ),
                tags$b("Turn off/on names of rows/columns"),
                # Display colnames?:
                checkboxInput(
                    inputId = "colnamesBubble",
                    label = "Show contrast names?",
                    value = TRUE
                ),
                hr(style = "border-top: 1px solid #000000;"), h4("Figure Settings"),
                numericInput(
                    inputId = "textSizeBubble",
                    label = "Figure Font Size", min = 4, max = 30,
                    value = 12, step = 0.5
                ),
                numericInput(
                    inputId = "figWidthBubble",
                    label = "Figure Width", min = 100, max = 2000,
                    value = 800, step = 10
                ),
                numericInput(
                    inputId = "figHeightBubble",
                    label = "Figure Height", min = 100, max = 2000,
                    value = 800, step = 10
                )
            )
        ),
        column(
            width = 9,
            box(
                title = "Bubble plots",
                width = NULL,
                solidHeader = TRUE,
                status = "success",
                h4("Download Settings"),
                splitLayout(
                  selectInput(inputId = "downloadFormatMultiBubble", label = "Format", choices = c("PDF", "SVG", "PNG"), selected = "PDF"),
                  selectInput(inputId = "dpiMultiBubble", label = "DPI", choices = c(72, 150, 300, 600, 1000), selected = 600)
                ),
                textInput(inputId = "filename_multiBubble", label = "Filename", value = "MultiBubble"),
                downloadButton(
                    outputId = "dlPlot_multiBubble",
                    label = "Download Bubble Plot"
                ),
                br(),
                plotOutput(outputId = "BubblePlot", inline = TRUE)
            ),
            box(
                title = "Bubble table",
                width = NULL,
                solidHeader = TRUE,
                status = "success",
                downloadButton(
                    outputId = "dlTable_multiBubble",
                    label = "Download Bubble Table (Excel)"
                ),
                br(),
                DT::dataTableOutput(outputId = "BubbleTable")
            )
        )
    )
)
