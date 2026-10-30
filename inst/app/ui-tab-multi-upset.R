tabItem(
  tabName = "multiUpsetPlotTab",
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
        tags$script(
          "
            Shiny.addCustomMessageHandler('txt', function (txt) {
                navigator.clipboard.writeText(txt);
            });
        "
        ),
        tags$p(
          "Upset plots are a way of visualising overlaps and intersections
          between or across multiple groups. They are similar to Venn diagrams,
          but condense much more information into a single figure. Here, we use
          upset plots to visualise overlapping DE features and pathways between
          multiple SUSHI-generated differential expression results."
        ),
        tags$p(
          "Upset plots have three modes. Please refer to the image at the
        following link, which explains the modes in pictures better than is
           possible in text:"
        ),
        helpText(
          a(
            "Click Here to View Image",
            href = "https://jokergoo.github.io/ComplexHeatmap-reference/book/08-upset_files/figure-html/unnamed-chunk-7-1.png",
            target = "_blank"
          )
        ),
        tags$p(
          "The 'Show single-group sets?' option allows you to either exclude
           (unticked) or include (ticked) the results for individual contrasts.
           If you leave it unticked, the upset plots will only show overlaps
           between two or more groups. Furthermore, when you click groups to see
           the genes/spanathways that overlap for one group, you will not see the
           results if the 'Show single-group sets?' button is unticked."
        ),
        tags$p(
          "Individual parameters for the overlapping DE features and pathways is
           possible. Not all settings affect all figures. For example, log fold
           change doesn't affect the pathway results, while GSEA does not use
           the raw p-values."
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        tabsetPanel(
          tabPanel(
            title = "Main",
            selectInput(
              inputId = "dataToUpset",
              label = "Select between plotting genes or pathways",
              choices = c("Genes", "ORA", "GSEA"),
              selected = "Genes",
              width = "50%"
            ),
            selectInput(
              inputId = "subdataToUpset",
              label = "Select the direction/ontology to plot",
              choices = c(),
              width = "50%"
            ),
            radioButtons(
              inputId = "upsetOrVenn",
              label = "Upset or Venn?",
              choices = c("Upset", "Venn"),
              selected = "Upset",
              inline = TRUE
            ),
            helpText("You can only make a Venn if you have loaded 4 contrasts or fewer."),
            hr(style = "border-top: 1px solid #000000;"),
            h4("DE Settings"),
            numericInput(
              inputId = "lfcUpset",
              label = "Log2 Fold Change Threshold:",
              value = 0,
              min = 0,
              max = 10,
              step = 0.25
            ),
            radioButtons(
              inputId = "pTypeUpset",
              choices = c("FDR", "Raw"),
              label = "P-Value:",
              selected = "FDR",
              inline = TRUE
            ),
            selectInput(
              inputId = "pThresholdUpset",
              choices = c(1.0, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001),
              label = "P-Value Threshold:",
              selected = 0.05
            )
          ),
          tabPanel(
            title = "Upset Settings",
            radioButtons(
              inputId = "upsetMode",
              label = "Upset Mode",
              choices = c("Distinct", "Intersect", "Union"),
              selected = "Distinct"
            ),
            helpText(
              a(
                "Click to see a diagram of the how the modes work",
                href = "https://jokergoo.github.io/ComplexHeatmap-reference/book/08-upset_files/figure-html/unnamed-chunk-7-1.png",
                target = "_blank"
              )
            ),
            numericInput(
              inputId = "upsetMinimum",
              label = "Minimum set size to be displayed",
              value = 1,
              min = 1,
              step = 1
            ),
            checkboxInput(inputId = "upsetSetMode", label = "Show single-group sets?", value = TRUE)
          ),
          tabPanel(
            title = "Figure Settings",
            radioButtons(
              inputId = "upsetColour",
              label = "Figure colour scheme",
              inline = TRUE,
              choiceNames = list(
                tags$span(style = "color:#000000;", "Greyscale"),
                HTML(
                  '<span style="font-weight: bold;"><span style="color:#E00000">C</span><span style="color:#E0E000">o</span><span style="color:#00E000">l</span><span style="color:#00E0E0">o</span><span style="color:#0000E0">u</span><span style="color:#E000E0">r</span></span>'
                )
              ),
              choiceValues = list("greyscale", "colour")
            ),
            numericInput(
              inputId = "textSizeUpset",
              label = "Figure Font Size",
              min = 4,
              max = 30,
              value = 12,
              step = 0.5
            ),
            numericInput(
              inputId = "figWidthUpset",
              label = "Figure Width",
              min = 100,
              max = 2000,
              value = 800,
              step = 10
            ),
            numericInput(
              inputId = "figHeightUpset",
              label = "Figure Height",
              min = 100,
              max = 2000,
              value = 400,
              step = 10
            )
          )
        )
      )
    ),
    column(
      width = 9,
      box(
        title = "Upset Plot",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        h4("Download Settings"),
        splitLayout(
          selectInput(
            inputId = "downloadFormatMultiUpset",
            label = "Format",
            choices = c("PDF", "SVG", "PNG"),
            selected = "PDF"
          ),
          selectInput(inputId = "dpiMultiUpset", label = "DPI", choices = c(72, 150, 300, 600, 1000), selected = 600)
        ),
        textInput(inputId = "filename_multiUpset", label = "Filename", value = "MultiUpset"),
        downloadButton(
          outputId = "dlPlot_multiUpset",
          label = "Download Upset/Venn Plot"
        ),
        br(),
        uiOutput(outputId = "upsetError", inline = TRUE),
        plotOutput(outputId = "UpsetPlot", inline = TRUE),
        br(),
        br(),
        uiOutput(outputId = "UpsetTableGenesUI"),
        actionButton("copy_link", "Copy to clipboard"),
        uiOutput("linkToEnrichr"),
        uiOutput(outputId = "UpsetIntersectGenesUI")
      )
    )
  )
)
