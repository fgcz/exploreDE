tabItem(
  tabName = "correlationsTab",
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
          "Pick two features, pick a count method, pick a correlation method... et voilà!"
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          tabPanel(
            title = "Main settings",
            selectizeInput(inputId = "correlationGene1", label = "Feature X", multiple = FALSE, choices = NULL, selected = NULL),
            selectizeInput(inputId = "correlationGene2", label = "Feature Y", multiple = FALSE, choices = NULL, selected = NULL),
            selectInput(inputId = "correlationColourBy", label = "Select column to colour by", choices = NULL, selected =  NULL),
            selectInput(inputId = "correlationShapeBy", label = "Select column to shape by", choices = NULL, selected =  NULL),
            selectInput(inputId = "correlationCounts", label = "Select count method to plot", choices = NULL, selected = NULL),
            helpText("Only applies to counts that are not already in log space."),
            checkboxInput(inputId = "correlationCountsLog", label = "Log2 Counts?", value = FALSE),
            selectInput(inputId = "correlationBatch", label = "Apply a batch-correction to the counts?", choices = NULL, selected = NULL, multiple = FALSE),
            
            hr(style = "border-top: 1px solid #000000;"), h4("Correlation Settings"),
            selectInput(inputId = "correlationMethod", label = "Correlation method", choices = c("Spearman", "Pearson", "Kendall"), selected = "Pearson", multiple = FALSE),
            fluidRow(
              column(6, radioButtons("correlationCorPosX", "Label position x", choices = c("Left", "Right"), selected = "Left", inline = TRUE)),
              column(6, radioButtons("correlationCorPosY", "Label position y", choices = c("Top", "Bottom"), selected = "Top", inline = TRUE))
            ),
          ),
          tabPanel(
            title = "Figure settings",
            splitLayout(
              sliderInput(inputId = "correlationPointSize", label = "Point size", min = 1, max = 10, value = 4, step = 0.5, ticks = FALSE, width = "85%"),
              sliderInput(inputId = "correlationPointStroke", label = "Point border", min = 0, max = 10, value = 0.6, step = 0.1, ticks = FALSE, width = "85%"),
              sliderInput(inputId = "correlationPointAlpha", label = "Point alpha", min = 0, max = 1, value = 1, step = 0.1, ticks = FALSE, width = "85%")
            ),
            splitLayout(
              sliderInput(inputId = "textSizeCorrelation", label = "Figure Font Size", min = 4, max = 30, value = 14, step = 0.5, ticks = FALSE, width = "85%"),
              sliderInput(inputId = "figWidthCorrelation", label = "Figure Width", min = 100, max = 2000, value = 600, step = 10, ticks = FALSE, width = "85%"),
              sliderInput(inputId = "figHeightCorrelation", label = "Figure Height", min = 100, max = 2000, value = 450, step = 10, ticks = FALSE, width = "85%")
            ),
            checkboxInput(inputId = "correlationShowNames", label = "Show sample names", value = FALSE),
            helpText("Increase the overlaps to show more labels"),
            splitLayout(
              sliderInput(inputId = "correlationMaxOverlaps", label = "Label overlaps", min = 1, max = 50, value = 10, step = 1, ticks = FALSE, width = "85%"),
              sliderInput(inputId = "correlationLabelPosX", label = "Nudge label x", min = 0, max = 20, value = 0, step = 1, ticks = FALSE, width = "85%"),
              sliderInput(inputId = "correlationLabelPosY", label = "Nudge label y", min = 0, max = 20, value = 0, step = 1, ticks = FALSE, width = "85%")
            ),
            uiOutput("correlationGroupBucket")
          )
        )
      )
    ),
    column(
      width = 8,
      box(
        title = "Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        downloadButton(outputId = "dlCorrelationButton", label = paste("Download Correlation (PDF)")),
        br(), br(),
        plotOutput(outputId = "correlationStatic", inline = TRUE)
      )
    )
  )
)