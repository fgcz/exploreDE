tabItem(
  tabName = "pcaPlotTab",
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
          "Principal component analysis plots are a form of dimensionalilty 
          reduction aimed at visualising high-dimensional data, e.g., the
          expression values of thousands of features in many samples,
          in a 2D (or 3D) space. Briefly, each dot represents a sample,
          and the closer dots are to one another on the plot, the more 
          similar those samples are to one another compared to the 
          other dots on the plot. In this case 'more similar' refers to
          the similarity in feature expression counts."
        ),
        tags$p(
          "The percentage of variance explained by each axis, or principal
           component, is relevant. Though axes 1 and 2 are the same length on 
           the plot, the variance explained by each might be very different. 
           This must be considered when comparing the distance between samples 
           along the axes. For example, if PC1 explains 60% of the variance,
           and PC2 explains 20%, the distance between samples along the x-axes,
            PC1, encompasses 40% more difference in feature expression than the 
           distance along the y-axis, PC2."),
        tags$p(
          "The different count transformation methods are explained in some
           detail in the boxplots tab.")
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          tabPanel(
            title = "Main settings",
            helpText("The PCA plot is completely separate from the DE test and will be identical in other live reports from the same dataset, it's just here for convenience."),
            numericInput(inputId = "pcaTopN", label = "Use top n features ranked by SD", value = 2000, min = 2, max = 1e5, step = 1, width = "75%"),
            hr(style = "border-top: 0.1px solid #000000;"),
            selectInput(inputId = "pcaX", label = "PC for x axis", choices = "PC1", selected = "PC1", width = "75%"),
            selectInput(inputId = "pcaY", label = "PC for y axis", choices = "PC2", selected = "PC2", width = "75%"),
            hr(style = "border-top: 0.1px solid #000000;"),
            selectInput(inputId = "pcaFactor1", label = "Colour by", choices = "", selected = "", width = "75%"),
            selectInput(inputId = "pcaFactor2", label = "Shape by", choices = "", selected = "", width = "75%"),
            hr(style = "border-top: 0.1px solid #000000;"),
            selectInput(inputId = "pcaCounts", label = "Counts to plot", choices = "", selected = "", width = "75%"),
            selectInput(inputId = "pcaBatch", label = "Batch-correct by", choices = "", multiple = FALSE, width = "75%"),
            checkboxInput(inputId = "pcaLog2", label = "Log2 counts?", value = TRUE),
            helpText("Log2 only applies to non-logged RNA counts, e.g., TPM, FPKM"),
            hr(style = "border-top: 0.1px solid #000000;"),
            tags$b("Show sample names on PCA?"),
            checkboxInput(inputId = "pcaShowNames", label = "Show sample names", value = FALSE),
            helpText("Not seeing labels? Try changing the max overlap in the figure settings tab."),
            tags$b("Keep PCA axes proportional to variance?"),
            checkboxInput(inputId = "pcaAxesProp", label = "Keep axes proportional", value = FALSE),
            tags$b("Add ellipses?"),
            checkboxInput(inputId = "pcaAddEllipses", label = "Add ellipses", value = TRUE),
            sliderInput(inputId = "pcaEllipsesAlpha", label = "Ellipses alpha", min = 0, max = 1, step = 0.1, value = 0.2, width = "65%"),
            tags$b("Centre and/or scale counts?"),
            checkboxInput(inputId = "pcaCentre", label = "Centre PCA?", value = TRUE),
            checkboxInput(inputId = "pcaScale", label = "Scale PCA?", value = FALSE),
            hr(style = "border-top: 0.1px solid #000000;"),
            checkboxGroupInput(inputId = "pcaGroups", label = "Select groups to plot", choices = "", selected = ""),
            uiOutput(outputId = "nrPeptidesPCAUI")
          ),
          tabPanel(
            title = "Figure settings",
            h4("Plot settings"),
            splitLayout(
              sliderInput(inputId = "figWidthPCA", label = "Width", min = 100, max = 2000, value = 800, step = 10, width = "85%"),
              sliderInput(inputId = "figHeightPCA", label = "Height", min = 100, max = 2000, value = 600, step = 10, width = "85%"),
              sliderInput(inputId = "textSizePCA", label = "Font Size", min = 4, max = 30, value = 12, step = 0.5, width = "85%", ticks = TRUE)
            ),
            splitLayout(
              sliderInput(inputId = "dotSizePCA", label = "Dot size", min = 1, max = 10, value = 5, step = 0.5, width = "85%", ticks = TRUE),
              sliderInput(inputId = "alphaPCA", label = "Dot alpha", min = 0.1, max = 1, value = 0.9, step = 0.1, width = "85%", ticks = TRUE),
              sliderInput(inputId = "dotBorderPCA", label = "Dot border", min = 0, max = 1, value = 0.2, step = 0.1, width = "85%"),
            ),
            br(),
            splitLayout(
              checkboxInput(inputId = "showAxesPCA", label = "Axes lines?", value = TRUE),
              checkboxInput(inputId = "boldPCA", label = "Bold?", value = TRUE),
              checkboxInput(inputId = "showLinesPCA", label = "Grid lines?", value = TRUE)
            ),
            hr(style = "border-top: 1px solid #000000;"), h4("Annotation settings"),
            numericInput(inputId = "pcaLabelMaxOverlap", label = "Number of max overlapping labels", min = 1, max = 1e4, value = 10, step = 1, width = "33%"),
            splitLayout(
              sliderInput(inputId = "geneLabelNudgePCAX", label = "Nudge Labels X", min = -10, max = 10, value = 0, step = 1, width = "85%", ticks = TRUE),
              sliderInput(inputId = "geneLabelNudgePCAY", label = "Nudge Labels Y", min = -10, max = 10, value = 0, step = 1, width = "85%", ticks = TRUE),
              sliderInput(inputId = "geneLabelSizePCA", label = "Label Size", min = 4, max = 30, value = 12, step = 0.5, width = "85%", ticks = TRUE)
            )
          ),
          tabPanel(
            title = "Download settings",
            selectInput(inputId = "downloadFormatPCA", label = "Select format", choices = c("PDF", "SVG", "PNG"), selected = "PDF"),
            selectInput(inputId = "dpiPCA", label = "PNG DPI", choices = c(72, 150, 300, 600, 1000), selected = 600),
            textInput(inputId = "filnamePCA", label = "Enter filename", value = "PCA")
          )
        ),
      ) # end of box
    ), # end of settings column
    column(
      width = 8,
      box(
        title = "PCA Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          tabPanel(
            title = "PCA Plot",
            textOutput("pcaDesign"),
            downloadButton(outputId = "dlPCAPlotButton", label = "Download PCA Plot"),
            downloadButton(outputId = "dlPCAPlotDFButton", HTML("Download PCA Coords Table (Excel)")),
            br(), br(),
            plotOutput(outputId = "pcaStatic", inline = TRUE, brush = "pcaBrush"),
            DT::dataTableOutput("pcaBrushTable"),
            uiOutput(outputId = "pcaNAWarning", inline = TRUE)
          ),
          tabPanel(
            title = "Paired Plot",
            downloadButton(outputId = "dlPCAPairButton", label = "Download Paired Plot"),
            br(),br(),
            plotOutput(outputId = "pcaPaired", inline = TRUE),
          )
        )
      ),
      box(
        title = "Scree Plot",
        width = 6, 
        solidHeader = TRUE,
        status = "primary",
        plotOutput("pcaScree", inline = TRUE),
        DT::dataTableOutput("pcaVars")
      ),
      box(
        title = "PCA Loadings",
        width = 6, 
        solidHeader = TRUE,
        status = "primary",
        DT::dataTableOutput("pcaLoadings"),
        br(), br(),
        downloadButton(outputId = "dlPCALoadDFButton", HTML("Download PCA Loadings Table (Excel)")),
      )
    )
  )
)