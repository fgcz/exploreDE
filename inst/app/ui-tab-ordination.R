tabItem(
  tabName = "ordinationTab",
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
          "Ordination methods are forms of dimensionality reduction aimed at
          visualising high-dimensional data (e.g., expression values of thousands
          of features in many samples) in a 2D space. Each dot represents a sample,
          and the closer dots are to one another on the plot, the more similar
          those samples are to one another."
        ),
        tags$p(
          tags$b("PCA"), " (Principal Component Analysis) works best with normally
          distributed data and is sensitive to outliers. ",
          tags$b("PCoA/MDS"), " (Principal Coordinates Analysis / Multidimensional Scaling)
          can use different distance metrics and is more flexible. ",
          tags$b("NMDS"), " (Non-metric Multidimensional Scaling) preserves rank-order
          distances and is robust to non-linear relationships. ",
          tags$b("GLM-PCA"), " uses generalized linear models and is specifically designed
          for count data (RNA-seq)."
        ),
        tags$p(
          "The percentage of variance explained by each axis (or stress for NMDS)
          indicates how well the 2D representation captures the original data structure.
          The different count transformation methods are explained in the boxplots tab."
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
            helpText("The ordination plot is completely separate from the DE test and will be identical in other live reports from the same dataset, it's just here for convenience."),

            # Method selector
            selectInput(
              inputId = "ordinationMethod",
              label = "Ordination method",
              choices = c("PCA" = "pca", "PCoA/MDS" = "pcoa", "NMDS" = "nmds", "GLM-PCA" = "glmpca"),
              selected = "pca",
              width = "75%"
            ),
            helpText("PCA: standard linear method. PCoA/MDS: distance-based. NMDS: non-metric. GLM-PCA: for count data."),
            hr(style = "border-top: 0.1px solid #000000;"),

            # Plot type selector (samples only or biplot)
            radioButtons(
              inputId = "ordinationPlotType",
              label = "Plot type",
              choices = c("Samples only" = "samples", "Biplot (samples + genes)" = "biplot"),
              selected = "samples",
              width = "75%"
            ),

            # Biplot-specific controls (conditionally shown)
            conditionalPanel(
              condition = "input.ordinationPlotType == 'biplot'",
              radioButtons(
                inputId = "biplotGeneSelectionMethod",
                label = "Gene selection method",
                choices = c(
                  "Top N genes (overall)" = "overall",
                  "Top N per PC" = "per_pc",
                  "Manual selection" = "manual"
                ),
                selected = "overall",
                width = "75%"
              ),

              # Controls for "overall" method
              conditionalPanel(
                condition = "input.biplotGeneSelectionMethod == 'overall'",
                numericInput(
                  inputId = "biplotNGenes",
                  label = "Number of top genes to show",
                  value = 10,
                  min = 1,
                  max = 100,
                  step = 1,
                  width = "75%"
                )
              ),

              # Controls for "per_pc" method
              conditionalPanel(
                condition = "input.biplotGeneSelectionMethod == 'per_pc'",
                numericInput(
                  inputId = "biplotNGenesPerPC",
                  label = "Number of genes per PC",
                  value = 5,
                  min = 1,
                  max = 50,
                  step = 1,
                  width = "75%"
                ),
                helpText("Selects top N genes from each selected PC based on absolute loading values")
              ),

              # Controls for "manual" method
              conditionalPanel(
                condition = "input.biplotGeneSelectionMethod == 'manual'",
                selectizeInput(
                  inputId = "biplotCustomGenes",
                  label = "Select genes:",
                  multiple = TRUE,
                  choices = NULL,
                  selected = NULL,
                  options = list(
                    placeholder = "Search and select genes",
                    plugins = list("remove_button", "drag_drop", "restore_on_backspace", "clear_button")
                  ),
                  width = "75%"
                ),
                textAreaInput(
                  inputId = "biplotCustomGenesText",
                  label = "Or paste gene list:",
                  placeholder = "gene1 gene2 gene3",
                  rows = 2,
                  width = "75%"
                ),
                helpText("Enter gene names separated by spaces, commas, or newlines")
              ),
              sliderInput(
                inputId = "biplotArrowScale",
                label = "Arrow length scale",
                min = 0.1,
                max = 3,
                value = 1,
                step = 0.1,
                width = "75%"
              ),

              # Single arrow color for overall and manual methods
              conditionalPanel(
                condition = "input.biplotGeneSelectionMethod != 'per_pc'",
                colourpicker::colourInput(
                  inputId = "biplotArrowColour",
                  label = "Arrow colour",
                  value = "#FF7F00",
                  showColour = "background"
                )
              ),

              # Two arrow colors for per_pc method
              conditionalPanel(
                condition = "input.biplotGeneSelectionMethod == 'per_pc'",
                colourpicker::colourInput(
                  inputId = "biplotArrowColourX",
                  label = "Arrow colour (PC X genes)",
                  value = "#FF7F00",
                  showColour = "background"
                ),
                colourpicker::colourInput(
                  inputId = "biplotArrowColourY",
                  label = "Arrow colour (PC Y genes)",
                  value = "#984EA3",
                  showColour = "background"
                )
              )
            ),
            hr(style = "border-top: 0.1px solid #000000;"),
            numericInput(inputId = "pcaTopN", label = "Use top n features ranked by SD", value = 2000, min = 2, max = 1e5, step = 1, width = "75%"),

            # Method-specific controls
            conditionalPanel(
              condition = "input.ordinationMethod == 'pcoa' || input.ordinationMethod == 'nmds'",
              selectInput(
                inputId = "pcoaDistance",
                label = "Distance metric",
                choices = c(
                  "Euclidean" = "euclidean",
                  "Bray-Curtis" = "bray",
                  "Canberra" = "canberra",
                  "Manhattan" = "manhattan",
                  "Jaccard" = "jaccard"
                ),
                selected = "euclidean",
                width = "75%"
              ),
              helpText("Euclidean: continuous data. Bray-Curtis: abundance data. Canberra: sensitive to small values.")
            ),
            conditionalPanel(
              condition = "input.ordinationMethod == 'glmpca'",
              selectInput(
                inputId = "glmpcaFamily",
                label = "GLM family",
                choices = c(
                  "Poisson" = "poi",
                  "Negative Binomial" = "nb",
                  "Binomial" = "binom"
                ),
                selected = "poi",
                width = "75%"
              ),
              helpText("Poisson: simple counts. Negative Binomial: overdispersed counts. Binomial: proportion data."),
              numericInput(
                inputId = "glmpcaMinibatch",
                label = "Minibatch size (for large datasets)",
                value = 100,
                min = 10,
                max = 1000,
                step = 10,
                width = "75%"
              )
            ),
            hr(style = "border-top: 0.1px solid #000000;"),
            selectInput(inputId = "pcaX", label = "Axis for x", choices = "PC1", selected = "PC1", width = "75%"),
            selectInput(inputId = "pcaY", label = "Axis for y", choices = "PC2", selected = "PC2", width = "75%"),
            hr(style = "border-top: 0.1px solid #000000;"),
            selectInput(inputId = "pcaFactor1", label = "Colour by", choices = "", selected = "", width = "75%"),
            selectInput(inputId = "pcaFactor2", label = "Shape by", choices = "", selected = "", width = "75%"),
            hr(style = "border-top: 0.1px solid #000000;"),
            selectInput(inputId = "pcaCounts", label = "Counts to plot", choices = "", selected = "", width = "75%"),
            selectInput(inputId = "pcaBatch", label = "Batch-correct by", choices = "", multiple = FALSE, width = "75%"),
            checkboxInput(inputId = "pcaLog2", label = "Log2 counts?", value = TRUE),
            helpText("Log2 only applies to non-logged RNA counts, e.g., TPM, FPKM"),
            hr(style = "border-top: 0.1px solid #000000;"),
            tags$b("Show sample names?"),
            checkboxInput(inputId = "pcaShowNames", label = "Show sample names", value = FALSE),
            helpText("Not seeing labels? Try changing the max overlap in the figure settings tab."),
            tags$b("Keep axes proportional to variance?"),
            checkboxInput(inputId = "pcaAxesProp", label = "Keep axes proportional", value = FALSE),
            tags$b("Add ellipses?"),
            checkboxInput(inputId = "pcaAddEllipses", label = "Add ellipses", value = TRUE),
            sliderInput(inputId = "pcaEllipsesAlpha", label = "Ellipses alpha", min = 0, max = 1, step = 0.1, value = 0.2, width = "65%"),

            # PCA-specific options
            conditionalPanel(
              condition = "input.ordinationMethod == 'pca'",
              tags$b("Centre and/or scale data?"),
              checkboxInput(inputId = "pcaCentre", label = "Centre data?", value = TRUE),
              checkboxInput(inputId = "pcaScale", label = "Scale data?", value = FALSE),
              helpText("Centering removes the mean, scaling divides by standard deviation.")
            ),
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
        title = "Ordination Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          tabPanel(
            title = "Ordination Plot",
            textOutput("pcaDesign"),
            downloadButton(outputId = "dlPCAPlotButton", label = "Download Plot"),
            downloadButton(outputId = "dlPCAPlotDFButton", HTML("Download Coords Table (Excel)")),
            br(), br(),
            shinycssloaders::withSpinner(
              plotOutput(outputId = "pcaStatic", brush = "pcaBrush", width = "100%", height = "600px"),
              type = 4,
              color = "#0275D8"
            ),
            uiOutput(outputId = "pcaNAWarning", inline = TRUE),
            br(), br(),
            DT::dataTableOutput("pcaBrushTable")
          ),
          tabPanel(
            title = "Paired Plot",
            downloadButton(outputId = "dlPCAPairButton", label = "Download Paired Plot"),
            br(), br(),
            plotOutput(outputId = "pcaPaired", inline = TRUE),
          )
        )
      ),
      box(
        title = "Variance Explained / Stress",
        width = 6,
        solidHeader = TRUE,
        status = "primary",
        plotOutput("pcaScree", inline = TRUE),
        DT::dataTableOutput("pcaVars")
      ),
      box(
        title = "Feature Loadings",
        width = 6,
        solidHeader = TRUE,
        status = "primary",
        DT::dataTableOutput("pcaLoadings"),
        br(), br(),
        downloadButton(outputId = "dlPCALoadDFButton", HTML("Download Loadings Table (Excel)")),
      )
    )
  )
)
