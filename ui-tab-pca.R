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
            numericInput(inputId = "pcaTopN", label = "Calculate PCA based on top n features ranked by standard deviation", value = 2000, min = 2, max = 1e5, step = 1),
            checkboxInput(inputId = "pcaCentre", label = "Centre PCA?", value = TRUE),
            checkboxInput(inputId = "pcaScale", label = "Scale PCA?", value = FALSE),
            selectInput(inputId = "pcaX", label = "PC for x axis", choices = "PC1", selected = "PC1"),
            selectInput(inputId = "pcaY", label = "PC for y axis", choices = "PC2", selected = "PC2"),
            selectInput(inputId = "pcaFactor1", label = "Select which factor to colour the samples by", choices = "", selected = ""),
            checkboxGroupInput(inputId = "pcaGroups", label = "Select groups to plot", choices = "", selected = ""),
            selectInput(inputId = "pcaFactor2", label = "Select which factor to shape the samples by", choices = "", selected = ""),
            selectInput(inputId = "pcaCounts", label = "Select count method to plot", choices = "", selected = ""),
            selectInput(inputId = "pcaBatch", label = "Apply a batch-correction to the counts?", choices = "", multiple = TRUE),
            ),
          tabPanel(
            title = "Figure settings",
            tags$b("Show sample names on PCA?"),
            checkboxInput(inputId = "pcaShowNames", label = "Show sample names", value = FALSE),
            tags$b("Keep PCA axes proportional to variance?"),
            checkboxInput(inputId = "pcaAxesProp", label = "Keep axes proportional", value = TRUE),
            checkboxInput(inputId = "showLinesPCA", label = "Show grid lines?", value = TRUE),
            checkboxInput(inputId = "showAxesPCA", label = "Show axes lines?", value = TRUE),
            checkboxInput(inputId = "boldPCA", label = "Use bold font?", value = TRUE),
            numericInput(inputId = "textSizePCA", label = "Figure Font Size", min = 4, max = 30, value = 12, step = 0.5),
            numericInput(inputId = "figWidthPCA", label = "Figure Width", min = 100, max = 2000, value = 800, step = 10),
            numericInput(inputId = "figHeightPCA", label = "Figure Height", min = 100, max = 2000, value = 600, step = 10)
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
        textOutput("pcaDesign"),
        downloadButton(outputId = "dlPCAPlotButton", label = "Download PCA Plot (PDF)"),
        downloadButton(outputId = "dlPCAPlotDFButton", HTML("Download PCA Coords Table (Excel)")),
        br(), br(),
        plotOutput(outputId = "pcaStatic", inline = TRUE, brush = "pcaBrush"),
        tableOutput("pcaBrushTable"),
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