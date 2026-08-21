# ORA tab -------------------------------------------------------------
tabItem(
  tabName = "oraTab",
  fluidRow(
    column(
      width = 5,
      box(
        title = "Info",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        collapsed = TRUE,
        tags$p(
          "The hypergeometric over-representation analysis (ORA) test gives an 
            estimate of whether a set of selected genes is enriched for 
            genes in specific categories. It takes as input a subset of 
            genes passing certain p-value/fold-change thresholds (i.e., 
            differentially expressed genes). It is recommended when the 
            difference between groups is large (e.g., 500+ genes above 
            p-value/fold-change thresholds)"
        ),
        tags$p(
          "The ORA results are based on Gene Ontology, and are divided first 
            into three ontologies: Biological Processes; Molecular Function, 
            and; Cellular Components. Within each ontology are three 
            directions: Both, where genes with both positive and negative
            log fold changes are input; Up, where genes with positive log
            fold changes only are input, and; Down, where genes with 
            negative log fold changes only are input."
        ),
        tags$p(
          "The main results table shows the Gene Ontology (GO) ID and
            description. The Gene Ratio shows the number of genes in the
            input set that are linked to that GO term, and the total size
            of the input gene set. The genes from each pathway can be pasted
            directly into the text boxes of other tabs, e.g., boxplots or 
            volcano plots. It is possible to search for specific genes to see 
            which pathways they are associated with by entering the gene name
            in the search box under the geneName column."
        ),
        tags$p(
          "By default, the network plot and barplot are absent. GO terms can be 
            selected from the results table by clicking on them. As you click 
            them, the row will be highlighted blue, and the GO term
            will appear in the network plot and barplot automatically. 
            Downloading the network plot will save a PDF version of the
            plot as it is currently displayed"
        ),
        tags$p(
          "To generate the list of DEGs that were input into the ORA test, and 
          that are visualised in any of the figures on the ORA tab, a raw p-value
          threshold of 0.01 and no log fold change threshold, were applied."
        )
      ),
      box(
        title = "Results Summary",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        collapsed = TRUE,
        tableOutput(outputId = "oraThreshOverview"),
        tableOutput(outputId = "oraOverview")
      ),
      box(
        title = "ORA Tables",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        textOutput("oraDesign"), br(),
        downloadButton(
          outputId = "dlTable_ORA",
          label = "Download ORA Results File"
        ),
        br(), br(),
        splitLayout(
          selectInput(inputId = "oraType", label = "Select GO to view", choices = c("BP", "MF", "CC"), selected = "BP", width = "85%"),
          selectInput(inputId = "oraDirection", label = "Select Direction to view", choices = c("upGenes", "downGenes", "bothGenes"), selected = "upGenes", width = "85%")
        ),
        DT::dataTableOutput(outputId = "selectedTable_ORA"),
        style = "overflow-y: scroll;"
      )),
    column(
      width = 7,
      box(
        title = "ORA Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        helpText("Click on pathways in the table on the left for them to be added to the figures here."),
        tabsetPanel(
          tabPanel(
            title = "Network plot",
            withSpinner(plotOutput(outputId = "cnetPlot_ORA", inline = TRUE), color = "#6E8D9E"),
            br(),
            downloadButton(outputId = "dlPlot_cnetORA", label = paste("Download Network Plot")),
            br(), 
            checkboxInput(inputId = "showGeneLabelsORA", label = "Show gene names?", value = TRUE),
            column(
              width = 6,
              h4("Network Plot Settings"),
              sliderInput(inputId = "scaleLimORA", label = "Scale colour limit (+/-)", min = 1, max = 10, value = 5, step = 0.5, width = "85%"),
              sliderInput(inputId = "nodeSizeORA", label = "Network Plot Node Size", min = 1, max = 30, value = 9, step = 0.5, width = "85%"),
              sliderInput(inputId = "nodeBorderORA", label = "Node border", min = 0, max = 2, value = 0.3, step = 0.1, width = "80%"),
            ),
            column(
              width = 6,
              h4("Label Settings"),
              sliderInput(inputId = "oraMaxOver", label = "Label max overlaps", min = 0, max = 100, value = 10, step = 1, width = "85%"),
              sliderInput(inputId = "oraDodgeY", label = "Label dodge y", min = -10, max = 10, value = 0, step = 1, width = "85%"),
              sliderInput(inputId = "oraDodgeX", label = "Label dodge x", min = -10, max = 10, value = 0, step = 1, width = "85%"),
              sliderInput(inputId = "oraLabelAlpha", label = "Label alpha", min = 0, max = 1, value = 0.7, step = 0.1, width = "85%")
            ),
          ),
          tabPanel(
            title = "Bar plot",
            withSpinner(plotOutput(outputId = "barPlot_ORA", inline = TRUE), color = "#6E8D9E"),
            br(),
            downloadButton(outputId = "dlPlot_barORA", label = paste("Download Bar Plot"))
          ),
          tabPanel(
            title = "Dot plot",
            withSpinner(plotOutput(outputId = "dotPlot_ORA", inline = TRUE), color = "#6E8D9E"),
            br(),
            downloadButton(outputId = "dlPlot_dotORA", label = paste("Download Dot Plot"))
          ),
          tabPanel(
            title = "Heatmap",
            withSpinner(plotOutput(outputId = "heatmap_ORA", inline = TRUE), color = "#6E8D9E"),
            br(),
            downloadButton(outputId = "dlPlot_heatmapORA", label = paste("Download Heatmap"))
          ),
          tabPanel(
            title = "Upset plot",
            withSpinner(plotOutput(outputId = "upset_ORA", inline = TRUE), color = "#6E8D9E"),
            br(),
            downloadButton(outputId = "dlPlot_upsetORA", label = paste("Download Upset Plot"))
          ),
          tabPanel(
            title = "Lollipop plot",
            withSpinner(plotOutput(outputId = "lollipopPlot_ORA", inline = TRUE), color = "#6E8D9E"),
            downloadButton(outputId = "dlPlot_lollipopORA", label = paste("Download Lollipop Plot")),
            br(),
            colourpicker::colourInput(inputId = "oraLollipopColour", label = "Dot/stem colour", value = "steelblue4", closeOnClick = TRUE, returnName = TRUE, width = "35%"),
            sliderInput(inputId = "oraLollipopScale", label = "Lollipop scale", min = 0.25, max = 5, value = 1, step = 0.25, width = "35%", ticks = FALSE)
          )
        )
      ),
      uiOutput("oraIsoformAggUI"),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        status = "primary",
        column(
          width = 4,
          h4("Plot sizes"),
          sliderInput(inputId = "textSizeORA", label = "Figure Font Size", min = 4, max = 30, value = 10, step = 0.5, width = "85%"),
          sliderInput(inputId = "figHeightORA", label = "Figure Height", min = 100, max = 2000, value = 400, step = 10, width = "85%"),
          sliderInput(inputId = "figWidthORA", label = "Figure Width", min = 100, max = 2000, value = 600, step = 10, width = "85%"),
        ),
        column(
          width = 4,
          h4("Download settings"),
          selectInput(inputId = "downloadFormatORA", label = "Select format", choices = c("PDF", "SVG", "PNG"), selected = "PDF", width = "85%"),
          selectInput(inputId = "dpiORA", label = "PNG DPI", choices = c(72, 150, 300, 600, 1000), selected = 600, width = "85%"),
          textInput(inputId = "filename_ORA", label = "Enter filename", value = "ORA", width = "85%")
        )
      ),
      box(
        title = "Gene Filtering (Visualisation Only)",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        status = "warning",
        helpText("IMPORTANT: These thresholds only affect visualisation (network plots and heatmaps). They do NOT change the ORA results, which were pre-calculated."),
        column(
          width = 6,
          h4("Threshold Settings"),
          selectInput(inputId = "oraVisPType", label = "P-Value Type:", choices = c("FDR", "Raw"), selected = "Raw", width = "85%"),
          numericInput(inputId = "oraVisPThreshold", label = "P-Value Threshold:", value = 0.05, min = 0, max = 1, step = 0.01, width = "85%"),
          numericInput(inputId = "oraVisLFC", label = "Log2 Fold Change Threshold:", value = 0, min = 0, max = 10, step = 0.25, width = "85%")
        ),
        column(
          width = 6,
          h4("Filter Mode"),
          radioButtons(
            inputId = "oraVisFilterMode",
            label = "How to apply filtering:",
            choices = c("Show all genes" = "none", "Filter labels only" = "labels", "Filter genes entirely" = "genes"),
            selected = "none"
          ),
          helpText("'Show all genes': Display all genes, label all genes."),
          helpText("'Filter labels only': Display all genes, but only label those meeting thresholds."),
          helpText("'Filter genes entirely': Only display/show genes meeting thresholds.")
        )
      )
    )
  )
)