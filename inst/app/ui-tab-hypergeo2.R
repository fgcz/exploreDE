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
        title = "ORA Tables",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        textOutput("oraDesign"), br(),
        downloadButton(
          outputId = "dlTableORA",
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
        title = "ORA Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        helpText("Click on pathways in the table on the left for them to be added to the figures here."),
        tabsetPanel(
          tabPanel(
            title = "Network plot",
            downloadButton(outputId = "dlCnetPlot_ORA", label = paste("Download Network Plot: PDF")),
            br(), br(),
            plotOutput(outputId = "cnetPlot_ORA", inline = TRUE),
            br(), 
            checkboxInput(inputId = "showGeneLabelsORA", label = "Show gene names?", value = TRUE)
          ),
          tabPanel(
            title = "Bar plot",
            downloadButton(outputId = "dlBarPlot_ORA", label = paste("Download Bar Plot: PDF")),
            br(), br(),
            plotOutput(outputId = "barPlot_ORA", inline = TRUE)
          ),
          tabPanel(
            title = "Dot plot",
            downloadButton(outputId = "dlDotPlot_ORA", label = paste("Download Dot Plot: PDF")),
            br(), br(),
            plotOutput(outputId = "dotPlot_ORA", inline = TRUE)
          ),
          tabPanel(
            title = "Heatmap",
            downloadButton(outputId = "dlHeatmap_ORA", label = paste("Download Heatmap: PDF")),
            br(), br(),
            plotOutput(outputId = "heatmap_ORA", inline = TRUE)
          ),
          tabPanel(
            title = "Upset plot",
            downloadButton(outputId = "dlUpset_ORA", label = paste("Download Upset Plot: PDF")),
            br(), br(),
            plotOutput(outputId = "upset_ORA", inline = TRUE)
          )
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = FALSE,
        status = "primary",
        column(
          width = 4,
          h4("Network Plot Settings"),
          sliderInput(inputId = "scaleLimORA", label = "Scale colour limit (+/-)", min = 1, max = 10, value = 5, step = 0.5, width = "85%"),
          sliderInput(inputId = "nodeSizeORA", label = "Network Plot Node Size", min = 1, max = 30, value = 9, step = 0.5, width = "85%"),
          sliderInput(inputId = "nodeBorderORA", label = "Node border", min = 0, max = 2, value = 0.3, step = 0.1, width = "80%"),
        ),
        column(
          width = 4,
          h4("Label Settings"),
          sliderInput(inputId = "oraMaxOver", label = "Label max overlaps", min = 0, max = 100, value = 10, step = 1, width = "85%"),
          sliderInput(inputId = "oraDodgeY", label = "Label dodge y", min = -10, max = 10, value = 0, step = 1, width = "85%"),
          sliderInput(inputId = "oraDodgeX", label = "Label dodge x", min = -10, max = 10, value = 0, step = 1, width = "85%"),
          sliderInput(inputId = "oraLabelAlpha", label = "Label alpha", min = 0, max = 1, value = 0.7, step = 0.1, width = "85%")
        ),
        column(
          width = 4,
          h4("Size Settings"),
          sliderInput(inputId = "textSizeORA", label = "Figure Font Size", min = 4, max = 30, value = 10, step = 0.5, width = "85%"),
          sliderInput(inputId = "figHeightORA", label = "Figure Height", min = 100, max = 2000, value = 400, step = 10, width = "85%"),
          sliderInput(inputId = "figWidthORA", label = "Figure Width", min = 100, max = 2000, value = 600, step = 10, width = "85%")
        )
      )
    )
  )
)