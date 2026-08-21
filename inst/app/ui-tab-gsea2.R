# GSEA tab -------------------------------------------------------------
tabItem(
  tabName = "gseaTab",
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
          "Gene Set Enrichment Analysis (GSEA) calculates an enrichment score
            for each annotation category (e.g., those in GO BP) by screening all
            the genes from a differential expression analysis and their
            associated fold-changes. It does not require a pre-selection based
            on p-value/fold-change. It is recommended when the difference
            between groups is small (i.e., applying thresholds would result in
            very few genes selected) or when combining results from different
            experiments."
        ),
        tags$p(
          "Similar to ORA, the GSEA is based on GO and split into the three
            ontologies: BP; MF, and; CC. However, since GSEA is ranked based on
            log fold change (LFC), there are no seperate results for specific
            LFC direction.
            Instead, direction of the pathway is represented by the enrichment
            score, where a positive enrichment score implies the pathway is
            comprised principally of genes with positive log fold change DEGs,
            and a negative enrichment score implies the pathway is mostly
            negative DEGs."
        ),
        tags$p(
          "Again, similarly to the ORA, all figures are absent by default until
            GO terms are clicked in the results table."
        ),
        h4("Running Score Plot"),
        tags$p(
          "When one pathway is selected, the green line shows the
          enrichment running score as a function of each gene within the pathway.
          The running enrichment score is based on the rank of each gene, i.e.,
          the more positive or negative the LFC of the gene, the higher the rank.
          So a higher rank gene has more of an impact on the score, e.g., genes
          at either end of the red/blue spectrum will cause the score to
          increase/decrease gretaer than genes in the middle.
          The black lines directly underneath denote each gene in the pathway.
          The position of the black line on the x-axis denotes the gene's rank.
          The histogram underneath shows how the distribution of the DEGs and
          their ranks.
          When multiple pathways are selected, the running score and gene rank
          lines become coloured by pathway."
        ),
        h4("Ridge Plot"),
        tags$p(
          "Analogous to the barplot in the ORA results, however here
          we can see each GO ID as a histogram showing the distribution of genes
          and their enrichment score. Peaks on the ridge plot reflect the number
          of genes accumulated at that enrichment score."
        ),
        h4("Upset plot"),
        tags$p(
          "These are like venn diagrams, that show the overlap
          between the DEGs in the pathways selected.
          The plot is split into two parts, a violin plot at the top, and a
          dot and line plot at the bottom.
          The dots and lines at the bottom show which group/s are being
          presented in the violin plot at the top. So, for example, where we see
          a line between two groups, known as a set, the violin plot above shows
          the LFCs of the overlapping DEGs between these two groups.
          The equivalent would be the overlapping section of a venn diagram.
          Where there is a dot and no line, then the DEGs are unique to that
          group. This is the same as the section of a venn that doesn't
          overlap. The sets are in size order, so largest sets on the left. It
          is hopefully a useful way to visualise the amount of overlap in genes
          between multiple GO terms."
        )
      ),
      box(
        title = "Results Summary",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        collapsed = TRUE,
        tableOutput(outputId = "gseaThreshOverview"),
        tableOutput(outputId = "gseaOverview")
      ),
      box(
        title = "GSEA Tables",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        textOutput("gseaDesign"), br(),
        downloadButton(
          outputId = "dlTable_GSEA",
          label = "Download GSEA Results File"
        ),
        br(), br(),
        selectInput(inputId = "gseaType", label = "Select result to view", choices = NA, selected = NA, width = "25%"),
        DT::dataTableOutput(outputId = "selectedTable_GSEA"),
        style = "overflow-y: scroll;"
      )),
    column(
      width = 7,
      box(
        title = "GSEA Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          tabPanel(
            title = "Network plot",
            withSpinner(plotOutput(outputId = "cnetPlot_GSEA", inline = TRUE), color = "#6E8D9E"),
            br(),
            downloadButton(outputId = "dlPlot_cnetGSEA", label = paste("Download Network Plot")),
            br(),
            checkboxInput(inputId = "showGeneLabelsGSEA", label = "Show gene names?", value = TRUE),
            column(
              width = 6,
              h4("Network Plot Settings"),
              sliderInput(inputId = "scaleLimGSEA", label = "Scale colour limit (+/-)", min = 1, max = 10, value = 5, step = 0.5, width = "85%"),
              sliderInput(inputId = "nodeSizeGSEA", label = "Network Plot Node Size", min = 1, max = 30, value = 9, step = 0.5, width = "85%"),
              sliderInput(inputId = "nodeBorderGSEA", label = "Node border", min = 0, max = 2, value = 0.3, step = 0.1, width = "80%"),
            ),
            column(
              width = 6,
              h4("Label Settings"),
              sliderInput(inputId = "gseaMaxOver", label = "Label max overlaps", min = 0, max = 100, value = 10, step = 1, width = "85%"),
              sliderInput(inputId = "gseaDodgeY", label = "Label dodge y", min = -10, max = 10, value = 0, step = 1, width = "85%"),
              sliderInput(inputId = "gseaDodgeX", label = "Label dodge x", min = -10, max = 10, value = 0, step = 1, width = "85%"),
              sliderInput(inputId = "gseaLabelAlpha", label = "Label alpha", min = 0, max = 1, value = 0.7, step = 0.1, width = "85%")
            ),
          ),
          tabPanel(
            title = "Heatmap",
            withSpinner(plotOutput(outputId = "heatmap_GSEA", inline = TRUE), color = "#6E8D9E"),
            br(),
            downloadButton(outputId = "dlPlot_heatmapGSEA", label = paste("Download Heatmap"))
          ),
          tabPanel(
            title = "Running Score",
            withSpinner(plotOutput(outputId = "runningScore_GSEA", inline = TRUE), color = "#6E8D9E"),
            br(),
            downloadButton(outputId = "dlPlot_runningScoreGSEA", label = paste("Download Running Score Plot"))
          ),
          tabPanel(
            title = "Ridge Plot",
            withSpinner(plotOutput(outputId = "ridgePlot_GSEA", inline = TRUE), color = "#6E8D9E"),
            br(),
            downloadButton(outputId = "dlPlot_ridgeGSEA", label = paste("Download Ridge Plot")),
            br(),
            selectInput(inputId = "ridgePlotColourByGSEA", label = "Colour plot by", choices = c("p.adjust", "enrichmentScore", "NES"), selected = "NES")
          ),
          tabPanel(
            title = "Upset plot",
            withSpinner(plotOutput(outputId = "upsetPlot_GSEA", inline = TRUE), color = "#6E8D9E"),
            br(),
            downloadButton(outputId = "dlPlot_upsetGSEA", label = paste("Download Upset Plot"))
          ),
          tabPanel(
            title = "Lollipop plot",
            withSpinner(plotOutput(outputId = "lollipopPlot_GSEA", inline = TRUE), color = "#6E8D9E"),
            downloadButton(outputId = "dlPlot_lollipopGSEA", label = paste("Download Lollipop Plot")),
            br(),
            sliderInput(inputId = "gseaLollipopScale", label = "Lollipop scale", min = 0.25, max = 5, value = 1, step = 0.25, width = "35%", ticks = FALSE)
          )
        )
      ),
      uiOutput("gseaIsoformAggUI"),
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
          sliderInput(inputId = "textSizeGSEA", label = "Figure Font Size", min = 4, max = 30, value = 10, step = 0.5, width = "85%"),
          sliderInput(inputId = "figHeightGSEA", label = "Figure Height", min = 100, max = 2000, value = 400, step = 10, width = "85%"),
          sliderInput(inputId = "figWidthGSEA", label = "Figure Width", min = 100, max = 2000, value = 600, step = 10, width = "85%"),
        ),
        column(
          width = 4,
          h4("Download settings"),
          selectInput(inputId = "downloadFormatGSEA", label = "Select format", choices = c("PDF", "SVG", "PNG"), selected = "PDF", width = "85%"),
          selectInput(inputId = "dpiGSEA", label = "PNG DPI", choices = c(72, 150, 300, 600, 1000), selected = 600, width = "85%"),
          textInput(inputId = "filename_GSEA", label = "Enter filename", value = "GSEA", width = "85%")
        )
      ),
      box(
        title = "Gene Filtering (Visualisation Only)",
        width = NULL,
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        status = "warning",
        helpText("IMPORTANT: These thresholds only affect visualisation (network plots and heatmaps). They do NOT change the GSEA results, which were pre-calculated."),
        column(
          width = 6,
          h4("Threshold Settings"),
          selectInput(inputId = "gseaVisPType", label = "P-Value Type:", choices = c("FDR", "Raw"), selected = "Raw", width = "85%"),
          numericInput(inputId = "gseaVisPThreshold", label = "P-Value Threshold:", value = 0.05, min = 0, max = 1, step = 0.01, width = "85%"),
          numericInput(inputId = "gseaVisLFC", label = "Log2 Fold Change Threshold:", value = 0, min = 0, max = 10, step = 0.25, width = "85%")
        ),
        column(
          width = 6,
          h4("Filter Mode"),
          radioButtons(
            inputId = "gseaVisFilterMode",
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
  ) # end of fluid row
) # end of tabItem: GSEA
