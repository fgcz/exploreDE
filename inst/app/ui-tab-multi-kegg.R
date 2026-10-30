tabItem(
  tabName = "multiKeggTab",
  fluidRow(
    column(
      width = 2,
      box(
        title = "Info",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        collapsible = TRUE,
        collapsed = TRUE,
        tags$p(
          "If your data is from human or mouse, you can plot specific KEGG pathways. Enter a KEGG pathway ID (e.g., hsa05200, or mmu05200), and that pathway will be downloaded.
          The log2 fold changes of significantly differentially expressed genes of your loaded contrasts will be overlayed onto each gene in the pathway, based on the settings selected.
          Each gene will be split into as many sections as you have contrasts, and each section will be coloured according to the log2FC of its respective contrast. The order is as you have in the table in the 'Import Contrasts' tab.
          Please note: KEGG orthologs can represent multiple different Ensembl/Entrez gene IDs. In these cases, the greatest absolute log2FC will be displayed automatically. Please consult the interactive pathway map on the KEGG
          website to see which orthologs are represented in a given gene box, if you see any strange inconsistencies.
          "
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        textInput(
          inputId = "KEGGID",
          label = "Input KEGG Pathway ID",
          placeholder = "e.g., 05200, hsa05200, mmu05200"
        ),
        numericInput(
          inputId = "lfcMultiKEGG",
          label = "Log2 Fold Change Threshold:",
          value = 0.5,
          min = 0,
          max = 10,
          step = 0.25
        ),
        radioButtons(
          inputId = "pTypeMultiKEGG",
          choices = c("FDR", "Raw"),
          label = "P-Value:",
          selected = "FDR"
        ),
        radioButtons(
          inputId = "pThresholdMultiKEGG",
          choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001),
          label = "P-Value Threshold:",
          selected = 0.05
        ),
        numericInput(
          inputId = "KEGGLimitColour",
          label = "Scale limit for colours (+/-)",
          value = 2,
          min = 1,
          max = 20,
          step = 0.5
        ),
      )
    ),
    column(
      width = 10,
      box(
        title = "KEGG Pathway Plot",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        h5("Gene box contrast order"),
        helpText(
          "Each gene box is divided to represent the different contrats loaded. From left to right, the order is:"
        ),
        tableOutput(outputId = "KEGGOrder"),
        hr(),
        textInput(inputId = "filename_multiKegg", label = NULL, value = "Multi_KEGG_Plot", placeholder = "Filename"),
        downloadButton("dlPlot_multiKegg", "Download KEGG Plot"),
        br(),
        br(),
        withSpinner(plotOutput(outputId = "KEGGout", inline = TRUE), color = "#00a65a")
      )
    )
  )
)
