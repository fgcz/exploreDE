tabItem(
  tabName = "keggTab",
  fluidRow(
    column(
      width = 3,
      box(
        title = "Info",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        collapsible = TRUE,
        collapsed = TRUE,
        tags$p(
          "Enter a KEGG pathway (such as hsa05200, mmu05200, etc.). The pathway 
          will be downloaded and gene boxes will be coloured if they are 
          differentially expressed. You can select the thresholds for defining a
          DE in the settings box. If a gene box is white, it's not DE. If a gene
          box is grey, it is DE, but the log2 fold change is close to zero. If 
          the gene box is red, the gene has a positive log2 fold change, and blue means 
          negative log2 fold change."
        ),
        tags$p(
          "The best way to view/save these images is to right click and save-as or
          open in a new tab. It's not possible to save them as PDFs, sorry about that!"
        ),
        tags$b("Please note:"),
        tags$p(
          "In any KEGG pathway, a single KEGG ortholog may actually represent 
           multiple gene/protein isoforms. In the case where multiple 
           genes/isoforms in your results map to the same KEGG ortholog, only the Ensembl gene 
           ID with the greatest log fold change will be mapped. If
           the different isoforms actually have diverging log fold changes, this 
           will not be represented on the KEGG pathway. The best way to identify if 
           your ortholog of interest has multiple genes/isoforms mapped to it is to 
           open the pathway on KEGG's website and hover over/click the ortholog of interest."
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        h4("KEGG Pathway plot settings"),
        textInput(
          inputId = "keggInput",
          placeholder = "hsa05200",
          value = NULL,
          label = "Enter the KEGG Pathway you wish to plot"
        ),
        numericInput(
          inputId = "lfcKEGG",
          label = "Log2 Fold Change Threshold:",
          value = 0.5,
          min = 0,
          max = 10,
          step = 0.25
        ),
        selectInput(
          inputId = "pTypeKEGG",
          choices = c("FDR", "Raw"),
          label = "P-Value:",
          selected = "FDR"
        ),
        selectInput(
          inputId = "pThresholdKEGG",
          choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001),
          label = "P-Value Threshold:",
          selected = 0.05
        ),
        numericInput(
          inputId = "keggLimitColour",
          label = "Scale limit for colours (+/-)",
          value = 1, min = 1, max = 20, step = 0.5
        )
      ) # end of box
    ), # end of settings column
    column(
      width = 9,
      box(
        title = "KEGG Pathway Plot",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        textOutput("keggPlotDesign"),
        br(),
        plotOutput("keggOutput", inline = TRUE)
      )
    )
  )
)