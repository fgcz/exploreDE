tabItem(
  tabName = "goTileTab",
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
          "Select one or more GO Terms to plot, and/or fetch a KEGG pathway. All the
          features in each pathway will be displayed as a tile. Cells are coloured by
          log2 fold change, and have stars denoting p-value. You can change whether the
          raw or adjusted p-values are displayed, and optionally limit the features shown
          to those passing a log2 fold change and p-value threshold."
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        uiOutput(outputId = "contrastSelectorUI_goTile"),
        h4("Tile Input settings"),
        # GO term tiles — only when the features carry GO annotation columns.
        conditionalPanel(
          condition = "output.goTileGoSupported",
          selectizeInput(
            inputId = "goTileInput",
            choices = NULL,
            multiple = TRUE,
            label = "Enter the GO Term(s) you wish to plot",
            options = list(
              placeholder = 'Select features',
              plugins = list('remove_button', 'drag_drop', 'restore_on_backspace', 'clear_button')
            )
          )
        ),
        # KEGG pathway tile — RNASeq + supported species (hsa/mmu/rno) only.
        conditionalPanel(
          condition = "output.goTileKeggSupported",
          tags$b("Or plot a KEGG pathway"),
          textInput(
            inputId = "goTileKeggInput",
            label = "Enter a KEGG pathway ID",
            placeholder = "e.g. 04110 or mmu04110",
            value = ""
          ),
          actionButton(
            inputId = "goTileKeggGo",
            label = "Fetch pathway genes",
            icon = icon("search")
          ),
          uiOutput("goTileKeggInfo"),
          helpText("The fetched pathway is shown as an additional tile alongside any GO terms.")
        ),
        hr(style = "border-top: 1px solid #000000;"),
        selectInput(
          inputId = "pTypeGoTile",
          choices = c("FDR", "Raw"),
          label = "P-Value:",
          selected = "FDR"
        ),
        numericInput(
          inputId = "tileLimitColour",
          label = "Scale limit for colours (+/-)",
          value = 4,
          min = 1,
          max = 20,
          step = 0.5
        ),
        checkboxInput(
          inputId = "degOnlyGoTile",
          label = "Only show significant features?",
          value = FALSE
        ),
        # Significance limits, revealed only when the box above is ticked
        # (mirrors the heatmap KEGG significance panel).
        conditionalPanel(
          condition = "input.degOnlyGoTile == true",
          numericInput(
            inputId = "lfcGoTile",
            label = "Log2 Fold Change Threshold:",
            value = 0,
            min = 0,
            max = 10,
            step = 0.25
          ),
          selectInput(
            inputId = "pThresholdGoTile",
            choices = c(1.0, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001),
            label = "P-Value Threshold:",
            selected = 0.05
          )
        ),
        hr(style = "border-top: 1px solid #000000;"),
        h4("Pathway Tile Plot display settings"),
        # Select whether to cluster the rows (features):
        checkboxInput(
          inputId = "clusterRowsGoTile",
          label = "Cluster features?",
          value = TRUE
        ),
        tags$b("Turn off/on names of rows/columns"),
        # Display colnames?:
        checkboxInput(
          inputId = "colnamesGoTile",
          label = "Show sample names?",
          value = TRUE
        ),
        # Display gene names?:
        checkboxInput(
          inputId = "geneNamesGoTile",
          label = "Show gene names?",
          value = TRUE
        ),
        # GoTile colour picker
        colourpicker::colourInput(
          inputId = "GoTileColourRed",
          label = "Tile Plot Red (+)",
          value = "#801717"
        ),
        colourpicker::colourInput(
          inputId = "GoTileColourWhite",
          label = "Tile Plot White (0)",
          value = "#FFFFFF"
        ),
        colourpicker::colourInput(
          inputId = "GoTileColourBlue",
          label = "Tile Plot Blue (-)",
          value = "#113D69"
        ),
        numericInput(
          inputId = "textSizeGoTile",
          label = "Figure Font Size",
          min = 4,
          max = 30,
          value = 12,
          step = 0.5
        ),
        numericInput(
          inputId = "figWidthGoTile",
          label = "Figure Width",
          min = 100,
          max = 2000,
          value = 400,
          step = 10
        ),
        numericInput(
          inputId = "figHeightGoTile",
          label = "Figure Height",
          min = 100,
          max = 2000,
          value = 1600,
          step = 10
        ),
        legend_offset_ui("goTile")
      ) # end of box
    ), # end of settings column
    column(
      width = 9,
      box(
        title = "Pathway Tile Plot",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        textOutput("goTilePlotDesign"),
        br(),
        h4("Download Settings"),
        splitLayout(
          selectInput(
            inputId = "downloadFormatGoTile",
            label = "Select format",
            choices = c("PDF", "SVG", "PNG"),
            selected = "PDF"
          ),
          selectInput(inputId = "dpiGoTile", label = "PNG DPI", choices = c(72, 150, 300, 600, 1000), selected = 600)
        ),
        textInput(inputId = "filename_goTile", label = "Enter filename", value = "Pathway_Tile"),
        downloadButton(outputId = "dlGoTilePlot", label = "Download Pathway Tile Plot"),
        br(),
        plotOutput("goTileOutput", inline = TRUE)
      ) # end of box
    ) # end of wide column
  ) # end of fluid row
) # end of tabItem: heatmaps
