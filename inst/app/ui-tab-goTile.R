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
          "Select a GO Term to plot. All the features in that pathway will be displayed.
          Cells are coloured by log2 fold change, and have stars denoting p-value. 
          You can change whether the raw or adjusted p-values are displayed. You can
          also only display features with a p-value of <= 0.05."
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        h4("GO Tile Input settings"),
        selectizeInput(
          inputId = "goTileInput",
          choices = c(),
          multiple = TRUE, 
          label = "Enter the GO Term you wish to plot"
        ),
        selectInput(
          inputId = "pTypeGoTile",
          choices = c("FDR", "Raw"),
          label = "P-Value:",
          selected = "FDR"
        ),
        numericInput(
          inputId = "tileLimitColour",
          label = "Scale limit for colours (+/-)",
          value = 4, min = 1, max = 20, step = 0.5
        ),
        checkboxInput(
          inputId = "degOnlyGoTile",
          label = "Only show features with significant p-values?",
          value = FALSE
        ),
        hr(style = "border-top: 1px solid #000000;"), 
        h4("GO Tile Plot display settings"),
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
          inputId = "GoTileColourRed", label = "Tile Plot Red (+)",
          value = "#801717"
        ),
        colourpicker::colourInput(
          inputId = "GoTileColourWhite", label = "Tile Plot White (0)",
          value = "#FFFFFF"
        ),
        colourpicker::colourInput(
          inputId = "GoTileColourBlue", label = "Tile Plot Blue (-)",
          value = "#113D69"
        ),
        numericInput(
          inputId = "textSizeGoTile",
          label = "Figure Font Size", min = 4, max = 30,
          value = 12, step = 0.5
        ),
        numericInput(
          inputId = "figWidthGoTile",
          label = "Figure Width", min = 100, max = 2000,
          value = 400, step = 10
        ),
        numericInput(
          inputId = "figHeightGoTile",
          label = "Figure Height", min = 100, max = 2000,
          value = 1600, step = 10
        )
      ) # end of box
    ), # end of settings column
    column(
      width = 9,
      box(
        title = "Go Tile Plot",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        textOutput("goTilePlotDesign"), br(),
        downloadButton(outputId = "dlGoTilePlot", label = "Download GO Tile Plot (PDF)"),
        br(),
        plotOutput("goTileOutput", inline = TRUE)
      ) # end of box
    ) # end of wide column
  ) # end of fluid row
) # end of tabItem: heatmaps