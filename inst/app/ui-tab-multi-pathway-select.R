tabItem(
  tabName = "multiPathwaySelectTab",
  fluidRow(
    column(
      width = 3,
      box(
        title = "Info",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        collapsible = TRUE,
        collapsed = TRUE,
        tags$p(
          "Hand-pick any GO term(s) and compare them across the loaded DE tests.
                    The full annotated GO term list is available, independent of significance.
                    The same selection drives both figures:"
        ),
        tags$ul(
          tags$li(
            tags$b("Lollipop plot: "),
            "x-axis is -log10 p-value, y-axis is the
                        pathway, and there is one lollipop per DE test, coloured by contrast."
          ),
          tags$li(
            tags$b("Network (cnetplot): "),
            "a gene-concept network of the selected
                        pathways. Each pathway node is a pie showing which contrasts it is enriched
                        in; gene nodes connect to the pathways that contain them."
          )
        ),
        tags$p(
          "Not sure which pathway to pick? Use ",
          tags$b("Find pathways by gene"),
          " to add
                    every tested pathway that contains a gene of interest."
        ),
        tags$p(
          tags$em(
            "Points/nodes are drawn from the ORA results, so a term only appears for
                    a contrast where it was tested. Non-significant terms have an adjusted p-value of
                    ~1; switch the P-value to Raw to spread them out."
          )
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        selectInput(
          inputId = "pwDirection",
          label = "Regulation direction",
          choices = c(
            "Up-Regulated" = "upGenes",
            "Down-Regulated" = "downGenes",
            "Both Directions" = "bothGenes"
          ),
          selected = "upGenes"
        ),
        selectInput(
          inputId = "pwPType",
          label = "P-value",
          choices = c("FDR", "Raw"),
          selected = "FDR"
        ),
        selectizeInput(
          inputId = "pathwaySelectInput",
          choices = NULL,
          multiple = TRUE,
          label = "Select the GO term(s) to plot",
          options = list(
            placeholder = "Select pathways",
            plugins = list("remove_button", "drag_drop", "restore_on_backspace", "clear_button")
          )
        ),
        hr(style = "border-top: 1px solid #000000;"),
        h4("Find pathways by gene"),
        selectizeInput(
          inputId = "pwGeneSearch",
          choices = NULL,
          multiple = TRUE,
          label = "Gene(s) of interest",
          options = list(
            placeholder = "e.g. TNF",
            plugins = list("remove_button", "restore_on_backspace", "clear_button")
          )
        ),
        actionButton(inputId = "pwGeneFindGo", label = "Add pathways containing these genes", icon = icon("search")),
        hr(style = "border-top: 1px solid #000000;"),
        h4("Figure Settings"),
        splitLayout(
          sliderInput(
            inputId = "textSizePwSelect",
            label = "Font Size",
            min = 4,
            max = 30,
            value = 12,
            step = 0.5,
            width = "85%",
            ticks = FALSE
          ),
          sliderInput(
            inputId = "figWidthPwSelect",
            label = "Width",
            min = 100,
            max = 2000,
            value = 800,
            step = 10,
            width = "85%",
            ticks = FALSE
          ),
          sliderInput(
            inputId = "figHeightPwSelect",
            label = "Height",
            min = 100,
            max = 2000,
            value = 800,
            step = 10,
            width = "85%",
            ticks = FALSE
          )
        ),
        hr(style = "border-top: 1px solid #000000;"),
        h4("Download Settings"),
        splitLayout(
          selectInput(
            inputId = "downloadFormatPwSelect",
            label = "Format",
            choices = c("PDF", "SVG", "PNG"),
            selected = "PDF"
          ),
          selectInput(inputId = "dpiPwSelect", label = "PNG DPI", choices = c(72, 150, 300, 600, 1000), selected = 600)
        )
      )
    ),
    column(
      width = 9,
      box(
        title = "Pathway selection plots",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        tabsetPanel(
          id = "pwSelectTabs",
          tabPanel(
            title = "Lollipop plot",
            br(),
            withSpinner(plotOutput(outputId = "pwLollipopPlot", inline = TRUE), color = "#00a65a"),
            br(),
            br(),
            sliderInput(
              inputId = "pwLollipopScale",
              label = "Lollipop scale",
              min = 0.25,
              max = 5,
              value = 1,
              step = 0.25,
              width = "50%",
              ticks = FALSE
            ),
            textInput(
              inputId = "filename_pwLollipop",
              label = NULL,
              value = "PathwaySelect_Lollipop",
              placeholder = "Filename",
              width = "30%"
            ),
            downloadButton(outputId = "dlPlot_pwLollipop", label = "Download This Plot")
          ),
          tabPanel(
            title = "Network (cnetplot)",
            br(),
            withSpinner(plotOutput(outputId = "pwCnetPlot", inline = TRUE), color = "#00a65a"),
            br(),
            br(),
            sliderInput(
              inputId = "pwCnetCatScale",
              label = "Pathway node scale",
              min = 0.25,
              max = 5,
              value = 1,
              step = 0.25,
              width = "35%",
              ticks = FALSE
            ),
            sliderInput(
              inputId = "pwCnetGeneScale",
              label = "Gene node scale",
              min = 0.25,
              max = 5,
              value = 1,
              step = 0.25,
              width = "35%",
              ticks = FALSE
            ),
            checkboxInput(inputId = "pwCnetGeneLabels", label = "Show gene labels?", value = TRUE),
            textInput(
              inputId = "filename_pwCnet",
              label = NULL,
              value = "PathwaySelect_Network",
              placeholder = "Filename",
              width = "35%"
            ),
            downloadButton(outputId = "dlPlot_pwCnet", label = "Download This Plot")
          ),
          tabPanel(
            title = "Data table",
            br(),
            downloadButton(outputId = "dlTable_pwBubble", label = "Download This Table (Excel)"),
            br(),
            br(),
            DT::dataTableOutput(outputId = "pwBubbleTable")
          )
        )
      )
    )
  )
)
