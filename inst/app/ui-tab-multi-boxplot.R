tabItem(
  tabName = "multiBoxplotTab",
  fluidRow(
    column(
      width = 4,
      box(
        title = "Info",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        collapsible = TRUE,
        collapsed = TRUE,
        tags$p(
          "Here, multiple features can be selected to generate boxplots
        of expression across selected sample groups. features can either be
        selected from the input feature drop down, which is in order of
        p-value and where specific features can be searched for by typing in the
        box, or typed/pasted in the text box."
        ),
        tags$p(
          "Sample groups can be removed and reordered by dragging and dropping
         the group names in the bucket list. If available, a second factor can
         be selected and will be used to either change the shape of the points
         (boxplot + ggplot), or change the facet (barplot). The boxplot prism theme
         doesn't seem to support different point shapes."
        ),
        tags$p(
          "Different transformation methods can be used to visualise the feature
         counts. Normalised + Log2 is the log2 of library-size normalised
         counts. FPKM, or fragments per kilobase per million, is the feature count
         following normalisation based on the mapping of fragments to features.
         VST, or variance stabilising transformed, counts are a version of the
         normalised + log2 but that has also been transformed to remove
         heteroscedasticity. The normalised counts have not been transformed."
        ),
        tags$p(
          "Turning the p-value stars on is still undergoing development and
         so may not work perfectly, and the plots might look a little different
         when this is toggled off or on. Currently, it will only pull in the fdr p-values,
         because plotting multiple features and multiple differential expression tests
         means the adjusted p-value is desired."
        ),
        tags$p(
          "It is possible to put all the plots on the same y-axis by checking
         the tick box. This makes it easy to visualise differences in feature
         counts across multiple features."
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        tabsetPanel(
          id = "multiBoxplotSettingsTabs",
          tabPanel(
            title = "Main settings",
            uiOutput("fb_search_mode_multiBoxplot"),
            selectizeInput(
              inputId = "multiBoxplotGenes",
              label = "Features for Boxplots:",
              multiple = TRUE,
              choices = NULL,
              selected = NULL,
              options = list(
                placeholder = "Select features",
                plugins = list("remove_button", "drag_drop", "restore_on_backspace", "clear_button")
              )
            ),
            helpText(
              "Drag and drop to rearrange feature orders. You can also search for features by typing in the box above."
            ),
            textAreaInput(
              inputId = "multiBoxplotGenesText",
              label = "Or, paste list of features:",
              placeholder = "feature1 feature2 feature3 feature4 ",
              cols = 1,
              rows = 2,
              resize = "vertical"
            ),
            uiOutput("multiBoxplotGTError"),
            selectInput(
              inputId = "multiBoxplotFactor1",
              label = "Select the main factor to plot by",
              choices = NULL,
              selected = NULL,
              multiple = FALSE
            ),
            selectInput(
              inputId = "multiBoxplotFactor2",
              label = "Select a second factor to plot by",
              choices = NULL,
              selected = NULL,
              multiple = FALSE
            ),
            helpText("Second factor determines point shape for the boxplot"),
            selectInput(
              inputId = "multiBoxplotCounts",
              label = "Select count method to plot",
              choices = NULL,
              selected = NULL,
              multiple = FALSE
            ),
            checkboxInput(inputId = "multiBoxplotZScore", label = "Transform to Z-Scores?", value = FALSE),
            helpText("Calculates per-gene Z-Scores. Works on any count method."),
            checkboxInput(inputId = "multiBoxplotCountsLog", label = "Log2 Counts?", value = FALSE),
            helpText("Only applies to counts that are not already on in log space."),
            selectInput(
              inputId = "multiBoxplotBatch",
              label = "Apply a batch-correction to the counts?",
              choices = NULL,
              multiple = FALSE
            ),
            selectInput(
              inputId = "multiBoxplotDesign",
              label = "Design factor to preserve during batch correction",
              choices = NULL,
              multiple = FALSE
            ),
            helpText(
              "The design factor tells removeBatchEffect which biological signal to protect. For two-factor DE this is locked to the primary grouping; for CountQC/Crispr you can pick any factor."
            ),
            hr(style = "border-top: 1px solid #000000;"),
            h4("Feature Bucket"),
            h5(
              "Features will be added to this bucket as you select them from the DE table tab and from the inputs in various tabs. You can use this bucket to quickly
              include/exclude these features from the plots as you need to by dragging and dropping them in order or into the exclude bucket."
            ),
            uiOutput("multiBoxplotGeneBucket1"),
            actionButton(inputId = "resetGeneBucketMultiBoxplot", label = "Empty the bucket?", icon = icon("bucket")),
            hr(style = "border-top: 1px solid #000000;"),
            uiOutput("multiBoxplotBucket")
          ),
          tabPanel(
            title = "Figure settings",
            h4("X Axis Settings"),
            checkboxInput(inputId = "multiBoxplotShowConditions", label = "Show condition labels", value = FALSE),
            checkboxInput(
              inputId = "multiBoxplotConditionFormat",
              label = "Change _./ etc. to line break",
              value = FALSE
            ),
            sliderInput(
              inputId = "multiBoxplotConditionAngle",
              label = "Condition label angle",
              min = 0,
              max = 90,
              value = 45,
              step = 45,
              width = "33%"
            ),
            hr(style = "border-top: 1px solid #000000;"),
            h4("Y Axis settings"),
            checkboxInput(inputId = "multiBoxplotYScaleFree", label = "Scale free y-axis?", value = FALSE),
            hr(style = "border-top: 1px solid #000000;"),
            h4("Point settings"),
            tags$b("Toggle to plot and/or dots, mean bar, violins, etc."),
            checkboxInput(inputId = "multiBoxplotPlotBorder", label = "Show border", value = TRUE),
            checkboxInput(inputId = "multiBoxplotShowPoint", label = "Show points", value = TRUE),
            checkboxInput(inputId = "multiBoxplotShowMeanBar", label = "Show mean bar", value = TRUE),
            checkboxInput(inputId = "multiBoxplotMeanBarFront", label = "Bring mean bar to front", value = FALSE),
            checkboxInput(inputId = "multiBoxplotShowBox", label = "Show boxes", value = FALSE),
            checkboxInput(inputId = "multiBoxplotShowViolin", label = "Show violins", value = FALSE),
            helpText("If you select both box and violin, you'll get neither!"),
            splitLayout(
              sliderInput(
                inputId = "multiBoxplotPointSize",
                label = "Point size",
                value = 4,
                min = 0.5,
                max = 20,
                step = 0.5,
                width = "80%",
                ticks = F
              ),
              sliderInput(
                inputId = "multiBoxplotPointDodge",
                label = "Point dodge",
                min = 1,
                max = 5,
                value = 2,
                step = 0.25,
                width = "80%",
                ticks = F
              ),
              sliderInput(
                inputId = "multiBoxplotPointAlpha",
                label = "Point alpha",
                min = 0.1,
                max = 1,
                value = 0.9,
                step = 0.1,
                width = "80%",
                ticks = F
              )
            ),
            splitLayout(
              sliderInput(
                inputId = "multiBoxplotPointBorder",
                label = "Point border",
                min = 0,
                max = 2,
                value = 0.5,
                step = 0.1,
                width = "80%",
                ticks = F
              ),
              sliderInput(
                inputId = "multiBoxplotBoxAlpha",
                label = "Box/violin alpha",
                min = 0.1,
                max = 1,
                value = 0.2,
                step = 0.1,
                width = "80%",
                ticks = F
              ),
              sliderInput(
                inputId = "multiBoxplotMeanLine",
                label = "Mean bar width",
                min = 0,
                max = 1,
                value = 0.5,
                step = 0.1,
                width = "80%",
                ticks = F
              )
            )
          ),
          tabPanel(
            title = "DE P-values",
            value = "dePValuesTab",
            tags$b("Show p-value stars on plot?"),
            checkboxInput(inputId = "multiBoxplotShowP", label = "Show me those stars!", value = TRUE),
            help("Adds adjusted p-values from the DE results to the plot for all available contrasts."),
            hr(),
            checkboxInput(inputId = "multiBoxplotColourP", label = "Colour P-values by Group?", value = FALSE),
            checkboxInput(inputId = "multiBoxplotPStack", label = "Stack P-values?", value = TRUE),
            sliderInput(
              inputId = "multiBoxplotPStackOffset",
              label = "Stack Offset",
              min = 0,
              max = 10,
              value = 0.7,
              step = 0.1,
              width = "66%",
              ticks = F
            ),
            hr(),
            splitLayout(
              sliderInput(
                inputId = "multiBoxplotShowPLabelSize",
                label = "P size",
                min = 0,
                max = 10,
                value = 3,
                step = 0.5,
                width = "80%",
                ticks = F
              ),
              sliderInput(
                inputId = "multiBoxplotShowPJust",
                label = "P just",
                value = 1,
                min = -1,
                max = 1,
                step = 0.1,
                width = "80%",
                ticks = F
              ),
              sliderInput(
                inputId = "multiBoxplotShowPDodge",
                label = "P height",
                min = 0,
                max = 10,
                value = 1,
                step = 0.5,
                width = "80%",
                ticks = F
              )
            ),
            splitLayout(
              sliderInput(
                inputId = "multiBoxplotShowPBracketSize",
                label = "Bracket size",
                min = 0,
                max = 10,
                value = 3,
                step = 0.5,
                width = "80%",
                ticks = F
              ),
              sliderInput(
                inputId = "multiBoxplotShowPTipSizeA",
                label = "Tip length A",
                min = 0,
                max = 25,
                value = 0,
                step = 0.5,
                width = "80%",
                ticks = F
              ),
              sliderInput(
                inputId = "multiBoxplotShowPTipSizeB",
                label = "Tip length B",
                min = 0,
                max = 25,
                value = 0,
                step = 0.5,
                width = "80%",
                ticks = F
              ),
              sliderInput(
                inputId = "multiBoxplotShowPBracketShorten",
                label = "Bracket shorten",
                min = 0,
                max = 1,
                value = 0,
                step = 0.05,
                width = "80%",
                ticks = F
              )
            )
          ),
          tabPanel(
            title = "Size & Colours",
            tags$b("Add dotted vertical lines between each group?"),
            checkboxInput(
              inputId = "multiBoxplotVertLines",
              label = "Add vertical lines",
              value = FALSE
            ),
            tags$b("Make figure greyscale?"),
            checkboxInput(
              inputId = "multiBoxplotGrey",
              label = "Make plots greyscale",
              value = FALSE
            ),
            numericInput(
              inputId = "multiBoxplotNCol",
              label = "Number of columns",
              min = 1,
              max = 10,
              value = 3,
              step = 1,
              width = "70%"
            ),
            splitLayout(
              sliderInput(
                inputId = "textSizeMultiBoxplot",
                label = "Figure Font Size",
                min = 4,
                max = 30,
                value = 12,
                step = 0.5,
                width = "85%",
                ticks = F
              ),
              sliderInput(
                inputId = "figWidthMultiBoxplot",
                label = "Figure Width",
                min = 100,
                max = 2000,
                value = 850,
                step = 10,
                width = "85%",
                ticks = F
              ),
              sliderInput(
                inputId = "figHeightMultiBoxplot",
                label = "Figure Height",
                min = 100,
                max = 2000,
                value = 600,
                step = 10,
                width = "85%",
                ticks = F
              )
            )
          ),
          tabPanel(
            title = "Download",
            downloadButton(outputId = "dlTable_multiBoxplotCounts", label = paste("Download Boxplot Counts")),
            helpText("You can download a small Excel file of the counts of the features you have selected."),
            checkboxInput(inputId = "includeDEStatsMultiBoxplot", label = "Include DE statistics?", value = FALSE),
            helpText("Adds log2FC, p-value, and FDR columns from all available contrasts."),
            hr(),
            h4("Plot Download Settings"),
            selectInput(
              inputId = "downloadFormatMultiBoxplot",
              label = "Select format",
              choices = c("PDF", "SVG", "PNG"),
              selected = "PDF"
            ),
            selectInput(
              inputId = "dpiMultiBoxplot",
              label = "PNG DPI",
              choices = c(72, 150, 300, 600, 1000),
              selected = 600
            ),
            textInput(inputId = "filename_multiBoxplot", label = "Enter filename", value = "MultiBoxplot"),
            br(),
            downloadButton(outputId = "dlPlot_multiBoxplot", label = paste("Download Boxplot")),
            downloadButton(outputId = "dlPlot_multiBarplot", label = paste("Download Barplot")),
            br(),
            br(),
            selectInput(
              inputId = "multiBoxplotFont",
              label = "Font",
              choices = c("serif", "sans", "mono"),
              selected = "sans"
            )
          )
        )
      )
    ),
    column(
      width = 8,
      box(
        title = "Info",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        textOutput("multiBoxplotDesign"),
        # uiOutput("multiBoxplotWarn", inline = T),
        tabsetPanel(
          id = "multiBoxbarplots",
          tabPanel(
            title = "Boxplot",
            withSpinner(
              plotOutput(outputId = "multiBoxplotStatic", inline = TRUE, brush = "multiBoxplotBrush"),
              color = "#00a65a"
            ),
            DT::dataTableOutput(outputId = "multiBoxplotTable")
          ),
          tabPanel(
            title = "Barplot",
            plotOutput(outputId = "multiBarplotStatic", inline = TRUE),
            br(),
            checkboxInput(inputId = "showDotsMultiBarplot", label = "Show dots?", value = TRUE),
            radioButtons(inputId = "multiBarplotSDorSE", label = "Error bars?", choices = c("None", "SE", "SD"))
          )
        )
      )
    )
  )
)
