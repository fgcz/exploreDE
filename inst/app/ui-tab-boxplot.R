tabItem(
  tabName = "boxplotTab",
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
          "Here, multiple features can be selected to generate boxplots
        of expression across selected sample groups. features can either be
        selected from the input feature drop down, which is in order of 
        p-value and where specific features can be searched for by typing in the
        box, or typed/pasted in the text box."
        ),
        tags$p("Sample groups can be removed and reordered by dragging and dropping
         the group names in the bucket list. If available, a second factor can
         be selected and will be used to either change the shape of the points
         (boxplot + ggplot), or change the facet (barplot). The boxplot prism theme 
         doesn't seem to support different point shapes."),
        tags$p("Different transformation methods can be used to visualise the feature 
         counts. Normalised + Log2 is the log2 of library-size normalised 
         counts. FPKM, or fragments per kilobase per million, is the feature count
         following normalisation based on the mapping of fragments to features. 
         VST, or variance stabilising transformed, counts are a version of the 
         normalised + log2 but that has also been transformed to remove 
         heteroscedasticity. The normalised counts have not been transformed."
        ),
        tags$p("Turning the p-value stars on is still undergoing development and
         so may not work perfectly, and the plots might look a little different 
         when this is toggled off or on. Currently, it will only pull in the fdr p-values,
         because plotting multiple features and multiple differential expression tests
         means the adjusted p-value is desired."), 
        tags$p("It is possible to put all the plots on the same y-axis by checking
         the tick box. This makes it easy to visualise differences in feature 
         counts across multiple features."
        )
      ),
      box(
        title = "Settings",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tabsetPanel(
          tabPanel(
            title = "Main settings",
            selectizeInput(
              inputId = "boxplotGenes", label = "Features for Boxplots:", multiple = TRUE, choices = NULL, selected = NULL, 
              options = list(placeholder = 'Select features', plugins = list('remove_button', 'drag_drop', 'restore_on_backspace', 'clear_button'))
            ),
            helpText("Drag and drop to rearragne feature orders. You can also search for features by typing in the box above."),
            textAreaInput(inputId = "boxplotGenesText", label = "Or, paste list of features:", placeholder = "feature1 feature2 feature3 feature4 ", cols = 1, rows = 2, resize = "vertical"),
            uiOutput("boxplotGTError"),
            selectInput(inputId = "boxplotFactor1", label = "Select the main factor to plot by", choices = NULL, selected = NULL, multiple = FALSE),
            selectInput(inputId = "boxplotFactor2", label = "Select a second factor to plot by", choices = NULL, selected = NULL, multiple = FALSE),
            helpText("Second factor determines point shape for the boxplot"),
            selectInput(inputId = "boxplotCounts", label = "Select count method to plot", choices = NULL, selected = NULL, multiple = FALSE),
            checkboxInput(inputId = "boxplotZScore", label = "Transform to Z-Scores?", value = FALSE),
            helpText("Calculates per-gene Z-Scores. Works on any count method."),
            checkboxInput(inputId = "boxplotCountsLog", label = "Log2 Counts?", value = FALSE),
            helpText("Only applies to counts that are not already on in log space."),
            selectInput(inputId = "boxplotBatch", label = "Apply a batch-correction to the counts?", choices = NULL, multiple = FALSE),
            hr(style = "border-top: 1px solid #000000;"), h4("Feature Bucket"),
            h5(
              "Features will be added to this bucket as you select them from the DE table tab and from the inputs in various tabs. You can use this bucket to quickly 
              include/exclude these features from the plots as you need to by dragging and dropping them in order or into the exclude bucket."),
            uiOutput("geneBucket1"),
            actionButton(inputId = "resetGeneBucketBoxplot", label = "Empty the bucket?", icon = icon("bucket")),
            hr(style = "border-top: 1px solid #000000;"),
            uiOutput("bucket")
          ),
          tabPanel(
            title = "Figure settings", 
            tags$b("X Axis Label Settings"),
            checkboxInput(inputId = "boxplotShowConditions", label = "Show condition labels", value = FALSE),
            checkboxInput(inputId = "boxplotConditionFormat", label = "Change _./ etc. to line break", value = FALSE),
            sliderInput(inputId = "boxplotConditionAngle", label = "Condition label angle", min = 0, max = 90, value = 45, step = 45, width = "33%"),
            tags$b("Toggle to plot and/or dots, mean bar, violins, etc."),
            checkboxInput(inputId = "boxplotShowPoint", label = "Show points", value = TRUE),
            checkboxInput(inputId = "boxplotShowMeanBar", label = "Show mean bar", value = TRUE),
            checkboxInput(inputId = "boxplotMeanBarFront", label = "Bring mean bar to front", value = FALSE),
            checkboxInput(inputId = "boxplotShowBox", label = "Show boxes", value = FALSE),
            checkboxInput(inputId = "boxplotShowViolin", label = "Show violins", value = FALSE),
            helpText("If you select both box and violin, you'll get neither!"),
            splitLayout(
              sliderInput(inputId = "boxplotPointSize", label = "Point size", value = 4, min = 0.5, max = 20, step = 0.5, width = "80%"),
              sliderInput(inputId = "boxplotPointDodge", label = "Point dodge", min = 1, max = 5, value = 2, step = 0.25, width = "80%"),
              sliderInput(inputId = "boxplotPointAlpha", label = "Point alpha", min = 0.1, max = 1, value = 0.9, step = 0.1, width = "80%"),
            ),
            splitLayout(
              sliderInput(inputId = "boxplotPointBorder", label = "Point border", min = 0, max = 2, value = 0.5, step = 0.1, width = "80%"),
              sliderInput(inputId = "boxplotBoxAlpha", label = "Box/violin alpha", min = 0.1, max = 1, value = 0.2, step = 0.1, width = "80%"),
              sliderInput(inputId = "boxplotMeanLine", label = "Mean bar width", min = 0, max = 1, value = 0.5, step = 0.1, width = "80%")
            ),
            hr(style = "border-top: 1px solid #000000;"), h4("DE P-value"),
            tags$b("Show p-value stars on plot?"),
            helpText("The highlight tool doesn't work when you have p-values displayed."),
            checkboxInput(inputId = "boxplotShowP", label = "Show me those stars!", value = FALSE),
            splitLayout(
              # sliderInput(inputId = "boxplotShowPHeight", label = "P height", value = 1, min = 0, max = 10, step = 1, width = "80%"),
              sliderInput(inputId = "boxplotShowPLabelSize", label = "P size", min = 0, max = 10, value = 3, step = 0.5, width = "80%"),
              sliderInput(inputId = "boxplotShowPDodge", label = "P height", min = 0, max = 10, value = 1, step = 0.5, width = "80%")
              ),
            splitLayout(
              sliderInput(inputId = "boxplotShowPBracketSize", label = "Bracket size", min = 0, max = 10, value = 3, step = 0.5, width = "80%"),
              sliderInput(inputId = "boxplotShowPTipSizeA", label = "Tip length A", min = 0, max = 25, value = 0, step = 0.5, width = "80%"),
              sliderInput(inputId = "boxplotShowPTipSizeB", label = "Tip length B", min = 0, max = 25, value = 0, step = 0.5, width = "80%")
            )
          ),
          tabPanel(
            title = "Size Options",
            # selectInput(
            #   inputId = "boxplotThemeChoice",
            #   label = "Prefer another plot theme?", 
            #   choices = c("Prism", "Minimal", "BW", "Classic", "Void"), 
            #   selected = "Prism"),
            tags$b("Add dotted vertical lines between each group?"),
            checkboxInput(
              inputId = "boxplotVertLines", 
              label = "Add vertical lines", 
              value = FALSE
            ),
            tags$b("Make figure greyscale?"),
            checkboxInput(
              inputId = "boxplotGrey",
              label = "Make plots greyscale",
              value = FALSE),
            numericInput(inputId = "boxplotNCol", label = "Number of columns", min = 1, max = 10, value = 3, step = 1),
            splitLayout(
              sliderInput(inputId = "textSizeBoxplot", label = "Figure Font Size", min = 4, max = 30, value = 12, step = 0.5, width = "85%"),
              sliderInput(inputId = "figWidthBoxplot", label = "Figure Width", min = 100, max = 2000, value = 850, step = 10, width = "85%"),
              sliderInput(inputId = "figHeightBoxplot", label = "Figure Height", min = 100, max = 2000, value = 600, step = 10, width = "85%")
            )
          ),
          tabPanel(
            title = "Download Settings",
            # selectInput(inputId = "boxplotDownloadFormat", label = "Select format", choices = c("PDF", "SVG", "PNG"), selected = "PDF"),
            # selectInput(inputId = "boxplotDPI", label = "PNG DPI", choices = c(72, 150, 300, 600, 1000), selected = 600),
            selectInput(inputId = "boxplotFont", label = "Font", choices = c("serif", "sans", "mono"), selected = "sans")
            # textInput(inputId = "boxplotFilename", label = "Enter filename", value = "Feature_Boxplot"),
          )
        )
      )
    ),
    column(
      width = 8,
      box(
        title = "Plots",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        textOutput("boxplotDesign"),
        downloadButton(
          outputId = "dlBoxplotButtonCounts",
          label = paste("Download Boxplot Counts")
        ), br(),
        tags$p(
          "You can download a small Excel file of the counts of the features you have selected. You can use these counts to quickly make figures in, e.g., Graphpad."
        ), br(),
        tabsetPanel(
          id = "boxbarplots",
          tabPanel(
            title = "Boxplot",
            br(),
            downloadButton(
              outputId = "dlBoxplotButton",
              label = paste("Download Boxplot (PDF)")
            ),
            br(), br(),
            withSpinner(
              plotOutput(
                outputId = "boxplotStatic", 
                inline = TRUE, 
                brush = "boxplotBrush"
              ),
              color="#19a7cf"
            ),
            DT::dataTableOutput(outputId = "boxplotBrushTable")
          ),
          tabPanel(
            title = "Barplot",
            br(),
            checkboxInput(inputId = "showDotsBarplot", label = "Show dots?", value = TRUE),
            downloadButton(
              outputId = "dlBarplotButton",
              label = paste("Download Barplot (PDF)")
            ),
            br(),
            plotOutput(
              outputId = "barplotStatic", 
              inline = TRUE
            )
          )
        )
      )
    )
  )
)