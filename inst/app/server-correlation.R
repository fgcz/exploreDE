output$correlationStatic <- renderPlot({
  NULL
})

# Update some inputs
updateSelectizeInput(session = session, inputId = "correlationGene1", label = "Feature X", choices = inputDataReactive()$genes, selected = "", server = TRUE)
updateSelectizeInput(session = session, inputId = "correlationGene2", label = "Feature Y", choices = inputDataReactive()$genes, selected = "", server = TRUE)
updateSelectInput(session = session, inputId = "correlationColourBy", label = "Select column to colour by", choices = inputDataReactive()$factors, selected = inputDataReactive()$groupingName)
updateSelectInput(session = session, inputId = "correlationShapeBy", label = "Select column to shape by", choices = c("None", inputDataReactive()$factors), selected = "None")
count_choices <- names(inputDataReactive()$countList)
default_count <- switch(inputDataReactive()$dataType,
  "RNASeq" = "Normalised + Log2",
  "proteomics" = if ("transformedData" %in% count_choices) "transformedData" else count_choices[1]
)
updateSelectInput(session = session, inputId = "correlationCounts", label = "Select count method to plot", choices = count_choices, selected = default_count)
# Automatically set batch correction if grouping2Name exists
batchDefault <- "None"
if (!is.null(inputDataReactive()$param$grouping2Name)) {
  if (inputDataReactive()$param$grouping2Name %in% inputDataReactive()$factors) {
    batchDefault <- inputDataReactive()$param$grouping2Name
  }
}
updateSelectInput(session = session, inputId = "correlationBatch", label = "Apply a batch-correction to the counts?", choices = c("None", inputDataReactive()$factors), selected = batchDefault)
updateSliderInput(session = session, inputId = "correlationCorPosX", label = "Label position x", min = 1, max = 100, value = 10, step = 1)
updateSliderInput(session = session, inputId = "correlationCorPosY", label = "Label position y", min = 1, max = 100, value = 10, step = 1)

observeEvent(
  {
    input$correlationColourBy
  },
  {
    output$correlationGroupBucket <- renderUI({
      req(input$tabs == "correlationsTab")
      bucket_list(
        header = "Drag and drop groups in order to be plotted (Ctrl/Cmd to select multiple)",
        group_name = "correlation_bucket_list",
        orientation = "horizontal",
        add_rank_list(
          text = "Include these groups",
          labels = as.list(unique(levels(as.factor(inputDataReactive()$dataset[[input$correlationColourBy]])))),
          input_id = "correlationKeepBucket",
          options = sortable::sortable_options(multiDrag = TRUE)
        ),
        add_rank_list(
          text = "Exclude these groups",
          labels = NULL,
          input_id = "correlationExcludeBucket",
          options = sortable::sortable_options(multiDrag = TRUE)
        )
      )
    })
    outputOptions(output, "correlationGroupBucket", suspendWhenHidden = FALSE)
  }
)

correlationPlot <- reactiveValues(plot = NULL)
observeEvent(
  {
    input$correlationMethod
    input$correlationGene1
    input$correlationGene2
    input$correlationColourBy
    input$correlationShapeBy
    input$correlationCounts
    input$correlationBatch
    input$correlationCorPosX
    input$correlationCorPosY
    input$correlationPointSize
    input$correlationPointStroke
    input$correlationPointAlpha
    input$textSizeCorrelation
    input$correlationKeepBucket
    input$correlationCountsLog
    input$correlationShowNames
    input$correlationMaxOverlaps
    input$correlationLabelPosX
    input$correlationLabelPosY
    # Add all color picker inputs as dependencies so plot updates when colors change
    lapply(seq_along(inputDataReactive()$factorLevels), function(i) {
      lapply(seq_along(inputDataReactive()$factorLevels[[i]]), function(j) {
        input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]], inputDataReactive()$factorLevels[[i]][[j]])]]
      })
    })
  },
  ignoreNULL = TRUE,
  ignoreInit = TRUE,
  {
    req(input$correlationGene1 != input$correlationGene2)
    req(nchar(input$correlationGene1) > 1 & nchar(input$correlationGene2) > 1)

    # prepare the metadat
    datasetCorrelation <- inputDataReactive()$dataset %>% data.frame(check.names = F)
    if (!"Name" %in% colnames(datasetCorrelation)) {
      datasetCorrelation <- datasetCorrelation %>% rownames_to_column("Name")
    }
    datasetCorrelation$Name <- factor(datasetCorrelation$Name, levels = datasetCorrelation$Name)
    datasetCorrelation <- datasetCorrelation[, c("Name", grepv(paste0("^", inputDataReactive()$factors, "$", collapse = "|"), colnames(datasetCorrelation)))]

    # get the counts
    countsCorrelation <- inputDataReactive()$countList[[input$correlationCounts]][c(input$correlationGene1, input$correlationGene2), ]

    # log2 if requested/necessary
    if (input$correlationCountsLog & !input$correlationCounts %in% c("Normalised + Log2", "VST", "transformedData")) {
      countsCorrelation <- log2(countsCorrelation + 1)
    }

    # remove batch effect if requested
    if (input$correlationBatch != "None") {
      countsCorrelation <- limma::removeBatchEffect(x = countsCorrelation, batch = inputDataReactive()$dataset[[input$correlationBatch]])
    }

    # prepare the counts
    countsCorrelation <- countsCorrelation %>%
      t() %>%
      data.frame(check.names = F) %>%
      rownames_to_column("Name") %>%
      left_join(datasetCorrelation, by = "Name")
    countsCorrelation <- countsCorrelation[which(countsCorrelation[[input$correlationColourBy]] %in% input$correlationKeepBucket), ]

    # get the colours
    coloursCorrelation <- setNames(lapply(input$correlationKeepBucket, function(k) {
      paste0(col2hex(input[[paste0("GroupColour", input$correlationColourBy, k)]]), "FF")
    }), input$correlationKeepBucket)

    # get the correct coefficient
    if (input$correlationMethod == "Spearman") {
      cor <- "rho"
    } else if (input$correlationMethod == "Kendall") {
      cor <- "tau"
    } else {
      cor <- "R"
    }

    # make the plot
    corr_plot1 <- ggplot(countsCorrelation, aes(x = .data[[input$correlationGene1]], y = .data[[input$correlationGene2]]))
    corr_plot1 <- corr_plot1 + geom_smooth(method = "lm", se = TRUE, colour = "black")
    corr_plot1 <- corr_plot1 + stat_cor(method = tolower(input$correlationMethod), cor.coef.name = cor, size = input$textSizeCorrelation / 3, show.legend = FALSE, label.x.npc = tolower(input$correlationCorPosX), label.y.npc = tolower(input$correlationCorPosY), label.args = list(fontface = "bold"))
    if (input$correlationShapeBy == "None") {
      corr_plot1 <- corr_plot1 + geom_point(aes(fill = .data[[input$correlationColourBy]]), pch = 21, size = input$correlationPointSize, stroke = input$correlationPointStroke, alpha = input$correlationPointAlpha)
    } else {
      corr_plot1 <- corr_plot1 + geom_point(aes(fill = .data[[input$correlationColourBy]], shape = .data[[input$correlationShapeBy]]), size = input$correlationPointSize, stroke = input$correlationPointStroke, alpha = input$correlationPointAlpha)
      corr_plot1 <- corr_plot1 + scale_shape_manual(values = c(rep(c(21, 22, 23, 24, 25), times = 10))[1:nlevels(as.factor(countsCorrelation[[input$correlationShapeBy]]))])
      corr_plot1 <- corr_plot1 + guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)), shape = guide_legend(override.aes = list(fill = "black", size = 3)))
    }
    if (input$correlationShowNames) {
      corr_plot1 <- corr_plot1 + geom_label_repel(
        aes(label = Name, fill = .data[[input$correlationColourBy]]),
        show.legend = FALSE,
        alpha = input$correlationPointAlpha,
        size = input$textSizeCorrelation / 3,
        nudge_x = input$correlationLabelPosX / 100,
        nudge_y = input$correlationLabelPosY / 100,
        segment.size = 0.5,
        segment.color = "black",
        max.overlaps = input$correlationMaxOverlaps
      )
    }
    corr_plot1 <- corr_plot1 + scale_fill_manual(breaks = names(coloursCorrelation), values = as.character(coloursCorrelation))
    corr_plot1 <- corr_plot1 + theme_prism(base_size = input$textSizeCorrelation)
    corr_plot1 <- corr_plot1 + theme(legend.text = element_text(size = input$textSizeCorrelation, color = "black", face = "bold"))
    corr_plot1
    correlationPlot$plot <- corr_plot1
    correlationPlot$counts <- countsCorrelation
  }
)

# view the plot
output$correlationStatic <- renderPlot(
  {
    req(!is.null(correlationPlot$plot))
    print(correlationPlot$plot)
  },
  width = function() {
    as.numeric(input$figWidthCorrelation)
  },
  height = function() {
    as.numeric(input$figHeightCorrelation)
  }
)
# download the plot
# download the plot
output$dlPlot_correlation <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = input$filename_correlation,
      design_name = inputDataReactive()$design,
      extension = tolower(input$downloadFormatCorrelation)
    )
  },
  content = function(file) {
    req(!is.null(correlationPlot$plot))
    export_plot_standard(
      plot_object = correlationPlot$plot,
      filename = file,
      format = tolower(input$downloadFormatCorrelation),
      width = as.numeric(input$figWidthCorrelation) / 75,
      height = as.numeric(input$figHeightCorrelation) / 75,
      dpi = as.numeric(input$dpiCorrelation)
    )
  }
)
# download the counts
# download the counts
output$dlTable_correlation <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_correlation, "_Counts"),
      design_name = inputDataReactive()$design,
      extension = "xlsx"
    )
  },
  content = function(file) {
    req(!is.null(correlationPlot$counts))
    export_excel_standard(
      df = correlationPlot$counts,
      filename = file,
      sheet_name = "Correlation Data",
      table_name = "CorrelationData"
    )
  }
)
