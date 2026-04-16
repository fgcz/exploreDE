# Create reactive for design to use across all outputs and download handlers
designReactive <- reactive({
  req(inputDataReactive())
  if (inputDataReactive()$dataType == "RNASeq") {
    return(inputDataReactive()$design)
  } else if (inputDataReactive()$dataType == "proteomics") {
    if (isTruthy(input$contrastSelected)) {
      return(input$contrastSelected)
    } else {
      return("No contrast selected")
    }
  } else if (inputDataReactive()$dataType %in% c("CountQC", "CrisprCountQC")) {
    return(inputDataReactive()$design)
  }
  return("Unknown")
})

output$correlationStatic <- renderPlot({
  NULL
})

# Update some inputs
updateSelectizeInput(session = session, inputId = "correlationGene1", label = "Feature X", choices = inputDataReactive()$genes, selected = "", server = TRUE)
updateSelectizeInput(session = session, inputId = "correlationGene2", label = "Feature Y", choices = inputDataReactive()$genes, selected = "", server = TRUE)

# Update feature choices when collapse guides toggle changes
observeEvent(input$correlationCollapseGuides, {
  if (inputDataReactive()$dataType == "CrisprCountQC") {
    if (isTRUE(input$correlationCollapseGuides)) {
      gene_choices <- sort(unique(inputDataReactive()$sgRNA_to_gene))
    } else {
      gene_choices <- inputDataReactive()$genes
    }
    updateSelectizeInput(session, "correlationGene1", choices = gene_choices, selected = "", server = TRUE)
    updateSelectizeInput(session, "correlationGene2", choices = gene_choices, selected = "", server = TRUE)
    correlationPlot$plot <- NULL
  }
}, ignoreInit = TRUE)
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

# Update feature choices when search mode changes
observeEvent(input$fb_search_by_correlation, ignoreInit = TRUE, {
  mode <- input$fb_search_by_correlation
  choices <- if (mode == "gene_name" && !is.null(inputDataReactive()$gene_choices)) {
    inputDataReactive()$gene_choices
  } else if (mode == "protein_Id" && !is.null(inputDataReactive()$protein_choices)) {
    inputDataReactive()$protein_choices
  } else {
    inputDataReactive()$genes
  }
  updateSelectizeInput(session, "correlationGene1", choices = choices, selected = "", server = TRUE)
  updateSelectizeInput(session, "correlationGene2", choices = choices, selected = "", server = TRUE)
  correlationPlot$plot <- NULL
})

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
    countsCorrelation <- inputDataReactive()$countList[[input$correlationCounts]]
    # Collapse sgRNA guides to gene level if requested
    if (inputDataReactive()$dataType == "CrisprCountQC" && isTRUE(input$correlationCollapseGuides)) {
      countsCorrelation <- collapse_guides_to_genes(countsCorrelation, inputDataReactive()$sgRNA_to_gene)
    }

    # Resolve gene/protein names to feature_ids if in alternative search mode
    corr_mode <- input$fb_search_by_correlation %||% "feature_id"
    corr_features <- list(input$correlationGene1, input$correlationGene2)
    corr_labels <- c(input$correlationGene1, input$correlationGene2)
    if (corr_mode != "feature_id") {
      corr_features <- lapply(corr_labels, function(sel) {
        resolve_features(sel, corr_mode)
      })
    }
    # Average across feature_ids when a gene/protein maps to multiple rows
    # Build axis labels that show "(mean of N features)" when averaging
    corr_axis_labels <- corr_labels
    corr_rows <- lapply(seq_along(corr_features), function(i) {
      fids <- corr_features[[i]]
      fids <- fids[fids %in% rownames(countsCorrelation)]
      req(length(fids) > 0)
      if (length(fids) == 1) {
        countsCorrelation[fids, , drop = FALSE]
      } else {
        corr_axis_labels[i] <<- paste0(corr_labels[i], " (mean of ", length(fids), " features)")
        showNotification(
          paste0(corr_labels[i], ": averaging ", length(fids), " features for correlation"),
          type = "warning", duration = 6
        )
        matrix(
          colMeans(countsCorrelation[fids, , drop = FALSE], na.rm = TRUE),
          nrow = 1, dimnames = list(NULL, colnames(countsCorrelation))
        )
      }
    })
    countsCorrelation <- do.call(rbind, corr_rows)
    rownames(countsCorrelation) <- corr_labels

    # Warn about sparsity — count non-NaN samples per feature
    n_samples <- ncol(countsCorrelation)
    for (i in seq_along(corr_labels)) {
      n_valid <- sum(!is.na(countsCorrelation[i, ]) & !is.nan(as.numeric(countsCorrelation[i, ])))
      if (n_valid < n_samples) {
        showNotification(
          paste0(corr_labels[i], ": data available for ", n_valid, " of ", n_samples, " samples"),
          type = if (n_valid < n_samples / 2) "error" else "warning",
          duration = 8
        )
      }
    }

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
    corr_plot1 <- corr_plot1 + labs(x = corr_axis_labels[1], y = corr_axis_labels[2])
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
      design_name = designReactive(),
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
      design_name = designReactive(),
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
