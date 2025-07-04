req(!is.null(inputDataReactive()$dataType))

if (inputDataReactive()$dataType == "proteomics") {
  updateNumericInput(session = session, inputId = "xLimVolcano", value = 2)
}

volcanoColoursDefault <- c(
  "Highlight" = "goldenrod3",
  "NotSignificant" = "grey60",
  "SignificantUp" = "lightcoral",
  "SignificantDown" = "steelblue",
  "FoldChangeUp" = "lightpink2",
  "FoldChangeDown" = "skyblue",
  "FoldChangeSignificantUp" = "firebrick3",
  "FoldChangeSignificantDown" = "deepskyblue4"
)

if (inputDataReactive()$dataType == "proteomics") {
  req(!is.null(inputDataReactive()$seqAnnoList))
  if("modelName" %in% colnames(inputDataReactive()$seqAnnoList[[1]])) {
    volcanoColoursDefault <- c(volcanoColoursDefault, "Imputed" = "burlywood3")
  }
}

output$volcanoColourPicker <- renderUI({
  lapply(names(volcanoColoursDefault), function(x) {
    colourpicker::colourInput(
      inputId = paste0("volcanoColour", x),
      label = x,
      value = as.character(volcanoColoursDefault[x]),
      closeOnClick = TRUE,
      returnName = TRUE
    )
  })
})
outputOptions(output, "volcanoColourPicker", suspendWhenHidden = FALSE)

updateSelectizeInput(
  session = session,
  inputId = "volcanoGenes",
  choices = inputDataReactive()$genes,
  selected = "",
  server = TRUE
)

observeEvent({
  genesReactive
}, {
  output$geneBucket2 <- renderUI({
    bucket_list(
      header = "Drag and drop features in order to be plotted",
      group_name = "bucket_list_group",
      orientation = "horizontal",
      add_rank_list(
        text = "Include these features in this order",
        labels = genesReactive$genes,
        input_id = "volcanoKeepBucketGenes"),
      add_rank_list(
        text = "Exclude these features",
        labels = NULL,
        input_id = "boxExcludeBucketGenes")
    )
  })
  outputOptions(output, "geneBucket2", suspendWhenHidden = FALSE)
})

observeEvent(input$resetGeneBucketVolcano, {
  genesReactive$genes <- NULL
  updateSelectizeInput(session, inputId = "volcanoGenes", choices = inputDataReactive()$genes, selected = NULL, server = T)
  updateTextAreaInput(session, inputId = "volcanoGenesText", value = NULL)
}, ignoreNULL = T, ignoreInit = T)

observeEvent(
  {
    input$pTypeVolcano
    input$contrastSelected
    inputDataReactive()
  }, ignoreInit = FALSE, ignoreNULL = TRUE,
  {
    if (inputDataReactive()$dataType == "RNASeq") {
      seqAnno <- inputDataReactive()$seqAnno
      design <- inputDataReactive()$design
      output$volcanoDesign <- renderText({
        design
      })
      seqAnnoFilt <- seqAnno[which(seqAnno$usedInTest),]
      updateTextInput(session = session, inputId = "filnameVolcano", value = design)
    }

    if (inputDataReactive()$dataType == "proteomics") {
      req(!is.null(input$contrastSelected))
      seqAnnoFilt <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]
      seqAnnoFilt <- seqAnnoFilt %>% dplyr::select(gene_name, log2Ratio, pValue, fdr)

      output$volcanoDesign <- renderText({
        input$contrastSelected
      })
      updateTextInput(session = session, inputId = "filnameVolcano", value = input$contrastSelected)
    }

    if (input$pTypeVolcano == "Raw") {
      pTypeVolcano <- "pValue"
    } else {
      pTypeVolcano <- "fdr"
    }

    updateNumericInput(
      session = session,
      inputId = "yLimVolcano",
      value = ceiling(max(-log10(seqAnnoFilt[[pTypeVolcano]]))*1.1)
    )
    updateNumericInput(
      session = session,
      inputId = "xLimVolcano",
      value = ceiling(max(abs(seqAnnoFilt[["log2Ratio"]]))*1.1)
    )
  }
)

if (inputDataReactive()$dataType == "proteomics") {
  req(!is.null(inputDataReactive()$seqAnnoList))

  if ("nrPeptides" %in% colnames(inputDataReactive()$seqAnnoList[[1]])) {
    output$nrPeptidesVolcanoUI <- renderUI({
      sliderInput(inputId = "nrPeptideVolcano", label = "Minimum number of peptides", min = 0, max = 20, value = 0, step = 1, width = "85%")
    })
    outputOptions(output, "nrPeptidesVolcanoUI", suspendWhenHidden = FALSE)
  } else {
    output$nrPeptidesVolcanoUI <- renderUI({ NULL })
  }

  if("modelName" %in% colnames(inputDataReactive()$seqAnnoList[[1]])) {
    output$showImputedVolcanoUI <- renderUI({
      checkboxInput(inputId = "showImputedVolcano", label = "Include imputed features?", value = TRUE)
    })
    output$highlightImputedVolcanoUI <- renderUI({
      checkboxInput(inputId = "highlightImputedVolcano", label = "Highlight imputed features?", value = FALSE)
    })
    outputOptions(output, "showImputedVolcanoUI", suspendWhenHidden = FALSE)
    outputOptions(output, "highlightImputedVolcanoUI", suspendWhenHidden = FALSE)
  } else {
    output$showImputedVolcanoUI <- renderUI({NULL})
    output$highlightImputedVolcanoUI <- renderUI({NULL})
  }
}

volcanoGenesReactive <- reactiveValues(Results = NULL)

observeEvent({
  input$volcanoBrush
  input$MABrush
  }, ignoreNULL = FALSE, ignoreInit = TRUE, {
    volcanoBrushGenes <- brushedPoints(volcanoResultsList()$volcanoTable, input$volcanoBrush)
    volcanoBrushGenes <- volcanoBrushGenes$gene_name
    maBrushGenes <- brushedPoints(volcanoResultsList()$volcanoTable, input$MABrush)
    maBrushGenes <- maBrushGenes$gene_name
    volcanoGenesReactive$volcanoBrushGenes = unique(c(volcanoBrushGenes, maBrushGenes))
})

observeEvent(input$lfcVolcano, ignoreInit = T, {
  showNotification("Volcano log2FC threshold changed")
})

# Volcano Plot Output ----
volcanoResultsList <- eventReactive(
  {
    input$lfcVolcano
    input$yLimVolcano
    input$xLimVolcano
    input$pTypeVolcano
    input$pThresholdVolcano
    input$volcanoShowGenes
    input$volcanoKeepBucketGenes
    input$contrastSelected
    input$nrPeptideVolcano
    input$showImputedVolcano
    input$highlightImputedVolcano
    lapply(names(volcanoColoursDefault), function(x) {
      input[[paste0("volcanoColour", x)]]
    })
  },
  ignoreInit = TRUE, ignoreNULL = FALSE,
  {
    if (inputDataReactive()$dataType == "RNASeq") {
      seqAnno <- inputDataReactive()$seqAnno
      design <- inputDataReactive()$design
      output$volcanoDesign <- renderText({
        design
      })
      seqAnnoFilt <- seqAnno[which(seqAnno$usedInTest),]
    }
    if (inputDataReactive()$dataType == "proteomics") {
      seqAnnoFilt <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]
      if (!is.null(input$nrPeptideVolcano)) {
        if (any(seqAnnoFilt$nrPeptides >= input$nrPeptideVolcano)) {
          seqAnnoFilt <- seqAnnoFilt[which(seqAnnoFilt$nrPeptides >= input$nrPeptideVolcano),]
        } else {
          shinyalert::shinyalert(title = "Oops!", text = "No features with this number of peptides!", type = "error", closeOnClickOutside = TRUE, showCancelButton = FALSE, showConfirmButton = TRUE, timer = 5000)
        }
      }
      if (!is.null(input$showImputedVolcano)) {
        if (!input$showImputedVolcano) {
          seqAnnoFilt <- seqAnnoFilt[grep("^Linear", seqAnnoFilt$modelName),]
        }
      }
      if (!is.null(input$highlightImputedVolcano)) {
        if (input$highlightImputedVolcano) {
          seqAnnoFilt <- seqAnnoFilt %>% dplyr::select(gene_name, description, log2Ratio, pValue, fdr, modelName)
        }
      } else {
        seqAnnoFilt <- seqAnnoFilt %>% dplyr::select(gene_name, description, log2Ratio, pValue, fdr)
      }

      output$volcanoDesign <- renderText({
        input$contrastSelected
      })
    }

    if (input$pTypeVolcano == "Raw") {
      pTypeVolcano <- "pValue"
    } else {
      pTypeVolcano <- "fdr"
    }

    # Get the colours
    volcanoColours <- unlist(setNames(lapply(names(volcanoColoursDefault), function(x) {
      input[[paste0("volcanoColour", x)]]
    }), names(volcanoColoursDefault)))

    # Create table of gene names, LFC, and p-value selected by user:
    volcanoTable <- seqAnnoFilt[, c("gene_name", "description", "log2Ratio", pTypeVolcano)]

    # Get the -log10 of the p-value:
    volcanoTable["log10p"] <- -log10(as.numeric(volcanoTable[[pTypeVolcano]]))

    # Get the -log10p of p-value threshold as a single value:
    log10pToFilter <- -log10(as.numeric(input$pThresholdVolcano))

    # Get mean log2 count for MA plot
    if (inputDataReactive()$dataType == "RNASeq") {
      rMs <- rowMeans(inputDataReactive()$countList$`Normalised + Log2`)
    } else {
      rMs <- rowMeans(inputDataReactive()$countList[[2]])
      names(rMs) <- names(rMs) %>% gsub("\\~.*", "", .)
    }
    volcanoTable$Log2_Mean <- rMs[match(volcanoTable$gene_name, names(rMs))]

    # Determine what significance Status each gene is and label accordingly
    ## Not significant:
    volcanoTable["Status"] <- "NotSignificant"
    ## Significant up only:
    volcanoTable[which(volcanoTable["log10p"] >= log10pToFilter & volcanoTable["log2Ratio"] <= as.numeric(input$lfcVolcano) & volcanoTable["log2Ratio"] > 0), "Status"] <- "SignificantUp"
    ## Significant down only:
    volcanoTable[which(volcanoTable["log10p"] >= log10pToFilter & volcanoTable["log2Ratio"] < 0 & volcanoTable["log2Ratio"] >= -as.numeric(input$lfcVolcano)), "Status"] <- "SignificantDown"
    ## LFC up only:
    volcanoTable[which(volcanoTable["log10p"] < log10pToFilter & volcanoTable["log2Ratio"] >= as.numeric(input$lfcVolcano)), "Status"] <- "FoldChangeUp"
    ## LFC down only:
    volcanoTable[which(volcanoTable["log10p"] < log10pToFilter & volcanoTable["log2Ratio"] <= -as.numeric(input$lfcVolcano)), "Status"] <- "FoldChangeDown"
    # LFC & p-value up:
    volcanoTable[which(volcanoTable["log10p"] >= log10pToFilter & volcanoTable["log2Ratio"] >= as.numeric(input$lfcVolcano)), "Status"] <- "FoldChangeSignificantUp"
    # LFC & p-value down:
    volcanoTable[which(volcanoTable["log10p"] >= log10pToFilter & volcanoTable["log2Ratio"] <= -as.numeric(input$lfcVolcano)), "Status"] <- "FoldChangeSignificantDown"
    
    if (!is.null(input$highlightImputedVolcano)) {
      if (input$highlightImputedVolcano) {
        volcanoTable[which(volcanoTable$gene_name %in% seqAnnoFilt$gene_name[grep("^Imputed", seqAnnoFilt$modelName)]), "Status"] <- "Imputed"
      }
    }

    # Set limits of data frame to x/y-axes limits
    volcanoTable$log10p_full <- volcanoTable$log10p
    volcanoTable$log2Ratio_full <- volcanoTable$log2Ratio
    volcanoTable$log10p[volcanoTable$log10p > as.numeric(input$yLimVolcano)] <- as.numeric(input$yLimVolcano)
    volcanoTable$log2Ratio <- shrinkToRange(volcanoTable$log2Ratio,
                                            theRange = c(-as.numeric(input$xLimVolcano), as.numeric(input$xLimVolcano)))

    # Summary table
    output$volcanoOverview <- function() {
      table(as.factor(volcanoTable$Status)) %>%
        kable(
          format = "html",
          caption = paste("DE Summary with Volcano Settings")
        ) %>%
        kable_styling(
          bootstrap_options = "striped",
          full_width = FALSE,
          position = "left"
        )
    }

    # Results table
    output$volcanoOverviewTable <- DT::renderDataTable({
      DT::datatable(
        data = volcanoTable[,c("gene_name", "description", "log2Ratio", pTypeVolcano, "log10p", "Status")],
        filter = "top",
        class = "cell-border stripe",
        rownames = FALSE,
        colnames = c("Feature Symbol", "Description", "Log2 Ratio", pTypeVolcano, paste("-log10", pTypeVolcano), "Significance")
      ) %>%
        DT::formatSignif(columns = c("log10p", pTypeVolcano, "log2Ratio"), digits = 3) %>%
        DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
        DT::formatStyle(columns = "log2Ratio", color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold") %>%
        DT::formatStyle(columns = "log10p", color = styleInterval(cuts = log10pToFilter, values = c("black", "green")), fontWeight = "bold") %>%
        DT::formatStyle(columns = "fdr", color = styleInterval(cuts = input$pThresholdVolcano, values = c("green", "black")), fontWeight = "bold") %>%
        DT::formatStyle(columns = "Status", color = styleEqual(levels = names(volcanoColours), values = as.character(volcanoColours)), fontWeight = "bold")
    })

    return(list(
      "volcanoTable" = volcanoTable,
      "volcanoColours" = volcanoColours,
      "design" = design
    ))
  }
)

observeEvent(
  {
    input$textSizeVolcano
    input$showLinesVolcano
    input$showAxesVolcano
    input$borderVolcano
    input$alphaVolcano
    input$dotSizeVolcano
    input$volcanoShowGenes
    input$volcanoKeepBucketGenes
    input$volcanoLabelAllUp
    input$volcanoLabelAllDown
    input$volcanoAnnotationHighlightColour
    input$volcanoLabelMaxOverlap
    input$showLinesVolcano
    input$geneLabelSizeVolcano
    input$geneLabelNudgeVolcanoX
    input$geneLabelNudgeVolcanoY
    input$volcanoPointBorder
    input$downloadFormatVolcano
    input$dpiVolcano
    input$filnameVolcano
    volcanoResultsList()
  }, ignoreInit = TRUE, ignoreNULL = FALSE,
  {

    req(!is.null(volcanoResultsList()$volcanoColours))
    volcanoTable <- volcanoResultsList()$volcanoTable
    volcanoColours <- volcanoResultsList()$volcanoColours
    design <- volcanoResultsList()$design

    volcanoTableFull <- volcanoTable
    volcanoTableFull$Label <- NA
    if (input$volcanoShowGenes & input$volcanoLabelAllUp) {
      volcanoTableFull$Label[which(volcanoTableFull$Status == "FoldChangeSignificantUp")] <- volcanoTableFull$gene_name[which(volcanoTableFull$Status == "FoldChangeSignificantUp")]
    }
    if (input$volcanoShowGenes & input$volcanoLabelAllDown) {
      volcanoTableFull$Label[which(volcanoTableFull$Status == "FoldChangeSignificantDown")] <- volcanoTableFull$gene_name[which(volcanoTableFull$Status == "FoldChangeSignificantDown")]
    }
    if (input$volcanoShowGenes & length(input$volcanoKeepBucketGenes) >= 1) {
      volcanoTableFull$Label[which(volcanoTableFull$gene_name %in% input$volcanoKeepBucketGenes)] <- volcanoTableFull$gene_name[which(volcanoTableFull$gene_name %in% input$volcanoKeepBucketGenes)]
      if (input$volcanoAnnotationHighlightColour) {
        volcanoTableFull$Status[which(volcanoTableFull$gene_name %in% input$volcanoKeepBucketGenes)] <- "Highlight"
      }
    }
    
    # Create the MA plot
    maplot <- ggplot(volcanoTableFull, aes(x = Log2_Mean, y = log2Ratio_full))
    maplot <- maplot +
      theme_prism(base_size = input$textSizeVolcano) +
      xlab("Log2 Normalised Mean") + ylab("Log2 Ratio")
    maplot <- maplot + geom_point(pch = 21, size = as.numeric(input$dotSizeVolcano), alpha = as.numeric(input$alphaVolcano), aes(fill = Status), stroke = input$volcanoPointBorder)
    maplot <- maplot + scale_fill_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
    if (input$volcanoShowGenes == TRUE) {
      maplot <- maplot +
        geom_label_repel(
          aes(label = Label, colour = Status),
          force = 2,
          nudge_y = as.numeric(input$geneLabelNudgeVolcanoY)/10,
          nudge_x = as.numeric(input$geneLabelNudgeVolcanoX)/10,
          size = (input$geneLabelSizeVolcano / 3),
          max.overlaps = input$volcanoLabelMaxOverlap,
          fontface = "bold",
          show.legend = F
      ) +
        scale_colour_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
    }
    maplot <- maplot + guides()
    figuresDataReactive$volcanoMA <- maplot
    
    
    volcanoTableFullMA <- volcanoTableFull %>% dplyr::select(gene_name, description, log2Ratio_full, Log2_Mean, log10p_full, Status, Label)
    output$MABrushTable <- renderDataTable({
      DT::datatable(
        data = volcanoTableFullMA[!is.na(volcanoTableFullMA$Label),],
        colnames = c("Feature name", "Description", "Log2 Ratio", "Log2 Mean", "-Log10p", "Status"),
        filter = "top",
        caption = "Click and drag over dots on the volcano plot to see those features in this table.",
        rownames = F) %>%
        DT::formatSignif(c("log2Ratio_full", "log10p_full", "Log2_Mean"), digits = 3) %>%
        DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
        DT::formatStyle(columns = c("log2Ratio_full"), color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold") %>%
        DT::formatStyle(columns = c("log10p_full"), color = styleInterval(cuts = -log10(as.numeric(input$pThresholdVolcano)), values = c("black", "green")), fontWeight = "bold") %>%
        DT::formatStyle(columns = "Status", color = styleEqual(levels = names(volcanoColours), values = as.character(volcanoColours)), fontWeight = "bold")
    })

    # Create the volcano plot
    volcanoStatic <- ggplot(
      data = volcanoTableFull, aes(x = log2Ratio, y = log10p, text = gene_name)
    )
    volcanoStatic <- volcanoStatic + ylim(0, input$yLimVolcano) + xlim(-input$xLimVolcano, input$xLimVolcano)
    volcanoStatic <- volcanoStatic + theme_prism(base_size = input$textSizeVolcano, border = input$borderVolcano)
    volcanoStatic <- volcanoStatic + geom_vline(xintercept = c(-as.numeric(input$lfcVolcano), as.numeric(input$lfcVolcano)), col = "black", linetype = "dashed")
    volcanoStatic <- volcanoStatic + geom_hline(yintercept = -log10(as.numeric(input$pThresholdVolcano)), col = "black", linetype = "dashed")
    if (input$volcanoPointBorder > 0) {
      volcanoStatic <- volcanoStatic + geom_point(aes(fill = Status), size = as.numeric(input$dotSizeVolcano), alpha = as.numeric(input$alphaVolcano), pch = 21, stroke = input$volcanoPointBorder)
      volcanoStatic <- volcanoStatic + scale_fill_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
      volcanoStatic <- volcanoStatic + guides(fill = guide_legend(override.aes = list(shape = 21,  size = 4, stroke = input$volcanoPointBorder)))
    } else if (input$volcanoPointBorder == 0) {
      volcanoStatic <- volcanoStatic + geom_point(aes(colour = Status), size = as.numeric(input$dotSizeVolcano)*0.8, alpha = as.numeric(input$alphaVolcano), pch = 19)
      volcanoStatic <- volcanoStatic + scale_colour_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
      volcanoStatic <- volcanoStatic + guides(colour = guide_legend(override.aes = list(shape = 19,  size = 4)))
    }
    if (input$volcanoShowGenes) {
      volcanoStatic <- volcanoStatic +
        geom_label_repel(
          aes(label = Label, colour = Status),
          force = 3,
          nudge_y = as.numeric(input$geneLabelNudgeVolcanoY)/10,
          nudge_x = as.numeric(input$geneLabelNudgeVolcanoX)/10,
          size = (input$geneLabelSizeVolcano / 3),
          max.overlaps = input$volcanoLabelMaxOverlap,
          fontface = "bold",
          show.legend = F
        ) + scale_colour_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
    }
    volcanoStatic <- volcanoStatic + labs(
      x = "Log2 Fold Change",
      y = "-log10 p-value",
      fill = "Significance"
    )

    if (!input$showAxesVolcano) {
      volcanoStatic <- volcanoStatic +
        theme(
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank()
        )
    }
    if (!input$showLinesVolcano) {
      volcanoStatic <- volcanoStatic +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
        )
    }
    figuresDataReactive$volcanoStatic <- volcanoStatic
    
    volcanoPlotly <- plot_ly(
      data = volcanoTableFull, 
      x = ~log2Ratio, 
      y = ~log10p, 
      type = 'scattergl', 
      mode = 'markers', 
      color = ~Status, 
      colors = volcanoColours,
      text = ~gene_name,
      hovertemplate = paste(
        "Gene: %{text}<br>",
        "log2FC: %{x:.2f}<br>",
        "-log10(p): %{y:.2f}<br>",
        "<extra></extra>"
      ),
      marker = list(
        size = input$dotSizeVolcano,
        color = "fill_colour",
        line = list(
          color = "black",
          width = input$volcanoPointBorder
        )
      )
    )
    figuresDataReactive$volcanoPlotly <- volcanoPlotly
    
    volcanoTableFull2 <- volcanoTableFull %>% dplyr::select(gene_name, description, log2Ratio, log2Ratio_full, log10p, log10p_full, Status, Label)
    output$volcanoBrushTable <- renderDataTable({
      DT::datatable(
        data = volcanoTableFull2[!is.na(volcanoTableFull2$Label),],
        colnames = c("Feature name", "Description", "Log2 Ratio (plot)", "Log2 Ratio (full)", "-Log10p (plot)", "-Log10p (full)", "Status"),
        filter = "top",
        caption = "Click and drag over dots on the volcano plot to see those features in this table.",
        rownames = F) %>%
        DT::formatSignif(c("log2Ratio", "log2Ratio_full", "log10p", "log10p_full"), digits = 3) %>%
        DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
        DT::formatStyle(columns = c("log2Ratio", "log2Ratio_full"), color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold") %>%
        DT::formatStyle(columns = c("log10p", "log10p_full"), color = styleInterval(cuts = -log10(as.numeric(input$pThresholdVolcano)), values = c("black", "green")), fontWeight = "bold") %>%
        DT::formatStyle(columns = "Status", color = styleEqual(levels = names(volcanoColours), values = as.character(volcanoColours)), fontWeight = "bold")
    })

    output$dlVolcanoDFButton <- downloadHandler(
      filename = function() {
        paste0(input$filnameVolcano, "_volcano.xlsx")
      },
      content = function(file) {
        openxlsx::write.xlsx(volcanoTable, file = file)
      }
    )
  }
)

# Render the volcano static plot
output$volcanoStatic <- renderPlot(
  {
    req(!is.null(figuresDataReactive$volcanoStatic))
    figuresDataReactive$volcanoStatic
  },
  width = function(){as.numeric(input$figWidthVolcano)},
  height = function(){as.numeric(input$figHeightVolcano)}
)
output$volcanoPlotly <- renderPlotly(
  {
    req(!is.null(figuresDataReactive$volcanoPlotly))
    volcanoPlotly <- figuresDataReactive$volcanoPlotly
    volcanoPlotly <- config(
      volcanoPlotly,
      displaylogo = FALSE,
      displayModeBar = TRUE,
      modeBarButtonsToRemove = c("toImage", "select2d", "lasso2d")
    )
    volcanoPlotly %>%
      layout(height = input$figHeightVolcano, width = input$figWidthVolcano)
  }
)
# Download button for the plot
output$dlVolcanoPlotButton <- downloadHandler(
  filename = function() {
    paste0(input$filnameVolcano, "_volcano.", tolower(input$downloadFormatVolcano))
  },
  content = function(file) {
    if (input$downloadFormatVolcano == "PDF") {
      pdf(file = file, width = as.numeric(input$figWidthVolcano/80), height = as.numeric(input$figHeightVolcano/80))
    } else if (input$downloadFormatVolcano == "SVG") {
      svg(file = file, width = as.numeric(input$figWidthVolcano/80), height = as.numeric(input$figHeightVolcano/80))
    } else if (input$downloadFormatVolcano == "PNG") {
      png(filename = file, width = as.numeric(input$figWidthVolcano/80), height = as.numeric(input$figHeightVolcano/80), units = "in", res = as.numeric(input$dpiVolcano))
    }
    print(figuresDataReactive$volcanoStatic)
    dev.off()
  }
)

output$MAPlot <- renderPlot({
  req(!is.null(figuresDataReactive$volcanoMA))
  figuresDataReactive$volcanoMA
}, height = as.numeric(input$figHeightVolcano)/1.5, width = as.numeric(input$figWidthVolcano)*1.2)
output$dlMAPlotButton <- downloadHandler(
  filename = function() {
    paste0(input$filnameVolcano, "_MA.pdf")
  },
  content = function(file) {
    ggsave(
      filename = file,
      plot = figuresDataReactive$volcanoMA,
      width = (as.numeric(input$figWidthVolcano) / 102),
      height = (as.numeric(input$figHeightVolcano)/ 115),
      limitsize = FALSE, units = "in"
    )
  }
)