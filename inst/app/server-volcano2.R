eds <- inputDataReactive()$.eds
req(eds)

# FDR-only producers (no raw p-value, e.g. SAINTexpress) lock the volcano
# significance axis to FDR rather than plotting an all-NA -log10(pValue).
if (!isTRUE(eds$has_pvalue())) {
  updateSelectInput(session = session, inputId = "pTypeVolcano", choices = c("FDR"), selected = "FDR")
}

if (eds$omics_class() == "proteomics") {
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

if (eds$omics_class() == "proteomics" && eds$has_de()) {
  if ("modelName" %in% colnames(eds$de(eds$primary_contrast()))) {
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
  inputId = "fb_select_volcano",
  choices = eds$feature_ids(),
  selected = "",
  server = TRUE
)

# Bucket render, reset handler removed — now in server-featureBucket.R

observeEvent(
  {
    input$pTypeVolcano
    input$contrastSelected
  },
  ignoreInit = TRUE,
  ignoreNULL = TRUE,
  {
    # Volcano plots are only applicable for RNASeq and proteomics data with DE results
    # CountQC data does not have differential expression results
    req(eds$data_type() == "RNASeq" || eds$omics_class() == "proteomics")
    req(eds$has_de())   # counts-only uploads carry no contrasts

    if (eds$data_type() == "RNASeq") {
      contrast    <- selectedContrast()
      design      <- eds$contrast_label(contrast)
      seqAnno     <- eds$de(contrast)
      output$volcanoDesign <- renderText({
        design
      })
      seqAnnoFilt <- seqAnno[which(seqAnno$usedInTest), ]
      updateTextInput(session = session, inputId = "filename_volcano", value = design)
    }

    if (eds$omics_class() == "proteomics") {
      req(!is.null(input$contrastSelected))
      req(input$contrastSelected %in% eds$contrasts())
      seqAnnoFilt <- eds$de(input$contrastSelected)
      seqAnnoFilt <- seqAnnoFilt %>% dplyr::select(feature_id, log2Ratio, pValue, fdr)

      output$volcanoDesign <- renderText({
        input$contrastSelected
      })
      updateTextInput(session = session, inputId = "filename_volcano", value = input$contrastSelected)
    }

    if (input$pTypeVolcano == "Raw") {
      pTypeVolcano <- "pValue"
    } else {
      pTypeVolcano <- "fdr"
    }

    updateNumericInput(
      session = session,
      inputId = "yLimVolcano",
      value = ceiling(max(-log10(seqAnnoFilt[[pTypeVolcano]])) * 1.1)
    )
    updateNumericInput(
      session = session,
      inputId = "xLimVolcano",
      value = ceiling(max(abs(seqAnnoFilt[["log2Ratio"]])) * 1.1)
    )
  }
)

if (eds$omics_class() == "proteomics" && eds$has_de()) {
  proto_de <- eds$de(eds$primary_contrast())

  if ("nrPeptides" %in% colnames(proto_de)) {
    output$nrPeptidesVolcanoUI <- renderUI({
      sliderInput(inputId = "nrPeptideVolcano", label = "Minimum number of peptides", min = 0, max = 20, value = 0, step = 1, width = "85%")
    })
    outputOptions(output, "nrPeptidesVolcanoUI", suspendWhenHidden = FALSE)
  } else {
    output$nrPeptidesVolcanoUI <- renderUI({
      NULL
    })
  }

  if ("modelName" %in% colnames(proto_de)) {
    output$showImputedVolcanoUI <- renderUI({
      checkboxInput(inputId = "showImputedVolcano", label = "Include imputed features?", value = TRUE)
    })
    output$highlightImputedVolcanoUI <- renderUI({
      checkboxInput(inputId = "highlightImputedVolcano", label = "Highlight imputed features?", value = FALSE)
    })
    outputOptions(output, "showImputedVolcanoUI", suspendWhenHidden = FALSE)
    outputOptions(output, "highlightImputedVolcanoUI", suspendWhenHidden = FALSE)
  } else {
    output$showImputedVolcanoUI <- renderUI({
      NULL
    })
    output$highlightImputedVolcanoUI <- renderUI({
      NULL
    })
  }
}

volcanoStagedGenes <- reactiveValues(brushGenes = character(0))

observeEvent(
  {
    input$volcanoBrush
    input$MABrush
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE,
  {
    volcanoBrushGenes <- brushedPoints(volcanoResultsList()$volcanoTable, input$volcanoBrush, xvar = "log2Ratio", yvar = "log10p")
    volcanoBrushGenes <- volcanoBrushGenes$feature_id
    maBrushGenes <- brushedPoints(volcanoResultsList()$volcanoTable, input$MABrush, xvar = "Log2_Mean", yvar = "log2Ratio_full")
    maBrushGenes <- maBrushGenes$feature_id
    volcanoStagedGenes$brushGenes <- unique(c(volcanoBrushGenes, maBrushGenes))
  }
)

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
    featureBucket$features
    input$fb_exclude_volcano
    input$contrastSelected
    input$nrPeptideVolcano
    input$showImputedVolcano
    input$highlightImputedVolcano
    lapply(names(volcanoColoursDefault), function(x) {
      input[[paste0("volcanoColour", x)]]
    })
  },
  ignoreInit = FALSE,
  ignoreNULL = FALSE,
  {
    # Volcano plots are only applicable for RNASeq and proteomics data with DE results
    # CountQC data does not have differential expression results
    req(eds$data_type() == "RNASeq" || eds$omics_class() == "proteomics")
    req(eds$has_de())   # counts-only uploads carry no contrasts

    if (eds$data_type() == "RNASeq") {
      contrast    <- selectedContrast()
      design      <- eds$contrast_label(contrast)
      seqAnno     <- eds$de(contrast)
      output$volcanoDesign <- renderText({
        design
      })
      seqAnnoFilt <- seqAnno[which(seqAnno$usedInTest), ]
    }
    if (eds$omics_class() == "proteomics") {
      req(!is.null(input$contrastSelected))
      req(input$contrastSelected %in% eds$contrasts())
      seqAnnoFilt <- eds$de(input$contrastSelected)
      if (!is.null(input$nrPeptideVolcano)) {
        if (any(seqAnnoFilt$nrPeptides >= input$nrPeptideVolcano)) {
          seqAnnoFilt <- seqAnnoFilt[which(seqAnnoFilt$nrPeptides >= input$nrPeptideVolcano), ]
        } else {
          shinyalert::shinyalert(title = "Oops!", text = "No features with this number of peptides!", type = "error", closeOnClickOutside = TRUE, showCancelButton = FALSE, showConfirmButton = TRUE, timer = 5000)
        }
      }
      if (!is.null(input$showImputedVolcano)) {
        if (!input$showImputedVolcano) {
          seqAnnoFilt <- seqAnnoFilt[grep("^Linear", seqAnnoFilt$modelName), ]
        }
      }
      if (!is.null(input$highlightImputedVolcano)) {
        if (input$highlightImputedVolcano) {
          seqAnnoFilt <- seqAnnoFilt %>% dplyr::select(feature_id, description, log2Ratio, pValue, fdr, modelName)
        }
      } else {
        seqAnnoFilt <- seqAnnoFilt %>% dplyr::select(feature_id, description, log2Ratio, pValue, fdr)
      }

      output$volcanoDesign <- renderText({
        input$contrastSelected
      })
    }

    if (length(input$fb_exclude_volcano) > 0) {
      seqAnnoFilt <- seqAnnoFilt[!seqAnnoFilt$feature_id %in% input$fb_exclude_volcano, ]
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
    volcanoTable <- seqAnnoFilt[, c("feature_id", "description", "log2Ratio", pTypeVolcano)]

    # Get the -log10 of the p-value:
    volcanoTable["log10p"] <- -log10(as.numeric(volcanoTable[[pTypeVolcano]]))

    # Get the -log10p of p-value threshold as a single value:
    log10pToFilter <- -log10(as.numeric(input$pThresholdVolcano))

    # Get mean log2 count for MA plot
    if (eds$data_type() == "RNASeq") {
      rMs <- rowMeans(eds$assay("normalised_log2"))
    } else {
      rMs <- rowMeans(eds$assay("intensities"))
    }
    volcanoTable$Log2_Mean <- rMs[match(volcanoTable$feature_id, names(rMs))]

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
        volcanoTable[which(volcanoTable$feature_id %in% seqAnnoFilt$feature_id[grep("^Imputed", seqAnnoFilt$modelName)]), "Status"] <- "Imputed"
      }
    }

    # Set limits of data frame to x/y-axes limits
    volcanoTable$log10p_full <- volcanoTable$log10p
    volcanoTable$log2Ratio_full <- volcanoTable$log2Ratio
    volcanoTable$log10p[volcanoTable$log10p > as.numeric(input$yLimVolcano)] <- as.numeric(input$yLimVolcano)
    volcanoTable$log2Ratio <- shrinkToRange(volcanoTable$log2Ratio,
      theRange = c(-as.numeric(input$xLimVolcano), as.numeric(input$xLimVolcano))
    )

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
        data = volcanoTable[, c("feature_id", "description", "log2Ratio", pTypeVolcano, "log10p", "Status")],
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
    featureBucket$features
    volcanoStagedGenes$brushGenes
    input$volcanoLabelAllUp
    input$volcanoLabelAllDown
    input$volcanoAnnotationHighlightColour
    input$volcanoLabelMaxOverlap
    input$showLinesVolcano
    input$geneLabelSizeVolcano
    input$geneLabelNudgeVolcanoX
    input$geneLabelNudgeVolcanoY
    input$volcanoPointBorder
    volcanoResultsList()
  },
  ignoreInit = TRUE,
  ignoreNULL = FALSE,
  {
    req(!is.null(volcanoResultsList()$volcanoColours))
    volcanoTable <- volcanoResultsList()$volcanoTable
    volcanoColours <- volcanoResultsList()$volcanoColours
    design <- volcanoResultsList()$design

    volcanoTableFull <- volcanoTable
    volcanoTableFull$Label <- NA
    if (input$volcanoShowGenes & input$volcanoLabelAllUp) {
      volcanoTableFull$Label[which(volcanoTableFull$Status == "FoldChangeSignificantUp")] <- volcanoTableFull$feature_id[which(volcanoTableFull$Status == "FoldChangeSignificantUp")]
    }
    if (input$volcanoShowGenes & input$volcanoLabelAllDown) {
      volcanoTableFull$Label[which(volcanoTableFull$Status == "FoldChangeSignificantDown")] <- volcanoTableFull$feature_id[which(volcanoTableFull$Status == "FoldChangeSignificantDown")]
    }
    # Annotate bucket genes AND staged (brushed but not yet added) genes
    allAnnotated <- unique(c(featureBucket$features, volcanoStagedGenes$brushGenes))
    if (input$volcanoShowGenes & length(allAnnotated) >= 1) {
      volcanoTableFull$Label[which(volcanoTableFull$feature_id %in% allAnnotated)] <- volcanoTableFull$feature_id[which(volcanoTableFull$feature_id %in% allAnnotated)]
      if (input$volcanoAnnotationHighlightColour) {
        volcanoTableFull$Status[which(volcanoTableFull$feature_id %in% featureBucket$features)] <- "Highlight"
      }
    }

    figuresDataReactive$volcanoTableFull <- volcanoTableFull

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
          nudge_y = as.numeric(input$geneLabelNudgeVolcanoY) / 10,
          nudge_x = as.numeric(input$geneLabelNudgeVolcanoX) / 10,
          size = (input$geneLabelSizeVolcano / 3),
          max.overlaps = input$volcanoLabelMaxOverlap,
          fontface = "bold",
          show.legend = F
        ) +
        scale_colour_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
    }
    maplot <- maplot + guides()
    figuresDataReactive$volcanoMA <- maplot

    # Create the volcano plot
    volcanoStatic <- ggplot(
      data = volcanoTableFull, aes(x = log2Ratio, y = log10p, text = feature_id)
    )
    volcanoStatic <- volcanoStatic + ylim(0, input$yLimVolcano) + xlim(-input$xLimVolcano, input$xLimVolcano)
    volcanoStatic <- volcanoStatic + theme_prism(base_size = input$textSizeVolcano, border = input$borderVolcano)
    volcanoStatic <- volcanoStatic + geom_vline(xintercept = c(-as.numeric(input$lfcVolcano), as.numeric(input$lfcVolcano)), col = "black", linetype = "dashed")
    volcanoStatic <- volcanoStatic + geom_hline(yintercept = -log10(as.numeric(input$pThresholdVolcano)), col = "black", linetype = "dashed")
    if (input$volcanoPointBorder > 0) {
      volcanoStatic <- volcanoStatic + geom_point(aes(fill = Status), size = as.numeric(input$dotSizeVolcano), alpha = as.numeric(input$alphaVolcano), pch = 21, stroke = input$volcanoPointBorder)
      volcanoStatic <- volcanoStatic + scale_fill_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
      volcanoStatic <- volcanoStatic + guides(fill = guide_legend(override.aes = list(shape = 21, size = 4, stroke = input$volcanoPointBorder)))
    } else if (input$volcanoPointBorder == 0) {
      volcanoStatic <- volcanoStatic + geom_point(aes(colour = Status), size = as.numeric(input$dotSizeVolcano) * 0.8, alpha = as.numeric(input$alphaVolcano), pch = 19)
      volcanoStatic <- volcanoStatic + scale_colour_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
      volcanoStatic <- volcanoStatic + guides(colour = guide_legend(override.aes = list(shape = 19, size = 4)))
    }
    if (input$volcanoShowGenes) {
      volcanoStatic <- volcanoStatic +
        geom_label_repel(
          aes(label = Label, colour = Status),
          force = 3,
          nudge_y = as.numeric(input$geneLabelNudgeVolcanoY) / 10,
          nudge_x = as.numeric(input$geneLabelNudgeVolcanoX) / 10,
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
      type = "scattergl",
      mode = "markers",
      color = ~Status,
      colors = volcanoColours,
      text = ~feature_id,
      hovertemplate = paste(
        "Gene: %{text}<br>",
        "log2FC: %{x:.2f}<br>",
        "-log10(p): %{y:.2f}<br>",
        "<extra></extra>"
      ),
      marker = list(
        size = input$dotSizeVolcano * 2,
        color = "fill_colour",
        line = list(
          color = "black",
          width = input$volcanoPointBorder
        )
      )
    )
    figuresDataReactive$volcanoPlotly <- volcanoPlotly

    output$dlTable_volcano <- downloadHandler(
      filename = function() {
        generate_filename(
          base_name = paste0(input$filename_volcano, "_Results"),
          design_name = NULL,
          extension = "xlsx"
        )
      },
      content = function(file) {
        export_excel_standard(
          df = volcanoTable,
          filename = file,
          sheet_name = "Volcano Results",
          table_name = "VolcanoResults"
        )
      }
    )
  }
)

output$volcanoBrushTable <- renderDataTable({
  req(figuresDataReactive$volcanoTableFull)
  volcanoTableFull <- figuresDataReactive$volcanoTableFull
  volcanoTableFull2 <- volcanoTableFull %>% dplyr::select(feature_id, description, log2Ratio, log2Ratio_full, log10p, log10p_full, Status)

  DT::datatable(
    data = brushedPoints(volcanoTableFull2, input$volcanoBrush),
    colnames = c("Feature name", "Description", "Log2 Ratio (plot)", "Log2 Ratio (full)", "-Log10p (plot)", "-Log10p (full)", "Status"),
    filter = "top",
    caption = "Click and drag over dots on the volcano plot to see those features in this table.",
    rownames = F
  ) %>%
    DT::formatSignif(c("log2Ratio", "log2Ratio_full", "log10p", "log10p_full"), digits = 3) %>%
    DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
    DT::formatStyle(columns = c("log2Ratio", "log2Ratio_full"), color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold") %>%
    DT::formatStyle(columns = c("log10p", "log10p_full"), color = styleInterval(cuts = -log10(as.numeric(input$pThresholdVolcano)), values = c("black", "green")), fontWeight = "bold") %>%
    DT::formatStyle(columns = "Status", color = styleEqual(levels = names(volcanoResultsList()$volcanoColours), values = as.character(volcanoResultsList()$volcanoColours)), fontWeight = "bold")
})

output$MABrushTable <- renderDataTable({
  req(figuresDataReactive$volcanoTableFull)
  volcanoTableFull <- figuresDataReactive$volcanoTableFull
  volcanoTableFullMA <- volcanoTableFull %>% dplyr::select(feature_id, description, log2Ratio_full, Log2_Mean, log10p_full, Status)

  DT::datatable(
    data = brushedPoints(volcanoTableFullMA, input$MABrush),
    colnames = c("Feature name", "Description", "Log2 Ratio", "Log2 Mean", "-Log10p", "Status"),
    filter = "top",
    caption = "Click and drag over dots on the volcano plot to see those features in this table.",
    rownames = F
  ) %>%
    DT::formatSignif(c("log2Ratio_full", "log10p_full", "Log2_Mean"), digits = 3) %>%
    DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
    DT::formatStyle(columns = c("log2Ratio_full"), color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold") %>%
    DT::formatStyle(columns = c("log10p_full"), color = styleInterval(cuts = -log10(as.numeric(input$pThresholdVolcano)), values = c("black", "green")), fontWeight = "bold") %>%
    DT::formatStyle(columns = "Status", color = styleEqual(levels = names(volcanoResultsList()$volcanoColours), values = as.character(volcanoResultsList()$volcanoColours)), fontWeight = "bold")
})
output$volcanoStatic <- renderPlot(
  {
    req(!is.null(figuresDataReactive$volcanoStatic))
    figuresDataReactive$volcanoStatic
  },
  width = function() {
    as.numeric(input$figWidthVolcano)
  },
  height = function() {
    as.numeric(input$figHeightVolcano)
  }
)
output$volcanoPlotly <- renderPlotly({
  req(!is.null(figuresDataReactive$volcanoPlotly))
  volcanoPlotly <- figuresDataReactive$volcanoPlotly
  volcanoPlotly <- plotly::config(
    volcanoPlotly,
    displaylogo = FALSE,
    displayModeBar = TRUE,
    modeBarButtonsToRemove = c("toImage", "select2d", "lasso2d")
  )
  volcanoPlotly %>%
    plotly::layout(height = input$figHeightVolcano, width = input$figWidthVolcano)
})
# Download button for the plot
# Download button for the plot
output$dlPlot_volcano <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = "volcano",
      design_name = input$filename_volcano,
      extension = tolower(input$downloadFormatVolcano)
    )
  },
  content = function(file) {
    export_plot_standard(
      plot_object = figuresDataReactive$volcanoStatic,
      filename = file,
      format = input$downloadFormatVolcano,
      dpi = input$dpiVolcano,
      width = as.numeric(input$figWidthVolcano / 80),
      height = as.numeric(input$figHeightVolcano / 80)
    )
  }
)

output$MAPlot <- renderPlot(
  {
    req(!is.null(figuresDataReactive$volcanoMA))
    figuresDataReactive$volcanoMA
  },
  height = as.numeric(input$figHeightVolcano) / 1.5,
  width = as.numeric(input$figWidthVolcano) * 1.2
)
output$dlPlot_MA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_volcano, "_MA"),
      design_name = NULL,
      extension = tolower(input$downloadFormatVolcano)
    )
  },
  content = function(file) {
    export_plot_standard(
      plot_object = figuresDataReactive$volcanoMA,
      filename = file,
      format = tolower(input$downloadFormatVolcano),
      width = as.numeric(input$figWidthVolcano / 80),
      height = as.numeric(input$figHeightVolcano / 80),
      dpi = as.numeric(input$dpiVolcano)
    )
  }
)
