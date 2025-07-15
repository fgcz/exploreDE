dataSummary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sqrt(sd(x[[col]]) / length(x[[col]])))
  }
  data_sum <- ddply(
    data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

zscore <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

# Update a bunch of input options ----
output$boxplotStatic <- renderPlot({NULL})

updateSelectInput(
  session = session, 
  inputId = "boxplotCounts", 
  choices = names(inputDataReactive()$countList), 
  selected = switch(
    inputDataReactive()$dataType, 
    "RNASeq" = "Normalised + Log2", 
    "proteomics" = "transformedData")
)

updateSelectizeInput(
  session = session,
  inputId = "boxplotGenes",
  choices = inputDataReactive()$genes,
  selected = "",
  server = TRUE
)

updateSelectInput(
  session = session,
  inputId = "boxplotFactor1",
  choices = inputDataReactive()$factors,
  selected = inputDataReactive()$factors[1]
)
updateSelectInput(
  session = session,
  inputId = "boxplotFactor2",
  choices = c("None", "Feature", inputDataReactive()$factors),
  selected = "None"
)

updateSelectInput(
  session = session,
  inputId = "boxplotBatch",
  choices = c("None", inputDataReactive()$factors),
  selected = "None"
)

observeEvent({
  input$boxplotFactor1
}, {
  output$bucket <- renderUI({
    bucket_list(
      header = "Drag and drop groups in order to be plotted",
      group_name = "bucket_list_group",
      orientation = "horizontal",
      add_rank_list(
        text = "Include these groups",
        labels = as.list(unique(levels(as.factor(inputDataReactive()$dataset[[input$boxplotFactor1]])))),
        input_id = "boxKeepBucket"),
      add_rank_list(
        text = "Exclude these groups",
        labels = NULL,
        input_id = "boxExcludeBucket")
    )  
  })
  outputOptions(output, "bucket", suspendWhenHidden = FALSE)
})

# Extract the appropriate data ----
if (inputDataReactive()$dataType == "RNASeq") {
  seqAnnoReactive <- reactiveValues("sa" = inputDataReactive()$seqAnno, "saF" = NULL)
  design <- inputDataReactive()$design
  designLevels <- inputDataReactive()$designLevels
  output$boxplotDesign <- renderText({inputDataReactive()$design})
} else if (inputDataReactive()$dataType == "proteomics") {
  seqAnnoReactive <- reactiveValues("sa" = NULL, "saF" = NULL)
  observeEvent(input$contrastSelected, {
    seqAnnoReactive$sa = inputDataReactive()$seqAnnoList[[input$contrastSelected]]
    seqAnnoReactive$sa <- seqAnnoReactive$sa %>% dplyr::select("gene_name", "log2Ratio", "pValue", "fdr")
    seqAnnoReactive$saF <- seqAnnoReactive$sa
    design <- input$contrastSelected
    output$boxplotDesign <- renderText({input$contrastSelected})
  }, ignoreInit = T, ignoreNULL = T)
}

# The main part of the gene bucket lives here ----
# It listens in on DEG table, volcano, and heatmap tabs, in addition to this tab
observeEvent({
  input$boxplotGenes
  input$boxplotGenesText
  input$degTable_rows_selected
  input$volcanoGenes
  input$volcanoGenesText
  input$heatmapGenes
  input$heatmapGenesText
  volcanoGenesReactive$volcanoBrushGenes
}, ignoreNULL = FALSE, ignoreInit = TRUE, {
  
  # Get all the genes that have been selected from the drop-downs of a respective tab
  genesList1 <- setNames(lapply(c("boxplotGenes", "volcanoGenes", "heatmapGenes"), function(geneList) {
    if(!is.null(input[[geneList]])) {
      input[[geneList]]
    }
  }), c("boxplotGenes", "volcanoGenes", "heatmapGenes"))
  
  # Get all the genes that have been typed into the textbox of a respective tab
  genesList2 <- setNames(lapply(c(
    "boxplotGenesText", "volcanoGenesText", "heatmapGenesText"
  ), function(geneText) {
    if(!is.null(input[[geneText]])) {
      strsplit(input[[geneText]], "\\ +|\\/+|\\,+|\n+|\r+|\n\r+|\t+")[[1]]
    }
  }), c("boxplotGenesText", "volcanoGenesText", "heatmapGenesText"))
  
  # Get the genes that were click-selected in the DE table 
  degMATable <- seqAnnoReactive$sa[order(seqAnnoReactive$sa$pValue), ]
  degMATable <- degMATable$gene_name[input$degTable_rows_selected]
  
  volcanoBrushGenes <- NULL
  if (!is.null(volcanoGenesReactive$volcanoBrushGenes)) {
    volcanoBrushGenes = volcanoGenesReactive$volcanoBrushGenes
  }
  
  genesUnlist <- unique(unlist(c(genesList1, genesList2, degMATable, volcanoBrushGenes)))
  genesUnlist <- genesUnlist[!genesUnlist == ""]
  
  genesUnlist <- genesUnlist[which(genesUnlist %in% seqAnnoReactive$sa$gene_name)]
  genesReactive$genes <- genesUnlist
})

observeEvent(input$resetGeneBucketBoxplot, {
  genesReactive$genes <- NULL
  updateSelectizeInput(session, inputId = "boxplotGenes", choices = inputDataReactive()$genes, selected = NULL, server = T)
  updateTextAreaInput(session, inputId = "boxplotGenesText", value = NULL)
  output$boxplotStatic <- NULL
}, ignoreNULL = T, ignoreInit = T)


observeEvent({genesReactive}, {
  output$geneBucket1 <- renderUI({
    bucket_list(
      header = "Drag and drop features in order to be plotted",
      group_name = "bucket_list_group",
      orientation = "horizontal",
      add_rank_list(
        text = "Include these features in this order",
        labels = genesReactive$genes,
        input_id = "keepBucketBoxplot"),
      add_rank_list(
        text = "Exclude these features",
        labels = NULL,
        input_id = "excludeBucketBoxplot")
    )
  })
  outputOptions(output, "geneBucket1", suspendWhenHidden = FALSE)
})

# Generate a reactive list of the count data to be plotted ----
# We do this so we can make aesthetic changes to the plot without having to recalculate the required boxplots and results 

debouncedBoxplotFactor1 <- reactive({
  input$keepBucketBoxplot
}) %>% debounce(200) 

boxplotCountsReactive <- eventReactive({
  debouncedBoxplotFactor1()
  input$boxplotFactor1
  input$boxKeepBucket
  input$boxplotFactor2
  input$boxplotBatch
  input$boxplotCounts
  input$contrastSelected
  input$boxplotShowPHeight
  input$boxplotCountsLog
  input$boxplotZScore
}, ignoreNULL = FALSE, ignoreInit = TRUE, {
  
  req(length(debouncedBoxplotFactor1()) >= 1)
  req(!is.null(input$boxplotFactor1))
  req(input$boxKeepBucket[[1]] %in% inputDataReactive()$dataset[[input$boxplotFactor1]])
  
  if (inputDataReactive()$dataType == "proteomics") {
    designLevels <- input$contrastSelected %>% strsplit("_vs_") %>% .[[1]]
  }
  
  genesToPlot <- debouncedBoxplotFactor1()
  if (length(genesToPlot) > 50) {
    genesToPlot <- genesToPlot[1:50]
  }
  datasetBoxplot <- inputDataReactive()$dataset
  
  # Get the count data
  countsBoxplot <- inputDataReactive()$countList[[input$boxplotCounts]]
  # Apply the following in this specific order as required: Log2, batch correct, Z-scale
  # If the user wants to log2, then do that
  if (input$boxplotCountsLog) {
    if (inputDataReactive()$dataType == "RNASeq") {
      if (input$boxplotCounts %in% c("TPM", "FPKM", "Normalised", "Raw")) {
        countsBoxplot <- log2(countsBoxplot+1)
      }
    }
  }
  # If the user wants to remove batch effect, then do that
  if (input$boxplotBatch != "None") {
    countsBoxplot <- limma::removeBatchEffect(countsBoxplot, datasetBoxplot[[input$boxplotBatch]])
  }
  # If the user wants to plot gene-wise Z-scores, then do that 
  if (input$boxplotZScore) {
    countsBoxplotz <- t(apply(countsBoxplot, 1, zscore))
    colnames(countsBoxplotz) <- colnames(countsBoxplot)
    countsBoxplot <- countsBoxplotz
    rm(countsBoxplotz)
  }
  # Clean up proteomics feature names 
  if (inputDataReactive()$dataType == "proteomics") {
    rownames(countsBoxplot) <- gsub("\\~.*", "", rownames(countsBoxplot))
  }
  countsBoxplot <- as.data.frame(t(countsBoxplot))
  countsBoxplot$names <- rownames(countsBoxplot)
  countsBoxplot <- countsBoxplot %>% dplyr::select(all_of(c("names", genesToPlot)))
  
  # Add the metadata
  countsBoxplot <- left_join(countsBoxplot, datasetBoxplot[,c("names", inputDataReactive()$factors)], by = "names")
  
  # Include only the conditions selected
  countsBoxplot <- countsBoxplot[countsBoxplot[[input$boxplotFactor1]] %in% input$boxKeepBucket, ]
  
  # Set the factor levels per the bucket list order
  countsBoxplot[[input$boxplotFactor1]] <- factor(
    countsBoxplot[[input$boxplotFactor1]],
    input$boxKeepBucket)
  
  # melt the table
  countsBoxplotMelt <- reshape2:::melt.data.frame(data = countsBoxplot, id.vars = c("names", inputDataReactive()$factors), value.name = "Counts", variable.name = "Feature")
  countsBoxplotMelt$Feature <- as.factor(x = countsBoxplotMelt$Feature)
  
  # Force the order of the genes as factor so the order the genes are
  # selected by the user is respected by the facet wrap
  countsBoxplotMelt$Feature <- factor(countsBoxplotMelt$Feature, genesToPlot)
  
  # Add vector of custom shapes for compatibility with ggprism
  if (!input$boxplotFactor2 == "None" & !input$boxplotFactor2 == "Feature") {
    shapes <- c(rep(c(15,16,17,18,8), times = 10))[1:nlevels(as.factor(countsBoxplot[[input$boxplotFactor2]]))]
    names(shapes) <- levels(as.factor(countsBoxplot[[input$boxplotFactor2]]))
    countsBoxplotMelt$shapes <- shapes[match(as.character(countsBoxplotMelt[[input$boxplotFactor2]]), names(shapes))]
  }
  
  # add the DE results 
  countsBoxplotMelt$fdr <- seqAnnoReactive$sa$fdr[match(countsBoxplotMelt$Feature, seqAnnoReactive$sa$gene_name)]
  countsBoxplotMelt$fdr[is.na(countsBoxplotMelt$fdr)] <- 1
  countsBoxplotMelt$log2Ratio <- seqAnnoReactive$sa$log2Ratio[match(countsBoxplotMelt$Feature, seqAnnoReactive$sa$gene_name)]
  countsBoxplotMelt$log2Ratio[is.na(countsBoxplotMelt$log2Ratio)] <- 0
  countsBoxplotMelt <- add_significance(
    countsBoxplotMelt,
    p.col = "fdr",
    output.col = "stars",
    cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
    symbols = c("****", "***", "**", "*", "ns")
  )
  countsBoxplotMelt$group1 <- designLevels[1]
  countsBoxplotMelt$group2 <- designLevels[2]
  countsBoxplotMelt$y.position <- lapply(genesToPlot, function(g) {
    rep(max(countsBoxplotMelt$Counts[countsBoxplotMelt$Feature == g], na.rm = T), 
        nrow(countsBoxplotMelt[countsBoxplotMelt$Feature == g, ]))
  }) %>% unlist()
  
  
  # Get a data frame for the barplots 
  countsBarplot <- dataSummary(
    data = countsBoxplotMelt[!is.na(countsBoxplotMelt$Counts),],
    varname = "Counts",
    groupnames = c(input$boxplotFactor1, "Feature"))
  countsBarplot[[input$boxplotFactor1]] <- factor(countsBarplot[[input$boxplotFactor1]], levels = input$boxKeepBucket)
  countsBarplot$Feature <- factor(countsBarplot$Feature, genesToPlot)
  for (f in c("fdr", "log2Ratio", "stars", "group1", "group2", "y.position")) {
    countsBarplot[[f]] <- countsBoxplotMelt[match(countsBarplot$Feature, countsBoxplotMelt$Feature), f]
  }
  
  return(list(
    "datasetBoxplot" = datasetBoxplot,
    "countsBoxplot" = countsBoxplot,
    "countsBoxplotMelt" = countsBoxplotMelt,
    "countsBarplot" = countsBarplot
  ))
  
})

observeEvent({
  boxplotCountsReactive()
  input$boxplotYLimLFC
  input$textSizeBoxplot
  input$boxplotNCol
  input$boxplotVertLines
  input$boxplotGrey
  input$boxplotShowP
  input$boxplotFreeX
  input$boxplotThemeChoice
  input$boxplotShowBox
  input$boxplotShowPoint
  input$boxplotShowViolin
  input$boxplotPoints
  input$showDotsBarplot
  input$contrastSelected
  input$boxplotPointSize
  input$boxplotPointAlpha
  input$boxplotMeanLine
  input$boxplotShowMeanBar
  input$boxplotBoxAlpha
  input$boxplotPointDodge
  input$boxplotShowPLabelSize
  input$boxplotShowPBracketSize
  input$boxplotShowPTipSizeA
  input$boxplotShowPTipSizeB
  input$boxplotShowPDodge
  input$boxplotPointBorder
  input$boxplotMeanBarFront
  input$boxplotShowConditions
  input$boxplotConditionAngle
  input$boxplotConditionFormat
  input$boxplotFont
  input$barplotSDorSE
  input$boxplotShowPJust
  input$boxplotYScaleFree
  input$boxplotPlotBorder
  lapply(seq_along(inputDataReactive()$factorLevels), function (i) {
    input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]])]]
  })
}, ignoreNULL = FALSE, ignoreInit = FALSE, {
  
  req(!is.null(boxplotCountsReactive()$countsBoxplotMelt))
  datasetBoxplot <- boxplotCountsReactive()$datasetBoxplot
  countsBoxplot <- boxplotCountsReactive()$countsBoxplot
  countsBoxplotMelt <- boxplotCountsReactive()$countsBoxplotMelt
  countsBarplot <- boxplotCountsReactive()$countsBarplot
  genesToPlot <- levels(as.factor(countsBoxplotMelt$variable))
  
  # Get the colours for the condition levels
  coloursBoxplot <- setNames(lapply(input$boxKeepBucket, function(k) {
    paste0(col2hex(input[[paste0("GroupColour", input$boxplotFactor1, k)]]), "FF")
  }), input$boxKeepBucket)
  
  if (inputDataReactive()$dataType == "proteomics") {
    designLevels <- input$contrastSelected %>% strsplit("_vs_") %>% .[[1]]
  }
  
  facetFree <- "free"
  if (input$boxplotYScaleFree) {
    countsBoxplotMelt$y.position <- max(countsBoxplotMelt$y.position)
    facetFree <- "free_x"
  }
  
  # Boxplot ----
  boxplot <- ggplot(data = countsBoxplotMelt, aes(x = .data[[input$boxplotFactor1]], y = Counts))
  if (input$boxplotShowBox & !input$boxplotShowViolin) {
    boxplot <- boxplot + geom_boxplot(aes(fill = .data[[input$boxplotFactor1]]), alpha = input$boxplotBoxAlpha, outlier.shape = NA, show.legend = FALSE)
  }
  if (!input$boxplotShowBox & input$boxplotShowViolin) {
    boxplot <- boxplot + geom_violin(aes(fill = .data[[input$boxplotFactor1]]), alpha = input$boxplotBoxAlpha, draw_quantiles = 0.5, show.legend = FALSE)
  }
  if (input$boxplotShowMeanBar & !input$boxplotMeanBarFront) {
    boxplot <- boxplot + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$boxplotMeanLine, show.legend = FALSE)
  }
  if (input$boxplotShowPoint) {
    if (input$boxplotFactor2 == "None" | input$boxplotFactor2 == "Feature") {
      boxplot <- boxplot + geom_beeswarm(aes(fill = .data[[input$boxplotFactor1]]), size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, shape = 21, cex = input$boxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9, stroke = input$boxplotPointBorder)
    } else {
      boxplot <- boxplot + geom_beeswarm(aes(fill = .data[[input$boxplotFactor1]], shape = .data[[input$boxplotFactor2]]), size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, cex = input$boxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9, stroke = input$boxplotPointBorder)
      boxplot <- boxplot + scale_shape_manual(values = c(rep(c(21,22,23,24,25), times = 10))[1:nlevels(as.factor(countsBoxplotMelt[[input$boxplotFactor2]]))]) 
      boxplot <- boxplot + guides(
        fill = guide_legend(override.aes = list(shape = 21)), 
        shape = guide_legend(override.aes = list(shape = c(rep(c(16,15,18,17,25), times = 10))[1:nlevels(as.factor(countsBoxplotMelt[[input$boxplotFactor2]]))]))
      )
    }
  }
  if (input$boxplotShowMeanBar & input$boxplotMeanBarFront) {
    boxplot <- boxplot + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$boxplotMeanLine, show.legend = FALSE)
  }
  if (!input$boxplotGrey) {
    boxplot <- boxplot + scale_fill_manual(breaks = names(coloursBoxplot), values = as.character(coloursBoxplot))
  } else {
    boxplot <- boxplot + scale_colour_grey()
    boxplot <- boxplot + scale_fill_grey()
  }
  if (input$boxplotConditionFormat) {
    xlabs <- levels(countsBoxplotMelt[[input$boxplotFactor1]])
    xlabs <- gsub("\\_|\\.|\\ |\\/", "\n", xlabs)
    boxplot <- boxplot + scale_x_discrete(labels = xlabs)
  }
  ylabToPlot <- paste(input$boxplotCounts, "Counts")
  if (input$boxplotCountsLog) {
    ylabToPlot <- paste(ylabToPlot, "(Log2)")
  }
  if (input$boxplotZScore) {
    ylabToPlot <- paste(ylabToPlot, "[Z-score]")
  }
  boxplot <- boxplot + labs(y = ylabToPlot)
  boxplot <- boxplot + facet_wrap(~Feature, scales = facetFree, ncol = input$boxplotNCol)
  if (input$boxplotVertLines) {
    boxplot <- boxplot + geom_vline(xintercept = seq_along(input$boxKeepBucket)[-length(seq_along(input$boxKeepBucket))]+0.5, linetype = "dashed", alpha = 0.7)
  }
  if (input$boxplotShowP) {
    boxplot <- boxplot + stat_pvalue_manual(
      data = countsBoxplotMelt, 
      y.position = countsBoxplotMelt$y.position*1.03,
      label = "stars", 
      size = input$boxplotShowPLabelSize*2,
      label.size = input$boxplotShowPLabelSize*2,
      bracket.size = input$boxplotShowPBracketSize/6,
      tip.length = c(input$boxplotShowPTipSizeA/50, input$boxplotShowPTipSizeB/50),
      bracket.nudge.y = input$boxplotShowPJust/3
      ) + coord_cartesian(clip = "off")
  }
  boxplot <- boxplot + theme_prism(base_size = input$textSizeBoxplot, axis_text_angle = input$boxplotConditionAngle, base_family = input$boxplotFont, border = input$boxplotPlotBorder)
  if (!input$boxplotShowConditions) {
    boxplot <- boxplot + theme(axis.text.x = element_blank())
  }
  figuresDataReactive$boxplot <- boxplot
  
  
  
  # Barplot ----
  countsBarplot <- countsBarplot %>% complete(!!sym(input$boxplotFactor1), Feature, fill = list(Counts = 0))
  barplot <- ggplot(data = countsBarplot, aes(x = .data[[input$boxplotFactor1]], y = Counts, fill = .data[[input$boxplotFactor1]]))
  barplot <- barplot + facet_wrap(~Feature, scales = "free", ncol = input$boxplotNCol)
  barplot <- barplot + geom_bar(stat = "identity", color = "black", position = position_dodge(), alpha = input$boxplotBoxAlpha, show.legend = FALSE)
  barplot <- barplot + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  if (input$barplotSDorSE %in% c("SE", "SD")) {
    barplot <- barplot + geom_errorbar(aes(ymin = Counts-.data[[tolower(input$barplotSDorSE)]], ymax = Counts+.data[[tolower(input$barplotSDorSE)]]), width = .2, position = position_dodge(.9))
  }
  if (input$showDotsBarplot) {
    barplot <- barplot + geom_beeswarm(data = countsBoxplotMelt, aes_string(x = input$boxplotFactor1, y = "Counts"), pch = 21, size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, cex = 3, corral = "gutter", corral.width = 0.9)
  }
  if (input$boxplotConditionFormat) {
    barplot <- barplot + scale_x_discrete(labels = xlabs)
  }
  barplot <- barplot + labs(y = ylabToPlot)
  if (!input$boxplotGrey) {
    barplot <- barplot + scale_fill_manual(breaks = names(coloursBoxplot), values = as.character(coloursBoxplot))
  } else {
    barplot <- barplot + scale_colour_grey()
    barplot <- barplot + scale_fill_grey()
  }
  barplot <- barplot + theme_prism(base_size = input$textSizeBoxplot, axis_text_angle = input$boxplotConditionAngle, base_family = input$boxplotFont)
  if (!input$boxplotShowConditions) {
    barplot <- barplot + theme(axis.text.x = element_blank())
  }
  if (input$boxplotShowP) {
    barplot <- barplot + stat_pvalue_manual(
      data = countsBarplot, 
      y.position = countsBarplot$y.position*1.03,
      label = "stars", 
      size = input$boxplotShowPLabelSize*2,
      label.size = input$boxplotShowPLabelSize*2,
      bracket.size = input$boxplotShowPBracketSize/6,
      tip.length = c(input$boxplotShowPTipSizeA/50, input$boxplotShowPTipSizeB/50),
      bracket.nudge.y = input$boxplotShowPJust/3
    ) + coord_cartesian(clip = "off")
  }
  figuresDataReactive$barplot <- barplot
  
  # Table ----
  output$boxplotTable <- renderDataTable({
    boxplot_brushed_df <- brushedPoints(
      countsBoxplotMelt,
      input$boxplotBrush,
      xvar = input$boxplotFactor1,
      yvar = "Counts"
    )
    DT::datatable(
      data = boxplot_brushed_df[, c("names", input$boxplotFactor1, "Feature", "Counts", "log2Ratio", "fdr"), drop = FALSE],
      filter = "top", 
      rownames = F) %>%
      DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
      formatSignif(c("Counts", "fdr", "log2Ratio"), digits = 3)
  })
  
  
})

lapply(c("boxplot", "barplot"), function(b) {
  output[[paste0(b, "Static")]] <- renderPlot({
    req(!is.null(figuresDataReactive[[b]]))
    figuresDataReactive[[b]]
  },
  width = function(){as.numeric(input$figWidthBoxplot)},
  height = function(){as.numeric(input$figHeightBoxplot)}
  )
  output[[paste0(b, "DL")]] <- downloadHandler(
    filename = function() {paste0(inputDataReactive()$design, b, ".pdf")},
    content = function(file) {
      req(!is.null(figuresDataReactive[[b]]))
      pdf(file = file, width = (input$figWidthBoxplot/85), height = (input$figHeightBoxplot/90))
      print(figuresDataReactive[[b]])
      dev.off()
    }
  )
})
output$dlBoxplotButtonCounts <- downloadHandler(
  filename = function() {paste0(inputDataReactive()$design, "BoxplotTable", ".xlsx")},
  content = function(file) {
    req(!is.null(boxplotCountsReactive()$countsBoxplotMelt))
    writexl::write_xlsx(
      boxplotCountsReactive()$countsBoxplotMelt[, c("names", "Condition", "Feature", "Counts", "log2Ratio", "fdr")],
      path = file)
  }
)