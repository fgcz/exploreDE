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
  choices = inputDataReactive()$factors,
  selected = NULL
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
})

# Extract the appropriate data ----
if (inputDataReactive()$dataType == "RNASeq") {
  seqAnnoReactive <- reactiveValues("sa" = inputDataReactive()$seqAnno, "saF" = NULL)
  # seqAnnoReactive$saF <- seqAnnoReactive$sa[which(seqAnnoReactive$sa$usedInTest),]
  design <- inputDataReactive()$design
  param <- inputDataReactive()$param
  designLevels <- inputDataReactive()$designLevels
  output$boxplotDesign <- renderText({design})
} else if (inputDataReactive()$dataType == "proteomics") {
  seqAnnoReactive <- reactiveValues("sa" = NULL, "saF" = NULL)
  observeEvent(input$contrastSelected, {
    seqAnnoReactive$sa = inputDataReactive()$seqAnnoList[[input$contrastSelected]]
    seqAnnoReactive$sa <- seqAnnoReactive$sa %>% dplyr::select("gene_name", "log2Ratio", "pValue", "fdr")
    seqAnnoReactive$saF <- seqAnnoReactive$sa
    design <- input$contrastSelected
    output$boxplotDesign <- renderText({design})
  }, ignoreInit = T, ignoreNULL = T)
}

# The main part of the gene bucket lives here ----
# It listens in on DEG table, volcano, and heatmap tabs, in addition to this tab
genesReactive <- reactiveValues(genes = NULL)
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
  # return(list(
  #   "genes" = genesUnlist
  # ))
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
})

# Generate a reactive list of the count data to be plotted ----
# We do this so we can make aesthetic changes to the plot without having to recalculate the required boxplots and results 
boxplotCountsReactive <- eventReactive({
  input$keepBucketBoxplot
  input$boxplotFactor1
  input$boxKeepBucket
  input$boxplotFactor2
  input$boxplotBatch
  input$boxplotCounts
  input$contrastSelected
  input$boxplotShowPHeight
  input$boxplotCountsLog
}, ignoreNULL = FALSE, ignoreInit = TRUE, {
  
  req(length(input$keepBucketBoxplot) >= 1)
  req(!is.null(input$boxplotFactor1))
  req(input$boxKeepBucket[[1]] %in% inputDataReactive()$dataset[[input$boxplotFactor1]])
  
  if (inputDataReactive()$dataType == "proteomics") {
    designLevels <- input$contrastSelected %>% strsplit("_vs_") %>% .[[1]]
  }
  
  genesToPlot <- input$keepBucketBoxplot
  if (length(genesToPlot) > 50) {
    genesToPlot <- genesToPlot[1:50]
  }
  datasetBoxplot <- inputDataReactive()$dataset
  
  # Get the count data
  countsBoxplot <- inputDataReactive()$countList[[input$boxplotCounts]]
  if (input$boxplotCountsLog) {
    if (inputDataReactive()$dataType == "RNASeq") {
      if (input$boxplotCounts %in% c("TPM", "FPKM", "Normalised", "Raw")) {
        countsBoxplot <- log2(countsBoxplot+1)
      }
    }
  }
  
  if (!is.null(input$boxplotBatch)) {
    for (i in seq_along(input$boxplotBatch)) {
      countsBoxplot <- limma::removeBatchEffect(countsBoxplot, datasetBoxplot[[input$boxplotBatch[i]]])
    }
  }
  if (inputDataReactive()$dataType == "proteomics") {
    rownames(countsBoxplot) <- gsub("\\~.*", "", rownames(countsBoxplot))
  }
  countsBoxplot <- as.data.frame(t(countsBoxplot))
  countsBoxplot$names <- rownames(countsBoxplot)
  countsBoxplot <- countsBoxplot %>% dplyr::select(names, genesToPlot)
  
  # Add the metadata
  countsBoxplot <- left_join(countsBoxplot, datasetBoxplot)
  
  # Include only the conditions selected
  countsBoxplot <- countsBoxplot[countsBoxplot[[input$boxplotFactor1]] %in% input$boxKeepBucket, ]
  
  # Set the factor levels per the bucket list order
  countsBoxplot[[input$boxplotFactor1]] <- factor(
    countsBoxplot[[input$boxplotFactor1]],
    input$boxKeepBucket)
  
  # melt the table
  countsBoxplotMelt <- countsBoxplot[, colnames(countsBoxplot) %in% c(inputDataReactive()$factorNames, "names", genesToPlot)]
  countsBoxplotMelt <- reshape2::melt(data = countsBoxplotMelt)
  countsBoxplotMelt$variable <- as.factor(x = countsBoxplotMelt$variable)
  
  # Force the order of the genes as factor so the order the genes are
  # selected by the user is respected by the facet wrap
  countsBoxplotMelt$variable <- factor(countsBoxplotMelt$variable, genesToPlot)
  
  # Add vector of custom shapes for compatibility with ggprism
  if (!input$boxplotFactor2 == "None" & !input$boxplotFactor2 == "Feature") {
    shapes <- c(rep(c(15,16,17,18,8), times = 10 ))[1:nlevels(as.factor(countsBoxplot[[input$boxplotFactor2]]))]
    names(shapes) <- levels(as.factor(countsBoxplot[[input$boxplotFactor2]]))
    countsBoxplotMelt$shapes <- shapes[match(as.character(countsBoxplotMelt[[input$boxplotFactor2]]), names(shapes))]
  }
  
  # Get list of p-values as *
  cbmList <- setNames(lapply(genesToPlot, function(g) {
    x <- countsBoxplotMelt[countsBoxplotMelt$variable == g, ]
    if (inputDataReactive()$dataType == "RNASeq") {
      if (param$featureLevel == "gene") {
        x$log2Ratio <- seqAnnoReactive$sa$log2Ratio[seqAnnoReactive$sa$gene_name == g]
        x$pValue <- seqAnnoReactive$sa$pValue[seqAnnoReactive$sa$gene_name == g]
        if (seqAnnoReactive$sa$usedInTest[seqAnnoReactive$sa$gene_name == g]) {
          x$fdr <- seqAnnoReactive$sa$fdr[seqAnnoReactive$sa$gene_name == g]
        } else {
          x$fdr <- 1
        }
      } else if (param$featureLevel == "isoform") {
        x$log2Ratio <- seqAnnoReactive$saF$log2Ratio[seqAnnoReactive$saF$gene_name == g]
        x$pValue <- seqAnnoReactive$saF$pValue[seqAnnoReactive$saF$gene_name == g]
        if (seqAnnoReactive$sa$usedInTest[seqAnnoReactive$sa$gene_name == g]) {
          x$fdr <- seqAnnoReactive$sa$fdr[seqAnnoReactive$sa$gene_name == g]
        } else {
          x$fdr <- 1
        }
      }
    } else if (inputDataReactive()$dataType == "proteomics") {
      x$log2Ratio <- seqAnnoReactive$saF$log2Ratio[which(seqAnnoReactive$saF$gene_name == g)]
      x$pValue <- seqAnnoReactive$saF$pValue[which(seqAnnoReactive$saF$gene_name == g)]
      x$fdr <- seqAnnoReactive$sa$fdr[seqAnnoReactive$sa$gene_name == g]
    }
    x <- add_significance(
      x,
      p.col = "fdr",
      output.col = "stars",
      cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    )
    x$yMax <- max(x$value, na.rm = T) # *(1 + input$boxplotShowPHeight/100)
    x$xMax <- designLevels[1]
    x$xMin <- designLevels[2]
    x
  }), genesToPlot)
  
  # Get a data frame for the barplots 
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
  
  countsBoxplotMeltForBarplot <- countsBoxplotMelt[!is.na(countsBoxplotMelt$value),]
  cbSum <- dataSummary(
    data = countsBoxplotMeltForBarplot,
    varname = "value",
    groupnames = c(input$boxplotFactor1, "variable"))
  cbSum[[input$boxplotFactor1]] <- factor(cbSum[[input$boxplotFactor1]], input$boxKeepBucket)
  cbSum$variable <- factor(cbSum$variable, genesToPlot)
  cbSumList <- setNames(lapply(unique(cbSum$variable), function(gene) {
    x <- cbSum[which(cbSum$variable == gene), ]
    x$variable <- droplevels(x$variable)
    if (inputDataReactive()$dataType == "RNASeq") {
      if (param$featureLevel == "gene") {
        x$log2Ratio <- seqAnnoReactive$sa$log2Ratio[match(x$variable, seqAnnoReactive$sa$gene_name)]
        x$pValue <- seqAnnoReactive$sa$pValue[match(x$variable, seqAnnoReactive$sa$gene_name)]
        if (seqAnnoReactive$sa$usedInTest[seqAnnoReactive$sa$gene_name == gene]) {
          x$fdr <- seqAnnoReactive$sa$fdr[seqAnnoReactive$sa$gene_name == gene]
        } else {
          x$fdr <- 1
        }
      } else if (param$featureLevel == "isoform") {
        x$log2Ratio <- seqAnnoReactive$saF$log2Ratio[match(x$variable, seqAnnoReactive$sa$gene_name)]
        x$pValue <- seqAnnoReactive$saF$pValue[match(x$variable, seqAnnoReactive$sa$gene_name)]
        if (seqAnnoReactive$sa$usedInTest[seqAnnoReactive$sa$gene_name == gene]) {
          x$fdr <- seqAnnoReactive$sa$fdr[seqAnnoReactive$sa$gene_name == gene]
        } else {
          x$fdr <- 1
        }
      }
    } else if (inputDataReactive()$dataType == "proteomics") {
      x$log2Ratio <- seqAnnoReactive$saF$log2Ratio[match(x$variable, seqAnnoReactive$sa$gene_name)]
      x$pValue <- seqAnnoReactive$saF$pValue[match(x$variable, seqAnnoReactive$sa$gene_name)]
      x$fdr <- seqAnnoReactive$sa$fdr[seqAnnoReactive$sa$gene_name == gene]
    }
    x <- add_significance(
      x,
      p.col = "fdr",
      output.col = "stars",
      cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    )
    x$yMax <- max(x$value, na.rm = T) # *1.1
    x$xMax <- designLevels[1]
    x$xMin <- designLevels[2]
    x
  }), unique(cbSum$variable))
  
  return(list(
    "datasetBoxplot" = datasetBoxplot,
    "countsBoxplot" = countsBoxplot,
    "countsBoxplotMelt" = countsBoxplotMelt,
    "cbmList" = cbmList,
    "cbSum" = cbSum,
    "cbSumList" = cbSumList
  ))
  
})

themes <- list(
  "Prism" = theme_prism(),
  "Minimal" = theme_minimal(),
  "BW" = theme_bw(),
  "Classic" = theme_classic(),
  "Void" = theme_void()
)

observeEvent({
  boxplotCountsReactive()
  # input$boxplotFactor1
  # input$boxKeepBucket
  # input$keepBucketBoxplot
  input$boxplotYLimLFC
  # input$figHeightBoxplot
  # input$figWidthBoxplot
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
  lapply(seq_along(inputDataReactive()$factorLevels), function (i) {
    input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]])]]
  })
}, ignoreNULL = FALSE, ignoreInit = FALSE, {
  
  req(!is.null(boxplotCountsReactive()$countsBoxplotMelt))
  # req(length(input$keepBucketBoxplot) >= 1)
  datasetBoxplot <- boxplotCountsReactive()$datasetBoxplot
  countsBoxplot <- boxplotCountsReactive()$countsBoxplot
  countsBoxplotMelt <- boxplotCountsReactive()$countsBoxplotMelt
  cbmList <- boxplotCountsReactive()$cbmList
  cbSum <- boxplotCountsReactive()$cbSum
  cbSumList <- boxplotCountsReactive()$cbSumList
  genesToPlot <- levels(as.factor(countsBoxplotMelt$variable))
  
  # Get the colours for the condition levels
  coloursBoxplot <- setNames(lapply(input$boxKeepBucket, function(k) {
    paste0(col2hex(input[[paste0("GroupColour", input$boxplotFactor1, "__", k)]]), "FF")
  }), input$boxKeepBucket)
  
  if (inputDataReactive()$dataType == "proteomics") {
    designLevels <- input$contrastSelected %>% strsplit("_vs_") %>% .[[1]]
  }
  
  # Make normal plot with no p-value stars
  if(!input$boxplotShowP) {
    g <- ggplot(data = countsBoxplotMelt, aes_string(x = input$boxplotFactor1, y = "value"))
    if (input$boxplotShowBox & !input$boxplotShowViolin) {
      g <- g + geom_boxplot(outlier.shape = NA, aes_string(fill = input$boxplotFactor1), alpha = input$boxplotBoxAlpha)
    }
    if (input$boxplotShowViolin & !input$boxplotShowBox) {
      g <- g + geom_violin(aes_string(fill = input$boxplotFactor1), alpha = input$boxplotBoxAlpha)
    }
    if (input$boxplotShowMeanBar & !input$boxplotMeanBarFront) {
      g <- g + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$boxplotMeanLine, show.legend = FALSE)
    }
    if (input$boxplotShowPoint) {
      if (input$boxplotFactor2 == "None" | input$boxplotFactor2 == "Feature") {
        g <- g + geom_beeswarm(aes_string(fill = input$boxplotFactor1), size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, shape = 21, cex = input$boxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9, stroke = input$boxplotPointBorder)
      } else {
        g <- g + geom_beeswarm(aes_string(fill = input$boxplotFactor1, shape = input$boxplotFactor2), size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, cex = input$boxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9, stroke = input$boxplotPointBorder)
        g <- g + scale_shape_manual(values = c(rep(c(21,22,23,24,25), times = 10))[1:nlevels(as.factor(countsBoxplotMelt[[input$boxplotFactor2]]))]) 
        g <- g + guides(
          fill = guide_legend(override.aes = list(shape = 21)), 
          shape = guide_legend(override.aes = list(shape = c(rep(c(16,15,18,17,25), times = 10))[1:nlevels(as.factor(countsBoxplotMelt[[input$boxplotFactor2]]))]))
        )
      }
    }
    if (input$boxplotShowMeanBar & input$boxplotMeanBarFront) {
      g <- g + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$boxplotMeanLine, show.legend = FALSE)
    }
    if (!input$boxplotGrey) {
      g <- g + scale_fill_manual(breaks = names(coloursBoxplot), values = as.character(coloursBoxplot))
    } else {
      g <- g + scale_colour_grey()
      g <- g + scale_fill_grey()
    }
    # g <- g + themes[[input$boxplotThemeChoice]]
    g <- g + theme_prism(base_size = input$textSizeBoxplot, axis_text_angle = input$boxplotConditionAngle, base_family = input$boxplotFont)
    if (input$boxplotConditionFormat) {
      xlabs <- levels(countsBoxplotMelt[[input$boxplotFactor1]])
      xlabs <- gsub("\\_|\\.|\\ |\\/", "\n", xlabs)
      g <- g + scale_x_discrete(labels = xlabs)
    }
    if (!input$boxplotShowConditions) {
      g <- g + theme(axis.text.x = element_blank())
    }
    if (input$boxplotCountsLog) {
      if (inputDataReactive()$dataType == "RNASeq") {
        if (input$boxplotCounts %in% c("TPM", "FPKM", "Normalised", "Raw")) {
          g <- g + labs(y = paste(input$boxplotCounts, "Counts (Log2)"))
        }
      }
    } else {
      g <- g + labs(y = paste(input$boxplotCounts, "Counts"))
    }
    g <- g + facet_wrap(~variable, scales = "free", ncol = input$boxplotNCol)
    if (input$boxplotVertLines) {
      g <- g + geom_vline(xintercept = seq_along(input$boxKeepBucket)[-length(seq_along(input$boxKeepBucket))]+0.5, linetype = "dashed", alpha = 0.7)
    }
    # g <- wrap_plots(g)
    figuresDataReactive$boxplotStatic <- g
  }
  
  # Make list of plots with p-value stars
  if(input$boxplotShowP) {
    plotList <- lapply(names(cbmList), function(gene) {
      g <- ggplot(data = cbmList[[gene]], aes_string(x = input$boxplotFactor1, y = "value"))
      if (input$boxplotShowBox & !input$boxplotShowViolin) {
        g <- g + geom_boxplot(outlier.shape = NA, aes_string(fill = input$boxplotFactor1), alpha = input$boxplotBoxAlpha)
      }
      if (input$boxplotShowViolin & !input$boxplotShowBox) {
        g <- g + geom_violin(aes_string(fill = input$boxplotFactor1), alpha = input$boxplotBoxAlpha)
      }
      if (input$boxplotShowMeanBar & !input$boxplotMeanBarFront) {
        g <- g + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$boxplotMeanLine, show.legend = FALSE)
      }
      if (input$boxplotShowPoint) {
        if (input$boxplotFactor2 == "None" | input$boxplotFactor2 == "Feature") {
          g <- g + geom_beeswarm(aes_string(fill = input$boxplotFactor1), size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, shape = 21, cex = input$boxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9, stroke = input$boxplotPointBorder)
        } else {
          g <- g + geom_beeswarm(aes_string(fill = input$boxplotFactor1, shape = input$boxplotFactor2), size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, cex = input$boxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9, stroke = input$boxplotPointBorder)
          g <- g + scale_shape_manual(values = c(rep(c(21,22,23,24,25), times = 10))[1:nlevels(as.factor(cbmList[[gene]][[input$boxplotFactor2]]))]) 
          g <- g + guides(
            fill = guide_legend(override.aes = list(shape = 21)), 
            shape = guide_legend(override.aes = list(shape = c(rep(c(16,15,18,17,25), times = 10))[1:nlevels(as.factor(cbmList[[gene]][[input$boxplotFactor2]]))]))
          )
        }
      }
      if (input$boxplotShowMeanBar & input$boxplotMeanBarFront) {
        g <- g + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$boxplotMeanLine, show.legend = FALSE)
      }
      if (!isTRUE(input$boxplotGrey)) {
        g <- g + scale_fill_manual(breaks = names(coloursBoxplot), values = as.character(coloursBoxplot))
      } else {
        g <- g + scale_colour_grey()
        g <- g + scale_fill_grey()
      }
      g <- g + coord_cartesian(
        clip = "off",
        ylim = c(
          min(cbmList[[gene]][["value"]]),
          max(cbmList[[gene]][["value"]]) + input$boxplotShowPDodge/2),
        xlim = c(1, length(unique(cbmList[[gene]][[input$boxplotFactor1]])))
      )
      if (input$boxplotShowP) {
        g <- g + add_pvalue(
          cbmList[[gene]], 
          xmin = "xMin", 
          xmax = "xMax", 
          label = "stars", 
          y.position = "yMax", 
          label.size = input$boxplotShowPLabelSize*2, 
          bracket.size = input$boxplotShowPBracketSize/6, 
          tip.length = c(input$boxplotShowPTipSizeA/50,  input$boxplotShowPTipSizeB/50),
          bracket.nudge.y = input$boxplotShowPDodge/2,
          vjust = -0.01
        )
      }
      g <- g + theme_prism(base_size = input$textSizeBoxplot, axis_text_angle = input$boxplotConditionAngle, base_family = input$boxplotFont)
      if (input$boxplotConditionFormat) {
        xlabs <- levels(cbmList[[gene]][[input$boxplotFactor1]])
        xlabs <- gsub("\\_|\\.|\\ |\\/", "\n", xlabs)
        g <- g + scale_x_discrete(labels = xlabs)
      }
      if (!input$boxplotShowConditions) {
        g <- g + theme(axis.text.x = element_blank())
      }
      g <- g + theme(axis.title = element_blank())
      if (input$boxplotCountsLog) {
        if (inputDataReactive()$dataType == "RNASeq") {
          if (input$boxplotCounts %in% c("TPM", "FPKM", "Normalised", "Raw")) {
            g <- g + labs(y = paste(input$boxplotCounts, "Counts (Log2)"))
          }
        }
      } else {
        g <- g + labs(y = paste(input$boxplotCounts, "Counts"))
      }
      g <- g + ggtitle(gene)
      if (input$boxplotVertLines) {
        g <- g + geom_vline(xintercept = seq_along(input$boxKeepBucket)[-length(seq_along(input$boxKeepBucket))]+0.5, linetype = "dashed", alpha = 0.7)
      }
      g
    })
    figuresDataReactive$boxplotStatic <- plotList
  }
  
  # Make a list of barplots
  plotListB <- lapply(names(cbmList), function(gene) {
    b <- ggplot(cbSumList[[gene]], aes_string(x = input$boxplotFactor1, y = "value", fill = input$boxplotFactor1))
    if (input$showDotsBarplot) {
      b <- b + geom_bar(stat="identity", color="black", position = position_dodge(), alpha = input$boxplotBoxAlpha)
      b <- b + geom_beeswarm(data = cbmList[[gene]], aes_string(x = input$boxplotFactor1, y = "value"), pch = 21, size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, cex = 3, corral = "gutter", corral.width = 0.9)
    } else {
      b <- b + geom_bar(stat="identity", color="black", position = position_dodge(), alpha = input$boxplotBoxAlpha)
    }
    if (input$boxplotCountsLog) {
      if (inputDataReactive()$dataType == "RNASeq") {
        if (input$boxplotCounts %in% c("TPM", "FPKM", "Normalised", "Raw")) {
          b <- b + labs(y = paste(input$boxplotCounts, "Counts (Log2)"), fill = input$boxplotFactor1)
        }
      }
    } else {
      b <- b + labs(y = paste(input$boxplotCounts, "Counts"), fill = input$boxplotFactor1)
    }
    if (!isTRUE(input$boxplotGrey)) {
      barColours <- colorRampPalette(brewer.pal(8, "Set2"))(length(names(cbSumList)))
      b <- b + scale_fill_manual(breaks = names(coloursBoxplot), values = as.character(coloursBoxplot))
    } else {
      b <- b + scale_fill_grey()
    }
    b <- b + coord_cartesian(
      clip = "off",
      ylim = c(
        min(0),
        max(cbmList[[gene]][["value"]]) + input$boxplotShowPDodge/2),
      xlim = c(1, length(unique(cbSumList[[gene]][[input$boxplotFactor1]])))
    )
    if (input$boxplotShowP) {
      b <- b + add_pvalue(
        cbSumList[[gene]], 
        xmin = "xMin", 
        xmax = "xMax", 
        label = "stars", 
        y.position = "yMax", 
        label.size = input$boxplotShowPLabelSize*2, 
        bracket.size = input$boxplotShowPBracketSize/6, 
        tip.length = c(input$boxplotShowPTipSizeA/50,  input$boxplotShowPTipSizeB/50),
        bracket.nudge.y = input$boxplotShowPDodge/2,
        vjust = -0.01
      )
    }
    b <- b + geom_errorbar(aes(ymin = value-sd, ymax = value+sd), width = .2, position = position_dodge(.9))
    b <- b + ggtitle(gene)
    b <- b + theme_prism(base_size = input$textSizeBoxplot, axis_text_angle = input$boxplotConditionAngle, base_family = input$boxplotFont)
    if (input$boxplotConditionFormat) {
      xlabs <- levels(cbmList[[gene]][[input$boxplotFactor1]])
      xlabs <- gsub("\\_|\\.|\\ |\\/", "\n", xlabs)
      b <- b + scale_x_discrete(labels = xlabs)
    }
    if (!input$boxplotShowConditions) {
      b <- b + theme(axis.text.x = element_blank())
    }
    b <- b + theme(axis.title = element_blank())
    b
  })
  figuresDataReactive$barplotStatic <- plotListB
  
  # if (length(genesToPlot) >= 1) {
  output$boxplotBrushTable <- renderTable({
    brushedPoints(countsBoxplotMelt, input$boxplotBrush)
  })
  output$boxplotBrushTable <- renderDataTable({
    DT::datatable(
      data = brushedPoints(countsBoxplotMelt, input$boxplotBrush), 
      filter = "top", 
      caption = "Click and drag over dots on the boxplot plot to see those samples in this table.", 
      rownames = F) %>%
      DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
      formatSignif("value", digits = 3)
  })
  
  countsBoxplotDl <- countsBoxplot[, c(genesToPlot, "names")] %>% left_join(datasetBoxplot[,c("names", inputDataReactive()$factorNames)])
  output$dlBoxplotButtonCounts <- downloadHandler(
    filename = function() { paste0("selected_counts_", design, ".xlsx") },
    content = function(file) { openxlsx::write.xlsx(countsBoxplotDl, file = file) }
  )
})


output$boxplotStatic <- renderPlot({
  req(!is.null(figuresDataReactive$boxplotStatic))
  if (!any(class(figuresDataReactive$boxplotStatic) == "list")) {
    figuresDataReactive$boxplotStatic
  } else {
    p <- ggpubr::ggarrange(plotlist = figuresDataReactive$boxplotStatic, common.legend = TRUE, legend = "right", ncol = input$boxplotNCol, nrow = ceiling(length(figuresDataReactive$boxplotStatic)/input$boxplotNCol))
    annotate_figure(p, left = textGrob(paste(input$boxplotCounts, "Counts"), rot = 90, vjust = 1, gp = gpar(cex = input$textSizeBoxplot/10)))
  }
},
width = function(){as.numeric(input$figWidthBoxplot)},
height = function(){as.numeric(input$figHeightBoxplot)}
)
output$barplotStatic <- renderPlot({
  req(!is.null(figuresDataReactive$barplotStatic))
  p <- ggpubr::ggarrange(plotlist = figuresDataReactive$barplotStatic, common.legend = TRUE, legend = "right", ncol = input$boxplotNCol, nrow = ceiling(length(figuresDataReactive$barplotStatic)/input$boxplotNCol))
  annotate_figure(p, left = textGrob(paste(input$boxplotCounts, "Counts"), rot = 90, vjust = 1, gp = gpar(cex = input$textSizeBoxplot/10)))
}, 
  width = function(){as.numeric(input$figWidthBoxplot)},
  height = function(){as.numeric(input$figHeightBoxplot)}
)
output$dlBarplotButton <- downloadHandler(
  filename = function() {paste0(design, "_barplot.pdf")},
  content = function(file) {
    pdf(file = file, width = (as.numeric(input$figWidthBoxplot)/95), height = (as.numeric(input$figHeightBoxplot)/100), onefile = FALSE)
    p <- ggpubr::ggarrange(plotlist = figuresDataReactive$barplotStatic, common.legend = TRUE, legend = "right", ncol = input$boxplotNCol, nrow = ceiling(length(figuresDataReactive$barplotStatic)/input$boxplotNCol))
    print(annotate_figure(p, left = textGrob(paste(input$boxplotCounts, "Counts"), rot = 90, vjust = 1, gp = gpar(cex = input$textSizeBoxplot/10))))
    dev.off()
  }
)
output$dlBoxplotButton <- downloadHandler(
  filename = function() {paste0(inputDataReactive()$design, "boxplot.pdf")},
  content = function(file) {
    req(!is.null(figuresDataReactive$boxplotStatic))
    if (!any(class(figuresDataReactive$boxplotStatic) == "list")) {
      pdf(file = file, width = function(){as.numeric(input$figWidthBoxplot)/85}, height = function(){as.numeric(input$figHeightBoxplot)/90})
      print(figuresDataReactive$boxplotStatic)
      dev.off()
    } else {
      pdf(file = file, width = function(){as.numeric(input$figWidthBoxplot)/65}, height = function(){as.numeric(input$figHeightBoxplot)/70}, onefile = FALSE)
      p <- ggpubr::ggarrange(plotlist = figuresDataReactive$boxplotStatic, common.legend = TRUE, legend = "right", ncol = input$boxplotNCol, nrow = ceiling(length(plotList)/input$boxplotNCol))
      print(annotate_figure(p, left = textGrob(paste(input$boxplotCounts, "Counts"), rot = 90, vjust = 1, gp = gpar(cex = input$textSizeBoxplot/10))))
      dev.off()
    }
  }
)