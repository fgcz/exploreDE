# Update a bunch of input options ----
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
  seqAnnoReactive$saF <- seqAnnoReactive$sa[which(seqAnnoReactive$sa$usedInTest),]
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
genesReactive <- eventReactive({
  input$boxplotGenes
  input$boxplotGenesText
  input$degTable_rows_selected
  input$volcanoGenes
  input$volcanoGenesText
  input$heatmapGenes
  input$heatmapGenesText
}, ignoreNULL = FALSE, ignoreInit = TRUE, {

  genesList1 <- setNames(lapply(c("boxplotGenes", "volcanoGenes", "heatmapGenes"), function(geneList) {
    if(!is.null(input[[geneList]])) {
      input[[geneList]]
    }
  }), c("boxplotGenes", "volcanoGenes", "heatmapGenes"))

  genesList2 <- setNames(lapply(c(
    "boxplotGenesText", "volcanoGenesText", "heatmapGenesText"
  ), function(geneText) {
    if(!is.null(input[[geneText]])) {
      strsplit(input[[geneText]], "\\ +|\\/+|\\,+|\n+|\r+|\n\r+|\t+")[[1]]
    }
  }), c("boxplotGenesText", "volcanoGenesText", "heatmapGenesText"))

  degMATable <- seqAnnoReactive$sa[order(seqAnnoReactive$sa$pValue), ]
  degMATable <- degMATable$gene_name[input$degTable_rows_selected]
  genesUnlist <- unique(unlist(c(genesList1, genesList2, degMATable)))
  genesUnlist <- genesUnlist[!genesUnlist == ""]
  return(list(
    "genes" = genesUnlist
  ))
})

observeEvent({genesReactive()}, {
  output$geneBucket1 <- renderUI({
    bucket_list(
      header = "Drag and drop features in order to be plotted",
      group_name = "bucket_list_group",
      orientation = "horizontal",
      add_rank_list(
        text = "Include these features in this order",
        labels = genesReactive()$genes,
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
  }, ignoreNULL = FALSE, ignoreInit = FALSE, {
    
    req(nchar(input$keepBucketBoxplot) >= 2)
    
    if (inputDataReactive()$dataType == "proteomics") {
      designLevels <- input$contrastSelected %>% strsplit("_vs_") %>% .[[1]]
    }
    
    genesToPlot <- input$keepBucketBoxplot
    datasetBoxplot <- inputDataReactive()$dataset

    # Get the count data
    countsBoxplot <- inputDataReactive()$countList[[input$boxplotCounts]]
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
          x$fdr <- seqAnnoReactive$sa$fdr[seqAnnoReactive$sa$gene_name == g]
        } else if (param$featureLevel == "isoform") {
          x$log2Ratio <- seqAnnoReactive$saF$log2Ratio[seqAnnoReactive$saF$gene_name == g]
          x$pValue <- seqAnnoReactive$saF$pValue[seqAnnoReactive$saF$gene_name == g]
          x$fdr <- seqAnnoReactive$saF$fdr[seqAnnoReactive$saF$gene_name == g]
        }
      } else if (inputDataReactive()$dataType == "proteomics") {
        x$log2Ratio <- seqAnnoReactive$saF$log2Ratio[which(seqAnnoReactive$saF$gene_name == g)]
        x$pValue <- seqAnnoReactive$saF$pValue[which(seqAnnoReactive$saF$gene_name == g)]
        x$fdr <- seqAnnoReactive$saF$fdr[which(seqAnnoReactive$saF$gene_name == g)]
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
          x$fdr <- seqAnnoReactive$sa$fdr[match(x$variable, seqAnnoReactive$sa$gene_name)]
        } else if (param$featureLevel == "isoform") {
          x$log2Ratio <- seqAnnoReactive$saF$log2Ratio[match(x$variable, seqAnnoReactive$sa$gene_name)]
          x$pValue <- seqAnnoReactive$saF$pValue[match(x$variable, seqAnnoReactive$sa$gene_name)]
          x$fdr <- seqAnnoReactive$saF$fdr[match(x$variable, seqAnnoReactive$sa$gene_name)]
        }
      } else if (inputDataReactive()$dataType == "proteomics") {
        x$log2Ratio <- seqAnnoReactive$saF$log2Ratio[match(x$variable, seqAnnoReactive$sa$gene_name)]
        x$pValue <- seqAnnoReactive$saF$pValue[match(x$variable, seqAnnoReactive$sa$gene_name)]
        x$fdr <- seqAnnoReactive$saF$fdr[match(x$variable, seqAnnoReactive$sa$gene_name)]
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
  # boxplotCountsReactive()
  input$boxplotFactor1
  input$boxKeepBucket
  input$boxplotYLimLFC
  input$figHeightBoxplot
  input$figWidthBoxplot
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
  lapply(seq_along(inputDataReactive()$factorLevels), function (i) {
    input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]])]]
  })
  input$keepBucketBoxplot
}, ignoreNULL = FALSE, ignoreInit = TRUE, {

  req(nchar(input$keepBucketBoxplot) >= 2)
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
    if (input$boxplotShowMeanBar) {
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
    if (!input$boxplotGrey) {
      g <- g + scale_fill_manual(breaks = names(coloursBoxplot), values = as.character(coloursBoxplot))
    } else {
      g <- g + scale_colour_grey()
      g <- g + scale_fill_grey()
    }
    g <- g + themes[[input$boxplotThemeChoice]]
    g <- g + labs(y = paste(input$boxplotCounts, "Counts"))
    g <- g + theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour = "grey20", size = input$textSizeBoxplot, angle = 0, hjust = 1, vjust = 0.5, face = "bold"),
      axis.title.x = element_text(colour = "grey20", size = input$textSizeBoxplot, angle = 0, hjust = .5, vjust = 0, face = "bold"),
      axis.title.y = element_text(colour = "grey20", size = input$textSizeBoxplot, angle = 90, hjust = .5, vjust = .5, face = "bold"),
      legend.text = element_text(colour = "black", size = input$textSizeBoxplot),
      legend.title = element_text(colour = "black", size = input$textSizeBoxplot),
      title = element_text(colour = "black", size = input$textSizeBoxplot),
      strip.text = element_text(size = input$textSizeBoxplot, face = "bold")
    )
    g <- g + facet_wrap(~variable, scales = "free", ncol = input$boxplotNCol)
    if (input$boxplotVertLines) {
      g <- g + geom_vline(xintercept = seq_along(input$boxKeepBucket)[-length(seq_along(input$boxKeepBucket))]+0.5, linetype = "dashed", alpha = 0.7)
    }
    
    output$boxplotStatic <- renderPlot({
      g
    }, width = as.numeric(input$figWidthBoxplot), height = as.numeric(input$figHeightBoxplot))
    output$dlBoxplotButton <- downloadHandler(
      filename = function() {paste0(design, "boxplot.pdf")},
      content = function(file) {
        pdf(file = file, width = (as.numeric(input$figWidthBoxplot)/85), height = (as.numeric(input$figHeightBoxplot)/90), )
        print(g)
        dev.off()
      }
    )
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
      if (input$boxplotShowMeanBar) {
        g <- g + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$boxplotMeanLine, show.legend = FALSE)
      }
      if (input$boxplotShowPoint) {
        if (input$boxplotFactor2 == "None" | input$boxplotFactor2 == "Feature") {
          g <- g + geom_beeswarm(aes_string(fill = input$boxplotFactor1), size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, shape = 21, cex = input$boxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9)
        } else {
          g <- g + geom_beeswarm(aes_string(fill = input$boxplotFactor1, shape = input$boxplotFactor2), size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, cex = input$boxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9)
          g <- g + scale_shape_manual(values = c(rep(c(21,22,23,24,25), times = 10))[1:nlevels(as.factor(cbmList[[gene]][[input$boxplotFactor2]]))]) 
          g <- g + guides(
            fill = guide_legend(override.aes = list(shape = 21)), 
            shape = guide_legend(override.aes = list(shape = c(rep(c(16,15,18,17,25), times = 10))[1:nlevels(as.factor(cbmList[[gene]][[input$boxplotFactor2]]))]))
          )
        }
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
      g <- g + themes[[input$boxplotThemeChoice]]
      g <- g + labs(y = paste(input$boxplotCounts, "Counts"))
      g <- g + ggtitle(gene)
      g <- g + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "grey20", size = input$textSizeBoxplot, angle = 0, hjust = 1, vjust = 0.5, face = "bold"),
        axis.title.y = element_blank(),
        legend.text = element_text(colour = "black", size = input$textSizeBoxplot),
        legend.title = element_text(colour = "black", size = input$textSizeBoxplot),
        title = element_text(colour = "black", size = input$textSizeBoxplot)
      )
      if (input$boxplotVertLines) {
        g <- g + geom_vline(xintercept = seq_along(input$boxKeepBucket)[-length(seq_along(input$boxKeepBucket))]+0.5, linetype = "dashed", alpha = 0.7)
      }
      g
    })
    
    output$boxplotStatic <- renderPlot({
      p <- ggpubr::ggarrange(plotlist = plotList, common.legend = TRUE, legend = "right", ncol = input$boxplotNCol, nrow = ceiling(length(plotList)/input$boxplotNCol))
      annotate_figure(p, left = textGrob(paste(input$boxplotCounts, "Counts"), rot = 90, vjust = 1, gp = gpar(cex = input$textSizeBoxplot/10)))
    }, width = as.numeric(input$figWidthBoxplot), height = as.numeric(input$figHeightBoxplot))
    output$dlBoxplotButton <- downloadHandler(
      filename = function() {paste0(design, "boxplot.pdf")},
      content = function(file) {
        pdf(file = file, width = (as.numeric(input$figWidthBoxplot)/65), height = (as.numeric(input$figHeightBoxplot)/70), onefile = FALSE)
        p <- ggpubr::ggarrange(plotlist = plotList, common.legend = TRUE, legend = "right", ncol = input$boxplotNCol, nrow = ceiling(length(plotList)/input$boxplotNCol))
        print(annotate_figure(p, left = textGrob(paste(input$boxplotCounts, "Counts"), rot = 90, vjust = 1, gp = gpar(cex = input$textSizeBoxplot/10))))
        dev.off()
      }
    )
  }

  # Make a list of barplots
  plotListB <- lapply(names(cbmList), function(gene) {
    b <- ggplot(cbSumList[[gene]], aes_string(x = input$boxplotFactor1, y = "value", fill = input$boxplotFactor1))
    if (input$showDotsBarplot) {
      b <- b + geom_bar(stat="identity", color="black", position = position_dodge(), alpha = input$boxplotBoxAlpha)
      b <- b + geom_jitter(data = cbmList[[gene]], aes_string(x = input$boxplotFactor1, y = "value"), pch = 21, size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, width = 0.3)
    } else {
      b <- b + geom_bar(stat="identity", color="black", position = position_dodge(), alpha = input$boxplotBoxAlpha)
    }
    b <- b + labs(y = paste(input$boxplotCounts, "Counts"), fill = input$boxplotFactor1)
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
    b <- b + themes[[input$boxplotThemeChoice]]
    b <- b + ggtitle(gene)
    b <- b + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour = "black", size = input$textSizeBoxplot, angle = 0, hjust = 1, vjust = 0.5, face = "bold"),
      legend.text = element_text(colour = "black", size = input$textSizeBoxplot),
      legend.title = element_text(colour = "black", size = input$textSizeBoxplot),
      title = element_text(colour = "black", size = input$textSizeBoxplot)
    )
    b
  })
  output$barplotStatic <- renderPlot({
    p <- ggpubr::ggarrange(plotlist = plotListB, common.legend = TRUE, legend = "right", ncol = input$boxplotNCol, nrow = ceiling(length(plotListB)/input$boxplotNCol))
    annotate_figure(p, left = textGrob(paste(input$boxplotCounts, "Counts"), rot = 90, vjust = 1, gp = gpar(cex = input$textSizeBoxplot/10)))
  }, width = as.numeric(input$figWidthBoxplot), height = as.numeric(input$figHeightBoxplot))
  output$dlBarplotButton <- downloadHandler(
    filename = function() {paste0(design, "_barplot.pdf")},
    content = function(file) {
      pdf(file = file, width = (as.numeric(input$figWidthBoxplot)/95), height = (as.numeric(input$figHeightBoxplot)/100), onefile = FALSE)
      p <- ggpubr::ggarrange(plotlist = plotListB, common.legend = TRUE, legend = "right", ncol = input$boxplotNCol, nrow = ceiling(length(plotListB)/input$boxplotNCol))
      print(annotate_figure(p, left = textGrob(paste(input$boxplotCounts, "Counts"), rot = 90, vjust = 1, gp = gpar(cex = input$textSizeBoxplot/10))))
      dev.off()
    }
  )

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
