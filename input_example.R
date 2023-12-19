input <- list()
names(inputDataReactive())

colourPaletteList <- list(
  "House colours" = c(
    "indianred3", "steelblue4", "chartreuse4", "grey30", 
    "goldenrod3", "firebrick4", "royalblue4", "mediumorchid3",
    "turquoise4", "darkolivegreen", "thistle4", "darkorange3", 
    "hotpink2", "burlywood3", "cadetblue4", "chocolate4", "firebrick"
  ),
  "Paired" = brewer.pal(12, "Paired"),
  "Set1" = brewer.pal(9, "Set1"),
  "Set2" = brewer.pal(8, "Set2"),
  "Set3" = brewer.pal(12, "Set3"),
  "Dark2" = brewer.pal(8, "Dark2")
)

input$colourPalette <- c("Dark2", "Paired")
for (i in seq_along(inputDataReactive()$factorLevels)) {
  input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]])]] <- rep(as.character(unlist(colourPaletteList[input$colourPalette])), times = 5)[i]
}

for (i in seq_along(1:50)) {
  input[[paste0("GroupColour", i)]] <- rep(RColorBrewer::brewer.pal(12, "Paired"), 10)[i]
}
input$selectedTable_ORA_rows_selected <- 1:5
input$textSizeORA <- 12
input$degTable_rows_selected <- c(12,68)
input$featureCalcPType <- "Raw"
input$featureCalcPValue <- 0.05
input$featureCalcLFC <- 0

input$pcaFactor1 <- inputDataReactive()$factors[1]
input$pcaGroups <- inputDataReactive()$designLevels
input$pcaCounts <- "Normalised + Log2"
input$pcaBatch <- NULL
input$pcaTopN <- 2000
input$pcaCentre = TRUE 
input$pcaScale = FALSE
input$textSizePCA <- 12
input$pcaX = "PC1"
input$pcaY = "PC2"
input$pcaFactor2 <- "None"
input$pcaAxesProp <- TRUE
input$pcaShowNames <- TRUE
input$showLinesPCA <- TRUE
input$showAxesPCA <- TRUE
input$boldPCA <- TRUE

input$boxplotFactor1 <- inputDataReactive()$factors[1]
input$keepBucketBoxplot <- inputDataReactive()$genes[1:5]
input$boxplotCounts <- names(inputDataReactive()$countList)[2]
input$boxplotBatch <- NULL
input$boxplotFactor2 <- "None"
input$boxKeepBucket <- unique(levels(as.factor(inputDataReactive()$dataset[[input$boxplotFactor1]])))
input$boxplotShowP <- FALSE
input$boxplotShowBox <- TRUE
input$boxplotShowViolin <- FALSE
input$boxplotShowMeanBar <- TRUE
input$boxplotShowPoint <- TRUE
input$boxplotPointDodge <- 2
input$boxplotBoxAlpha <- 0.2
input$boxplotPoints <- "Fancy"
input$textSizeBoxplot <- 12
input$boxplotGrey <- FALSE
input$boxplotThemeChoice <- "Prism"
input$boxplotNCol <- 3
input$boxplotVertLines <- FALSE 
input$figHeightBoxplot <- 800
input$figWidthBoxplot <- 800
input$showDotsBarplot <- TRUE
input$boxplotPointSize <- 5
input$boxplotPointAlpha <- 0.9
input$boxplotMeanLine <- 0.5
input$boxplotShowPHeight <- 1
input$boxplotShowPLabelSize <- 3
input$boxplotShowPBracketSize <- 3
input$boxplotShowPTipSizeA <- 0
input$boxplotShowPTipSizeB <- 1
input$boxplotShowPDodge <- 1

input$boxplotGenes <- inputDataReactive()$genes[1:5]
input$boxplotGenesText <- inputDataReactive()$genes[6:10]
input$boxplotFactor1 <- inputDataReactive()$factors[1]
input$boxplotFactor2 <- "None"
input$boxplotCounts <- names(inputDataReactive()$countList)[2]
input$boxplotBatch <- NULL
input$boxplotShowP <- TRUE
input$boxplotShowPoint <- TRUE
input$boxplotShowMeanBar <- TRUE
input$boxplotShowBox <- FALSE
input$boxplotShowViolin <- TRUE
input$boxplotPointSize <- 3
input$boxplotPointAlpha <- 0.9
input$boxplotBoxAlpha <- 0.2
input$boxplotMeanLine <- 0.5
# input$boxplotPointJitter
# input$boxplotMeanLineTop
# input$boxplotPoints
input$boxplotThemeChoice <- "Prism"
input$boxplotVertLines <- FALSE
input$boxplotGrey <- TRUE
input$boxplotNCol <- 3
input$textSizeBoxplot <- 12
input$figWidthBoxplot <- 600
input$figHeightBoxplot <- 900
input$showDotsBarplot <- TRUE


if (inputDataReactive()$dataType == "proteomics") {
  input$contrastSelected <- inputDataReactive()$contrasts[1]
  input$heatmapProteomicsColumnSelect <- inputDataReactive()$factors[1]
  input$downloadCount <- "transformed"
}

input$heatmapScale <- "Continuous"
input$pTypeHeatmap <- "Raw"
# input$hmKeepBucket <- unique(levels(as.factor(inputDataReactive()$dataset[,inputDataReactive()$factors[3]])))
input$hmKeepBucket <- unique(levels(as.factor(inputDataReactive()$dataset[,inputDataReactive()$factors[1]])))
input$heatmapCounts <- names(inputDataReactive()$countList)[2]

input$goTileInput <- "GO:0042254 ribosome biogenesis"
input$pTypeGoTile <- "FDR"
input$tileLimitColour <- 4
input$degOnlyGoTile <- FALSE
input$clusterRowsGoTile <- TRUE
input$colnamesGoTile <- TRUE
input$geneNamesGoTile <- TRUE
input$GoTileColourRed <- "#801717"
input$GoTileColourWhite <- "#FFFFFF"
input$GoTileColourBlue <- "#113D69"
input$textSizeGoTile <- 12
input$figWidthGoTile <- 400
input$figHeightGoTile <- 800

input$heatmapDPI <- 600
input$heatmapDownloadFormat <- "PDF"
input$heatmapFilename <- "both_regulated_heatmap"
input$heatmapLimitD <- 2
input$heatmapColourRed <- "#801717"
input$heatmapColourWhite <- "#FFFFFF"
input$heatmapColourBlue <- "#113D69"
input$pThresholdHeatmap <- 0.05
input$lfcHeatmap <- 0
input$heatmapGeneNumber <- 50
input$clusterColsHeatmap <- TRUE
input$clusterRowsHeatmap <- TRUE
input$colnamesHeatmap <- TRUE
input$geneNamesHeatmap <- TRUE
input$figHeightHeatmap <- 800
input$figWidthHeatmap <- 1000
input$textSizeHeatmap <- 12
input$heatmapBatch <- "None"
input$heatmapLimitCPalette <- "YlGnBu"
input$heatmapLimitCPaletteRev <- FALSE
input$heatmapLimitCHigh <- 1e5
input$heatmapLimitCMid <- 1e4
input$heatmapLimitCLow <- 1e3
input$heatmapFactor2 <- ""
input$showClusterColDend <- TRUE
input$showClusterRowDend <- TRUE
input$heatmapGoInput <- "GO:0042254 ribosome biogenesis"
input$heatmapGoExpr <- 0.05

input$pTypeVolcano <- "FDR"
input$lfcVolcano <- 0.25
input$pThresholdVolcano <- 0.05
input$textSizeVolcano <- 12
input$alphaVolcano <- 0.8
input$volcanoGenes <- inputDataReactive()$genes[1:5]
input$volcanoGenesText <- inputDataReactive()$genes[8:9]
input$volcanoShowGenes <- TRUE
input$showBorderVolcano <- TRUE
input$boxKeepBucketGenes <- inputDataReactive()$genes[1:5]
input$dotSizeVolcano <- 3
input$yLimVolcano <- 25
input$xLimVolcano <- 30
input$volcanoLabelAllUp <- TRUE
input$volcanoLabelAllDown <- TRUE
input$volcanoAnnotationHighlightColour <- FALSE
input$volcanoLabelMaxOverlap <- 50

input$oraType <- "BP"
input$oraDirection <- "upGenes"
input$gseaType <- "BP"
input$selectedTable_GSEA_rows_selected <- c(1,2,3,4)
input$textSizeGSEA <- 12
