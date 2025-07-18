if(inputDataReactive()$dataType == "RNASeq") {

  dataset <- inputDataReactive()$dataset
  factors <- inputDataReactive()$factors
  countList <- inputDataReactive()$countList
  factorLevels <- inputDataReactive()$factorLevels
  param <- inputDataReactive()$param
  design <- inputDataReactive()$design
  allPathways <- inputDataReactive()$allPathways
  
  output$goTilePlotDesign <- renderText({
    design
  })
  
  observeEvent(input$tabs, {
    if (input$tabs == "goTileTab") {
      if (is_empty(allPathways)) {
        showModal(modalDialog(title = "No GO Annotations", "Your genes have no GO annotations, so this tab is not functional for this dataset."))
      }
    }
  })
  
  if(!is_empty(allPathways)) {
  
    updateSelectizeInput(
      session = session,
      inputId = "goTileInput",
      choices = allPathways,
      selected = "",
      server = T
    )
    
    observeEvent(
      {
        input$figHeightGoTile
        input$figWidthGoTile
        input$textSizeGoTile
        input$goTileInput
        input$pTypeGoTile
        input$clusterRowsGoTile
        input$colnamesGoTile
        input$geneNamesGoTile
        input$GoTileColourRed
        input$GoTileColourWhite
        input$GoTileColourBlue
        input$tileLimitColour
        input$degOnlyGoTile
      },ignoreInit = TRUE,{
        
        figuresDataReactive$goTilePlot <- NULL
        if (input$pTypeGoTile == "Raw") {
          pTypeGoTile <- "pValue"
        } else {
          pTypeGoTile <- "fdr"
        }
        
        goToPlot <- gsub(" .*", "", input$goTileInput)
        dfToPlotList <- lapply(goToPlot, function(GO) {
          d1 <- inputDataReactive()$seqAnno[with(inputDataReactive()$seqAnno, grepl(GO, paste(`GO BP`, `GO MF`, `GO CC`))), ]
          d1 <- d1[d1$usedInTest == TRUE, ]
          if (input$degOnlyGoTile) {
            d1 <- d1[which(d1[[pTypeGoTile]] <= 0.05), ]
          }
          d1
        })
        names(dfToPlotList) <- goToPlot
        
        logRatioMatrixList <- lapply(goToPlot, function(GO) {
          m1 <- as.matrix(dfToPlotList[[GO]]$log2Ratio)
          colnames(m1) <- GO
          rownames(m1) <- dfToPlotList[[GO]]$gene_name
          m1
        }) %>% setNames(goToPlot)
        pValMatrixList <- lapply(goToPlot, function(GO) {
          m2 <- as.matrix(dfToPlotList[[GO]][[pTypeGoTile]])
          colnames(m2) <- GO
          rownames(m2) <- dfToPlotList[[GO]]$gene_name
          m2
        }) %>% setNames(goToPlot)
  
        cols <- colorRamp2(
          breaks = c(input$tileLimitColour, 0, -input$tileLimitColour),
          colors = c(
            paste0(input$GoTileColourRed, "FF"),
            paste0(input$GoTileColourWhite, "FF"),
            paste0(input$GoTileColourBlue, "FF")
          )
        )
  
        hmList <- lapply(goToPlot, function(GO) {
          m1 <- logRatioMatrixList[[GO]]
          m2 <- pValMatrixList[[GO]]
          
          Heatmap(
            matrix = m1,
            cell_fun = function(j, i, x, y, width, height, fill) {
              if (m2[i, j] <= 0.0001) {
                grid.text("****", x, y)
              } else if (m2[i, j] <= 0.001) {
                grid.text("***", x, y)
              } else if (m2[i, j] <= 0.01) {
                grid.text("**", x, y)
              } else if (m2[i, j] <= 0.05) {
                grid.text("*", x, y)
              }
            },
            name = "Log2FC",
            col = cols,
            cluster_rows = input$clusterRowsGoTile,
            show_column_names = input$colnamesGoTile,
            show_row_names = input$geneNamesGoTile,
            column_names_centered = TRUE,
            column_title = sub("\\s+", "\n", GO),
            column_title_gp = gpar(fontsize = ceiling(input$textSizeGoTile * 1.1), fontface = "bold"),
            row_names_gp = gpar(fontsize = input$textSizeGoTile),
            column_names_gp = gpar(fontsize = input$textSizeGoTile),
            heatmap_legend_param = list(
              title = "Log2FC", at = c(-input$tileLimitColour, 0, input$tileLimitColour),
              labels = c(-input$tileLimitColour, "0", input$tileLimitColour)
            ),
            width = ncol(m1) * unit(15, "mm"),
            height = nrow(m1) * unit(5, "mm"),
            rect_gp = gpar(col = "white", lwd = 0.5)
          )
        }) %>% setNames(goToPlot)
        
        figuresDataReactive$goTilePlot <- hmList
      }
    )
  }
} 

output$goTileOutput <- renderPlot({
  req(!is.null(figuresDataReactive$goTilePlot))
  req(length(figuresDataReactive$goTilePlot) > 0)
  
  grobList <- lapply(figuresDataReactive$goTilePlot, function(hm) {
    grid.grabExpr(draw(hm))
  })
  
  grid.arrange(grobs = grobList, ncol = length(grobList))
},
  width = function(){as.numeric(input$figWidthGoTile)},
  height = function(){as.numeric(input$figHeightGoTile)}
)

output$dlGoTilePlot <- downloadHandler(
  filename = function() {
    paste0(gsub(" .*", "", input$goTileInput), "_tile.pdf")
  },
  content = function(file) {
    req(!is.null(figuresDataReactive$goTilePlot))
    req(length(figuresDataReactive$goTilePlot) > 0)
    
    pdf(file = file,
        width = as.numeric(input$figWidthGoTile / 80),
        height = as.numeric(input$figHeightGoTile / 60))
    
    grobList <- lapply(figuresDataReactive$goTilePlot, function(hm) {
      grid.grabExpr(draw(hm))
    })
    
    grid.arrange(grobs = grobList, ncol = length(grobList))
    dev.off()
  }
)