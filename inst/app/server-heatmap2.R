req(!is.null(inputDataReactive()$dataType))

if (inputDataReactive()$dataType == "proteomics") {
  output$heatmapProteomicsColumnSelectUI <- renderUI({
    selectInput(inputId = "heatmapProteomicsColumnSelect", label = "Select main factor to plot by", choices = inputDataReactive()$factors, selected = inputDataReactive()$factors[1], multiple = F)
  })
  groupingNameHeatmap <- eventReactive(input$heatmapProteomicsColumnSelect, {
    return(list(
      "gn" = input$heatmapProteomicsColumnSelect
    ))
  })
} else if (inputDataReactive()$dataType == "RNASeq") {
  groupingNameHeatmap <- reactive(return(list("gn" = inputDataReactive()$param$groupingName)))
}

if (inputDataReactive()$dataType == "proteomics") {
  output$clusterWarningProteomics <- renderUI({
    helpText(
      "Clustering by row/column may break the heatmap for proteomics data, due 
      to the high number of NA values in the count table.")
  })
  updateCheckboxInput(session = session, inputId = "clusterColsHeatmap", value = FALSE)
  updateCheckboxInput(session = session, inputId = "clusterRowsHeatmap", value = FALSE)
  updateCheckboxInput(session = session, inputId = "heatmapLog2", value = FALSE)
}

updateSelectInput(
  session = session, 
  inputId = "heatmapCounts", 
  choices = names(inputDataReactive()$countList), 
  selected = switch(inputDataReactive()$dataType, "RNASeq" = "Normalised + Log2", "proteomics" = "transformedData")
)

updateSelectizeInput(
  session = session,
  inputId = "heatmapGenes",
  choices = inputDataReactive()$genes,
  selected = NULL,
  server = T, 
)

observeEvent({genesReactive}, {
  output$geneBucket3 <- renderUI({
    bucket_list(
      header = "Drag and drop features in order to be plotted",
      group_name = "bucket_list_group",
      orientation = "horizontal",
      add_rank_list(
        text = "Include these features in this order",
        labels = genesReactive$genes,
        input_id = "keepBucketHeatmap"),
      add_rank_list(
        text = "Exclude these features",
        labels = NULL,
        input_id = "excludeBucketHeatmap")
    )
  })
})

observeEvent(input$resetGeneBucketHeatmap, {
  genesReactive$genes <- NULL
  updateSelectizeInput(session, inputId = "heatmapGenes", choices = inputDataReactive()$genes, selected = NULL, server = T)
  updateTextAreaInput(session, inputId = "heatmapGenesText", value = NULL)
  output$heatmapCustom <- NULL
}, ignoreNULL = T, ignoreInit = T)

observeEvent(groupingNameHeatmap(), {
  output$heatmapBucket <- renderUI({
    bucket_list(
      header = "Drag and drop groups in order to be plotted",
      group_name = "hmBucket",
      orientation = "horizontal",
      add_rank_list(
        text = "Include these groups",
        labels = as.list(unique(levels(as.factor(inputDataReactive()$dataset[[groupingNameHeatmap()$gn]])))),
        input_id = "hmKeepBucket"),
      add_rank_list(
        text = "Exclude these groups",
        labels = NULL,
        input_id = "hmExcludeBucket")
    )
  })
  updateCheckboxGroupInput(
    session = session,
    inputId = "heatmapFactors",
    choices = inputDataReactive()$factors,
    selected = groupingNameHeatmap()$gn
  )
})

if (inputDataReactive()$dataType == "proteomics") {
  observeEvent(
    input$heatmapProteomicsColumnSelect,{
      updateCheckboxGroupInput(session = session, inputId = "heatmapFactors", selected = input$heatmapProteomicsColumnSelect)
    })
}

observeEvent(inputDataReactive(), {
  if (inputDataReactive()$dataType == "RNASeq") {
    req(!is.null(inputDataReactive()$allPathways))
    output$heatmapGOUI1 <- renderUI({
      selectInput(
        inputId = "heatmapGoInput",
        choices = inputDataReactive()$allPathways,
        selected = "",
        label = "Or, select the GO Term you wish to plot",
        selectize = TRUE
      )
    })
    output$heatmapGOUI2 <- renderUI({
      numericInput(
        inputId = "heatmapGoExpr",
        label = "Minimum median feature expression for GO plot",
        value = 0.05,
        min = 0,
        max = 1e100
      )
    })
  } else {
    output$heatmapGOUI1 <- renderUI({
      h5("Proteomics data doesn't support GO terms.")
    })
  }
})

updateSelectInput(session = session, inputId = "heatmapBatch", choices = c("None", inputDataReactive()$factors), selected = "None")

if (inputDataReactive()$dataType == "RNASeq") {
  output$heatmapGOUI <- renderUI({
    helpText("You can also plot the expression of features annotated to a given GO term (again, irrespective of p-value from the DE test).")
    selectizeInput(
      inputId = "heatmapGoInput",
      choices = c(),
      label = "Or, select the GO Term you wish to plot"
    )
    numericInput(
      inputId = "heatmapGoExpr",
      label = "Minimum median feature expression for GO plot",
      value = 0.05,
      min = 0,
      max = 1e100
    )
  })
}

observeEvent(inputDataReactive(), {
  if (inputDataReactive()$dataType == "proteomics") {
    updateSelectInput(session = session, inputId = "pTypeHeatmap", selected = "Raw")
  } else {
    updateSelectInput(session = session, inputId = "pTypeHeatmap", selected = "FDR")
  }
})

if (inputDataReactive()$dataType == "proteomics") {
  req(!is.null(inputDataReactive()$seqAnnoList))
  if ("nrPeptides" %in% colnames(inputDataReactive()$seqAnnoList[[1]])) {
    output$nrPeptidesHeatmapUI <- renderUI({
      sliderInput(inputId = "nrPeptideHeatmap", label = "Minimum number of peptides", min = 0, max = 20, value = 2, step = 1, width = "85%")
    })
  } else {
    output$nrPeptidesHeatmapUI <- renderUI({ NULL })
  }
}

observeEvent(input$heatmapZScore, {
  if (input$heatmapZScore == "Z-Score") {
    updateNumericInput(session = session, inputId = "heatmapAtLow", value = -2, min = -10, max = 10, step = 0.1)
    updateNumericInput(session = session, inputId = "heatmapAtMid", value = 0, min = -10, max = 10, step = 0.1)
    updateNumericInput(session = session, inputId = "heatmapAtHigh", value = 2, min = -10, max = 10, step = 0.1)
  } else if (input$heatmapZScore == "Centred") {
    updateNumericInput(session = session, inputId = "heatmapAtLow", value = -1, min = -10, max = 10, step = 0.1)
    updateNumericInput(session = session, inputId = "heatmapAtMid", value = 0, min = -10, max = 10, step = 0.1)
    updateNumericInput(session = session, inputId = "heatmapAtHigh", value = 1, min = -10, max = 10, step = 0.1)
  } else if (input$heatmapLog2 & input$heatmapZScore == "None") {
    updateNumericInput(session = session, inputId = "heatmapAtLow", value = 0, min = -10, max = 10, step = 0.1)
    updateNumericInput(session = session, inputId = "heatmapAtMid", value = 6, min = -10, max = 10, step = 0.1)
    updateNumericInput(session = session, inputId = "heatmapAtHigh", value = 12, min = -10, max = 10, step = 0.1)
  } else if (!input$heatmapLog2 & input$heatmapZScore == "None") {
    updateNumericInput(session = session, inputId = "heatmapAtLow", value = 0, min = 0, max = 1e6, step = 100)
    updateNumericInput(session = session, inputId = "heatmapAtMid", value = 1000, min = 0, max = 1e6, step = 100)
    updateNumericInput(session = session, inputId = "heatmapAtHigh", value = 10000, min = 0, max = 1e6, step = 100)
  }
})

observeEvent(
  {
    input$textSizeHeatmap
    input$heatmapBatch
    input$lfcHeatmap
    input$pTypeHeatmap
    input$pThresholdHeatmap
    input$heatmapLimitD
    input$heatmapLimitCHigh
    input$heatmapLimitCMid
    input$heatmapLimitCLow
    input$heatmapLimitCPalette
    input$heatmapLimitCPaletteRev
    input$heatmapGeneNumber
    input$hmKeepBucket
    input$heatmapFactors
    input$heatmapCounts
    input$clusterRowsHeatmap
    input$clusterColsHeatmap
    input$colnamesHeatmap
    input$geneNamesHeatmap
    input$heatmapColourRed
    input$heatmapColourWhite
    input$heatmapColourBlue
    input$showClusterColDend
    input$showClusterRowDend
    input$heatmapDPI
    input$heatmapDownloadFormat
    input$heatmapFilename
    input$nrPeptideHeatmap
    input$heatmapCustomTitle
    input$heatmapFeatureDirection
    input$heatmapLog2
    input$heatmapZScore
    input$heatmapAtLow
    input$heatmapAtMid
    input$heatmapAtHigh
    lapply(seq_along(inputDataReactive()$factorLevels), function (i) {
      input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]])]]
    })
    input$contrastSelected
    groupingNameHeatmap()
  }, ignoreNULL = FALSE, ignoreInit = TRUE, {
    tryCatch({
      if (!is.null(input$hmKeepBucket)) {
        
        if (inputDataReactive()$dataType == "proteomics") {
          design <- input$contrastSelected
          seqAnnoHeatmap <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]
          seqAnnoHeatmap[["usedInTest"]] = TRUE
          if (!is.null(input$nrPeptideHeatmap)) {
            seqAnnoHeatmap <- seqAnnoHeatmap[which(seqAnnoHeatmap$nrPeptides >= input$nrPeptideHeatmap),]
          }
        } else {
          seqAnnoHeatmap <- inputDataReactive()$seqAnno
          param <- inputDataReactive()$param
          design <- inputDataReactive()$design
        }
        
        # Get the desired p-value column
        if (input$pTypeHeatmap == "Raw") {
          pTypeHeatmap <- "pValue"
        } else {
          pTypeHeatmap <- "fdr"
        }
        
        # get the metadata for plotting
        datasetHeatmap <- inputDataReactive()$dataset
        datasetHeatmap <- datasetHeatmap[which(datasetHeatmap[[groupingNameHeatmap()$gn]] %in% input$hmKeepBucket), ]
        
        # Set the factor levels per the bucket list order
        datasetHeatmap[[groupingNameHeatmap()$gn]] <- factor(datasetHeatmap[[groupingNameHeatmap()$gn]], input$hmKeepBucket)
        datasetHeatmap <- datasetHeatmap[order(datasetHeatmap[[groupingNameHeatmap()$gn]]), ]
        
        # Subset count matrix:
        countsHeatmap <- inputDataReactive()$countList[[input$heatmapCounts]]
        countsHeatmap <- countsHeatmap[, rownames(datasetHeatmap)]
        
        # Apply the following in this specific order as required: Log2, batch correct, Z-scale
        # If the user wants to log2, then do that
        if (input$heatmapLog2 & input$heatmapCounts %in% c("FPKM", "TPM", "Normalised")) {
          countsHeatmap <- log2(countsHeatmap + 1)
        }
        # If the user wants to remove batch effect, then do that
        if (input$heatmapBatch != "None") {
          countsHeatmap <- limma::removeBatchEffect(countsHeatmap, datasetHeatmap[[input$heatmapBatch]])
        }
        # If the user wants to plot gene-wise Z-scores, then do that 
        if (input$heatmapZScore == "Z-Score") {
          countsHeatmapz <- t(apply(countsHeatmap, 1, zscore))
          colnames(countsHeatmapz) <- colnames(countsHeatmap)
          countsHeatmap <- countsHeatmapz
          rm(countsHeatmapz)
        }
        # If the user wants to centre the counts, then do that
        if (input$heatmapZScore == "Centred") {
          countsHeatmap <- sweep(countsHeatmap, 1, rowMeans(countsHeatmap, na.rm = T))
        }
        # Fix the proteomics rownames
        if (inputDataReactive()$dataType == "proteomics") {
          rownames(countsHeatmap) <- gsub("\\~.*", "", rownames(countsHeatmap))
        }
        countsHeatmap <- countsHeatmap[!is.na(rownames(countsHeatmap)), ]
        
        # Get the heatmap colour palette to use
        at <- c(input$heatmapAtLow, input$heatmapAtMid, input$heatmapAtHigh)
        if (input$heatmapLimitCPalette == "None") {
          heatmapColours <- colorRamp2(
            breaks = rev(at),
            colors = c(input$heatmapColourRed, input$heatmapColourWhite, input$heatmapColourBlue)
          )
        } else {
          if (input$heatmapLimitCPaletteRev) {
            heatmapColours <- colorRampPalette(brewer.pal(9, input$heatmapLimitCPalette))(64)
          } else {
            heatmapColours <- rev(colorRampPalette(brewer.pal(9, input$heatmapLimitCPalette))(64))
          }
        }
        
        # Top Annotation
        colourListHeatmap <- setNames(lapply(inputDataReactive()$factors, function(f) {
          setNames(lapply(unique(datasetHeatmap[[f]]), function(k) {
            paste0(col2hex(input[[paste0("GroupColour", f, "__", k)]]), "FF")
          }), unique(datasetHeatmap[[f]]))
        }), inputDataReactive()$factors)
        
        if (any(input$heatmapFactors != "")) {
          p1 <- paste(lapply(input$heatmapFactors, function(g) {
            paste0("'", g, "' = datasetHeatmap[['", g, "']]")
          }), collapse = ",")
          ha <- paste0(
            "HeatmapAnnotation("
            , p1, ", show_legend = T, annotation_legend_param = list(title_gp = gpar(fontsize = ", input$textSizeHeatmap, ", fontface = 'bold'), labels_gp = gpar(fontsize = ", input$textSizeHeatmap, ")), gp = gpar(col = 'black', lwd = 0.5))"
            )
          ha <- eval(parse(text = ha))
          for (f in c(input$heatmapFactors)) {
            for (l in levels(as.factor(datasetHeatmap[, f]))) {
              ha@anno_list[[f]]@color_mapping@full_col[[l]] <- colourListHeatmap[[f]][[l]]
              ha@anno_list[[f]]@fun@var_env$color_mapping@full_col[[l]] <- colourListHeatmap[[f]][[l]]
            }
          }
        }
        
        # Draw the heatmaps
        output[["heatmapDesignBoth Directions"]] <- renderText({design})
        # Get the features only in the specified DE direction
        geneDirections <- list(
          "Both" = seqAnnoHeatmap$gene_name[seqAnnoHeatmap[[pTypeHeatmap]] < input$pThresholdHeatmap & abs(seqAnnoHeatmap$log2Ratio) >= input$lfcHeatmap],
          "Up-regulated" = seqAnnoHeatmap$gene_name[seqAnnoHeatmap[[pTypeHeatmap]] < input$pThresholdHeatmap & seqAnnoHeatmap$log2Ratio >= input$lfcHeatmap],
          "Down-regulated" = seqAnnoHeatmap$gene_name[seqAnnoHeatmap[[pTypeHeatmap]] < input$pThresholdHeatmap & seqAnnoHeatmap$log2Ratio <= -input$lfcHeatmap]
        )
        countsHeatmapUpDown <- countsHeatmap[match(geneDirections[[input$heatmapFeatureDirection]], rownames(countsHeatmap)), ]
        countsHeatmapSig <- countsHeatmapUpDown
        if (as.numeric(input$heatmapGeneNumber) <= nrow(countsHeatmapSig)) {
          countsHeatmapSig <- countsHeatmapSig[1:as.numeric(input$heatmapGeneNumber), ]
        } 
        rownames(countsHeatmapSig) <- substr(rownames(countsHeatmapSig), 1, 15)
        
        hmName <- paste(input$heatmapCounts %>% gsub("\\+", "\n\\+", .))
        if (input$heatmapLog2) {
          hmName <- paste(hmName, "\n(Log2)")
        }
        if (input$heatmapZScore == "Z-Score") {
          hmName <- paste(hmName, "\n[Z Score]")
        } else if (input$heatmapZScore == "Centred") {
          hmName <- paste(hmName, "\n[Centred]")
        }
        
        heatmap <- ComplexHeatmap::Heatmap(
          matrix = countsHeatmapSig,
          top_annotation = ha,
          name = hmName,
          cluster_columns = input$clusterColsHeatmap,
          cluster_rows = input$clusterRowsHeatmap,
          show_column_names = input$colnamesHeatmap,
          show_row_names = input$geneNamesHeatmap,
          show_column_dend = input$showClusterColDend, 
          show_row_dend = input$showClusterRowDend, 
          rect_gp = gpar(col = "white", lwd = 0.5),
          use_raster = FALSE,
          row_dend_width = unit(1, "cm"),
          column_dend_height = unit(1, "cm"),
          col = heatmapColours, 
          row_names_gp = gpar(fontsize = input$textSizeHeatmap),
          column_names_gp = gpar(fontsize = input$textSizeHeatmap), 
          column_title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
          column_title = paste("Top", nrow(countsHeatmapSig), "DE Features:\n", design),
          heatmap_legend_param = list(labels_gp = gpar(fontsize = input$textSizeHeatmap), title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"))
        ) %>% draw(merge_legend = TRUE)
        figuresDataReactive[["heatmapDEFeatures"]] <- heatmap
        
        # Download the count file
        countsHeatmapSigDL <- as.data.frame(countsHeatmapSig) %>% rownames_to_column("gene_name")
        output[["dlHeatmapDFButtonDEFeatures"]] <- downloadHandler(
          filename = function() { paste0(design, "_DE_", "_heatmap_counts.xlsx") },
          content = function(file) { openxlsx::write.xlsx(countsHeatmapSigDL, file = file) }
        )
        # Show number of genes on plot (changes depending on settings applied)
        output[["hGNDE"]] <- renderText({
          paste0("Number of features per settings chosen: ", nrow(countsHeatmapSig))
        })
        
        # Make a custom heatmap from the genes selected
        observeEvent({
          input$keepBucketHeatmap
        }, {
          hmCustomMatrix <- countsHeatmap[input$keepBucketHeatmap, ]
          ch <- ComplexHeatmap::Heatmap(
            matrix = hmCustomMatrix,
            top_annotation = ha,
            name = input$heatmapCounts,
            cluster_columns = input$clusterColsHeatmap,
            cluster_rows = input$clusterRowsHeatmap,
            show_column_names = input$colnamesHeatmap,
            show_row_names = input$geneNamesHeatmap,
            col = heatmapColours,
            column_title = paste(input$heatmapCustomTitle),
            rect_gp = gpar(col = "white", lwd = 0.5)
          ) %>% draw(merge_legend = TRUE)
          figuresDataReactive$heatmapCustom <- ch
          
          output[[paste0("dlHeatmapDFButtonCustom")]] <- downloadHandler(
            filename = function() {
              paste0(design, "_Custom_heatmap.xlsx")
            },
            content = function(file) {
              openxlsx::write.xlsx((as.data.frame(hmCustomMatrix) %>% rownames_to_column("gene_name")), file = file)
            }
          )
        })
        
        # heatmapGoInput
        if (inputDataReactive()$dataType == "RNASeq") {
          observeEvent({
            input$heatmapGoInput
            input$heatmapGoExpr
          }, ignoreInit = T, {
            goToPlot <- gsub(" .*", "", input$heatmapGoInput)
            goGenes <- seqAnnoHeatmap$gene_name[with(seqAnnoHeatmap, grepl(goToPlot, paste(`GO BP`, `GO MF`, `GO CC`)))]
            
            if(!is.null(goGenes)) {
              hmGOMatrix <- countsHeatmap[rownames(countsHeatmap) %in% goGenes, ] %>% as.matrix()
              hmGOMatrix <- hmGOMatrix[abs(rowMedians(hmGOMatrix)) >= input$heatmapGoExpr, ]
              
              # Output heatmaps
              gh <- ComplexHeatmap::Heatmap(
                matrix = hmGOMatrix,
                top_annotation = ha,
                name = input$heatmapCounts,
                cluster_columns = input$clusterColsHeatmap,
                cluster_rows = input$clusterRowsHeatmap,
                show_column_names = input$colnamesHeatmap,
                show_row_names = input$geneNamesHeatmap,
                col = heatmapColours,
                column_title = sub("\\s+", "\n", input$heatmapGoInput),
                rect_gp = gpar(col = "white", lwd = 0.5)
              ) %>% draw(merge_legend = TRUE)
              figuresDataReactive$heatmapGO <- gh
              output[[paste0("dlHeatmapDFButtonGO")]] <- downloadHandler(
                filename = function() {
                  paste0(design, "_GO_heatmap.xlsx")
                },
                content = function(file) {
                  openxlsx::write.xlsx((as.data.frame(hmGOMatrix) %>% rownames_to_column("gene_name")), file = file)
                }
              )
            }
          })
        }
        
      }
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
)

# Output heatmaps
lapply(c("DEFeatures", "Custom", "GO"), function(sig) {
  output[[paste0("heatmap", sig)]] <- renderPlot(
    {
      req(!is.null(figuresDataReactive[[paste0("heatmap", sig)]]))
      draw(figuresDataReactive[[paste0("heatmap", sig)]], merge_legend = TRUE, padding = unit(c(10, 20, 10, 10), "mm"))
    },
    width = function(){as.numeric(input$figWidthHeatmap)},
    height = function(){as.numeric(input$figHeightHeatmap)}
  )
  output[[paste0("dlHeatmapButton", sig)]] <- downloadHandler(
    filename = function() {
      paste(input$heatmapFilename, sig, tolower(input$heatmapDownloadFormat), sep = ".")
    },
    content = function(file) {
      if (input$heatmapDownloadFormat == "PDF") {
        pdf(file = file, width = as.numeric(input$figWidthHeatmap/60), height = as.numeric(input$figHeightHeatmap/60))
      } else if (input$heatmapDownloadFormat == "SVG") {
        svg(file = file, width = as.numeric(input$figWidthHeatmap/60), height = as.numeric(input$figHeightHeatmap/60))
      } else if (input$heatmapDownloadFormat == "PNG") {
        png(filename = file, width = as.numeric(input$figWidthHeatmap/60), height = as.numeric(input$figHeightHeatmap/60), units = "in", res = as.numeric(input$heatmapDPI))
      }
      draw(figuresDataReactive[[paste0("heatmap", sig)]], merge_legend = TRUE, padding = unit(c(10, 20, 10, 10), "mm"))
      dev.off()
    }
  )
})