zscore <- function(x) {
  # Remove NAs for mean and sd calculation
  if (all(is.na(x))) {
    # If all values are NA, return all NAs
    return(rep(NA, length(x)))
  }
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    # If sd is NA (all NAs) or zero (no variance), return zeros (or NAs)
    return(rep(0, length(x)))
  }
  (x - m) / s
}

req(!is.null(inputDataReactive()$dataType))

if (inputDataReactive()$dataType == "proteomics") {
  output$heatmapProteomicsColumnSelectUI <- renderUI({
    selectInput(inputId = "heatmapProteomicsColumnSelect", label = "Select main factor to plot by", choices = inputDataReactive()$factors, selected = inputDataReactive()$factors[1], multiple = F)
  })
  outputOptions(output, "heatmapProteomicsColumnSelectUI", suspendWhenHidden = FALSE)
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
  outputOptions(output, "clusterWarningProteomics", suspendWhenHidden = FALSE)
  updateCheckboxInput(session = session, inputId = "clusterColsHeatmap", value = FALSE)
  updateCheckboxInput(session = session, inputId = "clusterRowsHeatmap", value = FALSE)
  updateCheckboxInput(session = session, inputId = "heatmapLog2", value = FALSE)
}

# Update inputs ----
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
  outputOptions(output, "geneBucket3", suspendWhenHidden = FALSE)
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
  outputOptions(output, "heatmapBucket", suspendWhenHidden = FALSE)
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
    updateSelectizeInput(
      session,
      inputId = "heatmapGoInput",
      choices = inputDataReactive()$allPathways,
      server = TRUE  # ← this is what activates server-side mode
    )
    output$heatmapGOUI2 <- renderUI({
      numericInput(
        inputId = "heatmapGoExpr",
        label = "Minimum median feature expression for GO plot",
        value = 0.05,
        min = 0,
        max = 1e100
      )
    })
    outputOptions(output, "heatmapGOUI2", suspendWhenHidden = FALSE)
  } else {
    output$heatmapGOUI1 <- renderUI({
      h5("Proteomics data doesn't support GO terms.")
    })
  }
})

updateSelectInput(session = session, inputId = "heatmapBatch", choices = c("None", inputDataReactive()$factors), selected = "None")

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
    outputOptions(output, "nrPeptidesHeatmapUI", suspendWhenHidden = FALSE)
  } else {
    output$nrPeptidesHeatmapUI <- renderUI({ NULL })
  }
}

observeEvent(input$heatmapZScore, {
  if (input$heatmapZScore %in% c("Z-Score", "Centred")) {
    updateNumericInput(session = session, inputId = "heatmapAtLow", value = -2, min = -10, max = 10, step = 0.1)
    updateNumericInput(session = session, inputId = "heatmapAtMid", value = 0, min = -10, max = 10, step = 0.1)
    updateNumericInput(session = session, inputId = "heatmapAtHigh", value = 2, min = -10, max = 10, step = 0.1)
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



# figure making ----
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
    input$keepBucketHeatmap
    input$heatmapCustomTitle
    input$heatmapGoInput
    input$heatmapGoExpr
  }, ignoreNULL = FALSE, ignoreInit = TRUE, {
    tryCatch({
      if (!is.null(input$hmKeepBucket)) {
        
        figuresDataReactive[["heatmapDEFeatures"]] <- NULL
        # Get some variables ----
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
        
        
        # Prepare the metadata for the top annotation ----
        datasetHeatmap <- inputDataReactive()$dataset
        datasetHeatmap <- datasetHeatmap[which(datasetHeatmap[[groupingNameHeatmap()$gn]] %in% input$hmKeepBucket), ]
        # Set the factor levels per the bucket list order
        datasetHeatmap[[groupingNameHeatmap()$gn]] <- factor(datasetHeatmap[[groupingNameHeatmap()$gn]], input$hmKeepBucket)
        datasetHeatmap <- datasetHeatmap[order(datasetHeatmap[[groupingNameHeatmap()$gn]]), ]
        
        
        # Prepare the selected count matrix ----
        countsHeatmap <- inputDataReactive()$countList[[input$heatmapCounts]]
        countsHeatmap <- countsHeatmap[, rownames(datasetHeatmap)]
        # Apply the following in this specific order as required: Log2, batch correct, Z-scale
        ## If the user wants to log2, then do that
        if (input$heatmapLog2 & input$heatmapCounts %in% c("TPM", "FPKM", "Normalised")) {
          countsHeatmap <- log2(countsHeatmap + 1)
        }
        ## If the user wants to remove batch effect, then do that
        if (input$heatmapBatch != "None") {
          countsHeatmap <- limma::removeBatchEffect(countsHeatmap, datasetHeatmap[[input$heatmapBatch]])
        }
        ## If the user wants to plot gene-wise Z-scores, then do that 
        if (input$heatmapZScore == "Z-Score") {
          countsHeatmapz <- t(apply(countsHeatmap, 1, zscore))
          colnames(countsHeatmapz) <- colnames(countsHeatmap)
          countsHeatmap <- countsHeatmapz
          rm(countsHeatmapz)
        }
        ## If the user wants to centre the counts, then do that
        if (input$heatmapZScore == "Centred") {
          countsHeatmap <- sweep(countsHeatmap, 1, rowMeans(countsHeatmap, na.rm = T))
        }
        ## Fix the proteomics rownames
        if (inputDataReactive()$dataType == "proteomics") {
          rownames(countsHeatmap) <- gsub("\\~.*", "", rownames(countsHeatmap))
        }
        
        
        # Prepare colours and generate top annotation ----
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
            paste0(col2hex(input[[paste0("GroupColour", f, k)]]), "FF")
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
        
        
        # DE Heatmaps ----
        # Get the features only in the specified DE direction
        geneDirections <- list(
          "Both" = seqAnnoHeatmap$gene_name[
            which(
              seqAnnoHeatmap[[pTypeHeatmap]] <= as.numeric(input$pThresholdHeatmap) & 
                abs(seqAnnoHeatmap$log2Ratio) >= as.numeric(input$lfcHeatmap)
              )
            ],
          "Up-regulated" = seqAnnoHeatmap$gene_name[
            which(
              seqAnnoHeatmap[[pTypeHeatmap]] <= as.numeric(input$pThresholdHeatmap) & 
                seqAnnoHeatmap$log2Ratio >= as.numeric(input$lfcHeatmap)
            )
          ],
          "Down-regulated" = seqAnnoHeatmap$gene_name[
            which(
              seqAnnoHeatmap[[pTypeHeatmap]] <= as.numeric(input$pThresholdHeatmap) & 
                seqAnnoHeatmap$log2Ratio <= -as.numeric(input$lfcHeatmap)
            )
          ]
        )
        
        ntp <- min(input$heatmapGeneNumber, length(geneDirections[[input$heatmapFeatureDirection]]))
        gtp <- rownames(countsHeatmap)[match(geneDirections[[input$heatmapFeatureDirection]][1:ntp], rownames(countsHeatmap))]
        rtp <- grep(paste0("^", gtp, "$", collapse = "|"), rownames(countsHeatmap))
        
        countsHeatmapSig <- countsHeatmap[rtp, ]
        countsHeatmapSig <- countsHeatmapSig[match(gtp, rownames(countsHeatmapSig)),]
        countsHeatmapSig <- countsHeatmapSig[!is.na(rownames(countsHeatmapSig)), ]
        output[["heatmapDesignBoth Directions"]] <- renderText({design})
        if (as.numeric(input$heatmapGeneNumber) <= nrow(countsHeatmapSig)) {
          countsHeatmapSig <- countsHeatmapSig[1:as.numeric(input$heatmapGeneNumber), ]
        } 
        rownames(countsHeatmapSig) <- substr(rownames(countsHeatmapSig), 1, 15)
        
        hmName <- paste(input$heatmapCounts %>% gsub("\\+", "\n\\+", .))
        if (input$heatmapLog2 & input$heatmapCounts %in% c("TPM", "FPKM", "Normalised")) {
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
        )
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
        
        
        
        # Custom heatmap ----
        req(input$keepBucketHeatmap != "")
        
        if (length(input$keepBucketHeatmap) == 1) {
          countsHeatmapCustom <- t(as.matrix(countsHeatmap[input$keepBucketHeatmap,]))
          rownames(countsHeatmapCustom) <- input$keepBucketHeatmap
        } else {
          countsHeatmapCustom <- countsHeatmap[input$keepBucketHeatmap, ]
        }
        req(nrow(countsHeatmapCustom) >= 1)
        ch <- ComplexHeatmap::Heatmap(
          matrix = countsHeatmapCustom,
          top_annotation = ha,
          name = input$heatmapCounts,
          cluster_columns = input$clusterColsHeatmap,
          cluster_rows = input$clusterRowsHeatmap,
          show_column_names = input$colnamesHeatmap,
          show_row_names = input$geneNamesHeatmap,
          col = heatmapColours,
          column_title = paste(input$heatmapCustomTitle),
          rect_gp = gpar(col = "white", lwd = 0.5)
        )
        figuresDataReactive$heatmapCustom <- ch
        
        output[[paste0("dlHeatmapDFButtonCustom")]] <- downloadHandler(
          filename = function() {
            paste0(design, "_Custom_heatmap.xlsx")
          },
          content = function(file) {
            openxlsx::write.xlsx((as.data.frame(countsHeatmapCustom) %>% rownames_to_column("gene_name")), file = file)
          }
        )
        
        
        # GO heatmap ----
        if (inputDataReactive()$dataType == "RNASeq") {
          goToPlot <- NULL
          goGenes <- NULL
          countsHeatmapGO <- NULL
          
          req(input$heatmapGoInput != "")
          
          goToPlot <- gsub(" .*", "", input$heatmapGoInput)
          goGenes <- seqAnnoHeatmap$gene_name[with(seqAnnoHeatmap, grepl(goToPlot, paste(`GO BP`, `GO MF`, `GO CC`)))]
          
          goGenes <- goGenes[which(goGenes %in% rownames(countsHeatmap))]
          if (length(goGenes) == 1) {
            countsHeatmapGO <- t(as.matrix(countsHeatmap[goGenes,]))
            rownames(countsHeatmapGO) <- goGenes
          } else {
            countsHeatmapGO <- countsHeatmap[goGenes, ]
            countsHeatmapGO <- countsHeatmapGO[abs(rowMedians(countsHeatmapGO)) >= input$heatmapGoExpr, ]
          }
          req(nrow(countsHeatmapGO) >= 1)
          req(ncol(countsHeatmapGO) == ha@anno_list[[1]]@fun@n)
          
          gh <- ComplexHeatmap::Heatmap(
            matrix = countsHeatmapGO,
            top_annotation = ha,
            name = input$heatmapCounts,
            cluster_columns = input$clusterColsHeatmap,
            cluster_rows = input$clusterRowsHeatmap,
            show_column_names = input$colnamesHeatmap,
            show_row_names = input$geneNamesHeatmap,
            col = heatmapColours,
            column_title = sub("\\s+", "\n", input$heatmapGoInput),
            rect_gp = gpar(col = "white", lwd = 0.5)
          )
          figuresDataReactive$heatmapGO <- gh
          output[[paste0("dlHeatmapDFButtonGO")]] <- downloadHandler(
            filename = function() {
              paste0(design, "_GO_heatmap.xlsx")
            },
            content = function(file) {
              openxlsx::write.xlsx((as.data.frame(countsHeatmapGO) %>% rownames_to_column("gene_name")), file = file)
            }
          )
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