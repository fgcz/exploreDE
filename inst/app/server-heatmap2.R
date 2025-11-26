req(inputDataReactive()$dataType)
defaultColumn <- switch(inputDataReactive()$dataType, "RNASeq" = "groupingName", "proteomics" = "factors")
output$heatmapColumnSelectUI <- renderUI({
  selectizeInput(
    inputId = "heatmapColumnSelect", 
    label = "Select factors to plot by", 
    choices = inputDataReactive()$factors, 
    selected = inputDataReactive()[[defaultColumn]][1], 
    multiple = T,
    options = list(placeholder = 'Select features', plugins = list('remove_button', 'drag_drop', 'restore_on_backspace', 'clear_button'))
  )
})
outputOptions(output, "heatmapColumnSelectUI", suspendWhenHidden = FALSE)

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
count_choices <- names(inputDataReactive()$countList)
default_count <- switch(inputDataReactive()$dataType,
                       "RNASeq" = "Normalised + Log2",
                       "proteomics" = if ("transformedData" %in% count_choices) "transformedData" else count_choices[1])
updateSelectInput(
  session = session,
  inputId = "heatmapCounts",
  choices = count_choices,
  selected = default_count
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

observeEvent({
  input$heatmapColumnSelect
}, {
  lapply(seq_along(input$heatmapColumnSelect), function(i) {
    output[[paste0("heatmapBucket", i)]] <- renderUI({
      bucket_list(
        header = input$heatmapColumnSelect[[i]],
        group_name = paste0("hmBucket", i),
        orientation = "horizontal",
        add_rank_list(
          text = "Include these groups",
          labels = as.list(unique(levels(as.factor(inputDataReactive()$dataset[[input$heatmapColumnSelect[[i]]]])))),
          input_id = paste0("hmKeepBucket", i)
        ),
        add_rank_list(
          text = "Exclude these groups",
          labels = NULL,
          input_id = paste0("hmExcludeBucket", i)
        )
      )
    })
    outputOptions(output, paste0("heatmapBucket", i), suspendWhenHidden = FALSE)
  })
})


observeEvent(inputDataReactive(), {
  if (inputDataReactive()$dataType == "RNASeq") {
    req(inputDataReactive()$allPathways)
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
  req(inputDataReactive()$seqAnnoList)
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

ht_opt$message = FALSE

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
    input$heatmapShowMeans
    input$heatmapAggMethod
    input$heatmapAtLow
    input$heatmapAtMid
    input$heatmapAtHigh
    input$contrastSelected
    input$keepBucketHeatmap
    input$heatmapCustomTitle
    input$heatmapGoInput
    input$heatmapGoExpr
    input$heatmapColumnSelect
    input$cellBorderSizeHeatmap1
    input$cellBorderColourHeatmap1
    input$cellBorderSizeHeatmap2
    input$cellBorderColourHeatmap2
    lapply(seq_along(inputDataReactive()$factorLevels), function (i) {
      input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]])]]
    })
    lapply(1:10, function(i) {
      input[[paste0("hmKeepBucket", i)]]
    })
  }, ignoreNULL = FALSE, ignoreInit = TRUE, {
    # tryCatch({
    if (!is.null(input$hmKeepBucket1)) {
      
      set.seed(2)
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
      for (i in seq_along(input$heatmapColumnSelect)) {
        datasetHeatmap <- datasetHeatmap[which(datasetHeatmap[[input$heatmapColumnSelect[[i]]]] %in% input[[paste0("hmKeepBucket", i)]]), ]
        datasetHeatmap[[input$heatmapColumnSelect[[i]]]] <- factor(datasetHeatmap[[input$heatmapColumnSelect[[i]]]], levels = input[[paste0("hmKeepBucket", i)]])
      }
      # Get the sample order
      datasetHeatmap <- datasetHeatmap %>% dplyr::arrange(across(all_of(input$heatmapColumnSelect)))
      
      
      # Prepare the selected count matrix ----
      countsHeatmap <- inputDataReactive()$countList[[input$heatmapCounts]]
      countsHeatmap <- countsHeatmap[, datasetHeatmap$Name]
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
        countsHeatmap <- t(scale(t(countsHeatmap)))
      }
      ## If the user wants to centre the counts, then do that
      if (input$heatmapZScore == "Centred") {
        countsHeatmap <- sweep(countsHeatmap, 1, rowMeans(countsHeatmap, na.rm = T))
      }
      ## Fix the proteomics rownames
      if (inputDataReactive()$dataType == "proteomics") {
        rownames(countsHeatmap) <- gsub("\\~.*", "", rownames(countsHeatmap))
      }

      ## Calculate per-condition aggregation if requested
      if (!is.null(input$heatmapShowMeans) && input$heatmapShowMeans) {
        # Get the grouping factor (first selected factor)
        grouping_factor <- input$heatmapColumnSelect[1]

        # Determine aggregation method
        agg_method <- if (!is.null(input$heatmapAggMethod)) input$heatmapAggMethod else "Mean"

        # Calculate aggregation for each condition
        condition_levels <- levels(as.factor(datasetHeatmap[[grouping_factor]]))

        # Store original rownames before aggregation
        original_rownames <- rownames(countsHeatmap)

        countsHeatmapAgg <- lapply(condition_levels, function(cond) {
          sample_names <- datasetHeatmap$Name[datasetHeatmap[[grouping_factor]] == cond]
          sample_data <- countsHeatmap[, sample_names, drop = FALSE]
          if (agg_method == "Mean") {
            rowMeans(sample_data, na.rm = TRUE)
          } else {
            rowMedians(as.matrix(sample_data), na.rm = TRUE)
          }
        })
        countsHeatmap <- do.call(cbind, countsHeatmapAgg)
        colnames(countsHeatmap) <- condition_levels
        rownames(countsHeatmap) <- original_rownames

        # Clean up NA/NaN/Inf values that can break clustering
        countsHeatmap[is.na(countsHeatmap) | is.nan(countsHeatmap) | is.infinite(countsHeatmap)] <- 0

        # Update datasetHeatmap to reflect the new structure (one row per condition)
        datasetHeatmapOriginal <- datasetHeatmap  # Keep original for colors
        datasetHeatmap <- data.frame(
          Name = condition_levels,
          row.names = NULL,
          stringsAsFactors = FALSE
        )
        # Add the grouping factor column
        datasetHeatmap[[grouping_factor]] <- factor(condition_levels, levels = condition_levels)
        # Add other factors with their first occurrence values (for annotation)
        for (f in input$heatmapColumnSelect[-1]) {
          if (f %in% colnames(datasetHeatmapOriginal)) {
            datasetHeatmap[[f]] <- sapply(condition_levels, function(cond) {
              first_val <- datasetHeatmapOriginal[[f]][datasetHeatmapOriginal[[grouping_factor]] == cond][1]
              as.character(first_val)
            })
            datasetHeatmap[[f]] <- factor(datasetHeatmap[[f]])
          }
        }
      }

      # Prepare colours and generate top annotation ----
      at <- c(input$heatmapAtLow, input$heatmapAtMid, input$heatmapAtHigh)
      if (input$heatmapLimitCPalette == "None") {
        heatmapColours <- colorRamp2(
          breaks = rev(at),
          colors = c(input$heatmapColourRed, input$heatmapColourWhite, input$heatmapColourBlue)
        )
      } else {
        palette_colours <- brewer.pal(9, input$heatmapLimitCPalette)
        selected_colours <- palette_colours[c(1, 5, 9)]
        if (input$heatmapLimitCPaletteRev) {
          heatmapColours <- colorRamp2(breaks = at, colors = selected_colours)
        } else {
          heatmapColours <- colorRamp2(at, rev(selected_colours))
        }
      }
      # Top Annotation
      colourListHeatmap <- lapply(input$heatmapColumnSelect, function(f) {
        lapply(levels(datasetHeatmap[[f]]), function(l) {
          paste0(col2hex(input[[paste0("GroupColour", f, l)]]), "FF")
        }) %>% setNames((levels(datasetHeatmap[[f]])))
      }) %>% setNames(input$heatmapColumnSelect)
      if (any(input$heatmapColumnSelect != "")) {
        p1 <- paste(lapply(input$heatmapColumnSelect, function(g) {
          paste0("'", g, "' = datasetHeatmap[['", g, "']]")
        }), collapse = ",")
        ha <- paste0(
          "HeatmapAnnotation("
          , p1, ", show_legend = T, annotation_legend_param = list(title_gp = gpar(fontsize = ", input$textSizeHeatmap, ", fontface = 'bold'), labels_gp = gpar(fontsize = ", input$textSizeHeatmap, ")), gp = gpar(col = '", input$cellBorderColourHeatmap1, "', lwd = ", input$cellBorderSizeHeatmap1/10,"))"
        )
        ha <- eval(parse(text = ha))
        for (f in c(input$heatmapColumnSelect)) {
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
      
      # Get the number of top genes to display
      ntp <- min(input$heatmapGeneNumber, length(geneDirections[[input$heatmapFeatureDirection]]))

      # Prepare heatmap name (used by all heatmaps)
      hmName <- paste(input$heatmapCounts %>% gsub("\\+", "\n\\+", .))
      if (input$heatmapLog2 & input$heatmapCounts %in% c("TPM", "FPKM", "Normalised")) {
        hmName <- paste(hmName, "\n(Log2)")
      }
      if (input$heatmapZScore == "Z-Score") {
        hmName <- paste(hmName, "\n[Z Score]")
      } else if (input$heatmapZScore == "Centred") {
        hmName <- paste(hmName, "\n[Centred]")
      }

      # Only generate DE Features heatmap if there are genes to display
      if (ntp > 0) {
        # Get the top gene names
        top_genes <- geneDirections[[input$heatmapFeatureDirection]][1:ntp]

      # Filter heatmap data to only include genes that exist in both lists
      # This combines the matching and filtering in one step
      countsHeatmapSig <- countsHeatmap[rownames(countsHeatmap) %in% top_genes, , drop = FALSE]

      # Reorder to match the original gene direction order
      gene_order <- top_genes[top_genes %in% rownames(countsHeatmapSig)]
      countsHeatmapSig <- countsHeatmapSig[gene_order, , drop = FALSE]

      output[["heatmapDesignBoth Directions"]] <- renderText({design})
      if (as.numeric(input$heatmapGeneNumber) <= nrow(countsHeatmapSig)) {
        countsHeatmapSig <- countsHeatmapSig[1:as.numeric(input$heatmapGeneNumber), ]
      }
      rownames(countsHeatmapSig) <- substr(rownames(countsHeatmapSig), 1, 15)
      
      # Determine aggregation label for title
      agg_label <- ""
      if (!is.null(input$heatmapShowMeans) && input$heatmapShowMeans) {
        agg_method <- if (!is.null(input$heatmapAggMethod)) input$heatmapAggMethod else "Mean"
        agg_label <- paste0("(Per-condition ", tolower(agg_method), "s)")
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
        rect_gp = gpar(col = input$cellBorderColourHeatmap2, lwd = input$cellBorderSizeHeatmap2/10),
        use_raster = FALSE,
        row_dend_width = unit(1, "cm"),
        column_dend_height = unit(1, "cm"),
        col = heatmapColours,
        row_names_gp = gpar(fontsize = input$textSizeHeatmap),
        column_names_gp = gpar(fontsize = input$textSizeHeatmap),
        column_title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
        column_title = paste("Top", nrow(countsHeatmapSig), "DE Features", agg_label, ":\n", design),
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
      } else {
        # When no DEGs are found (ntp = 0), provide informative placeholders
        figuresDataReactive[["heatmapDEFeatures"]] <- NULL

        output[["heatmapDesignBoth Directions"]] <- renderText({
          paste("No features found meeting current thresholds for:", design)
        })

        output[["hGNDE"]] <- renderText({
          "Number of features per settings chosen: 0"
        })

        output[["dlHeatmapDFButtonDEFeatures"]] <- downloadHandler(
          filename = function() { "no_data.txt" },
          content = function(file) {
            writeLines(
              paste0(
                "No differentially expressed features found for: ", design, "\n\n",
                "Current filter settings:\n",
                "- P-value type: ", input$pTypeHeatmap, "\n",
                "- P-value threshold: ", input$pThresholdHeatmap, "\n",
                "- Log2 fold-change threshold: ", input$lfcHeatmap, "\n",
                "- Feature direction: ", input$heatmapFeatureDirection, "\n\n",
                "Try adjusting these thresholds to find more features."
              ),
              file
            )
          }
        )
      } # End if (ntp > 0) - DE Features heatmap section


      # Custom heatmap ----
      if (length(input$keepBucketHeatmap) == 1) {
        countsHeatmapCustom <- t(as.matrix(countsHeatmap[input$keepBucketHeatmap,]))
        rownames(countsHeatmapCustom) <- input$keepBucketHeatmap
      } else {
        countsHeatmapCustom <- countsHeatmap[input$keepBucketHeatmap, ]
      }
      ch <- ComplexHeatmap::Heatmap(
        matrix = countsHeatmapCustom,
        top_annotation = ha,
        name = hmName,
        cluster_columns = input$clusterColsHeatmap,
        cluster_rows = input$clusterRowsHeatmap,
        show_column_names = input$colnamesHeatmap,
        show_row_names = input$geneNamesHeatmap,
        show_column_dend = input$showClusterColDend,
        show_row_dend = input$showClusterRowDend,
        rect_gp = gpar(col = input$cellBorderColourHeatmap2, lwd = input$cellBorderSizeHeatmap2/10),
        use_raster = FALSE,
        row_dend_width = unit(1, "cm"),
        column_dend_height = unit(1, "cm"),
        col = heatmapColours,
        row_names_gp = gpar(fontsize = input$textSizeHeatmap),
        column_names_gp = gpar(fontsize = input$textSizeHeatmap),
        column_title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
        column_title = paste(input$heatmapCustomTitle, agg_label),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = input$textSizeHeatmap), title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"))
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
        req(nrow(datasetHeatmap) == ncol(countsHeatmapGO))
        gh <- ComplexHeatmap::Heatmap(
          matrix = countsHeatmapGO,
          top_annotation = ha,
          name = hmName,
          cluster_columns = input$clusterColsHeatmap,
          cluster_rows = input$clusterRowsHeatmap,
          show_column_names = input$colnamesHeatmap,
          show_row_names = input$geneNamesHeatmap,
          show_column_dend = input$showClusterColDend,
          show_row_dend = input$showClusterRowDend,
          rect_gp = gpar(col = input$cellBorderColourHeatmap2, lwd = input$cellBorderSizeHeatmap2/10),
          use_raster = FALSE,
          row_dend_width = unit(1, "cm"),
          column_dend_height = unit(1, "cm"),
          col = heatmapColours,
          row_names_gp = gpar(fontsize = input$textSizeHeatmap),
          column_names_gp = gpar(fontsize = input$textSizeHeatmap),
          column_title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
          column_title = paste(sub("\\s+", "\n", input$heatmapGoInput), if (agg_label != "") paste0("\n", agg_label) else ""),
          heatmap_legend_param = list(labels_gp = gpar(fontsize = input$textSizeHeatmap), title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"))
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
    # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
)

# Output heatmaps
lapply(c("DEFeatures", "Custom", "GO"), function(sig) {
  output[[paste0("heatmap", sig)]] <- renderPlot(
    {
      req(figuresDataReactive[[paste0("heatmap", sig)]])
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
      # Safety check: ensure heatmap object exists before drawing
      heatmap_obj <- figuresDataReactive[[paste0("heatmap", sig)]]
      if (is.null(heatmap_obj)) {
        # Create informative placeholder file
        if (tolower(input$heatmapDownloadFormat) == "pdf") {
          pdf(file = file, width = 8, height = 6)
          plot.new()
          text(0.5, 0.5, paste("No data available for", sig, "heatmap\n\nPlease adjust filter settings or select features."), cex = 1.2)
          dev.off()
        } else {
          writeLines(paste("No data available for", sig, "heatmap. Please adjust filter settings or select features."), file)
        }
        return()
      }

      # Heatmap exists, proceed with normal download
      if (input$heatmapDownloadFormat == "PDF") {
        pdf(file = file, width = as.numeric(input$figWidthHeatmap/60), height = as.numeric(input$figHeightHeatmap/60))
      } else if (input$heatmapDownloadFormat == "SVG") {
        svg(file = file, width = as.numeric(input$figWidthHeatmap/60), height = as.numeric(input$figHeightHeatmap/60))
      } else if (input$heatmapDownloadFormat == "PNG") {
        png(filename = file, width = as.numeric(input$figWidthHeatmap/60), height = as.numeric(input$figHeightHeatmap/60), units = "in", res = as.numeric(input$heatmapDPI))
      }
      draw(heatmap_obj, merge_legend = TRUE, padding = unit(c(10, 20, 10, 10), "mm"))
      dev.off()
    }
  )
})