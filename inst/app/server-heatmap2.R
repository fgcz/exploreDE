observeEvent(inputDataReactive(), {
  req(inputDataReactive()$dataType)


  if (inputDataReactive()$dataType == "proteomics") {
    output$clusterWarningProteomics <- renderUI({
      helpText(
        "Clustering by row/column may break the heatmap for proteomics data, due
        to the high number of NA values in the count table."
      )
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
    "CountQC" = "Normalised + Log2",
    "proteomics" = if ("transformedData" %in% count_choices) "transformedData" else count_choices[1]
  )
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

  # Automatically set batch correction if grouping2Name exists
  batchDefault <- "None"
  if (!is.null(inputDataReactive()$param$grouping2Name)) {
    if (inputDataReactive()$param$grouping2Name %in% inputDataReactive()$factors) {
      batchDefault <- inputDataReactive()$param$grouping2Name
    }
  }
  updateSelectInput(session = session, inputId = "heatmapBatch", choices = c("None", inputDataReactive()$factors), selected = batchDefault)

  # Update p-value type based on data type
  if (inputDataReactive()$dataType == "proteomics") {
    updateSelectInput(session = session, inputId = "pTypeHeatmap", selected = "Raw")
  } else {
    updateSelectInput(session = session, inputId = "pTypeHeatmap", selected = "FDR")
  }

  # Proteomics-specific: peptide filter UI
  if (inputDataReactive()$dataType == "proteomics") {
    req(inputDataReactive()$seqAnnoList)
    if ("nrPeptides" %in% colnames(inputDataReactive()$seqAnnoList[[1]])) {
      output$nrPeptidesHeatmapUI <- renderUI({
        sliderInput(inputId = "nrPeptideHeatmap", label = "Minimum number of peptides", min = 0, max = 20, value = 2, step = 1, width = "85%")
      })
      outputOptions(output, "nrPeptidesHeatmapUI", suspendWhenHidden = FALSE)
    } else {
      output$nrPeptidesHeatmapUI <- renderUI({
        NULL
      })
    }
  }

  # RNASeq-specific: GO term selection
  if (inputDataReactive()$dataType %in% c("RNASeq", "CountQC")) {
    req(inputDataReactive()$allPathways)
    updateSelectizeInput(
      session,
      inputId = "heatmapGoInput",
      choices = inputDataReactive()$allPathways,
      server = TRUE # ← this is what activates server-side mode
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

output$heatmapColumnSelectUI <- renderUI({
  req(inputDataReactive())
  data <- inputDataReactive()
  req(data$dataType)

  # Calculate default selection
  defaultColumn <- switch(data$dataType,
    "RNASeq" = "groupingName",
    "proteomics" = "factors",
    "CountQC" = "factors"
  )

  defaultSelected <- NULL
  if (!is.null(data[[defaultColumn]]) && length(data[[defaultColumn]]) > 0) {
    defaultSelected <- data[[defaultColumn]][1]
  }

  # Get current selection and available choices
  choices <- data$factors
  currentSelection <- isolate(input$heatmapColumnSelect)

  # Validate current selection against new choices
  validSelection <- intersect(currentSelection, choices)

  # Determine selection to use (Sticky vs Default)
  # If we have a valid current selection, use it. Otherwise default.
  selected <- if (length(validSelection) > 0) validSelection else defaultSelected

  # Apply Mean/Median restriction
  maxItems <- NULL
  if (isTRUE(input$heatmapShowMeans)) {
    maxItems <- 1
    # Note: check length of 'selected' specifically
    if (length(selected) > 1) {
      selected <- selected[1]
      showNotification("Plotting mean/median supports only one factor. Excess factors removed.", type = "warning")
    }
  }

  selectizeInput(
    inputId = "heatmapColumnSelect",
    label = "Select factors to plot by",
    choices = choices,
    selected = selected,
    multiple = TRUE,
    options = list(
      placeholder = "Select features",
      plugins = list("remove_button", "drag_drop", "restore_on_backspace", "clear_button"),
      maxItems = maxItems
    )
  )
})
outputOptions(output, "heatmapColumnSelectUI", suspendWhenHidden = FALSE)

observeEvent(
  {
    genesReactive
  },
  {
    output$geneBucket3 <- renderUI({
      bucket_list(
        header = "Drag and drop features in order to be plotted (Ctrl/Cmd to select multiple)",
        group_name = "bucket_list_group",
        orientation = "horizontal",
        add_rank_list(
          text = "Include these features in this order",
          labels = genesReactive$genes,
          input_id = "keepBucketHeatmap",
          options = sortable::sortable_options(multiDrag = TRUE)
        ),
        add_rank_list(
          text = "Exclude these features",
          labels = NULL,
          input_id = "excludeBucketHeatmap",
          options = sortable::sortable_options(multiDrag = TRUE)
        )
      )
    })
    outputOptions(output, "geneBucket3", suspendWhenHidden = FALSE)
  }
)

observeEvent(input$resetGeneBucketHeatmap,
  {
    genesReactive$genes <- NULL
    updateSelectizeInput(session, inputId = "heatmapGenes", choices = inputDataReactive()$genes, selected = NULL, server = T)
    updateTextAreaInput(session, inputId = "heatmapGenesText", value = NULL)
    output$heatmapCustom <- NULL
  },
  ignoreNULL = T,
  ignoreInit = T
)

observeEvent(
  {
    input$heatmapColumnSelect
  },
  {
    # Clear all buckets first (up to 10)
    for (i in 1:10) {
      output[[paste0("heatmapBucket", i)]] <- renderUI({
        NULL
      })
    }

    # Only render buckets if there are factors selected
    if (!is.null(input$heatmapColumnSelect) && length(input$heatmapColumnSelect) > 0 &&
      !all(input$heatmapColumnSelect == "")) {
      lapply(seq_along(input$heatmapColumnSelect), function(i) {
        output[[paste0("heatmapBucket", i)]] <- renderUI({
          bucket_list(
            header = paste(input$heatmapColumnSelect[[i]], "(Drag/Drop levels - Ctrl/Cmd to select multiple)"),
            group_name = paste0("hmBucket", i),
            orientation = "horizontal",
            add_rank_list(
              text = "Include these groups",
              labels = as.list(unique(levels(as.factor(inputDataReactive()$dataset[[input$heatmapColumnSelect[[i]]]])))),
              input_id = paste0("hmKeepBucket", i),
              options = sortable::sortable_options(multiDrag = TRUE)
            ),
            add_rank_list(
              text = "Exclude these groups",
              labels = NULL,
              input_id = paste0("hmExcludeBucket", i),
              options = sortable::sortable_options(multiDrag = TRUE)
            )
          )
        })
        outputOptions(output, paste0("heatmapBucket", i), suspendWhenHidden = FALSE)
      })
    }
  }
)

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

# Update settings for CountQC data type
observeEvent(inputDataReactive(), {
  if (inputDataReactive()$dataType == "CountQC") {
    # For CountQC, hide gene names by default since there will be many features
    updateCheckboxInput(session = session, inputId = "geneNamesHeatmap", value = FALSE)
    # Set cell border size to 0 for CountQC (will have many features)
    updateSliderInput(session = session, inputId = "cellBorderSizeHeatmap2", value = 0)
  }
})

# Dynamic tab title based on data type
output$heatmapTabTitle <- renderUI({
  req(inputDataReactive())
  if (inputDataReactive()$dataType == "CountQC") {
    "SD Features"
  } else {
    "DE Features"
  }
})

# Helper function to select genes by standard deviation (for CountQC data)
select_genes_by_sd <- function(counts, top_n = 2000) {
  # Calculate SD for each gene
  gene_sds <- apply(counts, 1, sd, na.rm = TRUE)
  # Sort by SD (descending) and take top N
  top_genes <- names(head(sort(gene_sds, decreasing = TRUE), n = top_n))
  return(top_genes)
}

ht_opt$message <- FALSE

# figure making ----
# figure making ----
heatmapInputs <- reactive({
  list(
    input$textSizeHeatmap,
    input$heatmapBatch,
    input$lfcHeatmap,
    input$pTypeHeatmap,
    input$pThresholdHeatmap,
    input$heatmapLimitD,
    input$heatmapLimitCHigh,
    input$heatmapLimitCMid,
    input$heatmapLimitCLow,
    input$heatmapLimitCPalette,
    input$heatmapLimitCPaletteRev,
    input$heatmapGeneNumber,
    input$heatmapGeneNumberCountQC,
    input$heatmapCounts,
    input$clusterRowsHeatmap,
    input$clusterColsHeatmap,
    input$colnamesHeatmap,
    input$geneNamesHeatmap,
    input$heatmapColourRed,
    input$heatmapColourWhite,
    input$heatmapColourBlue,
    input$showClusterColDend,
    input$showClusterRowDend,
    input$heatmapDPI,
    input$heatmapDownloadFormat,
    input$heatmapFilename,
    input$nrPeptideHeatmap,
    input$heatmapFeatureDirection,
    input$heatmapLog2,
    input$heatmapZScore,
    input$heatmapShowMeans,
    input$heatmapAggMethod,
    input$heatmapAtLow,
    input$heatmapAtMid,
    input$heatmapAtHigh,
    input$contrastSelected,
    input$keepBucketHeatmap,
    input$heatmapCustomTitle,
    input$heatmapGoInput,
    input$heatmapGoExpr,
    input$heatmapColumnSelect,
    input$cellBorderSizeHeatmap1,
    input$cellBorderColourHeatmap1,
    input$cellBorderSizeHeatmap2,
    input$cellBorderColourHeatmap2,
    # Add all color picker inputs as dependencies so plot updates when colors change
    lapply(seq_along(inputDataReactive()$factorLevels), function(i) {
      lapply(seq_along(inputDataReactive()$factorLevels[[i]]), function(j) {
        input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]], inputDataReactive()$factorLevels[[i]][[j]])]]
      })
    }),
    lapply(1:10, function(i) {
      input[[paste0("hmKeepBucket", i)]]
    })
  )
})

# Debounce the inputs
heatmapInputsDebounced <- debounce(heatmapInputs, 1000)

observeEvent(
  {
    heatmapInputsDebounced()
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE,
  {
    # Check if user has removed all factors - if so, clear the plot
    if (is.null(input$heatmapColumnSelect) || length(input$heatmapColumnSelect) == 0 ||
      all(input$heatmapColumnSelect == "")) {
      figuresDataReactive[["heatmapDEFeatures"]] <- NULL
      figuresDataReactive[["heatmapCustom"]] <- NULL
      figuresDataReactive[["heatmapGO"]] <- NULL
      return(NULL)
    }

    # tryCatch({
    if (!is.null(input$hmKeepBucket1)) {
      set.seed(2)
      figuresDataReactive[["heatmapDEFeatures"]] <- NULL
      # Get some variables ----
      if (inputDataReactive()$dataType == "proteomics") {
        design <- input$contrastSelected
        seqAnnoHeatmap <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]
        seqAnnoHeatmap[["usedInTest"]] <- TRUE
        if (!is.null(input$nrPeptideHeatmap)) {
          seqAnnoHeatmap <- seqAnnoHeatmap[which(seqAnnoHeatmap$nrPeptides >= input$nrPeptideHeatmap), ]
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

      if (nrow(datasetHeatmap) == 0) {
        return(NULL)
      }

      # Prepare the selected count matrix ----
      countsHeatmap <- inputDataReactive()$countList[[input$heatmapCounts]]
      countsHeatmap <- countsHeatmap[, datasetHeatmap$Name]

      # Store original dataset BEFORE any aggregation for batch correction
      datasetHeatmapOriginal <- datasetHeatmap

      # Apply the following in this specific order: Batch -> Log2 -> Centre/Scale -> Aggregation
      ## 1. If the user wants to remove batch effect, then do that FIRST (use ORIGINAL dataset)
      if (input$heatmapBatch != "None") {
        req(input$heatmapBatch %in% colnames(datasetHeatmapOriginal))
        # Check that batch factor has more than one level after filtering
        batch_levels <- unique(datasetHeatmapOriginal[[input$heatmapBatch]])
        if (length(batch_levels) > 1) {
          countsHeatmap <- limma::removeBatchEffect(x = countsHeatmap, batch = datasetHeatmapOriginal[[input$heatmapBatch]])
        } else {
          showNotification(
            paste0("Batch correction skipped: '", input$heatmapBatch, "' has only ", length(batch_levels), " level in the selected samples"),
            type = "warning",
            duration = 5
          )
        }
      }
      ## 2. If the user wants to log2, then do that
      if (input$heatmapLog2 & input$heatmapCounts %in% c("TPM", "FPKM", "Normalised")) {
        countsHeatmap <- log2(countsHeatmap + 1)
      }
      ## Fix the proteomics rownames
      if (inputDataReactive()$dataType == "proteomics") {
        rownames(countsHeatmap) <- gsub("\\~.*", "", rownames(countsHeatmap))
      }
      ## 3. If the user wants to plot gene-wise Z-scores or centre, then do that
      if (input$heatmapZScore == "Z-Score") {
        countsHeatmap <- t(scale(t(countsHeatmap)))
      }
      if (input$heatmapZScore == "Centred") {
        countsHeatmap <- sweep(countsHeatmap, 1, rowMeans(countsHeatmap, na.rm = T))
      }

      ## 4. Calculate per-condition aggregation if requested (LAST step)
      if (!is.null(input$heatmapShowMeans) && input$heatmapShowMeans) {
        # Get the grouping factor (first selected factor)
        grouping_factor <- input$heatmapColumnSelect[1]

        # Determine aggregation method
        agg_method <- if (!is.null(input$heatmapAggMethod)) input$heatmapAggMethod else "Mean"

        # Get unique factor levels that actually exist in the data
        condition_levels <- lapply(unique(datasetHeatmap[[grouping_factor]]), function(x) {
          as.character(datasetHeatmap$Name[datasetHeatmap[[grouping_factor]] == x])
        }) %>% setNames(unique(datasetHeatmap[[grouping_factor]]))

        # Build a new aggregated matrix
        # Get the mean/median per factor level
        countsHeatmapAgg <- lapply(names(condition_levels), function(x) {
          sample_names <- condition_levels[[x]]
          sample_data <- countsHeatmap[, sample_names, drop = FALSE]
          if (agg_method == "Mean") {
            rowMeans(sample_data, na.rm = TRUE)
          } else {
            matrixStats::rowMedians(as.matrix(sample_data), na.rm = TRUE)
          }
        }) %>%
          setNames(names(condition_levels)) %>%
          data.frame()

        # Convert back to matrix and ensure proper formatting
        countsHeatmap <- as.matrix(countsHeatmapAgg)

        # Update datasetHeatmap to reflect the new structure (one row per condition)
        # Note: datasetHeatmapOriginal already stored above for batch correction
        datasetHeatmap <- data.frame(
          Name = names(condition_levels),
          row.names = NULL,
          stringsAsFactors = FALSE
        )
        # Add the grouping factor column
        datasetHeatmap[[grouping_factor]] <- factor(names(condition_levels), levels = names(condition_levels))

        # Check if additional factors can be annotated when showing means/medians
        # They can only be shown if they don't vary within each aggregated group
        for (f in input$heatmapColumnSelect[-1]) {
          if (f %in% colnames(datasetHeatmapOriginal)) {
            # Check if this factor varies within each condition group
            varies_within_groups <- any(sapply(condition_levels, function(cond) {
              vals <- unique(datasetHeatmapOriginal[[f]][datasetHeatmapOriginal[[grouping_factor]] == cond])
              length(vals) > 1
            }))

            if (varies_within_groups) {
              # Show warning and stop - cannot annotate with this factor
              showNotification(
                paste0("Cannot add '", f, "' to annotation: it has multiple values within '", grouping_factor, "' groups when showing aggregated data (mean/median). Please uncheck 'Show Mean/Median' or remove '", f, "' from the annotation."),
                type = "error",
                duration = 10
              )
              # Skip adding this factor to prevent crash
              next
            } else {
              # Factor is constant within groups - safe to annotate
              datasetHeatmap[[f]] <- sapply(condition_levels, function(cond) {
                first_val <- datasetHeatmapOriginal[[f]][datasetHeatmapOriginal[[grouping_factor]] == cond][1]
                as.character(first_val)
              })
              datasetHeatmap[[f]] <- factor(datasetHeatmap[[f]])
            }
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
          "HeatmapAnnotation(",
          p1, ", show_legend = T, annotation_legend_param = list(title_gp = gpar(fontsize = ", input$textSizeHeatmap, ", fontface = 'bold'), labels_gp = gpar(fontsize = ", input$textSizeHeatmap, ")), gp = gpar(col = '", input$cellBorderColourHeatmap1, "', lwd = ", input$cellBorderSizeHeatmap1 / 10, "))"
        )
        ha <- eval(parse(text = ha))
        # Only apply colors for factors that actually exist in datasetHeatmap
        for (f in c(input$heatmapColumnSelect)) {
          if (f %in% colnames(datasetHeatmap)) {
            for (l in levels(as.factor(datasetHeatmap[, f]))) {
              ha@anno_list[[f]]@color_mapping@full_col[[l]] <- colourListHeatmap[[f]][[l]]
              ha@anno_list[[f]]@fun@var_env$color_mapping@full_col[[l]] <- colourListHeatmap[[f]][[l]]
            }
          }
        }
      }

      # DE Heatmaps ----
      # Select genes based on data type
      if (inputDataReactive()$dataType == "CountQC") {
        # For CountQC: select by standard deviation
        top_n_countqc <- if (!is.null(input$heatmapGeneNumberCountQC)) input$heatmapGeneNumberCountQC else 2000
        top_genes <- select_genes_by_sd(
          counts = countsHeatmap,
          top_n = top_n_countqc
        )
        ntp <- length(top_genes)
      } else {
        # For RNASeq/proteomics: existing DE-based selection
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
        # Get the top gene names
        top_genes <- geneDirections[[input$heatmapFeatureDirection]][1:ntp]
      }

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

      # Determine aggregation label for title (used by all heatmap types)
      agg_label <- ""
      if (!is.null(input$heatmapShowMeans) && input$heatmapShowMeans) {
        agg_method <- if (!is.null(input$heatmapAggMethod)) input$heatmapAggMethod else "Mean"
        agg_label <- paste0("(Per-condition ", tolower(agg_method), "s)")
      }

      # Only generate DE Features heatmap if there are genes to display
      if (ntp > 0) {
        # Filter heatmap data to only include genes that exist in both lists
        # This combines the matching and filtering in one step
        countsHeatmapSig <- countsHeatmap[rownames(countsHeatmap) %in% top_genes, , drop = FALSE]

        # Reorder to match the original gene direction order
        gene_order <- top_genes[top_genes %in% rownames(countsHeatmapSig)]
        countsHeatmapSig <- countsHeatmapSig[gene_order, , drop = FALSE]

        output[["heatmapDesignBoth Directions"]] <- renderText({
          design
        })

        # Limit to requested number of features based on data type
        if (inputDataReactive()$dataType == "CountQC") {
          max_features <- min(as.numeric(input$heatmapGeneNumberCountQC), nrow(countsHeatmapSig))
          countsHeatmapSig <- countsHeatmapSig[1:max_features, ]
        } else {
          if (as.numeric(input$heatmapGeneNumber) <= nrow(countsHeatmapSig)) {
            countsHeatmapSig <- countsHeatmapSig[1:as.numeric(input$heatmapGeneNumber), ]
          }
        }
        rownames(countsHeatmapSig) <- substr(rownames(countsHeatmapSig), 1, 15)

        # Determine heatmap title based on data type
        heatmap_title <- if (inputDataReactive()$dataType == "CountQC") {
          paste("Top", nrow(countsHeatmapSig), "Features by SD", agg_label, ":\n", design)
        } else {
          paste("Top", nrow(countsHeatmapSig), "DE Features", agg_label, ":\n", design)
        }

        # For heatmaps with many features (>200), use no cell border to avoid clutter
        cell_border_size <- if (nrow(countsHeatmapSig) > 200) 0 else input$cellBorderSizeHeatmap2

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
          rect_gp = gpar(col = input$cellBorderColourHeatmap2, lwd = cell_border_size / 10),
          use_raster = FALSE,
          row_dend_width = unit(1, "cm"),
          column_dend_height = unit(1, "cm"),
          col = heatmapColours,
          row_names_gp = gpar(fontsize = input$textSizeHeatmap),
          column_names_gp = gpar(fontsize = input$textSizeHeatmap),
          column_title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
          column_title = heatmap_title,
          heatmap_legend_param = list(labels_gp = gpar(fontsize = input$textSizeHeatmap), title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"))
        )
        figuresDataReactive[["heatmapDEFeatures"]] <- heatmap

        # Download the count file
        countsHeatmapSigDL <- as.data.frame(countsHeatmapSig) %>% rownames_to_column("gene_name")
        output[["dlHeatmapDFButtonDEFeatures"]] <- downloadHandler(
          filename = function() {
            paste0(design, "_DE_", "_heatmap_counts.xlsx")
          },
          content = function(file) {
            openxlsx::write.xlsx(countsHeatmapSigDL, file = file)
          }
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
          filename = function() {
            "no_data.txt"
          },
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
        countsHeatmapCustom <- t(as.matrix(countsHeatmap[input$keepBucketHeatmap, ]))
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
        rect_gp = gpar(col = input$cellBorderColourHeatmap2, lwd = input$cellBorderSizeHeatmap2 / 10),
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
      if (inputDataReactive()$dataType %in% c("RNASeq", "CountQC")) {
        goToPlot <- NULL
        goGenes <- NULL
        countsHeatmapGO <- NULL

        goToPlot <- gsub(" .*", "", input$heatmapGoInput)
        goGenes <- seqAnnoHeatmap$gene_name[with(seqAnnoHeatmap, grepl(goToPlot, paste(`GO BP`, `GO MF`, `GO CC`)))]

        goGenes <- goGenes[which(goGenes %in% rownames(countsHeatmap))]
        if (length(goGenes) == 1) {
          countsHeatmapGO <- t(as.matrix(countsHeatmap[goGenes, ]))
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
          rect_gp = gpar(col = input$cellBorderColourHeatmap2, lwd = input$cellBorderSizeHeatmap2 / 10),
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
    width = function() {
      as.numeric(input$figWidthHeatmap)
    },
    height = function() {
      as.numeric(input$figHeightHeatmap)
    }
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
        pdf(file = file, width = as.numeric(input$figWidthHeatmap / 60), height = as.numeric(input$figHeightHeatmap / 60))
      } else if (input$heatmapDownloadFormat == "SVG") {
        svg(file = file, width = as.numeric(input$figWidthHeatmap / 60), height = as.numeric(input$figHeightHeatmap / 60))
      } else if (input$heatmapDownloadFormat == "PNG") {
        png(filename = file, width = as.numeric(input$figWidthHeatmap / 60), height = as.numeric(input$figHeightHeatmap / 60), units = "in", res = as.numeric(input$heatmapDPI))
      }
      draw(heatmap_obj, merge_legend = TRUE, padding = unit(c(10, 20, 10, 10), "mm"))
      dev.off()
    }
  )
})
