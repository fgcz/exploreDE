# GO heatmap: committed GO term (set on button click)
goHeatmapReactive <- reactiveValues(
  currentGoTerm = NULL,
  currentGoName = NULL,
  currentGenes = NULL
)

# KEGG heatmap reactive cache
keggHeatmapReactive <- reactiveValues(
  keggCode = NULL,
  cache = list(),
  currentGenes = NULL,
  currentPathwayName = NULL
)

observeEvent(inputDataReactive(), {
  req(inputDataReactive()$dataType)

  # Detect species for KEGG heatmap
  keggHeatmapReactive$keggCode <- tryCatch({
    spp <- ezRun::getSpecies(inputDataReactive()$param$refBuild)
    if (grepl("^Rattus", inputDataReactive()$param$refBuild)) spp <- "Rat"
    switch(spp, "Human" = "hsa", "Mouse" = "mmu", "Rat" = "rno", NULL)
  }, error = function(e) NULL)

  # Output flag for conditionalPanel in UI
  output$keggSpeciesSupported <- reactive({ !is.null(keggHeatmapReactive$keggCode) })
  outputOptions(output, "keggSpeciesSupported", suspendWhenHidden = FALSE)

  # Add KEGG heatmap tab dynamically for supported species
  if (!is.null(keggHeatmapReactive$keggCode)) {
    appendTab(
      inputId = "sig",
      tabPanel(
        title = "KEGG",
        column(
          width = 3,
          helpText(
            "Plot the expression of all features in a KEGG pathway.",
            "Unlike the KEGG Pathway tab, this shows individual gene expression",
            "rather than collapsing multiple genes to a single KEGG ortholog."
          ),
          textInput(
            inputId = "heatmapKeggInput",
            label = "Enter a KEGG pathway ID",
            placeholder = "e.g., hsa05200 or 05200",
            value = ""
          ),
          actionButton(
            inputId = "heatmapKeggGo",
            label = "Fetch pathway genes",
            icon = icon("search")
          ),
          uiOutput("heatmapKeggInfo"),
          hr(),
          h5("Subsetting"),
          conditionalPanel(
            condition = "output.test == 'RNASeq'",
            checkboxInput(
              inputId = "heatmapKeggDEOnly",
              label = "Show only DE genes",
              value = FALSE
            ),
            conditionalPanel(
              condition = "input.heatmapKeggDEOnly == true",
              numericInput(inputId = "heatmapKeggLFC", label = "Log2 Fold Change Threshold:", value = 0, min = 0, max = 10, step = 0.25),
              selectInput(inputId = "heatmapKeggPType", choices = c("FDR", "Raw"), label = "P-Value:", selected = "FDR"),
              selectInput(inputId = "heatmapKeggPThreshold", choices = c(1.0, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), label = "P-Value Threshold:", selected = 0.05)
            )
          ),
          numericInput(
            inputId = "heatmapKeggExpr",
            label = "Minimum median feature expression",
            value = 0,
            min = 0
          ),
          helpText(
            "Removes lowly-expressed genes. The threshold applies to the absolute",
            "median expression across all samples, after any normalisation, log2",
            "transform, or scaling selected in the shared heatmap settings.",
            "Increase to focus on more highly-expressed pathway members."
          )
        ),
        column(
          width = 9,
          withSpinner(
            plotOutput("heatmapKEGG", inline = TRUE),
            color = "#6E8D9E"
          )
        )
      )
    )
  }


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
    "CrisprCountQC" = "Normalised + Log2",
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
    inputId = "fb_select_heatmap",
    choices = inputDataReactive()$genes,
    selected = NULL,
    server = T,
  )

  # Populate exclude dropdown for CRISPR data
  if (inputDataReactive()$dataType == "CrisprCountQC") {
    updateSelectizeInput(session, "heatmapExcludeGuides", choices = inputDataReactive()$genes, selected = NULL, server = TRUE)
  }

  # Update feature choices when collapse guides toggle changes
  observeEvent(input$heatmapCollapseGuides, {
    if (inputDataReactive()$dataType == "CrisprCountQC" && isTRUE(input$heatmapCollapseGuides)) {
      gene_choices <- sort(unique(inputDataReactive()$sgRNA_to_gene))
    } else {
      gene_choices <- inputDataReactive()$genes
    }
    updateSelectizeInput(session, "fb_select_heatmap", choices = gene_choices, selected = "", server = TRUE)
    updateSelectizeInput(session, "heatmapExcludeGuides", choices = gene_choices, selected = NULL, server = TRUE)
  }, ignoreInit = TRUE)

  # Batch + design dropdowns — see server-boxplot.R for the decision table.
  .dt       <- inputDataReactive()$dataType
  .factors  <- inputDataReactive()$factors
  .grouping <- inputDataReactive()$groupingName
  .param    <- inputDataReactive()$param
  .twoFactor <- !is.null(.param$grouping2Name) && nzchar(.param$grouping2Name) &&
                .param$grouping2Name %in% .factors
  .deDesignCol <- switch(.dt, "RNASeq" = .grouping, "proteomics" = "G_", NA_character_)

  batchChoices <- c("None", .factors)
  batchDefault <- if (.twoFactor) .param$grouping2Name else "None"
  updateSelectInput(session = session, inputId = "heatmapBatch", choices = batchChoices, selected = batchDefault)

  designChoices <- if (.dt %in% c("RNASeq", "proteomics") && .twoFactor) {
    c("None", .deDesignCol)
  } else {
    c("None", .factors)
  }
  designDefault <- switch(.dt, "RNASeq" = .grouping, "proteomics" = "G_", "None")
  updateSelectInput(session = session, inputId = "heatmapDesign", choices = designChoices, selected = designDefault)

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
  if (inputDataReactive()$dataType %in% c("RNASeq", "CountQC", "CrisprCountQC")) {
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
    "CountQC" = "factors",
    "CrisprCountQC" = "factors"
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

# Bucket render and reset handler removed — now in server-featureBucket.R

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

# Update KEGG expression threshold step based on count scale
observeEvent({
  input$heatmapCounts
  input$heatmapLog2
  input$heatmapZScore
}, ignoreInit = TRUE, {
  isLogSpace <- grepl("Log2", input$heatmapCounts, ignore.case = TRUE) ||
    (input$heatmapLog2 && input$heatmapCounts %in% c("TPM", "FPKM", "Normalised"))
  isScaled <- input$heatmapZScore %in% c("Z-Score", "Centred")

  stepVal <- if (isScaled || isLogSpace) 0.1 else 1
  updateNumericInput(session = session, inputId = "heatmapKeggExpr", step = stepVal)
  updateNumericInput(session = session, inputId = "heatmapGoExpr", step = stepVal)
})

# Update settings for CountQC data type
observeEvent(inputDataReactive(), {
  if (inputDataReactive()$dataType == "CountQC") {
    # For CountQC, hide gene names by default since there will be many features
    updateCheckboxInput(session = session, inputId = "geneNamesHeatmap", value = FALSE)
    # Set cell border size to 0 for CountQC (will have many features)
    updateSliderInput(session = session, inputId = "cellBorderSizeHeatmap2", value = 0)
  } else if (inputDataReactive()$dataType == "CrisprCountQC") {
    # CRISPR data has fewer features, so show names by default
    updateCheckboxInput(session = session, inputId = "geneNamesHeatmap", value = TRUE)
  }
})

# Dynamic tab title based on data type
output$heatmapTabTitle <- renderUI({
  req(inputDataReactive())
  if (inputDataReactive()$dataType %in% c("CountQC", "CrisprCountQC")) {
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
    input$heatmapDesign,
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
    input$fb_keep_heatmap,
    input$heatmapCustomTitle,
    goHeatmapReactive$currentGenes,
    input$heatmapGoExpr,
    input$heatmapGoDEOnly,
    input$heatmapGoPType,
    input$heatmapGoPThreshold,
    input$heatmapGoLFC,
    keggHeatmapReactive$currentGenes,
    input$heatmapKeggExpr,
    input$heatmapKeggDEOnly,
    input$heatmapKeggPType,
    input$heatmapKeggPThreshold,
    input$heatmapKeggLFC,
    input$heatmapColumnSelect,
    input$cellBorderSizeHeatmap1,
    input$cellBorderColourHeatmap1,
    input$cellBorderSizeHeatmap2,
    input$cellBorderColourHeatmap2,
    input$heatmapCollapseGuides,
    input$heatmapExcludeGuides,
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
  ignoreInit = FALSE,
  {
    # Check if user has removed all factors - if so, clear the plot
    if (is.null(input$heatmapColumnSelect) || length(input$heatmapColumnSelect) == 0 ||
      all(input$heatmapColumnSelect == "")) {
      figuresDataReactive[["heatmapDEFeatures"]] <- NULL
      figuresDataReactive[["heatmapCustom"]] <- NULL
      figuresDataReactive[["heatmapGO"]] <- NULL
      figuresDataReactive[["heatmapKEGG"]] <- NULL
      return(NULL)
    }

    # tryCatch({
    req(input$hmKeepBucket1)
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
      # Collapse sgRNA guides to gene level if requested
      if (inputDataReactive()$dataType == "CrisprCountQC") {
        if (isTRUE(input$heatmapCollapseGuides)) {
          countsHeatmap <- collapse_guides_to_genes(countsHeatmap, inputDataReactive()$sgRNA_to_gene, exclude = input$heatmapExcludeGuides)
        } else if (length(input$heatmapExcludeGuides) > 0) {
          countsHeatmap <- countsHeatmap[!rownames(countsHeatmap) %in% input$heatmapExcludeGuides, , drop = FALSE]
        }
      }
      # Log2 and batch correction on FULL matrix (all genes, all samples) before subsetting
      is_logged <- input$heatmapCounts %in% c("Normalised + Log2", "VST") ||
        inputDataReactive()$dataType == "proteomics"
      if (input$heatmapLog2 & input$heatmapCounts %in% c("TPM", "FPKM", "Normalised")) {
        countsHeatmap <- log2(countsHeatmap + 1)
        is_logged <- TRUE
      }
      if (input$heatmapBatch != "None") {
        if (is_logged) {
          countsHeatmap <- apply_batch_correction(
            counts    = countsHeatmap,
            dataset   = inputDataReactive()$dataset,
            batchCol  = input$heatmapBatch,
            designCol = input$heatmapDesign
          )
        } else {
          showNotification("Batch correction requires log-scale data. Please select a logged count type or enable the log2 transform.",
            type = "warning", duration = 8)
        }
      }

      # Subset to selected samples
      countsHeatmap <- countsHeatmap[, datasetHeatmap$Name]
      datasetHeatmapOriginal <- datasetHeatmap
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
      if (inputDataReactive()$dataType %in% c("CountQC", "CrisprCountQC")) {
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
          "Both" = seqAnnoHeatmap$feature_id[
            which(
              seqAnnoHeatmap[[pTypeHeatmap]] <= as.numeric(input$pThresholdHeatmap) &
                abs(seqAnnoHeatmap$log2Ratio) >= as.numeric(input$lfcHeatmap)
            )
          ],
          "Up-regulated" = seqAnnoHeatmap$feature_id[
            which(
              seqAnnoHeatmap[[pTypeHeatmap]] <= as.numeric(input$pThresholdHeatmap) &
                seqAnnoHeatmap$log2Ratio >= as.numeric(input$lfcHeatmap)
            )
          ],
          "Down-regulated" = seqAnnoHeatmap$feature_id[
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
        if (inputDataReactive()$dataType %in% c("CountQC", "CrisprCountQC")) {
          max_features <- min(as.numeric(input$heatmapGeneNumberCountQC), nrow(countsHeatmapSig))
          countsHeatmapSig <- countsHeatmapSig[1:max_features, ]
        } else {
          if (as.numeric(input$heatmapGeneNumber) <= nrow(countsHeatmapSig)) {
            countsHeatmapSig <- countsHeatmapSig[1:as.numeric(input$heatmapGeneNumber), ]
          }
        }
        rownames(countsHeatmapSig) <- substr(rownames(countsHeatmapSig), 1, 15)

        # Determine heatmap title based on data type
        heatmap_title <- if (inputDataReactive()$dataType %in% c("CountQC", "CrisprCountQC")) {
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
        countsHeatmapSigDL <- as.data.frame(countsHeatmapSig) %>% rownames_to_column("feature_id")
        output[["dlTable_heatmapDE"]] <- downloadHandler(
          filename = function() {
            generate_filename(
              base_name = "DE_heatmap_counts",
              design_name = design,
              extension = "xlsx"
            )
          },
          content = function(file) {
            counts_export <- countsHeatmapSigDL

            # Determine which columns to include
            # Basic columns: feature_id + all sample columns
            basic_cols <- colnames(counts_export)

            # Add DE statistics if requested
            if (inputDataReactive()$dataType %in% c("RNASeq", "proteomics") &&
                !is.null(input$includeDEStatsHeatmap) &&
                input$includeDEStatsHeatmap) {

              if (inputDataReactive()$dataType == "RNASeq") {
                de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
                valid_de_cols <- intersect(de_cols, colnames(inputDataReactive()$seqAnno))
                if (length(valid_de_cols) > 1) {
                  de_data <- inputDataReactive()$seqAnno[, valid_de_cols, drop = FALSE]
                  counts_export <- counts_export %>% left_join(de_data, by = "feature_id")
                }
              } else if (!is.null(input$contrastSelected)) {
                de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
                valid_de_cols <- intersect(de_cols, colnames(inputDataReactive()$seqAnnoList[[input$contrastSelected]]))
                if (length(valid_de_cols) > 1) {
                  de_data <- inputDataReactive()$seqAnnoList[[input$contrastSelected]][, valid_de_cols, drop = FALSE]
                  counts_export <- counts_export %>% left_join(de_data, by = "feature_id")
                }
              }
            } else {
              # Export only basic columns (feature_id + sample columns)
              counts_export <- counts_export[, basic_cols, drop = FALSE]
            }

            export_excel_standard(
              df = counts_export,
              filename = file,
              sheet_name = "Heatmap Data",
              table_name = "HeatmapData"
            )
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
      if (length(input$fb_keep_heatmap) == 1) {
        countsHeatmapCustom <- t(as.matrix(countsHeatmap[input$fb_keep_heatmap, ]))
        rownames(countsHeatmapCustom) <- input$fb_keep_heatmap
      } else {
        countsHeatmapCustom <- countsHeatmap[input$fb_keep_heatmap, ]
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

      output[["dlTable_heatmapCustom"]] <- downloadHandler(
        filename = function() {
          generate_filename(
            base_name = "Custom_heatmap",
            design_name = design,
            extension = "xlsx"
          )
        },
        content = function(file) {
          counts_export <- as.data.frame(countsHeatmapCustom) %>% rownames_to_column("feature_id")

          # Determine which columns to include
          # Basic columns: feature_id + all sample columns
          basic_cols <- colnames(counts_export)

          # Add DE statistics if requested
          if (inputDataReactive()$dataType %in% c("RNASeq", "proteomics") &&
              !is.null(input$includeDEStatsHeatmap) &&
              input$includeDEStatsHeatmap) {

            if (inputDataReactive()$dataType == "RNASeq") {
              de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
              valid_de_cols <- intersect(de_cols, colnames(inputDataReactive()$seqAnno))
              if (length(valid_de_cols) > 1) {
                de_data <- inputDataReactive()$seqAnno[, valid_de_cols, drop = FALSE]
                counts_export <- counts_export %>% left_join(de_data, by = "feature_id")
              }
            } else if (!is.null(input$contrastSelected)) {
              de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
              valid_de_cols <- intersect(de_cols, colnames(inputDataReactive()$seqAnnoList[[input$contrastSelected]]))
              if (length(valid_de_cols) > 1) {
                de_data <- inputDataReactive()$seqAnnoList[[input$contrastSelected]][, valid_de_cols, drop = FALSE]
                counts_export <- counts_export %>% left_join(de_data, by = "feature_id")
              }
            }
          } else {
            # Export only basic columns (feature_id + sample columns)
            counts_export <- counts_export[, basic_cols, drop = FALSE]
          }

          export_excel_standard(
            df = counts_export,
            filename = file,
            sheet_name = "Custom Heatmap",
            table_name = "CustomHeatmap"
          )
        }
      )


      # GO heatmap ----
      if (inputDataReactive()$dataType %in% c("RNASeq", "CountQC", "CrisprCountQC")) {
        goGenes <- goHeatmapReactive$currentGenes
        countsHeatmapGO <- NULL
        goGenesTotal <- length(goGenes)

        if (!is.null(goGenes) && length(goGenes) > 0) {
          goGenes <- goGenes[which(goGenes %in% rownames(countsHeatmap))]
          goGenesTotal <- length(goGenes)

          # Apply DE filter if enabled
          if (!is.null(input$heatmapGoDEOnly) && input$heatmapGoDEOnly) {
            pCol <- if (input$heatmapGoPType == "Raw") "pValue" else "fdr"
            deFeatures <- seqAnnoHeatmap$feature_id[
              which(seqAnnoHeatmap[[pCol]] <= as.numeric(input$heatmapGoPThreshold) &
                    abs(seqAnnoHeatmap$log2Ratio) >= as.numeric(input$heatmapGoLFC))
            ]
            goGenes <- goGenes[goGenes %in% deFeatures]
          }

          if (length(goGenes) == 1) {
            countsHeatmapGO <- t(as.matrix(countsHeatmap[goGenes, ]))
            rownames(countsHeatmapGO) <- goGenes
          } else if (length(goGenes) > 1) {
            countsHeatmapGO <- countsHeatmap[goGenes, ]
            countsHeatmapGO <- countsHeatmapGO[abs(rowMedians(countsHeatmapGO)) >= input$heatmapGoExpr, ]
          }
        }

        if (!is.null(countsHeatmapGO) && nrow(countsHeatmapGO) > 0) {
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
          column_title = paste0(
            goHeatmapReactive$currentGoName %||% goHeatmapReactive$currentGoTerm %||% "",
            "\n", nrow(countsHeatmapGO), "/", goGenesTotal, " features plotted",
            if (agg_label != "") paste0("\n", agg_label) else ""
          ),
          heatmap_legend_param = list(labels_gp = gpar(fontsize = input$textSizeHeatmap), title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"))
        )
        figuresDataReactive$heatmapGO <- gh
        output[["dlTable_heatmapGO"]] <- downloadHandler(
          filename = function() {
            generate_filename(
              base_name = "GO_heatmap",
              design_name = design,
              extension = "xlsx"
            )
          },
          content = function(file) {
            counts_export <- as.data.frame(countsHeatmapGO) %>% rownames_to_column("feature_id")

            # Determine which columns to include
            # Basic columns: feature_id + all sample columns
            basic_cols <- colnames(counts_export)

            # Add DE statistics if requested
            if (inputDataReactive()$dataType %in% c("RNASeq", "proteomics") &&
                !is.null(input$includeDEStatsHeatmap) &&
                input$includeDEStatsHeatmap) {

              if (inputDataReactive()$dataType == "RNASeq") {
                de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
                valid_de_cols <- intersect(de_cols, colnames(inputDataReactive()$seqAnno))
                if (length(valid_de_cols) > 1) {
                  de_data <- inputDataReactive()$seqAnno[, valid_de_cols, drop = FALSE]
                  counts_export <- counts_export %>% left_join(de_data, by = "feature_id")
                }
              } else if (!is.null(input$contrastSelected)) {
                de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
                valid_de_cols <- intersect(de_cols, colnames(inputDataReactive()$seqAnnoList[[input$contrastSelected]]))
                if (length(valid_de_cols) > 1) {
                  de_data <- inputDataReactive()$seqAnnoList[[input$contrastSelected]][, valid_de_cols, drop = FALSE]
                  counts_export <- counts_export %>% left_join(de_data, by = "feature_id")
                }
              }
            } else {
              # Export only basic columns (feature_id + sample columns)
              counts_export <- counts_export[, basic_cols, drop = FALSE]
            }

            export_excel_standard(
              df = counts_export,
              filename = file,
              sheet_name = "GO Heatmap",
              table_name = "GOHeatmap"
            )
          }
        )
        } else {
          figuresDataReactive$heatmapGO <- NULL
        }
      }

      # KEGG heatmap ----
      if (inputDataReactive()$dataType %in% c("RNASeq", "CountQC", "CrisprCountQC") &&
          !is.null(keggHeatmapReactive$keggCode)) {
        keggGenes <- keggHeatmapReactive$currentGenes
        countsHeatmapKEGG <- NULL

        if (!is.null(keggGenes) && length(keggGenes) > 0) {
          keggGenes <- keggGenes[which(keggGenes %in% rownames(countsHeatmap))]

          # Apply DE filter if enabled
          if (!is.null(input$heatmapKeggDEOnly) && input$heatmapKeggDEOnly) {
            pCol <- if (input$heatmapKeggPType == "Raw") "pValue" else "fdr"
            deFeatures <- seqAnnoHeatmap$feature_id[
              which(seqAnnoHeatmap[[pCol]] <= as.numeric(input$heatmapKeggPThreshold) &
                    abs(seqAnnoHeatmap$log2Ratio) >= as.numeric(input$heatmapKeggLFC))
            ]
            keggGenes <- keggGenes[keggGenes %in% deFeatures]
          }

          if (length(keggGenes) == 1) {
            countsHeatmapKEGG <- t(as.matrix(countsHeatmap[keggGenes, ]))
            rownames(countsHeatmapKEGG) <- keggGenes
          } else if (length(keggGenes) > 1) {
            countsHeatmapKEGG <- countsHeatmap[keggGenes, ]
            if (!is.null(input$heatmapKeggExpr)) {
              countsHeatmapKEGG <- countsHeatmapKEGG[abs(rowMedians(countsHeatmapKEGG)) >= input$heatmapKeggExpr, ]
            }
          }

          keggGenesTotal <- length(keggHeatmapReactive$currentGenes)

          if (!is.null(countsHeatmapKEGG) && nrow(countsHeatmapKEGG) > 0) {
            req(nrow(datasetHeatmap) == ncol(countsHeatmapKEGG))

            keggTitle <- paste0(
              keggHeatmapReactive$currentPathwayName %||% "",
              "\n", nrow(countsHeatmapKEGG), "/", keggGenesTotal, " features plotted",
              if (agg_label != "") paste0("\n", agg_label) else ""
            )

            kh <- ComplexHeatmap::Heatmap(
              matrix = countsHeatmapKEGG,
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
              column_title = keggTitle,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = input$textSizeHeatmap),
                title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold")
              )
            )
            figuresDataReactive$heatmapKEGG <- kh

            output[["dlTable_heatmapKEGG"]] <- downloadHandler(
              filename = function() {
                generate_filename(
                  base_name = "KEGG_heatmap",
                  design_name = design,
                  extension = "xlsx"
                )
              },
              content = function(file) {
                counts_export <- as.data.frame(countsHeatmapKEGG) %>% rownames_to_column("feature_id")
                basic_cols <- colnames(counts_export)

                if (inputDataReactive()$dataType %in% c("RNASeq", "proteomics") &&
                    !is.null(input$includeDEStatsHeatmap) &&
                    input$includeDEStatsHeatmap) {
                  if (inputDataReactive()$dataType == "RNASeq") {
                    de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
                    valid_de_cols <- intersect(de_cols, colnames(inputDataReactive()$seqAnno))
                    if (length(valid_de_cols) > 1) {
                      de_data <- inputDataReactive()$seqAnno[, valid_de_cols, drop = FALSE]
                      counts_export <- counts_export %>% left_join(de_data, by = "feature_id")
                    }
                  } else if (!is.null(input$contrastSelected)) {
                    de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
                    valid_de_cols <- intersect(de_cols, colnames(inputDataReactive()$seqAnnoList[[input$contrastSelected]]))
                    if (length(valid_de_cols) > 1) {
                      de_data <- inputDataReactive()$seqAnnoList[[input$contrastSelected]][, valid_de_cols, drop = FALSE]
                      counts_export <- counts_export %>% left_join(de_data, by = "feature_id")
                    }
                  }
                } else {
                  counts_export <- counts_export[, basic_cols, drop = FALSE]
                }

                export_excel_standard(
                  df = counts_export,
                  filename = file,
                  sheet_name = "KEGG Heatmap",
                  table_name = "KEGGHeatmap"
                )
              }
            )
          } else {
            figuresDataReactive$heatmapKEGG <- NULL
          }
        }
      }
    # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
)

# GO heatmap: fetch GO term genes on button click
observeEvent(input$heatmapGoGo, {
  req(input$heatmapGoInput)
  req(inputDataReactive()$seqAnno)

  goId <- gsub(" .*", "", input$heatmapGoInput)
  goName <- sub("^\\S+\\s+", "", input$heatmapGoInput)
  seqAnno <- inputDataReactive()$seqAnno

  goGenes <- seqAnno$feature_id[with(seqAnno, grepl(goId, paste(`GO BP`, `GO MF`, `GO CC`)))]

  goHeatmapReactive$currentGoTerm <- goId
  goHeatmapReactive$currentGoName <- goName
  goHeatmapReactive$currentGenes <- goGenes
})

# GO heatmap: display GO term info
output$heatmapGoInfo <- renderUI({
  req(goHeatmapReactive$currentGoTerm)
  tagList(
    hr(),
    tags$b("GO term: "), goHeatmapReactive$currentGoName,
    br(),
    tags$b("Genes found in dataset: "), length(goHeatmapReactive$currentGenes)
  )
})

# KEGG heatmap: fetch pathway genes on button click
observeEvent(input$heatmapKeggGo, {
  req(input$heatmapKeggInput)
  req(keggHeatmapReactive$keggCode)

  keggCode <- keggHeatmapReactive$keggCode

  # Normalise the pathway ID input
  pathwayInput <- trimws(input$heatmapKeggInput)
  if (grepl("^[0-9]", pathwayInput)) {
    pathwayId <- paste0(keggCode, pathwayInput)
  } else {
    pathwayId <- pathwayInput
  }

  # Check cache
  if (pathwayId %in% names(keggHeatmapReactive$cache)) {
    cached <- keggHeatmapReactive$cache[[pathwayId]]
    keggHeatmapReactive$currentGenes <- cached$genes
    keggHeatmapReactive$currentPathwayName <- cached$name
    return()
  }

  # Fetch from KEGG API
  result <- withProgress(message = "Fetching KEGG pathway genes...", {
    tryCatch({
      # Get pathway name
      pathInfo <- KEGGREST::keggGet(pathwayId)
      pathwayName <- sub(" - .*$", "", pathInfo[[1]]$NAME)

      # Get gene list (returns ENTREZ IDs with species prefix)
      geneLinks <- KEGGREST::keggLink(keggCode, paste0("path:", pathwayId))
      entrezIds <- gsub(paste0(keggCode, ":"), "", geneLinks)

      # Match ENTREZ IDs to feature_ids in seqAnno
      seqAnno <- inputDataReactive()$seqAnno

      if ("entrez_id" %in% colnames(seqAnno)) {
        # Primary strategy: direct match on entrez_id column
        keggFeatureIds <- seqAnno$feature_id[seqAnno$entrez_id %in% entrezIds]
      } else {
        # Fallback: convert via pathview
        idMapping <- pathview::eg2id(entrezIds, category = c("SYMBOL", "ENSEMBL"), org = keggCode)
        geneSymbols <- na.omit(idMapping[, "SYMBOL"])
        ensemblIds <- na.omit(idMapping[, "ENSEMBL"])
        matchByName <- seqAnno$feature_id[seqAnno$gene_name %in% geneSymbols]
        matchById <- seqAnno$feature_id[seqAnno$gene_id %in% ensemblIds]
        keggFeatureIds <- unique(c(matchByName, matchById))
      }

      list(genes = keggFeatureIds, name = pathwayName)
    }, error = function(e) {
      showNotification(
        paste("Failed to fetch KEGG pathway:", e$message),
        type = "error", duration = 8
      )
      NULL
    })
  })

  if (!is.null(result)) {
    keggHeatmapReactive$cache[[pathwayId]] <- result
    keggHeatmapReactive$currentGenes <- result$genes
    keggHeatmapReactive$currentPathwayName <- result$name
  } else {
    keggHeatmapReactive$currentGenes <- NULL
    keggHeatmapReactive$currentPathwayName <- NULL
  }
})

# KEGG heatmap: display pathway info
output$heatmapKeggInfo <- renderUI({
  req(keggHeatmapReactive$currentPathwayName)
  tagList(
    hr(),
    tags$b("Pathway: "), keggHeatmapReactive$currentPathwayName,
    br(),
    tags$b("Genes found in dataset: "), length(keggHeatmapReactive$currentGenes)
  )
})

# Output heatmaps
lapply(c("DEFeatures", "Custom", "GO", "KEGG"), function(sig) {
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
  
  # Map 'sig' to button ID suffix
  # sig: "DEFeatures" -> ID: "dlPlot_heatmapDE"
  # sig: "Custom"     -> ID: "dlPlot_heatmapCustom"
  # sig: "GO"         -> ID: "dlPlot_heatmapGO"
  button_suffix <- switch(sig,
                          "DEFeatures" = "DE",
                          "Custom" = "Custom",
                          "GO" = "GO",
                          "KEGG" = "KEGG")
  
  output[[paste0("dlPlot_heatmap", button_suffix)]] <- downloadHandler(
    filename = function() {
      generate_filename(
        base_name = paste0("heatmap_", sig),
        design_name = input$heatmapFilename,
        extension = tolower(input$heatmapDownloadFormat)
      )
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
      # Using a custom draw wrapper because export_plot_standard expects a print/plot call, 
      # and ComplexHeatmap needs draw() with specific params.
      export_plot_standard(
        plot_object = function() {
          draw(heatmap_obj, merge_legend = TRUE, padding = unit(c(10, 20, 10, 10), "mm"))
        },
        filename = file,
        format = input$heatmapDownloadFormat,
        dpi = input$heatmapDPI,
        width = as.numeric(input$figWidthHeatmap / 60),
        height = as.numeric(input$figHeightHeatmap / 60)
      )
    }
  )
})
