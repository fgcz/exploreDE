eds <- inputDataReactive()$.eds
req(eds)

# Main reactive expression for ordination calculations
ordinationResults <- reactive({
  req(eds)

  # Comprehensive validation before proceeding
  req(
    !is.null(input$pcaFactor1),
    input$pcaFactor1 != "None",
    !is.null(input$pcaGroups),
    length(input$pcaGroups) > 0,
    input$pcaCounts %in% eds$assay_names(),
    !is.null(input$ordinationMethod)
  )

  # Validate selected groups exist in the dataset
  datasetOrd <- eds$samples()
  available_groups <- levels(as.factor(datasetOrd[[input$pcaFactor1]]))
  valid_groups <- input$pcaGroups[input$pcaGroups %in% available_groups]
  req(length(valid_groups) > 0)

  # Prepare data: start with full dataset and counts
  # GLM-PCA requires raw counts; other methods use the user-selected count type
  if (input$ordinationMethod == "glmpca") {
    if (!eds$has_assay("counts")) {
      showNotification("GLM-PCA requires raw count data which is not available", type = "error", duration = 5)
      return(NULL)
    }
    countsOrd <- eds$assay("counts")
  } else {
    countsOrd <- eds$assay(input$pcaCounts)
  }
  # Collapse sgRNA guides to gene level if requested
  if (eds$data_type() == "CrisprCountQC" && isTRUE(input$pcaCollapseGuides)) {
    countsOrd <- collapse_guides_to_genes(countsOrd, eds$sgrna_to_gene())
  }

  # Build cascade step map
  # Track whether counts are in log space (needed for batch correction)
  batch_col <- if (input$pcaBatch != "None") input$pcaBatch else NULL
  design_col <- input$pcaBatchDesign # read here so ordinationResults() tracks it as a dep even if the cascade step short-circuits
  is_logged <- eds$assay_is_log2(input$pcaCounts)

  # Features before and after the peptide filter, for the message below
  peptide_filter <- NULL
  step_map <- list(
    "Group subsetting" = function(cd) step_group_subset(cd$counts, cd$dataset, input$pcaFactor1, valid_groups),
    "Peptide filter" = function(cd) {
      if (eds$omics_class() == "proteomics" && length(eds$contrasts())) {
        res <- step_peptide_filter(cd$counts, cd$dataset, eds$de(eds$contrasts()[[1]]), input$nrPeptidePCA)
        peptide_filter <<- c(before = nrow(cd$counts), after = nrow(res$counts))
        res
      } else {
        cd
      }
    },
    "Log2 transform" = function(cd) {
      if (input$ordinationMethod == "glmpca") {
        return(cd) # GLM-PCA requires untransformed counts
      }
      prep <- prepare_expression(cd$counts, cd$dataset, is_log2 = is_logged, take_log2 = input$pcaLog2)
      is_logged <<- prep$is_log2
      list(counts = prep$counts, dataset = cd$dataset)
    },
    "Batch correction" = function(cd) {
      if (input$ordinationMethod == "glmpca") {
        return(cd) # GLM-PCA handles batch via X covariate, not removeBatchEffect
      }
      prep <- prepare_expression(
        cd$counts,
        cd$dataset,
        is_log2 = is_logged,
        batch = input$pcaBatch,
        design = design_col
      )
      list(counts = prep$counts, dataset = cd$dataset)
    },
    "Top N features (SD)" = function(cd) step_top_n_sd(cd$counts, cd$dataset, input$pcaTopN)
  )

  # Execute cascade in user-defined order
  cascade_order <- input$pcaCascadeOrder
  if (is.null(cascade_order) || length(cascade_order) == 0) {
    cascade_order <- c(
      "Peptide filter",
      "Log2 transform",
      "Batch correction",
      "Group subsetting",
      "Top N features (SD)"
    )
  }

  cd <- list(counts = countsOrd, dataset = datasetOrd)
  for (step_name in cascade_order) {
    if (step_name %in% names(step_map)) {
      cd <- step_map[[step_name]](cd)
    }
  }
  countsOrd <- cd$counts
  datasetOrd <- cd$dataset

  # Transpose for PCA/PCoA (GLM-PCA keeps genes x samples)
  if (input$ordinationMethod != "glmpca") {
    countsPrepared <- as.data.frame(t(countsOrd))
  } else {
    countsPrepared <- countsOrd
  }

  if (is.null(countsPrepared)) {
    showNotification(
      "Unable to prepare count data for ordination. Please check your data and settings.",
      type = "error",
      duration = 5
    )
    return(NULL)
  }

  # Enough features for ordination? Said in the plot area, naming the peptide
  # filter when it removed features.
  peptide_note <- if (!is.null(peptide_filter) && peptide_filter[["after"]] < peptide_filter[["before"]]) {
    paste0(
      " 'Minimum number of peptides' (",
      input$nrPeptidePCA,
      ") keeps ",
      peptide_filter[["after"]],
      " of ",
      peptide_filter[["before"]],
      " features."
    )
  }
  validate(need(
    nrow(countsOrd) >= 3,
    paste0("Ordination needs at least 3 features; these settings leave ", nrow(countsOrd), ".", peptide_note)
  ))

  # Calculate ordination based on selected method
  ord_result <- NULL

  if (input$ordinationMethod == "pca") {
    # Require PCA-specific inputs
    req(!is.null(input$pcaCentre), !is.null(input$pcaScale))

    ord_result <- calculate_pca(
      counts = countsPrepared,
      center = input$pcaCentre,
      scale = input$pcaScale
    )
  } else if (input$ordinationMethod == "pcoa") {
    # Require PCoA-specific inputs
    req(!is.null(input$pcoaDistance))

    ord_result <- calculate_pcoa(
      counts = countsPrepared,
      method = input$pcoaDistance,
      nmds = FALSE
    )
  } else if (input$ordinationMethod == "nmds") {
    # Require NMDS-specific inputs
    req(!is.null(input$pcoaDistance))

    ord_result <- calculate_pcoa(
      counts = countsPrepared,
      method = input$pcoaDistance,
      nmds = TRUE
    )
  } else if (input$ordinationMethod == "glmpca") {
    # Require GLM-PCA-specific inputs with strict validation
    req(
      !is.null(input$glmpcaFamily),
      is.character(input$glmpcaFamily),
      nchar(input$glmpcaFamily) > 0,
      input$glmpcaFamily %in% c("poi", "nb", "binom")
    )

    # GLM-PCA uses the cascade-filtered countsPrepared (genes x samples, raw counts)
    # The cascade already applied group subsetting, peptide filter, top N etc.
    # but skipped log2 and batch correction (handled here natively)
    if (is.null(countsPrepared) || nrow(countsPrepared) < 3 || ncol(countsPrepared) < 3) {
      showNotification("GLM-PCA: not enough data after preprocessing", type = "error", duration = 5)
      return(NULL)
    }

    # Build batch covariate matrix for GLM-PCA (instead of removeBatchEffect)
    X_covariates <- NULL
    if (!is.null(batch_col) && batch_col %in% colnames(datasetOrd)) {
      batch_factor <- as.factor(datasetOrd[[batch_col]])
      if (nlevels(batch_factor) > 1) {
        X_covariates <- model.matrix(~batch_factor)[, -1, drop = FALSE]
      }
    }

    # Show waiter during GLM-PCA computation
    waiter <- waiter::Waiter$new(
      id = "pcaStatic",
      html = tagList(
        waiter::spin_fading_circles(),
        h4("Computing GLM-PCA..."),
        p("This may take 1-2 minutes for large datasets")
      ),
      color = "#527085"
    )
    waiter$show()
    on.exit(waiter$hide(), add = TRUE)

    ord_result <- calculate_glmpca(
      counts = countsPrepared,
      family = input$glmpcaFamily,
      L = 2,
      X = X_covariates
    )
  }

  if (is.null(ord_result)) {
    return(NULL)
  }

  # Add dataset information to scores
  ord_result$scores$Name <- as.character(ord_result$scores$Name)
  ord_result$scores <- left_join(ord_result$scores, datasetOrd, by = "Name")
  ord_result$dataset <- datasetOrd
  ord_result$counts_used <- countsPrepared

  return(ord_result)
})

observe({
  # Build labelled choices for the count selector: key → human-readable label
  assay_keys <- eds$assay_names()
  count_choices <- setNames(assay_keys, vapply(assay_keys, eds$assay_label, character(1)))

  output$pcaDesign <- renderText({
    designReactive()
  })

  if (eds$data_type() == "RNASeq") {
    default_count <- if ("normalised_log2" %in% assay_keys) "normalised_log2" else assay_keys[1]
    updateSelectInput(session = session, inputId = "pcaCounts", choices = count_choices, selected = default_count)
  } else if (eds$omics_class() == "proteomics") {
    default_count <- if ("intensities" %in% assay_keys) "intensities" else assay_keys[1]
    updateSelectInput(session = session, inputId = "pcaCounts", choices = count_choices, selected = default_count)
  } else if (eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
    default_count <- if ("vst" %in% assay_keys) "vst" else assay_keys[1]
    updateSelectInput(session = session, inputId = "pcaCounts", choices = count_choices, selected = default_count)
  }
  # Batch + design dropdowns: batch_design_defaults() in R/utils_explore_dataset.R.
  # NB: the design input is `pcaBatchDesign` because `pcaDesign` is already a
  # textOutput on this tab (the DE design-formula string).
  bd <- batch_design_defaults(eds)
  updateSelectInput(session = session, inputId = "pcaBatch", choices = bd$batch_choices, selected = bd$batch_default)
  updateSelectInput(
    session = session,
    inputId = "pcaBatchDesign",
    choices = bd$design_choices,
    selected = bd$design_default
  )

  # Render cascade steps based on data type, log2 checkbox, and count method
  # Only available to staff (LDAP employee role)
  output$pcaCascadeUI <- renderUI({
    if (!identical(ldapRole(), "employee")) {
      return(helpText("Advanced preprocessing settings are available to FGCZ staff only."))
    }
    log2_checked <- !is.null(input$pcaLog2) && input$pcaLog2
    is_linear <- isTRUE(input$pcaCounts %in% eds$assay_names()) && !eds$assay_is_log2(input$pcaCounts)
    show_log2 <- log2_checked && is_linear

    steps <- c("Batch correction", "Group subsetting", "Top N features (SD)")
    if (eds$omics_class() == "proteomics") {
      steps <- c("Peptide filter", "Batch correction", "Group subsetting", "Top N features (SD)")
    }
    if (show_log2) {
      # Insert Log2 before Batch correction
      batch_pos <- which(steps == "Batch correction")
      steps <- append(steps, "Log2 transform", after = batch_pos - 1)
    }

    sortable::rank_list(
      text = "Preprocessing order",
      labels = as.list(steps),
      input_id = "pcaCascadeOrder",
      options = sortable::sortable_options(multiDrag = FALSE)
    )
  })

  lapply(c(1:2), function(i) {
    selected <- switch(i, "1" = eds$primary_grouping(), "2" = "None")
    updateSelectInput(
      session = session,
      inputId = paste0("pcaFactor", i),
      choices = c("None", eds$factor_names()),
      selected = selected
    )
  })
  observeEvent(
    {
      input$pcaFactor1
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE,
    {
      levs <- levels(as.factor(eds$samples()[[input$pcaFactor1]]))
      updateCheckboxGroupInput(
        session = session,
        inputId = "pcaGroups",
        choices = levs,
        # save-state restore: seed the saved group subset (valid= drops levels
        # that don't belong to the restored factor).
        selected = restore_default("pcaGroups", levs, valid = levs)
      )
    }
  )
  if (eds$omics_class() == "proteomics") {
    # Handle nrPeptides slider — check the primary contrast's DE table for nrPeptides column
    if (
      length(eds$contrasts()) > 0 &&
        "nrPeptides" %in% colnames(eds$de(eds$contrasts()[[1]]))
    ) {
      output$nrPeptidesPCAUI <- renderUI({
        sliderInput(
          inputId = "nrPeptidePCA",
          label = "Minimum number of peptides",
          min = 0,
          max = 20,
          value = 2,
          step = 1,
          width = "85%"
        )
      })
      outputOptions(output, "nrPeptidesPCAUI", suspendWhenHidden = FALSE)
    } else {
      output$nrPeptidesPCAUI <- renderUI({
        NULL
      })
    }

    # Update pcaTopN for proteomics data — cap at number of features in dataset
    n_features <- eds$n_features()
    default_value <- min(500, n_features)
    updateNumericInput(
      session = session,
      inputId = "pcaTopN",
      value = default_value,
      max = n_features
    )
  }
})

# Show notification when GLM-PCA is selected
observeEvent(input$ordinationMethod, {
  if (!is.null(input$ordinationMethod) && input$ordinationMethod == "glmpca") {
    showNotification(
      "GLM-PCA computation may take 1-2 minutes for large datasets. Please be patient.",
      type = "message",
      duration = 8,
      closeButton = TRUE
    )
  }
})

# Separate observer to update biplot gene choices only when ordination results change
# This prevents resetting the selection when user picks genes
observeEvent(
  list(ordinationResults(), input$fb_label_by_ordination),
  {
    ord_res <- ordinationResults()
    if (!is.null(ord_res) && !is.null(ord_res$loadings)) {
      # Values stay feature_ids; displayed names follow "Label features by"
      fids <- ord_res$loadings$gene
      updateSelectizeInput(
        session = session,
        inputId = "biplotCustomGenes",
        choices = setNames(fids, label_features(fids, input$fb_label_by_ordination)),
        selected = isolate(input$biplotCustomGenes),
        server = TRUE
      )
    }
  },
  ignoreNULL = FALSE
)

observeEvent(
  {
    input$pcaFactor1
    input$pcaX
    input$pcaY
    input$pcaTopN
    input$pcaCascadeOrder
    input$pcaGroups
    input$pcaFactor2
    input$pcaCounts
    input$pcaAxesProp
    input$pcaBatch
    input$pcaBatchDesign
    input$pcaCentre
    input$pcaScale
    input$textSizePCA
    input$showLinesPCA
    input$showAxesPCA
    input$boldPCA
    input$pcaShowNames
    input$pcaLog2
    input$dotSizePCA
    input$alphaPCA
    input$dotBorderPCA
    input$pcaLabelMaxOverlap
    input$geneLabelNudgePCAX
    input$geneLabelNudgePCAY
    input$geneLabelSizePCA
    input$nrPeptidePCA
    input$pcaAddEllipses
    input$pcaEllipsesAlpha
    # New ordination-specific inputs
    input$ordinationMethod
    input$ordinationPlotType
    input$pcoaDistance
    input$glmpcaFamily
    input$biplotGeneSelectionMethod
    input$biplotNGenes
    input$biplotNGenesPerPC
    input$biplotCustomGenes
    input$biplotCustomGenesText
    input$fb_label_by_ordination
    input$biplotArrowScale
    input$biplotArrowColour
    input$biplotArrowColourX
    input$biplotArrowColourY
    # Add all color picker + level label inputs as dependencies so plot updates when they change
    lapply(eds$factor_names(), function(fname) {
      lapply(eds$factor_levels(fname), function(lev) {
        list(input[[paste0("GroupColour", fname, lev)]], input[[paste0("GroupLabel", fname, lev)]])
      })
    })
  },
  ignoreInit = TRUE,
  ignoreNULL = FALSE,
  {
    # Get ordination results using the new reactive
    ord_res <- ordinationResults()

    # Silently stop if ordination results aren't ready yet (prevents brief error flash)
    req(ord_res)

    # Get valid groups and colors
    available_groups <- levels(as.factor(eds$samples()[[input$pcaFactor1]]))
    valid_groups <- input$pcaGroups[input$pcaGroups %in% available_groups]

    coloursPCA <- setNames(
      lapply(valid_groups, function(k) {
        paste0(col2hex(input[[paste0("GroupColour", input$pcaFactor1, k)]]), "FF")
      }),
      valid_groups
    )

    # Display labels from the Colours tab (fill factor, and shape factor).
    # If either is rich, both are, as they share the markdown legend text.
    pcaShapeFactor <- if (isTRUE(input$pcaFactor2 != "None")) input$pcaFactor2 else NULL
    pcaShapeLevels <- if (is.null(pcaShapeFactor)) NULL else levels(as.factor(eds$samples()[[pcaShapeFactor]]))
    pcaLevelLabels <- list(
      fill = level_labels(input, input$pcaFactor1, valid_groups),
      shape = if (is.null(pcaShapeFactor)) NULL else level_labels(input, pcaShapeFactor, pcaShapeLevels)
    )
    pcaLevelLabels$rich <- is_rich_labels(pcaLevelLabels$fill) || is_rich_labels(pcaLevelLabels$shape)
    if (pcaLevelLabels$rich) {
      pcaLevelLabels$fill <- level_labels(input, input$pcaFactor1, valid_groups, force_rich = TRUE)
      if (!is.null(pcaShapeFactor)) {
        pcaLevelLabels$shape <- level_labels(input, pcaShapeFactor, pcaShapeLevels, force_rich = TRUE)
      }
    }

    # Extract results from ordination
    pc_scores <- ord_res$scores
    pc_loadings <- ord_res$loadings

    # Create variance/stress table
    if (!is.null(ord_res$variance)) {
      # For PCA, PCoA (methods with variance explained)
      pc_eigenvalues <- tibble(
        PC = ord_res$component_names,
        variance = ord_res$variance
      ) %>%
        mutate(pct = (variance / sum(variance)) * 100) %>%
        mutate(pct_cum = cumsum(pct))
      pc_eigenvalues$PC <- factor(pc_eigenvalues$PC, levels = pc_eigenvalues$PC)
    } else if (!is.null(ord_res$stress)) {
      # For NMDS (has stress instead of variance)
      pc_eigenvalues <- tibble(
        PC = ord_res$component_names,
        stress = ord_res$stress
      )
    } else {
      # For GLM-PCA (no variance or stress metrics available)
      pc_eigenvalues <- tibble(
        PC = ord_res$component_names
      )
      pc_eigenvalues$PC <- factor(pc_eigenvalues$PC, levels = pc_eigenvalues$PC)
    }

    # Update axis selectors with smart defaults
    # If the old selection doesn't exist in new choices, reset to first two components
    new_choices <- as.character(pc_eigenvalues$PC)
    x_selected <- if (!is.null(input$pcaX) && input$pcaX %in% new_choices) input$pcaX else new_choices[1]
    y_selected <- if (!is.null(input$pcaY) && input$pcaY %in% new_choices) {
      input$pcaY
    } else {
      new_choices[min(2, length(new_choices))]
    }

    updateSelectInput(session = session, inputId = "pcaX", choices = new_choices, selected = x_selected)
    updateSelectInput(session = session, inputId = "pcaY", choices = new_choices, selected = y_selected)

    # Store results in reactive for safe access by download handlers and brush table
    figuresDataReactive$pc_loadings <- pc_loadings
    figuresDataReactive$pc_eigenvalues <- pc_eigenvalues
    figuresDataReactive$pc_scores <- pc_scores
    figuresDataReactive$ord_method <- ord_res$method
    figuresDataReactive$x_selected <- x_selected
    figuresDataReactive$y_selected <- y_selected
    figuresDataReactive$factors <- eds$factor_names()

    # Update plot type choices based on whether biplot is available
    if (!is.null(pc_loadings)) {
      # Biplot available - offer both options
      updateRadioButtons(
        session = session,
        inputId = "ordinationPlotType",
        choices = c("Samples only" = "samples", "Biplot (samples + genes)" = "biplot"),
        selected = if (!is.null(input$ordinationPlotType)) input$ordinationPlotType else "samples"
      )
    } else {
      # No loadings - biplot not available, force samples only
      updateRadioButtons(
        session = session,
        inputId = "ordinationPlotType",
        choices = c("Samples only" = "samples", "Biplot not available for this method" = "biplot_disabled"),
        selected = "samples"
      )
    }

    # Collect manual genes from dropdown and text area (for biplot)
    # Note: Gene choices are updated in a separate observer to prevent selection reset
    biplot_custom_genes <- NULL
    if (!is.null(input$biplotGeneSelectionMethod) && input$biplotGeneSelectionMethod == "manual") {
      # Get genes from dropdown
      genes_from_dropdown <- input$biplotCustomGenes

      # Get genes from text area (parse by various delimiters)
      genes_from_text <- NULL
      if (!is.null(input$biplotCustomGenesText) && nchar(input$biplotCustomGenesText) > 0) {
        genes_from_text <- strsplit(input$biplotCustomGenesText, "\\s+|/+|,+|\n+|\r+|\t+")[[1]]
        genes_from_text <- genes_from_text[genes_from_text != ""]
        # Proteomics: also accept pasted gene symbols / protein IDs
        genes_from_text <- unique(c(
          genes_from_text,
          fb_resolve(genes_from_text, "gene_name"),
          fb_resolve(genes_from_text, "protein_Id")
        ))
      }

      # Combine and deduplicate
      biplot_custom_genes <- unique(c(genes_from_dropdown, genes_from_text))

      # Validate against available genes
      if (!is.null(pc_loadings) && length(biplot_custom_genes) > 0) {
        biplot_custom_genes <- biplot_custom_genes[biplot_custom_genes %in% pc_loadings$gene]
      }
    }

    # Output the scree plot
    output$pcaScree <- renderPlot(
      {
        pc_eigenvalues2 <- pc_eigenvalues
        pc_eigenvalues2$PC <- gsub("PC|Axis|MDS|dim", "", pc_eigenvalues2$PC)
        pc_eigenvalues2$PC <- factor(pc_eigenvalues2$PC, levels = pc_eigenvalues2$PC)
        if (nrow(pc_eigenvalues2) > 15) {
          pc_eigenvalues2 <- pc_eigenvalues2[1:15, ]
        }

        if (!is.null(ord_res$variance)) {
          # Scree plot for variance-based methods (PCA, PCoA)
          pcaScree <- pc_eigenvalues2 %>%
            ggplot(aes(x = PC)) +
            geom_col(aes(y = pct)) +
            geom_line(aes(y = pct_cum, group = 1)) +
            geom_point(aes(y = pct_cum)) +
            labs(x = "Component", y = "Fraction variance explained") +
            theme_classic(base_size = as.numeric(input$textSizePCA))
        } else if (!is.null(ord_res$stress)) {
          # For NMDS, show stress as a single bar
          pcaScree <- ggplot(data = pc_eigenvalues2[1, ], aes(x = "Stress", y = stress)) +
            geom_col() +
            labs(x = "", y = "Stress") +
            theme_classic(base_size = as.numeric(input$textSizePCA))
        } else {
          # For GLM-PCA, show component list (no variance metrics available)
          pcaScree <- ggplot(data = pc_eigenvalues2, aes(x = PC, y = 1)) +
            geom_col(fill = "steelblue") +
            labs(x = "Component", y = "", title = "GLM-PCA Components") +
            theme_classic(base_size = as.numeric(input$textSizePCA)) +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        }
        pcaScree
      },
      width = 400,
      height = 400
    )

    # Output the variance table
    output$pcaVars <- DT::renderDataTable({
      if (!is.null(ord_res$variance)) {
        datatable(
          pc_eigenvalues[, c("PC", "pct", "pct_cum")],
          rownames = FALSE,
          colnames = c("Component", "Variance explained", "Cumulative variance")
        ) %>%
          formatRound(c("pct", "pct_cum"), digits = 3)
      } else if (!is.null(ord_res$stress)) {
        # For NMDS
        datatable(pc_eigenvalues[, c("PC", "stress")], rownames = FALSE, colnames = c("Component", "Stress"))
      } else {
        # For GLM-PCA (no variance metrics)
        datatable(pc_eigenvalues[, "PC", drop = FALSE], rownames = FALSE, colnames = c("Component"))
      }
    })

    # Output the loading table (only for methods with loadings)
    if (!is.null(pc_loadings) && nrow(pc_loadings) > 0) {
      # Use x_selected and y_selected (already validated) instead of input$pcaX/pcaY
      # This prevents crash when switching methods (input values lag behind updateSelectInput)
      top_genes <- pc_loadings %>%
        dplyr::select(all_of(c("gene", x_selected, y_selected))) %>%
        pivot_longer(matches("PC|Axis|MDS|dim|PCo"), names_to = "Component", values_to = "loading") %>%
        group_by(Component) %>%
        arrange(desc(abs(loading)))

      output$pcaLoadings <- DT::renderDataTable({
        datatable(top_genes, rownames = FALSE) %>% formatRound("loading", digits = 3)
      })
    } else {
      output$pcaLoadings <- DT::renderDataTable({
        datatable(data.frame(Message = "Loadings not available for this method"), rownames = FALSE)
      })
    }

    output$dlTable_ordinationLoadings <- downloadHandler(
      filename = function() {
        generate_filename(
          base_name = paste0(ord_res$method, "_loadings"),
          design_name = designReactive(),
          extension = "xlsx"
        )
      },
      content = function(file) {
        req(figuresDataReactive$pc_loadings)
        export_excel_standard(
          df = as.data.frame(figuresDataReactive$pc_loadings),
          filename = file,
          sheet_name = "Loadings",
          table_name = "Loadings"
        )
      }
    )

    # Generate the main ordination plot
    if (!is.null(input$ordinationPlotType) && input$ordinationPlotType == "biplot") {
      # Require biplot-specific inputs
      req(
        !is.null(input$biplotGeneSelectionMethod),
        !is.null(input$biplotArrowScale),
        !is.null(input$biplotArrowColour)
      )

      # Determine top_n based on selection method
      biplot_top_n <- if (input$biplotGeneSelectionMethod == "per_pc") {
        if (!is.null(input$biplotNGenesPerPC)) input$biplotNGenesPerPC else 5
      } else {
        if (!is.null(input$biplotNGenes)) input$biplotNGenes else 10
      }

      # Prepare biplot data (use validated axis selections)
      biplot_data <- prepare_biplot_data(
        ordination_result = ord_res,
        pc_x = x_selected,
        pc_y = y_selected,
        top_n = biplot_top_n,
        custom_genes = biplot_custom_genes,
        selection_method = input$biplotGeneSelectionMethod
      )

      if (!is.null(biplot_data)) {
        biplot_data$gene_data$label <- label_features(biplot_data$gene_data$gene, input$fb_label_by_ordination)
      }

      # If biplot_data is NULL (no loadings available), fall back to samples-only plot
      if (is.null(biplot_data)) {
        # Create temporary input with validated axis selections
        input_validated <- reactiveValuesToList(input)
        input_validated$pcaX <- x_selected
        input_validated$pcaY <- y_selected
        pcaPlot <- plot_ordination_samples(
          ordination_result = ord_res,
          pc_scores = pc_scores,
          input = input_validated,
          coloursPCA = coloursPCA,
          levelLabels = pcaLevelLabels
        )
      } else {
        # Determine arrow colors based on selection method
        arrow_colour <- if (input$biplotGeneSelectionMethod == "per_pc") {
          # For per_pc mode, pass both colors
          list(x = input$biplotArrowColourX, y = input$biplotArrowColourY)
        } else {
          # For other modes, use single color
          input$biplotArrowColour
        }

        # Create temporary input with validated axis selections
        input_validated <- reactiveValuesToList(input)
        input_validated$pcaX <- x_selected
        input_validated$pcaY <- y_selected
        pcaPlot <- plot_ordination_biplot(
          ordination_result = ord_res,
          pc_scores = pc_scores,
          input = input_validated,
          coloursPCA = coloursPCA,
          biplot_data = biplot_data,
          arrow_scale = input$biplotArrowScale,
          arrow_colour = arrow_colour,
          selection_method = input$biplotGeneSelectionMethod,
          levelLabels = pcaLevelLabels
        )
      }
    } else {
      # Samples only plot (use validated axis selections)
      input_validated <- reactiveValuesToList(input)
      input_validated$pcaX <- x_selected
      input_validated$pcaY <- y_selected
      pcaPlot <- plot_ordination_samples(
        ordination_result = ord_res,
        pc_scores = pc_scores,
        input = input_validated,
        coloursPCA = coloursPCA,
        levelLabels = pcaLevelLabels
      )
    }

    figuresDataReactive$pcaStatic <- pcaPlot

    # Paired Plot
    dim_cols <- grep("^(PC|PCo|NMDS|dim)\\d+", colnames(pc_scores))
    dim_cols <- head(dim_cols, 5) # max 5 dimensions
    req(length(dim_cols) >= 2)
    pairedPlot <- ggpairs(
      data = pc_scores,
      columns = dim_cols,
      aes(fill = .data[[input$pcaFactor1]], alpha = input$alphaPCA),
      upper = list(
        continuous = wrap(
          "points",
          alpha = input$alphaPCA,
          size = input$dotSizePCA,
          alpha = input$alphaPCA,
          pch = 21,
          stroke = input$dotBorderPCA
        )
      ),
      lower = list(
        continuous = wrap(
          "points",
          alpha = input$alphaPCA,
          size = input$dotSizePCA,
          alpha = input$alphaPCA,
          pch = 21,
          stroke = input$dotBorderPCA
        )
      ),
      diag = list(continuous = wrap("densityDiag", size = 0.3), discrete = "barDiag", na = "naDiag")
    ) +
      theme_bw(base_size = input$textSizePCA) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_colour_manual(values = as.character(coloursPCA), breaks = names(coloursPCA)) +
      scale_fill_manual(values = as.character(coloursPCA), breaks = names(coloursPCA))

    figuresDataReactive$pcaPaired <- pairedPlot

    output$dlTable_ordinationCoords <- downloadHandler(
      filename = function() {
        generate_filename(
          base_name = paste0(ord_res$method, "_coords"),
          design_name = designReactive(),
          extension = "xlsx"
        )
      },
      content = function(file) {
        export_excel_standard(
          df = figuresDataReactive$pc_scores,
          filename = file,
          sheet_name = "Coordinates",
          table_name = "Coordinates"
        )
      }
    )

    output$pcaNAWarning <- NULL
  }
)

output$pcaStatic <- renderPlot(
  {
    # Shows ordinationResults()' validation message, if any, in the plot area
    ordinationResults()
    req(figuresDataReactive$pcaStatic)
    figuresDataReactive$pcaStatic
  },
  width = function() {
    as.numeric(input$figWidthPCA)
  },
  height = function() {
    as.numeric(input$figHeightPCA)
  }
)

output$dlPlot_ordination <- downloadHandler(
  filename = function() {
    method_name <- if (!is.null(figuresDataReactive$ord_method)) figuresDataReactive$ord_method else "ordination"
    generate_filename(
      base_name = paste0(input$filename_ordination, "_", method_name),
      design_name = designReactive(),
      extension = tolower(input$downloadFormatPCA)
    )
  },
  content = function(file) {
    req(figuresDataReactive$pcaStatic)
    export_plot_standard(
      plot_object = figuresDataReactive$pcaStatic,
      filename = file,
      format = tolower(input$downloadFormatPCA),
      width = (input$figWidthPCA / 80),
      height = (input$figHeightPCA / 80),
      dpi = as.numeric(input$dpiPCA)
    )
  }
)

output$pcaPaired <- renderPlot(
  {
    req(figuresDataReactive$pcaPaired)
    figuresDataReactive$pcaPaired
  },
  width = function() {
    as.numeric(input$figWidthPCA)
  },
  height = function() {
    as.numeric(input$figHeightPCA)
  }
)
output$dlPlot_ordinationPairs <- downloadHandler(
  filename = function() {
    method_name <- if (!is.null(figuresDataReactive$ord_method)) figuresDataReactive$ord_method else "ordination"
    generate_filename(
      base_name = paste0(input$filename_ordination, "_", method_name, "_paired"),
      design_name = designReactive(),
      extension = tolower(input$downloadFormatPCA)
    )
  },
  content = function(file) {
    req(figuresDataReactive$pcaPaired)
    export_plot_standard(
      plot_object = figuresDataReactive$pcaPaired,
      filename = file,
      format = tolower(input$downloadFormatPCA),
      width = (input$figWidthPCA / 80),
      height = (input$figHeightPCA / 80),
      dpi = as.numeric(input$dpiPCA)
    )
  }
)

output$pcaBrushTable <- DT::renderDataTable({
  # Without scores the table would have no columns; DT's JavaScript fails on
  # that and the browser drops the rest of the update, including the plot's
  # validation message.
  pc_scores <- req(figuresDataReactive$pc_scores)
  xvar <- figuresDataReactive$x_selected
  yvar <- figuresDataReactive$y_selected
  brushed_df <- brushedPoints(pc_scores, input$pcaBrush, xvar = xvar, yvar = yvar)

  # Select columns to display
  cols_to_show <- c("Name", xvar, yvar)
  if (!is.null(input$pcaFactor1) && input$pcaFactor1 != "None") {
    cols_to_show <- c(cols_to_show, input$pcaFactor1)
  }
  cols_to_show <- cols_to_show[cols_to_show %in% colnames(brushed_df)]

  DT::datatable(
    data = brushed_df[, cols_to_show, drop = FALSE],
    filter = "top",
    caption = "Click and drag over dots on the PCA plot to see those samples in this table.",
    rownames = FALSE
  ) %>%
    DT::formatRound(columns = c(xvar, yvar), digits = 3)
})
