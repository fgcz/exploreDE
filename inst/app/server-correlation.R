# ---------------------------------------------------------------------------
# Phase 2 migration — complete. All reads go through inputDataReactive()$.eds.
# ---------------------------------------------------------------------------

# Create reactive for design to use across all outputs and download handlers
designReactive <- reactive({
  req(inputDataReactive())
  eds <- inputDataReactive()$.eds
  if (eds$has_de()) {
    contrast <- selectedContrast()
    eds$contrast_label(contrast)
  } else if (eds$omics_class() == "proteomics") {
    "No contrast selected"
  } else {
    eds$data_type() # QC fallback — no DE, use the data type name
  }
})

output$correlationStatic <- renderPlot({
  NULL
})

# Update some inputs
local({
  eds <- inputDataReactive()$.eds
  updateSelectizeInput(
    session = session,
    inputId = "correlationGene1",
    label = "Feature X",
    choices = eds$feature_ids(),
    selected = "",
    server = TRUE
  )
  updateSelectizeInput(
    session = session,
    inputId = "correlationGene2",
    label = "Feature Y",
    choices = eds$feature_ids(),
    selected = "",
    server = TRUE
  )
})

# Update feature choices when collapse guides toggle changes
observeEvent(
  input$correlationCollapseGuides,
  {
    if (isTRUE(restoring())) {
      return()
    } # save-state restore guard
    eds <- inputDataReactive()$.eds
    if (eds$data_type() == "CrisprCountQC") {
      if (isTRUE(input$correlationCollapseGuides)) {
        gene_choices <- sort(unique(eds$sgrna_to_gene()))
      } else {
        gene_choices <- eds$feature_ids()
      }
      updateSelectizeInput(session, "correlationGene1", choices = gene_choices, selected = "", server = TRUE)
      updateSelectizeInput(session, "correlationGene2", choices = gene_choices, selected = "", server = TRUE)
      correlationPlot$plot <- NULL
    }
  },
  ignoreInit = TRUE
)
local({
  eds <- inputDataReactive()$.eds
  updateSelectInput(
    session = session,
    inputId = "correlationColourBy",
    label = "Select column to colour by",
    choices = eds$factor_names(),
    selected = eds$primary_grouping()
  )
  updateSelectInput(
    session = session,
    inputId = "correlationShapeBy",
    label = "Select column to shape by",
    choices = c("None", eds$factor_names()),
    selected = "None"
  )
})
local({
  eds <- inputDataReactive()$.eds
  assay_keys <- eds$assay_names()
  count_choices <- setNames(assay_keys, vapply(assay_keys, eds$assay_label, character(1)))
  default_count <- if (eds$data_type() == "RNASeq" && "normalised_log2" %in% assay_keys) {
    "normalised_log2"
  } else if (eds$omics_class() == "proteomics" && "intensities" %in% assay_keys) {
    "intensities"
  } else {
    assay_keys[1]
  }
  updateSelectInput(
    session = session,
    inputId = "correlationCounts",
    label = "Select count method to plot",
    choices = count_choices,
    selected = default_count
  )
})
# Batch + design dropdowns — see server-boxplot.R for the decision table.
observeEvent(inputDataReactive(), {
  eds <- inputDataReactive()$.eds
  .dt <- eds$data_type()
  .oc <- eds$omics_class()
  .factors <- eds$factor_names()
  .grouping <- eds$primary_grouping()
  .secondary <- eds$secondary_grouping()
  .twoFactor <- !is.na(.secondary) && nzchar(.secondary) && .secondary %in% .factors
  # Proteomics grouping: SUSHI "G_" if present, else the primary grouping (h5ad).
  .protGroup <- if ("G_" %in% .factors) "G_" else .grouping
  .deDesignCol <- if (.dt == "RNASeq") {
    .grouping
  } else if (.oc == "proteomics") {
    .protGroup
  } else {
    NA_character_
  }

  batchChoices <- c("None", .factors)
  batchDefault <- if (.twoFactor) .secondary else "None"
  updateSelectInput(
    session = session,
    inputId = "correlationBatch",
    label = "Apply a batch-correction to the counts?",
    choices = batchChoices,
    selected = batchDefault
  )

  designChoices <- if ((.dt == "RNASeq" || .oc == "proteomics") && .twoFactor) {
    c("None", .deDesignCol)
  } else {
    c("None", .factors)
  }
  designDefault <- if (.dt == "RNASeq") {
    .grouping
  } else if (.oc == "proteomics") {
    .protGroup
  } else {
    "None"
  }
  updateSelectInput(
    session = session,
    inputId = "correlationDesign",
    label = "Design factor to preserve during batch correction",
    choices = designChoices,
    selected = designDefault
  )
})
updateSliderInput(
  session = session,
  inputId = "correlationCorPosX",
  label = "Label position x",
  min = 1,
  max = 100,
  value = 10,
  step = 1
)
updateSliderInput(
  session = session,
  inputId = "correlationCorPosY",
  label = "Label position y",
  min = 1,
  max = 100,
  value = 10,
  step = 1
)

observeEvent(
  {
    input$correlationColourBy
  },
  {
    output$correlationGroupBucket <- renderUI({
      req(input$tabs == "correlationsTab")
      eds <- inputDataReactive()$.eds
      .lvls <- unique(levels(as.factor(eds$samples()[[input$correlationColourBy]])))
      # save-state restore: seed include/exclude split + order from a saved session.
      .keep <- restore_default("correlationKeepBucket", .lvls, valid = .lvls)
      .excl <- restore_default("correlationExcludeBucket", character(0), valid = .lvls)
      .keep <- setdiff(.keep, .excl)
      bucket_list(
        header = "Drag and drop groups in order to be plotted (Ctrl/Cmd to select multiple)",
        group_name = "correlation_bucket_list",
        orientation = "horizontal",
        add_rank_list(
          text = "Include these groups",
          labels = as.list(.keep),
          input_id = "correlationKeepBucket",
          options = sortable::sortable_options(multiDrag = TRUE)
        ),
        add_rank_list(
          text = "Exclude these groups",
          labels = as.list(.excl),
          input_id = "correlationExcludeBucket",
          options = sortable::sortable_options(multiDrag = TRUE)
        )
      )
    })
    outputOptions(output, "correlationGroupBucket", suspendWhenHidden = FALSE)
  }
)

correlationPlot <- reactiveValues(plot = NULL)

# Update feature choices when search mode changes
observeEvent(input$fb_search_by_correlation, ignoreInit = TRUE, {
  if (isTRUE(restoring())) {
    return()
  } # save-state restore guard
  eds <- inputDataReactive()$.eds
  mode <- input$fb_search_by_correlation
  pc_gene <- eds$proteomics_choices("gene")$choices
  pc_protein <- eds$proteomics_choices("protein")$choices
  choices <- if (mode == "gene_name" && !is.null(pc_gene)) {
    pc_gene
  } else if (mode == "protein_Id" && !is.null(pc_protein)) {
    pc_protein
  } else {
    eds$feature_ids()
  }
  updateSelectizeInput(session, "correlationGene1", choices = choices, selected = "", server = TRUE)
  updateSelectizeInput(session, "correlationGene2", choices = choices, selected = "", server = TRUE)
  correlationPlot$plot <- NULL
})

observeEvent(
  {
    input$correlationMethod
    input$correlationGene1
    input$correlationGene2
    input$correlationColourBy
    input$correlationShapeBy
    input$correlationCounts
    input$correlationBatch
    input$correlationDesign
    input$correlationCorPosX
    input$correlationCorPosY
    input$correlationPointSize
    input$correlationPointStroke
    input$correlationPointAlpha
    input$textSizeCorrelation
    input$correlationKeepBucket
    input$correlationCountsLog
    input$correlationShowNames
    input$correlationMaxOverlaps
    input$correlationLabelPosX
    input$correlationLabelPosY
    # Add all color picker + level label inputs as dependencies so plot updates when they change
    local({
      eds_dep <- inputDataReactive()$.eds
      lapply(eds_dep$factor_names(), function(fname) {
        lapply(eds_dep$factor_levels(fname), function(lev) {
          list(input[[paste0("GroupColour", fname, lev)]], input[[paste0("GroupLabel", fname, lev)]])
        })
      })
    })
  },
  ignoreNULL = TRUE,
  ignoreInit = TRUE,
  {
    req(input$correlationGene1 != input$correlationGene2)
    req(nchar(input$correlationGene1) > 1 & nchar(input$correlationGene2) > 1)

    eds <- inputDataReactive()$.eds

    # prepare the metadata
    datasetCorrelation <- eds$samples() %>% data.frame(check.names = F)
    if (!"Name" %in% colnames(datasetCorrelation)) {
      datasetCorrelation <- datasetCorrelation %>% rownames_to_column("Name")
    }
    datasetCorrelation$Name <- factor(datasetCorrelation$Name, levels = datasetCorrelation$Name)
    datasetCorrelation <- datasetCorrelation[, c(
      "Name",
      grepv(paste0("^", eds$factor_names(), "$", collapse = "|"), colnames(datasetCorrelation))
    )]

    # get the counts (full matrix — all features, all samples)
    countsCorrelation <- eds$assay(input$correlationCounts)
    # Collapse sgRNA guides to gene level if requested
    if (eds$data_type() == "CrisprCountQC" && isTRUE(input$correlationCollapseGuides)) {
      countsCorrelation <- collapse_guides_to_genes(countsCorrelation, eds$sgrna_to_gene())
    }

    # Log2 and batch correction on full matrix BEFORE extracting genes
    is_logged <- input$correlationCounts %in% c("normalised_log2", "vst") || eds$omics_class() == "proteomics"
    if (input$correlationCountsLog & !is_logged) {
      countsCorrelation <- log2(countsCorrelation + 1)
      is_logged <- TRUE
    }
    if (input$correlationBatch != "None") {
      if (is_logged) {
        countsCorrelation <- apply_batch_correction(
          counts = countsCorrelation,
          dataset = eds$samples(),
          batchCol = input$correlationBatch,
          designCol = input$correlationDesign
        )
      } else {
        showNotification(
          "Batch correction requires log-scale data. Please select a logged count type or enable the log2 transform.",
          type = "warning",
          duration = 8
        )
      }
    }

    # Now extract the two genes of interest
    corr_mode <- input$fb_search_by_correlation %||% "feature_id"
    corr_features <- list(input$correlationGene1, input$correlationGene2)
    corr_labels <- c(input$correlationGene1, input$correlationGene2)
    if (corr_mode != "feature_id") {
      corr_features <- lapply(corr_labels, function(sel) {
        resolve_features(sel, corr_mode)
      })
    }
    # Average across feature_ids when a gene/protein maps to multiple rows
    corr_axis_labels <- corr_labels
    corr_rows <- lapply(seq_along(corr_features), function(i) {
      fids <- corr_features[[i]]
      fids <- fids[fids %in% rownames(countsCorrelation)]
      req(length(fids) > 0)
      if (length(fids) == 1) {
        countsCorrelation[fids, , drop = FALSE]
      } else {
        corr_axis_labels[i] <<- paste0(corr_labels[i], " (mean of ", length(fids), " features)")
        showNotification(
          paste0(corr_labels[i], ": averaging ", length(fids), " features for correlation"),
          type = "warning",
          duration = 6
        )
        matrix(
          colMeans(countsCorrelation[fids, , drop = FALSE], na.rm = TRUE),
          nrow = 1,
          dimnames = list(NULL, colnames(countsCorrelation))
        )
      }
    })
    countsCorrelation <- do.call(rbind, corr_rows)
    rownames(countsCorrelation) <- corr_labels

    # Warn about sparsity
    n_samples <- ncol(countsCorrelation)
    for (i in seq_along(corr_labels)) {
      n_valid <- sum(!is.na(countsCorrelation[i, ]) & !is.nan(as.numeric(countsCorrelation[i, ])))
      if (n_valid < n_samples) {
        showNotification(
          paste0(corr_labels[i], ": data available for ", n_valid, " of ", n_samples, " samples"),
          type = if (n_valid < n_samples / 2) "error" else "warning",
          duration = 8
        )
      }
    }

    # prepare the counts
    countsCorrelation <- countsCorrelation %>%
      t() %>%
      data.frame(check.names = F) %>%
      rownames_to_column("Name") %>%
      left_join(datasetCorrelation, by = "Name")
    countsCorrelation <- countsCorrelation[
      which(countsCorrelation[[input$correlationColourBy]] %in% input$correlationKeepBucket),
    ]

    # get the colours
    coloursCorrelation <- setNames(
      lapply(input$correlationKeepBucket, function(k) {
        paste0(col2hex(input[[paste0("GroupColour", input$correlationColourBy, k)]]), "FF")
      }),
      input$correlationKeepBucket
    )

    # Display labels from the Colours tab (fill factor, and shape factor).
    # If either is rich, both are, as they share the markdown legend text.
    corrShapeFactor <- if (isTRUE(input$correlationShapeBy != "None")) input$correlationShapeBy else NULL
    corrShapeLevels <- if (is.null(corrShapeFactor)) NULL else levels(as.factor(countsCorrelation[[corrShapeFactor]]))
    corrLabels <- level_labels(input, input$correlationColourBy, input$correlationKeepBucket)
    corrShapeLabels <- if (is.null(corrShapeFactor)) NULL else level_labels(input, corrShapeFactor, corrShapeLevels)
    corrRich <- is_rich_labels(corrLabels) || is_rich_labels(corrShapeLabels)
    if (corrRich) {
      corrLabels <- level_labels(input, input$correlationColourBy, input$correlationKeepBucket, force_rich = TRUE)
      if (!is.null(corrShapeFactor)) {
        corrShapeLabels <- level_labels(input, corrShapeFactor, corrShapeLevels, force_rich = TRUE)
      }
    }

    # get the correct coefficient
    if (input$correlationMethod == "Spearman") {
      cor <- "rho"
    } else if (input$correlationMethod == "Kendall") {
      cor <- "tau"
    } else {
      cor <- "R"
    }

    # make the plot
    corr_plot1 <- ggplot(
      countsCorrelation,
      aes(x = .data[[input$correlationGene1]], y = .data[[input$correlationGene2]])
    )
    corr_plot1 <- corr_plot1 + geom_smooth(method = "lm", se = TRUE, colour = "black")
    corr_plot1 <- corr_plot1 +
      stat_cor(
        method = tolower(input$correlationMethod),
        cor.coef.name = cor,
        size = input$textSizeCorrelation / 3,
        show.legend = FALSE,
        label.x.npc = tolower(input$correlationCorPosX),
        label.y.npc = tolower(input$correlationCorPosY),
        label.args = list(fontface = "bold")
      )
    if (input$correlationShapeBy == "None") {
      corr_plot1 <- corr_plot1 +
        geom_point(
          aes(fill = .data[[input$correlationColourBy]]),
          pch = 21,
          size = input$correlationPointSize,
          stroke = input$correlationPointStroke,
          alpha = input$correlationPointAlpha
        )
    } else {
      corr_plot1 <- corr_plot1 +
        geom_point(
          aes(fill = .data[[input$correlationColourBy]], shape = .data[[input$correlationShapeBy]]),
          size = input$correlationPointSize,
          stroke = input$correlationPointStroke,
          alpha = input$correlationPointAlpha
        )
      corr_plot1 <- corr_plot1 +
        scale_shape_manual(
          values = c(rep(c(21, 22, 23, 24, 25), times = 10))[
            1:nlevels(as.factor(countsCorrelation[[input$correlationShapeBy]]))
          ],
          labels = corrShapeLabels
        )
      corr_plot1 <- corr_plot1 +
        guides(
          fill = guide_legend(override.aes = list(shape = 21, size = 3)),
          shape = guide_legend(override.aes = list(fill = "black", size = 3))
        )
    }
    if (input$correlationShowNames) {
      corr_plot1 <- corr_plot1 +
        geom_label_repel(
          aes(label = Name, fill = .data[[input$correlationColourBy]]),
          show.legend = FALSE,
          alpha = input$correlationPointAlpha,
          size = input$textSizeCorrelation / 3,
          nudge_x = input$correlationLabelPosX / 100,
          nudge_y = input$correlationLabelPosY / 100,
          segment.size = 0.5,
          segment.color = "black",
          max.overlaps = input$correlationMaxOverlaps
        )
    }
    corr_plot1 <- corr_plot1 +
      scale_fill_manual(
        breaks = names(coloursCorrelation),
        values = as.character(coloursCorrelation),
        labels = corrLabels
      )
    corr_plot1 <- corr_plot1 + theme_prism(base_size = input$textSizeCorrelation)
    corr_plot1 <- corr_plot1 +
      theme(legend.text = element_text(size = input$textSizeCorrelation, color = "black", face = "bold"))
    if (corrRich) {
      corr_plot1 <- corr_plot1 + rich_label_theme(legend = TRUE)
    }
    corr_plot1 <- corr_plot1 + labs(x = corr_axis_labels[1], y = corr_axis_labels[2])
    corr_plot1
    correlationPlot$plot <- corr_plot1
    correlationPlot$counts <- countsCorrelation
  }
)

# view the plot
output$correlationStatic <- renderPlot(
  {
    req(!is.null(correlationPlot$plot))
    print(correlationPlot$plot)
  },
  width = function() {
    as.numeric(input$figWidthCorrelation)
  },
  height = function() {
    as.numeric(input$figHeightCorrelation)
  }
)
# download the plot
# download the plot
output$dlPlot_correlation <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = input$filename_correlation,
      design_name = designReactive(),
      extension = tolower(input$downloadFormatCorrelation)
    )
  },
  content = function(file) {
    req(!is.null(correlationPlot$plot))
    export_plot_standard(
      plot_object = correlationPlot$plot,
      filename = file,
      format = tolower(input$downloadFormatCorrelation),
      width = as.numeric(input$figWidthCorrelation) / 75,
      height = as.numeric(input$figHeightCorrelation) / 75,
      dpi = as.numeric(input$dpiCorrelation)
    )
  }
)
# download the counts
# download the counts
output$dlTable_correlation <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_correlation, "_Counts"),
      design_name = designReactive(),
      extension = "xlsx"
    )
  },
  content = function(file) {
    req(!is.null(correlationPlot$counts))
    export_excel_standard(
      df = correlationPlot$counts,
      filename = file,
      sheet_name = "Correlation Data",
      table_name = "CorrelationData"
    )
  }
)
