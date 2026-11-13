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

eds <- inputDataReactive()$.eds
req(eds)

# Sample ratios (Settings > Sample ratios):
# name -> list(name, numerator, denominator, values, assay, batch)
heatmapRatios <- reactiveValues(ratios = list())
HEATMAP_RATIO_MAX <- 20L
# The contrasts, offered on the Ratios heatmap next to the sample ratios as
# their log2 fold changes (values: contrast names; names: their labels)
contrastRatioChoices <- stats::setNames(eds$contrasts(), vapply(eds$contrasts(), eds$contrast_label, character(1)))

observe({
  # Detect species for KEGG heatmap
  keggHeatmapReactive$keggCode <- if (kegg_installed()) kegg_org_code(eds$ref_build())

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
              numericInput(
                inputId = "heatmapKeggLFC",
                label = "Log2 Fold Change Threshold:",
                value = 0,
                min = 0,
                max = 10,
                step = 0.25
              ),
              selectInput(
                inputId = "heatmapKeggPType",
                choices = c("FDR", "Raw"),
                label = "P-Value:",
                selected = "FDR"
              ),
              selectInput(
                inputId = "heatmapKeggPThreshold",
                choices = c(1.0, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001),
                label = "P-Value Threshold:",
                selected = 0.05
              )
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

  if (eds$omics_class() == "proteomics") {
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
  assay_keys <- eds$assay_names()
  count_choices <- setNames(assay_keys, vapply(assay_keys, eds$assay_label, character(1)))
  default_count <- if (eds$omics_class() == "proteomics") {
    if ("intensities" %in% assay_keys) "intensities" else assay_keys[1]
  } else {
    if ("normalised_log2" %in% assay_keys) "normalised_log2" else assay_keys[1]
  }
  updateSelectInput(
    session = session,
    inputId = "heatmapCounts",
    choices = count_choices,
    selected = default_count
  )

  updateSelectizeInput(
    session = session,
    inputId = "fb_select_heatmap",
    choices = eds$feature_ids(),
    selected = NULL,
    server = T,
  )

  for (id in c("heatmapRatioNumerator", "heatmapRatioDenominator")) {
    updateSelectizeInput(session = session, inputId = id, choices = eds$sample_ids(), selected = NULL, server = TRUE)
  }

  # Populate exclude dropdown for CRISPR data
  if (eds$data_type() == "CrisprCountQC") {
    updateSelectizeInput(session, "heatmapExcludeGuides", choices = eds$feature_ids(), selected = NULL, server = TRUE)
  }

  # Update feature choices when collapse guides toggle changes
  observeEvent(
    input$heatmapCollapseGuides,
    {
      if (isTRUE(restoring())) {
        return()
      } # save-state restore guard
      if (eds$data_type() == "CrisprCountQC" && identical(input$heatmapCollapseGuides, "gene")) {
        gene_choices <- sort(unique(eds$sgrna_to_gene()))
      } else {
        gene_choices <- eds$feature_ids()
      }
      updateSelectizeInput(session, "fb_select_heatmap", choices = gene_choices, selected = "", server = TRUE)
      updateTextAreaInput(session, "fb_text_heatmap", value = "")
      updateTextInput(session, "fb_pattern_heatmap", value = "")
      updateSelectizeInput(session, "heatmapExcludeGuides", choices = gene_choices, selected = NULL, server = TRUE)
    },
    ignoreInit = TRUE
  )

  # Batch + design dropdowns: batch_design_defaults() in R/utils_explore_dataset.R.
  bd <- batch_design_defaults(eds)
  updateSelectInput(
    session = session,
    inputId = "heatmapBatch",
    choices = bd$batch_choices,
    selected = bd$batch_default
  )
  updateSelectInput(
    session = session,
    inputId = "heatmapDesign",
    choices = bd$design_choices,
    selected = bd$design_default
  )

  # Update p-value type based on data type. FDR-only producers (no raw p-value,
  # e.g. SAINTexpress) are locked to FDR so the significance filter is not an
  # all-NA Raw axis; otherwise proteomics defaults to Raw, genomics to FDR.
  if (!isTRUE(eds$has_pvalue())) {
    updateSelectInput(session = session, inputId = "pTypeHeatmap", choices = c("FDR"), selected = "FDR")
  } else if (eds$omics_class() == "proteomics") {
    updateSelectInput(session = session, inputId = "pTypeHeatmap", selected = "Raw")
  } else {
    updateSelectInput(session = session, inputId = "pTypeHeatmap", selected = "FDR")
  }

  # Proteomics-specific: peptide filter UI
  if (eds$omics_class() == "proteomics") {
    if (
      length(eds$contrasts()) > 0 &&
        "nrPeptides" %in% colnames(eds$de(eds$contrasts()[[1]]))
    ) {
      output$nrPeptidesHeatmapUI <- renderUI({
        sliderInput(
          inputId = "nrPeptideHeatmap",
          label = "Minimum number of peptides",
          min = 0,
          max = 20,
          value = 2,
          step = 1,
          width = "85%"
        )
      })
      outputOptions(output, "nrPeptidesHeatmapUI", suspendWhenHidden = FALSE)
    } else {
      output$nrPeptidesHeatmapUI <- renderUI({
        NULL
      })
    }
  }

  # GO term selection. Like the KEGG sub-tab, the GO sub-tab only exists when
  # it can work: RNA-seq with DE (the fetch reads the DE table) and GO
  # annotation columns in the features.
  pathways <- if (eds$data_type() == "RNASeq" && eds$has_de()) eds$all_pathways() else character(0)
  if (!length(pathways)) {
    hideTab(inputId = "sig", target = "GO")
  } else {
    updateSelectizeInput(
      session,
      inputId = "heatmapGoInput",
      choices = pathways,
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
  }
})

output$heatmapColumnSelectUI <- renderUI({
  factors <- eds$factor_names()
  grouping <- eds$primary_grouping()

  # Calculate default selection
  defaultSelected <- if (eds$data_type() == "RNASeq" && !is.null(grouping) && !is.na(grouping)) {
    grouping
  } else if (length(factors) > 0) {
    factors[[1]]
  } else {
    NULL
  }

  # Get current selection and available choices
  choices <- factors
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
    if (
      !is.null(input$heatmapColumnSelect) &&
        length(input$heatmapColumnSelect) > 0 &&
        !all(input$heatmapColumnSelect == "")
    ) {
      lapply(seq_along(input$heatmapColumnSelect), function(i) {
        output[[paste0("heatmapBucket", i)]] <- renderUI({
          .lvls <- unique(levels(as.factor(eds$samples()[[input$heatmapColumnSelect[[i]]]])))
          # save-state restore: seed include/exclude split + order per factor.
          .keep <- restore_default(paste0("hmKeepBucket", i), .lvls, valid = .lvls)
          .excl <- restore_default(paste0("hmExcludeBucket", i), character(0), valid = .lvls)
          .keep <- setdiff(.keep, .excl)
          bucket_list(
            header = paste(input$heatmapColumnSelect[[i]], "(Drag/Drop levels - Ctrl/Cmd to select multiple)"),
            group_name = paste0("hmBucket", i),
            orientation = "horizontal",
            add_rank_list(
              text = "Include these groups",
              labels = as.list(.keep),
              input_id = paste0("hmKeepBucket", i),
              options = sortable::sortable_options(multiDrag = TRUE)
            ),
            add_rank_list(
              text = "Exclude these groups",
              labels = as.list(.excl),
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
  # save-state restore: restore_default yields the saved break during a restore
  # (this fires on Z-score/Log2 change, which a restore triggers).
  if (input$heatmapZScore %in% c("Z-Score", "Centred")) {
    updateNumericInput(
      session = session,
      inputId = "heatmapAtLow",
      value = restore_default("heatmapAtLow", -2),
      min = -10,
      max = 10,
      step = 0.1
    )
    updateNumericInput(
      session = session,
      inputId = "heatmapAtMid",
      value = restore_default("heatmapAtMid", 0),
      min = -10,
      max = 10,
      step = 0.1
    )
    updateNumericInput(
      session = session,
      inputId = "heatmapAtHigh",
      value = restore_default("heatmapAtHigh", 2),
      min = -10,
      max = 10,
      step = 0.1
    )
  } else if (input$heatmapLog2 & input$heatmapZScore == "None") {
    updateNumericInput(
      session = session,
      inputId = "heatmapAtLow",
      value = restore_default("heatmapAtLow", 0),
      min = -10,
      max = 10,
      step = 0.1
    )
    updateNumericInput(
      session = session,
      inputId = "heatmapAtMid",
      value = restore_default("heatmapAtMid", 6),
      min = -10,
      max = 10,
      step = 0.1
    )
    updateNumericInput(
      session = session,
      inputId = "heatmapAtHigh",
      value = restore_default("heatmapAtHigh", 12),
      min = -10,
      max = 10,
      step = 0.1
    )
  } else if (!input$heatmapLog2 & input$heatmapZScore == "None") {
    updateNumericInput(
      session = session,
      inputId = "heatmapAtLow",
      value = restore_default("heatmapAtLow", 0),
      min = 0,
      max = 1e6,
      step = 100
    )
    updateNumericInput(
      session = session,
      inputId = "heatmapAtMid",
      value = restore_default("heatmapAtMid", 1000),
      min = 0,
      max = 1e6,
      step = 100
    )
    updateNumericInput(
      session = session,
      inputId = "heatmapAtHigh",
      value = restore_default("heatmapAtHigh", 10000),
      min = 0,
      max = 1e6,
      step = 100
    )
  }
})

# Update KEGG expression threshold step based on count scale
observeEvent(
  {
    input$heatmapCounts
    input$heatmapLog2
    input$heatmapZScore
  },
  ignoreInit = TRUE,
  {
    isLogSpace <- input$heatmapCounts %in%
      c("normalised_log2", "vst") ||
      (input$heatmapLog2 && input$heatmapCounts %in% c("tpm", "fpkm", "normalised", "counts"))
    isScaled <- input$heatmapZScore %in% c("Z-Score", "Centred")

    stepVal <- if (isScaled || isLogSpace) 0.1 else 1
    updateNumericInput(session = session, inputId = "heatmapKeggExpr", step = stepVal)
    updateNumericInput(session = session, inputId = "heatmapGoExpr", step = stepVal)
  }
)

# Update settings for CountQC data type
observe({
  if (eds$data_type() == "CountQC") {
    # For CountQC, hide gene names by default since there will be many features
    updateCheckboxInput(session = session, inputId = "geneNamesHeatmap", value = FALSE)
    # Set cell border size to 0 for CountQC (will have many features)
    updateSliderInput(session = session, inputId = "cellBorderSizeHeatmap2", value = 0)
  } else if (eds$data_type() == "CrisprCountQC") {
    # CRISPR data has fewer features, so show names by default
    updateCheckboxInput(session = session, inputId = "geneNamesHeatmap", value = TRUE)
  }
})

# Dynamic tab title based on data type
output$heatmapTabTitle <- renderUI({
  if (eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
    "SD Features"
  } else {
    "DE Features"
  }
})

ht_opt$message <- FALSE

# Colour scale of the heatmaps on this tab, from the Figure settings colour
# and palette inputs; `at` overrides the three breaks (the Ratios heatmap
# passes symmetric ones, the count heatmaps use the Figure settings breaks).
.heatmap_colours <- function(at = NULL) {
  if (is.null(at)) {
    at <- c(input$heatmapAtLow, input$heatmapAtMid, input$heatmapAtHigh)
  }
  if (input$heatmapLimitCPalette == "None") {
    return(colorRamp2(
      breaks = rev(at),
      colors = c(input$heatmapColourRed, input$heatmapColourWhite, input$heatmapColourBlue)
    ))
  }
  palette_colours <- brewer.pal(9, input$heatmapLimitCPalette)
  selected_colours <- palette_colours[c(1, 5, 9)]
  if (input$heatmapLimitCPaletteRev) {
    colorRamp2(breaks = at, colors = selected_colours)
  } else {
    colorRamp2(at, rev(selected_colours))
  }
}

# Count matrices ----
# The count method of the Feature settings tab, all features and samples, with
# CRISPR guides collapsed to genes or left out as set there. Cached: a redraw
# does not compute the assay (VST, FPKM, ...) again. `collapsed`: rows are
# genes, not sgRNAs.
heatmapAssay <- reactive({
  req(input$heatmapCounts)
  counts <- eds$assay(input$heatmapCounts)
  collapsed <- eds$data_type() == "CrisprCountQC" && identical(input$heatmapCollapseGuides, "gene")
  if (eds$data_type() == "CrisprCountQC") {
    if (collapsed) {
      counts <- collapse_guides_to_genes(counts, eds$sgrna_to_gene(), exclude = input$heatmapExcludeGuides)
    } else if (length(input$heatmapExcludeGuides) > 0) {
      counts <- counts[!rownames(counts) %in% input$heatmapExcludeGuides, , drop = FALSE]
    }
  }
  list(counts = counts, collapsed = collapsed)
})

# heatmapAssay() with the log2 and batch settings applied; the heatmaps draw
# from it. Cached: colour, label and layout changes do not prepare it again.
heatmapPrepared <- reactive({
  prepare_expression(
    heatmapAssay()$counts,
    eds$samples(),
    is_log2 = eds$assay_is_log2(input$heatmapCounts),
    take_log2 = input$heatmapLog2,
    batch = input$heatmapBatch,
    design = input$heatmapDesign
  )
})

# Sample ratios ----
# The matrix a ratio is taken from: heatmapAssay() with the batch correction
# of the Feature settings tab, always on the log2 scale. Returns the matrix,
# whether guides were collapsed to genes, and the batch column applied ("None"
# when none was, or when the correction could not be applied and returned the
# matrix as it was).
heatmap_ratio_source <- function() {
  src <- heatmapAssay()
  prep <- prepare_expression(
    src$counts,
    eds$samples(),
    is_log2 = eds$assay_is_log2(input$heatmapCounts),
    take_log2 = TRUE,
    batch = input$heatmapBatch %||% "None",
    design = input$heatmapDesign %||% "None"
  )
  list(counts = prep$counts, collapsed = src$collapsed, batch = prep$batch)
}

# The ratio the add observer wants selected on the Ratios heatmap, picked up
# by the choices observer below when it next runs.
ratioSelectionPending <- reactiveVal(NULL)

observeEvent(input$heatmapRatioAdd, {
  name <- trimws(input$heatmapRatioName %||% "")
  num <- input$heatmapRatioNumerator
  den <- input$heatmapRatioDenominator
  if (!nzchar(name)) {
    showNotification("Give the ratio a name.", type = "warning")
    return()
  }
  if (name %in% names(heatmapRatios$ratios)) {
    showNotification(paste0("A ratio named '", name, "' already exists."), type = "warning")
    return()
  }
  if (name %in% c(contrastRatioChoices, names(contrastRatioChoices))) {
    showNotification(paste0("'", name, "' is a contrast of this dataset; choose another name."), type = "warning")
    return()
  }
  if (identical(name, "feature_id")) {
    showNotification("'feature_id' names the feature column of the download; choose another name.", type = "warning")
    return()
  }
  if (!length(num) || !length(den)) {
    showNotification("Choose at least one numerator and one denominator sample.", type = "warning")
    return()
  }
  if (length(intersect(num, den))) {
    showNotification("A sample cannot be in both the numerator and the denominator.", type = "warning")
    return()
  }
  if (length(heatmapRatios$ratios) >= HEATMAP_RATIO_MAX) {
    showNotification(
      paste("At most", HEATMAP_RATIO_MAX, "ratios can be defined; remove one first."),
      type = "warning"
    )
    return()
  }
  req(input$heatmapCounts %in% eds$assay_names())
  tryCatch(
    {
      src <- heatmap_ratio_source()
      heatmapRatios$ratios[[name]] <- list(
        name = name,
        numerator = num,
        denominator = den,
        values = sample_log2_ratio(src$counts, num, den),
        assay = eds$assay_label(input$heatmapCounts),
        collapsed = src$collapsed,
        batch = src$batch
      )
      ratioSelectionPending(name)
      updateTextInput(session = session, inputId = "heatmapRatioName", value = "")
      showNotification(paste0("Ratio '", name, "' added."), type = "message", duration = 4)
    },
    error = function(e) {
      showNotification(paste("Could not add the ratio:", conditionMessage(e)), type = "error", duration = 10)
    }
  )
})

observeEvent(input$heatmapRatioRemoveBtn, {
  drop <- input$heatmapRatioRemove
  req(length(drop) > 0)
  ratios <- heatmapRatios$ratios
  heatmapRatios$ratios <- ratios[setdiff(names(ratios), drop)]
})

output$heatmapRatioListUI <- renderUI({
  ratios <- heatmapRatios$ratios
  if (!length(ratios)) {
    return(helpText("No sample ratios defined yet."))
  }
  tagList(
    h5("Defined ratios"),
    lapply(ratios, function(r) {
      div(
        style = "margin-bottom: 6px;",
        tags$b(r$name),
        br(),
        tags$small(paste0(
          "mean(",
          paste(r$numerator, collapse = ", "),
          ") - mean(",
          paste(r$denominator, collapse = ", "),
          ") on ",
          r$assay,
          if (isTRUE(r$collapsed)) ", guides collapsed to genes",
          ", batch: ",
          r$batch
        ))
      )
    })
  )
})

# The ratios the Ratios heatmap can show: sample ratios first, then the
# contrasts' log2 fold changes. Keeps the current selection, plus the ratio
# the add observer just stored. Runs when the sample ratios change.
observe({
  custom <- names(heatmapRatios$ratios) %||% character(0)
  choices <- c(stats::setNames(custom, custom), contrastRatioChoices)
  current <- isolate(input$heatmapRatioPlot) %||% character(0)
  selected <- unique(c(intersect(current, choices), isolate(ratioSelectionPending())))
  ratioSelectionPending(NULL)
  updateSelectizeInput(
    session = session,
    inputId = "heatmapRatioPlot",
    choices = choices,
    selected = if (length(selected)) selected else NULL
  )
  updateSelectizeInput(session = session, inputId = "heatmapRatioRemove", choices = custom, selected = NULL)
})

# The Ratios heatmap does not depend on the contrast or the sample groups.
# Its matrix follows the ratio selection only (debounced); the heatmap object
# follows the matrix and the styling inputs, so a cosmetic change never
# recomputes the matrix. The matrix also feeds the download.
heatmapRatioSelection <- debounce(
  reactive({
    list(
      plot = input$heatmapRatioPlot,
      top_n = input$heatmapRatioTopN,
      rank_by = input$heatmapRatioRank,
      ratios = heatmapRatios$ratios
    )
  }),
  1000
)

heatmapRatioMatrix <- reactive({
  sel <- heatmapRatioSelection()
  custom <- sel$ratios
  ratioNames <- sel$plot %||% character(0)
  ratioNames <- ratioNames[ratioNames %in% c(names(custom), contrastRatioChoices)]
  if (!length(ratioNames)) {
    return(NULL)
  }
  ratioList <- lapply(ratioNames, function(nm) {
    if (nm %in% names(custom)) custom[[nm]]$values else contrast_log2_ratio(eds, nm)
  })
  names(ratioList) <- ifelse(
    ratioNames %in% names(custom),
    ratioNames,
    names(contrastRatioChoices)[match(ratioNames, contrastRatioChoices)]
  )
  dropped <- names(ratioList)[vapply(ratioList, is.null, logical(1))]
  if (length(dropped)) {
    showNotification(
      paste0(
        "No log2 fold change in the DE table of ",
        paste(dropped, collapse = ", "),
        ", so ",
        if (length(dropped) == 1) "it is" else "they are",
        " left out of the Ratios heatmap."
      ),
      type = "warning",
      duration = 8
    )
  }
  top_n <- sel$top_n
  top_n <- if (is.null(top_n) || is.na(top_n)) 50 else max(2, top_n)
  m <- ratio_heatmap_matrix(ratioList, top_n = top_n, rank_by = sel$rank_by %||% "variance")
  if (nrow(m) < 2) {
    # Ratios on genes (guides collapsed) and on sgRNAs share no features, for one
    showNotification(
      "The selected ratios share fewer than 2 features, so there is no Ratios heatmap.",
      type = "warning",
      duration = 8
    )
    return(NULL)
  }
  m
})

heatmapRatioStyle <- debounce(
  reactive({
    list(
      input$fb_label_by_heatmap,
      input$clusterColsHeatmap,
      input$clusterRowsHeatmap,
      input$geneNamesHeatmap,
      input$showClusterColDend,
      input$showClusterRowDend,
      input$cellBorderColourHeatmap2,
      input$cellBorderSizeHeatmap2,
      input$textSizeHeatmap,
      input$heatmapRatioLimit,
      input$heatmapColourRed,
      input$heatmapColourWhite,
      input$heatmapColourBlue,
      input$heatmapLimitCPalette,
      input$heatmapLimitCPaletteRev
    )
  }),
  1000
)

observeEvent(
  {
    heatmapRatioMatrix()
    heatmapRatioStyle()
  },
  ignoreNULL = FALSE,
  ignoreInit = FALSE,
  {
    ratioMatrix <- heatmapRatioMatrix()
    if (is.null(ratioMatrix)) {
      figuresDataReactive$heatmapRatios <- NULL
      return()
    }
    req(input$heatmapLimitCPalette)
    limit <- input$heatmapRatioLimit
    if (is.null(limit) || is.na(limit) || limit <= 0) {
      limit <- 2
    }
    rankLabel <- if (identical(input$heatmapRatioRank, "absolute")) {
      "largest absolute ratio"
    } else {
      "variance across ratios"
    }
    figuresDataReactive$heatmapRatios <- ComplexHeatmap::Heatmap(
      matrix = label_rows(ratioMatrix, input$fb_label_by_heatmap),
      name = "log2 ratio",
      cluster_columns = isTRUE(input$clusterColsHeatmap) && ncol(ratioMatrix) > 1,
      cluster_rows = isTRUE(input$clusterRowsHeatmap),
      show_column_names = TRUE,
      show_row_names = isTRUE(input$geneNamesHeatmap),
      show_column_dend = isTRUE(input$showClusterColDend),
      show_row_dend = isTRUE(input$showClusterRowDend),
      rect_gp = gpar(col = input$cellBorderColourHeatmap2, lwd = input$cellBorderSizeHeatmap2 / 10),
      use_raster = FALSE,
      row_dend_width = unit(1, "cm"),
      column_dend_height = unit(1, "cm"),
      col = .heatmap_colours(at = c(-limit, 0, limit)),
      row_names_gp = gpar(fontsize = input$textSizeHeatmap),
      column_names_gp = gpar(fontsize = input$textSizeHeatmap),
      column_title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
      column_title = paste("Top", nrow(ratioMatrix), "features by", rankLabel),
      heatmap_legend_param = list(
        labels_gp = gpar(fontsize = input$textSizeHeatmap),
        title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold")
      )
    )
  }
)

output[["dlTable_heatmapRatios"]] <- downloadHandler(
  filename = function() {
    generate_filename(base_name = "Ratio_heatmap", extension = "xlsx")
  },
  content = function(file) {
    m <- heatmapRatioMatrix()
    df <- if (is.null(m)) {
      data.frame(note = "No ratio heatmap is plotted.")
    } else {
      as.data.frame(m) %>% rownames_to_column("feature_id")
    }
    export_excel_standard(
      df = df,
      filename = file,
      sheet_name = "Ratio Heatmap",
      table_name = "RatioHeatmap"
    )
  }
)

# figure making ----
heatmapInputs <- reactive({
  list(
    input$textSizeHeatmap,
    input$heatmapBatch,
    input$heatmapDesign,
    input$lfcHeatmap,
    input$pTypeHeatmap,
    input$pThresholdHeatmap,
    input$heatmapLimitCPalette,
    input$heatmapLimitCPaletteRev,
    input$heatmapGeneNumber,
    input$heatmapGeneNumberCountQC,
    input$heatmapCounts,
    input$clusterRowsHeatmap,
    input$clusterColsHeatmap,
    input$colnamesHeatmap,
    input$geneNamesHeatmap,
    input$fb_label_by_heatmap,
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
    selectedContrast(),
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
    # Add all color picker + level label inputs as dependencies so plot updates when they change
    lapply(eds$factor_names(), function(fname) {
      lapply(eds$factor_levels(fname), function(lev) {
        list(input[[paste0("GroupColour", fname, lev)]], input[[paste0("GroupLabel", fname, lev)]])
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
    if (
      is.null(input$heatmapColumnSelect) ||
        length(input$heatmapColumnSelect) == 0 ||
        all(input$heatmapColumnSelect == "")
    ) {
      figuresDataReactive[["heatmapDEFeatures"]] <- NULL
      figuresDataReactive[["heatmapCustom"]] <- NULL
      figuresDataReactive[["heatmapGO"]] <- NULL
      figuresDataReactive[["heatmapKEGG"]] <- NULL
      for (sig in c("heatmapDEFeatures", "heatmapCustom", "heatmapGO", "heatmapKEGG")) {
        figuresDataReactive[[paste0(sig, "Table")]] <- NULL
      }
      return(NULL)
    }

    # tryCatch({
    req(input$hmKeepBucket1)
    figuresDataReactive[["heatmapDEFeatures"]] <- NULL
    # Get some variables ----
    peptide_note <- NULL # names the peptide filter when no features are left
    if (eds$omics_class() == "proteomics") {
      design <- selectedContrast()
      req(!is.null(design) && nzchar(design) && design %in% eds$contrasts())
      seqAnnoHeatmap <- eds$de(design)
      seqAnnoHeatmap[["usedInTest"]] <- TRUE
      if (!is.null(input$nrPeptideHeatmap) && "nrPeptides" %in% colnames(seqAnnoHeatmap)) {
        n_before <- nrow(seqAnnoHeatmap)
        seqAnnoHeatmap <- seqAnnoHeatmap[which(seqAnnoHeatmap$nrPeptides >= input$nrPeptideHeatmap), ]
        if (nrow(seqAnnoHeatmap) < n_before) {
          peptide_note <- paste0(
            "'Minimum number of peptides' (",
            input$nrPeptideHeatmap,
            ") keeps ",
            nrow(seqAnnoHeatmap),
            " of ",
            n_before,
            " features."
          )
        }
      }
    } else if (length(eds$contrasts()) > 0) {
      design <- eds$contrast_label(selectedContrast())
      seqAnnoHeatmap <- eds$de(selectedContrast())
    } else {
      # CountQC/CrisprCountQC — no DE results
      design <- eds$data_type()
      seqAnnoHeatmap <- NULL
    }

    # Get the desired p-value column
    if (input$pTypeHeatmap == "Raw") {
      pTypeHeatmap <- "pValue"
    } else {
      pTypeHeatmap <- "fdr"
    }

    # Prepare the metadata for the top annotation ----
    datasetHeatmap <- eds$samples()
    for (i in seq_along(input$heatmapColumnSelect)) {
      datasetHeatmap <- datasetHeatmap[
        which(datasetHeatmap[[input$heatmapColumnSelect[[i]]]] %in% input[[paste0("hmKeepBucket", i)]]),
      ]
      datasetHeatmap[[input$heatmapColumnSelect[[i]]]] <- factor(
        datasetHeatmap[[input$heatmapColumnSelect[[i]]]],
        levels = input[[paste0("hmKeepBucket", i)]]
      )
    }
    # Get the sample order
    datasetHeatmap <- datasetHeatmap %>% dplyr::arrange(across(all_of(input$heatmapColumnSelect)))

    if (nrow(datasetHeatmap) == 0) {
      return(NULL)
    }

    # The prepared count matrix (all features and samples; cached) ----
    countsHeatmap <- heatmapPrepared()$counts

    # Subset to selected samples
    countsHeatmap <- countsHeatmap[, datasetHeatmap$Name, drop = FALSE]
    datasetHeatmapOriginal <- datasetHeatmap
    ## 3. If the user wants to plot gene-wise Z-scores or centre, then do that
    if (input$heatmapZScore == "Z-Score") {
      countsHeatmap <- t(scale(t(countsHeatmap)))
    }
    if (input$heatmapZScore == "Centred") {
      countsHeatmap <- sweep(countsHeatmap, 1, rowMeans(countsHeatmap, na.rm = T))
    }

    # Column labels: sample names, unless aggregated (below)
    hmColumnLabels <- function(m) colnames(m)

    ## 4. Calculate per-condition aggregation if requested (LAST step)
    if (!is.null(input$heatmapShowMeans) && input$heatmapShowMeans) {
      # Get the grouping factor (first selected factor)
      grouping_factor <- input$heatmapColumnSelect[1]

      # Determine aggregation method
      agg_method <- if (!is.null(input$heatmapAggMethod)) input$heatmapAggMethod else "Mean"

      # Get unique factor levels that actually exist in the data
      condition_levels <- lapply(unique(datasetHeatmap[[grouping_factor]]), function(x) {
        as.character(datasetHeatmap$Name[datasetHeatmap[[grouping_factor]] == x])
      }) %>%
        setNames(unique(datasetHeatmap[[grouping_factor]]))

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

      # Aggregated columns are condition levels, but data.frame() above
      # mangled their names (make.names): label them with the Colours-tab
      # display label, falling back to the raw level.
      aggLevels <- names(condition_levels)
      aggLabels <- level_labels(input, grouping_factor, aggLevels)
      aggRich <- is_rich_labels(aggLabels)
      aggLabelMap <- setNames(aggLabels, make.names(aggLevels, unique = TRUE))
      hmColumnLabels <- function(m) {
        out <- unname(aggLabelMap[colnames(m)])
        fallback <- if (aggRich) escape_label_markdown(colnames(m)) else colnames(m)
        out <- ifelse(is.na(out), fallback, out)
        if (aggRich) ComplexHeatmap::gt_render(out) else out
      }

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
          ann <- aggregated_annotation(datasetHeatmapOriginal, grouping_factor, f, names(condition_levels))
          if (ann$varies) {
            # Show warning and stop - cannot annotate with this factor
            showNotification(
              paste0(
                "Cannot add '",
                f,
                "' to annotation: it has multiple values within '",
                grouping_factor,
                "' groups when showing aggregated data (mean/median). Please uncheck 'Show Mean/Median' or remove '",
                f,
                "' from the annotation."
              ),
              type = "error",
              duration = 10
            )
            # Skip adding this factor to prevent crash
            next
          } else {
            # Factor is constant within groups - safe to annotate
            datasetHeatmap[[f]] <- factor(unname(ann$values))
          }
        }
      }
    }

    # Prepare colours and generate top annotation ----
    heatmapColours <- .heatmap_colours()
    # Top Annotation
    colourListHeatmap <- lapply(input$heatmapColumnSelect, function(f) {
      vapply(
        levels(datasetHeatmap[[f]]),
        function(l) paste0(col2hex(input[[paste0("GroupColour", f, l)]]), "FF"),
        character(1)
      )
    }) %>%
      setNames(input$heatmapColumnSelect)
    if (any(input$heatmapColumnSelect != "")) {
      # Colours-tab display labels. Rich (formatted) labels get legends built
      # at draw time by draw_heatmap() (ComplexHeatmap mis-sizes gt_render
      # legend labels, see R/level_labels.R); otherwise the annotation's own
      # legends show the plain labels.
      annoFactors <- intersect(input$heatmapColumnSelect, colnames(datasetHeatmap))
      annoLevels <- lapply(setNames(annoFactors, annoFactors), function(g) levels(as.factor(datasetHeatmap[[g]])))
      annoRich <- any(vapply(
        annoFactors,
        function(g) is_rich_labels(level_labels(input, g, annoLevels[[g]])),
        logical(1)
      ))
      annoLegendParam <- lapply(annoFactors, function(g) {
        levs <- levels(as.factor(datasetHeatmap[[g]]))
        list(
          title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
          labels_gp = gpar(fontsize = input$textSizeHeatmap),
          at = levs,
          labels = unname(level_labels(input, g, levs, mode = "plain"))
        )
      }) %>%
        setNames(annoFactors)
      ha <- heatmap_top_annotation(
        datasetHeatmap,
        annoFactors,
        colourListHeatmap,
        legend_param = annoLegendParam,
        gp = gpar(col = input$cellBorderColourHeatmap1, lwd = input$cellBorderSizeHeatmap1 / 10),
        show_legend = !annoRich
      )
      if (annoRich) {
        attr(ha, "rich_legends") <- lapply(annoFactors, function(g) {
          heatmap_legend_spec(
            title = g,
            labels = level_labels(input, g, annoLevels[[g]], force_rich = TRUE),
            fill = unlist(colourListHeatmap[[g]]),
            title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
            labels_gp = gpar(fontsize = input$textSizeHeatmap)
          )
        })
      }
    }

    # DE Heatmaps ----
    # Select genes based on data type. Counts-only uploads (RNASeq/proteomics
    # with no DE) have no seqAnno to threshold on, so fall back to SD ranking
    # exactly like CountQC.
    if (eds$data_type() %in% c("CountQC", "CrisprCountQC") || !eds$has_de()) {
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

    # Prepare heatmap name (used by all heatmaps) — pretty label for the legend
    hmName <- eds$assay_label(input$heatmapCounts) %>% gsub("\\+", "\n\\+", .)
    if (input$heatmapLog2 && !eds$assay_is_log2(input$heatmapCounts)) {
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

    # The DE, Custom, GO and KEGG heatmaps differ only in matrix, title and
    # cell border size
    count_heatmap <- function(mat, title, border_size = input$cellBorderSizeHeatmap2) {
      ComplexHeatmap::Heatmap(
        matrix = mat,
        top_annotation = ha,
        name = hmName,
        cluster_columns = input$clusterColsHeatmap,
        cluster_rows = input$clusterRowsHeatmap,
        show_column_names = input$colnamesHeatmap,
        column_labels = hmColumnLabels(mat),
        show_row_names = input$geneNamesHeatmap,
        show_column_dend = input$showClusterColDend,
        show_row_dend = input$showClusterRowDend,
        rect_gp = gpar(col = input$cellBorderColourHeatmap2, lwd = border_size / 10),
        use_raster = FALSE,
        row_dend_width = unit(1, "cm"),
        column_dend_height = unit(1, "cm"),
        col = heatmapColours,
        row_names_gp = gpar(fontsize = input$textSizeHeatmap),
        column_names_gp = gpar(fontsize = input$textSizeHeatmap),
        column_title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
        column_title = title,
        heatmap_legend_param = list(
          labels_gp = gpar(fontsize = input$textSizeHeatmap),
          title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold")
        )
      )
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
      if (eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
        max_features <- min(as.numeric(input$heatmapGeneNumberCountQC), nrow(countsHeatmapSig))
        countsHeatmapSig <- countsHeatmapSig[seq_len(max_features), , drop = FALSE]
      } else {
        if (as.numeric(input$heatmapGeneNumber) <= nrow(countsHeatmapSig)) {
          countsHeatmapSig <- countsHeatmapSig[seq_len(as.numeric(input$heatmapGeneNumber)), , drop = FALSE]
        }
      }

      # Determine heatmap title based on data type
      heatmap_title <- if (eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
        paste("Top", nrow(countsHeatmapSig), "Features by SD", agg_label, ":\n", design)
      } else {
        paste("Top", nrow(countsHeatmapSig), "DE Features", agg_label, ":\n", design)
      }

      # For heatmaps with many features (>200), use no cell border to avoid clutter
      cell_border_size <- if (nrow(countsHeatmapSig) > 200) 0 else input$cellBorderSizeHeatmap2

      heatmap <- count_heatmap(
        # Row labels follow the "Label features by" radio (proteomics); the
        # matrix itself keeps feature_ids for the download table.
        label_rows(countsHeatmapSig, input$fb_label_by_heatmap, width = 15),
        heatmap_title,
        border_size = cell_border_size
      )
      figuresDataReactive[["heatmapDEFeatures"]] <- heatmap

      figuresDataReactive[["heatmapDEFeaturesTable"]] <- list(counts = countsHeatmapSig, design = design)
      # Show number of genes on plot (changes depending on settings applied)
      output[["hGNDE"]] <- renderText({
        paste0("Number of features per settings chosen: ", nrow(countsHeatmapSig))
      })
    } else {
      # When no DEGs are found (ntp = 0), provide informative placeholders
      figuresDataReactive[["heatmapDEFeatures"]] <- NULL
      figuresDataReactive[["heatmapDEFeaturesTable"]] <- NULL

      output[["heatmapDesignBoth Directions"]] <- renderText({
        paste("No features found meeting current thresholds for:", design, peptide_note)
      })

      output[["hGNDE"]] <- renderText({
        "Number of features per settings chosen: 0"
      })
    } # End if (ntp > 0) - DE Features heatmap section

    # Custom heatmap ----
    countsHeatmapCustom <- countsHeatmap[input$fb_keep_heatmap, , drop = FALSE]
    ch <- count_heatmap(
      label_rows(countsHeatmapCustom, input$fb_label_by_heatmap),
      paste(input$heatmapCustomTitle, agg_label)
    )
    figuresDataReactive$heatmapCustom <- ch

    figuresDataReactive[["heatmapCustomTable"]] <- list(counts = countsHeatmapCustom, design = design)

    # GO heatmap ----
    if (eds$data_type() %in% c("RNASeq", "CountQC", "CrisprCountQC")) {
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
            which(
              seqAnnoHeatmap[[pCol]] <= as.numeric(input$heatmapGoPThreshold) &
                abs(seqAnnoHeatmap$log2Ratio) >= as.numeric(input$heatmapGoLFC)
            )
          ]
          goGenes <- goGenes[goGenes %in% deFeatures]
        }

        countsHeatmapGO <- heatmap_term_matrix(countsHeatmap, goGenes, input$heatmapGoExpr)
      }

      if (!is.null(countsHeatmapGO) && nrow(countsHeatmapGO) > 0) {
        req(nrow(datasetHeatmap) == ncol(countsHeatmapGO))
        gh <- count_heatmap(
          countsHeatmapGO,
          paste0(
            goHeatmapReactive$currentGoName %||% goHeatmapReactive$currentGoTerm %||% "",
            "\n",
            nrow(countsHeatmapGO),
            "/",
            goGenesTotal,
            " features plotted",
            if (agg_label != "") paste0("\n", agg_label) else ""
          )
        )
        figuresDataReactive$heatmapGO <- gh
        figuresDataReactive[["heatmapGOTable"]] <- list(counts = countsHeatmapGO, design = design)
      } else {
        figuresDataReactive$heatmapGO <- NULL
        figuresDataReactive[["heatmapGOTable"]] <- NULL
      }
    }

    # KEGG heatmap ----
    if (eds$data_type() %in% c("RNASeq", "CountQC", "CrisprCountQC") && !is.null(keggHeatmapReactive$keggCode)) {
      keggGenes <- keggHeatmapReactive$currentGenes
      countsHeatmapKEGG <- NULL

      if (!is.null(keggGenes) && length(keggGenes) > 0) {
        keggGenes <- keggGenes[which(keggGenes %in% rownames(countsHeatmap))]

        # Apply DE filter if enabled
        if (!is.null(input$heatmapKeggDEOnly) && input$heatmapKeggDEOnly) {
          pCol <- if (input$heatmapKeggPType == "Raw") "pValue" else "fdr"
          deFeatures <- seqAnnoHeatmap$feature_id[
            which(
              seqAnnoHeatmap[[pCol]] <= as.numeric(input$heatmapKeggPThreshold) &
                abs(seqAnnoHeatmap$log2Ratio) >= as.numeric(input$heatmapKeggLFC)
            )
          ]
          keggGenes <- keggGenes[keggGenes %in% deFeatures]
        }

        countsHeatmapKEGG <- heatmap_term_matrix(countsHeatmap, keggGenes, input$heatmapKeggExpr)

        keggGenesTotal <- length(keggHeatmapReactive$currentGenes)

        if (!is.null(countsHeatmapKEGG) && nrow(countsHeatmapKEGG) > 0) {
          req(nrow(datasetHeatmap) == ncol(countsHeatmapKEGG))

          keggTitle <- paste0(
            keggHeatmapReactive$currentPathwayName %||% "",
            "\n",
            nrow(countsHeatmapKEGG),
            "/",
            keggGenesTotal,
            " features plotted",
            if (agg_label != "") paste0("\n", agg_label) else ""
          )

          kh <- count_heatmap(countsHeatmapKEGG, keggTitle)
          figuresDataReactive$heatmapKEGG <- kh

          figuresDataReactive[["heatmapKEGGTable"]] <- list(counts = countsHeatmapKEGG, design = design)
        } else {
          figuresDataReactive$heatmapKEGG <- NULL
          figuresDataReactive[["heatmapKEGGTable"]] <- NULL
        }
      }
    }
    # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
)

# GO heatmap: fetch GO term genes on button click
observeEvent(input$heatmapGoGo, {
  req(input$heatmapGoInput)
  req(length(eds$contrasts()) > 0)

  goId <- gsub(" .*", "", input$heatmapGoInput)
  goName <- sub("^\\S+\\s+", "", input$heatmapGoInput)
  # GO membership from features() (always carries the GO columns; .h5ad DE
  # tables don't). rownames(features()) are the DE table's feature_ids.
  feats <- eds$features()
  goGenes <- rownames(feats)[features_in_go_term(feats, goId)]

  goHeatmapReactive$currentGoTerm <- goId
  goHeatmapReactive$currentGoName <- goName
  goHeatmapReactive$currentGenes <- goGenes
})

# GO heatmap: display GO term info
output$heatmapGoInfo <- renderUI({
  req(goHeatmapReactive$currentGoTerm)
  tagList(
    hr(),
    tags$b("GO term: "),
    goHeatmapReactive$currentGoName,
    br(),
    tags$b("Genes found in dataset: "),
    length(goHeatmapReactive$currentGenes)
  )
})

# The table downloads of the DE, Custom, GO and KEGG heatmaps: the plotted
# values the drawing observer stores next to each figure (`<sig>Table`; none
# while a heatmap is empty), with the selected contrast's DE statistics when
# "Include DE statistics" is ticked (heatmap_export_table()).
heatmap_table_download <- function(sig, base_name, sheet_name, table_name) {
  downloadHandler(
    filename = function() {
      design <- figuresDataReactive[[paste0(sig, "Table")]]$design %||% ""
      generate_filename(base_name = base_name, design_name = design, extension = "xlsx")
    },
    content = function(file) {
      counts <- figuresDataReactive[[paste0(sig, "Table")]]$counts %||% matrix(numeric(0), 0, 0)
      de <- if (eds$has_de() && !is.null(selectedContrast())) eds$de(selectedContrast())
      export_excel_standard(
        df = heatmap_export_table(counts, de, input$includeDEStatsHeatmap),
        filename = file,
        sheet_name = sheet_name,
        table_name = table_name
      )
    }
  )
}
output[["dlTable_heatmapDE"]] <- heatmap_table_download(
  "heatmapDEFeatures",
  "DE_heatmap_counts",
  "Heatmap Data",
  "HeatmapData"
)
output[["dlTable_heatmapCustom"]] <- heatmap_table_download(
  "heatmapCustom",
  "Custom_heatmap",
  "Custom Heatmap",
  "CustomHeatmap"
)
output[["dlTable_heatmapGO"]] <- heatmap_table_download("heatmapGO", "GO_heatmap", "GO Heatmap", "GOHeatmap")
output[["dlTable_heatmapKEGG"]] <- heatmap_table_download("heatmapKEGG", "KEGG_heatmap", "KEGG Heatmap", "KEGGHeatmap")

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
  fetched <- withProgress(message = "Fetching KEGG pathway genes...", {
    tryCatch(
      kegg_fetch_pathway(pathwayId, keggCode),
      error = function(e) {
        showNotification(paste("Failed to fetch KEGG pathway:", e$message), type = "error", duration = 8)
        NULL
      }
    )
  })
  result <- if (!is.null(fetched)) {
    list(genes = kegg_map_genes(eds$de(eds$primary_contrast()), fetched$entrez, keggCode), name = fetched$name)
  }

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
    tags$b("Pathway: "),
    keggHeatmapReactive$currentPathwayName,
    br(),
    tags$b("Genes found in dataset: "),
    length(keggHeatmapReactive$currentGenes)
  )
})

# Output heatmaps
lapply(c("DEFeatures", "Custom", "GO", "KEGG", "Ratios"), function(sig) {
  output[[paste0("heatmap", sig)]] <- renderPlot(
    {
      req(figuresDataReactive[[paste0("heatmap", sig)]])
      draw_heatmap(
        figuresDataReactive[[paste0("heatmap", sig)]],
        merge_legend = TRUE,
        padding = unit(c(10, 20, 10, 10), "mm"),
        legend_dx = input$heatmapLegendDX,
        legend_dy = input$heatmapLegendDY
      )
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
  button_suffix <- switch(
    sig,
    "DEFeatures" = "DE",
    "Custom" = "Custom",
    "GO" = "GO",
    "KEGG" = "KEGG",
    "Ratios" = "Ratios"
  )

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
          text(
            0.5,
            0.5,
            paste("No data available for", sig, "heatmap\n\nPlease adjust filter settings or select features."),
            cex = 1.2
          )
          dev.off()
        } else {
          writeLines(
            paste("No data available for", sig, "heatmap. Please adjust filter settings or select features."),
            file
          )
        }
        return()
      }

      # Heatmap exists, proceed with normal download
      # Using a custom draw wrapper because export_plot_standard expects a print/plot call,
      # and ComplexHeatmap needs draw() with specific params.
      export_plot_standard(
        plot_object = function() {
          draw_heatmap(
            heatmap_obj,
            merge_legend = TRUE,
            padding = unit(c(10, 20, 10, 10), "mm"),
            legend_dx = input$heatmapLegendDX,
            legend_dy = input$heatmapLegendDY
          )
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
