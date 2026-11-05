# ---------------------------------------------------------------------------
# Phase 2 migration (2026-05-15) — partial.
# Second server-*.R file moved to the new ExploreDataset API. Same pattern
# as server-correlation.R.
# Migrated to eds$method() calls:
#   $factors, $dataset, $sgRNA_to_gene, $dataType, $countList, $design
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

# Update inputs
observeEvent(inputDataReactive(), {
  eds <- inputDataReactive()$.eds
  assay_keys <- eds$assay_names()
  count_choices <- setNames(assay_keys, vapply(assay_keys, eds$assay_label, character(1)))
  default_count <- if (eds$data_type() == "RNASeq" && "normalised_log2" %in% assay_keys) {
    "normalised_log2"
  } else if (eds$omics_class() == "proteomics" && "intensities" %in% assay_keys) {
    "intensities"
  } else if (eds$data_type() %in% c("CountQC", "CrisprCountQC") && "vst" %in% assay_keys) {
    "vst"
  } else {
    assay_keys[1]
  }

  updateSelectInput(
    session = session,
    inputId = "replicabilityCounts",
    choices = count_choices,
    selected = default_count
  )
})

output$replicabilityFactorSelectUI <- renderUI({
  eds <- inputDataReactive()$.eds
  fs <- eds$factor_names()
  selectizeInput(
    inputId = "replicabilityFactors",
    label = "Select factors for annotation",
    choices = fs,
    multiple = TRUE,
    selected = fs[1]
  )
})

# Update break values based on correlation data
observeEvent(replicabilityData(), {
  req(replicabilityData())

  # Get all correlation data
  all_corr <- replicabilityData()$all$corr
  top100_corr <- replicabilityData()$top100$corr

  # Use the min/mean/max from all genes
  corr_min <- min(all_corr, na.rm = TRUE)
  corr_mean <- mean(all_corr, na.rm = TRUE)
  corr_max <- max(all_corr, na.rm = TRUE)

  # Update the numeric inputs
  # save-state restore: yield the saved break while restoring (this observer
  # fires on every data recompute, which a restore triggers).
  updateNumericInput(session, "replicabilityAtLow", value = restore_default("replicabilityAtLow", round(corr_min, 2)))
  updateNumericInput(session, "replicabilityAtMid", value = restore_default("replicabilityAtMid", round(corr_mean, 2)))
  updateNumericInput(session, "replicabilityAtHigh", value = restore_default("replicabilityAtHigh", round(corr_max, 2)))
})

# Reactive for calculating correlation and clustering
replicabilityData <- reactive({
  req(input$replicabilityCounts)
  req(inputDataReactive())
  eds <- inputDataReactive()$.eds

  counts <- eds$assay(input$replicabilityCounts)
  # Collapse sgRNA guides to gene level if requested
  if (eds$data_type() == "CrisprCountQC" && isTRUE(input$replicabilityCollapseGuides)) {
    counts <- collapse_guides_to_genes(counts, eds$sgrna_to_gene())
  }

  # Apply Log2 if requested
  if (isTRUE(input$replicabilityCountsLog)) {
    counts <- log2(counts + 1)
  }

  # Helper to calculate corr and hclust
  calc_corr_hclust <- function(cnts) {
    # Correlation
    # Use pairwise.complete.obs to handle NAs (e.g. in proteomics)
    corr_mat <- cor(cnts, method = input$replicabilityCorrMethod, use = "pairwise.complete.obs")

    # Clustering
    # Transpose counts for dist (samples as rows)
    # Handle NAs in dist if necessary, but dist usually propagates NAs or requires imputation.
    # For now, we rely on standard dist behavior but could impute if needed.
    # Replacing NAs with 0 for clustering might be safer for proteomics if NAs are missing values.
    cnts_for_dist <- cnts
    cnts_for_dist[is.na(cnts_for_dist)] <- 0

    dist_mat <- dist(t(cnts_for_dist), method = input$replicabilityDistMethod)
    hclust_res <- hclust(dist_mat, method = input$replicabilityClustMethod)
    list(corr = corr_mat, hclust = hclust_res)
  }

  # 1. All Genes (filtered by variance)
  counts_all <- counts
  if (nrow(counts_all) > input$replicabilityTopGenes) {
    vars <- apply(counts_all, 1, var, na.rm = TRUE)
    counts_all <- counts_all[order(vars, decreasing = TRUE)[1:input$replicabilityTopGenes], ]
  }
  res_all <- calc_corr_hclust(counts_all)

  # 2. Top 100 Abundant Genes
  counts_top100 <- counts
  means <- rowMeans(counts_top100, na.rm = TRUE)
  # Take top 100 most abundant
  top_n <- min(100, nrow(counts_top100))
  counts_top100 <- counts_top100[order(means, decreasing = TRUE)[1:top_n], ]
  res_top100 <- calc_corr_hclust(counts_top100)

  list(
    all = list(corr = res_all$corr, hclust = res_all$hclust, dataset = eds$samples()),
    top100 = list(corr = res_top100$corr, hclust = res_top100$hclust, dataset = eds$samples())
  )
})

# Defined in R/helperFunctions.R; bound to this scope so they see this session's input.
make_rep_heatmap <- bind_to_scope(make_rep_heatmap)
make_rep_dendrogram <- bind_to_scope(make_rep_dendrogram)

# Render Heatmaps
output$replicabilityHeatmap <- renderPlot(
  {
    make_rep_heatmap(replicabilityData()$all, title_suffix = "All Features")
  },
  width = function() input$figWidthReplicability,
  height = function() input$figHeightReplicability
)

output$replicabilityHeatmapTop100 <- renderPlot(
  {
    make_rep_heatmap(replicabilityData()$top100, title_suffix = "Top 100 Abundant Features")
  },
  width = function() input$figWidthReplicability,
  height = function() input$figHeightReplicability
)

# Render Dendrograms
output$replicabilityDendrogram <- renderPlot(
  {
    make_rep_dendrogram(replicabilityData()$all, title_suffix = "All Features")
  },
  width = function() input$figWidthReplicability,
  height = function() input$figHeightReplicability
)

output$replicabilityDendrogramTop100 <- renderPlot(
  {
    make_rep_dendrogram(replicabilityData()$top100, title_suffix = "Top 100 Abundant Features")
  },
  width = function() input$figWidthReplicability,
  height = function() input$figHeightReplicability
)

# Download Handlers
output$replicabilityDownloadHeatmap <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_replicability, "_Heatmap"),
      design_name = designReactive(),
      suffix = gsub(" ", "", input$replicabilityDownloadSet),
      extension = tolower(input$downloadFormatReplicability)
    )
  },
  content = function(file) {
    req(replicabilityData())

    # Select dataset
    data_set <- if (input$replicabilityDownloadSet == "All Features") {
      replicabilityData()$all
    } else {
      replicabilityData()$top100
    }

    export_plot_standard(
      plot_object = function() {
        make_rep_heatmap(data_set, title_suffix = input$replicabilityDownloadSet)
      },
      filename = file,
      format = tolower(input$downloadFormatReplicability),
      width = as.numeric(input$figWidthReplicability / 72),
      height = as.numeric(input$figHeightReplicability / 72),
      dpi = as.numeric(input$dpiReplicability)
    )
  }
)

output$replicabilityDownloadDendrogram <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_replicability, "_Dendrogram"),
      design_name = designReactive(),
      suffix = gsub(" ", "", input$replicabilityDownloadSet),
      extension = tolower(input$downloadFormatReplicability)
    )
  },
  content = function(file) {
    req(replicabilityData())

    # Select dataset
    data_set <- if (input$replicabilityDownloadSet == "All Features") {
      replicabilityData()$all
    } else {
      replicabilityData()$top100
    }

    export_plot_standard(
      plot_object = function() {
        make_rep_dendrogram(data_set, title_suffix = input$replicabilityDownloadSet)
      },
      filename = file,
      format = tolower(input$downloadFormatReplicability),
      width = as.numeric(input$figWidthReplicability / 72),
      height = as.numeric(input$figHeightReplicability / 72),
      dpi = as.numeric(input$dpiReplicability)
    )
  }
)

output$replicabilityDownloadCorr <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_replicability, "_Correlation"),
      design_name = designReactive(),
      suffix = gsub(" ", "", input$replicabilityDownloadSet),
      extension = "xlsx"
    )
  },
  content = function(file) {
    req(replicabilityData())

    # Select dataset
    data_set <- if (input$replicabilityDownloadSet == "All Features") {
      replicabilityData()$all
    } else {
      replicabilityData()$top100
    }

    export_excel_standard(
      df = as.data.frame(data_set$corr) %>% rownames_to_column("Sample"),
      filename = file,
      sheet_name = "Correlation Matrix",
      table_name = "CorrelationMatrix"
    )
  }
)
