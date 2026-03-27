# Check if GSEA tab is accessed without results
observeEvent(input$tabs, {
  if (input$tabs == "gseaTab") {
    if (is.null(inputDataReactive()$goObjects$gsea) ||
      length(inputDataReactive()$goObjects$gsea) == 0) {
      showModal(modalDialog(
        title = "No GSEA Results",
        "No pathway analysis was performed for this dataset, so this tab is empty.",
        footer = modalButton("OK")
      ))
      return()
    }

    # Check if current selection has no results
    if (!is.null(input$gseaType) &&
      !is.null(inputDataReactive()$goObjects$gsea[[input$gseaType]])) {
      er <- inputDataReactive()$goObjects$gsea[[input$gseaType]]
      if (is.null(er) || nrow(er@result) == 0) {
        showModal(modalDialog(
          title = "No GSEA Results",
          "There are no results for this selection."
        ))
      }
    }
  }
})

# Initialize reactive values for GSEA plots
gseaPlotReactive <- reactiveValues(
  er = NULL,
  cnetplot = NULL,
  heatmap = NULL,
  runningscore = NULL,
  ridgeplot = NULL,
  upset = NULL
)


# Get data objects
gsea <- inputDataReactive()$goObjects$gsea
se <- inputDataReactive()$se
design <- inputDataReactive()$design
param <- inputDataReactive()$param

# Update GSEA type choices
updateSelectInput(
  session = session,
  inputId = "gseaType",
  choices = names(gsea),
  selected = names(gsea)[[1]]
)

# Render design output
output$gseaDesign <- renderText({
  design
})

# Set default FDR threshold if missing
if (is.null(param$fdrThreshGSEA)) {
  param$fdrThreshGSEA <- 0.05
}

# Create GSEA overview table
gsea.overview <- ezFrame("Number of Significant Pathways" = numeric())
for (i in seq_along(names(gsea))) {
  testName <- names(gsea)[i]
  gsea.overview[testName, "Number of Significant Pathways"] <- nrow(gsea[[i]]@result)
}

output$gseaOverview <- function() {
  gsea.overview %>%
    kable(format = "html", caption = "GSEA:") %>%
    kable_styling(
      bootstrap_options = "striped",
      full_width = FALSE,
      position = "left"
    )
}

# Create GSEA thresholds table
gseaThresholds <- ezFrame("GSEA Cut-Offs" = character())
gseaThresholds["Output GSEA Term FDR", "GSEA Cut-Offs"] <- param$fdrThreshGSEA

output$gseaThreshOverview <- function() {
  gseaThresholds %>%
    kable(format = "html", caption = "Cut-Offs:") %>%
    kable_styling(
      bootstrap_options = "striped",
      full_width = FALSE,
      position = "left"
    )
}

# Isoform-level aggregation UI (only shown for isoform data)
if (!is.null(param$featureLevel) && param$featureLevel == "isoform") {
  output$gseaIsoformAggUI <- renderUI({
    tagList(
      box(
        title = "Isoform Aggregation",
        width = NULL,
        solidHeader = TRUE,
        status = "info",
        collapsible = TRUE,
        collapsed = FALSE,
        tags$div(
          class = "callout callout-info",
          tags$p(
            tags$b("Note:"), "This dataset is at isoform level. Multiple transcripts per gene
            were passed to the enrichment analysis, each with their own log2 fold change. To
            display a single fold change per gene in the network plot, values are aggregated
            using the method selected below."
          )
        ),
        selectInput(
          inputId = "gseaIsoformAgg",
          label = "Aggregation method",
          choices = c("Max absolute" = "maxabs", "Mean" = "mean", "Median" = "median", "Max" = "max", "Min" = "min"),
          selected = "maxabs",
          width = "50%"
        )
      )
    )
  })
  outputOptions(output, "gseaIsoformAggUI", suspendWhenHidden = FALSE)
}

# Reactive values for current GSEA selection
gseaToUse <- reactiveValues(er = NULL)

# Main observer for GSEA type changes
observeEvent(
  input$gseaType,
  {
    gseaToUse$er <- NULL
    gseaToUse$er_symbol <- NULL
    er <- gsea[[input$gseaType]]

    # Check if results exist and show modal if empty
    if (input$tabs == "gseaTab" && (is.null(er) || nrow(er@result) == 0)) {
      showModal(modalDialog(
        title = "No GSEA Results",
        "There are no results for this selection."
      ))
      return()
    }

    # Process results if they exist
    if (!is.null(er) && nrow(er@result) > 0) {
      # Handle GO terms (BP, MF, CC) differently
      if (input$gseaType %in% c("BP", "MF", "CC")) {
        er@result$Label <- ""

        # Create labels for GO terms
        for (i in seq_along(er@result$Description)) {
          desc <- er@result$Description[i]
          label <- paste0(er@result$ID[i], "\n", desc)
          if (!is.na(label) && nchar(label) > 45) {
            label <- substr(label, 1, 45)
          }
          er@result$Label[i] <- label
        }
        colsToTable <- c("Label", "enrichmentScore", "NES", "qvalue", "geneName")
        # Convert gene IDs to gene names for core enrichment
        # Split all core enrichment strings and collect unique gene IDs
        all_core_genes <- er@result$core_enrichment %>%
          strsplit("/") %>%
          unlist() %>%
          unique()
        all_core_genes <- all_core_genes[nchar(all_core_genes) > 0]

        # Create lookup table for gene ID to gene name conversion
        gene_lookup <- setNames(
          metadata(se)$enrichInput$seqAnno$gene_name,
          metadata(se)$enrichInput$seqAnno$gene_id
        )
        er_symbol <- er
        # Convert core enrichment from gene IDs to gene names
        er_symbol@result$core_enrichment <- vapply(er_symbol@result$core_enrichment, function(genes_string) {
          gene_ids <- strsplit(genes_string, "/")[[1]]
          gene_names <- gene_lookup[gene_ids]
          gene_names <- gene_names[!is.na(gene_names) & nchar(gene_names) > 0]
          paste(gene_names, collapse = "/")
        }, character(1), USE.NAMES = FALSE)
      } else {
        # For other pathway databases (e.g., C2)
        er_symbol <- er
        colsToTable <- c("Description", "enrichmentScore", "NES", "qvalue", "core_enrichment")
      }
      # Prepare log2 ratio data for plotting (store raw, aggregate in plot observer)
      log2RatioGSEA <- metadata(se)$enrichInput$log2Ratio
      names(log2RatioGSEA) <- rowData(se)$gene_name[match(names(log2RatioGSEA), rowData(se)$gene_id)]
      log2RatioGSEA <- sort(log2RatioGSEA, decreasing = TRUE)

      # Update reactive values
      gseaToUse$er <- er
      gseaToUse$er_symbol <- er_symbol
      gseaToUse$log2RatioGSEA <- log2RatioGSEA
      gseaToUse$colsToTable <- colsToTable
    }


    # Clear existing plots
    for (x in names(gseaPlotReactive)) {
      gseaPlotReactive[[x]] <- NULL
    }

    # Reset table selection AFTER updating reactive values
    shinyjs::delay(100, {
      proxy <- DT::dataTableProxy("selectedTable_GSEA")
      DT::selectRows(proxy, NULL)
    })
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE
)

# Render the main GSEA table
output$selectedTable_GSEA <- DT::renderDataTable({
  req(gseaToUse$er, gseaToUse$colsToTable)

  DT::datatable(
    data = as.data.frame(gseaToUse$er)[, gseaToUse$colsToTable],
    filter = "top",
    class = "cell-border stripe",
    rownames = FALSE,
    caption = "Click pathways in this table to add them to the figures on the right"
  ) %>%
    DT::formatStyle(
      columns = colnames(.$x$data),
      `font-size` = "12px"
    ) %>%
    DT::formatSignif(
      columns = c("enrichmentScore", "NES", "qvalue"),
      digits = 3
    ) %>%
    DT::formatStyle(
      columns = c("enrichmentScore", "NES"),
      color = styleInterval(cuts = 0, values = c("blue", "darkorange")),
      fontWeight = "bold"
    ) %>%
    DT::formatStyle(
      columns = "qvalue",
      color = styleInterval(cuts = param$fdrThreshGSEA, values = c("green", "black")),
      fontWeight = "bold"
    )
})

# Download handler for GSEA results
# Download handler for GSEA results
output$dlTable_GSEA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = "GSEA_Results",
      design_name = design,
      extension = "xlsx"
    )
  },
  content = function(file) {
    sheets_list <- list()
    categories <- c("BP", "CC", "MF")
    for (category in categories) {
      if (!is.null(gsea[[category]])) {
        enrich_obj <- gsea[[category]]
        results_df <- as.data.frame(enrich_obj@result)
        sheet_name <- paste(category, sep = "_")
        sheets_list[[sheet_name]] <- results_df
      }
    }
    # Handle custom GSEA types if present (e.g., C2)
    other_types <- setdiff(names(gsea), categories)
    for (type in other_types) {
      enrich_obj <- gsea[[type]]
      results_df <- as.data.frame(enrich_obj@result)
      sheets_list[[type]] <- results_df
    }

    if (length(sheets_list) > 0) {
      write_xlsx(sheets_list, file)
    } else {
      # Fallback if no sheets were added
      write.csv(data.frame(Message = "No GSEA results available"), file)
    }
  }
)

# Observer for plot generation when rows are selected
observeEvent(
  {
    input$textSizeGSEA
    input$nodeSizeGSEA
    input$selectedTable_GSEA_rows_selected
    input$scaleLimGSEA
    input$gseaMaxOver
    input$gseaDodgeY
    input$gseaDodgeX
    input$nodeBorderGSEA
    input$ridgePlotColourByGSEA
    input$gseaLabelAlpha
    input$gseaVisPType
    input$gseaVisPThreshold
    input$gseaVisLFC
    input$gseaVisFilterMode
    input$gseaIsoformAgg
  },
  {
    req(gseaToUse$er_symbol)
    current_er <- isolate(gseaToUse$er)
    current_er_symbol <- isolate(gseaToUse$er_symbol)
    req(!is.null(input$selectedTable_GSEA_rows_selected))
    req(current_er)

    # Subset enrichment result for selected pathways
    erRP <- clusterProfiler::slice(current_er, input$selectedTable_GSEA_rows_selected)
    labelColumn <- ifelse(input$gseaType == "C2", "Description", "Label")

    # Apply gene filtering based on visualisation thresholds
    log2RatioGSEA <- gseaToUse$log2RatioGSEA
    # For isoform data, aggregate multiple transcripts per gene
    if (!is.null(param$featureLevel) && param$featureLevel == "isoform") {
      agg_method <- if (!is.null(input$gseaIsoformAgg)) input$gseaIsoformAgg else "maxabs"
      na_safe <- function(fn) function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else fn(x) }
      agg_fn <- switch(agg_method,
        "maxabs" = na_safe(function(x) x[which.max(abs(x))]),
        "mean"   = na_safe(mean),
        "median" = na_safe(median),
        "max"    = na_safe(max),
        "min"    = na_safe(min)
      )
      log2RatioGSEA <- vapply(split(log2RatioGSEA, names(log2RatioGSEA)), agg_fn, numeric(1))
      log2RatioGSEA <- sort(log2RatioGSEA, decreasing = TRUE)
    }
    log2RatioGSEA_full <- log2RatioGSEA # Keep full list for reference
    genes_pass_filter <- names(log2RatioGSEA) # Default: all genes pass

    if (!is.null(input$gseaVisFilterMode) && input$gseaVisFilterMode != "none") {
      # Get p-values for genes — try enrichInput path first, then fall back to direct path
      seqAnno <- metadata(se)$enrichInput$seqAnno
      if (is.null(seqAnno)) {
        seqAnno <- metadata(se)$seqAnno
      }

      # Determine which p-value column to use
      pval_col <- if (!is.null(input$gseaVisPType) && input$gseaVisPType == "FDR") "fdr" else "pValue"

      # Validate that the column exists
      if (!is.null(seqAnno) && pval_col %in% colnames(seqAnno) && "gene_name" %in% colnames(seqAnno)) {
        # Match genes by name
        gene_match_idx <- match(names(log2RatioGSEA), seqAnno$gene_name)
        gene_match_idx <- gene_match_idx[!is.na(gene_match_idx)]

        if (length(gene_match_idx) > 0) {
          # Get p-values for matched genes
          gene_pvals <- seqAnno[[pval_col]][gene_match_idx]
          names(gene_pvals) <- seqAnno$gene_name[gene_match_idx]

          # Filter genes based on thresholds
          genes_pass_filter <- names(log2RatioGSEA)[
            !is.na(gene_pvals[names(log2RatioGSEA)]) &
              gene_pvals[names(log2RatioGSEA)] <= input$gseaVisPThreshold &
              abs(log2RatioGSEA) >= input$gseaVisLFC
          ]
          genes_pass_filter <- genes_pass_filter[!is.na(genes_pass_filter)]

          # Apply filter mode
          if (input$gseaVisFilterMode == "genes" && length(genes_pass_filter) > 0) {
            # Filter genes entirely - only keep genes meeting thresholds
            log2RatioGSEA <- log2RatioGSEA[names(log2RatioGSEA) %in% genes_pass_filter]
          }
        }
      }
      # For "labels" mode, keep all genes in log2RatioGSEA but use genes_pass_filter for labeling
    }

    # In "genes" mode, also remove filtered genes from the enrichment result
    # so cnetplot doesn't render them as grey nodes
    if (!is.null(input$gseaVisFilterMode) && input$gseaVisFilterMode == "genes" && length(genes_pass_filter) > 0) {
      current_er_symbol@result$core_enrichment <- sapply(strsplit(current_er_symbol@result$core_enrichment, "/"), function(g) {
        g <- g[g %in% genes_pass_filter]
        paste(g, collapse = "/")
      })
    }

    # Network plot
    cn <- enrichplot::cnetplot(
      x = current_er_symbol,
      node_label = "none",
      foldChange = log2RatioGSEA,
      showCategory = current_er_symbol@result$Description[input$selectedTable_GSEA_rows_selected]
    )

    # Apply colour scale limits
    cn$data$foldChange[!is.na(cn$data$foldChange) & cn$data$foldChange > input$scaleLimGSEA] <- input$scaleLimGSEA
    cn$data$foldChange[!is.na(cn$data$foldChange) & cn$data$foldChange < -input$scaleLimGSEA] <- -input$scaleLimGSEA

    # Add custom styling
    cn <- cn +
      geom_point(data = cn$data, aes(x = x, y = y, size = size), shape = 21, stroke = input$nodeBorderGSEA) +
      scale_colour_gradient2(
        name = "fold change",
        low = "deepskyblue4",
        mid = "white",
        high = "firebrick",
        midpoint = 0
      )

    # Prepare data for labels
    cd1 <- cn$data
    cd1$name[!is.na(cd1$foldChange)] <- NA
    cd1$Label <- NA
    cd1$Label[seq_along(input$selectedTable_GSEA_rows_selected)] <-
      current_er_symbol@result[[labelColumn]][input$selectedTable_GSEA_rows_selected]

    cd2 <- cn$data
    cd2$name[is.na(cd2$foldChange)] <- NA
    # Apply gene label filtering if in "labels" mode
    if (!is.null(input$gseaVisFilterMode) && input$gseaVisFilterMode == "labels") {
      cd2$name[!(cd2$name %in% genes_pass_filter)] <- NA
    }

    # Add gene labels if requested
    if (input$showGeneLabelsGSEA) {
      cn <- cn + geom_text_repel(
        data = cd2,
        aes(label = name, x = x, y = y),
        size = input$textSizeGSEA / 3,
        max.overlaps = input$gseaMaxOver,
        nudge_x = input$gseaDodgeX / 10,
        nudge_y = input$gseaDodgeY / 10,
        fontface = "bold",
        color = "black",
        bg.color = "white",
        bg.r = 0.15,
        na.rm = TRUE
      )
    }

    # Add pathway labels
    cn <- cn + geom_label_repel(
      data = cd1,
      aes(label = Label, x = x, y = y, fill = Label),
      size = input$textSizeGSEA / 3,
      max.overlaps = input$gseaMaxOver,
      nudge_x = input$gseaDodgeX / 10,
      nudge_y = input$gseaDodgeY / 10,
      fontface = "bold",
      show.legend = FALSE,
      fill = alpha("white", input$gseaLabelAlpha),
      na.rm = TRUE
    )

    gseaPlotReactive$cnetplot <- cn

    # Ridge plot
    rp <- enrichplot::ridgeplot(x = erRP, fill = input$ridgePlotColourByGSEA)
    rp <- rp +
      scale_y_discrete(labels = erRP@result$Label[match(
        levels(as.factor(rp$plot_env$gs2val.df$category)),
        erRP@result$Description
      )]) +
      theme(axis.text.y = element_text(vjust = -0.01, size = input$textSizeGSEA)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_fill_gradient2(
        low = "dodgerblue4",
        high = "firebrick4",
        mid = "white",
        midpoint = 0
      )

    gseaPlotReactive$ridgeplot <- rp

    # Heatmap - use filtered genes when filtering is active
    # Prepare matrix for ComplexHeatmap
    sel_idx <- input$selectedTable_GSEA_rows_selected
    sel_er <- current_er_symbol@result[sel_idx, ]

    # Extract genes for each pathway
    pathway_genes <- lapply(strsplit(sel_er$core_enrichment, "/"), function(g) g[nchar(g) > 0])
    names(pathway_genes) <- sel_er$Description

    # Identify all unique genes across these pathways that also pass the filter
    all_genes_in_pathways <- unique(unlist(pathway_genes))
    genes_to_show <- intersect(all_genes_in_pathways, genes_pass_filter)

    if (length(genes_to_show) > 0) {
      # Create matrix
      mat <- matrix(NA, nrow = length(genes_to_show), ncol = length(sel_idx))
      rownames(mat) <- genes_to_show
      colnames(mat) <- sel_er$Description

      # Fill matrix
      for (i in seq_along(pathway_genes)) {
        p_genes <- pathway_genes[[i]]
        # Genes in this pathway that are to be shown
        common_genes <- intersect(p_genes, genes_to_show)
        if (length(common_genes) > 0) {
          # Get LFCs
          mat[common_genes, i] <- log2RatioGSEA_full[common_genes]
        }
      }

      # Create Heatmap
      limit <- input$scaleLimGSEA
      col_fun <- circlize::colorRamp2(
        c(-limit, 0, limit),
        c("deepskyblue4", "white", "firebrick")
      )

      # Custom clustering function to handle NAs (genes not in pathway)
      clust_fun <- function(m) {
        m[is.na(m)] <- 0
        hclust(dist(m))
      }

      hp <- ComplexHeatmap::Heatmap(
        mat,
        name = "Log2FC",
        col = col_fun,
        na_col = "grey95",
        cluster_rows = clust_fun,
        cluster_columns = clust_fun,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = input$textSizeGSEA),
        column_names_gp = grid::gpar(fontsize = input$textSizeGSEA),
        border = TRUE,
        rect_gp = grid::gpar(col = "black", lwd = 0.5)
      )

      gseaPlotReactive$heatmap <- hp
    } else {
      gseaPlotReactive$heatmap <- NULL
    }

    # GSEA running score plot
    rs <- enrichplot::gseaplot2(
      x = current_er,
      geneSetID = input$selectedTable_GSEA_rows_selected,
      color = RColorBrewer::brewer.pal(max(3, length(input$selectedTable_GSEA_rows_selected)), "Set2")[seq_along(input$selectedTable_GSEA_rows_selected)]
    )

    gseaPlotReactive$runningscore <- rs

    # UpSet plot
    geneSetsUpset <- setNames(
      lapply(input$selectedTable_GSEA_rows_selected, function(i) {
        colToUse <- if (input$gseaType %in% c("BP", "MF", "CC")) "geneName" else "core_enrichment"
        g <- current_er@result[[colToUse]][i] %>%
          strsplit("/") %>%
          .[[1]]
        g[nchar(g) > 1]
      }),
      current_er@result[[labelColumn]][input$selectedTable_GSEA_rows_selected]
    )

    m1 <- ComplexHeatmap::make_comb_mat(geneSetsUpset)

    up <- ComplexHeatmap::UpSet(
      m = m1,
      comb_order = order(comb_size(m1), decreasing = TRUE),
      top_annotation = upset_top_annotation(
        m1,
        add_numbers = TRUE,
        annotation_name_rot = 90,
        annotation_name_offset = unit(10, "mm")
      ),
      left_annotation = upset_left_annotation(m1, add_numbers = TRUE)
    )

    gseaPlotReactive$upset <- up
  }
)

# Render plot outputs ----

# Network plot
output$cnetPlot_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["cnetplot"]]))
    gseaPlotReactive[["cnetplot"]]
  },
  width = function() as.numeric(input$figWidthGSEA),
  height = function() as.numeric(input$figHeightGSEA)
)

output$dlPlot_cnetGSEA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_GSEA, "_", input$gseaType, "_network"),
      design_name = design,
      extension = tolower(input$downloadFormatGSEA)
    )
  },
  content = function(file) {
    req(gseaPlotReactive[["cnetplot"]])
    export_plot_standard(
      plot_object = gseaPlotReactive[["cnetplot"]],
      filename = file,
      format = tolower(input$downloadFormatGSEA),
      width = as.numeric(input$figWidthGSEA) / 50,
      height = as.numeric(input$figHeightGSEA) / 50,
      dpi = as.numeric(input$dpiGSEA)
    )
  }
)

# Heatmap
output$heatmap_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["heatmap"]]))
    ComplexHeatmap::draw(gseaPlotReactive[["heatmap"]])
  },
  width = function() as.numeric(input$figWidthGSEA),
  height = function() as.numeric(input$figHeightGSEA * 1.5)
)

output$dlPlot_heatmapGSEA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_GSEA, "_", input$gseaType, "_heatmap"),
      design_name = design,
      extension = tolower(input$downloadFormatGSEA)
    )
  },
  content = function(file) {
    req(gseaPlotReactive[["heatmap"]])
    export_plot_standard(
      plot_object = gseaPlotReactive[["heatmap"]],
      filename = file,
      format = tolower(input$downloadFormatGSEA),
      width = as.numeric(input$figWidthGSEA) / 25,
      height = as.numeric(input$figHeightGSEA) / 50,
      dpi = as.numeric(input$dpiGSEA)
    )
  }
)

# Running score
output$runningScore_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["runningscore"]]))
    gseaPlotReactive[["runningscore"]]
  },
  width = function() as.numeric(input$figWidthGSEA),
  height = function() as.numeric(input$figHeightGSEA)
)

output$dlPlot_runningScoreGSEA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_GSEA, "_", input$gseaType, "_runningscore"),
      design_name = design,
      extension = tolower(input$downloadFormatGSEA)
    )
  },
  content = function(file) {
    req(gseaPlotReactive[["runningscore"]])
    export_plot_standard(
      plot_object = gseaPlotReactive[["runningscore"]],
      filename = file,
      format = tolower(input$downloadFormatGSEA),
      width = as.numeric(input$figWidthGSEA) / 50,
      height = as.numeric(input$figHeightGSEA) / 50,
      dpi = as.numeric(input$dpiGSEA)
    )
  }
)

# Ridge plot
output$ridgePlot_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["ridgeplot"]]))
    gseaPlotReactive[["ridgeplot"]]
  },
  width = function() as.numeric(input$figWidthGSEA),
  height = function() as.numeric(input$figHeightGSEA)
)

output$dlPlot_ridgeGSEA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_GSEA, "_", input$gseaType, "_ridgeplot"),
      design_name = design,
      extension = tolower(input$downloadFormatGSEA)
    )
  },
  content = function(file) {
    req(gseaPlotReactive[["ridgeplot"]])
    export_plot_standard(
      plot_object = gseaPlotReactive[["ridgeplot"]],
      filename = file,
      format = tolower(input$downloadFormatGSEA),
      width = as.numeric(input$figWidthGSEA) / 50,
      height = as.numeric(input$figHeightGSEA) / 50,
      dpi = as.numeric(input$dpiGSEA)
    )
  }
)

# UpSet plot
output$upsetPlot_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["upset"]]))
    draw(gseaPlotReactive[["upset"]], padding = unit(c(5, 5, 5, 20), "mm"))
  },
  width = function() as.numeric(input$figWidthGSEA) * 1.3,
  height = function() as.numeric(input$figHeightGSEA)
)

output$dlPlot_upsetGSEA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_GSEA, "_", input$gseaType, "_upset"),
      design_name = design,
      extension = tolower(input$downloadFormatGSEA)
    )
  },
  content = function(file) {
    req(gseaPlotReactive[["upset"]])
    export_plot_standard(
      plot_object = gseaPlotReactive[["upset"]],
      filename = file,
      format = tolower(input$downloadFormatGSEA),
      width = as.numeric(input$figWidthGSEA) / 70,
      height = as.numeric(input$figHeightGSEA) / 90,
      dpi = as.numeric(input$dpiGSEA)
    )
  }
)
