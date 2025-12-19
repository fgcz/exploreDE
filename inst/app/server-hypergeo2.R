observeEvent(input$tabs, {
  if (input$tabs == "oraTab") {
    if (is_empty(inputDataReactive()$goObjects$ora)) {
      showModal(modalDialog(title = "No ORA Results", "No pathway analysis was performed so this tab is empty."))
    }
  }
})

oraPlotReactive <- reactiveValues(
  "er" = NULL,
  "cnetplot" = NULL,
  "barplot" = NULL,
  "dotplot" = NULL,
  "heatmap" = NULL,
  "upset" = NULL
)

if (inputDataReactive()$dataType == "RNASeq" && !is_empty(inputDataReactive()$goObjects$ora)) {
  # Get data
  ora <- inputDataReactive()$goObjects$ora
  se <- inputDataReactive()$se
  design <- inputDataReactive()$design
  param <- inputDataReactive()$param

  output$oraDesign <- renderText({
    design
  })

  # Some older datasets seem to be missing this parameter, so we will just assume it for now
  if (is.null(param$fdrThreshORA)) {
    param$fdrThreshORA <- 0.05
  }

  # Download full ORA results file
  # Download full ORA results file
  output$dlTable_ORA <- downloadHandler(
    filename = function() {
      generate_filename(
        base_name = paste0("ORA_Results"),
        design_name = design,
        extension = "xlsx"
      )
    },
    content = function(file) {
      sheets_list <- list()
      categories <- c("BP", "CC", "MF")
      gene_sets <- c("upGenes", "downGenes", "bothGenes")
      for (category in categories) {
        for (gene_set in gene_sets) {
          enrich_obj <- ora[[category]][[gene_set]]
          results_df <- as.data.frame(enrich_obj@result)
          sheet_name <- paste(category, gene_set, sep = "_")
          sheets_list[[sheet_name]] <- results_df
        }
      }
      write_xlsx(sheets_list, file)
    }
  )

  # ORA Overview Table
  ora.overview <- ezFrame(
    "Number of Significant Pathways" = numeric()
  )
  for (goDomain in names(ora)) {
    for (resName in names(ora[[goDomain]])) {
      testName <- paste0(goDomain, ": ", resName)
      ora.overview[testName, "Number of Significant Pathways"] <-
        nrow(ora[[goDomain]][[resName]]@result[
          as.numeric(
            ora[[goDomain]][[resName]]@result$p.adjust
          ) < param$fdrThreshORA,
        ])
    }
  }
  output$oraOverview <- function() {
    ora.overview %>%
      kable(
        format = "html",
        caption =
          paste("ORA GO:")
      ) %>%
      kable_styling(
        bootstrap_options = "striped",
        full_width = FALSE,
        position = "left"
      )
  }
  # ORA thresholds
  oraThresholds <- ezFrame(
    "ORA Cut-Offs" = character()
  )
  oraThresholds["Input DEG P-value", "ORA Cut-Offs"] <- param$pValThreshGO
  oraThresholds["Input DEG Log2 Ratio", "ORA Cut-Offs"] <- param$log2RatioThreshGO
  oraThresholds["Input up-regulated DEG total", "ORA Cut-Offs"] <- length(metadata(se)$enrichInput$selections$upGenes)
  oraThresholds["Input down-regulated DEG total", "ORA Cut-Offs"] <- length(metadata(se)$enrichInput$selections$downGenes)
  oraThresholds["Output ORA Term FDR", "ORA Cut-Offs"] <- param$fdrThreshORA
  output$oraThreshOverview <- function() {
    oraThresholds %>%
      kable(
        format = "html",
        caption =
          paste("Cut-Offs:")
      ) %>%
      kable_styling(
        bootstrap_options = "striped",
        full_width = FALSE,
        position = "left"
      )
  }

  # Generate main ORA table ----
  observeEvent(
    {
      input$oraType
      input$oraDirection
    },
    {
      er <- ora[[input$oraType]][[input$oraDirection]]
      # if (input$tabs == "oraTab") {
      if (is.null(er) || nrow(er@result) == 0) {
        showModal(modalDialog(title = "No ORA Results", "There are no results for this selection."))
      }
      # }
      if (!is.null(er) && nrow(er@result) != 0) {
        er@result$Label <- ""
        for (i in seq_along(er@result$Description)) {
          desc <- er@result$Description[i]
          label <- paste0(er@result$ID[i], "\n", desc)
          if (!is.na(label) && nchar(label) > 45) {
            label <- substr(label, 1, 45)
          }
          er@result$Label[i] <- label
        }
        output$selectedTable_ORA <- DT::renderDataTable({
          oraDfToShow <- as.data.frame(er)[, c("Label", "GeneRatio", "geneName")]
          oraDfToShow$geneNumber <- sapply(strsplit(as.character(oraDfToShow$GeneRatio), "/"), function(x) as.numeric(x[1]))

          DT::datatable(
            data = oraDfToShow,
            filter = "top", class = "cell-border stripe",
            rownames = FALSE, caption = "Click pathways in this table to add them to the figures on the right",
            options = list(
              columnDefs = list(
                list(visible = FALSE, targets = c(grep("geneNumber", colnames(oraDfToShow)) - 1)),
                list(orderData = c(grep("geneNumber", colnames(oraDfToShow)) - 1), targets = c(grep("GeneRatio", colnames(oraDfToShow)) - 1))
              )
            ),
          ) %>%
            DT::formatStyle(
              columns = colnames(.$x$data), `font-size` = "12px"
            )
        })
      }
    }
  )

  observeEvent(
    {
      input$textSizeORA
      input$nodeSizeORA
      input$selectedTable_ORA_rows_selected
      input$showGeneLabelsORA
      input$scaleLimORA
      input$oraMaxOver
      input$oraDodgeY
      input$oraDodgeX
      input$nodeBorderORA
      input$oraLabelAlpha
      input$oraVisPType
      input$oraVisPThreshold
      input$oraVisLFC
      input$oraVisFilterMode
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE,
    {
      req(!is.null(input$selectedTable_ORA_rows_selected))

      er <- ora[[input$oraType]][[input$oraDirection]]
      req(nrow(er@result) >= 2)

      er@result$Label <- ""
      for (i in seq_along(er@result$Description)) {
        desc <- er@result$Description[i]
        label <- paste0(er@result$ID[i], "\n", desc)
        if (!is.na(label) && nchar(label) > 45) {
          label <- substr(label, 1, 45)
        }
        er@result$Label[i] <- label
      }

      # Get log2FCs of genes for network plot
      sigGeneIds <- metadata(se)$enrichInput$selections[[input$oraDirection]]
      log2RatioORA <- metadata(se)$enrichInput$seqAnno$log2Ratio[metadata(se)$enrichInput$seqAnno$gene_id %in% sigGeneIds]
      names(log2RatioORA) <- metadata(se)$enrichInput$seqAnno$gene_name[metadata(se)$enrichInput$seqAnno$gene_id %in% sigGeneIds]
      log2RatioORA <- sort(log2RatioORA, decreasing = T)

      # Apply gene filtering based on visualisation thresholds
      log2RatioORA_full <- log2RatioORA # Keep full list for reference
      genes_pass_filter <- names(log2RatioORA) # Default: all genes pass

      cat("\n=== ORA Filtering Debug ===\n")
      cat("Filter mode:", input$oraVisFilterMode, "\n")
      cat("P-value type:", input$oraVisPType, "| Threshold:", input$oraVisPThreshold, "\n")
      cat("LFC threshold:", input$oraVisLFC, "\n")
      cat("Total genes before filter:", length(log2RatioORA), "\n")

      if (!is.null(input$oraVisFilterMode) && input$oraVisFilterMode != "none") {
        # Get p-values for genes
        seqAnno <- metadata(se)$enrichInput$seqAnno

        # Determine which p-value column to use
        pval_col <- if (!is.null(input$oraVisPType) && input$oraVisPType == "FDR") "fdr" else "pValue"

        cat("Looking for column:", pval_col, "\n")
        cat("seqAnno columns:", paste(colnames(seqAnno), collapse = ", "), "\n")

        # Validate that the column exists
        if (!is.null(seqAnno) && pval_col %in% colnames(seqAnno) && "gene_name" %in% colnames(seqAnno)) {
          # Match genes by name
          gene_match_idx <- match(names(log2RatioORA), seqAnno$gene_name)
          gene_match_idx <- gene_match_idx[!is.na(gene_match_idx)]

          cat("Matched genes:", length(gene_match_idx), "\n")

          if (length(gene_match_idx) > 0) {
            # Get p-values for matched genes
            gene_pvals <- seqAnno[[pval_col]][gene_match_idx]
            names(gene_pvals) <- seqAnno$gene_name[gene_match_idx]

            # Filter genes based on thresholds
            genes_pass_filter <- names(log2RatioORA)[
              !is.na(gene_pvals[names(log2RatioORA)]) &
                gene_pvals[names(log2RatioORA)] <= input$oraVisPThreshold &
                abs(log2RatioORA) >= input$oraVisLFC
            ]
            genes_pass_filter <- genes_pass_filter[!is.na(genes_pass_filter)]

            cat("Genes passing filter:", length(genes_pass_filter), "\n")

            # Apply filter mode
            if (input$oraVisFilterMode == "genes" && length(genes_pass_filter) > 0) {
              # Filter genes entirely - only keep genes meeting thresholds
              log2RatioORA <- log2RatioORA[names(log2RatioORA) %in% genes_pass_filter]
              cat("After filtering log2RatioORA:", length(log2RatioORA), "genes\n")
            }
          }
        } else {
          cat("Validation failed - column or seqAnno not found\n")
        }
        # For "labels" mode, keep all genes in log2RatioORA but use genes_pass_filter for labeling
      }
      cat("=== End Debug ===\n\n")

      er@result$geneID <- er@result$geneName
      nl <- "category"
      if (input$showGeneLabelsORA) {
        nl <- "all"
      }

      oraPlotReactive$er <- er

      # network plot ----
      cn <- cnetplot(
        x = er,
        foldChange = log2RatioORA,
        showCategory = er@result$Description[input$selectedTable_ORA_rows_selected],
        node_label = "none",
        cex.params = list(category_label = (input$textSizeORA / 15), gene_label = (input$textSizeORA / 18), gene_node = input$nodeSizeORA / 10, category_node = input$nodeSizeORA / 8)
      )
      # Set the colour limits to the user-specified values
      cn$data$foldChange[!is.na(cn$data$foldChange) & cn$data$foldChange > input$scaleLimORA] <- input$scaleLimORA
      cn$data$foldChange[!is.na(cn$data$foldChange) & cn$data$foldChange < -input$scaleLimORA] <- -input$scaleLimORA
      # Add our own geom_points on top
      cn <- cn + geom_point(data = cn$data, aes(x = x, y = y, size = size), shape = 21, stroke = input$nodeBorderORA)
      # Get two copies of the plot data table so we can get gene and node labels and modify them for plotting
      cd1 <- cn$data
      cd1$name[!is.na(cd1$foldChange)] <- NA
      cd1$Label <- NA
      cd1$Label[1:length(input$selectedTable_ORA_rows_selected)] <- er@result$Label[input$selectedTable_ORA_rows_selected]
      cd2 <- cn$data
      cd2$name[is.na(cd2$foldChange)] <- NA
      # Apply gene label filtering if in "labels" mode
      if (!is.null(input$oraVisFilterMode) && input$oraVisFilterMode == "labels") {
        cd2$name[!(cd2$name %in% genes_pass_filter)] <- NA
      }
      # Some reason, enrichplot changed the default such that up is blue and down is red. We prefer the opposite...
      if (input$oraDirection == "upGenes") {
        cn <- cn + scale_color_gradientn(name = "fold change", colours = c("white", "firebrick"))
      } else if (input$oraDirection == "downGenes") {
        cn <- cn + scale_color_gradientn(name = "fold change", colours = c("deepskyblue4", "white"))
      } else if (input$oraDirection == "bothGenes") {
        cn <- cn + scale_colour_gradient2(name = "fold change", low = "deepskyblue4", mid = "white", high = "firebrick", midpoint = 0)
      }
      # Add gene and node labels in that order
      if (input$showGeneLabelsORA) {
        cn <- cn +
          geom_text_repel(
            data = cd2, aes(label = name, x = x, y = y), size = input$textSizeORA / 3, max.overlaps = input$oraMaxOver, nudge_x = input$oraDodgeX / 10, nudge_y = input$oraDodgeY / 10,
            fontface = "bold", color = "black", bg.color = "white", bg.r = .15, na.rm = T
          )
      }
      cn <- cn + geom_label_repel(
        data = cd1, aes(label = Label, x = x, y = y, fill = Label), size = input$textSizeORA / 3, max.overlaps = input$oraMaxOver, nudge_x = input$oraDodgeX / 10, nudge_y = input$oraDodgeY / 10, fontface = "bold", show.legend = F, fill = alpha(c("white"), input$oraLabelAlpha), na.rm = T
      )
      oraPlotReactive$cnetplot <- cn

      # Prepare matrix for ComplexHeatmap
      sel_idx <- input$selectedTable_ORA_rows_selected
      sel_er <- er@result[sel_idx, ]

      pathway_genes <- strsplit(sel_er$geneName, "/")
      names(pathway_genes) <- sel_er$Description

      all_genes_in_pathways <- unique(unlist(pathway_genes))
      genes_to_show <- intersect(all_genes_in_pathways, genes_pass_filter)

      if (length(genes_to_show) > 0) {
        mat <- matrix(NA, nrow = length(genes_to_show), ncol = length(sel_idx))
        rownames(mat) <- genes_to_show
        colnames(mat) <- sel_er$Description

        for (i in seq_along(pathway_genes)) {
          p_genes <- pathway_genes[[i]]
          common <- intersect(p_genes, genes_to_show)
          if (length(common) > 0) {
            mat[common, i] <- log2RatioORA_full[common]
          }
        }

        limit <- input$scaleLimORA
        if (input$oraDirection == "upGenes") {
          col_fun <- circlize::colorRamp2(c(0, limit), c("white", "firebrick"))
        } else if (input$oraDirection == "downGenes") {
          col_fun <- circlize::colorRamp2(c(-limit, 0), c("deepskyblue4", "white"))
        } else {
          col_fun <- circlize::colorRamp2(c(-limit, 0, limit), c("deepskyblue4", "white", "firebrick"))
        }

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
          row_names_gp = grid::gpar(fontsize = input$textSizeORA),
          column_names_gp = grid::gpar(fontsize = input$textSizeORA),
          border = TRUE,
          rect_gp = grid::gpar(col = "black", lwd = 0.5)
        )
        oraPlotReactive$heatmap <- hp
      } else {
        oraPlotReactive$heatmap <- NULL
      }

      bp <- barplot(
        er,
        showCategory = er@result$Description[input$selectedTable_ORA_rows_selected],
        title = paste(input$oraType, input$oraDirection)
      )
      bp <- bp +
        scale_y_discrete(
          labels = rev(
            er@result$Label[
              match(bp$plot_env$df$ID, er@result$ID)
            ]
          )
        ) +
        theme_prism() +
        theme(
          axis.text.y = element_text(
            vjust = -0.01, size = input$textSizeORA
          )
        )
      oraPlotReactive$barplot <- bp

      dp <- enrichplot::dotplot(
        er,
        x = "GeneRatio",
        showCategory = er@result$Description[input$selectedTable_ORA_rows_selected],
        title = paste(input$oraType, input$oraDirection)
      )
      dp <- dp +
        scale_y_discrete(
          labels = rev(
            er@result$Label[
              match(dp$plot_env$df$ID, er@result$ID)
            ]
          )
        ) +
        theme_prism() +
        theme(
          axis.text.y = element_text(
            vjust = -0.01, size = input$textSizeORA
          )
        )
      oraPlotReactive$dotplot <- dp

      geneSetsUpset <- setNames(lapply(input$selectedTable_ORA_rows_selected, function(i) {
        g <- er@result$geneName[i] %>%
          strsplit("/") %>%
          .[[1]]
        g <- g[nchar(g) > 1]
        g
      }), er@result$Label[input$selectedTable_ORA_rows_selected])
      m1 <- ComplexHeatmap::make_comb_mat(geneSetsUpset)
      up <- ComplexHeatmap::UpSet(
        m = m1,
        comb_order = order(comb_size(m1), decreasing = T),
        top_annotation = upset_top_annotation(m1, add_numbers = TRUE, annotation_name_rot = 90, annotation_name_offset = unit(10, "mm")),
        left_annotation = upset_left_annotation(m1, add_numbers = TRUE)
      )
      oraPlotReactive$upset <- up
    }
  )
}
# Network plot
output$cnetPlot_ORA <- renderPlot(
  {
    req(!is.null(oraPlotReactive[["cnetplot"]]))
    oraPlotReactive[["cnetplot"]]
  },
  width = function() as.numeric(input$figWidthORA),
  height = function() as.numeric(input$figHeightORA)
)
output$dlPlot_cnetORA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_ORA, "_", input$oraType, "_network"),
      design_name = design,
      extension = tolower(input$downloadFormatORA)
    )
  },
  content = function(file) {
    req(oraPlotReactive[["cnetplot"]])
    export_plot_standard(
      plot_object = oraPlotReactive[["cnetplot"]],
      filename = file,
      format = tolower(input$downloadFormatORA),
      width = (as.numeric(input$figWidthORA) / 80),
      height = (as.numeric(input$figHeightORA) / 80),
      dpi = as.numeric(input$dpiORA)
    )
  }
)

# bar plot
output$barPlot_ORA <- renderPlot(
  {
    req(!is.null(oraPlotReactive[["barplot"]]))
    oraPlotReactive[["barplot"]]
  },
  width = function() as.numeric(input$figWidthORA),
  height = function() as.numeric(input$figHeightORA)
)
output$dlPlot_barORA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_ORA, "_", input$oraType, "_bar"),
      design_name = design,
      extension = tolower(input$downloadFormatORA)
    )
  },
  content = function(file) {
    req(oraPlotReactive[["barplot"]])
    export_plot_standard(
      plot_object = oraPlotReactive[["barplot"]],
      filename = file,
      format = tolower(input$downloadFormatORA),
      width = (as.numeric(input$figWidthORA) / 50),
      height = (as.numeric(input$figHeightORA) / 50),
      dpi = as.numeric(input$dpiORA)
    )
  }
)

# Dot plot
output$dotPlot_ORA <- renderPlot(
  {
    req(!is.null(oraPlotReactive[["dotplot"]]))
    oraPlotReactive[["dotplot"]]
  },
  width = function() as.numeric(input$figWidthORA),
  height = function() as.numeric(input$figHeightORA)
)
output$dlPlot_dotORA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_ORA, "_", input$oraType, "_dot"),
      design_name = design,
      extension = tolower(input$downloadFormatORA)
    )
  },
  content = function(file) {
    req(oraPlotReactive[["dotplot"]])
    export_plot_standard(
      plot_object = oraPlotReactive[["dotplot"]],
      filename = file,
      format = tolower(input$downloadFormatORA),
      width = (as.numeric(input$figWidthORA) / 50),
      height = (as.numeric(input$figHeightORA) / 50),
      dpi = as.numeric(input$dpiORA)
    )
  }
)

# Heatmap
output$heatmap_ORA <- renderPlot(
  {
    req(!is.null(oraPlotReactive[["heatmap"]]))
    ComplexHeatmap::draw(oraPlotReactive[["heatmap"]])
  },
  width = function() as.numeric(input$figWidthORA / 2),
  height = function() as.numeric(input$figHeightORA * 1.5)
)
output$dlPlot_heatmapORA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_ORA, "_", input$oraType, "_heatmap"),
      design_name = design,
      extension = tolower(input$downloadFormatORA)
    )
  },
  content = function(file) {
    req(oraPlotReactive[["heatmap"]])
    export_plot_standard(
      plot_object = oraPlotReactive[["heatmap"]],
      filename = file,
      format = tolower(input$downloadFormatORA),
      width = (as.numeric(input$figWidthORA) / 25),
      height = (as.numeric(input$figHeightORA) / 50),
      dpi = as.numeric(input$dpiORA)
    )
  }
)

# upset
output$upset_ORA <- renderPlot(
  {
    req(!is.null(oraPlotReactive[["upset"]]))
    draw(oraPlotReactive[["upset"]], padding = unit(c(5, 5, 5, 20), "mm"))
  },
  width = function() as.numeric(input$figWidthORA) * 1.3,
  height = function() as.numeric(input$figHeightORA)
)
output$dlPlot_upsetORA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_ORA, "_", input$oraType, "_upset"),
      design_name = design,
      extension = tolower(input$downloadFormatORA)
    )
  },
  content = function(file) {
    req(oraPlotReactive[["upset"]])
    export_plot_standard(
      plot_object = oraPlotReactive[["upset"]],
      filename = file,
      format = tolower(input$downloadFormatORA),
      width = (as.numeric(input$figWidthORA) / 70),
      height = (as.numeric(input$figHeightORA) / 90),
      dpi = as.numeric(input$dpiORA)
    )
  }
)
