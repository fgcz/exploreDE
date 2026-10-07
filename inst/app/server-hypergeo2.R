eds <- inputDataReactive()$.eds
req(eds)
# Enrichment reads go through the backend-agnostic accessors so an uploaded
# .h5ad surfaces ORA like a SUSHI SE. eds$enrichment() returns a FLAT list keyed
# BP/upGenes... (SE) or ora/BP/upGenes... (h5ad); strip the prefix and rebuild
# the ontology x direction nesting this module indexes as ora[[type]][[dir]].
contrast <- if (length(eds$contrasts()) > 0) eds$primary_contrast() else NULL

# A cell in the ora[[type]][[dir]] structure below is EITHER a real
# clusterProfiler S4 object (the native path -- existing behaviour, unchanged)
# or a plain tidy data.frame (a BYO dataset whose pathway data was never a
# clusterProfiler object at all -- see plans/exploredelight/enrichment_detach_plan.md). These
# two helpers let the rest of this module read either uniformly; only the
# genuinely S4-only plots (cnetplot/barplot/dotplot, which dispatch on the
# class) need to branch on `.is_native()`.
.is_native <- function(x) inherits(x, c("enrichResult", "gseaResult"))
.res_df    <- function(x) if (.is_native(x)) x@result else x

.ora_raw <- function() {
  if (is.null(contrast)) return(NULL)
  nested <- list()
  flat <- eds$enrichment(contrast, "ora")
  if (length(flat)) {
    names(flat) <- sub("^ora/", "", names(flat))
    for (k in names(flat)) {
      parts <- strsplit(k, "/", fixed = TRUE)[[1]]
      if (length(parts) >= 2) { dom <- parts[1]; dir <- parts[length(parts)] }
      else { dom <- "GO"; dir <- parts[1] }
      nested[[dom]][[dir]] <- flat[[k]]
    }
  }
  # Fill any gaps from the tidy accessor, split by its `ontology` column the
  # same way -- this is the ONLY source of data for a dataset whose pathway
  # entries were never clusterProfiler objects (plain data.frames: nothing to
  # rebuild objects from). Objects always win where both exist (cnetplot below
  # needs one), so this only ADDS cells, never overwrites one already set from
  # `flat` above.
  tab <- eds$enrichment_table(contrast, "ora")
  if (is.data.frame(tab) && nrow(tab) && "ontology" %in% colnames(tab)) {
    for (k in unique(tab$ontology)) {
      parts <- strsplit(k, "/", fixed = TRUE)[[1]]
      if (length(parts) >= 2) { dom <- parts[1]; dir <- parts[length(parts)] }
      else { dom <- "GO"; dir <- parts[1] }
      if (is.null(nested[[dom]][[dir]]))
        nested[[dom]][[dir]] <- tab[tab$ontology == k, , drop = FALSE]
    }
  }
  if (!length(nested)) return(NULL)
  nested
}

observeEvent(input$tabs, {
  if (input$tabs == "oraTab") {
    if (is_empty(.ora_raw())) {
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
  "upset" = NULL,
  "lollipop" = NULL
)

# Contrast-dependent state. (Re)built by the observeEvent(selectedContrast())
# block below (higher priority than the oraType observer that reads it) and read
# by the observers/outputs. Single-contrast data fires it once.
ora     <- NULL
ei      <- list(seq_anno = NULL, ranked_l2fc = NULL, selections = list())
param   <- list()
design  <- ""
sel     <- list()
.id2sym <- function(ids) ids

if (eds$data_type() == "RNASeq" && !is_empty(.ora_raw())) {
 observeEvent(selectedContrast(), {
  # Get data through the accessor layer (backend-agnostic).
  contrast <<- selectedContrast()
  ora    <<- .ora_raw()
  ei     <<- eds$enrichment_input(contrast)
  design <<- if (!is.null(contrast)) eds$contrast_label(contrast) else ""
  # ORA thresholds / featureLevel: SUSHI SE metadata$param, else defaults.
  p <- if (eds$backend_kind() == "se") S4Vectors::metadata(eds$backend())$param else list()
  if (is.null(p)) p <- list()
  if (is.null(p$featureLevel)) p$featureLevel <- eds$feature_id_kind()
  if (is.null(p$fdrThreshORA)) p$fdrThreshORA <- 0.05
  param <<- p

  # Significant-gene selections, with bothGenes synthesised when the backend
  # only provides up/down (uploaded h5ad).
  s <- ei$selections %||% list()
  if (is.null(s$bothGenes)) s$bothGenes <- union(s$upGenes, s$downGenes)
  sel <<- s

  # Backend-agnostic gene-id -> symbol remap (no-op when feature_ids are already
  # symbols, e.g. an uploaded h5ad), so labels never blank out.
  .id2sym <<- local({
    sa <- ei$seq_anno
    gid <- gsym <- NULL
    if (!is.null(sa)) {
      gid  <- if ("gene_id" %in% colnames(sa)) sa$gene_id else if ("label_gene_id" %in% colnames(sa)) sa$label_gene_id else NULL
      gsym <- if ("gene_name" %in% colnames(sa)) sa$gene_name else if ("label_gene_symbol" %in% colnames(sa)) sa$label_gene_symbol else NULL
    }
    if (is.null(gid) || is.null(gsym)) {
      function(ids) ids
    } else {
      map <- setNames(as.character(gsym), as.character(gid))
      function(ids) { s <- unname(map[as.character(ids)]); ifelse(is.na(s) | !nzchar(s), ids, s) }
    }
  })

  output$oraDesign <- renderText({
    design
  })

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
          if (is.null(enrich_obj)) next
          sheet_name <- paste(category, gene_set, sep = "_")
          sheets_list[[sheet_name]] <- .res_df(enrich_obj)
        }
      }
      write_xlsx(sheets_list, file)
    }
  )

  # ORA Overview Table
  ora.overview <- data.frame(check.names = FALSE,
    "Number of Significant Pathways" = numeric()
  )
  for (goDomain in names(ora)) {
    for (resName in names(ora[[goDomain]])) {
      testName <- paste0(goDomain, ": ", resName)
      rdf <- .res_df(ora[[goDomain]][[resName]])
      ora.overview[testName, "Number of Significant Pathways"] <-
        nrow(rdf[as.numeric(rdf$p.adjust) < param$fdrThreshORA, ])
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
  oraThresholds <- data.frame(check.names = FALSE,
    "ORA Cut-Offs" = character()
  )
  oraThresholds["Input DE P-value", "ORA Cut-Offs"] <- param$pValThreshGO %||% NA
  oraThresholds["Input DE Log2 Ratio", "ORA Cut-Offs"] <- param$log2RatioThreshGO %||% NA
  oraThresholds["Input up-regulated DE total", "ORA Cut-Offs"] <- length(sel$upGenes)
  oraThresholds["Input down-regulated DE total", "ORA Cut-Offs"] <- length(sel$downGenes)
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

  # Isoform-level aggregation UI (only shown for isoform data)
  if (!is.null(param$featureLevel) && param$featureLevel == "isoform") {
    output$oraIsoformAggUI <- renderUI({
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
            inputId = "oraIsoformAgg",
            label = "Aggregation method",
            choices = c("Max absolute" = "maxabs", "Mean" = "mean", "Median" = "median", "Max" = "max", "Min" = "min"),
            selected = "maxabs",
            width = "50%"
          )
        )
      )
    })
    outputOptions(output, "oraIsoformAggUI", suspendWhenHidden = FALSE)
  }
 }, ignoreNULL = FALSE, priority = 10)

  # Generate main ORA table (also refires on contrast switch, after the
  # higher-priority setup block above refreshes ora/ei/sel/.id2sym).
  observeEvent(
    {
      input$oraType
      input$oraDirection
      selectedContrast()
    },
    {
      # Empty any stale figures and clear the row selection so a contrast /
      # ontology / direction switch resets the plots (mirrors the GSEA tab),
      # rather than leaving the previous pathway's network plot on screen until
      # a new row is clicked. The plot observer below only refires on a fresh
      # row selection, so without this the old plots linger.
      for (x in names(oraPlotReactive)) {
        oraPlotReactive[[x]] <- NULL
      }
      shinyjs::delay(100, {
        proxy <- DT::dataTableProxy("selectedTable_ORA")
        DT::selectRows(proxy, NULL)
      })

      er  <- ora[[input$oraType]][[input$oraDirection]]
      rdf <- .res_df(er)
      # if (input$tabs == "oraTab") {
      if (is.null(er) || nrow(rdf) == 0) {
        showModal(modalDialog(title = "No ORA Results", "There are no results for this selection."))
      }
      # }
      if (!is.null(er) && nrow(rdf) != 0) {
        rdf$Label <- ""
        for (i in seq_along(rdf$Description)) {
          desc <- rdf$Description[i]
          label <- paste0(rdf$ID[i], "\n", desc)
          if (!is.na(label) && nchar(label) > 45) {
            label <- substr(label, 1, 45)
          }
          rdf$Label[i] <- label
        }
        # Rebuild geneName from geneID via the backend-agnostic remap.
        rdf$geneName <- sapply(strsplit(rdf$geneID, "/"), function(ids) {
          paste(.id2sym(ids), collapse = "/")
        })
        output$selectedTable_ORA <- DT::renderDataTable({
          oraDfToShow <- rdf[, c("Label", "GeneRatio", "p.adjust", "geneName")]
          colnames(oraDfToShow)[colnames(oraDfToShow) == "p.adjust"] <- "FDR"
          oraDfToShow$geneNumber <- sapply(strsplit(as.character(oraDfToShow$GeneRatio), "/"), function(x) as.numeric(x[1]))

          DT::datatable(
            data = oraDfToShow,
            filter = "top", class = "cell-border stripe",
            rownames = FALSE,
            caption = paste0(
              "Click pathways in this table to add them to the figures on the right. ",
              "Only pathways with FDR < ", param$fdrThreshORA, " (green) are plotted; ",
              "clicking any other pathway will not add it to the figures."
            ),
            options = list(
              columnDefs = list(
                list(visible = FALSE, targets = c(grep("geneNumber", colnames(oraDfToShow)) - 1)),
                list(orderData = c(grep("geneNumber", colnames(oraDfToShow)) - 1), targets = c(grep("GeneRatio", colnames(oraDfToShow)) - 1))
              )
            ),
          ) %>%
            DT::formatStyle(
              columns = colnames(.$x$data), `font-size` = "12px"
            ) %>%
            DT::formatSignif(columns = "FDR", digits = 3) %>%
            DT::formatStyle(
              columns = "FDR",
              color = styleInterval(cuts = param$fdrThreshORA, values = c("green", "black")),
              fontWeight = "bold"
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
      input$oraIsoformAgg
      input$oraLollipopScale
      input$oraLollipopColour
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE,
    {
      req(!is.null(input$selectedTable_ORA_rows_selected))

      er  <- ora[[input$oraType]][[input$oraDirection]]
      req(!is.null(er))
      rdf <- .res_df(er)
      req(nrow(rdf) >= 2)

      rdf$Label <- ""
      for (i in seq_along(rdf$Description)) {
        desc <- rdf$Description[i]
        label <- paste0(rdf$ID[i], "\n", desc)
        if (!is.na(label) && nchar(label) > 45) {
          label <- substr(label, 1, 45)
        }
        rdf$Label[i] <- label
      }

      # Only pathways passing the term FDR cut-off (green in the table) are
      # plotted -- enrichplot already silently drops the rest from the native
      # plots, so apply the same rule to every figure and say so.
      sel_rows <- input$selectedTable_ORA_rows_selected
      sel_rows <- sel_rows[!is.na(rdf$p.adjust[sel_rows]) & rdf$p.adjust[sel_rows] <= param$fdrThreshORA]
      n_dropped <- length(input$selectedTable_ORA_rows_selected) - length(sel_rows)
      if (n_dropped > 0) {
        showNotification(
          paste0(n_dropped, " selected pathway(s) with FDR > ", param$fdrThreshORA,
                 " not added to the figures."),
          id = "oraFdrDrop", type = "warning"
        )
      }
      if (!length(sel_rows)) {
        for (x in names(oraPlotReactive)) oraPlotReactive[[x]] <- NULL
        return()
      }

      # Get log2FCs of genes for network plot
      sigGeneIds <- sel[[input$oraDirection]]
      log2RatioORA <- ei$ranked_l2fc[names(ei$ranked_l2fc) %in% sigGeneIds]
      gene_ids_ora <- names(log2RatioORA)
      # Relabel to gene symbols where possible (no-op when already symbols).
      names(log2RatioORA) <- .id2sym(gene_ids_ora)
      # For isoform data, aggregate multiple transcripts per gene
      if (!is.null(param$featureLevel) && param$featureLevel == "isoform") {
        agg_method <- if (!is.null(input$oraIsoformAgg)) input$oraIsoformAgg else "maxabs"
        na_safe <- function(fn) function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else fn(x) }
        agg_fn <- switch(agg_method,
          "maxabs" = na_safe(function(x) x[which.max(abs(x))]),
          "mean"   = na_safe(mean),
          "median" = na_safe(median),
          "max"    = na_safe(max),
          "min"    = na_safe(min)
        )
        log2RatioORA <- vapply(split(log2RatioORA, names(log2RatioORA)), agg_fn, numeric(1))
      }
      log2RatioORA <- sort(log2RatioORA, decreasing = T)

      # Apply gene filtering based on visualisation thresholds
      log2RatioORA_full <- log2RatioORA # Keep full list for reference
      genes_pass_filter <- names(log2RatioORA) # Default: all genes pass

      if (!is.null(input$oraVisFilterMode) && input$oraVisFilterMode != "none") {
        # Per-gene p-values (from the enrichment input's annotation; absent for
        # uploads -> filter degrades off).
        seqAnno <- ei$seq_anno

        # Determine which p-value column to use
        pval_col <- if (!is.null(input$oraVisPType) && input$oraVisPType == "FDR") "fdr" else "pValue"

        # Validate that the column exists
        if (!is.null(seqAnno) && pval_col %in% colnames(seqAnno) && "gene_name" %in% colnames(seqAnno)) {
          # Match genes by name
          gene_match_idx <- match(names(log2RatioORA), seqAnno$gene_name)
          gene_match_idx <- gene_match_idx[!is.na(gene_match_idx)]

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

            # Apply filter mode
            if (input$oraVisFilterMode == "genes") {
              if (length(genes_pass_filter) == 0) {
                showNotification("No genes meet the current filter thresholds. Lower the Log2FC or P-Value threshold.", type = "warning")
                return()
              }
              log2RatioORA <- log2RatioORA[names(log2RatioORA) %in% genes_pass_filter]
            }
          }
        }
        # For "labels" mode, keep all genes in log2RatioORA but use genes_pass_filter for labeling
      }

      # Rebuild geneName from geneID via the backend-agnostic remap.
      rdf$geneName <- sapply(strsplit(rdf$geneID, "/"), function(ids) {
        paste(.id2sym(ids), collapse = "/")
      })

      # In "genes" mode, also remove filtered genes from the enrichment result
      # so cnetplot doesn't render them as grey nodes
      if (!is.null(input$oraVisFilterMode) && input$oraVisFilterMode == "genes") {
        rdf$geneName <- sapply(strsplit(rdf$geneName, "/"), function(g) {
          g <- g[g %in% genes_pass_filter]
          paste(g, collapse = "/")
        })
        # Check that selected pathways still have genes after filtering
        sel_genes <- rdf$geneName[sel_rows]
        if (!any(nchar(sel_genes) > 0)) {
          showNotification("No genes in the selected pathways meet the current filter thresholds. Lower the Log2FC or P-Value threshold.", type = "warning")
          return()
        }
      }

      rdf$geneID <- rdf$geneName
      nl <- "category"
      if (input$showGeneLabelsORA) {
        nl <- "all"
      }

      # cnetplot needs the real clusterProfiler object (its geneSets slot, not
      # just the result table -- see plans/exploredelight/enrichment_detach_plan.md §9), so it
      # only runs on the native path. Sync the augmented table back into the
      # S4 object's result slot first so cnetplot sees the same Label/geneName
      # this reactive just built -- this is otherwise a no-op for the native
      # path (identical to mutating er@result directly, just expressed via rdf
      # so the table/heatmap/upset/lollipop code below can stay one code path
      # for both the native and tidy-only cases).
      if (.is_native(er)) {
        er@result <- rdf
        oraPlotReactive$er <- er

        # network plot ----
        cn <- cnetplot(
          x = er,
          foldChange = log2RatioORA,
          showCategory = er@result$Description[sel_rows],
          node_label = "none"
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
        cd1$Label[1:length(sel_rows)] <- er@result$Label[sel_rows]
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
      } else {
        oraPlotReactive$er <- rdf
        oraPlotReactive$cnetplot <- NULL
        showNotification("Network plot needs the original enrichment object; not available for this dataset.", type = "message")
      }

      # Prepare matrix for ComplexHeatmap
      sel_idx <- sel_rows
      sel_er <- rdf[sel_idx, ]

      pathway_genes <- lapply(strsplit(sel_er$geneName, "/"), function(g) g[nchar(g) > 0])
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

      if (.is_native(er)) {
        bp <- barplot(
          er,
          showCategory = er@result$Description[sel_rows],
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
          showCategory = er@result$Description[sel_rows],
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
      } else {
        # barplot()/dotplot() dispatch on the enrichResult class -- with no
        # native S4 object (tidy-only pathway data), rebuild the same view by
        # hand from the result table alone (both only ever needed `result`;
        # see plans/exploredelight/enrichment_detach_plan.md).
        sel_df <- rdf[sel_rows, ]
        sel_df$Label <- factor(sel_df$Label, levels = rev(sel_df$Label[order(sel_df$Count)]))

        bp <- ggplot(sel_df, aes(x = Count, y = Label, fill = p.adjust)) +
          geom_col() +
          scale_fill_continuous(low = "red", high = "blue", name = "p.adjust") +
          labs(x = "Count", y = NULL, title = paste(input$oraType, input$oraDirection)) +
          theme_prism() +
          theme(axis.text.y = element_text(vjust = -0.01, size = input$textSizeORA))
        oraPlotReactive$barplot <- bp

        sel_df$GeneRatioNum <- sapply(strsplit(as.character(sel_df$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
        dp <- ggplot(sel_df, aes(x = GeneRatioNum, y = Label, size = Count, colour = p.adjust)) +
          geom_point() +
          scale_colour_continuous(low = "red", high = "blue", name = "p.adjust") +
          labs(x = "GeneRatio", y = NULL, title = paste(input$oraType, input$oraDirection)) +
          theme_prism() +
          theme(axis.text.y = element_text(vjust = -0.01, size = input$textSizeORA))
        oraPlotReactive$dotplot <- dp
      }

      lolli_df <- rdf[sel_rows, , drop = FALSE]
      lolli_df$negLogP <- -log10(lolli_df$p.adjust)
      lolli_df$Label <- factor(lolli_df$Label, levels = lolli_df$Label[order(lolli_df$negLogP)])
      scl <- input$oraLollipopScale
      ts <- input$textSizeORA
      lp <- ggplot(lolli_df, aes(y = Label, x = negLogP)) +
        geom_segment(aes(x = 0, xend = negLogP, y = Label, yend = Label), linewidth = scl) +
        geom_point(aes(size = Count), colour = input$oraLollipopColour) +
        scale_size(range = c(ts / 3, ts) * scl, name = "Gene Count") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(x = "-log10(FDR)", y = NULL, title = paste(input$oraType, input$oraDirection)) +
        theme_prism(base_size = ts)
      oraPlotReactive$lollipop <- lp

      geneSetsUpset <- setNames(lapply(sel_rows, function(i) {
        g <- rdf$geneName[i] %>%
          strsplit("/") %>%
          .[[1]]
        g <- unique(g[nchar(g) > 1])
        g
      }), rdf$Label[sel_rows])
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
    draw_heatmap(oraPlotReactive[["heatmap"]],
                 legend_dx = input$ORALegendDX, legend_dy = input$ORALegendDY)
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
      plot_object = function() draw_heatmap(oraPlotReactive[["heatmap"]],
                                            legend_dx = input$ORALegendDX, legend_dy = input$ORALegendDY),
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

# Lollipop plot
output$lollipopPlot_ORA <- renderPlot(
  {
    req(!is.null(oraPlotReactive[["lollipop"]]))
    oraPlotReactive[["lollipop"]]
  },
  width = function() as.numeric(input$figWidthORA),
  height = function() as.numeric(input$figHeightORA)
)
output$dlPlot_lollipopORA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_ORA, "_", input$oraType, "_lollipop"),
      design_name = design,
      extension = tolower(input$downloadFormatORA)
    )
  },
  content = function(file) {
    req(oraPlotReactive[["lollipop"]])
    export_plot_standard(
      plot_object = oraPlotReactive[["lollipop"]],
      filename = file,
      format = tolower(input$downloadFormatORA),
      width = (as.numeric(input$figWidthORA) / 50),
      height = (as.numeric(input$figHeightORA) / 50),
      dpi = as.numeric(input$dpiORA)
    )
  }
)
