eds <- inputDataReactive()$.eds
req(eds)
# Enrichment reads go through the backend-agnostic accessors so an uploaded
# .h5ad surfaces GSEA just like a SUSHI SE. eds$enrichment() returns a FLAT list
# whose keys are BP/MF/CC for the SE backend and gsea/BP... for h5ad — strip the
# prefix so both normalise to the BP/MF/CC the rest of this module indexes by.
contrast <- if (length(eds$contrasts()) > 0) eds$primary_contrast() else NULL

# See server-hypergeo2.R for the identical rationale: a cell in `gsea` below is
# EITHER a real clusterProfiler S4 object (native path, unchanged) or a plain
# tidy data.frame (a BYO dataset with no clusterProfiler object at all).
.is_native <- function(x) inherits(x, c("enrichResult", "gseaResult"))
.res_df <- function(x) if (.is_native(x)) x@result else x

# FDR threshold for "ora" / "gsea": the SUSHI parameter (fdrThreshORA /
# fdrThreshGSEA) for an SE; for an h5ad, the per-method value that
# build_explore_h5ad() recorded in uns$enrichment_filter; otherwise 0.05.
.enr_fdr_thresh <- function(kind) {
  p <- if (eds$backend_kind() == "se") S4Vectors::metadata(eds$backend())$param else NULL
  v <- switch(kind, ora = p$fdrThreshORA, gsea = p$fdrThreshGSEA)
  if (is.null(v)) {
    f <- eds$enrichment_filter()
    if (is.list(f)) v <- f[[kind]]
  }
  if (is.null(v) || is.na(v)) 0.05 else as.numeric(v)
}
# Keep only the terms with p.adjust <= thr. ezRun stores every tested term (the
# full table goes to its Excel export); the app, like the HTML report, only
# ever shows the significant ones.
.enr_keep_sig <- function(x, thr) {
  df <- .res_df(x)
  if (!is.data.frame(df) || !nrow(df) || !"p.adjust" %in% colnames(df)) {
    return(x)
  }
  keep <- !is.na(df$p.adjust) & df$p.adjust <= thr
  if (.is_native(x)) {
    x@result <- df[keep, , drop = FALSE]
    x
  } else {
    df[keep, , drop = FALSE]
  }
}

.gsea_raw <- function() {
  if (is.null(contrast)) {
    return(NULL)
  }
  g <- list()
  flat <- eds$enrichment(contrast, "gsea")
  if (length(flat)) {
    names(flat) <- sub("^gsea/", "", names(flat))
    g <- flat
  }
  # Fill any gaps from the tidy accessor -- the only source at all for a
  # dataset whose GSEA entries were never clusterProfiler objects -- keyed by
  # its `ontology` column the same way enrichment() names its leaves. Native
  # always wins where both exist (it alone carries geneList/geneSets, needed
  # for cnetplot/ridgeplot/gseaplot2 below).
  tab <- eds$enrichment_table(contrast, "gsea")
  if (is.data.frame(tab) && nrow(tab) && "ontology" %in% colnames(tab)) {
    for (k in unique(tab$ontology)) {
      if (is.null(g[[k]])) g[[k]] <- tab[tab$ontology == k, , drop = FALSE]
    }
  }
  g <- lapply(g, .enr_keep_sig, thr = .enr_fdr_thresh("gsea"))
  if (!length(g)) {
    return(NULL)
  }
  g
}

# Check if GSEA tab is accessed without results
observeEvent(input$tabs, {
  if (input$tabs == "gseaTab") {
    g <- .gsea_raw()
    if (is.null(g) || length(g) == 0) {
      showModal(modalDialog(
        title = "No GSEA Results",
        "No pathway analysis was performed for this dataset, so this tab is empty.",
        footer = modalButton("OK")
      ))
      return()
    }

    # Check if current selection has no results
    if (!is.null(input$gseaType) && !is.null(g[[input$gseaType]])) {
      er <- g[[input$gseaType]]
      if (is.null(er) || nrow(.res_df(er)) == 0) {
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
  upset = NULL,
  lollipop = NULL
)
gseaToUse <- reactiveValues(er = NULL)

# Contrast-dependent state. (Re)built by the observeEvent(selectedContrast())
# block below and read by the observers/outputs. Kept as plain module vars so the
# deep observers see the latest without being re-registered (no stacking). For a
# single-contrast dataset selectedContrast() is constant, so this fires once.
gsea <- NULL
ei <- list(seq_anno = NULL, ranked_l2fc = NULL, selections = list())
param <- list()
design <- ""
.id2sym <- function(ids) ids

# Higher priority than the gseaType observer below so the module vars are fresh
# before that observer reads them on a contrast switch.
observeEvent(
  selectedContrast(),
  {
    contrast <<- selectedContrast()
    gsea <<- .gsea_raw()
    ei <<- if (!is.null(contrast)) {
      eds$enrichment_input(contrast)
    } else {
      list(seq_anno = NULL, ranked_l2fc = NULL, selections = list())
    }
    design <<- if (!is.null(contrast)) eds$contrast_label(contrast) else ""

    # GSEA thresholds / featureLevel: SUSHI SE metadata$param, else defaults.
    p <- if (eds$backend_kind() == "se") S4Vectors::metadata(eds$backend())$param else list()
    if (is.null(p)) {
      p <- list()
    }
    if (is.null(p$featureLevel)) {
      p$featureLevel <- eds$feature_id_kind()
    }
    p$fdrThreshGSEA <- .enr_fdr_thresh("gsea")
    param <<- p

    # Backend-agnostic gene-id -> symbol remap (no-op when feature_ids are already
    # symbols, e.g. an uploaded h5ad), so labels never blank out.
    .id2sym <<- local({
      sa <- ei$seq_anno
      gid <- gsym <- NULL
      if (!is.null(sa)) {
        gid <- if ("gene_id" %in% colnames(sa)) {
          sa$gene_id
        } else if ("label_gene_id" %in% colnames(sa)) {
          sa$label_gene_id
        } else {
          NULL
        }
        gsym <- if ("gene_name" %in% colnames(sa)) {
          sa$gene_name
        } else if ("label_gene_symbol" %in% colnames(sa)) {
          sa$label_gene_symbol
        } else {
          NULL
        }
      }
      if (is.null(gid) || is.null(gsym)) {
        function(ids) ids
      } else {
        map <- setNames(as.character(gsym), as.character(gid))
        function(ids) {
          s <- unname(map[as.character(ids)])
          ifelse(is.na(s) | !nzchar(s), ids, s)
        }
      }
    })

    # Update GSEA type choices for this contrast
    updateSelectInput(
      session = session,
      inputId = "gseaType",
      choices = names(gsea),
      selected = if (length(gsea)) names(gsea)[[1]] else NULL
    )

    output$gseaDesign <- renderText({
      design
    })

    # GSEA overview table
    gsea.overview <- data.frame(check.names = FALSE, "Number of Significant Pathways" = numeric())
    for (i in seq_along(names(gsea))) {
      testName <- names(gsea)[i]
      gsea.overview[testName, "Number of Significant Pathways"] <- nrow(.res_df(gsea[[i]]))
    }
    output$gseaOverview <- function() {
      gsea.overview %>%
        kable(format = "html", caption = "GSEA:") %>%
        kable_styling(bootstrap_options = "striped", full_width = FALSE, position = "left")
    }

    # GSEA thresholds table
    gseaThresholds <- data.frame(check.names = FALSE, "GSEA Cut-Offs" = character())
    gseaThresholds["Output GSEA Term FDR", "GSEA Cut-Offs"] <- param$fdrThreshGSEA
    output$gseaThreshOverview <- function() {
      gseaThresholds %>%
        kable(format = "html", caption = "Cut-Offs:") %>%
        kable_styling(bootstrap_options = "striped", full_width = FALSE, position = "left")
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
                tags$b("Note:"),
                "This dataset is at isoform level. Multiple transcripts per gene
              were passed to the enrichment analysis, each with their own log2 fold change. To
              display a single fold change per gene in the network plot, values are aggregated
              using the method selected below."
              )
            ),
            selectInput(
              inputId = "gseaIsoformAgg",
              label = "Aggregation method",
              choices = c(
                "Max absolute" = "maxabs",
                "Mean" = "mean",
                "Median" = "median",
                "Max" = "max",
                "Min" = "min"
              ),
              selected = "maxabs",
              width = "50%"
            )
          )
        )
      })
      outputOptions(output, "gseaIsoformAggUI", suspendWhenHidden = FALSE)
    }
  },
  ignoreNULL = FALSE,
  priority = 10
)

# Main observer for GSEA type changes (also refires on contrast switch, after the
# higher-priority setup block above refreshes gsea/ei/.id2sym).
observeEvent(
  {
    input$gseaType
    selectedContrast()
  },
  {
    gseaToUse$er <- NULL
    gseaToUse$er_symbol <- NULL
    gseaToUse$rdf <- NULL
    gseaToUse$rdf_symbol <- NULL
    gseaToUse$is_native <- FALSE
    er <- gsea[[input$gseaType]]
    rdf <- .res_df(er)

    # Check if results exist and show modal if empty
    if (input$tabs == "gseaTab" && (is.null(er) || nrow(rdf) == 0)) {
      showModal(modalDialog(
        title = "No GSEA Results",
        "There are no results for this selection."
      ))
      return()
    }

    # Process results if they exist
    if (!is.null(er) && nrow(rdf) > 0) {
      # Handle GO terms (BP, MF, CC) differently
      if (input$gseaType %in% c("BP", "MF", "CC")) {
        rdf$Label <- ""

        # Create labels for GO terms
        for (i in seq_along(rdf$Description)) {
          desc <- rdf$Description[i]
          label <- paste0(rdf$ID[i], "\n", desc)
          if (!is.na(label) && nchar(label) > 45) {
            label <- substr(label, 1, 45)
          }
          rdf$Label[i] <- label
        }
        colsToTable <- c("Label", "enrichmentScore", "NES", "p.adjust", "geneName")
        # Rebuild geneName from core_enrichment via the backend-agnostic remap
        # (falls back to the raw ID where no symbol is known).
        rdf$geneName <- sapply(strsplit(rdf$core_enrichment, "/"), function(ids) {
          paste(.id2sym(ids), collapse = "/")
        })
        rdf_symbol <- rdf
        # Convert core enrichment from gene IDs to gene names
        rdf_symbol$core_enrichment <- vapply(
          rdf_symbol$core_enrichment,
          function(genes_string) {
            gene_ids <- strsplit(genes_string, "/")[[1]]
            paste(.id2sym(gene_ids), collapse = "/")
          },
          character(1),
          USE.NAMES = FALSE
        )
      } else {
        # For other pathway databases (e.g., C2)
        rdf_symbol <- rdf
        colsToTable <- c("Description", "enrichmentScore", "NES", "p.adjust", "core_enrichment")
      }
      # Prepare log2 ratio data for plotting (store raw, aggregate in plot observer)
      log2RatioGSEA <- ei$ranked_l2fc
      gene_ids_gsea <- names(log2RatioGSEA)
      # Relabel to gene symbols where possible (no-op when already symbols).
      names(log2RatioGSEA) <- .id2sym(gene_ids_gsea)
      log2RatioGSEA <- sort(log2RatioGSEA, decreasing = TRUE)

      # Update reactive values. gseaplot2/ridgeplot/cnetplot/slice need the
      # real S4 object (its geneList/geneSets, not just the result table --
      # see plans/exploredelight/enrichment_detach_plan.md §9), so er/er_symbol stay NULL for
      # tidy-only pathway data; everything else reads rdf/rdf_symbol either way.
      if (.is_native(er)) {
        er@result <- rdf
        er_symbol <- er
        er_symbol@result <- rdf_symbol
        gseaToUse$er <- er
        gseaToUse$er_symbol <- er_symbol
        gseaToUse$is_native <- TRUE
      }
      gseaToUse$rdf <- rdf
      gseaToUse$rdf_symbol <- rdf_symbol
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
  req(gseaToUse$rdf, gseaToUse$colsToTable)

  # Show FDR (p.adjust, which the enrichment was thresholded on) in place of qvalue.
  gseaDfToShow <- gseaToUse$rdf[, gseaToUse$colsToTable]
  colnames(gseaDfToShow)[colnames(gseaDfToShow) == "p.adjust"] <- "FDR"

  DT::datatable(
    data = gseaDfToShow,
    filter = "top",
    class = "cell-border stripe",
    rownames = FALSE,
    caption = paste0(
      "Click pathways in this table to add them to the figures on the right. ",
      "Only pathways with FDR <= ",
      param$fdrThreshGSEA,
      " are listed."
    )
  ) %>%
    DT::formatStyle(
      columns = colnames(.$x$data),
      `font-size` = "12px"
    ) %>%
    DT::formatSignif(
      columns = c("enrichmentScore", "NES", "FDR"),
      digits = 3
    ) %>%
    DT::formatStyle(
      columns = c("enrichmentScore", "NES"),
      color = styleInterval(cuts = 0, values = c("blue", "darkorange")),
      fontWeight = "bold"
    ) %>%
    DT::formatStyle(
      columns = "FDR",
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
        sheet_name <- paste(category, sep = "_")
        sheets_list[[sheet_name]] <- .res_df(gsea[[category]])
      }
    }
    # Handle custom GSEA types if present (e.g., C2)
    other_types <- setdiff(names(gsea), categories)
    for (type in other_types) {
      sheets_list[[type]] <- .res_df(gsea[[type]])
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
    input$gseaLollipopScale
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE,
  {
    req(gseaToUse$rdf_symbol)
    current_rdf <- isolate(gseaToUse$rdf)
    current_rdf_symbol <- isolate(gseaToUse$rdf_symbol)
    current_er <- isolate(gseaToUse$er)
    current_er_symbol <- isolate(gseaToUse$er_symbol)
    is_native <- isolate(gseaToUse$is_native)
    req(!is.null(input$selectedTable_GSEA_rows_selected))
    req(current_rdf)

    # Only pathways passing the term FDR cut-off (green in the table) are
    # plotted -- enrichplot already silently drops the rest from the native
    # plots, so apply the same rule to every figure and say so.
    sel_rows <- input$selectedTable_GSEA_rows_selected
    sel_rows <- sel_rows[!is.na(current_rdf$p.adjust[sel_rows]) & current_rdf$p.adjust[sel_rows] <= param$fdrThreshGSEA]
    n_dropped <- length(input$selectedTable_GSEA_rows_selected) - length(sel_rows)
    if (n_dropped > 0) {
      showNotification(
        paste0(n_dropped, " selected pathway(s) with FDR > ", param$fdrThreshGSEA, " not added to the figures."),
        id = "gseaFdrDrop",
        type = "warning"
      )
    }
    if (!length(sel_rows)) {
      for (x in names(gseaPlotReactive)) {
        gseaPlotReactive[[x]] <- NULL
      }
      return()
    }

    # Row-subset for the selected pathways, as a plain data.frame either way.
    # For the native path this is identical in content to
    # clusterProfiler::slice(current_er, ...)@result (used below wherever the
    # real S4 object -- and its geneList/geneSets -- is still needed).
    erRP_df <- current_rdf[sel_rows, , drop = FALSE]
    labelColumn <- ifelse(input$gseaType == "C2", "Description", "Label")

    # Apply gene filtering based on visualisation thresholds
    log2RatioGSEA <- gseaToUse$log2RatioGSEA
    # For isoform data, aggregate multiple transcripts per gene
    if (!is.null(param$featureLevel) && param$featureLevel == "isoform") {
      agg_method <- if (!is.null(input$gseaIsoformAgg)) input$gseaIsoformAgg else "maxabs"
      na_safe <- function(fn) {
        function(x) {
          x <- x[!is.na(x)]
          if (length(x) == 0) NA_real_ else fn(x)
        }
      }
      agg_fn <- switch(
        agg_method,
        "maxabs" = na_safe(function(x) x[which.max(abs(x))]),
        "mean" = na_safe(mean),
        "median" = na_safe(median),
        "max" = na_safe(max),
        "min" = na_safe(min)
      )
      log2RatioGSEA <- vapply(split(log2RatioGSEA, names(log2RatioGSEA)), agg_fn, numeric(1))
      log2RatioGSEA <- sort(log2RatioGSEA, decreasing = TRUE)
    }
    log2RatioGSEA_full <- log2RatioGSEA # Keep full list for reference
    genes_pass_filter <- names(log2RatioGSEA) # Default: all genes pass

    if (!is.null(input$gseaVisFilterMode) && input$gseaVisFilterMode != "none") {
      # Per-gene p-values for optional visualisation filtering, from the
      # enrichment input's annotation; absent for uploads -> filter degrades off.
      seqAnno <- ei$seq_anno

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
          if (input$gseaVisFilterMode == "genes") {
            if (length(genes_pass_filter) == 0) {
              showNotification(
                "No genes meet the current filter thresholds. Lower the Log2FC or P-Value threshold.",
                type = "warning"
              )
              return()
            }
            log2RatioGSEA <- log2RatioGSEA[names(log2RatioGSEA) %in% genes_pass_filter]
          }
        }
      }
      # For "labels" mode, keep all genes in log2RatioGSEA but use genes_pass_filter for labeling
    }

    # In "genes" mode, also remove filtered genes from the enrichment result
    # so cnetplot doesn't render them as grey nodes
    if (!is.null(input$gseaVisFilterMode) && input$gseaVisFilterMode == "genes") {
      current_rdf_symbol$core_enrichment <- sapply(strsplit(current_rdf_symbol$core_enrichment, "/"), function(g) {
        g <- g[g %in% genes_pass_filter]
        paste(g, collapse = "/")
      })
      # Check that selected pathways still have genes after filtering
      sel_genes <- current_rdf_symbol$core_enrichment[sel_rows]
      if (!any(nchar(sel_genes) > 0)) {
        showNotification(
          "No genes in the selected pathways meet the current filter thresholds. Lower the Log2FC or P-Value threshold.",
          type = "warning"
        )
        return()
      }
      if (is_native) current_er_symbol@result <- current_rdf_symbol
    }

    # Network plot and ridge plot both need the real S4 object -- geneSets for
    # cnetplot, geneList+geneSets for ridgeplot -- so both are native-only (see
    # plans/exploredelight/enrichment_detach_plan.md §9). Skip with a notification for
    # tidy-only pathway data rather than crashing.
    if (is_native) {
      cn <- enrichplot::cnetplot(
        x = current_er_symbol,
        node_label = "none",
        foldChange = log2RatioGSEA,
        showCategory = current_er_symbol@result$Description[sel_rows]
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
      cd1$Label[seq_along(sel_rows)] <-
        current_er_symbol@result[[labelColumn]][sel_rows]

      cd2 <- cn$data
      cd2$name[is.na(cd2$foldChange)] <- NA
      # Apply gene label filtering if in "labels" mode
      if (!is.null(input$gseaVisFilterMode) && input$gseaVisFilterMode == "labels") {
        cd2$name[!(cd2$name %in% genes_pass_filter)] <- NA
      }

      # Add gene labels if requested
      if (input$showGeneLabelsGSEA) {
        cn <- cn +
          geom_text_repel(
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
      cn <- cn +
        geom_label_repel(
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
      erRP <- clusterProfiler::slice(current_er, sel_rows)
      rp <- enrichplot::ridgeplot(x = erRP, fill = input$ridgePlotColourByGSEA)
      rp <- rp +
        scale_y_discrete(
          labels = erRP@result$Label[match(
            levels(as.factor(rp$plot_env$gs2val.df$category)),
            erRP@result$Description
          )]
        ) +
        theme(axis.text.y = element_text(vjust = -0.01, size = input$textSizeGSEA)) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_fill_gradient2(
          low = "dodgerblue4",
          high = "firebrick4",
          mid = "white",
          midpoint = 0
        )

      gseaPlotReactive$ridgeplot <- rp
    } else {
      gseaPlotReactive$cnetplot <- NULL
      gseaPlotReactive$ridgeplot <- NULL
      showNotification(
        "Network and ridge plots need the original enrichment object; not available for this dataset.",
        type = "message"
      )
    }

    # Lollipop plot (-log10 FDR, coloured by NES)
    lolli_df <- erRP_df
    lolli_df$negLogP <- -log10(lolli_df$p.adjust)
    lolli_df[[labelColumn]] <- factor(
      lolli_df[[labelColumn]],
      levels = lolli_df[[labelColumn]][order(lolli_df$negLogP)]
    )
    scl <- input$gseaLollipopScale
    ts <- input$textSizeGSEA
    lp <- ggplot(lolli_df, aes(y = .data[[labelColumn]], x = negLogP)) +
      geom_segment(aes(x = 0, xend = negLogP, y = .data[[labelColumn]], yend = .data[[labelColumn]]), linewidth = scl) +
      geom_point(aes(size = setSize, colour = NES)) +
      scale_colour_gradient2(low = "royalblue4", mid = "white", high = "firebrick4", midpoint = 0) +
      scale_size(range = c(ts / 3, ts) * scl, name = "Gene Set Size") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x = "-log10(FDR)", y = NULL, title = input$gseaType) +
      theme_prism(base_size = ts)
    gseaPlotReactive$lollipop <- lp

    # Heatmap - use filtered genes when filtering is active
    # Prepare matrix for ComplexHeatmap
    sel_idx <- sel_rows
    sel_er <- current_rdf_symbol[sel_idx, ]

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

    # GSEA running score plot -- needs the real S4 object's full geneList, so
    # native-only (see plans/exploredelight/enrichment_detach_plan.md §9).
    if (is_native) {
      rs <- enrichplot::gseaplot2(
        x = current_er,
        geneSetID = sel_rows,
        color = RColorBrewer::brewer.pal(max(3, length(sel_rows)), "Set2")[seq_along(sel_rows)]
      )
      gseaPlotReactive$runningscore <- rs
    } else {
      gseaPlotReactive$runningscore <- NULL
      showNotification(
        "Running-score plot needs the original enrichment object; not available for this dataset.",
        type = "message"
      )
    }

    # UpSet plot
    geneSetsUpset <- setNames(
      lapply(sel_rows, function(i) {
        colToUse <- if (input$gseaType %in% c("BP", "MF", "CC")) "geneName" else "core_enrichment"
        g <- current_rdf[[colToUse]][i] %>%
          strsplit("/") %>%
          .[[1]]
        unique(g[nchar(g) > 1])
      }),
      current_rdf[[labelColumn]][sel_rows]
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
    draw_heatmap(gseaPlotReactive[["heatmap"]], legend_dx = input$GSEALegendDX, legend_dy = input$GSEALegendDY)
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
      plot_object = function() {
        draw_heatmap(gseaPlotReactive[["heatmap"]], legend_dx = input$GSEALegendDX, legend_dy = input$GSEALegendDY)
      },
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

output$lollipopPlot_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["lollipop"]]))
    gseaPlotReactive[["lollipop"]]
  },
  width = function() as.numeric(input$figWidthGSEA),
  height = function() as.numeric(input$figHeightGSEA)
)

output$dlPlot_lollipopGSEA <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_GSEA, "_", input$gseaType, "_lollipop"),
      design_name = design,
      extension = tolower(input$downloadFormatGSEA)
    )
  },
  content = function(file) {
    req(gseaPlotReactive[["lollipop"]])
    export_plot_standard(
      plot_object = gseaPlotReactive[["lollipop"]],
      filename = file,
      format = tolower(input$downloadFormatGSEA),
      width = as.numeric(input$figWidthGSEA) / 50,
      height = as.numeric(input$figHeightGSEA) / 50,
      dpi = as.numeric(input$dpiGSEA)
    )
  }
)
