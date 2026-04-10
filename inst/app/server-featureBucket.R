# =============================================================================
# Central Feature Bucket — Single Source of Truth
# =============================================================================
# Sourced BEFORE all tab server files. Provides the shared feature bucket
# reactive, API functions, centralised bucket UI rendering, and reset handling.
#
# Architecture:
#   - featureBucket$features: canonical set of features (insertion order)
#   - Each tab renders its own bucket_list from featureBucket$features
#   - Each tab's consumer reads input$fb_keep_<tab> for LOCAL ordering/exclusion
#   - Reorder/exclude is LOCAL per tab (no cross-tab sync — avoids stale-input loops)
#   - Adding features (selectize, textarea, DEG table, volcano brush) goes through
#     bucket_add() which updates the central store → all tabs re-render

# --- Central reactive store ---
featureBucket <- reactiveValues(
  features = character(0),
  version = 0L
)

# --- API functions ---

bucket_add <- function(new_features, source = "unknown") {
  valid_ids <- seqAnnoReactive$sa$feature_id
  # Also accept gene-level names for CrisprCountQC collapsed mode
  if (!is.null(inputDataReactive()$sgRNA_to_gene)) {
    valid_ids <- c(valid_ids, unique(inputDataReactive()$sgRNA_to_gene))
  }
  valid <- new_features[new_features %in% valid_ids & new_features != ""]
  current <- isolate(featureBucket$features)
  to_add <- setdiff(valid, current)
  if (length(to_add) > 0) {
    featureBucket$features <- c(current, to_add)
    featureBucket$version <- isolate(featureBucket$version) + 1L
  }
}

bucket_remove <- function(to_remove, source = "unknown") {
  current <- isolate(featureBucket$features)
  remaining <- current[!current %in% to_remove]
  if (length(remaining) != length(current)) {
    featureBucket$features <- remaining
    featureBucket$version <- isolate(featureBucket$version) + 1L
  }
}

bucket_clear <- function(source = "unknown") {
  featureBucket$features <- character(0)
  featureBucket$version <- isolate(featureBucket$version) + 1L
}

parse_gene_text <- function(text) {
  if (is.null(text) || text == "") return(character(0))
  genes <- strsplit(text, "\\ +|\\/+|\\,+|\n+|\r+|\n\r+|\t+")[[1]]
  genes[genes != ""]
}

# Resolve gene symbols / protein IDs to feature_ids when not in feature_id mode
resolve_features <- function(input_genes, search_mode = "feature_id") {
  if (search_mode == "gene_name" && !is.null(inputDataReactive()$gene_lookup)) {
    expanded <- unlist(inputDataReactive()$gene_lookup[input_genes], use.names = FALSE)
    direct <- input_genes[input_genes %in% seqAnnoReactive$sa$feature_id]
    return(unique(c(expanded, direct)))
  }
  if (search_mode == "protein_Id" && !is.null(inputDataReactive()$protein_lookup)) {
    expanded <- unlist(inputDataReactive()$protein_lookup[input_genes], use.names = FALSE)
    direct <- input_genes[input_genes %in% seqAnnoReactive$sa$feature_id]
    return(unique(c(expanded, direct)))
  }
  return(input_genes)
}

# --- Centralised bucket UI rendering (one per participating tab) ---
fb_tabs <- c("boxplot", "volcano", "heatmap", "deg")

lapply(fb_tabs, function(tab) {
  observe({
    features <- featureBucket$features
    output[[paste0("fb_bucket_", tab)]] <- renderUI({
      if (length(features) == 0) {
        div(
          class = "text-muted text-center", style = "padding: 20px;",
          icon("hand-pointer"),
          " Select features from any tab to populate the bucket."
        )
      } else {
        bucket_list(
          header = "Drag and drop features (Ctrl/Cmd to select multiple)",
          group_name = paste0("fb_group_", tab),
          orientation = "horizontal",
          add_rank_list(
            text = "Include these features in this order",
            labels = features,
            input_id = paste0("fb_keep_", tab),
            options = sortable::sortable_options(multiDrag = TRUE)
          ),
          add_rank_list(
            text = "Exclude these features",
            labels = NULL,
            input_id = paste0("fb_exclude_", tab),
            options = sortable::sortable_options(multiDrag = TRUE)
          )
        )
      }
    })
  })
})

# NOTE: No reorder/exclude sync observers. Each tab's drag-drop order and
# exclusions are LOCAL. This avoids infinite loops caused by stale input
# values when bucket_list re-renders (input$fb_keep_* fires with the old
# value before the client receives the new UI, which triggers bucket_replace,
# which drops features, which triggers re-add from DEG table, etc.).
# Consumers read input$fb_keep_<tab> for their local order.

# --- Per-tab input observers (contributors) ---

# Boxplot tab
observeEvent(
  { input$fb_select_boxplot; input$fb_text_boxplot },
  ignoreNULL = FALSE, ignoreInit = TRUE,
  {
    search_mode <- input$fb_search_by_boxplot %||% "feature_id"
    genes <- c(input$fb_select_boxplot %||% character(0),
               parse_gene_text(input$fb_text_boxplot))
    genes <- resolve_features(genes, search_mode)
    bucket_add(genes, source = "boxplot")
  }
)

# Volcano tab (selectize + text only; brush is handled via staging area)
observeEvent(
  { input$fb_select_volcano; input$fb_text_volcano },
  ignoreNULL = FALSE, ignoreInit = TRUE,
  {
    search_mode <- input$fb_search_by_volcano %||% "feature_id"
    genes <- c(input$fb_select_volcano %||% character(0),
               parse_gene_text(input$fb_text_volcano))
    genes <- resolve_features(genes, search_mode)
    bucket_add(genes, source = "volcano")
  }
)

# Heatmap tab
observeEvent(
  { input$fb_select_heatmap; input$fb_text_heatmap },
  ignoreNULL = FALSE, ignoreInit = TRUE,
  {
    search_mode <- input$fb_search_by_heatmap %||% "feature_id"
    genes <- c(input$fb_select_heatmap %||% character(0),
               parse_gene_text(input$fb_text_heatmap))
    genes <- resolve_features(genes, search_mode)
    bucket_add(genes, source = "heatmap")
  }
)

# DEG table tab (row selection)
observeEvent(input$degTable_rows_selected, ignoreNULL = FALSE, ignoreInit = TRUE, {
  sa <- seqAnnoReactive$sa
  if (!inputDataReactive()$dataType %in% c("CountQC", "CrisprCountQC")) {
    sa <- sa[order(sa$pValue), ]
  }
  selected <- sa$feature_id[input$degTable_rows_selected]
  bucket_add(selected, source = "deg_table")
})

# --- Volcano brush staging area handler ---
observeEvent(input$fb_add_brushed, ignoreInit = TRUE, {
  req(length(volcanoStagedGenes$brushGenes) > 0)
  bucket_add(volcanoStagedGenes$brushGenes, source = "volcano_brush")
  showNotification(
    paste0("Added ", length(volcanoStagedGenes$brushGenes), " features to bucket"),
    type = "message", duration = 3
  )
})

# --- Centralised reset handler ---
fb_input_tabs <- c("boxplot", "volcano", "heatmap")
fb_reset_tabs <- c("boxplot", "volcano", "heatmap", "deg")

lapply(fb_reset_tabs, function(tab) {
  observeEvent(input[[paste0("fb_reset_", tab)]], ignoreInit = TRUE, {
    bucket_clear(source = paste0("reset_", tab))
    # Clear ALL selectize and textarea inputs across ALL tabs
    for (t in fb_input_tabs) {
      updateSelectizeInput(session, paste0("fb_select_", t),
                           choices = inputDataReactive()$genes, selected = NULL, server = TRUE)
      updateTextAreaInput(session, paste0("fb_text_", t), value = "")
    }
    # Clear volcano staged genes
    volcanoStagedGenes$brushGenes <- character(0)
    # Null out figure reactives to clear all bucket-driven plots
    figuresDataReactive$boxplot <- NULL
    figuresDataReactive$barplot <- NULL
    figuresDataReactive$heatmapCustom <- NULL
    figuresDataReactive$volcanoStatic <- NULL
    figuresDataReactive$volcanoMA <- NULL
    figuresDataReactive$volcanoPlotly <- NULL
  })
})

# =============================================================================
# Phase 8: Gene Symbol Search (sub-protein proteomics)
# =============================================================================

lapply(fb_input_tabs, function(tab) {
  output[[paste0("fb_search_mode_", tab)]] <- renderUI({
    # Only show radio when at least one alternative lookup exists
    has_gene <- !is.null(inputDataReactive()$gene_lookup)
    has_protein <- !is.null(inputDataReactive()$protein_lookup)
    req(has_gene || has_protein)

    choices <- c("Feature ID" = "feature_id")
    if (has_protein) choices <- c(choices, "Protein ID" = "protein_Id")
    if (has_gene) choices <- c(choices, "Gene Symbol" = "gene_name")

    radioButtons(
      inputId = paste0("fb_search_by_", tab),
      label = "Search by:",
      choices = choices,
      selected = "feature_id",
      inline = TRUE
    )
  })
})

lapply(fb_input_tabs, function(tab) {
  observeEvent(input[[paste0("fb_search_by_", tab)]], ignoreInit = TRUE, {
    mode <- input[[paste0("fb_search_by_", tab)]]

    # Clear everything — search mode change invalidates bucket contents
    bucket_clear(source = paste0("search_mode_", tab))
    for (t in fb_input_tabs) {
      updateTextAreaInput(session, paste0("fb_text_", t), value = "")
    }
    volcanoStagedGenes$brushGenes <- character(0)
    # Clear bucket-driven plots (boxplot, heatmap) but NOT volcano/MA
    # — volcano renders all points independently; only labels come from bucket
    figuresDataReactive$boxplot <- NULL
    figuresDataReactive$barplot <- NULL
    figuresDataReactive$heatmapCustom <- NULL

    # Update selectize choices for the new mode across all tabs
    choices <- if (mode == "gene_name") {
      sort(names(inputDataReactive()$gene_lookup))
    } else if (mode == "protein_Id") {
      sort(names(inputDataReactive()$protein_lookup))
    } else {
      inputDataReactive()$genes
    }
    for (t in fb_input_tabs) {
      updateSelectizeInput(session, paste0("fb_select_", t),
                           choices = choices, selected = NULL, server = TRUE)
    }
  })
})
