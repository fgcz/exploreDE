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

eds <- inputDataReactive()$.eds
req(eds)

# --- Central reactive store ---
featureBucket <- reactiveValues(
  features = character(0),
  version = 0L
)

# --- API functions ---

bucket_add <- function(new_features, source = "unknown") {
  # Counts-only uploads have no DE table (sa is NULL); the valid feature
  # universe is then simply every feature in the dataset.
  valid_ids <- if (!is.null(seqAnnoReactive$sa)) {
    seqAnnoReactive$sa$feature_id
  } else {
    eds$feature_ids()
  }
  # Also accept gene-level names for CrisprCountQC collapsed mode
  sg <- eds$sgrna_to_gene()
  if (!is.null(sg)) {
    valid_ids <- c(valid_ids, unique(sg))
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
  gene_lookup    <- eds$proteomics_choices("gene")$lookup
  protein_lookup <- eds$proteomics_choices("protein")$lookup
  if (search_mode == "gene_name" && !is.null(gene_lookup)) {
    expanded <- unlist(gene_lookup[input_genes], use.names = FALSE)
    direct <- input_genes[input_genes %in% seqAnnoReactive$sa$feature_id]
    return(unique(c(expanded, direct)))
  }
  if (search_mode == "protein_Id" && !is.null(protein_lookup)) {
    expanded <- unlist(protein_lookup[input_genes], use.names = FALSE)
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
        # Render an empty bucket_list so that input$fb_keep_<tab> is actively
        # set to character(0). A plain div() would leave the old input value
        # cached, causing the debounced reactive to miss changes when the same
        # feature is re-added after a bucket clear (e.g. search mode switch).
        hidden_bucket <- bucket_list(
          header = NULL,
          group_name = paste0("fb_group_", tab),
          orientation = "horizontal",
          add_rank_list(
            text = NULL,
            labels = NULL,
            input_id = paste0("fb_keep_", tab)
          )
        )
        # bucket_list() returns a tagList, which tagAppendAttributes() would
        # print as text; hide it with a wrapper instead.
        hidden_bucket <- div(style = "display:none;", hidden_bucket)
        tagList(
          div(
            class = "text-muted text-center", style = "padding: 20px;",
            icon("hand-pointer"),
            " Select features from any tab to populate the bucket."
          ),
          hidden_bucket
        )
      } else {
        # Long IDs (e.g. phospho sites with several peptides) are cut off with
        # an ellipsis; hovering shows the full ID. Named labels make sortable
        # return the feature_id itself (data-rank-id), not the visible text.
        labels <- setNames(lapply(features, function(f) tags$span(title = f, f)), features)
        div(
          class = "fb-bucket",
          tags$style(HTML(
            ".fb-bucket .rank-list-item { max-width: 16em; overflow: hidden;
               text-overflow: ellipsis; white-space: nowrap; }"
          )),
          bucket_list(
            header = "Drag and drop features (Ctrl/Cmd to select multiple)",
            group_name = paste0("fb_group_", tab),
            orientation = "horizontal",
            add_rank_list(
              text = "Include these features in this order",
              labels = labels,
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

# Boxplot, volcano and heatmap tabs: selectize + text box (the volcano brush is
# handled via the staging area). Each tab remembers what its inputs last
# contributed, so a feature removed from the selectize (its "x") or deleted
# from the text box also leaves the bucket, not just the input.
lapply(c("boxplot", "volcano", "heatmap"), function(tab) {
  sel_id  <- paste0("fb_select_", tab)
  text_id <- paste0("fb_text_", tab)
  contributed <- character(0)
  observeEvent(
    { input[[sel_id]]; input[[text_id]] },
    ignoreNULL = FALSE, ignoreInit = TRUE,
    {
      search_mode <- input[[paste0("fb_search_by_", tab)]] %||% "feature_id"
      genes <- c(input[[sel_id]] %||% character(0),
                 parse_gene_text(input[[text_id]]))
      genes <- resolve_features(genes, search_mode)
      bucket_remove(setdiff(contributed, genes), source = tab)
      bucket_add(genes, source = tab)
      contributed <<- genes
    }
  )
})

# DEG table tab (row selection)
observeEvent(input$degTable_rows_selected, ignoreNULL = FALSE, ignoreInit = TRUE, {
  sa <- seqAnnoReactive$sa
  if (!eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
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
    # Clear ALL selectize and textarea inputs across ALL tabs,
    # respecting the current search mode so choices stay consistent
    mode <- input[[paste0("fb_search_by_", tab)]] %||% "feature_id"
    gene_ch    <- eds$proteomics_choices("gene")$choices
    protein_ch <- eds$proteomics_choices("protein")$choices
    choices <- if (mode == "gene_name" && !is.null(gene_ch)) {
      gene_ch
    } else if (mode == "protein_Id" && !is.null(protein_ch)) {
      protein_ch
    } else {
      eds$feature_ids()
    }
    for (t in fb_input_tabs) {
      updateSelectizeInput(session, paste0("fb_select_", t),
                           choices = choices, selected = NULL, server = TRUE)
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
# Feature label mode (proteomics): which ID to print on heatmap rows / biplot
# arrows. Inverts the proteomics_choices lookups (label -> feature_ids) into a
# feature_id -> label map, so it works for both the SE and h5ad backends.
# =============================================================================

feature_label_map <- function(mode = "feature_id") {
  kind <- switch(mode, gene_name = "gene", protein_Id = "protein", NULL)
  if (is.null(kind)) return(NULL)
  lookup <- eds$proteomics_choices(kind)$lookup
  if (is.null(lookup)) return(NULL)
  setNames(rep(names(lookup), lengths(lookup)), unlist(lookup, use.names = FALSE))
}

# Map feature_ids to display labels; falls back to the feature_id when a
# feature has no label in the chosen mode.
label_features <- function(ids, mode = "feature_id") {
  map <- feature_label_map(mode %||% "feature_id")
  if (is.null(map)) return(ids)
  labs <- unname(map[ids])
  ifelse(is.na(labs) | !nzchar(labs), ids, labs)
}

# Plotting copy of a matrix with rownames swapped for display labels. Relabels
# rownames rather than passing Heatmap(row_labels = ...) because ComplexHeatmap
# errors on any explicit row_labels for a 0-row matrix (e.g. empty custom bucket).
label_rows <- function(mat, mode = "feature_id", width = NULL) {
  if (!NROW(mat)) return(mat)
  labs <- label_features(rownames(mat), mode)
  if (!is.null(width)) labs <- substr(labs, 1, width)
  rownames(mat) <- labs
  mat
}

lapply(c("heatmap", "ordination", "volcano", "multiTile"), function(tab) {
  output[[paste0("fb_label_mode_", tab)]] <- renderUI({
    has_gene    <- !is.null(eds$proteomics_choices("gene")$lookup)
    has_protein <- !is.null(eds$proteomics_choices("protein")$lookup)
    req(has_gene || has_protein)

    choices <- c("Feature ID" = "feature_id")
    if (has_protein) choices <- c(choices, "Protein ID" = "protein_Id")
    if (has_gene) choices <- c(choices, "Gene Symbol" = "gene_name")

    radioButtons(
      inputId = paste0("fb_label_by_", tab),
      label = "Label features by:",
      choices = choices,
      selected = if (has_gene) "gene_name" else "feature_id",
      inline = TRUE
    )
  })
  # Render even while its settings tab is hidden: plots read the radio's value,
  # and an unrendered radio is NULL, so the Gene Symbol default would silently
  # fall back to feature IDs until someone opened that tab.
  outputOptions(output, paste0("fb_label_mode_", tab), suspendWhenHidden = FALSE)
})

# =============================================================================
# Phase 8: Gene Symbol Search (sub-protein proteomics)
# =============================================================================

# Selectize choices for a "Search by" mode; `default` is used in feature_id mode
# (or when the chosen lookup is absent).
search_mode_choices <- function(mode, default) {
  gene_ch    <- eds$proteomics_choices("gene")$choices
  protein_ch <- eds$proteomics_choices("protein")$choices
  if (identical(mode, "gene_name") && !is.null(gene_ch)) return(gene_ch)
  if (identical(mode, "protein_Id") && !is.null(protein_ch)) return(protein_ch)
  default
}

# Multi DE tabs manage their own selectize / clearing (see server-multi-*.R).
lapply(c(fb_input_tabs, "correlation", "multiBoxplot", "multiTile"), function(tab) {
  output[[paste0("fb_search_mode_", tab)]] <- renderUI({
    # Only show radio when at least one alternative lookup exists
    has_gene    <- !is.null(eds$proteomics_choices("gene")$lookup)
    has_protein <- !is.null(eds$proteomics_choices("protein")$lookup)
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
    if (isTRUE(restoring())) return()  # save-state restore guard
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
    gene_ch    <- eds$proteomics_choices("gene")$choices
    protein_ch <- eds$proteomics_choices("protein")$choices
    choices <- if (mode == "gene_name" && !is.null(gene_ch)) {
      gene_ch
    } else if (mode == "protein_Id" && !is.null(protein_ch)) {
      protein_ch
    } else {
      eds$feature_ids()
    }
    for (t in fb_input_tabs) {
      updateSelectizeInput(session, paste0("fb_select_", t),
                           choices = choices, selected = NULL, server = TRUE)
    }
  })
})
