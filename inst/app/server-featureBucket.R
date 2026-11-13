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
#     fb_add() which updates the central store → all tabs re-render

eds <- inputDataReactive()$.eds
req(eds)

# --- Central reactive store ---
featureBucket <- reactiveValues(
  features = character(0),
  version = 0L
)

# "Search by" mode (Feature ID / Gene Symbol / Protein ID) of the bucket tabs
# (fb_input_tabs). ONE app-wide value, like activeContrast: A switch on any tab
# clears the shared bucket and re-points every tab's selectize, so all their
# radios have to show the same mode.
fbSearchMode <- reactiveVal("feature_id")

# --- API functions ---
# The session's bucket and lookups around the pure helpers bucket_add(),
# bucket_remove() and fb_resolve() (R/). seqAnnoReactive is set by
# server-boxplot.R, sourced after this file; it is read when they run.

# Ids that may enter the bucket: the DE table's, or every feature of a
# counts-only upload (no DE table); for CRISPR also the genes guides collapse to.
fb_known_ids <- function() {
  ids <- if (!is.null(seqAnnoReactive$sa)) seqAnnoReactive$sa$feature_id else eds$feature_ids()
  c(ids, unique(eds$sgrna_to_gene()))
}

# Store a new bucket; the version counts the changes.
fb_set <- function(features) {
  if (!identical(features, isolate(featureBucket$features))) {
    featureBucket$features <- features
    featureBucket$version <- isolate(featureBucket$version) + 1L
  }
}

fb_add <- function(new) fb_set(bucket_add(isolate(featureBucket$features), new, fb_known_ids()))
fb_remove <- function(ids) fb_set(bucket_remove(isolate(featureBucket$features), ids))

# Empty the bucket; the version changes even when it was empty already.
fb_clear <- function() {
  featureBucket$features <- character(0)
  featureBucket$version <- isolate(featureBucket$version) + 1L
}

fb_resolve <- function(tokens, mode = "feature_id") {
  lookups <- list(gene = eds$proteomics_choices("gene")$lookup, protein = eds$proteomics_choices("protein")$lookup)
  resolve_features(tokens, mode, lookups, seqAnnoReactive$sa$feature_id)
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
            class = "text-muted text-center",
            style = "padding: 20px;",
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

# Boxplot, volcano and heatmap tabs: selectize + text box + pattern box (the
# volcano brush is handled via the staging area). Each tab remembers what its
# inputs last contributed, so a feature removed from the selectize (its "x"),
# deleted from the text box or no longer matched by the pattern also leaves
# the bucket, not just the input. The pattern is matched against the names
# the tab's picker offers (feature IDs, gene symbols or protein IDs; gene
# names when CRISPR guides are collapsed), so its matches resolve to feature
# IDs like typed names do.
lapply(c("boxplot", "volcano", "heatmap"), function(tab) {
  sel_id <- paste0("fb_select_", tab)
  text_id <- paste0("fb_text_", tab)
  collapse_id <- paste0(tab, "CollapseGuides")
  # The pattern box sends on blur where shiny supports it; the debounce covers
  # older shiny, where it sends every keystroke.
  pattern_val <- debounce(reactive(input[[paste0("fb_pattern_", tab)]]), 600)
  regex_val <- debounce(reactive(input[[paste0("fb_pattern_regex_", tab)]]), 600)
  pattern_info <- reactiveVal("")
  contributed <- character(0)
  observeEvent(
    {
      input[[sel_id]]
      input[[text_id]]
      pattern_val()
      regex_val()
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE,
    {
      search_mode <- input[[paste0("fb_search_by_", tab)]] %||% "feature_id"
      pattern <- trimws(pattern_val() %||% "")
      matched <- if (nzchar(pattern)) {
        candidates <- if (eds$data_type() == "CrisprCountQC" && identical(input[[collapse_id]], "gene")) {
          sort(unique(eds$sgrna_to_gene()))
        } else {
          search_mode_choices(search_mode, eds$feature_ids())
        }
        match_feature_pattern(pattern, candidates, regex = isTRUE(regex_val()))
      } else {
        character(0)
      }
      # The cap counts features: a protein ID or gene symbol can stand for
      # many peptides or sites.
      pattern_ids <- cap_pattern_features(matched, function(x) fb_resolve(x, search_mode))
      pattern_info(fb_pattern_message(pattern, matched, pattern_ids))
      genes <- c(input[[sel_id]] %||% character(0), parse_gene_text(input[[text_id]]))
      genes <- unique(c(fb_resolve(genes, search_mode), pattern_ids))
      fb_remove(setdiff(contributed, genes))
      fb_add(genes)
      contributed <<- genes
    }
  )
  output[[paste0("fb_pattern_info_", tab)]] <- renderText(pattern_info())
})

# Empty the selectize, text and pattern boxes of every bucket tab, offering
# `choices` in the selectize (the reset and the search-mode handlers below).
fb_clear_inputs <- function(choices) {
  for (t in fb_input_tabs) {
    updateSelectizeInput(session, paste0("fb_select_", t), choices = choices, selected = NULL, server = TRUE)
    updateTextAreaInput(session, paste0("fb_text_", t), value = "")
    updateTextInput(session, paste0("fb_pattern_", t), value = "")
  }
}

# DEG table tab (row selection)
observeEvent(input$degTable_rows_selected, ignoreNULL = FALSE, ignoreInit = TRUE, {
  sa <- seqAnnoReactive$sa
  if (!eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
    sa <- sa[order(sa$pValue), ]
  }
  selected <- sa$feature_id[input$degTable_rows_selected]
  fb_add(selected)
})

# --- Volcano brush staging area handler ---
observeEvent(input$fb_add_brushed, ignoreInit = TRUE, {
  req(length(volcanoStagedGenes$brushGenes) > 0)
  fb_add(volcanoStagedGenes$brushGenes)
  showNotification(
    paste0("Added ", length(volcanoStagedGenes$brushGenes), " features to bucket"),
    type = "message",
    duration = 3
  )
})

# --- Centralised reset handler ---
fb_input_tabs <- c("boxplot", "volcano", "heatmap")
fb_reset_tabs <- c("boxplot", "volcano", "heatmap", "deg")

lapply(fb_reset_tabs, function(tab) {
  observeEvent(input[[paste0("fb_reset_", tab)]], ignoreInit = TRUE, {
    fb_clear()
    # Clear the inputs of ALL tabs, respecting the current search mode so
    # choices stay consistent
    fb_clear_inputs(search_mode_choices(fbSearchMode(), eds$feature_ids()))
    # Clear volcano staged genes
    volcanoStagedGenes$brushGenes <- character(0)
    # Boxplot / barplot / custom heatmap are NOT nulled here: They follow their
    # own tab's fb_keep_<tab>, which empties when that tab next shows the
    # bucket. Nulling them from here, while their tab is hidden, left them
    # blank for good if the bucket was refilled with the same features
    # (fb_keep_<tab> then never changes, so nothing redraws them).
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
  if (is.null(kind)) {
    return(NULL)
  }
  lookup <- eds$proteomics_choices(kind)$lookup
  if (is.null(lookup)) {
    return(NULL)
  }
  setNames(rep(names(lookup), lengths(lookup)), unlist(lookup, use.names = FALSE))
}

# Map feature_ids to display labels; falls back to the feature_id when a
# feature has no label in the chosen mode.
label_features <- function(ids, mode = "feature_id") {
  map <- feature_label_map(mode %||% "feature_id")
  if (is.null(map)) {
    return(ids)
  }
  labs <- unname(map[ids])
  ifelse(is.na(labs) | !nzchar(labs), ids, labs)
}

# Plotting copy of a matrix with rownames swapped for display labels. Relabels
# rownames rather than passing Heatmap(row_labels = ...) because ComplexHeatmap
# errors on any explicit row_labels for a 0-row matrix (e.g. empty custom bucket).
label_rows <- function(mat, mode = "feature_id", width = NULL) {
  if (!NROW(mat)) {
    return(mat)
  }
  labs <- label_features(rownames(mat), mode)
  if (!is.null(width)) {
    labs <- substr(labs, 1, width)
  }
  rownames(mat) <- labs
  mat
}

lapply(c("heatmap", "ordination", "volcano", "multiTile"), function(tab) {
  output[[paste0("fb_label_mode_", tab)]] <- renderUI({
    # Isoform feature_ids already carry the gene symbol (Gene_ENST...), and
    # labelling by symbol alone would make a gene's transcripts indistinguishable.
    req(!identical(eds$feature_id_kind(), "isoform"))
    has_gene <- !is.null(eds$proteomics_choices("gene")$lookup)
    has_protein <- !is.null(eds$proteomics_choices("protein")$lookup)
    req(has_gene || has_protein)

    choices <- c("Feature ID" = "feature_id")
    if (has_protein) {
      choices <- c(choices, "Protein ID" = "protein_Id")
    }
    if (has_gene) {
      choices <- c(choices, "Gene Symbol" = "gene_name")
    }

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
  gene_ch <- eds$proteomics_choices("gene")$choices
  protein_ch <- eds$proteomics_choices("protein")$choices
  if (identical(mode, "gene_name") && !is.null(gene_ch)) {
    return(gene_ch)
  }
  if (identical(mode, "protein_Id") && !is.null(protein_ch)) {
    return(protein_ch)
  }
  default
}

# Multi DE tabs manage their own selectize / clearing (see server-multi-*.R).
lapply(c(fb_input_tabs, "correlation", "multiBoxplot", "multiTile"), function(tab) {
  output[[paste0("fb_search_mode_", tab)]] <- renderUI({
    # Only show radio when at least one alternative lookup exists
    has_gene <- !is.null(eds$proteomics_choices("gene")$lookup)
    has_protein <- !is.null(eds$proteomics_choices("protein")$lookup)
    req(has_gene || has_protein)

    choices <- c("Feature ID" = "feature_id")
    if (has_protein) {
      choices <- c(choices, "Protein ID" = "protein_Id")
    }
    if (has_gene) {
      choices <- c(choices, "Gene Symbol" = "gene_name")
    }

    radioButtons(
      inputId = paste0("fb_search_by_", tab),
      label = "Search by:",
      choices = choices,
      # Bucket tabs open on the app-wide mode; a radio first drawn later (e.g.
      # Heatmap > Custom) must not reset it to Feature ID.
      selected = if (tab %in% fb_input_tabs) isolate(fbSearchMode()) else "feature_id",
      inline = TRUE
    )
  })
})

lapply(fb_input_tabs, function(tab) {
  observeEvent(input[[paste0("fb_search_by_", tab)]], ignoreInit = TRUE, {
    mode <- input[[paste0("fb_search_by_", tab)]]
    if (isTRUE(restoring())) {
      fbSearchMode(mode)
      return()
    } # save-state restore guard
    # Not a switch: The radio was just drawn, or was synced from another tab.
    # (Treating its first draw as a switch emptied the bucket and blanked the
    # boxplot whenever Heatmap > Custom was opened.)
    if (identical(mode, fbSearchMode())) {
      return()
    }
    fbSearchMode(mode)
    for (t in setdiff(fb_input_tabs, tab)) {
      updateRadioButtons(session, paste0("fb_search_by_", t), selected = mode)
    }

    # Clear everything — search mode change invalidates bucket contents.
    # Bucket-driven plots (boxplot, barplot, custom heatmap) follow their own
    # tab's bucket, so they are not nulled here (see the reset handler above).
    fb_clear()
    # Empty every tab's inputs, with the selectize choices of the new mode
    fb_clear_inputs(search_mode_choices(mode, eds$feature_ids()))
    volcanoStagedGenes$brushGenes <- character(0)
  })
})
