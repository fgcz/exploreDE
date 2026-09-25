## Count QC tab
## Plots (RNASeq / CountQC / CrisprCountQC only; the tab itself is hidden
## for proteomics in the sidebar):
##   - Total reads per sample
##   - Features above threshold (% y-axis)
##   - Dominant genes boxplot with family highlighter
##   - Intra-group scatter with 3 colour-by flavours
##     (GC / Width flavours require gc / featWidth in seqAnno)

eds <- inputDataReactive()$.eds
req(eds)

# ---- Sync Settings tabset to the Results tabset -----------------------------
observeEvent(input$countQcPlotsTabs, ignoreInit = TRUE, {
  mapping <- c(
    "Total Reads"               = "Bars",
    "Features Above Threshold"  = "Bars",
    "Dominant Genes"            = "Dominant",
    "Intra-group Scatter"       = "Intra-group"
  )
  target <- mapping[[input$countQcPlotsTabs]]
  if (!is.null(target)) {
    updateTabsetPanel(session, "countQcSettingsTabs", selected = target)
  }
})

# ---- Safe input accessors ----------------------------------------------------
cqW <- function() { w <- input$countQcFigWidth;  if (is.null(w)) 800 else w }
cqH <- function() { h <- input$countQcFigHeight; if (is.null(h)) 500 else h }
cqT <- function() { t <- input$countQcTextSize;  if (is.null(t)) 12  else t }

# ---- Gene-family presets -----------------------------------------------------
countQcFamilyPresets <- list(
  "Mitochondrial"    = "^MT-",
  "Ribosomal (cyto)" = "^RP[LS]",
  "Ribosomal (mito)" = "^MRP[LS]",
  "Haemoglobin"      = "^HB[ABEGMQZ]",
  "Immunoglobulin"   = "^IG[HKL][VDJC]",
  "TCR"              = "^TR[ABGD][VDJC]",
  "Heat shock"       = "^HSP",
  "Histone"          = "^HIST"
)

# ---- Populate the "colour bars by factor" selector ---------------------------
observe({
  factors <- eds$factor_names()
  req(length(factors) > 0)
  grouping <- eds$primary_grouping()
  default <- if (!is.null(grouping) && !is.na(grouping) && grouping %in% factors) grouping else factors[[1]]
  updateSelectInput(session, "countQcBarColourBy",
                    choices = factors, selected = default)
})

# ---- Sample / group reorder bucket for bar charts ---------------------------
# When the user switches factor, rebuild the include/exclude lists from that
# factor's levels. Same pattern used by the Boxplots and Heatmaps tabs.
observeEvent(input$countQcBarColourBy, ignoreInit = FALSE, {
  req(input$countQcBarColourBy)
  fc <- input$countQcBarColourBy
  dataset <- eds$samples()
  req(fc %in% colnames(dataset))
  lvls <- levels(as.factor(dataset[[fc]]))
  output$countQcBarBucketUI <- renderUI({
    # save-state restore: seed include/exclude split + order from a saved session.
    .keep <- restore_default("countQcBarKeepBucket", lvls, valid = lvls)
    .excl <- restore_default("countQcBarExcludeBucket", character(0), valid = lvls)
    .keep <- setdiff(.keep, .excl)
    sortable::bucket_list(
      header = "Drag levels to include / exclude / reorder (Ctrl/Cmd for multi-select)",
      group_name = "countQcBarBucket",
      orientation = "horizontal",
      sortable::add_rank_list(
        text = "Include these groups",
        labels = as.list(.keep),
        input_id = "countQcBarKeepBucket",
        options = sortable::sortable_options(multiDrag = TRUE)
      ),
      sortable::add_rank_list(
        text = "Exclude these groups",
        labels = as.list(.excl),
        input_id = "countQcBarExcludeBucket",
        options = sortable::sortable_options(multiDrag = TRUE)
      )
    )
  })
})

# ---- Facet selector + bucket for bar charts ---------------------------------
# Populate facet-by choices with "(none)" + all factors. The observer rebuilds
# the facet bucket whenever the user picks a new facet factor.
observe({
  factors <- eds$factor_names()
  req(length(factors) > 0)
  updateSelectInput(session, "countQcBarFacetBy",
                    choices = c("(none)", factors),
                    selected = isolate(input$countQcBarFacetBy) %||% "(none)")
})

observeEvent(input$countQcBarFacetBy, ignoreInit = FALSE, {
  fc <- input$countQcBarFacetBy
  if (is.null(fc) || fc == "(none)") {
    output$countQcBarFacetBucketUI <- renderUI({ NULL })
    return()
  }
  dataset <- eds$samples()
  req(fc %in% colnames(dataset))
  lvls <- levels(as.factor(dataset[[fc]]))
  output$countQcBarFacetBucketUI <- renderUI({
    # save-state restore: seed include/exclude split + order from a saved session.
    .keep <- restore_default("countQcBarFacetKeepBucket", lvls, valid = lvls)
    .excl <- restore_default("countQcBarFacetExcludeBucket", character(0), valid = lvls)
    .keep <- setdiff(.keep, .excl)
    sortable::bucket_list(
      header = "Facet levels to include / exclude / reorder",
      group_name = "countQcBarFacetBucket",
      orientation = "horizontal",
      sortable::add_rank_list(
        text = "Include these facets",
        labels = as.list(.keep),
        input_id = "countQcBarFacetKeepBucket",
        options = sortable::sortable_options(multiDrag = TRUE)
      ),
      sortable::add_rank_list(
        text = "Exclude these facets",
        labels = as.list(.excl),
        input_id = "countQcBarFacetExcludeBucket",
        options = sortable::sortable_options(multiDrag = TRUE)
      )
    )
  })
})

# ---- Populate count-type selector (intra-group) ------------------------------
observe({
  assay_keys <- eds$assay_names()
  req(length(assay_keys) > 0)
  count_choices <- setNames(assay_keys, vapply(assay_keys, eds$assay_label, character(1)))
  default <- if ("normalised" %in% assay_keys) "normalised" else assay_keys[[1]]
  updateSelectInput(session, "countQcIntraCountType",
                    choices = count_choices, selected = default)
})

# ---- Populate highlight selectize -------------------------------------------
observe({
  updateSelectizeInput(
    session = session, inputId = "countQcHighlightSelect",
    choices = eds$feature_ids(), selected = character(0),
    server = TRUE
  )
})

# ---- Intra-group: factor and level selectors --------------------------------
output$countQcIntraFactorUI <- renderUI({
  factors <- eds$factor_names()
  req(length(factors) > 0)
  grouping <- eds$primary_grouping()
  default <- if (!is.null(grouping) && !is.na(grouping) && grouping %in% factors) grouping else factors[[1]]
  selectInput("countQcIntraFactor", "Group by factor",
              choices = factors, selected = default)
})

output$countQcIntraGroupUI <- renderUI({
  req(input$countQcIntraFactor)
  dataset <- eds$samples()
  req(input$countQcIntraFactor %in% colnames(dataset))
  levels_avail <- unique(as.character(dataset[[input$countQcIntraFactor]]))
  selectInput("countQcIntraGroup", "Group level to display",
              choices = levels_avail,
              selected = restore_default("countQcIntraGroup", levels_avail[[1]], valid = levels_avail))
})

# ==============================================================================
# 1. Bar charts — total reads + features above threshold
# ==============================================================================

# Visibility: only RNASeq / CountQC (proteomics has no raw counts + presentFlag)
countQcBarVisible <- reactive({
  eds$data_type() %in% c("RNASeq", "CountQC", "CrisprCountQC")
})

# Bar / Dominant plot tabs live in a shared tabsetPanel; each render function
# uses validate() to show a friendly message when the data type doesn't apply.

# Shared bar-data reactive
# Filters samples to the current "Include" bucket (if populated) and orders
# samples by the user's chosen group order, then by original sample order
# within each group.
countQcBarData <- reactive({
  req(countQcBarVisible())
  # Prefer the unfiltered SE (rawData.rds) when available so the denominator
  # is the full transcriptome — matches the ezRun CountQC.Rmd "Genomic
  # features with reads above threshold" plot. Falls back to the main SE
  # for RNASeq (deResult.rds already carries all features) and for
  # CountQC/CrisprCountQC outputs that don't ship rawData.rds.
  # seUnfiltered is still attached as a legacy gStoreList field by initInputData;
  # falls back to eds$backend() (the filtered SE) when absent.
  seU <- inputDataReactive()$seUnfiltered
  if (!is.null(seU)) {
    req("counts" %in% assayNames(seU))
    raw_counts   <- assay(seU, "counts")
    present_flag <- if ("presentFlag" %in% assayNames(seU)) assay(seU, "presentFlag") else (raw_counts > 0) * 1L
  } else if (eds$backend_kind() == "se") {
    # SUSHI SE without a separate rawData.rds (e.g. RNASeq deResult): use the
    # filtered SE directly, exactly as before.
    se <- eds$backend()
    req("counts" %in% assayNames(se))
    raw_counts   <- assay(se, "counts")
    present_flag <- if ("presentFlag" %in% assayNames(se)) assay(se, "presentFlag") else (raw_counts > 0) * 1L
  } else {
    # Uploaded h5ad: no unfiltered SE and no presentFlag assay — use the
    # dataset's own counts via the accessor.
    req(eds$has_assay("counts"))
    raw_counts   <- eds$assay("counts")
    present_flag <- (raw_counts > 0) * 1L
  }

  dataset <- eds$samples()
  factorCol <- input$countQcBarColourBy %||% eds$primary_grouping()
  req(factorCol %in% colnames(dataset))

  # Per-sample group from dataset
  group_vec <- as.character(dataset[[factorCol]][match(colnames(raw_counts), rownames(dataset))])

  # Order of levels — from the Include bucket if present, else the dataset's natural order
  level_order <- input$countQcBarKeepBucket
  if (is.null(level_order) || !length(level_order)) {
    level_order <- levels(as.factor(dataset[[factorCol]]))
  }

  # Optional facet factor
  facetCol <- input$countQcBarFacetBy %||% "(none)"
  use_facet <- facetCol != "(none)" && facetCol %in% colnames(dataset) && facetCol != factorCol
  if (use_facet) {
    facet_vec <- as.character(dataset[[facetCol]][match(colnames(raw_counts), rownames(dataset))])
    facet_order <- input$countQcBarFacetKeepBucket
    if (is.null(facet_order) || !length(facet_order)) {
      facet_order <- levels(as.factor(dataset[[facetCol]]))
    }
  } else {
    facet_vec <- rep(NA_character_, length(group_vec))
    facet_order <- NA_character_
  }

  keep_samples <- group_vec %in% level_order &
                  (!use_facet | facet_vec %in% facet_order)
  raw_counts   <- raw_counts[, keep_samples, drop = FALSE]
  present_flag <- present_flag[, keep_samples, drop = FALSE]
  group_vec    <- group_vec[keep_samples]
  facet_vec    <- facet_vec[keep_samples]
  # Sort samples: by facet order (if set), then by group order
  if (use_facet) {
    sort_idx <- order(match(facet_vec, facet_order),
                      match(group_vec, level_order),
                      seq_along(group_vec))
  } else {
    sort_idx <- order(match(group_vec, level_order), seq_along(group_vec))
  }
  raw_counts   <- raw_counts[, sort_idx, drop = FALSE]
  present_flag <- present_flag[, sort_idx, drop = FALSE]
  group_vec    <- group_vec[sort_idx]
  facet_vec    <- facet_vec[sort_idx]

  n_features <- nrow(raw_counts)
  present_count <- colSums(present_flag)

  df <- tibble::tibble(
    Sample = factor(colnames(raw_counts), levels = colnames(raw_counts)),
    Total_M = colSums(raw_counts) / 1e6,
    Present = present_count,
    Pct_num = 100 * present_count / n_features,
    Group = factor(group_vec, levels = level_order)
  )
  if (use_facet) {
    df$FacetGroup <- factor(facet_vec, levels = facet_order)
  }
  attr(df, "use_facet") <- use_facet
  df
})

# Colours come from the app-wide GroupColour* inputs populated by the Colours
# tab. Falls back to black if a level's colour hasn't been instantiated yet.
countQcGroupColours <- reactive({
  req(input$countQcBarColourBy)
  factorCol <- input$countQcBarColourBy
  df <- countQcBarData()
  lvls <- levels(df$Group)
  cols <- vapply(lvls, function(lv) {
    v <- input[[paste0("GroupColour", factorCol, lv)]]
    if (is.null(v) || !nzchar(v)) "black" else v
  }, character(1))
  setNames(cols, lvls)
})

# Display labels from the Colours tab for the colour-by factor and the facet
# factor. If either is rich, both are (shared markdown theme).
countQcLevelLabels <- function(df) {
  groupCol <- input$countQcBarColourBy %||% eds$primary_grouping()
  facetCol <- input$countQcBarFacetBy
  has_facet <- "FacetGroup" %in% colnames(df)
  make <- function(force) list(
    group = level_labels(input, groupCol, levels(df$Group), force_rich = force),
    facet = if (has_facet) level_labels(input, facetCol, levels(df$FacetGroup), force_rich = force) else NULL
  )
  labs <- make(FALSE)
  rich <- is_rich_labels(labs$group) || is_rich_labels(labs$facet)
  if (rich) labs <- make(TRUE)
  labs$rich <- rich
  labs
}

maybe_facet <- function(p, df, lbl = NULL) {
  if (isTRUE(attr(df, "use_facet")) && "FacetGroup" %in% colnames(df)) {
    p <- p + facet_grid(cols = vars(FacetGroup),
                        scales = "free_x", space = "free_x",
                        labeller = if (is.null(lbl$facet)) "label_value" else as_labeller(lbl$facet))
  }
  if (isTRUE(lbl$rich)) p <- p + rich_label_theme(legend = TRUE, strip = TRUE)
  p
}

build_total_reads_plot <- function() {
  validate(need(isTRUE(try(countQcBarVisible(), silent = TRUE)),
                "Total-reads plot not applicable for this data type."))
  df <- countQcBarData()
  lbl <- countQcLevelLabels(df)
  p <- ggplot(df, aes(x = Sample, y = Total_M, fill = Group)) +
    geom_col(colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = countQcGroupColours(), labels = lbl$group) +
    ggprism::theme_prism(base_size = cqT(), axis_text_angle = 45) +
    labs(title = "Total Reads per Sample",
         y = "Counts (millions)", x = NULL)
  maybe_facet(p, df, lbl)
}

build_features_above_plot <- function() {
  validate(need(isTRUE(try(countQcBarVisible(), silent = TRUE)),
                "Features-above-threshold plot not applicable for this data type."))
  df <- countQcBarData()
  lbl <- countQcLevelLabels(df)
  p <- ggplot(df, aes(x = Sample, y = Pct_num, fill = Group)) +
    geom_col(colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = countQcGroupColours(), labels = lbl$group) +
    scale_y_continuous(
      labels = function(v) paste0(v, "%"),
      limits = c(0, 105), expand = expansion(mult = c(0, 0.02))
    ) +
    ggprism::theme_prism(base_size = cqT(), axis_text_angle = 45) +
    labs(title = "Features Above Threshold",
         y = "% features present", x = NULL)
  if (isTRUE(input$countQcBarShowPctLabel)) {
    p <- p + geom_text(
      aes(label = paste0(signif(Pct_num, 3), "%")),
      vjust = -0.3,
      size  = input$countQcBarPctLabelSize %||% 3.5
    )
  }
  maybe_facet(p, df, lbl)
}

output$countQcTotalReads <- renderPlot(
  build_total_reads_plot(), width = cqW, height = cqH
)

output$countQcFeaturesAbove <- renderPlot(
  build_features_above_plot(), width = cqW, height = cqH
)

# ==============================================================================
# 2. Dominant genes boxplot (RNASeq / CountQC)
# ==============================================================================

countQcDomVisible <- reactive(countQcBarVisible())

# (The dominant-genes plot is a tab inside the shared Results tabsetPanel.)

countQcDominantData <- reactive({
  req(countQcDomVisible())
  req(eds$has_assay("counts"))
  counts <- eds$assay("counts")
  # Keep rowData row-aligned to counts (features() is p-value sorted).
  rd <- eds$features()[rownames(counts), , drop = FALSE]
  labels <- if ("gene_name" %in% colnames(rd) && any(nzchar(rd$gene_name))) {
    ifelse(is.na(rd$gene_name) | !nzchar(rd$gene_name),
           rownames(counts), rd$gene_name)
  } else {
    rownames(counts)
  }
  labels <- make.unique(labels)
  col_sums <- colSums(counts)
  col_sums[col_sums == 0] <- 1  # guard against zero-col (empty sample)
  frac <- sweep(counts, 2, col_sums, FUN = "/")
  rownames(frac) <- labels
  df <- as.data.frame(frac)
  df$Gene <- labels
  long <- tidyr::pivot_longer(df, cols = -Gene, names_to = "Sample", values_to = "Fraction")
  long$Sample <- factor(long$Sample, levels = colnames(counts))
  long
})

# Classify a gene name to a family preset / custom / "Other"
classify_family <- function(gene, presets_used, user_regex) {
  for (fam in names(presets_used)) {
    if (grepl(presets_used[[fam]], gene, ignore.case = TRUE)) return(fam)
  }
  if (nzchar(user_regex %||% "") &&
      grepl(user_regex, gene, ignore.case = TRUE)) return("Custom")
  "Other"
}

build_dominant_plot <- function() {
  validate(need(isTRUE(try(countQcDomVisible(), silent = TRUE)),
                "Dominant-genes plot not applicable for this data type."))
  df <- countQcDominantData()
  n_top <- input$countQcDomTopN %||% 3
  n_top <- max(0L, as.integer(n_top))
  colour_mode <- input$countQcDomColourBy %||% "Rank"
  fam_col <- input$countQcColDomFamily %||% "#1965B0"

  # Top-N layer (chosen by expression fraction)
  topN <- if (n_top > 0) {
    df |>
      dplyr::group_by(Sample) |>
      dplyr::slice_max(Fraction, n = n_top, with_ties = FALSE) |>
      dplyr::mutate(rank = factor(dplyr::row_number(), levels = seq_len(n_top))) |>
      dplyr::ungroup()
  } else {
    df[0, ][, c("Gene", "Sample", "Fraction"), drop = FALSE] |>
      dplyr::mutate(rank = factor(character(0)))
  }

  # Family matching: resolve the regex used for both colouring (Family mode)
  # and the always-on family-triangle layer.
  preset_names <- input$countQcDomPresets
  presets_used <- countQcFamilyPresets[preset_names]
  user_regex   <- input$countQcDomPattern %||% ""
  regex_parts  <- c(unname(unlist(presets_used)),
                    if (nzchar(user_regex)) user_regex else NULL)
  hit_regex    <- if (length(regex_parts)) paste(regex_parts, collapse = "|") else NA_character_

  family <- if (!is.na(hit_regex)) {
    df[grepl(hit_regex, df$Gene, ignore.case = TRUE), , drop = FALSE]
  } else {
    df[0, ]
  }

  p <- ggplot(df, aes(x = Sample, y = Fraction)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90")

  # --- Top-N point + label layer with per-mode colouring -------------------
  if (nrow(topN) > 0) {
    if (colour_mode == "Rank") {
      # Interpolate three rank-colour anchors across N ranks so N > 3 works.
      base_cols <- c(input$countQcColDomRank1 %||% "red",
                     input$countQcColDomRank2 %||% "orange",
                     input$countQcColDomRank3 %||% "gold")
      rank_palette <- if (n_top == 1) base_cols[1]
                      else if (n_top == 2) base_cols[c(1, 3)]
                      else grDevices::colorRampPalette(base_cols)(n_top)
      names(rank_palette) <- as.character(seq_len(n_top))

      p <- p +
        geom_point(data = topN, aes(x = Sample, y = Fraction, colour = rank),
                   size = 2, position = position_dodge(width = 0.4),
                   inherit.aes = FALSE) +
        ggrepel::geom_text_repel(
          data = topN, aes(x = Sample, y = Fraction, label = Gene, colour = rank),
          size = 2.8, nudge_y = 0.02, direction = "y",
          segment.size = 0.2, segment.alpha = 0.4,
          max.overlaps = Inf, inherit.aes = FALSE
        ) +
        scale_colour_manual(values = rank_palette, name = "Top-gene rank") +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
    } else {
      # Family mode — classify each top-N gene and colour by membership
      topN$family <- vapply(as.character(topN$Gene),
                            classify_family, FUN.VALUE = character(1),
                            presets_used, user_regex)
      lvls <- c(names(presets_used),
                if (nzchar(user_regex)) "Custom" else NULL,
                "Other")
      topN$family <- factor(topN$family, levels = lvls)
      n_fams <- length(lvls)
      palette_base <- if (n_fams <= 8) {
        RColorBrewer::brewer.pal(max(3, n_fams), "Set1")[seq_len(n_fams)]
      } else {
        scales::hue_pal()(n_fams)
      }
      fam_palette <- setNames(palette_base, lvls)
      if ("Other" %in% lvls) fam_palette["Other"] <- fam_col

      p <- p +
        geom_point(data = topN, aes(x = Sample, y = Fraction, colour = family),
                   size = 2, position = position_dodge(width = 0.4),
                   inherit.aes = FALSE) +
        ggrepel::geom_text_repel(
          data = topN, aes(x = Sample, y = Fraction, label = Gene, colour = family),
          size = 2.8, nudge_y = 0.02, direction = "y",
          segment.size = 0.2, segment.alpha = 0.4,
          max.overlaps = Inf, inherit.aes = FALSE
        ) +
        scale_colour_manual(values = fam_palette, name = "Gene family", drop = FALSE) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
    }
  }

  # --- Always-on family-triangle layer (unmatched-or-matched based on
  # regex-matched data in `family`). Colour: single fam_col so it doesn't
  # collide with the top-N scale above.
  if (nrow(family) > 0) {
    max_labels <- input$countQcDomFamilyLabelsPerSample %||% 3
    max_labels <- max(0L, as.integer(max_labels))
    family_labelled <- if (max_labels > 0) {
      family |>
        dplyr::group_by(Sample) |>
        dplyr::slice_max(Fraction, n = max_labels, with_ties = FALSE) |>
        dplyr::ungroup()
    } else {
      family[0, , drop = FALSE]
    }

    p <- p +
      geom_point(data = family, aes(x = Sample, y = Fraction),
                 inherit.aes = FALSE, colour = fam_col,
                 shape = 17, size = 1.4, alpha = 0.55)
    if (nrow(family_labelled) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = family_labelled, aes(x = Sample, y = Fraction, label = Gene),
        inherit.aes = FALSE, colour = fam_col,
        size = 2.6, segment.alpha = 0.3, max.overlaps = Inf
      )
    }
  }

  p +
    ggprism::theme_prism(base_size = cqT(), axis_text_angle = 90) +
    labs(y = "Fraction of total expression",
         title = "Dominant genes per sample")
}

output$countQcDominant <- renderPlot(
  build_dominant_plot(), width = cqW, height = cqH
)

# ==============================================================================
# 3. Intra-group scatter (all data types; GC/Width require gc/featWidth)
# ==============================================================================

.countQc_seqAnno_or_null <- function() {
  if (length(eds$contrasts()) > 0) eds$de(eds$primary_contrast()) else NULL
}
countQcHasGC <- reactive({
  sa <- .countQc_seqAnno_or_null(); if (is.null(sa)) return(FALSE)
  "gc" %in% colnames(sa) && !all(is.na(sa$gc))
})
countQcHasWidth <- reactive({
  sa <- .countQc_seqAnno_or_null(); if (is.null(sa)) return(FALSE)
  "featWidth" %in% colnames(sa) && !all(is.na(sa$featWidth))
})

# Hide GC/Width radio choices when the columns are absent
observe({
  choices <- "Standard"
  if (isTRUE(try(countQcHasGC(), silent = TRUE)))    choices <- c(choices, "GC content")
  if (isTRUE(try(countQcHasWidth(), silent = TRUE))) choices <- c(choices, "Gene width")
  current <- isolate(input$countQcScatterColourBy)
  selected <- if (!is.null(current) && current %in% choices) current else "Standard"
  updateRadioButtons(session, "countQcScatterColourBy",
                     choices = choices, selected = selected, inline = TRUE)
})

# Resolve highlight tokens -> feature_ids
countQcHighlightIds <- reactive({
  tokens <- c(input$countQcHighlightSelect, parse_gene_text(input$countQcHighlightText))
  tokens <- unique(tokens[nzchar(tokens)])
  if (!length(tokens)) return(character(0))
  sa <- .countQc_seqAnno_or_null()
  if (is.null(sa)) return(character(0))
  hits <- sa$feature_id[sa$feature_id %in% tokens]
  if ("gene_name" %in% colnames(sa)) {
    hits <- unique(c(hits, sa$feature_id[sa$gene_name %in% tokens]))
  }
  hits
})

# Intra-group data + types, given the user's colour-by mode
countQcIntraData <- reactive({
  req(input$countQcIntraFactor)
  req(input$countQcIntraGroup)
  req(input$countQcIntraCountType)
  dataset  <- eds$samples()
  # SE-only escape hatch (presentFlag/isPresent assays); NULL for uploads.
  se       <- if (eds$backend_kind() == "se") eds$backend() else NULL
  dataType <- eds$data_type()
  # Proteomics: use the currently selected contrast's DE table
  seqAnno <- if (eds$omics_class() == "proteomics") {
    contr <- input$contrastSelected %||% eds$contrasts()[[1]]
    eds$de(contr)
  } else if (length(eds$contrasts()) > 0) {
    eds$de(eds$primary_contrast())
  } else NULL
  req(!is.null(seqAnno))

  req(input$countQcIntraFactor %in% colnames(dataset))
  fc <- input$countQcIntraFactor
  grp <- input$countQcIntraGroup
  idx <- which(as.character(dataset[[fc]]) == grp)
  validate(need(length(idx) >= 2, "Need at least 2 samples in the selected group."))

  cts <- eds$assay(input$countQcIntraCountType)
  req(!is.null(cts))
  # Apply background offset if the count type looks linear
  offset <- if (input$countQcIntraCountType %in% c("normalised", "counts", "fpkm", "tpm")) {
    eds$background_expression()
  } else 10
  signal <- cts + offset

  # Align signal to seqAnno by feature_id (with isoform fallback). Keeps rd
  # and signal in matching row order so the types vector and plotted x/y
  # stay aligned.
  aln <- align_signal_to_seqanno(signal, seqAnno)
  signal <- aln$signal; rd <- aln$rd
  req(nrow(rd) > 0)

  samples <- colnames(signal)[idx]
  # geometric mean over samples in the group
  log_sig <- log2(pmax(signal[, samples, drop = FALSE], .Machine$double.eps))
  group_geomean <- 2^rowMeans(log_sig)

  # Build types vector per colour-by mode (via shared helper)
  mode <- input$countQcScatterColourBy %||% "Standard"
  present_mask <- if (dataType == "RNASeq" && "usedInTest" %in% colnames(rd)) {
    rd$usedInTest
  } else if (!is.null(se) && "presentFlag" %in% assayNames(se)) {
    pf <- rowMeans(assay(se, "presentFlag")) >= 0.5
    if (!is.null(names(pf))) pf <- pf[match(rd$feature_id, names(pf))]
    pf[is.na(pf)] <- FALSE
    pf
  } else if (!is.null(se) && "isPresent" %in% assayNames(se)) {
    ip <- rowMeans(assay(se, "isPresent")) >= 0.5
    if (!is.null(names(ip))) ip <- ip[match(rd$feature_id, names(ip))]
    ip[is.na(ip)] <- FALSE
    ip
  } else rep(TRUE, nrow(rd))

  types_info <- build_qc_scatter_types(
    rd, mode, present_mask,
    p_thresh  = input$countQcScatterPval   %||% 0.05,
    lfc_thresh = input$countQcScatterLFC   %||% 1,
    gc_lo = input$countQcScatterGCLow  %||% 0.4,
    gc_hi = input$countQcScatterGCHigh %||% 0.6,
    w_lo  = input$countQcScatterWidthLow  %||% 500,
    w_hi  = input$countQcScatterWidthHigh %||% 5000,
    std_cols = list(
      Absent       = input$countQcColAbsent  %||% "grey70",
      Present      = input$countQcColPresent %||% "black",
      Significants = input$countQcColSig     %||% "#A6CEE3"
    ),
    gc_cols = list(
      low  = input$countQcColGCLow  %||% "#1F78B4",
      mid  = input$countQcColGCMid  %||% "black",
      high = input$countQcColGCHigh %||% "#FF7F00"
    ),
    width_cols = list(
      low  = input$countQcColWidthLow  %||% "#B15928",
      mid  = input$countQcColWidthMid  %||% "black",
      high = input$countQcColWidthHigh %||% "#33A02C"
    )
  )

  list(
    signal = signal,
    samples = samples,
    group = grp,
    group_geomean = group_geomean,
    types_vec = types_info$types,
    types_colours = types_info$colours,
    draw_order = types_info$draw_order,
    feature_id = rd$feature_id,
    label = if ("gene_name" %in% colnames(rd) && any(nzchar(rd$gene_name))) {
      ifelse(is.na(rd$gene_name) | !nzchar(rd$gene_name),
             rd$feature_id, rd$gene_name)
    } else rd$feature_id
  )
})

build_intra_scatter_plot <- function() {
  d <- countQcIntraData()
  hl_ids <- countQcHighlightIds()

  guide_col    <- input$countQcColGuide        %||% "blue"
  guide_fold   <- input$countQcGuideFold       %||% 2
  show_id      <- input$countQcShowIdentity    %||% TRUE
  show_fold    <- input$countQcShowFoldBands   %||% TRUE

  panel_plots <- lapply(d$samples, function(s) {
    p <- make_qc_scatter(
      x = d$group_geomean,
      y = d$signal[, s],
      types_vec = d$types_vec,
      types_colours = d$types_colours,
      xlab = paste("Avg of", d$group),
      ylab = s,
      title = s,
      point_size = input$countQcScatterPointSize %||% 0.6,
      alpha = input$countQcScatterAlpha %||% 0.4,
      base_size = max(cqT() - 2, 6),
      draw_order = d$draw_order,
      guide_colour    = guide_col,
      guide_fold      = guide_fold,
      show_identity   = show_id,
      show_fold_bands = show_fold
    ) + theme(legend.position = "none")

    if (length(hl_ids)) {
      hl <- data.frame(
        feature_id = d$feature_id,
        x = d$group_geomean, y = d$signal[, s],
        label = d$label, stringsAsFactors = FALSE
      )
      hl <- hl[hl$feature_id %in% hl_ids, , drop = FALSE]
      if (nrow(hl)) {
        col <- input$countQcColHighlight %||% "goldenrod3"
        p <- p + geom_point(
          data = hl, aes(x = x, y = y), inherit.aes = FALSE,
          colour = col, size = input$countQcHighlightSize %||% 3
        )
        if (isTRUE(input$countQcHighlightLabel)) {
          p <- p + ggrepel::geom_text_repel(
            data = hl, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, colour = col,
            size = input$countQcHighlightLabelSize %||% 3.5,
            min.segment.length = 0, max.overlaps = Inf,
            segment.alpha = 0.5
          )
        }
      }
    }
    p
  })

  # Keep the legend from the first panel so users can still see it
  legend_plot <- make_qc_scatter(
    x = d$group_geomean, y = d$signal[, d$samples[[1]]],
    types_vec = d$types_vec,
    types_colours = d$types_colours,
    xlab = "", ylab = "", title = NULL,
    point_size = input$countQcScatterPointSize %||% 0.6,
    alpha = input$countQcScatterAlpha %||% 0.4,
    base_size = cqT(),
    draw_order = d$draw_order,
    guide_colour    = guide_col,
    guide_fold      = guide_fold,
    show_identity   = show_id,
    show_fold_bands = show_fold
  )
  legend_grob <- cowplot::get_legend(legend_plot)

  ncol <- min(length(panel_plots), input$countQcIntraNcol %||% 4)
  grid <- cowplot::plot_grid(plotlist = panel_plots, ncol = ncol, align = "hv")
  cowplot::plot_grid(grid, legend_grob, ncol = 2, rel_widths = c(1, 0.15))
}

intraPlotHeight <- function() {
  d <- try(countQcIntraData(), silent = TRUE)
  if (inherits(d, "try-error") || is.null(d)) return(cqH())
  ncol <- min(length(d$samples), input$countQcIntraNcol %||% 4)
  nrow <- max(1, ceiling(length(d$samples) / ncol))
  max(cqH(), nrow * 320)
}
intraPlotWidth <- function() {
  d <- try(countQcIntraData(), silent = TRUE)
  if (inherits(d, "try-error") || is.null(d)) return(cqW())
  ncol <- min(length(d$samples), input$countQcIntraNcol %||% 4)
  max(cqW(), ncol * 300)
}

output$countQcIntraScatter <- renderPlot(
  build_intra_scatter_plot(),
  width  = intraPlotWidth,
  height = intraPlotHeight
)

# ==============================================================================
# Downloads
# ==============================================================================

countQcDownloadFile <- function(plot_reactive, suffix, width_in = 8, height_in = 5) {
  downloadHandler(
    filename = function() {
      ext <- tolower(input$countQcDlFormat %||% "pdf")
      paste0(input$countQcDlName %||% "Count_QC", "_", suffix, ".", ext)
    },
    content = function(file) {
      fmt <- tolower(input$countQcDlFormat %||% "pdf")
      if (fmt == "png") {
        ggsave(file, plot = plot_reactive(), device = "png",
               width = width_in, height = height_in,
               dpi = as.numeric(input$countQcDlDpi %||% 300), bg = "white")
      } else {
        plt <- plot_reactive()
        if (fmt == "pdf") plt <- prepare_plot_for_pdf(plt)
        ggsave(file, plot = plt, device = fmt,
               width = width_in, height = height_in, bg = "white")
      }
    }
  )
}

output$dlPlot_countQcTotalReads   <- countQcDownloadFile(
  reactive(build_total_reads_plot()), "total_reads",
  width_in = cqW() / 96, height_in = cqH() / 96)
output$dlPlot_countQcFeaturesAbove <- countQcDownloadFile(
  reactive(build_features_above_plot()), "features_above",
  width_in = cqW() / 96, height_in = cqH() / 96)
output$dlPlot_countQcDominant     <- countQcDownloadFile(
  reactive(build_dominant_plot()), "dominant_genes",
  width_in = cqW() / 96, height_in = cqH() / 96)
output$dlPlot_countQcIntra        <- countQcDownloadFile(
  reactive(build_intra_scatter_plot()), "intra_group",
  width_in = intraPlotWidth() / 96, height_in = intraPlotHeight() / 96)
