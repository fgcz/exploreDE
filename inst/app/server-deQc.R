## DE QC tab — RNASeq + proteomics
## Plots:
##   - Between-group scatter (Standard / GC / Width colouring)
##   - P-value histogram
##   - GC + Length LOESS bias (RNASeq only; seqAnno needs gc / featWidth)

req(!is.null(inputDataReactive()$dataType))
deQcVisible <- reactive({
  inputDataReactive()$dataType %in% c("RNASeq", "proteomics")
})

# ---- Hide inapplicable plot tabs per dataType --------------------------------
# GC/Length LOESS bias plots require log2Ratio + gc/featWidth, which
# proteomics seqAnnoList doesn't carry. Scatter colour-by radio also restricts
# to Standard when gc/featWidth columns are absent.
observe({
  req(inputDataReactive()$dataType)
  dataType <- inputDataReactive()$dataType
  is_rna <- dataType == "RNASeq"
  if (is_rna) {
    showTab("deQcPlotsTabs", "GC Bias",     session = session)
    showTab("deQcPlotsTabs", "Length Bias", session = session)
  } else {
    hideTab("deQcPlotsTabs", "GC Bias",     session = session)
    hideTab("deQcPlotsTabs", "Length Bias", session = session)
  }
  # Adjust scatter colour-by choices based on what's available
  choices <- "Standard"
  if (is_rna) {
    if (isTRUE(try(deQcHasGC(),    silent = TRUE))) choices <- c(choices, "GC content")
    if (isTRUE(try(deQcHasWidth(), silent = TRUE))) choices <- c(choices, "Gene width")
  }
  current <- isolate(input$deQcScatterColourBy)
  selected <- if (!is.null(current) && current %in% choices) current else "Standard"
  updateRadioButtons(session, "deQcScatterColourBy",
                     choices = choices, selected = selected, inline = TRUE)
})

# ---- Sync Settings tabset to the Results tabset -----------------------------
# Clicking a plot tab jumps Settings to the matching sub-tab so the controls
# for that plot are right there.
observeEvent(input$deQcPlotsTabs, ignoreInit = TRUE, {
  mapping <- c(
    "Between-group Scatter" = "Between-group",
    "P-value Distribution"  = "P-value hist",
    "GC Bias"               = "Bias",
    "Length Bias"           = "Bias"
  )
  target <- mapping[[input$deQcPlotsTabs]]
  if (!is.null(target)) {
    updateTabsetPanel(session, "deQcSettingsTabs", selected = target)
  }
})

# ---- Populate highlight selectize with available feature_ids -----------------
observe({
  req(deQcVisible())
  updateSelectizeInput(
    session = session,
    inputId = "deQcHighlightSelect",
    choices = inputDataReactive()$genes,
    selected = character(0),
    server = TRUE
  )
})

# ---- Proteomics-only selectors -----------------------------------------------
# Contrast (reuses the app-wide input$contrastSelected); count-type + group
# pair for building the scatter. Choices embedded in the renderUI so the
# selectInputs materialise with sensible defaults the moment the UI is drawn
# (rather than relying on update* observers that race the plot render).
output$deQcProteomicsSelectorsUI <- renderUI({
  if (!isTRUE(inputDataReactive()$dataType == "proteomics")) return(NULL)
  ch_counts <- names(inputDataReactive()$countList) %||% character(0)
  default_count <- if ("Normalised" %in% ch_counts) "Normalised" else ch_counts[[1]]

  contr <- input$contrastSelected
  levels_avail <- if (!is.null(contr) && nzchar(contr)) {
    v <- strsplit(contr, "_vs_")[[1]]
    if (length(v) == 2) v else character(0)
  } else character(0)
  if (!length(levels_avail)) {
    dataset <- inputDataReactive()$dataset
    f1 <- inputDataReactive()$factors[[1]]
    if (!is.null(f1) && f1 %in% colnames(dataset)) {
      levels_avail <- unique(as.character(dataset[[f1]]))
    }
  }
  default_sample <- if (length(levels_avail)) levels_avail[[1]] else NULL
  default_ref    <- if (length(levels_avail) >= 2) levels_avail[[2]] else default_sample

  tagList(
    selectInput("deQcProtCountType", "Counts to use",
                choices = ch_counts, selected = default_count),
    splitLayout(
      selectInput("deQcProtSampleGroup", "Sample group (y)",
                  choices = levels_avail, selected = default_sample, width = "85%"),
      selectInput("deQcProtRefGroup", "Reference group (x)",
                  choices = levels_avail, selected = default_ref, width = "85%"),
      cellArgs = list(style = "overflow: visible;")
    ),
    hr()
  )
})
outputOptions(output, "deQcProteomicsSelectorsUI", suspendWhenHidden = FALSE)

# ---- Resolve highlight selections -> feature_ids -----------------------------
deQcHighlightIds <- reactive({
  tokens <- c(input$deQcHighlightSelect, parse_gene_text(input$deQcHighlightText))
  tokens <- unique(tokens[nzchar(tokens)])
  if (!length(tokens)) return(character(0))
  sa <- inputDataReactive()$seqAnno
  # Match on feature_id first, then gene_name
  hits <- sa$feature_id[sa$feature_id %in% tokens]
  if ("gene_name" %in% colnames(sa)) {
    hits <- unique(c(hits, sa$feature_id[sa$gene_name %in% tokens]))
  }
  hits
})

# ---- Reactive: between-group scatter data ------------------------------------
# Returns ref/sample mean expression vectors + a per-colour-mode types/colours
# descriptor. Works for RNASeq (param-driven groups) and proteomics
# (user-picked groups + contrast's seqAnnoList).
deQcScatterData <- reactive({
  req(deQcVisible())
  dataType <- inputDataReactive()$dataType
  countList <- inputDataReactive()$countList
  dataset <- inputDataReactive()$dataset
  param <- inputDataReactive()$param
  req(!is.null(countList))

  if (dataType == "RNASeq") {
    seqAnno <- inputDataReactive()$seqAnno
    req(!is.null(seqAnno))
    grouping <- param$groupingName
    req(!is.null(param$sampleGroup), !is.null(param$refGroup),
        grouping %in% colnames(dataset))
    signal <- countList[["Normalised"]] + (param$backgroundExpression %||% 10)
    rd <- seqAnno
    sample_label <- param$sampleGroup
    ref_label    <- param$refGroup
    sample_idx <- which(dataset[[grouping]] == param$sampleGroup)
    ref_idx    <- which(dataset[[grouping]] == param$refGroup)
  } else {
    # proteomics
    contr <- input$contrastSelected
    req(!is.null(contr), nzchar(contr))
    sa_list <- inputDataReactive()$seqAnnoList
    req(!is.null(sa_list), contr %in% names(sa_list))
    rd <- sa_list[[contr]]
    signal <- countList[[input$deQcProtCountType %||% names(countList)[[1]]]]
    req(!is.null(signal))
    grouping <- inputDataReactive()$factors[[1]]
    # prefer the G_ column used by proteomics tabs if it exists
    if ("G_" %in% colnames(dataset)) grouping <- "G_"
    req(grouping %in% colnames(dataset))
    sample_label <- input$deQcProtSampleGroup %||% NA
    ref_label    <- input$deQcProtRefGroup    %||% NA
    validate(need(!is.na(sample_label) & !is.na(ref_label) &
                    sample_label != ref_label,
                  "Pick two distinct groups (sample vs reference)."))
    sample_idx <- which(as.character(dataset[[grouping]]) == sample_label)
    ref_idx    <- which(as.character(dataset[[grouping]]) == ref_label)
  }
  req(length(sample_idx) > 0, length(ref_idx) > 0)

  # Align signal to rd by feature_id. For isoform-level data, init replaces
  # special chars in rownames(signal) ([:/()-] + space → "."), so a direct
  # match can leave NAs; the helper applies the same substitution as a
  # fallback and drops any features that can't be matched either way. This
  # keeps rd and signal in the same row order — critical for the types
  # vector (from rd) and x/y (from signal) to line up.
  aln <- align_signal_to_seqanno(signal, rd)
  signal <- aln$signal; rd <- aln$rd
  req(nrow(rd) > 0)

  sample_mean <- rowMeans(signal[, sample_idx, drop = FALSE])
  ref_mean    <- rowMeans(signal[, ref_idx, drop = FALSE])

  # Present mask: RNASeq = usedInTest; proteomics has no such flag
  present_mask <- if (dataType == "RNASeq" && "usedInTest" %in% colnames(rd)) {
    rd$usedInTest
  } else rep(TRUE, nrow(rd))

  mode <- input$deQcScatterColourBy %||% "Standard"
  # Auto-fallback: if mode needs gc/featWidth but column is missing
  if (mode == "GC content" && (!"gc" %in% colnames(rd) || all(is.na(rd$gc)))) mode <- "Standard"
  if (mode == "Gene width" && (!"featWidth" %in% colnames(rd) || all(is.na(rd$featWidth)))) mode <- "Standard"

  types_info <- build_qc_scatter_types(
    rd, mode, present_mask,
    p_thresh   = input$deQcScatterPval  %||% 0.05,
    lfc_thresh = input$deQcScatterLFC   %||% 1,
    gc_lo = input$deQcScatterGCLow  %||% 0.4,
    gc_hi = input$deQcScatterGCHigh %||% 0.6,
    w_lo  = input$deQcScatterWidthLow  %||% 500,
    w_hi  = input$deQcScatterWidthHigh %||% 5000,
    std_cols = list(
      Absent       = input$deQcColAbsent  %||% "grey70",
      Present      = input$deQcColPresent %||% "black",
      Significants = input$deQcColSig     %||% "#A6CEE3"
    ),
    gc_cols = list(
      low  = input$deQcColGCLow  %||% "#1F78B4",
      mid  = input$deQcColGCMid  %||% "black",
      high = input$deQcColGCHigh %||% "#FF7F00"
    ),
    width_cols = list(
      low  = input$deQcColWidthLow  %||% "#B15928",
      mid  = input$deQcColWidthMid  %||% "black",
      high = input$deQcColWidthHigh %||% "#33A02C"
    )
  )

  labels <- if ("gene_name" %in% colnames(rd) && any(nzchar(rd$gene_name))) {
    ifelse(is.na(rd$gene_name) | !nzchar(rd$gene_name),
           rd$feature_id, rd$gene_name)
  } else rd$feature_id

  # Table for brushedPoints() — log_x / log_y are needed because brush coords
  # on scale_x_log10() / scale_y_log10() axes come back in log-transformed space.
  plotDf <- data.frame(
    feature_id = rd$feature_id,
    gene_name  = labels,
    x = ref_mean, y = sample_mean,
    log_x = suppressWarnings(log10(ref_mean)),
    log_y = suppressWarnings(log10(sample_mean)),
    types = types_info$types,
    log2Ratio = if ("log2Ratio" %in% colnames(rd)) rd$log2Ratio else NA_real_,
    pValue    = if ("pValue"    %in% colnames(rd)) rd$pValue    else NA_real_,
    fdr       = if ("fdr"       %in% colnames(rd)) rd$fdr       else NA_real_,
    stringsAsFactors = FALSE
  )

  list(
    x = ref_mean, y = sample_mean,
    types_vec = types_info$types,
    types_colours = types_info$colours,
    draw_order = types_info$draw_order,
    xlab = ref_label, ylab = sample_label,
    feature_id = rd$feature_id,
    label = labels,
    plotDf = plotDf
  )
})

# ---- Output: design / comparison label ---------------------------------------
output$deQcScatterDesign <- renderText({
  req(deQcVisible())
  dataType <- inputDataReactive()$dataType
  if (dataType == "RNASeq") {
    inputDataReactive()$design %||% ""
  } else {
    input$contrastSelected %||% ""
  }
})

# ---- Output: between-group scatter -------------------------------------------
# Safe input accessors with defaults — non-initial-tab inputs can be NULL
# before the user activates their tab. Keeps plots rendering eagerly.
deQcW <- function() { w <- input$deQcFigWidth;  if (is.null(w)) 700 else w }
deQcH <- function() { h <- input$deQcFigHeight; if (is.null(h)) 500 else h }
deQcT <- function() { t <- input$deQcTextSize;  if (is.null(t)) 12  else t }

# Overlay highlight points + labels onto a ggplot.
# `lookup` is a data frame with columns: feature_id, x, y, label
add_highlight_layer <- function(p, lookup) {
  ids <- deQcHighlightIds()
  if (!length(ids)) return(p)
  hl <- lookup[lookup$feature_id %in% ids, , drop = FALSE]
  if (!nrow(hl)) return(p)
  col <- input$deQcColHighlight %||% "goldenrod3"
  sz  <- input$deQcHighlightSize %||% 3
  p <- p + ggplot2::geom_point(
    data = hl, ggplot2::aes(x = x, y = y),
    inherit.aes = FALSE, colour = col, size = sz
  )
  if (isTRUE(input$deQcHighlightLabel)) {
    p <- p + ggrepel::geom_text_repel(
      data = hl, ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE, colour = col,
      size = input$deQcHighlightLabelSize %||% 4,
      min.segment.length = 0, max.overlaps = Inf,
      segment.alpha = 0.5
    )
  }
  p
}

build_between_group_plot <- function() {
  d <- deQcScatterData()
  p <- make_qc_scatter(
    x = d$x, y = d$y,
    types_vec = d$types_vec,
    types_colours = d$types_colours,
    xlab = d$xlab, ylab = d$ylab,
    title = "Between-group expression comparison",
    point_size = input$deQcScatterPointSize %||% 1,
    alpha = input$deQcScatterAlpha %||% 0.6,
    base_size = deQcT(),
    draw_order = d$draw_order,
    guide_colour    = input$deQcColGuide       %||% "blue",
    guide_fold      = input$deQcGuideFold      %||% 2,
    show_identity   = input$deQcShowIdentity   %||% TRUE,
    show_fold_bands = input$deQcShowFoldBands  %||% TRUE
  )
  lookup <- data.frame(
    feature_id = d$feature_id, x = d$x, y = d$y,
    label = d$label, stringsAsFactors = FALSE
  )
  add_highlight_layer(p, lookup)
}

output$betweenGroupScatter <- renderPlot(
  build_between_group_plot(),
  width = deQcW, height = deQcH
)

# Brush-driven table of selected features.
# Manual filter (rather than brushedPoints) because the plot uses
# scale_x_log10()/scale_y_log10(); the brush box dimensions can come back in
# either log-transformed or raw space depending on the Shiny/ggplot build, so
# we evaluate both and keep whichever yields hits.
output$deQcScatterBrushTable <- DT::renderDataTable({
  brush <- input$deQcScatterBrush
  df <- deQcScatterData()$plotDf
  empty <- df[0, , drop = FALSE]
  make_dt <- function(d) {
    cols_show <- c("feature_id", "gene_name", "x", "y",
                   "log2Ratio", "pValue", "fdr", "types")
    cols_show <- intersect(cols_show, colnames(d))
    d <- d[, cols_show, drop = FALSE]
    num_cols <- vapply(d, is.numeric, logical(1))
    d[num_cols] <- lapply(d[num_cols], function(v) signif(v, 4))
    DT::datatable(d, rownames = FALSE,
                  options = list(pageLength = 8, scrollX = TRUE,
                                 dom = "tip", order = list()))
  }
  if (is.null(brush) || !nrow(df)) return(make_dt(empty))

  keep_log <- !is.na(df$log_x) & !is.na(df$log_y) &
              df$log_x >= brush$xmin & df$log_x <= brush$xmax &
              df$log_y >= brush$ymin & df$log_y <= brush$ymax
  keep_raw <- !is.na(df$x) & !is.na(df$y) &
              df$x >= brush$xmin & df$x <= brush$xmax &
              df$y >= brush$ymin & df$y <= brush$ymax
  picked <- if (sum(keep_log, na.rm = TRUE) > 0) {
    df[which(keep_log), , drop = FALSE]
  } else {
    df[which(keep_raw), , drop = FALSE]
  }
  make_dt(picked)
})

# ---- Output: p-value histogram -----------------------------------------------
# Build the p-value histogram data frame for whichever dataType applies.
deQcPhistDf <- reactive({
  req(deQcVisible())
  pcol <- input$deQcPhistType %||% "pValue"
  dataType <- inputDataReactive()$dataType
  if (dataType == "RNASeq") {
    sa <- inputDataReactive()$seqAnno
    req(pcol %in% colnames(sa))
    df <- tibble::tibble(
      pval = sa[[pcol]],
      Category = ifelse(sa$usedInTest, "Used in test", "Not expressed")
    )
  } else {
    contr <- input$contrastSelected
    req(!is.null(contr), nzchar(contr))
    sa <- inputDataReactive()$seqAnnoList[[contr]]
    req(!is.null(sa), pcol %in% colnames(sa))
    # modelName-aware categorisation (Linear vs Imputed) — falls back
    # to a single category when modelName isn't present.
    Category <- if ("modelName" %in% colnames(sa)) {
      ifelse(grepl("^Linear", sa$modelName), "Linear model", "Imputed")
    } else rep("All features", nrow(sa))
    df <- tibble::tibble(pval = sa[[pcol]], Category = Category)
  }
  df[!is.na(df$pval), , drop = FALSE]
})

output$pvalHistogram <- renderPlot({
  df <- deQcPhistDf()
  req(nrow(df) > 0)
  pcol <- input$deQcPhistType %||% "pValue"
  nbins <- input$deQcPhistBins %||% 200
  # expected-uniform line only meaningful on the "used"/"linear" category
  n_used <- sum(df$Category %in% c("Used in test", "Linear model", "All features"))

  col_used    <- input$deQcColPhistUsed    %||% "#1965B0"
  col_notused <- input$deQcColPhistNotUsed %||% "#FF7B00"

  ggplot(df, aes(x = pval, fill = Category)) +
    geom_histogram(bins = nbins, boundary = 0, colour = "white", linewidth = 0.15) +
    geom_hline(yintercept = n_used / nbins, linetype = "dashed",
               colour = "black", linewidth = 0.6) +
    scale_fill_manual(values = c(
      "Used in test"  = col_used,
      "Not expressed" = col_notused,
      "Linear model"  = col_used,
      "Imputed"       = col_notused,
      "All features"  = col_used
    )) +
    ggprism::theme_prism(base_size = deQcT()) +
    labs(
      x = if (pcol == "pValue") "P-value" else "FDR",
      y = "Count",
      title = "P-value distribution",
      subtitle = "Dashed line = expected uniform count under the null"
    )
},
width  = deQcW,
height = deQcH)

# ---- Bias plots (GC, length) -------------------------------------------------
# Which bias plots are supported depends on seqAnno columns.
deQcHasGC <- reactive({
  if (!isTRUE(inputDataReactive()$dataType == "RNASeq")) return(FALSE)
  sa <- inputDataReactive()$seqAnno
  "gc" %in% colnames(sa) && !all(is.na(sa$gc))
})
deQcHasWidth <- reactive({
  if (!isTRUE(inputDataReactive()$dataType == "RNASeq")) return(FALSE)
  sa <- inputDataReactive()$seqAnno
  "featWidth" %in% colnames(sa) && !all(is.na(sa$featWidth))
})

# Bias tabs are always present in the UI; each plot uses validate/need to
# show a friendly message when its required column is missing.

deQcBiasData <- reactive({
  validate(need(isTRUE(inputDataReactive()$dataType == "RNASeq"),
                "LOESS bias plots are only available for RNASeq DE results."))
  sa <- inputDataReactive()$seqAnno
  req("log2Ratio" %in% colnames(sa), "pValue" %in% colnames(sa))
  df <- sa[which(sa$usedInTest & !is.na(sa$pValue) & !is.na(sa$log2Ratio)), , drop = FALSE]
  req(nrow(df) > 0)
  df$types <- ifelse(
    abs(df$log2Ratio) >= (input$deQcBiasLFC %||% 1) &
      df$pValue <= (input$deQcBiasPval %||% 0.05),
    "Significants", "Present"
  )
  df$types <- factor(df$types, levels = c("Present", "Significants"))
  df
})

build_bias_plot <- function(df, xvar, xlab, title) {
  bias_cols <- c(
    "Present"      = input$deQcColBiasPresent %||% "black",
    "Significants" = input$deQcColBiasSig     %||% "#A6CEE3"
  )
  p <- ggplot(df, aes(x = .data[[xvar]], y = log2Ratio, colour = types)) +
    geom_point(size = input$deQcBiasPointSize %||% 0.8,
               alpha = input$deQcBiasAlpha %||% 0.3) +
    scale_color_manual(values = bias_cols) +
    geom_hline(yintercept = 0, colour = "blue") +
    ggprism::theme_prism(base_size = deQcT()) +
    theme(legend.title = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = xlab, y = "Log2 Fold Change", title = title)
  if (isTRUE(input$deQcBiasLoess)) {
    p <- p + geom_smooth(aes(group = 1), method = "loess", se = TRUE,
                         colour = "#DC050C", fill = "#DC050C",
                         alpha = 0.15, linewidth = 0.7, formula = y ~ x)
  }
  p
}

bias_labels <- function(df) {
  if ("gene_name" %in% colnames(df) && any(nzchar(df$gene_name))) {
    ifelse(is.na(df$gene_name) | !nzchar(df$gene_name), df$feature_id, df$gene_name)
  } else {
    df$feature_id
  }
}

build_gc_bias_plot <- function() {
  validate(need(isTRUE(try(deQcHasGC(), silent = TRUE)),
                "GC bias plot not available: the 'gc' column is missing from seqAnno."))
  df <- deQcBiasData()
  req("gc" %in% colnames(df))
  df <- df[!is.na(df$gc), , drop = FALSE]
  p <- build_bias_plot(df, xvar = "gc", xlab = "GC Content",
                       title = "GC content vs fold change")
  lookup <- data.frame(
    feature_id = df$feature_id, x = df$gc, y = df$log2Ratio,
    label = bias_labels(df), stringsAsFactors = FALSE
  )
  add_highlight_layer(p, lookup)
}

build_length_bias_plot <- function() {
  validate(need(isTRUE(try(deQcHasWidth(), silent = TRUE)),
                "Length bias plot not available: the 'featWidth' column is missing from seqAnno."))
  df <- deQcBiasData()
  req("featWidth" %in% colnames(df))
  df <- df[!is.na(df$featWidth) & df$featWidth > 0, , drop = FALSE]
  df$log10featWidth <- log10(df$featWidth)
  p <- build_bias_plot(df, xvar = "log10featWidth",
                       xlab = "Log10 Feature Width (bp)",
                       title = "Gene length vs fold change")
  lookup <- data.frame(
    feature_id = df$feature_id, x = df$log10featWidth, y = df$log2Ratio,
    label = bias_labels(df), stringsAsFactors = FALSE
  )
  add_highlight_layer(p, lookup)
}

output$gcBiasPlot <- renderPlot(
  build_gc_bias_plot(), width = deQcW, height = deQcH
)

output$lengthBiasPlot <- renderPlot(
  build_length_bias_plot(), width = deQcW, height = deQcH
)

# Brush tables for the bias plots. Both axes are linear so we use the raw
# columns (gc / log10featWidth / log2Ratio) directly.
deQcBiasBrushTable <- function(xvar) {
  brush_obj <- if (xvar == "gc") input$deQcGCBrush else input$deQcLengthBrush
  if (is.null(brush_obj)) return(data.frame()[0, ])
  df <- deQcBiasData()
  if (!nrow(df)) return(df[0, ])
  if (xvar == "log10featWidth") {
    df <- df[!is.na(df$featWidth) & df$featWidth > 0, , drop = FALSE]
    df$log10featWidth <- log10(df$featWidth)
  } else {
    df <- df[!is.na(df$gc), , drop = FALSE]
  }
  picked <- brushedPoints(df, brush_obj, xvar = xvar, yvar = "log2Ratio")
  cols_show <- c("feature_id", "gene_name", xvar, "log2Ratio", "pValue", "fdr", "types")
  cols_show <- intersect(cols_show, colnames(picked))
  picked <- picked[, cols_show, drop = FALSE]
  num_cols <- vapply(picked, is.numeric, logical(1))
  picked[num_cols] <- lapply(picked[num_cols], function(v) signif(v, 4))
  picked
}

output$deQcGCBrushTable <- DT::renderDataTable({
  DT::datatable(
    deQcBiasBrushTable("gc"), rownames = FALSE,
    options = list(pageLength = 8, scrollX = TRUE, dom = "tip", order = list())
  )
})

output$deQcLengthBrushTable <- DT::renderDataTable({
  DT::datatable(
    deQcBiasBrushTable("log10featWidth"), rownames = FALSE,
    options = list(pageLength = 8, scrollX = TRUE, dom = "tip", order = list())
  )
})

# ---- Downloads ---------------------------------------------------------------
deQcDownloadFile <- function(plot_reactive, suffix) {
  downloadHandler(
    filename = function() {
      ext <- tolower(input$deQcDlFormat)
      paste0(input$deQcDlName, "_", suffix, ".", ext)
    },
    content = function(file) {
      w_in <- input$deQcFigWidth / 96
      h_in <- input$deQcFigHeight / 96
      fmt  <- tolower(input$deQcDlFormat)
      if (fmt == "png") {
        ggsave(file, plot = plot_reactive(), device = "png",
               width = w_in, height = h_in,
               dpi = as.numeric(input$deQcDlDpi), bg = "white")
      } else {
        ggsave(file, plot = plot_reactive(), device = fmt,
               width = w_in, height = h_in, bg = "white")
      }
    }
  )
}

pvalHistogramPlot <- reactive({
  df <- deQcPhistDf()
  pcol <- input$deQcPhistType %||% "pValue"
  nbins <- input$deQcPhistBins %||% 200
  n_used <- sum(df$Category %in% c("Used in test", "Linear model", "All features"))
  col_used    <- input$deQcColPhistUsed    %||% "#1965B0"
  col_notused <- input$deQcColPhistNotUsed %||% "#FF7B00"
  ggplot(df, aes(x = pval, fill = Category)) +
    geom_histogram(bins = nbins, boundary = 0, colour = "white", linewidth = 0.15) +
    geom_hline(yintercept = n_used / nbins, linetype = "dashed",
               colour = "black", linewidth = 0.6) +
    scale_fill_manual(values = c(
      "Used in test"  = col_used, "Not expressed" = col_notused,
      "Linear model"  = col_used, "Imputed"       = col_notused,
      "All features"  = col_used
    )) +
    ggprism::theme_prism(base_size = deQcT()) +
    labs(
      x = if (pcol == "pValue") "P-value" else "FDR",
      y = "Count",
      title = "P-value distribution",
      subtitle = "Dashed line = expected uniform count under the null"
    )
})

output$dlPlot_betweenGroupScatter <- deQcDownloadFile(reactive(build_between_group_plot()), "between_group_scatter")
output$dlPlot_pvalHistogram       <- deQcDownloadFile(pvalHistogramPlot,                    "pvalue_histogram")
output$dlPlot_gcBias              <- deQcDownloadFile(reactive(build_gc_bias_plot()),      "gc_bias")
output$dlPlot_lengthBias          <- deQcDownloadFile(reactive(build_length_bias_plot()),  "length_bias")
