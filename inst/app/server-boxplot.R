dataSummary <- function(data, varname, groupnames) {
  summary_func <- function(x, col) {
    c(
      mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE),
      se = sd(x[[col]], na.rm = TRUE) / sqrt(sum(!is.na(x[[col]])))
    )
  }
  data_sum <- ddply(
    data, groupnames,
    .fun = summary_func, varname
  )
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

eds <- inputDataReactive()$.eds
req(eds)

# Create reactive for design to use across all outputs and download handlers
designReactive <- reactive({
  if (eds$data_type() == "RNASeq") {
    return(eds$contrast_label(eds$primary_contrast()))
  } else if (eds$omics_class() == "proteomics") {
    if (isTruthy(input$contrastSelected)) {
      return(input$contrastSelected)
    } else {
      return("No contrast selected")
    }
  } else if (eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
    return(eds$data_type())
  }
  return("Unknown")
})

# Update a bunch of input options ----
output$boxplotStatic <- renderPlot({
  req(figuresDataReactive$boxplot)
  figuresDataReactive$boxplot
})

.assay_keys    <- eds$assay_names()
.count_choices <- setNames(.assay_keys, vapply(.assay_keys, eds$assay_label, character(1)))
.default_count <- if (eds$omics_class() == "proteomics") {
  if ("intensities" %in% .assay_keys) "intensities" else .assay_keys[1]
} else {
  if ("normalised_log2" %in% .assay_keys) "normalised_log2" else .assay_keys[1]
}
updateSelectInput(
  session = session,
  inputId = "boxplotCounts",
  choices = .count_choices,
  selected = .default_count
)

updateSelectizeInput(
  session = session,
  inputId = "fb_select_boxplot",
  choices = eds$feature_ids(),
  selected = "",
  server = TRUE
)

# Update feature choices when collapse guides toggle changes
observeEvent(input$boxplotCollapseGuides, {
  if (isTRUE(restoring())) return()  # save-state restore guard
  if (eds$data_type() == "CrisprCountQC") {
    if (isTRUE(input$boxplotCollapseGuides)) {
      gene_choices <- sort(unique(eds$sgrna_to_gene()))
    } else {
      gene_choices <- eds$feature_ids()
    }
    # Clear everything — feature level has changed
    featureBucket$features <- character(0)
    featureBucket$version <- isolate(featureBucket$version) + 1L
    updateSelectizeInput(session, "fb_select_boxplot", choices = gene_choices, selected = NULL, server = TRUE)
    updateTextAreaInput(session, "fb_text_boxplot", value = "")
    # Clear plots
    figuresDataReactive$boxplot <- NULL
    figuresDataReactive$barplot <- NULL
  }
}, ignoreInit = TRUE)

updateSelectInput(
  session = session,
  inputId = "boxplotFactor1",
  choices = eds$factor_names(),
  selected = eds$primary_grouping()
)
updateSelectInput(
  session = session,
  inputId = "boxplotFactor2",
  choices = c("None", "Feature", eds$factor_names()),
  selected = "None"
)

# Batch + design dropdowns. For two-factor DE (param$grouping2Name set),
# auto-select batch = grouping2Name and lock design to the DE design column
# so the tab loads in the analysis's own batch-corrected view. Otherwise
# batch defaults to None and any factor can be chosen as design (free
# exploration for CountQC/CrisprCountQC, groupingName pre-selected for DE).
observe({
  .dt        <- eds$data_type()
  .factors   <- eds$factor_names()
  .grouping  <- eds$primary_grouping()
  .secondary <- eds$secondary_grouping()
  .twoFactor <- !is.na(.secondary) && nzchar(.secondary) && .secondary %in% .factors
  # Proteomics grouping: SUSHI data carries the prolfqua "G_" column; an uploaded
  # h5ad won't, so fall back to the dataset's primary grouping factor.
  .protGroup <- if ("G_" %in% .factors) "G_" else .grouping
  .deDesignCol <- if (.dt == "RNASeq") .grouping
                  else if (eds$omics_class() == "proteomics") .protGroup
                  else NA_character_

  batchChoices <- c("None", .factors)
  batchDefault <- if (.twoFactor) .secondary else "None"
  updateSelectInput(
    session = session,
    inputId = "boxplotBatch",
    choices = batchChoices,
    selected = batchDefault
  )

  designChoices <- if ((.dt == "RNASeq" || eds$omics_class() == "proteomics") && .twoFactor) {
    c("None", .deDesignCol)
  } else {
    c("None", .factors)
  }
  designDefault <- if (.dt == "RNASeq") .grouping
                   else if (eds$omics_class() == "proteomics") .protGroup
                   else "None"
  updateSelectInput(
    session = session,
    inputId = "boxplotDesign",
    choices = designChoices,
    selected = designDefault
  )
})

observeEvent(
  {
    input$boxplotFactor1
  },
  {
    output$bucket <- renderUI({
      .lvls <- unique(levels(as.factor(eds$samples()[[input$boxplotFactor1]])))
      # save-state restore: seed the include/exclude split + order from a saved
      # session (valid= drops levels that don't belong to the current factor).
      .keep <- restore_default("boxKeepBucket", .lvls, valid = .lvls)
      .excl <- restore_default("boxExcludeBucket", character(0), valid = .lvls)
      .keep <- setdiff(.keep, .excl)
      bucket_list(
        header = "Drag and drop groups in order to be plotted (Ctrl/Cmd to select multiple)",
        group_name = "bucket_list_group",
        orientation = "horizontal",
        add_rank_list(
          text = "Include these groups",
          labels = as.list(.keep),
          input_id = "boxKeepBucket",
          options = sortable::sortable_options(multiDrag = TRUE)
        ),
        add_rank_list(
          text = "Exclude these groups",
          labels = as.list(.excl),
          input_id = "boxExcludeBucket",
          options = sortable::sortable_options(multiDrag = TRUE)
        )
      )
    })
    outputOptions(output, "bucket", suspendWhenHidden = FALSE)
  }
)

# Extract the appropriate data ----
if (eds$data_type() == "RNASeq") {
  .contrast <- eds$primary_contrast()
  seqAnnoReactive <- reactiveValues("sa" = eds$de(.contrast), "saF" = NULL)
  output$boxplotDesign <- renderText({
    eds$contrast_label(.contrast)
  })
} else if (eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
  # CountQC has no DE results — leave seqAnno NULL; every downstream consumer
  # gates on dataType, so this is never queried.
  seqAnnoReactive <- reactiveValues("sa" = NULL, "saF" = NULL)
} else if (eds$omics_class() == "proteomics") {
  # Initialize seqAnno with the first available contrast immediately so that
  # bucket_add() validation works before the user explicitly changes the selector.
  .first_contrast <- if (length(eds$contrasts())) eds$contrasts()[[1]] else NULL
  .first_sa <- if (!is.null(.first_contrast)) {
    eds$de(.first_contrast) %>%
      dplyr::select("feature_id", "log2Ratio", "pValue", "fdr")
  } else NULL
  seqAnnoReactive <- reactiveValues("sa" = .first_sa, "saF" = .first_sa)
  observeEvent(input$contrastSelected,
    {
      seqAnnoReactive$sa <- eds$de(input$contrastSelected) %>%
        dplyr::select("feature_id", "log2Ratio", "pValue", "fdr")
      seqAnnoReactive$saF <- seqAnnoReactive$sa
      output$boxplotDesign <- renderText({
        input$contrastSelected
      })
    },
    ignoreInit = T,
    ignoreNULL = T
  )
}

# God observer, reset handler, and bucket render removed — now in server-featureBucket.R

# observeEvent({
#   figuresDataReactive
#   }, {
#     if (!is.null(figuresDataReactive$boxplotStatic)) {
#       output$boxplotWarn <- renderUI({NULL})
#     } else {
#       output$boxplotWarn <- renderUI({
#         div(
#           style = "text-align: center; padding: 50px; color: #666; font-size: 18px;",
#           icon("info-circle", style = "font-size: 48px; color: #3c8dbc; margin-bottom: 20px;"),
#           h3("Please Select Features", style = "color: #3c8dbc;"),
#           p("Select one or more features from the 'Features for Boxplots' dropdown above to generate your plots."),
#           p("You can also paste a list of features in the text area provided.", style = "font-style: italic;")
#         )
#       })
#     }
#   })

# Generate a reactive list of the count data to be plotted ----
# We do this so we can make aesthetic changes to the plot without having to recalculate the required boxplots and results

debouncedFeatures <- reactive({ input$fb_keep_boxplot }) %>% debounce(200)

boxplotCountsReactive <- eventReactive(
  {
    debouncedFeatures()
    input$boxplotFactor1
    input$boxKeepBucket
    input$boxplotFactor2
    input$boxplotBatch
    input$boxplotDesign
    input$boxplotCounts
    input$contrastSelected
    input$boxplotShowPHeight
    input$boxplotCountsLog
    input$boxplotZScore
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE,
  {
    req(length(debouncedFeatures()) >= 1)
    req(!is.null(input$boxplotFactor1))
    req(input$boxKeepBucket[[1]] %in% eds$samples()[[input$boxplotFactor1]])

    # Resolve the contrast's test/reference group labels (used for stat_pvalue_manual
    # brackets below). RNA-seq uses primary contrast; proteomics tracks the picker.
    designLevels <- NULL
    if (eds$omics_class() == "proteomics") {
      if (input$contrastSelected %in% eds$contrasts()) {
        cm <- eds$contrast_meta(input$contrastSelected)
        designLevels <- c(cm$test_group, cm$reference_group)
      }
    } else if (eds$data_type() == "RNASeq" && length(eds$contrasts())) {
      cm <- eds$contrast_meta(eds$primary_contrast())
      designLevels <- c(cm$test_group, cm$reference_group)
    }

    genesToPlot <- debouncedFeatures()
    if (length(genesToPlot) > 50) {
      genesToPlot <- genesToPlot[1:50]
    }
    datasetBoxplot <- eds$samples()
    datasetBoxplot$Name <- factor(datasetBoxplot$Name, levels = datasetBoxplot$Name)

    # Get the count data
    countsBoxplot <- eds$assay(input$boxplotCounts)
    # Collapse sgRNA guides to gene level if requested
    if (eds$data_type() == "CrisprCountQC" && isTRUE(input$boxplotCollapseGuides)) {
      countsBoxplot <- collapse_guides_to_genes(countsBoxplot, eds$sgrna_to_gene())
      # Map any sgRNA-level feature names to gene names
      sgRNA_map <- eds$sgrna_to_gene()
      genesToPlot <- unique(ifelse(genesToPlot %in% names(sgRNA_map), sgRNA_map[genesToPlot], genesToPlot))
      genesToPlot <- genesToPlot[genesToPlot %in% rownames(countsBoxplot)]
    }

    # Apply the following in this specific order as required: Log2, batch correct, Z-scale
    # If the user wants to log2, then do that
    is_logged <- input$boxplotCounts %in% c("normalised_log2", "vst") ||
      eds$omics_class() == "proteomics"
    if (input$boxplotCountsLog) {
      if (eds$data_type() %in% c("RNASeq", "CountQC", "CrisprCountQC")) {
        if (input$boxplotCounts %in% c("tpm", "fpkm", "normalised", "counts")) {
          countsBoxplot <- log2(countsBoxplot + 1)
          is_logged <- TRUE
        }
      }
    }
    # If the user wants to remove batch effect, then do that
    if (input$boxplotBatch != "None") {
      if (is_logged) {
        countsBoxplot <- apply_batch_correction(
          counts    = countsBoxplot,
          dataset   = datasetBoxplot,
          batchCol  = input$boxplotBatch,
          designCol = input$boxplotDesign
        )
      } else {
        showNotification("Batch correction requires log-scale data. Please select a logged count type or enable the log2 transform.",
          type = "warning", duration = 8)
      }
    }
    # If the user wants to plot gene-wise Z-scores, then do that
    if (input$boxplotZScore) {
      countsBoxplot <- t(scale(t(countsBoxplot)))
    }
    # Proteomics rowname remapping is handled inside eds$assay() via .fid_lookup.
    countsBoxplot <- as.data.frame(t(countsBoxplot))
    countsBoxplot$Name <- rownames(countsBoxplot)
    countsBoxplot <- countsBoxplot %>% dplyr::select(all_of(c("Name", genesToPlot)))

    # Add the metadata
    countsBoxplot <- left_join(countsBoxplot, datasetBoxplot[, c("Name", eds$factor_names())], by = "Name")

    # Include only the conditions selected
    countsBoxplot <- countsBoxplot[countsBoxplot[[input$boxplotFactor1]] %in% input$boxKeepBucket, ]

    # Set the factor levels per the bucket list order
    countsBoxplot[[input$boxplotFactor1]] <- factor(
      countsBoxplot[[input$boxplotFactor1]],
      input$boxKeepBucket
    )
    countsBoxplot[["Name"]] <- factor(countsBoxplot[["Name"]], levels = countsBoxplot[["Name"]])

    # melt the table
    countsBoxplotMelt <- reshape2:::melt.data.frame(data = countsBoxplot, id.vars = c("Name", eds$factor_names()), value.name = "Counts", variable.name = "Feature")
    countsBoxplotMelt$Feature <- as.factor(x = countsBoxplotMelt$Feature)

    # Force the order of the genes as factor so the order the genes are
    # selected by the user is respected by the facet wrap
    countsBoxplotMelt$Feature <- factor(countsBoxplotMelt$Feature, genesToPlot)

    # Add vector of custom shapes for compatibility with ggprism
    if (!input$boxplotFactor2 == "None" & !input$boxplotFactor2 == "Feature") {
      shapes <- c(rep(c(15, 16, 17, 18, 8), times = 10))[1:nlevels(as.factor(countsBoxplot[[input$boxplotFactor2]]))]
      names(shapes) <- levels(as.factor(countsBoxplot[[input$boxplotFactor2]]))
      countsBoxplotMelt$shapes <- shapes[match(as.character(countsBoxplotMelt[[input$boxplotFactor2]]), names(shapes))]
    }

    # add the DE results
    if (!eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
      countsBoxplotMelt$fdr <- seqAnnoReactive$sa$fdr[match(countsBoxplotMelt$Feature, seqAnnoReactive$sa$feature_id)]
      countsBoxplotMelt$fdr[is.na(countsBoxplotMelt$fdr)] <- 1
      countsBoxplotMelt$log2Ratio <- seqAnnoReactive$sa$log2Ratio[match(countsBoxplotMelt$Feature, seqAnnoReactive$sa$feature_id)]
      countsBoxplotMelt$log2Ratio[is.na(countsBoxplotMelt$log2Ratio)] <- 0
      countsBoxplotMelt <- add_significance(
        countsBoxplotMelt,
        p.col = "fdr",
        output.col = "stars",
        cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
        symbols = c("****", "***", "**", "*", "ns")
      )
      countsBoxplotMelt$group1 <- designLevels[1]
      countsBoxplotMelt$group2 <- designLevels[2]
      countsBoxplotMelt$y.position <- lapply(genesToPlot, function(g) {
        rep(
          max(countsBoxplotMelt$Counts[countsBoxplotMelt$Feature == g], na.rm = T),
          nrow(countsBoxplotMelt[countsBoxplotMelt$Feature == g, ])
        )
      }) %>% unlist()
    }

    # Get a data frame for the barplots
    countsBarplot <- dataSummary(
      data = countsBoxplotMelt[!is.na(countsBoxplotMelt$Counts), ],
      varname = "Counts",
      groupnames = c(input$boxplotFactor1, "Feature")
    )
    countsBarplot[[input$boxplotFactor1]] <- factor(countsBarplot[[input$boxplotFactor1]], levels = input$boxKeepBucket)
    countsBarplot$Feature <- factor(countsBarplot$Feature, genesToPlot)
    if (!eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
      for (f in c("fdr", "log2Ratio", "stars", "group1", "group2", "y.position")) {
        countsBarplot[[f]] <- countsBoxplotMelt[match(countsBarplot$Feature, countsBoxplotMelt$Feature), f]
      }
    }

    return(list(
      "datasetBoxplot" = datasetBoxplot,
      "countsBoxplot" = countsBoxplot,
      "countsBoxplotMelt" = countsBoxplotMelt,
      "countsBarplot" = countsBarplot
    ))
  }
)

observeEvent(
  {
    boxplotCountsReactive()
    input$boxplotYLimLFC
    input$textSizeBoxplot
    input$boxplotNCol
    input$boxplotVertLines
    input$boxplotGrey
    input$boxplotShowP
    input$boxplotFreeX
    input$boxplotThemeChoice
    input$boxplotShowBox
    input$boxplotShowPoint
    input$boxplotShowViolin
    input$boxplotPoints
    input$showDotsBarplot
    input$contrastSelected
    input$boxplotPointSize
    input$boxplotPointAlpha
    input$boxplotMeanLine
    input$boxplotShowMeanBar
    input$boxplotBoxAlpha
    input$boxplotPointDodge
    input$boxplotShowPLabelSize
    input$boxplotShowPBracketSize
    input$boxplotShowPTipSizeA
    input$boxplotShowPTipSizeB
    input$boxplotShowPDodge
    input$boxplotPointBorder
    input$boxplotMeanBarFront
    input$boxplotShowConditions
    input$boxplotConditionAngle
    input$boxplotConditionFormat
    input$boxplotFont
    input$barplotSDorSE
    input$boxplotShowPJust
    input$boxplotYScaleFree
    input$boxplotPlotBorder
    # Add all color picker inputs as dependencies so plot updates when colors change
    lapply(eds$factor_names(), function(fname) {
      lapply(eds$factor_levels(fname), function(lev) {
        input[[paste0("GroupColour", fname, lev)]]
      })
    })
  },
  ignoreNULL = FALSE,
  ignoreInit = FALSE,
  {
    req(!is.null(boxplotCountsReactive()$countsBoxplotMelt))
    # If bucket is empty (e.g., just cleared by collapse toggle), clear plots
    if (length(input$fb_keep_boxplot) < 1) {
      figuresDataReactive$boxplot <- NULL
      figuresDataReactive$barplot <- NULL
      return()
    }
    datasetBoxplot <- boxplotCountsReactive()$datasetBoxplot
    countsBoxplot <- boxplotCountsReactive()$countsBoxplot
    countsBoxplotMelt <- boxplotCountsReactive()$countsBoxplotMelt
    countsBarplot <- boxplotCountsReactive()$countsBarplot
    genesToPlot <- levels(as.factor(countsBoxplotMelt$variable))

    # Get the colours for the condition levels
    coloursBoxplot <- setNames(lapply(input$boxKeepBucket, function(k) {
      paste0(col2hex(input[[paste0("GroupColour", input$boxplotFactor1, k)]]), "FF")
    }), input$boxKeepBucket)

    if (eds$omics_class() == "proteomics" && input$contrastSelected %in% eds$contrasts()) {
      cm <- eds$contrast_meta(input$contrastSelected)
      designLevels <- c(cm$test_group, cm$reference_group)
    }

    facetFree <- "free"
    if (input$boxplotYScaleFree) {
      countsBoxplotMelt$y.position <- max(countsBoxplotMelt$y.position)
      facetFree <- "free_x"
    }

    # Build a facet labeller that annotates feature_ids with gene symbols
    # (proteomics only, when gene_name info exists in the seqAnno)
    facet_labeller <- "label_value"
    if (eds$omics_class() == "proteomics" && !is.null(input$contrastSelected) &&
        input$contrastSelected %in% eds$contrasts()) {
      sa_full <- eds$de(input$contrastSelected)
      if ("gene_name" %in% colnames(sa_full)) {
        gene_map <- setNames(sa_full$gene_name, sa_full$feature_id)
        facet_labeller <- labeller(Feature = function(fids) {
          gnames <- gene_map[fids]
          ifelse(!is.na(gnames) & gnames != fids,
                 paste0(gnames, "\n", fids),
                 fids)
        })
      }
    }

    # Boxplot ----
    boxplot <- ggplot(data = countsBoxplotMelt, aes(x = .data[[input$boxplotFactor1]], y = Counts))
    if (input$boxplotShowBox & !input$boxplotShowViolin) {
      boxplot <- boxplot + geom_boxplot(aes(fill = .data[[input$boxplotFactor1]]), alpha = input$boxplotBoxAlpha, outlier.shape = NA, show.legend = FALSE)
    }
    if (!input$boxplotShowBox & input$boxplotShowViolin) {
      boxplot <- boxplot + geom_violin(aes(fill = .data[[input$boxplotFactor1]]), alpha = input$boxplotBoxAlpha, draw_quantiles = 0.5, show.legend = FALSE)
    }
    if (input$boxplotShowMeanBar & !input$boxplotMeanBarFront) {
      boxplot <- boxplot + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$boxplotMeanLine, show.legend = FALSE)
    }
    if (input$boxplotShowPoint) {
      if (input$boxplotFactor2 == "None" | input$boxplotFactor2 == "Feature") {
        boxplot <- boxplot + geom_beeswarm(aes(fill = .data[[input$boxplotFactor1]]), size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, shape = 21, cex = input$boxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9, stroke = input$boxplotPointBorder)
      } else {
        boxplot <- boxplot + geom_beeswarm(aes(fill = .data[[input$boxplotFactor1]], shape = .data[[input$boxplotFactor2]]), size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, cex = input$boxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9, stroke = input$boxplotPointBorder)
        boxplot <- boxplot + scale_shape_manual(values = c(rep(c(21, 22, 23, 24, 25), times = 10))[1:nlevels(as.factor(countsBoxplotMelt[[input$boxplotFactor2]]))])
        boxplot <- boxplot + guides(
          fill = guide_legend(override.aes = list(shape = 21)),
          shape = guide_legend(override.aes = list(shape = c(rep(c(16, 15, 18, 17, 25), times = 10))[1:nlevels(as.factor(countsBoxplotMelt[[input$boxplotFactor2]]))]))
        )
      }
    }
    if (input$boxplotShowMeanBar & input$boxplotMeanBarFront) {
      boxplot <- boxplot + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$boxplotMeanLine, show.legend = FALSE)
    }
    if (!input$boxplotGrey) {
      boxplot <- boxplot + scale_fill_manual(breaks = names(coloursBoxplot), values = as.character(coloursBoxplot))
    } else {
      boxplot <- boxplot + scale_colour_grey()
      boxplot <- boxplot + scale_fill_grey()
    }
    if (input$boxplotConditionFormat) {
      xlabs <- levels(countsBoxplotMelt[[input$boxplotFactor1]])
      xlabs <- gsub("\\_|\\.|\\ |\\/", "\n", xlabs)
      boxplot <- boxplot + scale_x_discrete(labels = xlabs)
    }
    ylabToPlot <- paste(input$boxplotCounts, "Counts")
    if (input$boxplotCountsLog) {
      ylabToPlot <- paste(ylabToPlot, "(Log2)")
    }
    if (input$boxplotZScore) {
      ylabToPlot <- paste(ylabToPlot, "[Z-score]")
    }
    boxplot <- boxplot + labs(y = ylabToPlot)
    boxplot <- boxplot + facet_wrap(~Feature, scales = facetFree, ncol = input$boxplotNCol, labeller = facet_labeller)
    if (input$boxplotVertLines) {
      boxplot <- boxplot + geom_vline(xintercept = seq_along(input$boxKeepBucket)[-length(seq_along(input$boxKeepBucket))] + 0.5, linetype = "dashed", alpha = 0.7)
    }
    if (input$boxplotShowP) {
      boxplot <- boxplot + stat_pvalue_manual(
        data = countsBoxplotMelt,
        y.position = countsBoxplotMelt$y.position * 1.03,
        label = "stars",
        size = input$boxplotShowPLabelSize * 2,
        label.size = input$boxplotShowPLabelSize * 2,
        bracket.size = input$boxplotShowPBracketSize / 6,
        tip.length = c(input$boxplotShowPTipSizeA / 50, input$boxplotShowPTipSizeB / 50),
        bracket.nudge.y = input$boxplotShowPJust / 3
      ) + coord_cartesian(clip = "off")
    }
    boxplot <- boxplot + theme_prism(base_size = input$textSizeBoxplot, axis_text_angle = input$boxplotConditionAngle, base_family = input$boxplotFont, border = input$boxplotPlotBorder)
    if (!input$boxplotShowConditions) {
      boxplot <- boxplot + theme(axis.text.x = element_blank())
    }
    figuresDataReactive$boxplot <- boxplot


    # Barplot ----
    countsBarplot <- countsBarplot %>% complete(!!sym(input$boxplotFactor1), Feature, fill = list(Counts = 0))
    barplot <- ggplot(data = countsBarplot, aes(x = .data[[input$boxplotFactor1]], y = Counts, fill = .data[[input$boxplotFactor1]]))
    barplot <- barplot + facet_wrap(~Feature, scales = "free", ncol = input$boxplotNCol, labeller = facet_labeller)
    barplot <- barplot + geom_bar(stat = "identity", color = "black", position = position_dodge(), alpha = input$boxplotBoxAlpha, show.legend = FALSE)
    barplot <- barplot + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
    if (input$barplotSDorSE %in% c("SE", "SD")) {
      barplot <- barplot + geom_errorbar(aes(ymin = Counts - .data[[tolower(input$barplotSDorSE)]], ymax = Counts + .data[[tolower(input$barplotSDorSE)]]), width = .2, position = position_dodge(.9))
    }
    if (input$showDotsBarplot) {
      barplot <- barplot + geom_beeswarm(data = countsBoxplotMelt, aes_string(x = input$boxplotFactor1, y = "Counts"), pch = 21, size = input$boxplotPointSize, alpha = input$boxplotPointAlpha, cex = 3, corral = "gutter", corral.width = 0.9)
    }
    if (input$boxplotConditionFormat) {
      barplot <- barplot + scale_x_discrete(labels = xlabs)
    }
    barplot <- barplot + labs(y = ylabToPlot)
    if (!input$boxplotGrey) {
      barplot <- barplot + scale_fill_manual(breaks = names(coloursBoxplot), values = as.character(coloursBoxplot))
    } else {
      barplot <- barplot + scale_colour_grey()
      barplot <- barplot + scale_fill_grey()
    }
    barplot <- barplot + theme_prism(base_size = input$textSizeBoxplot, axis_text_angle = input$boxplotConditionAngle, base_family = input$boxplotFont)
    if (!input$boxplotShowConditions) {
      barplot <- barplot + theme(axis.text.x = element_blank())
    }
    if (input$boxplotShowP) {
      barplot <- barplot + stat_pvalue_manual(
        data = countsBarplot,
        y.position = countsBarplot$y.position * 1.03,
        label = "stars",
        size = input$boxplotShowPLabelSize * 2,
        label.size = input$boxplotShowPLabelSize * 2,
        bracket.size = input$boxplotShowPBracketSize / 6,
        tip.length = c(input$boxplotShowPTipSizeA / 50, input$boxplotShowPTipSizeB / 50),
        bracket.nudge.y = input$boxplotShowPJust / 3
      ) + coord_cartesian(clip = "off")
    }
    figuresDataReactive$barplot <- barplot
  }
)

output$boxplotTable <- renderDataTable({
  req(boxplotCountsReactive()$countsBoxplotMelt)
  countsBoxplotMelt <- boxplotCountsReactive()$countsBoxplotMelt
  boxplot_brushed_df <- brushedPoints(
    countsBoxplotMelt,
    input$boxplotBrush,
    xvar = input$boxplotFactor1,
    yvar = "Counts"
  )
  if (!eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
    colsToShowHere <- c("Name", input$boxplotFactor1, "Feature", "Counts", "log2Ratio", "fdr")
    colsToFormat <- c("Counts", "fdr", "log2Ratio")
  } else {
    colsToShowHere <- c("Name", input$boxplotFactor1, "Feature", "Counts")
    colsToFormat <- c("Counts")
  }
  DT::datatable(
    data = boxplot_brushed_df[, colsToShowHere, drop = FALSE],
    filter = "top",
    rownames = F
  ) %>%
    DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
    formatSignif(colsToFormat, digits = 3)
})

lapply(c("boxplot", "barplot"), function(b) {
  output[[paste0(b, "Static")]] <- renderPlot(
    {
      req(!is.null(figuresDataReactive[[b]]))
      figuresDataReactive[[b]]
    },
    width = function() {
      as.numeric(input$figWidthBoxplot)
    },
    height = function() {
      as.numeric(input$figHeightBoxplot)
    }
  )
  
  # Map internal name 'b' to new button ID suffix
  button_id <- paste0("dlPlot_", b)
  
  output[[button_id]] <- downloadHandler(
    filename = function() {
      plot_type_suffix <- if (b == "boxplot") "_Boxplot" else "_Barplot"
      generate_filename(
        base_name = paste0(input$filename_boxplot, plot_type_suffix),
        design_name = designReactive(),
        extension = tolower(input$downloadFormatBoxplot)
      )
    },
    content = function(file) {
      req(!is.null(figuresDataReactive[[b]]))
      export_plot_standard(
        plot_object = figuresDataReactive[[b]],
        filename = file,
        format = tolower(input$downloadFormatBoxplot),
        width = as.numeric(input$figWidthBoxplot) / 75,
        height = as.numeric(input$figHeightBoxplot) / 75,
        dpi = as.numeric(input$dpiBoxplot)
      )
    }
  )
})
output$dlTable_boxplot <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_boxplot, "_Counts"),
      design_name = designReactive(),
      extension = "xlsx"
    )
  },
  content = function(file) {
    req(!is.null(boxplotCountsReactive()$countsBoxplotMelt))
    counts_data <- boxplotCountsReactive()$countsBoxplotMelt

    # Determine which columns to include
    # Basic columns: Name, factor columns, Feature, Counts
    basic_cols <- c("Name", eds$factor_names(), "Feature", "Counts")

    if ((eds$data_type() == "RNASeq" || eds$omics_class() == "proteomics") &&
        !is.null(input$includeDEStatsBoxplot) &&
        input$includeDEStatsBoxplot) {
      # Add DE statistics: log2Ratio, pValue, fdr
      # First add pValue if it doesn't exist (it's not in countsBoxplotMelt by default)
      if (eds$data_type() == "RNASeq") {
        seqAnno_dl <- eds$de(eds$primary_contrast())
        if ("pValue" %in% colnames(seqAnno_dl)) {
          counts_data$pValue <- seqAnno_dl$pValue[match(counts_data$Feature, seqAnno_dl$feature_id)]
        }
      } else if (!is.null(input$contrastSelected) &&
                 input$contrastSelected %in% eds$contrasts()) {
        seqAnno_dl <- eds$de(input$contrastSelected)
        if ("pValue" %in% colnames(seqAnno_dl)) {
          counts_data$pValue <- seqAnno_dl$pValue[match(counts_data$Feature, seqAnno_dl$feature_id)]
        }
      }

      # Select columns: basic + DE stats (log2Ratio and fdr already exist, pValue just added)
      de_cols <- c("log2Ratio", "pValue", "fdr")
      available_de_cols <- intersect(de_cols, colnames(counts_data))
      cols_to_export <- c(basic_cols, available_de_cols)
    } else {
      # Export only basic columns
      cols_to_export <- basic_cols
    }

    # Filter to selected columns
    counts_data <- counts_data[, intersect(cols_to_export, colnames(counts_data)), drop = FALSE]

    export_excel_standard(
      df = counts_data,
      filename = file,
      sheet_name = "Boxplot Data",
      table_name = "BoxplotData"
    )
  }
)

# Show/Hide DE P-values tab based on data type
observe({
  if (eds$data_type() == "RNASeq" || eds$omics_class() == "proteomics") {
    showTab(inputId = "boxplotSettingsTabs", target = "dePValuesTab")
  } else {
    hideTab(inputId = "boxplotSettingsTabs", target = "dePValuesTab")
  }
})
