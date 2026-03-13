
dataSummary <- function(data, varname, groupnames) {
  summary_func <- function(x, col) {
    c(
      mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE),
      se = sqrt(sd(x[[col]]) / length(x[[col]]))
    )
  }
  data_sum <- plyr::ddply(
    data, groupnames,
    .fun = summary_func, varname
  )
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# Proxy Data Reactive for Single-Dataset Emulation
multiBoxplotDataReactive <- reactive({
    req(multiDeReactive$results)
    req(length(multiDeReactive$results) >= 1)
    # Return the first result as the primary dataset source
    # This tricks the boxplot logic into thinking it has a standard single input
    multiDeReactive$results[[1]]
})


# Update a bunch of input options ----
output$multiBoxplotStatic <- renderPlot({
  NULL
})

observe({
    req(multiBoxplotDataReactive())
    d <- multiBoxplotDataReactive()
    
    updateSelectInput(
      session = session,
      inputId = "multiBoxplotCounts",
      choices = names(d$countList),
      selected = switch(d$dataType,
        "RNASeq" = "Normalised + Log2",
        "CountQC" = "Normalised + Log2",
        "proteomics" = "transformedData",
        "Normalised + Log2" # Fallback
      )
    )
    
    updateSelectizeInput(
      session = session,
      inputId = "multiBoxplotGenes",
      choices = d$genes,
      selected = "",
      server = TRUE
    )
    
    updateSelectInput(
      session = session,
      inputId = "multiBoxplotFactor1",
      choices = d$factors,
      selected = d$groupingName
    )
    updateSelectInput(
      session = session,
      inputId = "multiBoxplotFactor2",
      choices = c("None", "Feature", d$factors),
      selected = "None"
    )
    
    # Automatically set batch correction if grouping2Name exists
    batchDefault <- "None"
    if (!is.null(d$param$grouping2Name) && !is.null(d$factors)) {
      if (d$param$grouping2Name %in% d$factors) {
        batchDefault <- d$param$grouping2Name
      }
    }
    updateSelectInput(
      session = session,
      inputId = "multiBoxplotBatch",
      choices = c("None", d$factors),
      selected = batchDefault
    )
})

observeEvent(
  {
    input$multiBoxplotFactor1
  },
  {
    output$multiBoxplotBucket <- renderUI({
      req(multiBoxplotDataReactive())
      d <- multiBoxplotDataReactive()
      req(input$multiBoxplotFactor1 %in% colnames(d$dataset))
      
      bucket_list(
        header = "Drag and drop groups in order to be plotted (Ctrl/Cmd to select multiple)",
        group_name = "multi_gemini_bucket_list_group",
        orientation = "horizontal",
        add_rank_list(
          text = "Include these groups",
          labels = as.list(unique(levels(as.factor(d$dataset[[input$multiBoxplotFactor1]])))),
          input_id = "multiBoxKeepBucket",
          options = sortable::sortable_options(multiDrag = TRUE)
        ),
        add_rank_list(
          text = "Exclude these groups",
          labels = NULL,
          input_id = "multiBoxExcludeBucket",
          options = sortable::sortable_options(multiDrag = TRUE)
        )
      )
    })
    outputOptions(output, "multiBoxplotBucket", suspendWhenHidden = FALSE)
  }
)

# Extract the appropriate data ----
# Using proxy, mostly for "design" text or simple access
observe({
    req(multiBoxplotDataReactive())
    d <- multiBoxplotDataReactive()
    output$multiBoxplotDesign <- renderText({
       # Just show info about number of contrasts loaded
       paste("Plots based on", length(multiDeReactive$results), "loaded contrasts.")
    })
})
# Note: Legacy seqAnnoReactive block removed as we handle multi-seqAnno differently.


# The main part of the gene bucket lives here ----
multiGenesReactive <- reactiveValues(genes = NULL)

observeEvent(
  {
    input$multiBoxplotGenes
    input$multiBoxplotGenesText
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE,
  {
    d <- multiBoxplotDataReactive()
    req(d)
    
    # Get all the genes that have been selected from the drop-downs
    genesList1 <- if (!is.null(input$multiBoxplotGenes)) input$multiBoxplotGenes else character(0)

    # Get all the genes that have been typed into the textbox
    genesList2 <- if (!is.null(input$multiBoxplotGenesText)) {
        strsplit(input$multiBoxplotGenesText, "\\ +|\\/+|\\,+|\n+|\r+|\n\r+|\t+")[[1]]
    } else {
        character(0)
    }

    genesUnlist <- unique(c(genesList1, genesList2))
    genesUnlist <- genesUnlist[genesUnlist != ""]
    # Verify against first contrast genes
    genesUnlist <- genesUnlist[genesUnlist %in% d$genes]
    
    multiGenesReactive$genes <- genesUnlist
  }
)

observeEvent(input$resetGeneBucketMultiBoxplot,
  {
    multiGenesReactive$genes <- NULL
    updateSelectizeInput(session, inputId = "multiBoxplotGenes", selected = NULL)
    updateTextAreaInput(session, inputId = "multiBoxplotGenesText", value = NULL)
    output$multiBoxplotStatic <- NULL
  },
  ignoreNULL = T,
  ignoreInit = T
)


observeEvent(
  {
    multiGenesReactive$genes
  },
  {
    output$multiBoxplotGeneBucket1 <- renderUI({
      bucket_list(
        header = "Drag and drop features in order to be plotted (Ctrl/Cmd to select multiple)",
        group_name = "multi_gemini_feature_bucket_group",
        orientation = "horizontal",
        add_rank_list(
          text = "Include these features in this order",
          labels = multiGenesReactive$genes,
          input_id = "multiBoxFeatureKeepBucket",
          options = sortable::sortable_options(multiDrag = TRUE)
        ),
        add_rank_list(
          text = "Exclude these features",
          labels = NULL,
          input_id = "multiBoxFeatureExcludeBucket",
          options = sortable::sortable_options(multiDrag = TRUE)
        )
      )
    })
    outputOptions(output, "multiBoxplotGeneBucket1", suspendWhenHidden = FALSE)
  }
)


# Define catalyst colour palette (from exploreDE) ----


# Create colour pickers for each group ----




# Generate a reactive list of the count data to be plotted ----
debouncedMultiBoxplotFactor1 <- reactive({
  input$multiBoxFeatureKeepBucket
}) %>% debounce(200)

multiBoxplotCountsReactive <- eventReactive(
  {
    debouncedMultiBoxplotFactor1()
    input$multiBoxplotFactor1
    input$multiBoxKeepBucket
    input$multiBoxplotFactor2
    input$multiBoxplotBatch
    input$multiBoxplotCounts
    input$multiBoxplotCountsLog
    input$multiBoxplotZScore
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE,
  {
    req(length(debouncedMultiBoxplotFactor1()) >= 1)
    req(!is.null(input$multiBoxplotFactor1))
    
    # DATA HYGIENE CHECKING
    validate(
        need(length(debouncedMultiBoxplotFactor1()) >= 1, "No genes selected. Please select genes to plot."),
        need(input$multiBoxplotFactor1 %in% colnames(multiBoxplotDataReactive()$dataset), "Selected grouping factor is not present in the dataset.")
    )
    
    d <- multiBoxplotDataReactive() # Proxy
    req(input$multiBoxKeepBucket[[1]] %in% d$dataset[[input$multiBoxplotFactor1]])

    genesToPlot <- debouncedMultiBoxplotFactor1()
    if (length(genesToPlot) > 50) {
      genesToPlot <- genesToPlot[1:50]
    }

    datasetBoxplot <- d$dataset
    # STRICT Order
    datasetBoxplot$Name <- factor(datasetBoxplot$Name, levels = datasetBoxplot$Name)

    # Get the count data
    countsBoxplot <- d$countList[[input$multiBoxplotCounts]]

    # Transformations
    if (input$multiBoxplotCountsLog) {
       if (d$dataType %in% c("RNASeq", "CountQC")) {
         if (input$multiBoxplotCounts %in% c("TPM", "FPKM", "Normalised", "Raw")) {
            countsBoxplot <- log2(countsBoxplot + 1)
         }
       }
    }
    if (input$multiBoxplotBatch != "None") {
      countsBoxplot <- limma::removeBatchEffect(countsBoxplot, datasetBoxplot[[input$multiBoxplotBatch]])
    }
    if (input$multiBoxplotZScore) {
      countsBoxplot <- t(scale(t(countsBoxplot)))
    }
    if (d$dataType == "proteomics") {
      rownames(countsBoxplot) <- gsub("\\~.*", "", rownames(countsBoxplot))
    }

    countsBoxplot <- as.data.frame(t(countsBoxplot))
    countsBoxplot$Name <- rownames(countsBoxplot)
    countsBoxplot <- countsBoxplot %>% dplyr::select(all_of(c("Name", genesToPlot))) # strictly select only needed

    countsBoxplot <- left_join(countsBoxplot, datasetBoxplot[, c("Name", d$factors)], by = "Name")
    countsBoxplot <- countsBoxplot[countsBoxplot[[input$multiBoxplotFactor1]] %in% input$multiBoxKeepBucket, ]

    # Set factor levels strictly per bucket list
    countsBoxplot[[input$multiBoxplotFactor1]] <- factor(
      countsBoxplot[[input$multiBoxplotFactor1]],
      input$multiBoxKeepBucket
    )
    countsBoxplot[["Name"]] <- factor(countsBoxplot[["Name"]], levels = countsBoxplot[["Name"]])

    countsBoxplotMelt <- reshape2::melt(data = countsBoxplot, id.vars = c("Name", d$factors), value.name = "Counts", variable.name = "Feature")
    countsBoxplotMelt$Feature <- as.factor(x = countsBoxplotMelt$Feature)

    # Strictly force feature order
    countsBoxplotMelt$Feature <- factor(countsBoxplotMelt$Feature, genesToPlot)

    if (!input$multiBoxplotFactor2 == "None" & !input$multiBoxplotFactor2 == "Feature") {
      shapes <- c(rep(c(15, 16, 17, 18, 8), times = 10))[1:nlevels(as.factor(countsBoxplot[[input$multiBoxplotFactor2]]))]
      names(shapes) <- levels(as.factor(countsBoxplot[[input$multiBoxplotFactor2]]))
      countsBoxplotMelt$shapes <- shapes[match(as.character(countsBoxplotMelt[[input$multiBoxplotFactor2]]), names(shapes))]
    }

    # NEW MULTI-CONTRAST P-VALUE LOGIC
    # Replaces single-contrast logic from boxplot.R
    pvalueDF <- list()
    for (nm in names(multiDeReactive$results)) {
        res <- multiDeReactive$results[[nm]]
        if (!is.null(res$seqAnno)) {
            # Find comparison groups
            group1 <- NULL
            group2 <- NULL

            # Helper to parse param$comparison or design string
            designStr <- if (!is.null(res$param$comparison)) res$param$comparison else if (!is.null(res$design)) res$design else ""

            if (grepl("--over--", designStr)) {
                parts <- strsplit(designStr, "--over--")[[1]]
                group1 <- parts[1]
                group2 <- parts[2]
            } else if (grepl("_vs_", designStr)) {
                 parts <- strsplit(designStr, "_vs_")[[1]]
                 group1 <- parts[1]
                 group2 <- parts[2]
            } else if (grepl(" vs ", designStr)) {
                 parts <- strsplit(designStr, " vs ")[[1]]
                 group1 <- parts[1]
                 group2 <- parts[2]
            }

            # If we successfully parsed groups, and genes match
            if (!is.null(group1) && !is.null(group2)) {

                 # Filtering Step: Only show P-values if the contrast groups exist in the currently plotted factor
                 currentFactorLevels <- unique(datasetBoxplot[[input$multiBoxplotFactor1]])
                 if (!all(c(group1, group2) %in% currentFactorLevels)) {
                     next
                 }

                 sa <- res$seqAnno
                 # Ensure unique Gene ID matching if possible, or use gene_name
                 # Assuming feature level gene_name
                 matched <- sa[sa$feature_id %in% genesToPlot, ]

                 if (nrow(matched) > 0) {
                     df <- data.frame(
                         Feature = matched$feature_id,
                         fdr = matched$fdr,
                         log2Ratio = matched$log2Ratio,
                         Contrast = nm,
                         group1 = group1,
                         group2 = group2,
                         stringsAsFactors = FALSE
                     )
                     pvalueDF[[nm]] <- df
                 }
            }
        }
    }
    pvalueDF <- do.call(rbind, pvalueDF)

    if (!is.null(pvalueDF)) {
         pvalueDF$Feature <- factor(pvalueDF$Feature, levels = genesToPlot)
         pvalueDF <- rstatix::add_significance(
            pvalueDF,
            p.col = "fdr",
            output.col = "stars",
            cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
            symbols = c("****", "***", "**", "*", "ns")
         )

         # Calculate y.positions based on max counts
         # We need max count per feature
         maxC <- tapply(countsBoxplotMelt$Counts, countsBoxplotMelt$Feature, max, na.rm=TRUE)

         pvalueDF$y.position <- 0
         for (i in 1:nrow(pvalueDF)) {
             ft <- as.character(pvalueDF$Feature[i])
             if (ft %in% names(maxC)) {
                 pvalueDF$y.position[i] <- maxC[ft]
             }
         }
    }

    # Barplot Summary
    countsBarplot <- dataSummary(
      data = countsBoxplotMelt[!is.na(countsBoxplotMelt$Counts), ],
      varname = "Counts",
      groupnames = c(input$multiBoxplotFactor1, "Feature")
    )
    countsBarplot[[input$multiBoxplotFactor1]] <- factor(countsBarplot[[input$multiBoxplotFactor1]], levels = input$multiBoxKeepBucket)
    countsBarplot$Feature <- factor(countsBarplot$Feature, genesToPlot)

    return(list(
      "datasetBoxplot" = datasetBoxplot,
      "countsBoxplot" = countsBoxplot,
      "countsBoxplotMelt" = countsBoxplotMelt,
      "countsBarplot" = countsBarplot,
      "pvalueDF" = pvalueDF
    ))
  }
)

observeEvent(
  {
    multiBoxplotCountsReactive()
    input$multiBoxplotYScaleFree
    input$textSizeMultiBoxplot
    input$multiBoxplotNCol
    input$multiBoxplotVertLines
    input$multiBoxplotGrey
    input$multiBoxplotShowP
    # ...
    input$multiBoxplotShowBox
    input$multiBoxplotShowPoint
    input$multiBoxplotShowViolin
    input$showDotsMultiBarplot
    input$multiBoxplotPointSize
    input$multiBoxplotPointAlpha
    input$multiBoxplotMeanLine
    input$multiBoxplotShowMeanBar
    input$multiBoxplotBoxAlpha
    input$multiBoxplotPointDodge
    input$multiBoxplotShowPLabelSize
    input$multiBoxplotShowPBracketSize
    input$multiBoxplotShowPTipSizeA
    input$multiBoxplotShowPTipSizeB
    input$multiBoxplotShowPJust
    input$multiBoxplotShowPDodge
    input$multiBoxplotShowPBracketShorten
    input$multiBoxplotPointBorder
    input$multiBoxplotMeanBarFront
    input$multiBoxplotShowConditions
    input$multiBoxplotConditionAngle
    input$multiBoxplotConditionFormat
    input$multiBarplotSDorSE
    input$multiBoxplotPStack
    input$multiBoxplotPStackOffset
    input$multiBoxplotFont
    input$multiBoxplotPlotBorder
    input$multiBoxplotColourP
    # Colour dependencies
    # Watch for changes in Contrast Colours (for P-values)
    lapply(names(multiDeReactive$results), function(n) {
      input[[paste0("multiGroupColour", gsub(" ", "_", n))]]
    })
    
    # Watch for changes in Group Colours (for Boxplot/Dots)
    # We need to access the current dataset levels to know which inputs to watch
    if (!is.null(multiBoxplotCountsReactive())) {
        d_melt <- multiBoxplotCountsReactive()$datasetBoxplot
        levels_to_watch <- unique(d_melt[[input$multiBoxplotFactor1]])
        lapply(levels_to_watch, function(lvl) {
            input[[paste0("GroupColour", input$multiBoxplotFactor1, lvl)]]
        })
    }
  },
  ignoreNULL = FALSE,
  ignoreInit = FALSE,
  {
    req(!is.null(multiBoxplotCountsReactive()$countsBoxplotMelt))
    d <- multiBoxplotCountsReactive()
    datasetBoxplot <- d$datasetBoxplot
    countsBoxplot <- d$countsBoxplot
    countsBoxplotMelt <- d$countsBoxplotMelt
    countsBarplot <- d$countsBarplot
    pvalueDF <- d$pvalueDF

    genesToPlot <- levels(as.factor(countsBoxplotMelt$Feature))

    # Fetch Group Colours
    groupLevels <- levels(as.factor(countsBoxplotMelt[[input$multiBoxplotFactor1]]))
    cols <- setNames(lapply(groupLevels, function(lvl) {

        # Check Main App Group Colours (from Colours Tab)
        # ID format from server-summary.R: paste0("GroupColour", FactorName, LevelName)
        # We need to construct this ID using the selected factor and the current level.
        # Note: server-summary.R iterates factorLevels. We assume input$multiBoxplotFactor1 is the FactorName.
        inputId <- paste0("GroupColour", input$multiBoxplotFactor1, lvl)
        col <- input[[inputId]]

        if (is.null(col)) {
             # Fallback if input not ready
             "grey50"
        } else {
             paste0(gplots::col2hex(col), "FF")
        }
    }), groupLevels)

    cols <- unlist(cols)


    facetFree <- "free"
    if (input$multiBoxplotYScaleFree) {
       facetFree <- if (input$multiBoxplotYScaleFree) "free" else "fixed"
    }

    # Boxplot ----
    boxplot <- ggplot(data = countsBoxplotMelt, aes(x = .data[[input$multiBoxplotFactor1]], y = Counts))
    if (input$multiBoxplotShowBox & !input$multiBoxplotShowViolin) {
      boxplot <- boxplot + geom_boxplot(aes(fill = .data[[input$multiBoxplotFactor1]]), alpha = input$multiBoxplotBoxAlpha, outlier.shape = NA, show.legend = FALSE)
    }
    if (!input$multiBoxplotShowBox & input$multiBoxplotShowViolin) {
      boxplot <- boxplot + geom_violin(aes(fill = .data[[input$multiBoxplotFactor1]]), alpha = input$multiBoxplotBoxAlpha, draw_quantiles = 0.5, show.legend = FALSE)
    }
    if (input$multiBoxplotShowMeanBar & !input$multiBoxplotMeanBarFront) {
      boxplot <- boxplot + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$multiBoxplotMeanLine, show.legend = FALSE)
    }
    if (input$multiBoxplotShowPoint) {
       if (input$multiBoxplotFactor2 == "None" | input$multiBoxplotFactor2 == "Feature") {
        boxplot <- boxplot + geom_beeswarm(aes(fill = .data[[input$multiBoxplotFactor1]]), size = input$multiBoxplotPointSize, alpha = input$multiBoxplotPointAlpha, shape = 21, cex = input$multiBoxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9, stroke = input$multiBoxplotPointBorder)
       } else {
        boxplot <- boxplot + geom_beeswarm(aes(fill = .data[[input$multiBoxplotFactor1]], shape = .data[[input$multiBoxplotFactor2]]), size = input$multiBoxplotPointSize, alpha = input$multiBoxplotPointAlpha, cex = input$multiBoxplotPointDodge, method = "swarm", corral = "none", corral.width = 0.9, stroke = input$multiBoxplotPointBorder)
        boxplot <- boxplot + scale_shape_manual(values = c(rep(c(21, 22, 23, 24, 25), times = 10))[1:nlevels(as.factor(countsBoxplotMelt[[input$multiBoxplotFactor2]]))])
       }
    }
    if (input$multiBoxplotShowMeanBar & input$multiBoxplotMeanBarFront) {
      boxplot <- boxplot + stat_summary(geom = "crossbar", fun = mean, position = position_dodge(0.9), colour = "black", size = 0.4, width = input$multiBoxplotMeanLine, show.legend = FALSE)
    }

    if (input$multiBoxplotGrey) {
      boxplot <- boxplot + scale_colour_grey() + scale_fill_grey()
    } else {
      boxplot <- boxplot + scale_fill_manual(values = cols)
    }

    if (input$multiBoxplotConditionFormat) {
      xlabs <- levels(countsBoxplotMelt[[input$multiBoxplotFactor1]])
      xlabs <- gsub("\\_|\\.|\\ |\\/", "\n", xlabs)
      boxplot <- boxplot + scale_x_discrete(labels = xlabs)
    }

    ylabToPlot <- paste(input$multiBoxplotCounts, "Counts")
    if (input$multiBoxplotCountsLog) ylabToPlot <- paste(ylabToPlot, "(Log2)")
    if (input$multiBoxplotZScore) ylabToPlot <- paste(ylabToPlot, "[Z-score]")
    boxplot <- boxplot + labs(y = ylabToPlot)

    boxplot <- boxplot + facet_wrap(~Feature, scales = facetFree, ncol = input$multiBoxplotNCol)

    if (input$multiBoxplotVertLines) {
      boxplot <- boxplot + geom_vline(xintercept = seq_along(input$multiBoxKeepBucket)[-length(seq_along(input$multiBoxKeepBucket))] + 0.5, linetype = "dashed", alpha = 0.7)
    }

    # MULTI-CONTRAST P-VALUE PLOTTING
    if (input$multiBoxplotShowP && !is.null(pvalueDF) && nrow(pvalueDF) > 0) {

        # Get unique contrasts
        contrasts <- unique(pvalueDF$Contrast)

        for (i in seq_along(contrasts)) {
            con <- contrasts[i]
            dat <- pvalueDF[pvalueDF$Contrast == con, ]

            # Stacking offset
            # If stack is ON, add offset * i.
            # Base Y is already set to max(count).
            offset <- 0
            if (input$multiBoxplotPStack) {
                offset <- (i - 1) * (max(countsBoxplotMelt$Counts, na.rm=TRUE) * (input$multiBoxplotPStackOffset / 10))
            }

            # Check if y.position exists and add offset
            if (!is.null(dat$y.position)) {
                # P Height (Dodge) controls the buffer above the boxplot default y
                buffer <- if (input$multiBoxplotShowPDodge > 0) (input$multiBoxplotShowPDodge / 20) else 0
                dat$y.position <- dat$y.position + (dat$y.position * buffer) + offset
            }

            # Fetch colour from Data Selection inputs
            colInput <- input[[paste0("multiGroupColour", gsub(" ", "_", con))]]
            if (is.null(colInput)) colInput <- "black"

            # User option to force black P-values
            if (!input$multiBoxplotColourP) {
                colInput <- "black"
            }

            # Manual Bracket Shortening
            # Map group1/group2 to x-coordinates
            # We must use the exact levels used in the plot
            x_levels <- levels(as.factor(countsBoxplotMelt[[input$multiBoxplotFactor1]]))
            dat$x1 <- match(dat$group1, x_levels)
            dat$x2 <- match(dat$group2, x_levels)

            # Apply shortening symmetrically
            # Logic: Center = (x1 + x2) / 2. Original Width = abs(x1 - x2).
            # New Width = Original Width * (1 - shorten_fraction)
            # xmin = Center - (New Width / 2)
            # xmax = Center + (New Width / 2)
            if (input$multiBoxplotShowPBracketShorten > 0) {
                center <- (dat$x1 + dat$x2) / 2
                width <- abs(dat$x1 - dat$x2)
                new_width <- width * (1 - input$multiBoxplotShowPBracketShorten)
                dat$xmin <- center - (new_width / 2)
                dat$xmax <- center + (new_width / 2)

                # Use xmin/xmax mapping
                boxplot <- boxplot + stat_pvalue_manual(
                    data = dat,
                    y.position = "y.position",
                    xmin = "xmin",
                    xmax = "xmax",
                    label = "stars",
                    size = input$multiBoxplotShowPLabelSize * 2,
                    label.size = input$multiBoxplotShowPLabelSize * 2,
                    bracket.size = input$multiBoxplotShowPBracketSize / 6,
                    tip.length = c(input$multiBoxplotShowPTipSizeA / 50, input$multiBoxplotShowPTipSizeB / 50),
                    color = colInput,
                    bracket.nudge.y = input$multiBoxplotShowPJust
                )
            } else {
                # Standard plotting if shortening is 0 (to stay safe)
                boxplot <- boxplot + stat_pvalue_manual(
                    data = dat,
                    y.position = "y.position",
                    label = "stars",
                    size = input$multiBoxplotShowPLabelSize * 2,
                    label.size = input$multiBoxplotShowPLabelSize * 2,
                    bracket.size = input$multiBoxplotShowPBracketSize / 6,
                    tip.length = c(input$multiBoxplotShowPTipSizeA / 50, input$multiBoxplotShowPTipSizeB / 50),
                    color = colInput,
                    bracket.nudge.y = input$multiBoxplotShowPJust
                )
            }
        }
    }

    boxplot <- boxplot + theme_prism(base_size = input$textSizeMultiBoxplot, axis_text_angle = input$multiBoxplotConditionAngle, base_family = input$multiBoxplotFont, border = input$multiBoxplotPlotBorder)

    if (!input$multiBoxplotShowConditions) {
      boxplot <- boxplot + theme(axis.text.x = element_blank())
    }

    multiBoxplotFiguresReactive$boxplot <- boxplot


    # Barplot logic can follow similar pattern if needed (omitted for brevity, focus on boxplot)
    # Just copying minimal barplot
    countsBarplot <- countsBarplot %>% tidyr::complete(!!sym(input$multiBoxplotFactor1), Feature, fill = list(Counts = 0))
    barplot <- ggplot(data = countsBarplot, aes(x = .data[[input$multiBoxplotFactor1]], y = Counts, fill = .data[[input$multiBoxplotFactor1]]))
    barplot <- barplot + facet_wrap(~Feature, scales = "free", ncol = input$multiBoxplotNCol)
    barplot <- barplot + geom_bar(stat = "identity", color = "black", position = position_dodge(), alpha = input$multiBoxplotBoxAlpha, show.legend = FALSE)
    barplot <- barplot + theme_prism(base_size = input$textSizeMultiBoxplot)

    multiBoxplotFiguresReactive$barplot <- barplot
  }
)

multiBoxplotFiguresReactive <- reactiveValues(boxplot = NULL, barplot = NULL)

lapply(c("boxplot", "barplot"), function(b) {
  output[[paste0("multi", toupper(substr(b, 1, 1)), substr(b, 2, nchar(b)), "Static")]] <- renderPlot(
    {
      req(!is.null(multiBoxplotFiguresReactive[[b]]))
      # Additional hygiene check at render time
      validate(
        need(length(multiGenesReactive$genes) > 0, "Please select at least one gene to plot.")
      )
      multiBoxplotFiguresReactive[[b]]
    },
    width = function() {
      as.numeric(input$figWidthMultiBoxplot)
    },
    height = function() {
      as.numeric(input$figHeightMultiBoxplot)
    }
  )
})

# Download Handlers
# Download Handlers
output$dlTable_multiBoxplotCounts <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_multiBoxplot, "_Counts"),
      design_name = "Multi",
      extension = "xlsx"
    )
  },
  content = function(file) {
    req(multiBoxplotCountsReactive()$countsBoxplot)
    d <- multiBoxplotCountsReactive()
    
    data_to_export <- d$countsBoxplot
    
    # If stats are requested
    if (input$includeDEStatsMultiBoxplot && !is.null(d$pvalueDF)) {
       wb <- openxlsx::createWorkbook()
       addWorksheet(wb, "Counts")
       writeDataTable(wb, "Counts", data_to_export, tableStyle = "TableStyleMedium9", withFilter = TRUE)
       
       addWorksheet(wb, "Statistics")
       writeDataTable(wb, "Statistics", d$pvalueDF, tableStyle = "TableStyleMedium9", withFilter = TRUE)
       
       saveWorkbook(wb, file, overwrite = TRUE)
    } else {
      export_excel_standard(
        df = data_to_export,
        filename = file,
        sheet_name = "Counts",
        table_name = "Counts"
      )
    }
  }
)

output$dlPlot_multiBoxplot <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_multiBoxplot, "_Boxplot"),
      design_name = "Multi",
      extension = tolower(input$downloadFormatMultiBoxplot)
    )
  },
  content = function(file) {
    req(multiBoxplotFiguresReactive$boxplot)
    export_plot_standard(
      plot_object = multiBoxplotFiguresReactive$boxplot,
      filename = file,
      format = tolower(input$downloadFormatMultiBoxplot),
      width = as.numeric(input$figWidthMultiBoxplot) / 75,
      height = as.numeric(input$figHeightMultiBoxplot) / 75,
      dpi = as.numeric(input$dpiMultiBoxplot)
    )
  }
)

output$dlPlot_multiBarplot <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_multiBoxplot, "_Barplot"),
      design_name = "Multi",
      extension = tolower(input$downloadFormatMultiBoxplot)
    )
  },
  content = function(file) {
    req(multiBoxplotFiguresReactive$barplot)
    export_plot_standard(
      plot_object = multiBoxplotFiguresReactive$barplot,
      filename = file,
      format = tolower(input$downloadFormatMultiBoxplot),
      width = as.numeric(input$figWidthMultiBoxplot) / 75,
      height = as.numeric(input$figHeightMultiBoxplot) / 75,
      dpi = as.numeric(input$dpiMultiBoxplot)
    )
  }
)
