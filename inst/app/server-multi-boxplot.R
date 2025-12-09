# Multi-Contrast Boxplot Server Logic
# This module creates boxplots showing expression data from multiple contrasts

# Define catalyst colour palette (from exploreDE) ----
catalystCols <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
    "#FF7B00", "#FDC362", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
    "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
    "#66A61E", "#E6F5C9", "#F781BF", "#FF00FF", "#FDDAEC", "#CCCCCC",
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"
)
catalystCols <- c(catalystCols[c(seq(1, 30, by = 2), seq(2, 30, by = 2))])
cc2 <- catalystCols
cc2 <- colorspace::darken(cc2, 0.4)
catalystCols <- c(catalystCols, cc2)
catalystCols <- paste0(catalystCols, "FF")

# Check compatibility of loaded contrasts ----
multiBoxplotCompatibilityReactive <- reactive({
    req(multiDeReactive$results)
    results <- multiDeReactive$results

    if (length(results) < 2) {
        return(list(
            "compatible" = FALSE,
            "message" = "You need at least 2 contrasts loaded to use this feature.",
            "samples" = NULL,
            "factors" = NULL,
            "param" = NULL
        ))
    }

    # Extract sample information and groupingName from each contrast
    sampleInfo <- lapply(results, function(r) {
        param <- r$param
        dataset <- r$dataset

        # param$groupingName is the column name (single string like "Condition")
        groupingName <- param$groupingName

        # Fallback if groupingName is missing (e.g. some Proteomics data)
        if (is.null(groupingName) && length(dataset) > 0) {
            # Try to guess or use first factor
            factors <- colnames(dataset)[grep("Sample|Name|names|channel|raw\\.file", colnames(dataset), invert = T, ignore.case = T)]
            if (length(factors) > 0) {
                groupingName <- factors[1]
            } else {
                # Fallback: Try to find a likely grouping column
                potential_groups <- colnames(dataset)[grep("condition|group|treatment|genotype|batch", colnames(dataset), ignore.case = TRUE)]
                if (length(potential_groups) > 0) {
                    groupingName <- potential_groups[1]
                } else {
                    groupingName <- "Name"
                }
            }
        }

        # param$grouping is a named vector: sample_name -> group_value
        groupingVector <- param$grouping

        # Get the dataset column data (sorted by sample Name for comparison)
        groupingData <- if (!is.null(groupingName) && groupingName %in% colnames(dataset)) {
            dataset[[groupingName]][order(dataset$Name)]
        } else {
            NULL
        }

        list(
            "samples" = sort(dataset$Name), # Sample names from dataset
            "nsamples" = nrow(dataset),
            "groupingName" = groupingName, # Column name string
            "groupingVector" = groupingVector, # Named vector from params
            "groupingData" = groupingData # Actual column values from dataset (sorted)
        )
    })

    # Check if all contrasts have the same samples
    firstSamples <- sampleInfo[[1]]$samples
    allSamesamples <- all(sapply(sampleInfo, function(x) {
        identical(sort(x$samples), sort(firstSamples))
    }))

    # Check if all contrasts use the same groupingName (column name)
    firstGroupingName <- sampleInfo[[1]]$groupingName

    # Validate groupingName is a single string
    if (is.null(firstGroupingName) || length(firstGroupingName) != 1 || !is.character(firstGroupingName)) {
        return(list(
            "compatible" = FALSE,
            "message" = "✗ Could not determine grouping factor name from first contrast",
            "samples" = NULL,
            "groupingName" = NULL,
            "dataset" = NULL
        ))
    }

    # Check all have same groupingName (string comparison)
    allSameGroupingName <- all(sapply(seq_along(sampleInfo), function(i) {
        thisName <- sampleInfo[[i]]$groupingName
        if (is.null(thisName) || length(thisName) != 1) {
            return(FALSE)
        }
        match <- identical(thisName, firstGroupingName)
        match
    }))

    # Check if the grouping column data is identical across all contrasts
    # Convert to character to avoid factor level mismatches
    firstGroupingData <- sampleInfo[[1]]$groupingData
    firstGroupingDataChar <- if (!is.null(firstGroupingData)) as.character(firstGroupingData) else NULL

    allSameGroupingData <- all(sapply(seq_along(sampleInfo), function(i) {
        if (is.null(sampleInfo[[i]]$groupingData) || is.null(firstGroupingData)) {
            return(FALSE)
        }

        # Convert to character for comparison to avoid factor level differences
        thisGroupingDataChar <- as.character(sampleInfo[[i]]$groupingData)
        match <- identical(thisGroupingDataChar, firstGroupingDataChar)
        match
    }))

    # Check if count dimensions match
    countInfo <- lapply(results, function(r) {
        se <- r$se
        list(
            "ngenes" = nrow(se),
            "nsamples" = ncol(se)
        )
    })

    allSameCounts <- all(sapply(countInfo, function(x) {
        x$ngenes == countInfo[[1]]$ngenes && x$nsamples == countInfo[[1]]$nsamples
    }))

    # All checks must pass
    compatible <- allSamesamples && allSameGroupingName && allSameGroupingData && allSameCounts

    if (compatible) {
        message <- paste0(
            "✓ All ", length(results), " contrasts are compatible! ",
            "They share ", length(firstSamples), " samples, ",
            countInfo[[1]]$ngenes, " genes, and use identical grouping (",
            firstGroupingName, " with same group assignments). You can proceed with multi-contrast boxplots."
        )
    } else {
        issues <- c()
        if (!allSamesamples) issues <- c(issues, "Different samples across contrasts")
        if (!allSameGroupingName) issues <- c(issues, "Different grouping factor names")
        if (!allSameGroupingData) issues <- c(issues, "Different group assignments (same factor name but different values per sample)")
        if (!allSameCounts) issues <- c(issues, "Different count matrix dimensions")
        message <- paste0(
            "✗ Contrasts are NOT compatible for multi-contrast boxplots. Issues: ",
            paste(issues, collapse = "; ")
        )
    }

    return(list(
        "compatible" = compatible,
        "message" = message,
        "samples" = if (compatible) firstSamples else NULL,
        "groupingName" = if (compatible) firstGroupingName else NULL,
        "dataset" = if (compatible) results[[1]]$dataset else NULL
    ))
})

# Render compatibility message ----
output$boxplotCompatibility <- renderUI({
    compat <- multiBoxplotCompatibilityReactive()

    if (compat$compatible) {
        tags$div(
            style = "padding: 10px;",
            tags$span(
                style = "color: #00AA00; font-size: 18px; font-weight: bold;",
                icon("check-circle"),
                " ", compat$message
            )
        )
    } else {
        tags$div(
            style = "padding: 10px;",
            tags$span(
                style = "color: #DD0000; font-size: 18px; font-weight: bold;",
                icon("times-circle"),
                " ", compat$message
            ),
            tags$p(
                style = "margin-top: 10px;",
                "Multi-contrast boxplots require all loaded contrasts to come from the same experiment
        with identical samples and count matrices. Please check your data."
            )
        )
    }
})

# Update UI elements based on compatibility ----
observe({
    req(multiDeReactive$results)
    req(length(multiDeReactive$results) > 0)

    compat <- multiBoxplotCompatibilityReactive()

    if (compat$compatible) {
        # Get gene list from first contrast
        genes <- multiDeReactive$results[[1]]$seqAnno$gene_name

        updateSelectizeInput(
            session = session,
            inputId = "multiBoxplotGenes",
            choices = genes,
            selected = "",
            server = TRUE
        )

        # Update count options based on data type
        # Assuming all loaded data is of the same type (enforced by multi-data-selection)
        firstType <- multiDeReactive$results[[1]]$dataType

        # Get count list from first result
        if (!is.null(multiDeReactive$results[[1]]$countList)) {
            choices <- names(multiDeReactive$results[[1]]$countList)
        } else {
            # Build count list from SE assays if countList not provided
            se <- multiDeReactive$results[[1]]$se
            availableAssays <- assayNames(se)
            choices <- c()

            if ("counts" %in% availableAssays) {
                choices <- c(choices, "Raw")
            }
            if ("xNorm" %in% availableAssays || "signal" %in% availableAssays || "counts" %in% availableAssays) {
                choices <- c(choices, "Normalised", "Normalised + Log2")
            }
            if ("counts" %in% availableAssays) {
                choices <- c(choices, "TPM", "FPKM", "VST")
            }
            if ("transformedData" %in% availableAssays) {
                choices <- c(choices, "transformedData")
            }
        }

        # Select default based on data type
        if (firstType == "proteomics") {
            selected <- "transformedData"
            if (!"transformedData" %in% choices && length(choices) > 0) selected <- choices[1]
        } else {
            if ("Normalised + Log2" %in% choices) {
                selected <- "Normalised + Log2"
            } else if (length(choices) > 0) {
                selected <- choices[1]
            }
        }

        updateSelectInput(session, "multiBoxplotCounts", choices = choices, selected = selected)
    }
})

# Create bucket list for sample selection ----
observe({
    req(multiBoxplotCompatibilityReactive()$compatible)

    compat <- multiBoxplotCompatibilityReactive()
    groupingName <- compat$groupingName

    req(!is.null(groupingName))
    req(groupingName %in% colnames(compat$dataset))

    groupLevels <- unique(compat$dataset[[groupingName]])

    output$multiBoxplotBucket <- renderUI({
        tagList(
            h5(paste0("Grouping factor: ", groupingName)),
            bucket_list(
                header = "Drag groups to include/exclude from plot (Ctrl/Cmd to select multiple)",
                group_name = "multi_boxplot_bucket_group",
                orientation = "horizontal",
                add_rank_list(
                    text = "Include these groups",
                    labels = as.list(groupLevels),
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
        )
    })
    outputOptions(output, "multiBoxplotBucket", suspendWhenHidden = FALSE)
})

# Create colour pickers for each group ----
output$multiBoxplotGroupColours <- renderUI({
    req(multiBoxplotCompatibilityReactive()$compatible)

    compat <- multiBoxplotCompatibilityReactive()
    groupingName <- compat$groupingName

    req(!is.null(groupingName))
    req(groupingName %in% colnames(compat$dataset))

    groupLevels <- unique(compat$dataset[[groupingName]])

    # Create colour pickers for each group with catalyst defaults
    colourPickers <- lapply(seq_along(groupLevels), function(i) {
        group <- groupLevels[i]
        inputId <- paste0("multiBoxplotGroupColour_", gsub(" |/|\\.", "_", group))

        # Get catalyst colour for this index (cycle if more groups than colours)
        defaultColour <- sub("FF$", "", catalystCols[(i - 1) %% length(catalystCols) + 1])

        colourpicker::colourInput(
            inputId = inputId,
            label = group,
            value = defaultColour,
            showColour = "both",
            palette = "limited",
            allowedCols = sub("FF$", "", catalystCols),
            returnName = FALSE
        )
    })

    do.call(tagList, colourPickers)
})

# Main reactive: prepare boxplot data ----
multiBoxplotDataReactive <- reactive({
    req(multiBoxplotCompatibilityReactive()$compatible)
    req(multiDeReactive$results)
    req(input$multiBoxKeepBucket)

    # Check if any genes are selected (either via dropdown or text input)
    hasGenesDropdown <- !is.null(input$multiBoxplotGenes) && length(input$multiBoxplotGenes) > 0
    hasGenesText <- !is.null(input$multiBoxplotGenesText) && nzchar(trimws(input$multiBoxplotGenesText))

    if (!hasGenesDropdown && !hasGenesText) {
        return(NULL)
    }

    results <- multiDeReactive$results
    compat <- multiBoxplotCompatibilityReactive()
    groupingName <- compat$groupingName

    req(!is.null(groupingName))

    # Get gene list
    genesList <- c(
        input$multiBoxplotGenes,
        if (!is.null(input$multiBoxplotGenesText)) {
            strsplit(input$multiBoxplotGenesText, "\\ +|\\/+|\\,+|\n+|\r+|\n\r+|\t+")[[1]]
        }
    )
    genesList <- unique(genesList[genesList != ""])

    # Limit to max genes
    if (length(genesList) > input$multiBoxplotMaxGenes) {
        genesList <- genesList[1:input$multiBoxplotMaxGenes]
    }

    req(length(genesList) >= 1)

    # Get counts from first contrast (they should all be identical)
    se <- results[[1]]$se
    param <- results[[1]]$param

    # Check if countList is already provided (e.g. for Proteomics)
    if (!is.null(results[[1]]$countList)) {
        countList <- results[[1]]$countList
    } else {
        # Build count list like exploreDE from SE
        availableAssays <- assayNames(se)
        countList <- list()

        # Raw counts
        if ("counts" %in% availableAssays) {
            countList[["Raw"]] <- assay(se, "counts")
        }

        # Normalized counts
        if ("xNorm" %in% availableAssays) {
            countList[["Normalised"]] <- assay(se, "xNorm")
        } else if ("signal" %in% availableAssays) {
            countList[["Normalised"]] <- assay(se, "signal")
        } else if ("counts" %in% availableAssays) {
            countList[["Normalised"]] <- assay(se, "counts")
        }

        # Normalised + Log2 (default)
        if (!is.null(countList[["Normalised"]])) {
            backgroundExp <- if (!is.null(param$backgroundExpression)) param$backgroundExpression else 1
            countList[["Normalised + Log2"]] <- log2(countList[["Normalised"]] + backgroundExp)
        }

        # TPM and FPKM (only if not proteomics)
        if ("counts" %in% availableAssays) {
            countList[["FPKM"]] <- getRpkm(se)
            countList[["TPM"]] <- getTpm(se)
        }

        # VST transformation
        if ("counts" %in% availableAssays) {
            forVST <- as.matrix(assay(se, "counts"))
            mode(forVST) <- "integer"
            vstCountMatrix <- DESeq2::vst(forVST, blind = FALSE)
            countList[["VST"]] <- vstCountMatrix
        }

        # Proteomics transformed data (fallback if not in countList)
        if ("transformedData" %in% availableAssays) {
            countList[["transformedData"]] <- assay(se, "transformedData")
        }
    }

    # Select the counts to use based on user choice
    # Update the choices in the UI first!
    updateSelectInput(session, "multiBoxplotCounts", choices = names(countList), selected = input$multiBoxplotCounts)

    if (input$multiBoxplotCounts %in% names(countList)) {
        counts <- countList[[input$multiBoxplotCounts]]
    } else {
        # Fallback
        if ("transformedData" %in% names(countList)) {
            counts <- countList[["transformedData"]]
        } else if ("Normalised + Log2" %in% names(countList)) {
            counts <- countList[["Normalised + Log2"]]
        } else if (length(countList) > 0) {
            counts <- countList[[1]]
        } else {
            return(NULL) # No counts available
        }
    }

    # Apply log2 transformation if requested (and data is not already log-transformed)
    if (input$multiBoxplotCountsLog && !grepl("Log2|VST|transformedData", input$multiBoxplotCounts)) {
        if (input$multiBoxplotCounts %in% c("TPM", "FPKM", "Normalised", "Raw")) {
            backgroundExp <- if (!is.null(param$backgroundExpression)) param$backgroundExpression else 1
            counts <- log2(counts + backgroundExp)
        }
    }

    # Match gene IDs to gene names from seqAnno
    seqAnno <- results[[1]]$seqAnno

    # Reorder counts to match seqAnno and use gene names as rownames
    if (param$featureLevel == "gene") {
        # Check if count rownames are gene_ids or gene_names
        matchByGeneId <- sum(rownames(counts) %in% seqAnno$gene_id)
        matchByGeneName <- sum(rownames(counts) %in% seqAnno$gene_name)

        if (matchByGeneId > matchByGeneName) {
            # Rownames are gene IDs - convert to gene names
            geneOrder <- intersect(seqAnno$gene_id, rownames(counts))
            counts <- counts[geneOrder, , drop = FALSE]
            rownames(counts) <- seqAnno$gene_name[match(geneOrder, seqAnno$gene_id)]
        } else {
            # Rownames are already gene names - just reorder to match seqAnno
            geneOrder <- intersect(seqAnno$gene_name, rownames(counts))
            counts <- counts[geneOrder, , drop = FALSE]
        }
    } else if (param$featureLevel == "isoform") {
        geneOrder <- intersect(seqAnno$transcript_id, rownames(counts))
        counts <- counts[geneOrder, , drop = FALSE]
        rownames(counts) <- seqAnno$gene_name[match(geneOrder, seqAnno$transcript_id)]
    }

    # Apply Z-score if requested
    if (input$multiBoxplotZScore) {
        counts <- t(scale(t(counts)))
    }

    # Get dataset metadata
    dataset <- compat$dataset

    # Keep only genes that exist in count matrix, but maintain the input order
    genesListOrdered <- genesList[genesList %in% rownames(counts)]

    if (length(genesListOrdered) == 0) {
        return(NULL)
    }

    req(length(genesListOrdered) >= 1)

    # Extract counts in the order specified by user
    counts <- counts[genesListOrdered, , drop = FALSE]

    # Update genesList to the ordered, filtered version
    genesList <- genesListOrdered

    # Transpose and merge with metadata
    countsDF <- as.data.frame(t(counts))
    countsDF$Name <- rownames(countsDF)

    countsDF <- left_join(countsDF, dataset, by = "Name")

    # Filter by selected groups
    countsDF <- countsDF[countsDF[[groupingName]] %in% input$multiBoxKeepBucket, ]

    # Set factor levels per bucket order
    countsDF[[groupingName]] <- factor(
        countsDF[[groupingName]],
        levels = input$multiBoxKeepBucket
    )

    # Get all metadata column names
    metadataCols <- colnames(dataset)[!colnames(dataset) %in% "Name"]

    # Melt the data
    countsMelt <- reshape2::melt(
        countsDF,
        id.vars = c("Name", metadataCols),
        value.name = "Counts",
        variable.name = "Feature"
    )
    countsMelt$Feature <- factor(countsMelt$Feature, levels = genesList)

    # Add DE results from all contrasts
    pvalueData <- list()
    for (i in seq_along(results)) {
        contrastName <- names(results)[i]
        seqAnno <- results[[i]]$seqAnno
        param <- results[[i]]$param

        # Extract the actual groups being compared from param
        # These are always levels in the grouping factor column
        refGroup <- param$refGroup
        sampleGroup <- param$sampleGroup

        if (is.null(refGroup) || is.null(sampleGroup)) {
            next
        }

        # Match genes to DE results
        for (gene in genesList) {
            geneIdx <- which(seqAnno$gene_name == gene)
            if (length(geneIdx) > 0) {
                pvalueData[[length(pvalueData) + 1]] <- data.frame(
                    Feature = gene,
                    Contrast = contrastName,
                    group1 = sampleGroup, # Sample group (numerator)
                    group2 = refGroup, # Reference group (denominator)
                    pValue = seqAnno$pValue[geneIdx],
                    fdr = seqAnno$fdr[geneIdx],
                    log2Ratio = seqAnno$log2Ratio[geneIdx],
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    pvalueDF <- if (length(pvalueData) > 0) do.call(rbind, pvalueData) else NULL

    # Add significance stars
    if (!is.null(pvalueDF) && nrow(pvalueDF) > 0) {
        pvalueDF$pval <- pvalueDF[[input$multiBoxplotPType]]
        pvalueDF <- pvalueDF %>%
            mutate(
                stars = case_when(
                    pval <= 0.0001 ~ "****",
                    pval <= 0.001 ~ "***",
                    pval <= 0.01 ~ "**",
                    pval <= 0.05 ~ "*",
                    TRUE ~ "ns"
                )
            )

        # Calculate y.position and x positions for brackets (stacked by contrast)
        # Get all group levels to calculate numeric positions
        groupLevels <- levels(countsDF[[groupingName]])

        yposData <- list()
        for (gene in genesList) {
            geneData <- countsMelt[countsMelt$Feature == gene, ]
            maxY <- max(geneData$Counts, na.rm = TRUE)
            minY <- min(geneData$Counts, na.rm = TRUE)
            yRange <- maxY - minY

            # Get contrasts for this gene
            geneContrasts <- pvalueDF$Contrast[pvalueDF$Feature == gene]
            nContrasts <- length(geneContrasts)

            # Base position (controlled by multiBoxplotPYOffset slider)
            baseYPos <- maxY + (yRange * input$multiBoxplotPYOffset / 100)

            for (i in seq_along(geneContrasts)) {
                # If stacking is enabled, offset each bracket vertically
                # If not stacking, all brackets at same height
                stackOffset <- if (input$multiBoxplotPStack) {
                    (i - 1) * (yRange * input$multiBoxplotPOffset / 100)
                } else {
                    0
                }

                # Get the group names for this contrast
                contrastData <- pvalueDF[pvalueDF$Contrast == geneContrasts[i] & pvalueDF$Feature == gene, ]
                group1_name <- contrastData$group1[1]
                group2_name <- contrastData$group2[1]

                # Find numeric positions of these groups (1-indexed)
                xmin_pos <- which(groupLevels == group1_name)
                xmax_pos <- which(groupLevels == group2_name)

                # Ensure xmin < xmax
                if (xmin_pos > xmax_pos) {
                    temp <- xmin_pos
                    xmin_pos <- xmax_pos
                    xmax_pos <- temp
                }

                # Apply bracket width adjustment
                # When width is 100%, use full span. When less, shrink towards center
                widthFactor <- input$multiBoxplotPWidth / 100
                bracketCenter <- (xmin_pos + xmax_pos) / 2
                bracketHalfSpan <- (xmax_pos - xmin_pos) / 2 * widthFactor

                xmin_adjusted <- bracketCenter - bracketHalfSpan
                xmax_adjusted <- bracketCenter + bracketHalfSpan

                yposData[[length(yposData) + 1]] <- data.frame(
                    Feature = gene,
                    Contrast = geneContrasts[i],
                    y.position = baseYPos + stackOffset,
                    xmin = xmin_adjusted,
                    xmax = xmax_adjusted,
                    stringsAsFactors = FALSE
                )
            }
        }
        yposDF <- do.call(rbind, yposData)

        # Merge y positions and x positions into pvalue data
        pvalueDF <- left_join(pvalueDF, yposDF, by = c("Feature", "Contrast"))
    }

    return(list(
        "countsMelt" = countsMelt,
        "pvalueDF" = pvalueDF,
        "genesList" = genesList
    ))
})

# Render the boxplot ----
output$boxplotPlot <- renderPlot(
    {
        req(multiBoxplotDataReactive())
        req(multiDeReactive$results)

        plotData <- multiBoxplotDataReactive()
        countsMelt <- plotData$countsMelt
        pvalueDF <- plotData$pvalueDF
        results <- multiDeReactive$results
        compat <- multiBoxplotCompatibilityReactive()
        groupingName <- compat$groupingName

        req(!is.null(groupingName))

        # Get colours for groups from colour pickers
        groupLevels <- levels(countsMelt[[groupingName]])
        colours <- setNames(
            sapply(groupLevels, function(k) {
                inputId <- paste0("multiBoxplotGroupColour_", gsub(" |/|\\.", "_", k))
                col <- input[[inputId]]
                if (!is.null(col)) {
                    paste0(gplots::col2hex(col), "FF")
                } else {
                    # Fallback to catalyst colours if colour picker not ready
                    idx <- match(k, groupLevels)
                    catalystCols[(idx - 1) %% length(catalystCols) + 1]
                }
            }),
            groupLevels
        )

        # Get contrast colours for p-value brackets
        contrastColours <- setNames(
            lapply(names(results), function(k) {
                col <- input[[paste0("multiGroupColour", gsub(" ", "_", k))]]
                if (!is.null(col)) {
                    paste0(gplots::col2hex(col))
                } else {
                    RColorBrewer::brewer.pal(max(3, length(results)), "Dark2")[match(k, names(results))]
                }
            }),
            names(results)
        )

        # Create base plot
        p <- ggplot(countsMelt, aes(x = .data[[groupingName]], y = Counts))

        # Add boxplot
        if (input$multiBoxplotShowBox && !input$multiBoxplotShowViolin) {
            p <- p + geom_boxplot(
                aes(fill = .data[[groupingName]]),
                alpha = input$multiBoxplotBoxAlpha,
                outlier.shape = NA,
                show.legend = FALSE
            )
        }

        # Add violin plot
        if (!input$multiBoxplotShowBox && input$multiBoxplotShowViolin) {
            p <- p + geom_violin(
                aes(fill = .data[[groupingName]]),
                alpha = input$multiBoxplotBoxAlpha,
                draw_quantiles = 0.5,
                show.legend = FALSE
            )
        }

        # Add mean bar (behind other elements)
        if (input$multiBoxplotShowMeanBar && !input$multiBoxplotMeanBarFront) {
            p <- p + stat_summary(
                geom = "crossbar",
                fun = mean,
                position = position_dodge(0.9),
                colour = "black",
                size = 0.4,
                width = input$multiBoxplotMeanLine,
                show.legend = FALSE
            )
        }

        # Add points (using beeswarm for better distribution)
        if (input$multiBoxplotShowPoint) {
            p <- p + geom_beeswarm(
                aes(fill = .data[[groupingName]]),
                size = input$multiBoxplotPointSize,
                alpha = input$multiBoxplotPointAlpha,
                shape = 21,
                cex = 2,
                method = "swarm",
                corral = "none",
                corral.width = 0.9,
                stroke = 0.5
            )
        }

        # Add mean bar (in front of other elements)
        if (input$multiBoxplotShowMeanBar && input$multiBoxplotMeanBarFront) {
            p <- p + stat_summary(
                geom = "crossbar",
                fun = mean,
                position = position_dodge(0.9),
                colour = "black",
                size = 0.4,
                width = input$multiBoxplotMeanLine,
                show.legend = FALSE
            )
        }

        # Apply colours
        p <- p + scale_fill_manual(values = colours)

        # Add p-value brackets if requested
        if (input$multiBoxplotShowP && !is.null(pvalueDF) && nrow(pvalueDF) > 0) {
            # Add brackets for each contrast
            # Use black if requested, otherwise use contrast colours
            for (contrast in names(results)) {
                contrastData <- pvalueDF[pvalueDF$Contrast == contrast, ]
                if (nrow(contrastData) > 0) {
                    bracketColour <- if (input$multiBoxplotPBlack) {
                        "black"
                    } else {
                        as.character(contrastColours[[contrast]])
                    }

                    p <- p + stat_pvalue_manual(
                        data = contrastData,
                        label = "stars",
                        y.position = "y.position",
                        xmin = "xmin",
                        xmax = "xmax",
                        color = bracketColour,
                        size = input$multiBoxplotPLabelSize,
                        label.size = input$multiBoxplotPLabelSize,
                        bracket.size = input$multiBoxplotPBracketSize,
                        tip.length = input$multiBoxplotPTipSize,
                        inherit.aes = FALSE
                    )
                }
            }
            p <- p + coord_cartesian(clip = "off")
        }

        # Facet by gene
        p <- p + facet_wrap(~Feature, scales = "free_y", ncol = input$multiBoxplotNCol)

        # Labels
        ylab <- paste(input$multiBoxplotCounts, "Counts")
        if (input$multiBoxplotCountsLog && !grepl("Log2|VST", input$multiBoxplotCounts)) {
            ylab <- paste(ylab, "(Log2)")
        }
        if (input$multiBoxplotZScore) {
            ylab <- paste(ylab, "[Z-score]")
        }
        p <- p + labs(y = ylab, x = groupingName)

        # Theme - match exploreDE style
        p <- p + theme_prism(base_size = input$multiTextSizeBoxplot, axis_text_angle = input$multiBoxplotConditionAngle)

        print(p)
    },
    width = function() input$multiFigWidthBoxplot,
    height = function() input$multiFigHeightBoxplot
)

# Download handlers ----
output$dlBoxplot <- downloadHandler(
    filename = function() {
        paste0("multi_contrast_boxplot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
        plotData <- multiBoxplotDataReactive()
        req(plotData)

        # Recreate plot (copy logic from renderPlot)
        countsMelt <- plotData$countsMelt
        pvalueDF <- plotData$pvalueDF
        results <- multiDeReactive$results
        compat <- multiBoxplotCompatibilityReactive()
        groupingName <- compat$groupingName

        # Get colours for groups from colour pickers
        groupLevels <- levels(countsMelt[[groupingName]])
        colours <- setNames(
            sapply(groupLevels, function(k) {
                inputId <- paste0("multiBoxplotGroupColour_", gsub(" |/|\\.", "_", k))
                col <- input[[inputId]]
                if (!is.null(col)) {
                    paste0(gplots::col2hex(col), "FF")
                } else {
                    # Fallback to catalyst colours if colour picker not ready
                    idx <- match(k, groupLevels)
                    catalystCols[(idx - 1) %% length(catalystCols) + 1]
                }
            }),
            groupLevels
        )

        # Get contrast colours for p-value brackets
        contrastColours <- setNames(
            lapply(names(results), function(k) {
                col <- input[[paste0("multiGroupColour", gsub(" ", "_", k))]]
                if (!is.null(col)) {
                    paste0(gplots::col2hex(col))
                } else {
                    RColorBrewer::brewer.pal(max(3, length(results)), "Dark2")[match(k, names(results))]
                }
            }),
            names(results)
        )

        # Create base plot
        p <- ggplot(countsMelt, aes(x = .data[[groupingName]], y = Counts))

        # Add boxplot
        if (input$multiBoxplotShowBox && !input$multiBoxplotShowViolin) {
            p <- p + geom_boxplot(
                aes(fill = .data[[groupingName]]),
                alpha = input$multiBoxplotBoxAlpha,
                outlier.shape = NA,
                show.legend = FALSE
            )
        }

        # Add violin
        if (!input$multiBoxplotShowBox && input$multiBoxplotShowViolin) {
            p <- p + geom_violin(
                aes(fill = .data[[groupingName]]),
                alpha = input$multiBoxplotBoxAlpha,
                draw_quantiles = 0.5,
                show.legend = FALSE
            )
        }

        # Add mean bar (behind other elements)
        if (input$multiBoxplotShowMeanBar && !input$multiBoxplotMeanBarFront) {
            p <- p + stat_summary(
                geom = "crossbar",
                fun = mean,
                position = position_dodge(0.9),
                colour = "black",
                size = 0.4,
                width = input$multiBoxplotMeanLine,
                show.legend = FALSE
            )
        }

        # Add points (using beeswarm for better distribution)
        if (input$multiBoxplotShowPoint) {
            p <- p + geom_beeswarm(
                aes(fill = .data[[groupingName]]),
                size = input$multiBoxplotPointSize,
                alpha = input$multiBoxplotPointAlpha,
                shape = 21,
                cex = 2,
                method = "swarm",
                corral = "none",
                corral.width = 0.9,
                stroke = 0.5
            )
        }

        # Add mean bar (in front of other elements)
        if (input$multiBoxplotShowMeanBar && input$multiBoxplotMeanBarFront) {
            p <- p + stat_summary(
                geom = "crossbar",
                fun = mean,
                position = position_dodge(0.9),
                colour = "black",
                size = 0.4,
                width = input$multiBoxplotMeanLine,
                show.legend = FALSE
            )
        }

        # Apply colours
        p <- p + scale_fill_manual(values = colours)

        # Add p-value brackets if requested
        if (input$multiBoxplotShowP && !is.null(pvalueDF) && nrow(pvalueDF) > 0) {
            # Add brackets for each contrast
            # Use black if requested, otherwise use contrast colours
            for (contrast in names(results)) {
                contrastData <- pvalueDF[pvalueDF$Contrast == contrast, ]
                if (nrow(contrastData) > 0) {
                    bracketColour <- if (input$multiBoxplotPBlack) {
                        "black"
                    } else {
                        as.character(contrastColours[[contrast]])
                    }

                    p <- p + stat_pvalue_manual(
                        data = contrastData,
                        label = "stars",
                        y.position = "y.position",
                        xmin = "xmin",
                        xmax = "xmax",
                        color = bracketColour,
                        size = input$multiBoxplotPLabelSize,
                        label.size = input$multiBoxplotPLabelSize,
                        bracket.size = input$multiBoxplotPBracketSize,
                        tip.length = input$multiBoxplotPTipSize,
                        inherit.aes = FALSE
                    )
                }
            }
            p <- p + coord_cartesian(clip = "off")
        }

        # Facet by gene
        p <- p + facet_wrap(~Feature, scales = "free_y", ncol = input$multiBoxplotNCol)

        # Labels
        ylab <- paste(input$multiBoxplotCounts, "Counts")
        if (input$multiBoxplotCountsLog && !grepl("Log2|VST", input$multiBoxplotCounts)) {
            ylab <- paste(ylab, "(Log2)")
        }
        if (input$multiBoxplotZScore) {
            ylab <- paste(ylab, "[Z-score]")
        }
        p <- p + labs(y = ylab, x = groupingName)

        # Theme - match exploreDE style
        p <- p + theme_prism(base_size = input$multiTextSizeBoxplot, axis_text_angle = input$multiBoxplotConditionAngle)

        pdf(file, width = input$multiFigWidthBoxplot / 72, height = input$multiFigHeightBoxplot / 72)
        print(p)
        dev.off()
    }
)

# Download Excel table ----
output$dlBoxplotTable <- downloadHandler(
    filename = function() {
        paste0("multi_contrast_boxplot_data_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
        req(multiBoxplotDataReactive())
        plotData <- multiBoxplotDataReactive()
        writexl::write_xlsx(plotData$countsMelt, file)
    }
)

# Render data table ----
output$boxplotTable <- renderDataTable({
    req(multiBoxplotDataReactive())
    req(multiDeReactive$results)

    plotData <- multiBoxplotDataReactive()
    results <- multiDeReactive$results
    genesList <- plotData$genesList
    countsMelt <- plotData$countsMelt
    compat <- multiBoxplotCompatibilityReactive()
    groupingName <- compat$groupingName

    # Build summary table: one row per gene, columns for each contrast
    tableData <- lapply(genesList, function(gene) {
        # Start with gene name
        row <- data.frame(gene_name = gene, stringsAsFactors = FALSE)

        # Add DE stats for each contrast
        for (contrastName in names(results)) {
            seqAnno <- results[[contrastName]]$seqAnno

            # Determine which column to use for matching
            matchCol <- if ("gene_name" %in% colnames(seqAnno)) {
                "gene_name"
            } else if ("gene_id" %in% colnames(seqAnno)) {
                "gene_id"
            } else {
                NULL
            }

            # Find this gene in seqAnno
            geneRow <- if (!is.null(matchCol)) {
                seqAnno[seqAnno[[matchCol]] == gene, ]
            } else {
                data.frame()
            }

            if (nrow(geneRow) > 0) {
                # Extract log2FC
                lfcCol <- grep("^(log2?[._ ]?(fold)?[._ ]?(change|ratio|fc)|logfc|log2Ratio)$",
                    colnames(seqAnno),
                    value = TRUE, ignore.case = TRUE
                )
                row[[paste0(contrastName, "_log2FC")]] <- if (length(lfcCol) > 0) geneRow[[lfcCol[1]]] else NA

                # Extract p-value
                pvalCol <- grep("^(pvalue|p[._]?val|pval)$", colnames(seqAnno), value = TRUE, ignore.case = TRUE)
                row[[paste0(contrastName, "_pvalue")]] <- if (length(pvalCol) > 0) geneRow[[pvalCol[1]]] else NA

                # Extract padj
                padjCol <- grep("^(padj|fdr|p[._]?adj|adj[._]?p)$", colnames(seqAnno), value = TRUE, ignore.case = TRUE)
                row[[paste0(contrastName, "_padj")]] <- if (length(padjCol) > 0) geneRow[[padjCol[1]]] else NA
            } else {
                row[[paste0(contrastName, "_log2FC")]] <- NA
                row[[paste0(contrastName, "_pvalue")]] <- NA
                row[[paste0(contrastName, "_padj")]] <- NA
            }
        }

        # Calculate mean counts per group for this gene
        geneCounts <- countsMelt[countsMelt$Feature == gene, ]
        meansByGroup <- aggregate(Counts ~ get(groupingName), data = geneCounts, FUN = mean)
        colnames(meansByGroup) <- c("Group", "MeanCount")

        # Add mean counts as columns
        for (i in seq_len(nrow(meansByGroup))) {
            row[[paste0("Mean_", meansByGroup$Group[i])]] <- meansByGroup$MeanCount[i]
        }

        return(row)
    })

    tableData <- bind_rows(tableData)

    # Find column names for formatting
    log2fc_cols <- grep("_log2FC$", colnames(tableData), value = TRUE)
    pval_cols <- grep("_pvalue$", colnames(tableData), value = TRUE)
    padj_cols <- grep("_padj$", colnames(tableData), value = TRUE)
    mean_cols <- grep("^Mean_", colnames(tableData), value = TRUE)

    # Render table
    dt <- DT::datatable(
        tableData,
        filter = "top",
        class = "cell-border stripe",
        rownames = FALSE,
        options = list(
            scrollX = TRUE,
            pageLength = 25,
            autoWidth = TRUE
        )
    )

    # Format numeric columns
    if (length(log2fc_cols) > 0) {
        dt <- dt %>% DT::formatSignif(columns = log2fc_cols, digits = 3)
    }
    if (length(pval_cols) > 0) {
        dt <- dt %>% DT::formatSignif(columns = pval_cols, digits = 4)
    }
    if (length(padj_cols) > 0) {
        dt <- dt %>% DT::formatSignif(columns = padj_cols, digits = 4)
    }
    if (length(mean_cols) > 0) {
        dt <- dt %>% DT::formatSignif(columns = mean_cols, digits = 2)
    }

    # Apply styling
    if (length(pval_cols) > 0) {
        dt <- dt %>% DT::formatStyle(
            columns = pval_cols,
            color = styleInterval(cuts = 0.05, values = c("green", "black")),
            fontWeight = "bold"
        )
    }
    if (length(padj_cols) > 0) {
        dt <- dt %>% DT::formatStyle(
            columns = padj_cols,
            color = styleInterval(cuts = 0.05, values = c("green", "black")),
            fontWeight = "bold"
        )
    }
    if (length(log2fc_cols) > 0) {
        dt <- dt %>% DT::formatStyle(
            columns = log2fc_cols,
            color = styleInterval(cuts = 0, values = c("blue", "orange")),
            fontWeight = "bold"
        )
    }

    dt
})
