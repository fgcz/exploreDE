# Update inputs
observeEvent(inputDataReactive(), {
    count_choices <- names(inputDataReactive()$countList)
    default_count <- switch(inputDataReactive()$dataType,
        "RNASeq" = "Normalised + Log2",
        "proteomics" = if ("transformedData" %in% count_choices) "transformedData" else count_choices[1],
        "CountQC" = "VST"
    ) # Default to VST for CountQC if available, or whatever is best

    if (inputDataReactive()$dataType == "CountQC" && "VST" %in% count_choices) {
        default_count <- "VST"
    }

    updateSelectInput(session = session, inputId = "replicabilityCounts", choices = count_choices, selected = default_count)
})

output$replicabilityFactorSelectUI <- renderUI({
    selectizeInput(
        inputId = "replicabilityFactors",
        label = "Select factors for annotation",
        choices = inputDataReactive()$factors,
        multiple = TRUE,
        selected = inputDataReactive()$factors[1]
    )
})

# Update break values based on correlation data
observeEvent(replicabilityData(), {
    req(replicabilityData())

    # Get all correlation data
    all_corr <- replicabilityData()$all$corr
    top100_corr <- replicabilityData()$top100$corr

    # Use the min/mean/max from all genes
    corr_min <- min(all_corr, na.rm = TRUE)
    corr_mean <- mean(all_corr, na.rm = TRUE)
    corr_max <- max(all_corr, na.rm = TRUE)

    # Update the numeric inputs
    updateNumericInput(session, "replicabilityAtLow", value = round(corr_min, 2))
    updateNumericInput(session, "replicabilityAtMid", value = round(corr_mean, 2))
    updateNumericInput(session, "replicabilityAtHigh", value = round(corr_max, 2))
})

# Reactive for calculating correlation and clustering
replicabilityData <- reactive({
    req(input$replicabilityCounts)
    req(inputDataReactive())

    counts <- inputDataReactive()$countList[[input$replicabilityCounts]]

    # Apply Log2 if requested
    if (isTRUE(input$replicabilityCountsLog)) {
        counts <- log2(counts + 1)
    }

    # Helper to calculate corr and hclust
    calc_corr_hclust <- function(cnts) {
        # Correlation
        # Use pairwise.complete.obs to handle NAs (e.g. in proteomics)
        corr_mat <- cor(cnts, method = input$replicabilityCorrMethod, use = "pairwise.complete.obs")

        # Clustering
        # Transpose counts for dist (samples as rows)
        # Handle NAs in dist if necessary, but dist usually propagates NAs or requires imputation.
        # For now, we rely on standard dist behavior but could impute if needed.
        # Replacing NAs with 0 for clustering might be safer for proteomics if NAs are missing values.
        cnts_for_dist <- cnts
        cnts_for_dist[is.na(cnts_for_dist)] <- 0

        dist_mat <- dist(t(cnts_for_dist), method = input$replicabilityDistMethod)
        hclust_res <- hclust(dist_mat, method = input$replicabilityClustMethod)
        list(corr = corr_mat, hclust = hclust_res)
    }

    # 1. All Genes (filtered by variance)
    counts_all <- counts
    if (nrow(counts_all) > input$replicabilityTopGenes) {
        vars <- apply(counts_all, 1, var, na.rm = TRUE)
        counts_all <- counts_all[order(vars, decreasing = TRUE)[1:input$replicabilityTopGenes], ]
    }
    res_all <- calc_corr_hclust(counts_all)

    # 2. Top 100 Abundant Genes
    counts_top100 <- counts
    means <- rowMeans(counts_top100, na.rm = TRUE)
    # Take top 100 most abundant
    top_n <- min(100, nrow(counts_top100))
    counts_top100 <- counts_top100[order(means, decreasing = TRUE)[1:top_n], ]
    res_top100 <- calc_corr_hclust(counts_top100)

    list(
        all = list(corr = res_all$corr, hclust = res_all$hclust, dataset = inputDataReactive()$dataset),
        top100 = list(corr = res_top100$corr, hclust = res_top100$hclust, dataset = inputDataReactive()$dataset)
    )
})

# Render Heatmap
# Helper function for Heatmap
make_rep_heatmap <- function(data, title_suffix = "") {
    req(data)

    # Prepare annotations
    ha_top <- NULL
    ha_left <- NULL
    row_names_gp <- gpar(fontsize = input$replicabilityTextSize)
    column_names_gp <- gpar(fontsize = input$replicabilityTextSize)

    if (!is.null(input$replicabilityFactors)) {
        dataset <- data$dataset
        sample_names <- rownames(data$corr)
        dataset <- dataset[match(sample_names, dataset$Name), ]

        # Get colors
        colourList <- lapply(input$replicabilityFactors, function(f) {
            levs <- levels(as.factor(dataset[[f]]))
            colors <- vapply(levs, function(l) {
                col_input <- input[[paste0("GroupColour", f, l)]]
                if (is.null(col_input)) {
                    return("#808080FF")
                }
                tryCatch(paste0(col2hex(col_input), "FF"), error = function(e) "#808080FF")
            }, character(1))
            names(colors) <- levs
            return(colors)
        }) %>% setNames(input$replicabilityFactors)

        # Color labels based on the first selected factor
        factor_col <- input$replicabilityFactors[1]
        factor_levels <- levels(as.factor(dataset[[factor_col]]))
        colors <- vapply(factor_levels, function(l) {
            col_input <- input[[paste0("GroupColour", factor_col, l)]]
            if (is.null(col_input)) {
                return("#808080FF")
            }
            tryCatch(paste0(col2hex(col_input), "FF"), error = function(e) "#808080FF")
        }, character(1))

        sample_colors <- colors[as.character(dataset[[factor_col]])]
        names(sample_colors) <- dataset$Name

        row_names_gp <- gpar(fontsize = input$replicabilityTextSize, col = sample_colors)
        column_names_gp <- gpar(fontsize = input$replicabilityTextSize, col = sample_colors)

        anno_df <- dataset[, input$replicabilityFactors, drop = FALSE]

        ha_top <- HeatmapAnnotation(
            df = anno_df,
            col = colourList,
            show_legend = TRUE,
            annotation_name_gp = gpar(fontsize = input$replicabilityTextSize, fontface = "bold"),
            gp = gpar(col = input$replicabilityBorderColor1, lwd = input$replicabilityBorderSize1 / 10)
        )

        ha_left <- rowAnnotation(
            df = anno_df,
            col = colourList,
            show_legend = FALSE,
            annotation_name_gp = gpar(fontsize = input$replicabilityTextSize, fontface = "bold"),
            gp = gpar(col = input$replicabilityBorderColor1, lwd = input$replicabilityBorderSize1 / 10)
        )
    }

    # Get color breaks
    at <- c(input$replicabilityAtLow, input$replicabilityAtMid, input$replicabilityAtHigh)

    # Determine color scheme - use palette if selected, otherwise use manual colors
    if (!is.null(input$replicabilityPalette) && input$replicabilityPalette != "None") {
        palette_colours <- brewer.pal(9, input$replicabilityPalette)
        selected_colours <- palette_colours[c(1, 5, 9)]
        if (input$replicabilityPaletteRev) {
            heatmap_colours <- colorRamp2(at, selected_colours)
        } else {
            heatmap_colours <- colorRamp2(at, rev(selected_colours))
        }
    } else {
        heatmap_colours <- colorRamp2(at, c(input$replicabilityColorLow, input$replicabilityColorMid, input$replicabilityColorHigh))
    }

    hm <- ComplexHeatmap::Heatmap(
        data$corr,
        name = "Correlation",
        col = heatmap_colours,
        top_annotation = ha_top,
        left_annotation = ha_left,
        show_row_names = input$replicabilityShowNames,
        show_column_names = input$replicabilityShowNames,
        row_names_gp = row_names_gp,
        column_names_gp = column_names_gp,
        heatmap_legend_param = list(title_gp = gpar(fontsize = input$replicabilityTextSize, fontface = "bold")),
        rect_gp = gpar(col = input$replicabilityBorderColor2, lwd = input$replicabilityBorderSize2 / 10),
        column_title = paste0("Correlation Heatmap: ", title_suffix),
        column_title_gp = gpar(fontsize = input$replicabilityTextSize + 2, fontface = "bold")
    )

    draw(hm, merge_legend = TRUE, padding = unit(c(10, 10, 10, 10), "mm"))
}

# Helper function for Dendrogram
make_rep_dendrogram <- function(data, title_suffix = "") {
    req(data)
    dend <- as.dendrogram(data$hclust)

    if (!is.null(input$replicabilityFactors)) {
        dataset <- data$dataset
        factor_col <- input$replicabilityFactors[1]
        factor_levels <- levels(as.factor(dataset[[factor_col]]))
        colors <- vapply(factor_levels, function(l) {
            col_input <- input[[paste0("GroupColour", factor_col, l)]]
            if (is.null(col_input)) {
                return("#808080FF")
            }
            tryCatch(paste0(col2hex(col_input), "FF"), error = function(e) "#808080FF")
        }, character(1))

        sample_colors <- colors[as.character(dataset[[factor_col]])]
        names(sample_colors) <- dataset$Name
        labels_colors(dend) <- sample_colors[labels(dend)]

        colors_mat <- as.data.frame(dataset[, input$replicabilityFactors, drop = FALSE])
        for (f in input$replicabilityFactors) {
            colors_mat[[f]] <- vapply(colors_mat[[f]], function(l) {
                col_input <- input[[paste0("GroupColour", f, l)]]
                if (is.null(col_input)) {
                    return("#808080FF")
                }
                tryCatch(paste0(col2hex(col_input), "FF"), error = function(e) "#808080FF")
            }, character(1))
        }
        colors_mat <- colors_mat[match(labels(dend), dataset$Name), , drop = FALSE]
    }

    par(mar = c(20, 8, 6, 4))
    plot(dend,
        main = paste0("Cluster Dendrogram (", input$replicabilityDistMethod, "/", input$replicabilityClustMethod, "): ", title_suffix),
        xlab = "", sub = "", cex = input$replicabilityTextSize / 10
    )

    if (!is.null(input$replicabilityFactors)) {
        colored_bars(colors = colors_mat, dend = dend, sort_by_labels_order = FALSE)
    }
}

# Render Heatmaps
output$replicabilityHeatmap <- renderPlot(
    {
        make_rep_heatmap(replicabilityData()$all, title_suffix = "All Genes")
    },
    width = function() input$figWidthReplicability,
    height = function() input$figHeightReplicability
)

output$replicabilityHeatmapTop100 <- renderPlot(
    {
        make_rep_heatmap(replicabilityData()$top100, title_suffix = "Top 100 Abundant Genes")
    },
    width = function() input$figWidthReplicability,
    height = function() input$figHeightReplicability
)

# Render Dendrograms
output$replicabilityDendrogram <- renderPlot(
    {
        make_rep_dendrogram(replicabilityData()$all, title_suffix = "All Genes")
    },
    width = function() input$figWidthReplicability,
    height = function() input$figHeightReplicability
)

output$replicabilityDendrogramTop100 <- renderPlot(
    {
        make_rep_dendrogram(replicabilityData()$top100, title_suffix = "Top 100 Abundant Genes")
    },
    width = function() input$figWidthReplicability,
    height = function() input$figHeightReplicability
)

# Download Handlers
output$replicabilityDownloadHeatmap <- downloadHandler(
    filename = function() {
        paste0("Replicability_Heatmap.pdf")
    },
    content = function(file) {
        req(replicabilityData())
        data <- replicabilityData()

        # Re-create heatmap logic for download
        # (Ideally this would be a reactive, but for now we duplicate logic for simplicity)
        # Prepare annotations
        ha_top <- NULL
        ha_left <- NULL
        row_names_gp <- gpar(fontsize = input$replicabilityTextSize)
        column_names_gp <- gpar(fontsize = input$replicabilityTextSize)

        if (!is.null(input$replicabilityFactors)) {
            dataset <- data$dataset
            sample_names <- rownames(data$corr)
            dataset <- dataset[match(sample_names, dataset$Name), ]

            colourList <- lapply(input$replicabilityFactors, function(f) {
                levs <- levels(as.factor(dataset[[f]]))
                colors <- vapply(levs, function(l) {
                    col_input <- input[[paste0("GroupColour", f, l)]]
                    if (is.null(col_input)) {
                        return("#808080FF")
                    }
                    tryCatch(paste0(col2hex(col_input), "FF"), error = function(e) "#808080FF")
                }, character(1))
                names(colors) <- levs
                return(colors)
            }) %>% setNames(input$replicabilityFactors)

            # Color labels
            factor_col <- input$replicabilityFactors[1]
            factor_levels <- levels(dataset[[factor_col]])
            colors <- sapply(factor_levels, function(l) {
                paste0(col2hex(input[[paste0("GroupColour", factor_col, l)]]), "FF")
            })
            sample_colors <- colors[as.character(dataset[[factor_col]])]
            names(sample_colors) <- dataset$Name
            row_names_gp <- gpar(fontsize = input$replicabilityTextSize, col = sample_colors)
            column_names_gp <- gpar(fontsize = input$replicabilityTextSize, col = sample_colors)

            anno_df <- dataset[, input$replicabilityFactors, drop = FALSE]

            ha_top <- HeatmapAnnotation(
                df = anno_df,
                col = colourList,
                show_legend = TRUE,
                annotation_name_gp = gpar(fontsize = input$replicabilityTextSize, fontface = "bold"),
                gp = gpar(col = input$replicabilityBorderColor1, lwd = input$replicabilityBorderSize1 / 10)
            )
            ha_left <- rowAnnotation(
                df = anno_df,
                col = colourList,
                show_legend = FALSE,
                annotation_name_gp = gpar(fontsize = input$replicabilityTextSize, fontface = "bold"),
                gp = gpar(col = input$replicabilityBorderColor1, lwd = input$replicabilityBorderSize1 / 10)
            )
        }

        # Get color breaks
        at <- c(input$replicabilityAtLow, input$replicabilityAtMid, input$replicabilityAtHigh)

        # Determine color scheme - use palette if selected, otherwise use manual colors
        if (!is.null(input$replicabilityPalette) && input$replicabilityPalette != "None") {
            palette_colours <- brewer.pal(9, input$replicabilityPalette)
            selected_colours <- palette_colours[c(1, 5, 9)]
            if (input$replicabilityPaletteRev) {
                heatmap_colours <- colorRamp2(at, selected_colours)
            } else {
                heatmap_colours <- colorRamp2(at, rev(selected_colours))
            }
        } else {
            heatmap_colours <- colorRamp2(at, c(input$replicabilityColorLow, input$replicabilityColorMid, input$replicabilityColorHigh))
        }

        hm <- ComplexHeatmap::Heatmap(
            data$corr,
            name = "Correlation",
            col = heatmap_colours,
            top_annotation = ha_top,
            left_annotation = ha_left,
            show_row_names = input$replicabilityShowNames,
            show_column_names = input$replicabilityShowNames,
            row_names_gp = row_names_gp,
            column_names_gp = column_names_gp,
            heatmap_legend_param = list(title_gp = gpar(fontsize = input$replicabilityTextSize, fontface = "bold")),
            rect_gp = gpar(col = input$replicabilityBorderColor2, lwd = input$replicabilityBorderSize2 / 10)
        )

        pdf(file, width = input$figWidthReplicability / 72, height = input$figHeightReplicability / 72)
        draw(hm, merge_legend = TRUE, padding = unit(c(10, 10, 10, 10), "mm"))
        dev.off()
    }
)

output$replicabilityDownloadDendrogram <- downloadHandler(
    filename = function() {
        paste0("Replicability_Dendrogram.pdf")
    },
    content = function(file) {
        req(replicabilityData())
        data <- replicabilityData()

        pdf(file, width = input$figWidthReplicability / 72, height = input$figHeightReplicability / 72)

        dend <- as.dendrogram(data$hclust)
        if (!is.null(input$replicabilityFactors)) {
            dataset <- data$dataset
            factor_col <- input$replicabilityFactors[1]
            factor_levels <- levels(dataset[[factor_col]])
            colors <- sapply(factor_levels, function(l) {
                paste0(col2hex(input[[paste0("GroupColour", factor_col, l)]]), "FF")
            })
            sample_colors <- colors[as.character(dataset[[factor_col]])]
            names(sample_colors) <- dataset$Name
            labels_colors(dend) <- sample_colors[labels(dend)]

            colors_mat <- as.data.frame(dataset[, input$replicabilityFactors, drop = FALSE])
            for (f in input$replicabilityFactors) {
                colors_mat[[f]] <- sapply(colors_mat[[f]], function(l) {
                    paste0(col2hex(input[[paste0("GroupColour", f, l)]]), "FF")
                })
            }
            colors_mat <- colors_mat[match(labels(dend), dataset$Name), , drop = FALSE]
        }

        par(mar = c(20, 8, 6, 4))
        plot(dend,
            main = paste("Cluster Dendrogram (", input$replicabilityDistMethod, "/", input$replicabilityClustMethod, ")"),
            xlab = "", sub = "", cex = input$replicabilityTextSize / 10
        )

        if (!is.null(input$replicabilityFactors)) {
            colored_bars(colors = colors_mat, dend = dend, sort_by_labels_order = FALSE)
        }

        dev.off()
    }
)

output$replicabilityDownloadCorr <- downloadHandler(
    filename = function() {
        paste0("Replicability_Correlation_Matrix.xlsx")
    },
    content = function(file) {
        req(replicabilityData())
        data <- replicabilityData()
        write.xlsx(as.data.frame(data$corr), file, rowNames = TRUE)
    }
)
