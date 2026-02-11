# Multi-Contrast Tile Plot Server Logic

tilePlots <- reactiveValues(top = NULL, custom = NULL, go = NULL, mat1 = NULL, df3 = NULL, custom_df3 = NULL, go_df3 = NULL)

multiMatReactive <- eventReactive(multiDeReactive$results, {
    req(multiDeReactive$results)
    results <- multiDeReactive$results

    dfList <- list()
    for (n in seq_along(names(results))) {
        # Check if seqAnno exists directly (standard format) or need to extract from seqAnnoList (proteomics format)
        if (!is.null(results[[n]]$seqAnno)) {
            # Standard format (genomics or properly formatted proteomics)
            seqAnno <- results[[n]]$seqAnno
        } else if (!is.null(results[[n]]$seqAnnoList)) {
            # Proteomics format with seqAnnoList - extract the first/only contrast
            seqAnno <- results[[n]]$seqAnnoList[[1]]
        } else {
            next # Skip if no data
        }

        df1 <- data.frame(
            "Gene" = seqAnno$gene_name,
            "LFC" = seqAnno$log2Ratio,
            "pValue" = seqAnno$pValue,
            "FDR" = seqAnno$fdr
        )
        colnames(df1)[2] <- names(results)[n]
        colnames(df1)[3] <- paste(names(results)[n], "pValue", sep = "_")
        colnames(df1)[4] <- paste(names(results)[n], "fdr", sep = "_")
        dfList[[names(results)[n]]] <- df1
    }
    # Make a big df of all the results together
    if (length(dfList) == 0) {
        return(NULL)
    } else if (length(dfList) == 1) {
        df2 <- dfList[[1]]
    } else {
        df2 <- merge(dfList[[1]], dfList[[2]], all = TRUE)
        if (length(dfList) >= 3) {
            for (i in 3:length(dfList)) {
                df2 <- merge(df2, dfList[[i]], all = TRUE)
            }
        }
    }
    for (i in grep("Gene|pValue|fdr", colnames(df2), invert = T)) {
        df2[, i][is.na(df2[, i])] <- 0
    }

    # Get a vector of all possible pathways for the selection drop down
    req(multiDeReactive$results)
    results <- multiDeReactive$results

    dfList <- list()
    for (n in seq_along(names(results))) {
        # Get the appropriate data table (degResultsTable for genomics, seqAnno/seqAnnoList for proteomics)
        degTable <- NULL
        if (!is.null(results[[n]]$degResultsTable)) {
            degTable <- results[[n]]$degResultsTable
        } else if (!is.null(results[[n]]$seqAnno)) {
            degTable <- results[[n]]$seqAnno
        } else if (!is.null(results[[n]]$seqAnnoList)) {
            degTable <- results[[n]]$seqAnnoList[[1]]
        }

        if (is.null(degTable)) {
            next # Skip if no data
        }

        # Check if required columns exist
        if (!"log2Ratio" %in% colnames(degTable)) {
            # Try to find alias
            if ("logFC" %in% colnames(degTable)) {
                degTable$log2Ratio <- degTable$logFC
            }
        }

        # Ensure gene_name exists
        if (is.null(degTable$gene_name)) {
            if ("Gene" %in% colnames(degTable)) {
                degTable$gene_name <- degTable$Gene
            } else {
                degTable$gene_name <- rownames(degTable)
            }
        }

        # Ensure vectors are not NULL to prevent data.frame from dropping columns
        if (is.null(degTable$gene_name)) degTable$gene_name <- character(0)
        if (is.null(degTable$log2Ratio)) degTable$log2Ratio <- numeric(0)
        if (is.null(degTable$pValue)) degTable$pValue <- numeric(0)
        if (is.null(degTable$fdr)) degTable$fdr <- numeric(0)

        df1 <- data.frame(
            "Gene" = degTable$gene_name,
            "LFC" = degTable$log2Ratio,
            "pValue" = degTable$pValue,
            "FDR" = degTable$fdr
        )
        colnames(df1)[2] <- names(results)[n]
        colnames(df1)[3] <- paste(names(results)[n], "pValue", sep = "_")
        colnames(df1)[4] <- paste(names(results)[n], "fdr", sep = "_")
        dfList[[names(results)[n]]] <- df1
    }
    # Make a big df of all the results together
    if (length(dfList) == 0) {
        return(NULL)
    } else if (length(dfList) == 1) {
        df2 <- dfList[[1]]
    } else {
        df2 <- merge(dfList[[1]], dfList[[2]], all = TRUE)
        if (length(dfList) >= 3) {
            for (i in 3:length(dfList)) {
                df2 <- merge(df2, dfList[[i]], all = TRUE)
            }
        }
    }
    for (i in grep("Gene|pValue|fdr", colnames(df2), invert = T)) {
        df2[, i][is.na(df2[, i])] <- 0
    }

    # Get a vector of all possible pathways for the selection drop down
    # PROTEOMICS FIX: Check if GO columns exist
    allPathways <- character(0)

    # Get the first available data table to extract GO terms
    firstTable <- NULL
    if (!is.null(results[[1]]$degResultsTable)) {
        firstTable <- results[[1]]$degResultsTable
    } else if (!is.null(results[[1]]$seqAnno)) {
        firstTable <- results[[1]]$seqAnno
    } else if (!is.null(results[[1]]$seqAnnoList)) {
        firstTable <- results[[1]]$seqAnnoList[[1]]
    }

    if (!is.null(firstTable) && !is.null(firstTable$`GO BP`)) {
        tryCatch(
            {
                allPathways <- c(
                    unique(unlist(strsplit(x = c(firstTable$`GO BP`), "; "))),
                    unique(unlist(strsplit(x = c(firstTable$`GO MF`), "; "))),
                    unique(unlist(strsplit(x = c(firstTable$`GO CC`), "; ")))
                )
                # Only run go2term if we have pathways
                if (length(allPathways) > 0) {
                    goTerms <- clusterProfiler::go2term(allPathways)
                    allPathways <- paste(goTerms[, 1], goTerms[, 2])
                }
            },
            error = function(e) {
                warning("Could not process GO terms: ", e$message)
            }
        )
    }

    # Split the combined LFC and p-vals into separate matrices and store in a list
    df3 <- df2

    # Ensure Gene column exists
    if (!"Gene" %in% colnames(df3)) {
        warning("Gene column not found in df2")
        return(NULL)
    }

    rownames(df3) <- NULL

    # Handle duplicates before column_to_rownames
    if (any(duplicated(df3$Gene))) {
        df3 <- df3[!duplicated(df3$Gene), ]
    }

    df3 <- column_to_rownames(df3, "Gene")
    df3[is.na(df3)] <- 1
    matLFC <- as.matrix(df3[, names(results), drop = FALSE])
    matRaw <- as.matrix(df3[, grep("pValue", colnames(df3)), drop = FALSE])
    matFDR <- as.matrix(df3[, grep("fdr", colnames(df3)), drop = FALSE])
    matList <- list("LFC" = matLFC, "pValue" = matRaw, "fdr" = matFDR)

    return(list(
        "allPathways" = allPathways,
        "matList" = matList,
        "dfList" = dfList,
        "df2" = df2,
        "df3" = df3
    ))
})

# 3. Tile Plot Logic

observe({
    req(multiMatReactive()$matList)
    updateSelectizeInput(session = session, inputId = "tileGenes", choices = multiMatReactive()$df2$Gene, selected = "", server = T)
    updateSelectizeInput(session = session, inputId = "tileGenesPathway", choices = multiMatReactive()$allPathways, selected = "", server = T)

    # Hide GO Pathway tab if no pathways available
    if (length(multiMatReactive()$allPathways) == 0) {
        hideTab(inputId = "multiTileTabs", target = "GO Pathway Genes")
    } else {
        showTab(inputId = "multiTileTabs", target = "GO Pathway Genes")
    }

    # Use updateSliderInput for sliders!
    # And ensure value is not below min
    valH <- (input$tileDEGNumber * 60)
    valW <- round(ncol(multiMatReactive()$matList$LFC) * 1.3) * 100

    updateSliderInput(session = session, inputId = "figHeightTile", value = max(100, valH))
    # updateSliderInput(session = session, inputId = "figWidthTile", value = max(100, valW))
})

observeEvent(
    {
        input$tileSortBy
        input$tileDEGNumber
        input$lfcTile
        input$pTypeTile
        input$pThresholdTile

        input$clusterRowsTile
        input$clusterColsTile

        input$colnamesTile

        input$tileGenes
        input$tileGenesText

        input$tileColourBlue
        input$tileColourWhite
        input$tileColourRed
        input$tileAtLow
        input$tileAtMid
        input$tileAtHigh
        input$tileLimitCPalette
        input$tileLimitCPaletteRev

        input$cellBorderColourTile1
        input$cellBorderSizeTile1
        input$cellBorderColourTile2
        input$cellBorderSizeTile2

        input$textSizeTile
        input$figWidthHeatmapTile
        input$figHeightTile

        lapply(seq_along(multiDeReactive$results), function(i) {
            input[[(paste0("multiGroupColour", names(multiDeReactive$results)[[i]]) %>% gsub(" ", "_", .))]]
        })
    },
    ignoreNULL = TRUE,
    ignoreInit = FALSE,
    {
        req(multiDeReactive$results)
        tryCatch({
        results <- multiDeReactive$results
        matList <- multiMatReactive()$matList
        dfList <- multiMatReactive()$dfList
        df2 <- multiMatReactive()$df2

        # get the latest colours chosen
        # Default colors from the palette (same as used in multiColourPickerUI)
        catalystCols <- c(
            "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
            "#FF7B00", "#FDC362", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
            "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
            "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
            "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
        )
        catalystCols <- c(catalystCols[c(seq(1, 30, by = 2), seq(2, 30, by = 2))])

        contrastColours <- setNames(lapply(seq_along(names(multiDeReactive$results)), function(i) {
            k <- names(multiDeReactive$results)[i]
            k_clean <- gsub(" ", "_", k)
            inputId <- paste0("multiGroupColour", k_clean)
            col <- input[[inputId]]
            if (is.null(col)) {
                # Use default palette color
                defaultCol <- rep(catalystCols, times = 5)[i]
                return(paste0(defaultCol, "FF"))
            }
            paste0(gplots::col2hex(col), "FF")
        }), names(multiDeReactive$results))
        coloursToUse <- unlist(contrastColours)

        # Determine whether the user wants to display raw or adjusted p-values
        if (input$pTypeTile == "Raw") {
            pTypeTile <- "pValue"
        } else {
            pTypeTile <- "fdr"
        }

        # Get the LFCs and the chosen p-values
        mat1 <- matList$LFC
        mat2 <- matList[[pTypeTile]]

        # Get the gene names from each contrast per the settings chosen
        highlightList <- list()
        for (i in seq_along(dfList)) {
            # Get the table for a contrast
            h <- dfList[[i]]
            # Filter the table based on the chosen settings
            h <- h[which(abs(h[, 2]) >= as.numeric(input$lfcTile) & h[, grep(pTypeTile, colnames(h))] <= as.numeric(input$pThresholdTile)), ]

            # If the user wants the top n by p-value, sort
            if (input$tileSortBy == "P-value") {
                h <- h[order(h[, grep(pTypeTile, colnames(h))]), ]
                genes <- h$Gene[1:min(nrow(h), input$tileDEGNumber)]
                # If the user wants the top n by LFC, sort
                # Two settings for LFC: By absolute value, or by top n + and -
            } else if (input$tileSortBy == "Absolute Log2FC") {
                h <- h[order(abs(h[, 2]), decreasing = T), ]
                genes <- h$Gene[1:min(nrow(h), input$tileDEGNumber)]
            } else if (input$tileSortBy == "Top n +/- Log2FC") {
                h <- h[order(h[, 2], decreasing = T), ]
                genes <- h$Gene[1:min(nrow(h), floor(input$tileDEGNumber / 2))]
                h <- h[order(h[, 2], decreasing = F), ]
                genes <- c(genes, h$Gene[1:min(nrow(h), floor(input$tileDEGNumber / 2))])
            }

            # Assign colours if genes found
            if (length(genes) > 0) {
                names(genes) <- rep(coloursToUse[names(dfList)[i]], length(genes))
                highlightList[[names(dfList)[i]]] <- genes
            }
        }

        # Take all the genes selected
        # req(length(highlightList) > 1) # This req might be too strict if only one contrast has DEGs

        if (length(highlightList) == 0) {
            tilePlots$top <- NULL # Clear previous plot
            tilePlots$mat1 <- NULL # Clear matrix
            tilePlots$df3 <- NULL # Clear table
            return()
        }

        highlighters <- unlist(highlightList, use.names = T)
        names(highlighters) <- gsub(".*.\\#", "#", names(highlighters))
        dupedHighlighters <- as.character(highlighters)[duplicated(as.character(highlighters))]

        # We want to know which gene came from which contrast, so we'll highlight the row labels accordingly
        highlightCols <- highlighters
        names(highlightCols)[which(as.character(highlightCols) %in% dupedHighlighters)] <- "#000000"
        highlightCols <- highlightCols[which(!is.na(highlightCols))]
        highlightCols <- highlightCols[!duplicated(highlightCols)]
        fontsizes <- rep(input$textSizeTile, length(highlightCols))
        fontfaces <- rep("bold", length(highlightCols))

        # Filter the lfc and p-values
        mat1 <- mat1[which(rownames(mat1) %in% highlightCols), , drop = FALSE]
        mat2 <- mat2[which(rownames(mat2) %in% highlightCols), , drop = FALSE]
        df3 <- df2[df2$Gene %in% highlightCols, ]
        # Rename LFC columns but not the Gene column
        lfc_cols <- grep("pValue|fdr|Gene", colnames(df3), invert = TRUE)
        colnames(df3)[lfc_cols] <- paste0(colnames(df3)[lfc_cols], "_log2Ratio")

        if (nrow(mat1) == 0) {
            tilePlots$top <- NULL
            tilePlots$mat1 <- NULL
            tilePlots$df3 <- NULL
            return()
        }

        # Reorder highlightCols to match mat1 rows
        highlightCols <- highlightCols[match(rownames(mat1), as.character(highlightCols))]
        fontcolors <- names(highlightCols)

        # Ensure valid colors
        fontcolors[is.na(fontcolors)] <- "#000000"

        # Explicit character conversion for anno_text (ComplexHeatmap requires character input)
        txt_labels <- as.character(highlightCols)
        rowAnno <- ComplexHeatmap::rowAnnotation(rows = anno_text(txt_labels, gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

        # Get some colours for the tiles
        at <- c(input$tileAtLow, input$tileAtMid, input$tileAtHigh)
        if (input$tileLimitCPalette == "None") {
            tileColours <- colorRamp2(
                breaks = rev(at),
                colors = c(input$tileColourRed, input$tileColourWhite, input$tileColourBlue)
            )
        } else {
            palette_colours <- brewer.pal(9, input$tileLimitCPalette)
            selected_colours <- palette_colours[c(1, 5, 9)]
            if (input$tileLimitCPaletteRev) {
                tileColours <- colorRamp2(breaks = at, colors = selected_colours)
            } else {
                tileColours <- colorRamp2(at, rev(selected_colours))
            }
        }

        # Make a nice top annotation thing
        ha <- HeatmapAnnotation(
            "Tests" = names(results),
            col = list(
                "Tests" = coloursToUse
            ),
            show_legend = T,
            annotation_legend_param = list(
                title_gp = gpar(fontsize = input$textSizeTile, fontface = "bold"),
                labels_gp = gpar(fontsize = input$textSizeTile)
            ),
            gp = gpar(col = input$cellBorderColourTile1, lwd = input$cellBorderSizeTile1 / 10)
        )

        colLabels <- colnames(mat1)

        hm <- Heatmap(
            matrix = mat1,
            cell_fun = function(j, i, x, y, width, height, fill) {
                gb <- textGrob("*")
                gb_w <- convertWidth(grobWidth(gb), "mm")
                gb_h <- convertHeight(grobHeight(gb), "mm")
                if (mat2[i, j] <= 0.0001) {
                    grid.text("****", x, y - gb_h * 0.5 + gb_w * 0.4, gp = gpar(fontsize = round(input$textSizeTile * 1.7), fontfamily = "HersheySerif"))
                } else if (mat2[i, j] <= 0.001) {
                    grid.text("***", x, y - gb_h * 0.5 + gb_w * 0.4, gp = gpar(fontsize = round(input$textSizeTile * 1.7), fontfamily = "HersheySerif"))
                } else if (mat2[i, j] <= 0.01) {
                    grid.text("**", x, y - gb_h * 0.5 + gb_w * 0.4, gp = gpar(fontsize = round(input$textSizeTile * 1.7), fontfamily = "HersheySerif"))
                } else if (mat2[i, j] <= 0.05) {
                    grid.text("*", x, y - gb_h * 0.5 + gb_w * 0.4, gp = gpar(fontsize = round(input$textSizeTile * 1.7), fontfamily = "HersheySerif"))
                }
            },
            name = "Log2FC",
            col = tileColours,
            right_annotation = rowAnno,
            cluster_columns = input$clusterColsTile,
            cluster_rows = input$clusterRowsTile,
            show_column_names = input$colnamesTile,
            show_row_names = F,
            column_names_centered = TRUE,
            column_split = factor(
                rep(c(names(results)), 1),
                levels = c(names(results))
            ),
            border = TRUE,
            column_gap = unit(1, "mm"),
            column_title = paste("Top", input$tileDEGNumber, "DEGs"),
            top_annotation = ha,
            column_labels = colLabels,
            rect_gp = gpar(col = input$cellBorderColourTile2, lwd = input$cellBorderSizeTile2 / 10),
            row_dend_width = unit(1, "cm"),
            column_dend_height = unit(1, "cm"),
            row_names_gp = gpar(fontsize = input$textSizeTile),
            column_names_gp = gpar(fontsize = input$textSizeTile),
            column_title_gp = gpar(fontsize = input$textSizeTile, fontface = "bold"),
            heatmap_legend_param = list(labels_gp = gpar(fontsize = input$textSizeTile), title_gp = gpar(fontsize = input$textSizeTile, fontface = "bold"))
        )

        # Store the data in reactive values
        tilePlots$top <- hm
        tilePlots$mat1 <- mat1
        tilePlots$df3 <- df3
        }, error = function(e) {
            if (inherits(e, "shiny.silent.error")) stop(e)
            msg <- paste0("Tile plot build error: ", e$message,
                          "\nCall: ", deparse(conditionCall(e)))
            showNotification(msg, type = "error", duration = 30)
            tilePlots$top <- NULL
            tilePlots$mat1 <- NULL
            tilePlots$df3 <- NULL
        })
    },
)

# Output heatmaps
output$TopNTileMap <- renderPlot(
    {
        if (is.null(tilePlots$mat1) || is.null(tilePlots$top)) {
            par(mar = c(0, 0, 0, 0))
            plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(1, 1, "No genes meet the selected criteria.", cex = 1.5)
            return()
        }
        validate(need(nrow(tilePlots$mat1) > 0, "No genes meet the selected criteria (LFC/P-value thresholds)."))
        tryCatch(
            draw(tilePlots$top, merge_legend = TRUE),
            error = function(e) {
                if (inherits(e, "shiny.silent.error")) stop(e)
                msg <- paste0("Top tile render error: ", e$message,
                              "\nCall: ", deparse(conditionCall(e)))
                showNotification(msg, type = "error", duration = 30)
                plot.new(); text(0.5, 0.5, msg, cex = 0.8)
            }
        )
    },
    width = function() as.numeric(input$figWidthTile),
    height = function() as.numeric(input$figHeightTile)
)

# download button for the heatmaps
# download button for the heatmaps
output$dlPlot_tileTop <- downloadHandler(
    filename = function() {
        generate_filename(
            base_name = input$filename_tileTop,
            design_name = "Multi",
            extension = tolower(input$downloadFormatMultiTile)
        )
    },
    content = function(file) {
        req(tilePlots$top)
        export_plot_standard(
            plot_object = tilePlots$top,
            filename = file,
            format = tolower(input$downloadFormatMultiTile),
            width = as.numeric(input$figWidthTile) / 75,
            height = as.numeric(input$figHeightTile) / 75,
            dpi = as.numeric(input$dpiMultiTile)
        )
    }
)

output$dlTable_tileTopRes <- downloadHandler(
    filename = function() {
        generate_filename(
            base_name = paste0(input$filename_tileTop, "_Log2FCs"),
            design_name = "Multi",
            extension = "xlsx"
        )
    },
    content = function(file) {
        req(tilePlots$mat1)
        mat1s <- tilePlots$mat1 %>%
            data.frame(check.names = F) %>%
            rownames_to_column("Feature")
        export_excel_standard(
            df = mat1s,
            filename = file,
            sheet_name = "Log2FC Matrix",
            table_name = "Log2FCMatrix"
        )
    }
)

# Output the results table
output$TopNTileTopTable <- renderDataTable({
    req(tilePlots$df3)
    df3 <- as.data.frame(tilePlots$df3)
    DT::datatable(
        data = df3,
        filter = "top",
        class = "cell-border stripe",
        rownames = FALSE,
        options = list(scrollX = TRUE, autoWidth = TRUE)
    ) %>%
        DT::formatStyle(
            columns = colnames(df3)[grep("log2Ratio", colnames(df3))],
            color = styleInterval(cuts = 0, values = c("blue", "darkorange")),
            fontWeight = "bold"
        ) %>%
        DT::formatSignif(
            columns = colnames(df3)[grep("pValue|fdr|log2Ratio", colnames(df3))],
            digits = 4
        ) %>%
        DT::formatStyle(
            columns = colnames(df3)[grep("pValue|fdr", colnames(df3))],
            color = styleInterval(
                cuts = 0.05,
                values = c("green", "black")
            ),
            fontWeight = "bold"
        )
})

# Download the results table
# Download the results table
output$dlTable_tileTop <- downloadHandler(
    filename = function() {
        generate_filename(
            base_name = paste0(input$filename_tileTop, "_Results"),
            design_name = "Multi",
            extension = "xlsx"
        )
    },
    content = function(file) {
        req(tilePlots$df3)
        export_excel_standard(
            df = tilePlots$df3,
            filename = file,
            sheet_name = "Top DEGs Results",
            table_name = "TopDEGsResults"
        )
    }
)

observeEvent(
    {
        input$tileGenes
        input$tileGenesText
        input$tileGenesPathway
        input$GOTooMany
        input$GOTooManyN
    },
    {
        customList <- list(
            "Custom" = list(),
            "GO" = list()
        )

        hGT <- NULL
        if (!is.null(input$tileGenesText)) {
            hGT <- strsplit(input$tileGenesText, "\\ +|\\/+|\\,+|\n+|\r+|\n\r+|\t+")[[1]]
        }
        customList[["Custom"]][["highlightersC"]] <- c(input$tileGenes, hGT)
        customList[["Custom"]][["plotTitle"]] <- paste0("Custom Gene Tile Plot")

        if (nchar(input$tileGenesPathway) >= 1) {
            # PROTEOMICS FIX: Check if GO columns exist before grep
            # Get the first available data table
            firstTable <- NULL
            if (!is.null(multiDeReactive$results[[1]]$degResultsTable)) {
                firstTable <- multiDeReactive$results[[1]]$degResultsTable
            } else if (!is.null(multiDeReactive$results[[1]]$seqAnno)) {
                firstTable <- multiDeReactive$results[[1]]$seqAnno
            } else if (!is.null(multiDeReactive$results[[1]]$seqAnnoList)) {
                firstTable <- multiDeReactive$results[[1]]$seqAnnoList[[1]]
            }

            if (!is.null(firstTable) && !is.null(firstTable$`GO BP`)) {
                customList[["GO"]][["highlightersC"]] <- firstTable$gene_name[grep(
                    (input$tileGenesPathway %>% strsplit(" ") %>% .[[1]] %>% .[1]),
                    paste(
                        firstTable$`GO BP`,
                        firstTable$`GO MF`,
                        firstTable$`GO CC`
                    )
                )]
                customList[["GO"]][["plotTitle"]] <- paste0(input$tileGenesPathway, " Tile Plot")
            }
        }

        lapply(c("Custom", "GO"), function(plotty) {
            if (length(customList[[plotty]][["highlightersC"]]) >= 1) { # Changed from 3 to 1 to allow smaller selections
                # Determine p-value type
                if (input$pTypeTile == "Raw") {
                    pTypeTile <- "pValue"
                } else {
                    pTypeTile <- "fdr"
                }

                df3 <- multiMatReactive()$df2[multiMatReactive()$df2$Gene %in% customList[[plotty]][["highlightersC"]], ]
                # Rename LFC columns but not the Gene column
                lfc_cols <- grep("pValue|fdr|Gene", colnames(df3), invert = TRUE)
                colnames(df3)[lfc_cols] <- paste0(colnames(df3)[lfc_cols], "_log2Ratio")

                # Some GO Pathways have hundreds or thousands of genes, and plotting them all might be ugly, so we can filter the table in two ways
                # This can either be done separately or additionally
                # 1. we can make sure any gene is DE in at least one contrast
                if (input$GOTooMany) {
                    df3 <- df3[apply(df3[grep(pTypeTile, colnames(df3))] <= as.numeric(input$pThresholdTile), 1, any), ]
                }

                # 2. we can select the top N genes from each contrast for that pathway
                if (!is.na(input$GOTooManyN)) {
                    tooManyHighlighters <- list()
                    for (n in names(multiDeReactive$results)) {
                        df4 <- df3[order(df3[, grep(paste0(n, "_", pTypeTile), colnames(df3))]), ]
                        tooManyHighlighters[[n]] <- df4$Gene[1:min(nrow(df4), input$GOTooManyN)]
                    }
                    tooManyHighlighters <- unique(unlist(tooManyHighlighters))
                    df3 <- df3[which(df3$Gene %in% tooManyHighlighters), ]
                }

                # If both options are selected, the table will be subset in the order of the code, so 1. and then 2.

                mat1 <- multiMatReactive()$matList$LFC
                mat1 <- mat1[which(rownames(mat1) %in% df3$Gene), , drop = FALSE]

                mat2 <- multiMatReactive()$matList[[pTypeTile]]
                mat2 <- mat2[which(rownames(mat2) %in% df3$Gene), , drop = FALSE]

                if (nrow(mat1) > 0) {
                    # Re-calculate rowAnno and ha based on current state for custom plots
                    # This is a simplified version, assuming coloursToUse, at, tileColours, pTypeTile are available from the main observeEvent
                    # For a robust solution, these might need to be passed or re-derived.
                    # For now, let's assume they are in scope or globally accessible.

                    # Get the latest colours chosen (re-deriving for safety)
                    contrastColours <- setNames(lapply((names(multiDeReactive$results) %>% gsub(" ", "_", .)), function(k) {
                        col <- input[[paste0("multiGroupColour", k)]]
                        if (is.null(col)) {
                            return("#000000FF")
                        }
                        paste0(gplots::col2hex(col), "FF")
                    }), names(multiDeReactive$results))
                    coloursToUse <- unlist(contrastColours)

                    # Get some colours for the tiles (re-deriving for safety)
                    at <- c(input$tileAtLow, input$tileAtMid, input$tileAtHigh)
                    if (input$tileLimitCPalette == "None") {
                        tileColours <- colorRamp2(
                            breaks = rev(at),
                            colors = c(input$tileColourRed, input$tileColourWhite, input$tileColourBlue)
                        )
                    } else {
                        palette_colours <- brewer.pal(9, input$tileLimitCPalette)
                        selected_colours <- palette_colours[c(1, 5, 9)]
                        if (input$tileLimitCPaletteRev) {
                            tileColours <- colorRamp2(breaks = at, colors = selected_colours)
                        } else {
                            tileColours <- colorRamp2(at, rev(selected_colours))
                        }
                    }

                    # Make a nice top annotation thing (re-deriving for safety)
                    ha <- HeatmapAnnotation(
                        "Tests" = names(multiDeReactive$results),
                        col = list(
                            "Tests" = coloursToUse
                        ),
                        show_legend = T,
                        annotation_legend_param = list(
                            title_gp = gpar(fontsize = input$textSizeTile, fontface = "bold"),
                            labels_gp = gpar(fontsize = input$textSizeTile)
                        ),
                        gp = gpar(col = input$cellBorderColourTile1, lwd = input$cellBorderSizeTile1 / 10)
                    )

                    # For custom plots, we don't have the same highlightCols logic for rowAnno,
                    # so we'll just show row names directly.
                    # If a rowAnno is desired, it would need to be constructed here based on customList[[plotty]][["highlightersC"]]
                    # For now, setting rowAnno to NULL or omitting it will show row names if show_row_names = TRUE

                    hm <- Heatmap(
                        matrix = mat1,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                            if (mat2[i, j] <= 0.0001) {
                                grid.text("****", x, y)
                            } else if (mat2[i, j] <= 0.001) {
                                grid.text("***", x, y)
                            } else if (mat2[i, j] <= 0.01) {
                                grid.text("**", x, y)
                            } else if (mat2[i, j] <= 0.05) {
                                grid.text("*", x, y)
                            }
                        },
                        name = "Log2FC",
                        col = tileColours,
                        # right_annotation = rowAnno[which(rownames(matList$LFC) %in% df3$Gene)], # This line was problematic, removed or adapted
                        cluster_columns = input$clusterColsTile,
                        cluster_rows = input$clusterRowsTile,
                        show_column_names = input$colnamesTile,
                        show_row_names = TRUE, # Changed to TRUE for custom plots
                        column_names_centered = TRUE,
                        column_split = factor(
                            rep(c(names(multiDeReactive$results)), 1),
                            levels = c(names(multiDeReactive$results))
                        ),
                        border = TRUE,
                        column_gap = unit(1, "mm"),
                        column_title = customList[[plotty]][["plotTitle"]],
                        top_annotation = ha,
                        # column_labels = colLabels, # colLabels might not be defined here, assuming it's from main observeEvent
                        rect_gp = gpar(col = input$cellBorderColourTile2, lwd = input$cellBorderSizeTile2 / 10),
                        row_dend_width = unit(1, "cm"),
                        column_dend_height = unit(1, "cm"),
                        row_names_gp = gpar(fontsize = input$textSizeTile),
                        column_names_gp = gpar(fontsize = input$textSizeTile),
                        column_title_gp = gpar(fontsize = input$textSizeTile, fontface = "bold"),
                        heatmap_legend_param = list(labels_gp = gpar(fontsize = input$textSizeTile), title_gp = gpar(fontsize = input$textSizeTile, fontface = "bold"))
                    )

                    # Store in reactive values
                    if (plotty == "Custom") {
                        tilePlots$custom <- hm
                        tilePlots$custom_df3 <- df3
                    } else {
                        tilePlots$go <- hm
                        tilePlots$go_df3 <- df3
                    }
                } else {
                    # If no genes remain after filtering, clear the plot and table
                    if (plotty == "Custom") {
                        tilePlots$custom <- NULL
                        tilePlots$custom_df3 <- NULL
                    } else {
                        tilePlots$go <- NULL
                        tilePlots$go_df3 <- NULL
                    }
                }
            } else {
                # If no highlighters, clear the plot and table
                if (plotty == "Custom") {
                    tilePlots$custom <- NULL
                    tilePlots$custom_df3 <- NULL
                } else {
                    tilePlots$go <- NULL
                    tilePlots$go_df3 <- NULL
                }
            }
        })
    }
)

# Output Custom/GO heatmaps and tables (Top Level)
lapply(c("Custom", "GO"), function(plotty) {
    output[[paste0(plotty, "TileMap")]] <- renderPlot(
        {
            hm <- if (plotty == "Custom") tilePlots$custom else tilePlots$go
            validate(need(!is.null(hm), "No genes selected or meeting criteria for this plot."))
            tryCatch(
                draw(hm, merge_legend = TRUE),
                error = function(e) {
                    if (inherits(e, "shiny.silent.error")) stop(e)
                    showNotification(paste(plotty, "tile plot render error:", e$message), type = "error", duration = 15)
                    plot.new(); text(0.5, 0.5, paste("Render error:", e$message), cex = 1.2)
                }
            )
        },
        width = function() as.numeric(input$figWidthTile),
        height = function() as.numeric(input$figHeightTile)
    )

    output[[paste0("dlPlot_tile", plotty)]] <- downloadHandler(
        filename = function() {
            generate_filename(
                base_name = input[[paste0("filename_tile", plotty)]],
                design_name = "Multi",
                extension = tolower(input$downloadFormatMultiTile)
            )
        },
        content = function(file) {
            if (plotty == "Custom") req(tilePlots$custom) else req(tilePlots$go)
            hm <- if (plotty == "Custom") tilePlots$custom else tilePlots$go
            export_plot_standard(
                plot_object = hm,
                filename = file,
                format = tolower(input$downloadFormatMultiTile),
                width = as.numeric(input$figWidthTile) / 75,
                height = as.numeric(input$figHeightTile) / 75,
                dpi = as.numeric(input$dpiMultiTile)
            )
        }
    )

    output[[paste0("TopNTile", plotty, "Table")]] <- renderDataTable({
        if (plotty == "Custom") req(tilePlots$custom_df3) else req(tilePlots$go_df3)
        df3 <- if (plotty == "Custom") tilePlots$custom_df3 else tilePlots$go_df3

        DT::datatable(
            data = df3,
            filter = "top",
            class = "cell-border stripe",
            rownames = FALSE,
            options = list(scrollX = TRUE, autoWidth = TRUE)
        ) %>%
            DT::formatStyle(
                columns = colnames(df3)[grep("log2Ratio", colnames(df3))],
                color = styleInterval(cuts = 0, values = c("blue", "darkorange")),
                fontWeight = "bold"
            ) %>%
            DT::formatSignif(
                columns = colnames(df3)[grep("pValue|fdr|log2Ratio", colnames(df3))],
                digits = 4
            ) %>%
            DT::formatStyle(
                columns = colnames(df3)[grep("pValue|fdr", colnames(df3))],
                color = styleInterval(
                    cuts = 0.05,
                    values = c("green", "black")
                ),
                fontWeight = "bold"
            )
    })

    output[[paste0("dlTable_tile", plotty)]] <- downloadHandler(
        filename = function() {
                generate_filename(
                base_name = paste0(input[[paste0("filename_tile", plotty)]], "_Results"),
                design_name = "Multi",
                extension = "xlsx"
            )
        },
        content = function(file) {
            if (plotty == "Custom") req(tilePlots$custom_df3) else req(tilePlots$go_df3)
            df3 <- if (plotty == "Custom") tilePlots$custom_df3 else tilePlots$go_df3
            export_excel_standard(
                df = df3,
                filename = file,
                sheet_name = paste0(plotty, " Results"),
                table_name = paste0(plotty, "Results")
            )
        }
    )
})
