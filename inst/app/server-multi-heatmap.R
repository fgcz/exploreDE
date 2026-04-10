# Server logic for Multi Pathway Heatmap

# Define data structure options ----
dataToPlotHeatmap <- list(
    "ORA" = c(
        "BP: Up-Regulated" = "BPupGenes",
        "BP: Down-Regulated" = "BPdownGenes",
        "BP: Both Directions" = "BPbothGenes",
        "MF: Up-Regulated" = "MFupGenes",
        "MF: Down-Regulated" = "MFdownGenes",
        "MF: Both Directions" = "MFbothGenes",
        "CC: Up-Regulated" = "CCupGenes",
        "CC: Down-Regulated" = "CCdownGenes",
        "CC: Both Directions" = "CCbothGenes"
    ),
    "GSEA" = c(
        "BP" = "BP",
        "MF" = "MF",
        "CC" = "CC"
    )
)

# Initialize subdataset selector ----
updateSelectInput(
    session = session,
    inputId = "subdataToHeatmap",
    choices = dataToPlotHeatmap[[1]],
    selected = dataToPlotHeatmap[[1]][1]
)

# Update subdataset choices when main data type changes ----
observeEvent(input$dataToHeatmap, {
    updateSelectInput(
        session = session,
        inputId = "subdataToHeatmap",
        choices = dataToPlotHeatmap[[input$dataToHeatmap]],
        selected = dataToPlotHeatmap[[input$dataToHeatmap]][1]
    )
})

# Show NES option only for GSEA ----
output$heatmapNESOption <- renderUI({
    if (input$dataToHeatmap == "GSEA") {
        checkboxInput(
            inputId = "heatmapUseNES",
            label = "Use NES instead of p-value/FDR?",
            value = TRUE
        )
    }
})

# Main reactive: prepare heatmap data ----
multiPathwayHeatmapDataReactive <- reactive({
    req(multiDeReactive$results)
    req(length(multiDeReactive$results) > 0)
    req(input$dataToHeatmap)
    req(input$subdataToHeatmap)

    results <- multiDeReactive$results

    dataToUse <- switch(input$dataToHeatmap,
        "ORA" = "oraList",
        "GSEA" = "gseaList"
    )

    pathwayToUse <- switch(input$dataToHeatmap,
        "ORA" = "ora",
        "GSEA" = "gsea"
    )

    # Check if data exists
    req(any(length(upsetGeneList[[dataToUse]][[input$subdataToHeatmap]]) >= 2))

    geneListHeatmap <- upsetGeneList[[dataToUse]][[input$subdataToHeatmap]]

    # Extract data frames for each contrast ----
    dfListHeatmap <- setNames(lapply(names(geneListHeatmap), function(con) {
        if (input$dataToHeatmap == "ORA") {
            # Check existence
            if (is.null(results[[con]][[pathwayToUse]][[substr(input$subdataToHeatmap, 1, 2)]][[substr(input$subdataToHeatmap, 3, nchar(input$subdataToHeatmap))]])) {
                return(NULL)
            }

            r1 <- as.data.frame(results[[con]][[pathwayToUse]][[substr(input$subdataToHeatmap, 1, 2)]][[substr(input$subdataToHeatmap, 3, nchar(input$subdataToHeatmap))]]@result)
        } else if (input$dataToHeatmap == "GSEA") {
            if (is.null(results[[con]][[pathwayToUse]][[substr(input$subdataToHeatmap, 1, 2)]])) {
                return(NULL)
            }

            r1 <- as.data.frame(results[[con]][[pathwayToUse]][[substr(input$subdataToHeatmap, 1, 2)]]@result)
        }
        r2 <- r1[which(r1$Label %in% geneListHeatmap[[con]]), ]
        return(r2)
    }), names(geneListHeatmap))

    # Remove NULLs
    dfListHeatmap <- dfListHeatmap[!sapply(dfListHeatmap, is.null)]

    # Get top pathways from each contrast ----
    highlightListHeatmap <- setNames(lapply(names(dfListHeatmap), function(con) {
        x <- dfListHeatmap[[con]]
        if (nrow(x) > input$heatmapPathwayNumber) {
            x$Label[1:input$heatmapPathwayNumber]
        } else {
            x$Label
        }
    }), names(dfListHeatmap))

    allHighlightedPathways <- unique(unlist(highlightListHeatmap))

    req(length(allHighlightedPathways) > 0)

    # Build matrix for heatmap ----
    # Determine what values to extract
    useNES <- !is.null(input$heatmapUseNES) && input$heatmapUseNES && input$dataToHeatmap == "GSEA"

    if (useNES) {
        valueCol <- "NES"
    } else {
        valueCol <- input$heatmapColourBy
    }

    # Create matrix: rows = pathways, cols = contrasts
    heatmapMat <- matrix(
        NA,
        nrow = length(allHighlightedPathways),
        ncol = length(dfListHeatmap),
        dimnames = list(allHighlightedPathways, names(dfListHeatmap))
    )

    # Fill matrix with values
    for (con in names(dfListHeatmap)) {
        df <- dfListHeatmap[[con]]
        for (pathway in allHighlightedPathways) {
            if (pathway %in% df$Label) {
                val <- df[df$Label == pathway, valueCol]
                if (length(val) > 0 && !is.na(val)) {
                    if (useNES) {
                        # Use NES directly
                        heatmapMat[pathway, con] <- val
                    } else {
                        # Use -log10 transformation for p-values
                        heatmapMat[pathway, con] <- -log10(val + 1e-300) # Add small value to avoid log(0)
                    }
                }
            }
        }
    }

    # Replace remaining NAs with 0
    heatmapMat[is.na(heatmapMat)] <- 0

    # Create full table for download
    fullTable <- do.call(rbind, lapply(names(dfListHeatmap), function(con) {
        df <- dfListHeatmap[[con]]
        df <- df[df$Label %in% allHighlightedPathways, ]
        if (nrow(df) > 0) {
            df$Contrast <- con
            df
        } else {
            NULL
        }
    }))

    return(list(
        "heatmapMat" = heatmapMat,
        "fullTable" = fullTable,
        "useNES" = useNES,
        "valueCol" = valueCol
    ))
})

multiPathwayHeatmapObjectReactive <- reactive({
    req(multiPathwayHeatmapDataReactive())
    req(multiDeReactive$results)

    plotData <- multiPathwayHeatmapDataReactive()
    heatmapMat <- plotData$heatmapMat
    useNES <- plotData$useNES
    results <- multiDeReactive$results

    if (nrow(heatmapMat) == 0 || ncol(heatmapMat) == 0) return(NULL)

    # Get contrast colours for top annotation
    contrastColours <- setNames(
        lapply(names(results), function(k) {
            col <- input[[paste0("multiGroupColour", gsub(" ", "_", k))]]
            if (!is.null(col)) {
                tryCatch(col2hex(col), error = function(e) "#000000")
            } else {
                RColorBrewer::brewer.pal(max(3, length(results)), "Dark2")[match(k, names(results))]
            }
        }),
        names(results)
    )

    contrastColours <- contrastColours[colnames(heatmapMat)]

    topAnnotation <- HeatmapAnnotation(
        Contrast = colnames(heatmapMat),
        col = list(Contrast = unlist(contrastColours)),
        annotation_name_gp = gpar(fontsize = input$textSizeHeatmap),
        annotation_legend_param = list(
            Contrast = list(
                title = "Contrast",
                title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
                labels_gp = gpar(fontsize = input$textSizeHeatmap)
            )
        ),
        show_legend = TRUE,
        border = TRUE,
        gp = gpar(col = "black", lwd = 1)
    )

    if (useNES) {
        maxVal <- max(abs(range(heatmapMat, na.rm = TRUE)))
        colFun <- colorRamp2(
            c(-maxVal, 0, maxVal),
            c("#1965B0", "#FFFFFF", "#DC050C")
        )
        legendTitle <- "NES"
    } else {
        colFun <- colorRamp2(
            c(min(heatmapMat, na.rm = TRUE), max(heatmapMat, na.rm = TRUE)),
            c(input$heatmapColourLow, input$heatmapColourHigh)
        )
        if (plotData$valueCol == "pvalue") {
            legendTitle <- "-log10(p-value)"
        } else {
            legendTitle <- "-log10(FDR)"
        }
    }

    ht <- Heatmap(
        heatmapMat,
        name = legendTitle,
        col = colFun,
        top_annotation = topAnnotation,
        cluster_rows = input$clusterRowsHeatmap,
        cluster_columns = input$clusterColsMultiHeatmap,
        show_row_names = input$rownamesHeatmap,
        show_column_names = input$colnamesHeatmap,
        row_names_gp = gpar(fontsize = input$textSizeHeatmap),
        column_names_gp = gpar(fontsize = input$textSizeHeatmap),
        column_names_rot = 45,
        row_names_side = "right",
        row_names_max_width = max_text_width(rownames(heatmapMat), gp = gpar(fontsize = input$textSizeHeatmap)),
        column_names_max_height = unit(8, "cm"),
        heatmap_legend_param = list(
            title_gp = gpar(fontsize = input$textSizeHeatmap, fontface = "bold"),
            labels_gp = gpar(fontsize = input$textSizeHeatmap),
            legend_height = unit(6, "cm")
        ),
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
            val <- heatmapMat[i, j]
            if (useNES) {
                if (!is.na(val) && abs(val) > 0.01) {
                    grid.text(sprintf("%.2f", val), x, y, gp = gpar(fontsize = input$textSizeHeatmap * 0.8))
                }
            } else {
                if (!is.na(val) && val > 0) {
                    grid.text(sprintf("%.2f", val), x, y, gp = gpar(fontsize = input$textSizeHeatmap * 0.8))
                }
            }
        }
    )
    ht
})

# Render heatmap ----
output$pathwayHeatmapPlot <- renderPlot(
    {
        req(multiPathwayHeatmapObjectReactive())
        ht <- multiPathwayHeatmapObjectReactive()
        draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", padding = unit(c(2, 2, 2, 2), "cm"))
    },
    width = function() input$figWidthHeatmap,
    height = function() input$figHeightHeatmap
)

# Render data table ----
output$pathwayHeatmapTable <- renderDataTable({
    req(multiPathwayHeatmapDataReactive())

    fullTable <- multiPathwayHeatmapDataReactive()$fullTable

    # Select relevant columns
    if (input$dataToHeatmap == "ORA") {
        cols <- c("Contrast", "Description", "Count", "GeneRatio", "pvalue", "p.adjust")
    } else {
        cols <- c("Contrast", "Description", "setSize", "NES", "pvalue", "p.adjust")
    }

    tableData <- fullTable[, intersect(cols, colnames(fullTable))]

    dt <- DT::datatable(
        tableData,
        filter = "top",
        class = "cell-border stripe",
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 25)
    ) %>%
        DT::formatRound(columns = c("pvalue", "p.adjust"), digits = 3)

    # Add colour formatting for p-values (green if < 0.05, black otherwise)
    dt <- dt %>%
        DT::formatStyle(
            "pvalue",
            color = DT::styleInterval(0.05, c("green", "black"))
        ) %>%
        DT::formatStyle(
            "p.adjust",
            color = DT::styleInterval(0.05, c("green", "black"))
        )

    # Add colour formatting for NES if GSEA (red for positive, blue for negative)
    if (input$dataToHeatmap == "GSEA" && "NES" %in% colnames(tableData)) {
        dt <- dt %>%
            DT::formatRound(columns = "NES", digits = 3) %>%
            DT::formatStyle(
                "NES",
                color = DT::styleInterval(0, c("blue", "red"))
            )
    }

    dt
})

# Download PDF ----
# Download PDF ----
output$dlPlot_multiHeatmap <- downloadHandler(
    filename = function() {
        generate_filename(
            base_name = paste0(input$filename_multiHeatmap, "_Heatmap"),
            design_name = "Multi",
            extension = tolower(input$downloadFormatMultiHeatmap)
        )
    },
    content = function(file) {
        req(multiPathwayHeatmapObjectReactive())
        ht <- multiPathwayHeatmapObjectReactive()
        export_plot_standard(
            plot_object = ht,
            filename = file,
            format = tolower(input$downloadFormatMultiHeatmap),
            width = as.numeric(input$figWidthHeatmap) / 72,
            height = as.numeric(input$figHeightHeatmap) / 72,
            dpi = as.numeric(input$dpiMultiHeatmap)
        )
    }
)

# Download Excel table ----
output$dlTable_multiHeatmap <- downloadHandler(
    filename = function() {
        generate_filename(
            base_name = paste0(input$filename_multiHeatmap, "_Data"),
            design_name = "Multi",
            extension = "xlsx"
        )
    },
    content = function(file) {
        req(multiPathwayHeatmapDataReactive())
        fullTable <- multiPathwayHeatmapDataReactive()$fullTable
        export_excel_standard(
            df = fullTable,
            filename = file,
            sheet_name = "Pathway Data",
            table_name = "PathwayData"
        )
    }
)
