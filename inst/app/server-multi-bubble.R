# Server logic for Multi Bubble Plot

dataToPlotBubble <- list(
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

updateSelectInput(
    session = session,
    inputId = "subdataToBubble",
    choices = dataToPlotBubble[[1]],
    selected = dataToPlotBubble[[1]][1]
)

observeEvent(input$dataToBubble, {
    updateSelectInput(
        session = session,
        inputId = "subdataToBubble",
        choices = dataToPlotBubble[[input$dataToBubble]],
        selected = dataToPlotBubble[[input$dataToBubble]][1]
    )
})

observeEvent(
    {
        input$dataToBubble
        input$subdataToBubble
        input$bubblePathwayNumber
        input$colourByBubble
        input$colnamesBubble
        input$textSizeBubble
        input$figWidthBubble
        input$figHeightBubble
        lapply(seq_along(multiDeReactive$results), function(i) {
            input[[(paste0("multiGroupColour", names(multiDeReactive$results)[[i]]) %>% gsub(" ", "_", .))]]
        })
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE,
    {
        req(input$dataToBubble)
        dataToUseBubble <- switch(input$dataToBubble,
            "ORA" = "oraList",
            "GSEA" = "gseaList"
        )
        pathwayToUseBubble <- switch(input$dataToBubble,
            "ORA" = "ora",
            "GSEA" = "gsea"
        )
        columnsToUse <- switch(input$dataToBubble,
            "ORA" = c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "geneName", "Label", "Contrast"),
            "GSEA" = c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment", "geneName", "Label", "Contrast")
        )

        req(any(length(upsetGeneList[[dataToUseBubble]][[input$subdataToBubble]]) >= 2))
        req(!is.null(multiDeReactive$results))
        results <- multiDeReactive$results
        geneListBubble <- upsetGeneList[[dataToUseBubble]][[input$subdataToBubble]]

        # get the latest colours chosen
        contrastColoursBubble <- list()
        for (i in seq_along(results)) {
            contrastColoursBubble[[names(results)[i]]] <- paste0(col2hex(input[[paste0("multiGroupColour", names(results)[i] %>% gsub(" ", "_", .))]]))
        }
        coloursToUseBubble <- unlist(contrastColoursBubble)

        contrastColoursBubble <- setNames(lapply((names(multiDeReactive$results) %>% gsub(" ", "_", .)), function(k) {
            paste0(col2hex(input[[paste0("multiGroupColour", k)]]), "FF")
        }), names(multiDeReactive$results))
        coloursToUseBubble <- unlist(contrastColoursBubble)
        # req(grepl("^#", coloursToUseBubble[[1]]))

        # get a list of df
        dfListBubble <- setNames(lapply(names(geneListBubble), function(con) {
            if (input$dataToBubble == "ORA") {
                # Check if result exists
                if (is.null(results[[con]][[pathwayToUseBubble]][[substr(input$subdataToBubble, 1, 2)]][[substr(input$subdataToBubble, 3, nchar(input$subdataToBubble))]])) {
                    return(NULL)
                }

                r1 <- as.data.frame(results[[con]][[pathwayToUseBubble]][[substr(input$subdataToBubble, 1, 2)]][[substr(input$subdataToBubble, 3, nchar(input$subdataToBubble))]]@result)
            } else if (input$dataToBubble == "GSEA") {
                if (is.null(results[[con]][[pathwayToUseBubble]][[substr(input$subdataToBubble, 1, 2)]])) {
                    return(NULL)
                }

                r1 <- as.data.frame(results[[con]][[pathwayToUseBubble]][[substr(input$subdataToBubble, 1, 2)]]@result)
            }
            r2 <- r1[which(r1$Label %in% geneListBubble[[con]]), ]
            return(r2)
        }), names(geneListBubble))

        # Remove NULLs
        dfListBubble <- dfListBubble[!sapply(dfListBubble, is.null)]

        # get a list of pathways to extract/plot/highlight
        highlightListBubble <- setNames(lapply(names(dfListBubble), function(con) {
            x <- dfListBubble[[con]]
            if (nrow(x) > input$bubblePathwayNumber) {
                x$Label[1:input$bubblePathwayNumber]
            } else {
                x$Label
            }
        }), names(dfListBubble))
        highlightersBubble <- unlist(highlightListBubble)
        names(highlightersBubble) <- rep(names(highlightListBubble), lengths(highlightListBubble))
        # names(highlightersBubble) <- gsub("^\\d+|\\d+$", "", names(highlightersBubble))

        # we want to know which pathway came from which contrast, so we'll highlight the row labels accordingly
        highlightColsBubble <- highlightersBubble
        names(highlightColsBubble) <- as.character(coloursToUseBubble)[match(names(highlightColsBubble), names(coloursToUseBubble))]
        names(highlightColsBubble)[duplicated(highlightColsBubble) | duplicated(highlightColsBubble, fromLast = TRUE)] <- "#000000"
        highlightColsBubble <- highlightColsBubble[!duplicated(highlightColsBubble)]

        # get a list of df of only the highlighted pathways
        plotDFBubbleListBubble <- setNames(lapply(names(dfListBubble), function(con) {
            x <- dfListBubble[[con]]
            x$Contrast <- con
            x <- x[which(x$Label %in% highlightersBubble), columnsToUse]
        }), names(dfListBubble))

        # make that list a single df for plotting
        plotDFBubble <- rbindlist(plotDFBubbleListBubble)
        plotDFBubble$Contrast <- factor(plotDFBubble$Contrast, levels = names(multiDeReactive$results))
        plotDFBubble$Colours <- names(highlightColsBubble)[match(plotDFBubble$Label, as.character(highlightColsBubble))]

        # make the plot
        if (input$dataToBubble == "ORA") {
            bubblePlot <- ggplot(plotDFBubble, aes(x = Contrast, y = Label, fill = Contrast, size = Count))
            bubblePlot <- bubblePlot + scale_fill_manual(values = coloursToUseBubble)
        } else if (input$dataToBubble == "GSEA") {
            # cols <- colorRamp2(breaks = c(-3, 0, 3), colors = c("#113D69FF", "#FFFFFFFF", "#801717FF"))
            # if(sum(plotDFBubble$NES < 0) == nrow(plotDFBubble) | sum(plotDFBubble$NES > 0) == nrow(plotDFBubble)) {
            #   col_fun <- cols(c(min(plotDFBubble$NES), mean(plotDFBubble$NES), max(plotDFBubble$NES)))
            # } else {
            #   col_fun <- cols(c(min(plotDFBubble$NES), 0, max(plotDFBubble$NES)))
            # }
            bubblePlot <- ggplot(plotDFBubble, aes(x = Contrast, y = Label, fill = NES, size = setSize))
            bubblePlot <- bubblePlot + scale_fill_gradient2(low = "royalblue4", mid = "white", high = "firebrick4", midpoint = 0)
        }
        bubblePlot <- bubblePlot + geom_point(alpha = 0.8, pch = 21)
        bubblePlot <- bubblePlot + scale_size(range = c(input$textSizeBubble / 2, input$textSizeBubble * 1.5), name = "Gene Count")
        bubblePlot <- bubblePlot + theme_classic()
        highlightColsBubble <- highlightColsBubble[match(sort(unique(plotDFBubble$Label)), as.character(highlightColsBubble))]
        bubblePlot <- bubblePlot + theme(
            axis.text.x = element_text(size = input$textSizeBubble, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = input$textSizeBubble, colour = names(highlightColsBubble), face = "bold"),
            axis.title.x = element_blank(),
            legend.title = element_text(size = input$textSizeBubble, face = "bold"),
            legend.text = element_text(size = input$textSizeBubble, face = "bold")
        )
        if (input$dataToBubble == "ORA") {
            bubblePlot <- bubblePlot + guides(fill = guide_legend(override.aes = list(size = input$textSizeBubble)))
        }
        if (!input$colnamesBubble) {
            bubblePlot <- bubblePlot + theme(axis.text.x = element_blank())
        }

        # Output the results table
        columnsToTableBubble <- switch(input$dataToBubble,
            "ORA" = c("ID", "Description", "Contrast", "Count", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneName"),
            "GSEA" = c("ID", "Description", "Contrast", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "geneName")
        )
        output$BubbleTable <- renderDataTable({
            renderTable <- 'DT::datatable(
        data = (plotDFBubble %>% dplyr::select(all_of(columnsToTableBubble))),
        filter = "top",
        class = "cell-border stripe",
        rownames = FALSE,
        options = list(scrollX = TRUE, autoWidth = TRUE)
        ) %>%
        DT::formatSignif(columns = c("pvalue", "p.adjust"), digits = 4) %>%
        DT::formatStyle(columns = c("pvalue", "p.adjust"), color = styleInterval(cuts = 0.05, values = c("green", "black")), fontWeight = "bold")'
            if (input$dataToBubble == "ORA") {
                eval(parse(text = renderTable))
            } else {
                eval(parse(text = paste0(
                    renderTable,
                    '%>% DT::formatStyle(columns = c("enrichmentScore", "NES"), color = styleInterval(cuts = 0, values = c("blue", "orange")), fontWeight = "bold") %>%
          DT::formatSignif(columns = c("enrichmentScore", "NES"), digits = 4)'
                )))
            }
        })
        # Download the results table
        output$dlBubbleTable <- downloadHandler(
            filename = function() {
                paste0("Top_", input$bubblePathwayNumber, "_", input$dataToBubble, "_Bubble.xlsx")
            },
            content = function(file) {
                openxlsx::write.xlsx((plotDFBubble %>% dplyr::select(all_of(columnsToTableBubble))), file = file)
            }
        )

        output$BubblePlot <- renderPlot(
            {
                print(bubblePlot)
            },
            width = as.numeric(input$figWidthBubble),
            height = as.numeric(input$figHeightBubble)
        )

        output$dlBubblePlot <- downloadHandler(
            filename = function() {
                paste0("bubble.pdf")
            },
            content = function(file) {
                pdf(
                    file = file,
                    width = (as.numeric(input$figWidthBubble) / 70),
                    height = (as.numeric(input$figHeightBubble) / 80),
                )
                print(bubblePlot)
                dev.off()
            }
        )
    }
)
