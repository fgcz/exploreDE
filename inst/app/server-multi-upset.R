# Server logic for Multi Upset Plot

# Update dataToUpset choices based on data type (hide ORA/GSEA for proteomics)
observe({
    req(multiDeReactive$results)
    req(length(multiDeReactive$results) > 0)

    # Check if any result is proteomics
    dataTypes <- sapply(multiDeReactive$results, function(x) x$dataType)
    isProteomics <- any(dataTypes == "proteomics")

    if (isProteomics) {
        # Only show "Genes" for proteomics
        updateSelectInput(session, "dataToUpset",
                         choices = c("Genes"),
                         selected = "Genes")
    } else {
        # Show all options for genomics
        updateSelectInput(session, "dataToUpset",
                         choices = c("Genes", "ORA", "GSEA"),
                         selected = "Genes")
    }
})

# 1. Generate Gene Lists (Ported from multiDEG/server-import2.R)
multiUpsetReactive <- reactiveValues(plot = NULL)
upsetGeneList <- reactiveValues(results = NULL)
upsetGeneListToCopy <- reactiveVal(NULL)

observeEvent(input$copy_link, {
    req(upsetGeneListToCopy())
    session$sendCustomMessage("txt", upsetGeneListToCopy())
})

observeEvent(
    {
        input$lfcUpset
        input$pTypeUpset
        input$pThresholdUpset
        multiDeReactive$results
    },
    {
        req(multiDeReactive$results)
        req(length(multiDeReactive$results) > 0)

        req(input$pTypeUpset)
        if (input$pTypeUpset == "Raw") {
            pTypeUpset <- "pValue"
            pTypeUpset2 <- "pvalue"
        } else {
            pTypeUpset <- "fdr"
            pTypeUpset2 <- "p.adjust"
        }

        results <- multiDeReactive$results

        geneLists <- list("upGenes" = list(), "downGenes" = list(), "bothGenes" = list())

        for (i in seq_along(results)) {
            # In multiDEG, results[[i]] is an object with seqAnno, etc.
            # Our multiDeReactive$results has the same structure.

            seqAnno1 <- results[[i]]$seqAnno

            # Filter usedInTest if available (it is in exploreDE output)
            if ("usedInTest" %in% colnames(seqAnno1)) {
                seqAnno1 <- seqAnno1[seqAnno1$usedInTest == TRUE, ]
            }

            # Find p-value column
            # Note: exploreDE's convert_genomics_se standardizes columns to pValue/fdr/log2Ratio
            # But let's be robust.
            p <- colnames(seqAnno1)[grep(paste(c(pTypeUpset, pTypeUpset2), collapse = "|"), colnames(seqAnno1), ignore.case = TRUE)][1]

            # If not found, try standard names
            if (is.na(p)) {
                if (pTypeUpset == "pValue") {
                    p <- "pValue"
                } else {
                    p <- "fdr"
                }
            }

            # Calculate lists
            up <- seqAnno1$gene_name[which(as.numeric(seqAnno1[[p]]) <= as.numeric(input$pThresholdUpset) & as.numeric(seqAnno1$log2Ratio) >= as.numeric(input$lfcUpset))]
            down <- seqAnno1$gene_name[which(as.numeric(seqAnno1[[p]]) <= as.numeric(input$pThresholdUpset) & as.numeric(seqAnno1$log2Ratio) <= -as.numeric(input$lfcUpset))]
            both <- seqAnno1$gene_name[which(as.numeric(seqAnno1[[p]]) <= as.numeric(input$pThresholdUpset) & abs(as.numeric(seqAnno1$log2Ratio)) >= as.numeric(input$lfcUpset))]

            geneLists[["upGenes"]][[names(results)[i]]] <- up
            geneLists[["downGenes"]][[names(results)[i]]] <- down
            geneLists[["bothGenes"]][[names(results)[i]]] <- both
        }

        # Pathway Lists (ORA/GSEA)
        # Only if runGO was true for the first dataset (or check per dataset)
        # multiDEG assumes all have it if the first has it.
        # We should check per dataset.

        # Make a big df of all the results together
        dfList <- list()
        for (n in seq_along(names(results))) {
            df1 <- data.frame(
                "Gene" = results[[n]]$seqAnno$gene_name,
                "LFC" = results[[n]]$seqAnno$log2Ratio,
                "pValue" = results[[n]]$seqAnno$pValue,
                "FDR" = results[[n]]$seqAnno$fdr
            )
            colnames(df1)[2] <- names(results)[n]
            colnames(df1)[3] <- paste(names(results)[n], "pValue", sep = "_")
            colnames(df1)[4] <- paste(names(results)[n], "fdr", sep = "_")
            dfList[[names(results)[n]]] <- df1
        }

        if (length(dfList) == 1) {
            df2 <- dfList[[1]]
        } else {
            df2 <- merge(dfList[[1]], dfList[[2]], all = TRUE)
            if (length(results) >= 3) {
                for (i in 3:length(results)) {
                    df2 <- merge(df2, dfList[[i]], all = TRUE)
                }
            }
        }
        oraList <- list()
        gseaList <- list()

        # Initialize lists
        for (ont in c("BP", "MF", "CC")) {
            gseaList[[ont]] <- list()
            for (sig in c("upGenes", "downGenes", "bothGenes")) {
                oraList[[paste0(ont, sig)]] <- list()
            }
        }

        # Populate
        for (d in names(results)) {
            if (!is.null(results[[d]]$param$runGO) && results[[d]]$param$runGO) {
                # GSEA
                for (ont in c("BP", "MF", "CC")) {
                    if (!is.null(results[[d]]$gsea[[ont]])) {
                        res <- results[[d]]$gsea[[ont]]@result
                        col <- if (pTypeUpset == "fdr") "p.adjust" else "pvalue"

                        if (col %in% colnames(res)) {
                            filtered <- res$Label[which(as.numeric(res[[col]]) <= as.numeric(input$pThresholdUpset))]
                            gseaList[[ont]][[d]] <- filtered
                        }
                    }
                }

                # ORA
                for (ont in c("BP", "MF", "CC")) {
                    for (sig in c("bothGenes", "upGenes", "downGenes")) {
                        if (!is.null(results[[d]]$ora[[ont]][[sig]])) {
                            res <- results[[d]]$ora[[ont]][[sig]]@result
                            col <- if (pTypeUpset == "fdr") "p.adjust" else "pvalue"

                            if (col %in% colnames(res)) {
                                filtered <- res$Label[which(as.numeric(res[[col]]) <= as.numeric(input$pThresholdUpset))]
                                oraList[[paste0(ont, sig)]][[d]] <- filtered
                            }
                        }
                    }
                }
            }
        }

        upsetGeneList$genesList <- geneLists
        upsetGeneList$oraList <- oraList
        upsetGeneList$gseaList <- gseaList
    }
)


# 2. Upset Plot Logic (Ported from multiDEG/server-upset4.R)
dataToPlot <- list(
    "Genes" = c(
        "Up-Regulated" = "upGenes",
        "Down-regulated" = "downGenes",
        "Both Directions" = "bothGenes"
    ),
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
    inputId = "subdataToUpset",
    choices = dataToPlot[[1]],
    selected = dataToPlot[[1]][1]
)

observeEvent(input$dataToUpset, {
    updateSelectInput(
        session = session,
        inputId = "subdataToUpset",
        choices = dataToPlot[[input$dataToUpset]],
        selected = dataToPlot[[input$dataToUpset]][1]
    )
})

observe({
    req(input$dataToUpset)
    dataToUse <- switch(input$dataToUpset,
        "Genes" = "genesList",
        "ORA" = "oraList",
        "GSEA" = "gseaList"
    )
    req(upsetGeneList[[dataToUse]])
    req(multiDeReactive$results)
    results <- multiDeReactive$results
    gl1 <- upsetGeneList[[dataToUse]][[input$subdataToUpset]]
    req(length(gl1) > 1)
    output$UpsetTableGenesUI <- renderUI({
        checkboxGroupInput(
            inputId = "UpsetTableGroupsGenes",
            label = "Select which contrasts you wish to see the overlapping genes/pathways between",
            choiceNames = names(gl1),
            choiceValues = seq_along(names(gl1))
        )
    })
})

observeEvent(
    {
        input$dataToUpset
        input$subdataToUpset
        input$upsetOrVenn
        input$upsetMode
        input$upsetMinimum
        input$upsetSetMode
        input$lfcUpset
        input$pTypeUpset
        input$pThresholdUpset
        lapply(seq_along(multiDeReactive$results), function(i) {
            input[[(paste0("multiGroupColour", names(multiDeReactive$results)[[i]]) %>% gsub(" ", "_", .))]]
        })
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE,
    {
        updateCheckboxGroupInput(session = session, inputId = "UpsetTableGroupsGenes", selected = character(0))
    }
)

observeEvent(
    {
        input$dataToUpset
        input$subdataToUpset
        input$upsetOrVenn
        input$upsetMode
        input$upsetMinimum
        input$upsetSetMode
        input$textSizeUpset
        input$figWidthUpset
        input$figHeightUpset
        input$upsetColour
        input$lfcUpset
        input$pTypeUpset
        input$pThresholdUpset
        input$UpsetTableGroupsGenes
        upsetGeneList$genesList
        lapply(seq_along(multiDeReactive$results), function(i) {
            input[[(paste0("multiGroupColour", names(multiDeReactive$results)[[i]]) %>% gsub(" ", "_", .))]]
        })
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE,
    {
        req(input$dataToUpset)
        tryCatch({
        dataToUse <- switch(input$dataToUpset,
            "Genes" = "genesList",
            "ORA" = "oraList",
            "GSEA" = "gseaList"
        )
        req(length(upsetGeneList[[dataToUse]]) >= 2)
        results <- multiDeReactive$results
        gl1 <- upsetGeneList[[dataToUse]][[input$subdataToUpset]]
        req(!is.null(gl1))

        # Colours: Use multiGroupColour
        colours <- setNames(lapply((names(gl1) %>% gsub(" ", "_", .)), function(k) {
            paste0(col2hex(input[[paste0("multiGroupColour", k)]]), "FF")
        }), names(gl1))
        colours <- unlist(colours)
        greyscale <- col2hex(c("grey87", "grey99"))

        # Makin' venns like its 2012 ----
        if (input$upsetOrVenn == "Venn" & length(results) <= 3) {
            multiUpsetReactive$plot <- ggvenn::ggvenn(gl1, fill_color = as.character(colours), set_name_size = input$textSizeUpset * 0.5, text_size = input$textSizeUpset * 0.4)
            output$UpsetPlot <- renderPlot(
                {
                    tryCatch(
                        multiUpsetReactive$plot,
                        error = function(e) {
                            if (inherits(e, "shiny.silent.error")) stop(e)
                            msg <- paste0("Venn render error: ", e$message,
                                          "\nCall: ", deparse(conditionCall(e)))
                            showNotification(msg, type = "error", duration = 30)
                            plot.new(); text(0.5, 0.5, msg, cex = 0.8)
                        }
                    )
                },
                width = as.numeric(input$figWidthUpset),
                height = as.numeric(input$figHeightUpset) * 2
            )

            if (length(input$UpsetTableGroupsGenes) >= 1) {
                m1 <- make_comb_mat(gl1, mode = "distinct")
                combViewList <- list()
                for (i in seq_along(gl1)) {
                    if (i %in% as.numeric(input$UpsetTableGroupsGenes)) {
                        combViewList[[i]] <- 1
                    } else {
                        combViewList[[i]] <- 0
                    }
                }
                combView <- paste0(unlist(combViewList), collapse = "")
                if (input$dataToUpset == "Genes") {
                    ec <- extract_comb(m1, as.character(combView))
                    if (length(ec) >= 1) {
                        output$UpsetIntersectGenesUI <- renderUI({
                            HTML(paste0(ec, collapse = "<br>"))
                        })
                        observeEvent(input$copy_link, {
                            text <- paste0(ec, collapse = "\n")
                            session$sendCustomMessage("txt", text)
                        })
                        jsCall <- paste0('enrich({list: "', paste(ec, collapse = "\\n"), '", popup: true});')
                        output$linkToEnrichr <- renderUI({
                            HTML(paste0("<a href='javascript:void(0)' onclick='", jsCall, "'>Analyse at Enrichr website</a>"))
                        })
                    } else {
                        paste("No genes in this selection.\n")
                    }
                } else {
                    output$UpsetIntersectGenesUI <- function() {
                        ec <- extract_comb(m1, as.character(combView))
                        if (length(ec) >= 1) {
                            kable(x = ec, format = "html", col.names = "Overlapping Pathways") %>%
                                kable_styling(bootstrap_options = "striped", full_width = FALSE, position = "left")
                        } else {
                            paste("No pathways in this selection.\n")
                        }
                    }
                }
            } else {
                output$UpsetIntersectGenesUI <- NULL
            }
        }

        # Makin' upsets like its 2022 ----
        if (input$upsetOrVenn == "Upset") {
            output$UpsetPlot <- NULL
            output$upsetError <- NULL
            if (!length(which(as.numeric(lengths(gl1)) >= 1)) >= 2) {
                output$upsetError <- renderUI({
                    h5("There are no overlapping DEGs/pathways based on current settings/input data.")
                })
            }
            if (length(which(as.numeric(lengths(gl1)) >= 1)) >= 2) {
                m1 <- make_comb_mat(gl1, mode = tolower(input$upsetMode))
                if (input$upsetSetMode == FALSE) {
                    m1 <- m1[comb_degree(m1) >= 2]
                }
                m1 <- m1[comb_size(m1) >= input$upsetMinimum]
                ss1 <- set_size(m1)
                cs1 <- comb_size(m1)
                g <- reactive({
                    cs1
                })
                combColsToUse <- switch(input$upsetColour,
                    # "colour" = brewer.pal(10, "RdYlBu")[comb_degree(m1)[rev(order(comb_size(m1)))]],
                    "colour" = "black",
                    "greyscale" = "black"
                )
                colsToUse <- switch(input$upsetColour,
                    "colour" = list(colorspace::lighten(colours, amount = 0.2)[rev(order(set_size(m1)))], colorspace::lighten(colours, amount = 0.2)[]),
                    "greyscale" = c(rep("#F0F0F0", nrow(m1)))
                )

                ht <- switch(input$upsetColour,
                    "colour" = {
                        ComplexHeatmap::UpSet(
                            m = m1,
                            column_title = paste(input$dataToUpset, names(dataToPlot[[input$dataToUpset]])[which(dataToPlot[[input$dataToUpset]] %in% input$subdataToUpset)]),
                            set_order = order(ss1),
                            comb_order = rev(order(comb_size(m1))),
                            row_title = "Contrast",
                            pt_size = unit((input$textSizeUpset / 3), "mm"),
                            lwd = unit((input$textSizeUpset / 4), "mm"),
                            comb_col = combColsToUse,
                            bg_pt_col = "gray92",
                            bg_col = colsToUse[[1]],
                            row_names_gp = gpar(fontsize = input$textSizeUpset),
                            top_annotation = upset_top_annotation(
                                m1,
                                add_numbers = TRUE,
                                numbers_gp = gpar(fontsize = input$textSizeUpset * 0.8),
                                axis_param = list(gp = gpar(fontsize = input$textSizeUpset * 0.8)),
                                annotation_name_gp = gpar(fontsize = input$textSizeUpset), annotation_name_rot = 90
                            ),
                            left_annotation = upset_left_annotation(
                                m1,
                                add_numbers = TRUE,
                                numbers_gp = gpar(fontsize = input$textSizeUpset * 0.8),
                                gp = gpar(fontsize = input$textSizeUpset * 0.8, fill = colsToUse[[2]]), axis_param = list(gp = gpar(fontsize = input$textSizeUpset * 0.8))
                            )
                        )
                    },
                    "greyscale" = {
                        ComplexHeatmap::UpSet(
                            m = m1,
                            column_title = paste(input$dataToUpset, names(dataToPlot[[input$dataToUpset]])[which(dataToPlot[[input$dataToUpset]] %in% input$subdataToUpset)]),
                            set_order = order(ss1),
                            comb_order = rev(order(comb_size(m1))),
                            row_title = "Contrast",
                            pt_size = unit((input$textSizeUpset / 3), "mm"),
                            lwd = unit((input$textSizeUpset / 4), "mm"),
                            top_annotation = upset_top_annotation(
                                m1,
                                add_numbers = TRUE,
                                numbers_gp = gpar(fontsize = input$textSizeUpset * 0.8),
                                axis_param = list(gp = gpar(fontsize = input$textSizeUpset * 0.8)),
                                annotation_name_gp = gpar(fontsize = input$textSizeUpset * 0.8), annotation_name_rot = 90
                            ),
                            left_annotation = upset_left_annotation(
                                m1,
                                add_numbers = TRUE,
                                numbers_gp = gpar(fontsize = input$textSizeUpset * 0.8),
                                axis_param = list(gp = gpar(fontsize = input$textSizeUpset * 0.8))
                            )
                        )
                    }
                )

                multiUpsetReactive$plot <- ht
                output$UpsetPlot <- renderPlot(
                    {
                        validate(
                            need(length(g()) > 0, "There are no overlapping DEGs/pathways based on current settings/input data.")
                        )
                        tryCatch(
                            draw(ht),
                            error = function(e) {
                                if (inherits(e, "shiny.silent.error")) stop(e)
                                msg <- paste0("Upset render error: ", e$message,
                                              "\nCall: ", deparse(conditionCall(e)))
                                showNotification(msg, type = "error", duration = 30)
                                plot.new(); text(0.5, 0.5, msg, cex = 0.8)
                            }
                        )
                    },
                    width = as.numeric(input$figWidthUpset),
                    height = as.numeric(input$figHeightUpset)
                ) # end of figure

                if (length(input$UpsetTableGroupsGenes) >= 1) {
                    combViewList <- list()
                    for (i in seq_along(gl1)) {
                        if (i %in% as.numeric(input$UpsetTableGroupsGenes)) {
                            combViewList[[i]] <- 1
                        } else {
                            combViewList[[i]] <- 0
                        }
                    }
                    combView <- paste0(unlist(combViewList), collapse = "")
                    if (input$dataToUpset == "Genes") {
                        ec <- extract_comb(m1, as.character(combView))
                        if (length(ec) >= 1) {
                            output$UpsetIntersectGenesUI <- renderUI({
                                HTML(paste0(ec, collapse = "<br>"))
                            })
                            # Update reactive val for clipboard
                            upsetGeneListToCopy(paste0(ec, collapse = "\n"))

                            jsCall <- paste0('enrich({list: "', paste(ec, collapse = "\\n"), '", popup: true});')
                            output$linkToEnrichr <- renderUI({
                                HTML(paste0("<a href='javascript:void(0)' onclick='", jsCall, "'>Analyse at Enrichr website</a>"))
                            })
                        } else {
                            paste("No genes in this selection.\n")
                        }
                    } else {
                        output$UpsetIntersectGenesUI <- function() {
                            ec <- extract_comb(m1, as.character(combView))
                            if (length(ec) >= 1) {
                                kable(x = ec, format = "html", col.names = "Overlapping Pathways") %>%
                                    kable_styling(bootstrap_options = "striped", full_width = FALSE, position = "left")
                            } else {
                                paste("No pathways in this selection.\n")
                            }
                        }
                    }
                } else {
                    output$UpsetIntersectGenesUI <- NULL
                }
            }
        }
        }, error = function(e) {
            if (inherits(e, "shiny.silent.error")) stop(e)
            msg <- paste0("Upset/Venn build error: ", e$message,
                          "\nCall: ", deparse(conditionCall(e)))
            showNotification(msg, type = "error", duration = 30)
            output$UpsetPlot <- renderPlot({
                plot.new(); text(0.5, 0.5, msg, cex = 0.8)
            })
        })
    }
)


output$dlPlot_multiUpset <- downloadHandler(
    filename = function() {
        generate_filename(
            base_name = paste0(input$filename_multiUpset, "_Plot"),
            design_name = "Multi",
            extension = tolower(input$downloadFormatMultiUpset)
        )
    },
    content = function(file) {
        req(multiUpsetReactive$plot)
        plot_obj <- multiUpsetReactive$plot
        
        # Determine if it's Upset or Venn to handle drawing method?
        # export_plot_standard uses print() by default, or special handling for ComplexHeatmap if specified?
        # R logic in export_plot_standard checks inherits(plot_object, "Heatmap") or "HeatmapList".
        # ComplexHeatmap::UpSet returns a HeatmapList object usually.
        # ggvenn returns a ggplot object.
        
        export_plot_standard(
            plot_object = plot_obj,
            filename = file,
            format = tolower(input$downloadFormatMultiUpset),
            width = as.numeric(input$figWidthUpset) / 75,
            height = as.numeric(input$figHeightUpset) / 75,
            dpi = as.numeric(input$dpiMultiUpset)
        )
    }
)
