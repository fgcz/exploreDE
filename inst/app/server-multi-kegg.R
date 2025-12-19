# Server logic for Multi KEGG Plot

observeEvent(
    {
        input$KEGGID
        input$lfcKEGG
        input$pTypeKEGG
        input$pThresholdKEGG
        input$KEGGLimitColour
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE,
    {
        req(multiDeReactive$results)
        req(length(multiDeReactive$results) > 0)
        results <- multiDeReactive$results

        # Determine species from first contrast
        # exploreDE stores refBuild in param
        refBuild <- results[[1]]$param$refBuild
        if (is.null(refBuild) && !is.null(results[[1]]$param$ezRef)) {
            refBuild <- results[[1]]$param$ezRef@refBuild
        }

        sp <- "other"
        if (!is.null(refBuild)) {
            # Simple check for human/mouse
            if (grepl("Homo_sapiens|GRCh|hg", refBuild, ignore.case = TRUE)) {
                sp <- "hsa"
            } else if (grepl("Mus_musculus|GRCm|mm", refBuild, ignore.case = TRUE)) {
                sp <- "mmu"
            }
        }

        if (sp == "other") {
            showModal(modalDialog(
                title = "Your data is not human or mouse.",
                paste0("KEGG Pathways are only compatible with human and mouse data, unfortunately."),
                easyClose = TRUE,
                footer = NULL
            ))
        }

        if (sp != "other") {
            req(input$KEGGID)
            if (nchar(input$KEGGID) >= 5) { # Allow short IDs like 05200
                if (input$pTypeKEGG == "Raw") {
                    pTypeKEGG <- "pValue"
                } else {
                    pTypeKEGG <- "fdr"
                }

                pathviewList <- list()
                for (i in seq_along(results)) {
                    df1 <- results[[i]]$seqAnno
                    # Ensure we have entrez_id
                    if (!"entrez_id" %in% colnames(df1)) {
                        # Try to map if missing? For now assume it's there or fail gracefully
                        next
                    }

                    df1 <- df1[order(df1$pValue), ]
                    # df1 <- df1[df1$usedInTest == TRUE,] # exploreDE might not have this column, skip

                    # Filter
                    df1 <- df1[which(df1[[pTypeKEGG]] <= as.numeric(input$pThresholdKEGG) & abs(df1$log2Ratio) >= as.numeric(input$lfcKEGG)), ]

                    gL1 <- data.frame("Gene" = as.character(df1$entrez_id), "log2Ratio" = df1$log2Ratio)
                    colnames(gL1)[2] <- names(results)[i]
                    gL1 <- gL1[order(gL1[, 2], decreasing = TRUE), ]
                    gL1 <- gL1[which(!gL1$Gene == ""), ]
                    gL1 <- gL1[which(!is.na(gL1$Gene)), ]
                    pathviewList[[names(results)[i]]] <- gL1
                }

                if (length(pathviewList) > 0) {
                    pvDF <- merge(pathviewList[[1]], pathviewList[[2]], all = TRUE)
                    if (length(results) >= 3) {
                        for (i in 3:length(results)) {
                            if (names(results)[i] %in% names(pathviewList)) {
                                pvDF <- merge(pvDF, pathviewList[[names(results)[i]]], all = TRUE)
                            }
                        }
                    }
                    pvDF <- pvDF[!duplicated(pvDF$Gene), ]
                    rownames(pvDF) <- NULL
                    pvDF <- column_to_rownames(pvDF, "Gene")
                    pvDF <- as.matrix(pvDF)

                    outfile <- paste0(input$KEGGID, ".pathview.multi.png")

                    output$KEGGout <- renderImage(
                        {
                            # Use temporary directory to avoid Docker permissions issues
                            temp_dir <- tempdir()
                            old_wd <- getwd()

                            tryCatch({
                                # Change to temp directory where we have write permissions
                                setwd(temp_dir)

                                # Check if pathview is loaded
                                if (!requireNamespace("pathview", quietly = TRUE)) {
                                    stop("Package 'pathview' is not installed.")
                                }

                                pathview::pathview(
                                    gene.data = pvDF,
                                    gene.idtype = "ENTREZ",
                                    pathway.id = input$KEGGID,
                                    key.pos = "bottomright",
                                    species = sp,
                                    low = list("gene" = "dodgerblue3"),
                                    mid = list("gene" = "lightgrey"),
                                    high = list("gene" = "indianred3"),
                                    kegg.native = TRUE,
                                    limit = list(gene = input$KEGGLimitColour),
                                    map.null = FALSE,
                                    bins = 12
                                )

                                # Get the full path to the generated file
                                # pathview generates file with .png extension appended to outfile name usually
                                # actually pathview generates [pathway.id].pathview.multi.png

                                generated_file <- paste0(input$KEGGID, ".pathview.multi.png")
                                output_file <- file.path(temp_dir, generated_file)

                                if (file.exists(output_file)) {
                                    list(src = output_file, contentType = "image/png")
                                } else {
                                    # Fallback or error
                                    list(src = "", alt = "Plot generation failed")
                                }
                            }, error = function(e) {
                                cat("KEGG plot error:", conditionMessage(e), "\n")
                                list(src = "", alt = paste("Error:", conditionMessage(e)))
                            }, finally = {
                                # Always restore working directory
                                setwd(old_wd)
                            })
                        },
                        deleteFile = TRUE
                    )
                }

                output$KEGGOrder <- function() {
                    d1 <- lapply(1:length(results), function(i) {
                        names(results)[i]
                    }) %>% data.frame()
                    colnames(d1) <- paste("Contrast", 1:length(results))
                    d1 %>%
                        kable("html", align = "l") %>%
                        kable_styling(
                            bootstrap_options = "striped",
                            full_width = FALSE,
                            position = "left"
                        )
                }
            }
        }
    }
)

output$dlPlot_multiKegg <- downloadHandler(
    filename = function() {
        generate_filename(
            base_name = input$filename_multiKegg,
            design_name = "Multi",
            extension = "png"
        )
    },
    content = function(file) {
        req(multiDeReactive$results)
        req(input$KEGGID)
        
        # Determine species
        results <- multiDeReactive$results
        refBuild <- results[[1]]$param$refBuild
        if (is.null(refBuild) && !is.null(results[[1]]$param$ezRef)) {
            refBuild <- results[[1]]$param$ezRef@refBuild
        }
        
        sp <- "other"
        if (!is.null(refBuild)) {
            if (grepl("Homo_sapiens|GRCh|hg", refBuild, ignore.case = TRUE)) sp <- "hsa"
            else if (grepl("Mus_musculus|GRCm|mm", refBuild, ignore.case = TRUE)) sp <- "mmu"
        }
        
        req(sp != "other")
        
        # Prepare data matrix (this logic is duplicated from renderImage, but safer to re-run for download)
        if (input$pTypeKEGG == "Raw") pTypeKEGG <- "pValue" else pTypeKEGG <- "fdr"
        
        pathviewList <- list()
        for (i in seq_along(results)) {
            df1 <- results[[i]]$seqAnno
            if (!"entrez_id" %in% colnames(df1)) next
            
            df1 <- df1[which(df1[[pTypeKEGG]] <= as.numeric(input$pThresholdKEGG) & abs(df1$log2Ratio) >= as.numeric(input$lfcKEGG)), ]
            
            gL1 <- data.frame("Gene" = as.character(df1$entrez_id), "log2Ratio" = df1$log2Ratio)
            colnames(gL1)[2] <- names(results)[i]
            gL1 <- gL1[which(!is.na(gL1$Gene) & gL1$Gene != ""), ]
            pathviewList[[names(results)[i]]] <- gL1
        }
        
        req(length(pathviewList) > 0)
        
        pvDF <- merge(pathviewList[[1]], pathviewList[[2]], all = TRUE)
        if (length(results) >= 3) {
            for (i in 3:length(results)) {
                if (names(results)[i] %in% names(pathviewList)) {
                    pvDF <- merge(pvDF, pathviewList[[names(results)[i]]], all = TRUE)
                }
            }
        }
        pvDF <- pvDF[!duplicated(pvDF$Gene), ]
        rownames(pvDF) <- NULL
        pvDF <- tibble::column_to_rownames(pvDF, "Gene")
        pvDF <- as.matrix(pvDF)
        
        # Generate plot in temp dir
        temp_dir <- tempdir()
        old_wd <- getwd()
        setwd(temp_dir)
        
        tryCatch({
            pathview::pathview(
                gene.data = pvDF,
                gene.idtype = "ENTREZ",
                pathway.id = input$KEGGID,
                key.pos = "bottomright",
                species = sp,
                low = list("gene" = "dodgerblue3"),
                mid = list("gene" = "lightgrey"),
                high = list("gene" = "indianred3"),
                kegg.native = TRUE,
                limit = list(gene = input$KEGGLimitColour),
                map.null = FALSE,
                bins = 12
            )
            
            generated_file <- paste0(input$KEGGID, ".pathview.multi.png")
            
            if (file.exists(generated_file)) {
                file.copy(generated_file, file)
            }
        }, finally = {
            setwd(old_wd)
        })
    }
)
