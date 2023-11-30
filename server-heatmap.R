req(!is.null(inputDataReactive()$dataType))

if (inputDataReactive()$dataType == "proteomics") {
  output$heatmapProteomicsColumnSelectUI <- renderUI({
    selectInput(inputId = "heatmapProteomicsColumnSelect", label = "Select main factor to plot by", choices = inputDataReactive()$factors, selected = inputDataReactive()$factors[1], multiple = F)
  })
  groupingNameHeatmap <- eventReactive(input$heatmapProteomicsColumnSelect, {
    return(list(
      "gn" = input$heatmapProteomicsColumnSelect
    ))
  })
} else if (inputDataReactive()$dataType == "RNASeq") {
  groupingNameHeatmap <- reactive(return(list("gn" = inputDataReactive()$param$groupingName)))
}

if (inputDataReactive()$dataType == "proteomics") {
  output$clusterWarningProteomics <- renderUI({
    helpText(
      "Clustering by row/column may break the heatmap for proteomics data, due 
      to the high number of NA values in the count table.")
  })
  updateCheckboxInput(session = session, inputId = "clusterColsHeatmap", value = FALSE)
  updateCheckboxInput(session = session, inputId = "clusterRowsHeatmap", value = FALSE)
}

updateSelectInput(
  session = session, 
  inputId = "heatmapCounts", 
  choices = names(inputDataReactive()$countList), 
  selected = switch(inputDataReactive()$dataType, "RNASeq" = "Normalised + Log2", "proteomics" = "transformedData")
  )

updateSelectizeInput(
  session = session,
  inputId = "heatmapGenes",
  choices = inputDataReactive()$genes,
  selected = "",
  server = T
)

observeEvent({genesReactive()}, {
  output$geneBucket3 <- renderUI({
    bucket_list(
      header = "Drag and drop features in order to be plotted",
      group_name = "bucket_list_group",
      orientation = "horizontal",
      add_rank_list(
        text = "Include these features in this order",
        labels = genesReactive()$genes,
        input_id = "keepBucketHeatmap"),
      add_rank_list(
        text = "Exclude these features",
        labels = NULL,
        input_id = "excludeBucketHeatmap")
    )
  })
})

observeEvent(groupingNameHeatmap(), {
  output$heatmapBucket <- renderUI({
    bucket_list(
      header = "Drag and drop groups in order to be plotted",
      group_name = "hmBucket",
      orientation = "horizontal",
      add_rank_list(
        text = "Include these groups",
        labels = as.list(unique(levels(as.factor(inputDataReactive()$dataset[[groupingNameHeatmap()$gn]])))),
        input_id = "hmKeepBucket"),
      add_rank_list(
        text = "Exclude these groups",
        labels = NULL,
        input_id = "hmExcludeBucket")
    )
  })
})

updateCheckboxGroupInput(
  session = session,
  inputId = "heatmapFactor2",
  choices = inputDataReactive()$factors,
  selected = ""
)

if (inputDataReactive()$dataType == "proteomics") {
  observeEvent(
    input$heatmapProteomicsColumnSelect,{
      updateCheckboxGroupInput(session = session, inputId = "heatmapFactor2", selected = input$heatmapProteomicsColumnSelect)
  })
}

if (inputDataReactive()$dataType == "RNASeq") {
  updateSelectizeInput(
    session = session,
    inputId = "heatmapGoInput",
    choices = inputDataReactive()$allPathways,
    selected = "",
    server = T
  )
}

updateSelectInput(
  session = session,
  inputId = "heatmapBatch",
  choices = inputDataReactive()$factors,
  selected = NULL
)

if (inputDataReactive()$dataType == "RNASeq") {
  output$heatmapGOUI <- renderUI({
    helpText("You can also plot the expression of features annotated to a given GO term (again, irrespective of p-value from the DE test).")
    selectizeInput(
      inputId = "heatmapGoInput",
      choices = c(),
      label = "Or, select the GO Term you wish to plot"
    )
    numericInput(
      inputId = "heatmapGoExpr",
      label = "Minimum median feature expression for GO plot",
      value = 0.05,
      min = 0,
      max = 1e100
    )
  })
}

observeEvent(inputDataReactive(), {
  if (inputDataReactive()$dataType == "proteomics") {
    updateSelectInput(session = session, inputId = "pTypeHeatmap", selected = "Raw")
  } else {
    updateSelectInput(session = session, inputId = "pTypeHeatmap", selected = "FDR")
  }
})

observeEvent(
  {
    input$figHeightHeatmap
    input$figWidthHeatmap
    input$textSizeHeatmap
    input$heatmapBatch
    input$lfcHeatmap
    input$pTypeHeatmap
    input$pThresholdHeatmap
    input$heatmapScale
    input$heatmapLimitD
    input$heatmapLimitCHigh
    input$heatmapLimitCMid
    input$heatmapLimitCLow
    input$heatmapLimitCPalette
    input$heatmapLimitCPaletteRev
    input$heatmapGeneNumber
    input$hmKeepBucket
    input$heatmapFactor2
    input$heatmapCounts
    input$clusterRowsHeatmap
    input$clusterColsHeatmap
    input$colnamesHeatmap
    input$geneNamesHeatmap
    input$heatmapColourRed
    input$heatmapColourWhite
    input$heatmapColourBlue
    input$showClusterColDend
    input$showClusterRowDend
    lapply(seq_along(inputDataReactive()$factorLevels), function (i) {
      input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]])]]
    })
    input$contrastSelected
    groupingNameHeatmap()
  }, ignoreNULL = FALSE, ignoreInit = TRUE, {
    tryCatch({
      if (!is.null(input$hmKeepBucket)) {
        
        if (inputDataReactive()$dataType == "proteomics") {
          design <- input$contrastSelected
          seqAnnoHeatmap <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]
          seqAnnoHeatmap[["usedInTest"]] = TRUE
        } else {
          seqAnnoHeatmap <- inputDataReactive()$seqAnno
          param <- inputDataReactive()$param
          design <- inputDataReactive()$design
        }
        
        output$heatmapDesign <- renderText({design})

        # Get the desired p-value column
        if (input$pTypeHeatmap == "Raw") {
          pTypeHeatmap <- "pValue"
        } else {
          pTypeHeatmap <- "fdr"
        }

        # get the metadata for plotting
        datasetHeatmap <- inputDataReactive()$dataset
        datasetHeatmap <- datasetHeatmap[which(datasetHeatmap[[groupingNameHeatmap()$gn]] %in% input$hmKeepBucket), ]
        
        # Set the factor levels per the bucket list order
        datasetHeatmap[[groupingNameHeatmap()$gn]] <- factor(datasetHeatmap[[groupingNameHeatmap()$gn]], input$hmKeepBucket)
        datasetHeatmap <- datasetHeatmap[order(datasetHeatmap[[groupingNameHeatmap()$gn]]), ]

        # Subset count matrix:
        countsHeatmap <- inputDataReactive()$countList[[input$heatmapCounts]]
        if (!is.null(input$heatmapBatch)) {
          for (i in seq_along(input$heatmapBatch)) {
            countsHeatmap <- limma::removeBatchEffect(countsHeatmap, datasetHeatmap[[input$heatmapBatch[i]]])
          }
        }
        countsHeatmap <- countsHeatmap[, rownames(datasetHeatmap)]
        if (input$heatmapScale == "Diverging") {
          if (input$heatmapCounts %in% c("FPKM", "TPM", "Normalised")) {
            countsHeatmap <- log2(countsHeatmap + param$backgroundExpression)
          }
          # Subtract the gene row mean from each sample
          countsHeatmap <- sweep(countsHeatmap, 1, rowMeans(countsHeatmap, na.rm = T))
          heatmapColours <- colorRamp2(
            breaks = c(input$heatmapLimitD, 0, -input$heatmapLimitD),
            colors = c(input$heatmapColourRed, input$heatmapColourWhite, input$heatmapColourBlue)
          )
          at <- c(-input$heatmapLimitD, 0, input$heatmapLimitD)
        } else if (input$heatmapScale == "Continuous") {
          countsHeatmap[countsHeatmap > input$heatmapLimitCHigh] <- input$heatmapLimitCHigh
          if (input$heatmapLimitCPaletteRev) {
            heatmapColours <- colorRampPalette(brewer.pal(9, input$heatmapLimitCPalette))(64)
          } else {
            heatmapColours <- rev(colorRampPalette(brewer.pal(9, input$heatmapLimitCPalette))(64))
          }
          at <- c(
            as.numeric(input$heatmapLimitCLow),
            as.numeric(input$heatmapLimitCMid), 
            as.numeric(input$heatmapLimitCHigh)
          )
        }
        if (inputDataReactive()$dataType == "proteomics") {
          rownames(countsHeatmap) <- gsub("\\~.*", "", rownames(countsHeatmap))
        }
        countsHeatmap <- countsHeatmap[match(seqAnnoHeatmap$gene_name, rownames(countsHeatmap)), ]
        countsHeatmap <- countsHeatmap[!is.na(rownames(countsHeatmap)), ]

        geneDirections <- list(
          "Both" = seqAnnoHeatmap$gene_name[which(
            as.numeric(seqAnnoHeatmap[[pTypeHeatmap]]) <= as.numeric(input$pThresholdHeatmap) &
              as.numeric(abs(seqAnnoHeatmap[["log2Ratio"]])) >= as.numeric(input$lfcHeatmap) &
              seqAnnoHeatmap[["usedInTest"]] == TRUE
          )],
          "Up" = seqAnnoHeatmap$gene_name[which(
            as.numeric(seqAnnoHeatmap[[pTypeHeatmap]]) <= as.numeric(input$pThresholdHeatmap) &
              as.numeric(seqAnnoHeatmap[["log2Ratio"]]) >= as.numeric(input$lfcHeatmap) &
              seqAnnoHeatmap[["usedInTest"]] == TRUE
          )],
          "Down" = seqAnnoHeatmap$gene_name[which(
            as.numeric(seqAnnoHeatmap[[pTypeHeatmap]]) <= as.numeric(input$pThresholdHeatmap) &
              as.numeric(seqAnnoHeatmap[["log2Ratio"]]) <= -as.numeric(input$lfcHeatmap) &
              seqAnnoHeatmap[["usedInTest"]] == TRUE
          )]
        )


        # Top Annotation
        # I realise this is the ugliest code in the world and I will refactor it
        # but just getting it to work at all was a minor miracle
        # HA = heatmap annotation
        ha <- paste0("HeatmapAnnotation('", groupingNameHeatmap()$gn, "' = datasetHeatmap[['", groupingNameHeatmap()$gn,"']])")
        if (!is.null(input$heatmapFactor2)) {
          p1 <- paste(lapply(input$heatmapFactor2, function(g) {
            paste0("'", g, "' = datasetHeatmap[['", g, "']]")
          }), collapse = ",")
          ha <- paste0("HeatmapAnnotation('", groupingNameHeatmap()$gn, "' = datasetHeatmap[['", groupingNameHeatmap()$gn,"']],", p1,")")
        }
        ha <- eval(parse(text = ha))
        
        colourListHeatmap <- setNames(lapply(inputDataReactive()$factors, function(f) {
          setNames(lapply(unique(datasetHeatmap[[f]]), function(k) {
            paste0(col2hex(input[[paste0("GroupColour", f, "__", k)]]), "FF")
          }), levels(datasetHeatmap[[f]]))
        }), inputDataReactive()$factors)
        
        for (f in c(groupingNameHeatmap()$gn, input$heatmapFactor2)) {
          for (l in levels(as.factor(datasetHeatmap[,f]))) {
            ha@anno_list[[f]]@color_mapping@full_col[[l]] <- colourListHeatmap[[f]][[l]]
            ha@anno_list[[f]]@fun@var_env$color_mapping@full_col[[l]] <- colourListHeatmap[[f]][[l]]
          }
        }
        # Draw the heatmaps
        lapply(c("Both", "Up", "Down"), function(sig) {
          countsHeatmapSig <- countsHeatmap
          countsHeatmapSig <- countsHeatmapSig[geneDirections[[sig]], ]
          
          # rowSumsVector <- rowSums(countsHeatmapSig, na.rm = T)
          # zeroGeneMeans <- rownames(countsHeatmapSig)[which(rowSums(countsHeatmapSig, na.rm = T) == 0)]
          # if (length(zeroGeneMeans) >= 1) {
          #   output[[paste0(sig, "HeatmapZeroMeanMessage")]] <- renderText({
          #     c("The following features have been excluded because they have no mean calculated: ", zeroGeneMeans)
          #   })
          # }
          # countsHeatmapSig <- countsHeatmapSig[which(!rownames(countsHeatmapSig) %in% zeroGeneMeans),]
          if (as.numeric(input$heatmapGeneNumber) <= length(geneDirections[[sig]])) {
            countsHeatmapSig <- countsHeatmapSig[1:as.numeric(input$heatmapGeneNumber), ]
          } 
          
          heatmap <- ComplexHeatmap::Heatmap(
            matrix = countsHeatmapSig,
            top_annotation = ha,
            name = input$heatmapCounts %>% gsub("\\+", "\n\\+", .),
            cluster_columns = input$clusterColsHeatmap,
            cluster_rows = input$clusterRowsHeatmap,
            show_column_names = input$colnamesHeatmap,
            show_row_names = input$geneNamesHeatmap,
            show_column_dend = input$showClusterColDend, 
            show_row_dend = input$showClusterRowDend, 
            use_raster = TRUE,
            row_dend_width = unit(3, "cm"),
            column_dend_height = unit(2, "cm"),
            col = heatmapColours,
            # heatmap_legend_param = list(at = at, legend_width = unit(14, "mm")),
            column_title = paste(
              "Top", nrow(countsHeatmapSig), sig, "Regulated DEGs:"
              # param$sampleGroup, "vs", param$refGroup
              )
          )

          # Output heatmaps
          output[[paste0("heatmap", sig)]] <- renderPlot(
            {
              draw(heatmap, merge_legend = TRUE)
            },
            width = as.numeric(input$figWidthHeatmap),
            height = as.numeric(input$figHeightHeatmap)
          )

          # download button for the heatmaps
          output[[paste0("dlHeatmap", sig, "Button")]] <- downloadHandler(
            filename = function() {
              paste0(design, "_", sig, "_heatmap.pdf")
            },
            content = function(file) {
              pdf(file = file, width = as.numeric(input$figWidthHeatmap/70), height = as.numeric(input$figHeightHeatmap/70))
              draw(heatmap, merge_legend = TRUE)
              dev.off()
            }
          )
          output[[paste0("dlHeatmap", sig, "ButtonPNG")]] <- downloadHandler(
            filename = function() {
              paste0(design, "_", sig, "_heatmap.png")
            },
            content = function(file) {
              png(filename = file, width = as.numeric(input$figWidthHeatmap)/60, height = as.numeric(input$figHeightHeatmap)/60, units = "in", res = 600)
              draw(heatmap, merge_legend = TRUE)
              dev.off()
            }
          )
          output[[paste0("dlHeatmap", sig, "DFButton")]] <- downloadHandler(
            filename = function() {
              paste0(design, "_", sig, "_heatmap.xlsx")
            },
            content = function(file) {
              openxlsx::write.xlsx((as.data.frame(countsHeatmapSig) %>% rownames_to_column("gene_name")), file = file)
            }
          )

          # Show number of genes on plot (changes depending on settings applied)
          output[[paste0("hGN", sig)]] <- renderText({
            paste0("Number of features per settings chosen: ", nrow(heatmapResults[[sig]]$matrix))
          })

        }) # end of lapply

        # Make a custom heatmap from the genes selected
        observeEvent({
          input$keepBucketHeatmap
        }, {
          
          hmCustomMatrix <- inputDataReactive()$countList[[input$heatmapCounts]]
          if (!is.null(input$heatmapBatch)) {
            for (i in seq_along(input$heatmapBatch)) {
              hmCustomMatrix <- limma::removeBatchEffect(hmCustomMatrix, datasetHeatmap[[input$heatmapBatch[i]]])
            }
          }
          hmCustomMatrix <- hmCustomMatrix[, rownames(datasetHeatmap)]
          if (input$heatmapScale == "Diverging") {
            if (input$heatmapCounts %in% c("FPKM", "TPM", "Normalised")) {
              hmCustomMatrix <- log2(hmCustomMatrix + param$backgroundExpression)
            }
            hmCustomMatrix <- sweep(hmCustomMatrix, 1 , rowMeans(hmCustomMatrix, na.rm = T))
          }

          if(!is.null(input$keepBucketHeatmap)) {
            hmCustomMatrix <- hmCustomMatrix[rownames(hmCustomMatrix) %in% input$keepBucketHeatmap, ] %>% as.matrix()
            zeroGeneMeansCustom <- rownames(hmCustomMatrix)[which(rowSums(hmCustomMatrix, na.rm = T) == 0)]
            # if (length(zeroGeneMeansCustom) >= 1) {
            #   output[[paste0("CustomHeatmapZeroMeanMessage")]] <- renderText({
            #     c("The following features have been excluded because they have no mean calculated: ", zeroGeneMeansCustom)
            #   })
            # }
            hmCustomMatrix <- hmCustomMatrix[which(!rownames(hmCustomMatrix) %in% zeroGeneMeansCustom),]

            if (input$heatmapScale == "Diverging") {
              hmCustomColours <- colorRamp2(
                breaks = c(
                  input$heatmapLimitD,
                  0,
                  -input$heatmapLimitD
                ),
                colors = c(
                  input$heatmapColourRed,
                  input$heatmapColourWhite,
                  input$heatmapColourBlue
                )
              )
            } else if (input$heatmapScale == "Continuous") {
              hmCustomMatrix[hmCustomMatrix > input$heatmapLimitC] <- input$heatmapLimitC
              hmCustomColours <- rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(64))
            }
            
            # Output heatmaps
            output$heatmapCustom <- renderPlot(
              {
                ch <- ComplexHeatmap::Heatmap(
                  matrix = hmCustomMatrix,
                  top_annotation = ha,
                  name = input$heatmapCounts,
                  cluster_columns = input$clusterColsHeatmap,
                  cluster_rows = input$clusterRowsHeatmap,
                  show_column_names = input$colnamesHeatmap,
                  show_row_names = input$geneNamesHeatmap,
                  col = hmCustomColours,
                  column_title = paste("Custom heatmap")
                )
                draw(ch, merge_legend = TRUE)
              },
              width = as.numeric(input$figWidthHeatmap),
              height = as.numeric(input$figHeightHeatmap)
            )
            output[[paste0("dlHeatmapCustomButton")]] <- downloadHandler(
              filename = function() {
                paste0(design, "_Custom_heatmap.pdf")
              },
              content = function(file) {
                pdf(file = file, width = as.numeric(input$figWidthHeatmap)/70, height = as.numeric(input$figHeightHeatmap)/70)
                draw(ComplexHeatmap::Heatmap(
                  matrix = hmCustomMatrix,
                  top_annotation = ha,
                  name = input$heatmapCounts,
                  cluster_columns = input$clusterColsHeatmap,
                  cluster_rows = input$clusterRowsHeatmap,
                  show_column_names = input$colnamesHeatmap,
                  col = hmCustomColours,
                  column_title = paste("Custom heatmap")
                ), merge_legend = TRUE)
                dev.off()
              }
            )
          }
        })
        output[[paste0("dlHeatmapCustomDFButton")]] <- downloadHandler(
          filename = function() {
            paste0(design, "_Custom_heatmap.xlsx")
          },
          content = function(file) {
            openxlsx::write.xlsx((as.data.frame(hmCustomMatrix) %>% rownames_to_column("gene_name")), file = file)
          }
        )

        # heatmapGoInput
        if (inputDataReactive()$dataType == "RNASeq") {
          observeEvent({
            input$heatmapGoInput
            input$heatmapGoExpr
          }, ignoreInit = T, {
            withProgress(message = "Making GO Heatmap. Can take a few minutes!", {
              goToPlot <- gsub(" .*", "", input$heatmapGoInput)
              goGenes <- seqAnnoHeatmap$gene_name[with(seqAnnoHeatmap, grepl(goToPlot, paste(`GO BP`, `GO MF`, `GO CC`)))]

              if(!is.null(goGenes)) {
                hmGOMatrix <- countsHeatmap[rownames(countsHeatmap) %in% goGenes, ] %>% as.matrix()

                if (input$heatmapScale == "Diverging") {
                  hmGOMatrix <- hmGOMatrix - rowMeans(hmGOMatrix)
                  if (input$heatmapCounts %in% c("FPKM", "TPM", "Normalised")) {
                    hmGOMatrix <- log2(hmGOMatrix+param$backgroundExpression)
                  }
                  hmCustomColours <- colorRamp2(
                    breaks = c(
                      input$heatmapLimitD,
                      0,
                      -input$heatmapLimitD
                    ),
                    colors = c(
                      input$heatmapColourRed,
                      input$heatmapColourWhite,
                      input$heatmapColourBlue
                    )
                  )
                } else if (input$heatmapScale == "Continuous") {
                  hmGOMatrix[hmGOMatrix > input$heatmapLimitC] <- input$heatmapLimitC
                  hmCustomColours <- rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(64))
                }
                hmGOMatrix <- hmGOMatrix[abs(rowMedians(hmGOMatrix)) >= input$heatmapGoExpr, ]

                # Output heatmaps
                output$heatmapGO <- renderPlot(
                  {
                    gh <- ComplexHeatmap::Heatmap(
                      matrix = hmGOMatrix,
                      top_annotation = ha,
                      name = input$heatmapCounts,
                      cluster_columns = input$clusterColsHeatmap,
                      cluster_rows = input$clusterRowsHeatmap,
                      show_column_names = input$colnamesHeatmap,
                      show_row_names = input$geneNamesHeatmap,
                      col = hmCustomColours,
                      column_title = sub("\\s+", "\n", input$heatmapGoInput)
                    )
                    draw(gh, merge_legend = TRUE)
                  },
                  width = as.numeric(input$figWidthHeatmap),
                  height = as.numeric(input$figHeightHeatmap)
                )
                output[[paste0("dlHeatmapGOButton")]] <- downloadHandler(
                  filename = function() {
                    paste0(design, "_GO_heatmap.pdf")
                  },
                  content = function(file) {
                    pdf(file = file, width = as.numeric(input$figWidthHeatmap)/70, height = as.numeric(input$figHeightHeatmap)/70)
                    draw(ComplexHeatmap::Heatmap(
                      matrix = hmGOMatrix,
                      top_annotation = ha,
                      name = input$heatmapCounts,
                      cluster_columns = input$clusterColsHeatmap,
                      cluster_rows = input$clusterRowsHeatmap,
                      show_column_names = input$colnamesHeatmap,
                      col = hmCustomColours,
                      column_title = sub("\\s+", "\n", input$heatmapGoInput)
                    ), merge_legend = TRUE)
                    dev.off()
                  }
                )
                output[[paste0("dlHeatmapGODFButton")]] <- downloadHandler(
                  filename = function() {
                    paste0(design, "_GO_heatmap.xlsx")
                  },
                  content = function(file) {
                    openxlsx::write.xlsx((as.data.frame(hmGOMatrix) %>% rownames_to_column("gene_name")), file = file)
                  }
                )
              }
            })
          })
        }

      } # here
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
) # End of heatmaps observeEvent
