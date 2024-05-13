if(inputDataReactive()$dataType == "RNASeq") {
  param <- inputDataReactive()$param
  
  observeEvent(input$tabs, {
    if (input$tabs == "oraTab") {
      if (is_empty(inputDataReactive()$ora)) {
        showModal(modalDialog(title = "No ORA Results", "No pathway analysis was performed so this tab is empty."))
      }
    }
  })
  
  if(isTRUE(param$runGO & param$featureLevel != "isoform")) {
    ora <- inputDataReactive()$ora
    oraHTML <- inputDataReactive()$oraHTML
    se <- inputDataReactive()$se
    design <- inputDataReactive()$design
    
    output$oraDesign <- renderText({design})
    
    # Download full ORA results file
    output$dlTableORA <- downloadHandler(
      filename = function() {
        basename(oraHTML)
      },
      content = function(file) {
        file.copy(from = oraHTML, to = file)
      }
    )
    
    # ORA Overview Table
    ora.overview <- ezFrame(
      "Number of Significant Pathways" = numeric()
    )
    for (goDomain in names(ora)) {
      for (resName in names(ora[[goDomain]])) {
        testName <- paste0(goDomain, ": ", resName)
        ora.overview[testName, "Number of Significant Pathways"] <-
          nrow(ora[[goDomain]][[resName]]@result[
            as.numeric(
              ora[[goDomain]][[resName]]@result$p.adjust
            ) < param$fdrThreshORA,
          ])
      }
    }
    output$oraOverview <- function() {
      ora.overview %>%
        kable(
          format = "html",
          caption =
            paste("ORA GO:")
        ) %>%
        kable_styling(
          bootstrap_options = "striped",
          full_width = FALSE,
          position = "left"
        )
    }
    # ORA thresholds
    oraThresholds <- ezFrame(
      "ORA Cut-Offs" = character()
    )
    oraThresholds["Input DEG P-value", "ORA Cut-Offs"] <- param$pValThreshGO
    oraThresholds["Input DEG Log2 Ratio", "ORA Cut-Offs"] <- param$log2RatioThreshGO
    oraThresholds["Input up-regulated DEG total", "ORA Cut-Offs"] <- length(metadata(se)$enrichInput$selections$upGenes)
    oraThresholds["Input down-regulated DEG total", "ORA Cut-Offs"] <- length(metadata(se)$enrichInput$selections$downGenes)
    oraThresholds["Output ORA Term FDR", "ORA Cut-Offs"] <- param$fdrThreshORA
    output$oraThreshOverview <- function() {
      oraThresholds %>%
        kable(
          format = "html",
          caption =
            paste("Cut-Offs:")
        ) %>%
        kable_styling(
          bootstrap_options = "striped",
          full_width = FALSE,
          position = "left"
        )
    }
    
    observeEvent({
      input$oraType
      input$oraDirection
      }, {
      er <- ora[[input$oraType]][[input$oraDirection]]
      if (input$tabs == "oraTab") {
        if (is.null(er) || nrow(er@result) == 0) {
          showModal(modalDialog(title = "No ORA Results", "There are no results for this selection."))
        }
      }
      if (!is.null(er) && nrow(er@result) != 0) {
        er@result$Label <- ""
        for (i in seq_along(er@result$Description)) {
          desc <- er@result$Description[i]
          label <- paste0(er@result$ID[i], "\n", desc)
          if (!is.na(label) && nchar(label) > 45) {
            label <- substr(label, 1, 45)
          }
          er@result$Label[i] <- label
        }
        output$selectedTable_ORA <- DT::renderDataTable({
          oraDfToShow <- as.data.frame(er)[, c("Label", "GeneRatio", "geneName")]
          oraDfToShow$geneNumber <- sapply(strsplit(as.character(oraDfToShow$GeneRatio), "/"), function(x) as.numeric(x[1]))
          
          DT::datatable(
            data = oraDfToShow,
            filter = "top", class = "cell-border stripe",
            rownames = FALSE, caption = "Click pathways in this table to add them to the figures on the right",
            options = list(
              columnDefs = list(
                # Hide the numeric column used for sorting
                list(visible = FALSE, targets = c(grep("geneNumber", colnames(oraDfToShow))-1)),
                # Specify the column you want to display and make sortable based on the hidden one
                list(orderData = c(grep("geneNumber", colnames(oraDfToShow))-1), targets = c(grep("GeneRatio", colnames(oraDfToShow))-1))
              )
            ),
          ) %>%
            DT::formatStyle(
              columns = colnames(.$x$data), `font-size` = "12px"
            )
        })
      }
    })
    
    oraPlotList <- eventReactive({
      input$textSizeORA
      input$nodeSizeORA
      input$selectedTable_ORA_rows_selected
      input$showGeneLabelsORA
    }, {
      
      req(!is.null(input$selectedTable_ORA_rows_selected))
      
      er <- ora[[input$oraType]][[input$oraDirection]]
      if (!is.null(er) && nrow(er@result) >=2) {
        er@result$Label <- ""
        for (i in seq_along(er@result$Description)) {
          desc <- er@result$Description[i]
          label <- paste0(er@result$ID[i], "\n", desc)
          if (!is.na(label) && nchar(label) > 45) {
            label <- substr(label, 1, 45)
          }
          er@result$Label[i] <- label
        }
        
        # Get log2FCs of genes for network plot
        sigGeneIds <- metadata(se)$enrichInput$selections[[input$oraDirection]]
        log2RatioORA <- metadata(se)$enrichInput$seqAnno$log2Ratio[metadata(se)$enrichInput$seqAnno$gene_id %in% sigGeneIds]
        names(log2RatioORA) <- metadata(se)$enrichInput$seqAnno$gene_name[metadata(se)$enrichInput$seqAnno$gene_id %in% sigGeneIds]
        log2RatioORA <- sort(log2RatioORA, decreasing = T)
        er@result$geneID <- er@result$geneName
        nl <- "category"
        if (input$showGeneLabelsORA) {
          nl <- "all"
        }
        
        if (!is.null(er) && nrow(er@result) >= 2) {
          cn <- enrichplot::cnetplot(
            x = er,
            color.params = list(foldChange = log2RatioORA),
            showCategory = er@result$Description[input$selectedTable_ORA_rows_selected],
            node_label = nl,
            cex.params = list(category_label = (input$textSizeORA / 15), gene_label = (input$textSizeORA / 18), gene_node = (input$textSizeORA / 16), category_node = (input$textSizeORA / 13))
          )
          if (input$oraDirection == "upGenes") {
            cn <- cn + scale_color_gradientn(name = "fold change", colours = c("white", "firebrick3"))
          } else if (input$oraDirection == "downGenes") {
            cn <- cn + scale_color_gradientn(name = "fold change", colours = c("deepskyblue4", "white"))
          } else if (input$oraDirection == "bothGenes") {
            cn <- cn + scale_colour_gradient2(name = "fold change", low = "deepskyblue4", mid = "white", high = "firebrick3", midpoint = 0)
          }
          hp <- enrichplot::heatplot(
            x = er,
            foldChange = log2RatioORA,
            showCategory = er@result$Description[input$selectedTable_ORA_rows_selected],
          ) + scale_colour_distiller(palette = "RdBu")
          bp <- barplot(
            er,
            showCategory = er@result$Description[input$selectedTable_ORA_rows_selected],
            title = paste(input$oraType, input$oraDirection)
            )
          bp <- bp +
            scale_y_discrete(
              labels = rev(
                er@result$Label[
                  match(bp$plot_env$df$ID, er@result$ID)])) +
            theme(
              axis.text.y = element_text(
                vjust = -0.01, size = input$textSizeORA))
          dp <- enrichplot::dotplot(
            er,
            x = "GeneRatio",
            showCategory = er@result$Description[input$selectedTable_ORA_rows_selected],
            title = paste(input$oraType, input$oraDirection)
          )
          dp <- dp +
            scale_y_discrete(
              labels = rev(
                er@result$Label[
                  match(dp$plot_env$df$ID, er@result$ID)])) +
            theme(
              axis.text.y = element_text(
                vjust = -0.01, size = input$textSizeORA))
        }
        
        return(list(
          "cnetplot" = cn,
          "barplot" = bp,
          "dotplot" = dp,
          "heatmap" = hp
        ))
      }
    })
    
    # Custom ORA results output
    observeEvent(
      {
        input$figHeightORA
        input$figWidthORA
        oraPlotList()
      }, ignoreNULL = FALSE, {
        
        # Network plot
        output$cnetPlot_ORA <- renderPlot(
          {
            oraPlotList()[["cnetplot"]]
          },
          width = as.numeric(input$figWidthORA),
          height = as.numeric(input$figHeightORA)
        )
        output$dlCnetPlot_ORA <- downloadHandler(
          filename = function() {
            paste0(design, "ORA_", input$oraType, "_network.pdf")
          },
          content = function(file) {
            pdf(file = file, width = (as.numeric(input$figWidthORA)/50), height = (as.numeric(input$figHeightORA)/50), )
            print(oraPlotList()[["cnetplot"]])
            dev.off()
          }
        )

        # bar plot
        output$barPlot_ORA <- renderPlot(
          {
            oraPlotList()[["barplot"]]
          },
          width = as.numeric(input$figWidthORA),
          height = as.numeric(input$figHeightORA)
        )
        output$dlBarPlot_ORA <- downloadHandler(
          filename = function() {
            paste0(design, "ORA_", input$oraType, "_bar.pdf")
          },
          content = function(file) {
            pdf(file = file, width = (as.numeric(input$figWidthORA)/50), height = (as.numeric(input$figHeightORA)/50), )
            print(oraPlotList()[["barplot"]])
            dev.off()
          }
        )

        # Dot plot
        output$dotPlot_ORA <- renderPlot(
          {
            oraPlotList()[["dotplot"]]
          },
          width = as.numeric(input$figWidthORA),
          height = as.numeric(input$figHeightORA)
        )
        output$dlDotPlot_ORA <- downloadHandler(
          filename = function() {
            paste0(design, "ORA_", input$oraType, "_dot.pdf")
          },
          content = function(file) {
            pdf(file = file, width = (as.numeric(input$figWidthORA)/50), height = (as.numeric(input$figHeightORA)/50), )
            print(oraPlotList()[["dotplot"]])
            dev.off()
          }
        )

        # Heatmap
        output$heatmap_ORA <- renderPlot(
          {
            oraPlotList()[["heatmap"]]
          },
          width = as.numeric(input$figWidthORA*2),
          height = as.numeric(input$figHeightORA)
        )
        output$dlHeatmap_ORA <- downloadHandler(
          filename = function() {
            paste0(design, "ORA_", input$oraType, "_heatmap.pdf")
          },
          content = function(file) {
            pdf(file = file, width = (as.numeric(input$figWidthORA)/25), height = (as.numeric(input$figHeightORA)/50), )
            print(oraPlotList()[["heatmap"]])
            dev.off()
          }
        )

      }
    )
  }
} else {
  observeEvent(input$tabs, {
    if (input$tabs == "oraTab") {
      showModal(modalDialog(title = "No ORA Results", "Unfortunately, our proteomics results don't include any pathway analyses! We hope to change this very soon."))
    }
  })
}