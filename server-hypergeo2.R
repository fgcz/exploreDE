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
    
    # Some older datasets seem to be missing this parameter, so we will just assume it for now
    if (is.null(param$fdrThreshORA)) {
      param$fdrThreshORA <- 0.05
    }
    
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
      input$scaleLimORA
      input$oraMaxOver
      input$oraDodgeY
      input$oraDodgeX
      input$nodeBorderORA
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
          # Draw the main plot 
          cn <- cnetplot(
            x = er,
            foldChange = log2RatioORA,
            showCategory = er@result$Description[input$selectedTable_ORA_rows_selected],
            node_label = "none",
            cex.params = list(category_label = (input$textSizeORA / 15), gene_label = (input$textSizeORA / 18), gene_node = input$nodeSizeORA/10, category_node = input$nodeSizeORA/8)
          )
          # Set the colour limits to the user-specified values 
          cn$data$foldChange[!is.na(cn$data$foldChange) & cn$data$foldChange > input$scaleLimORA] <- input$scaleLimORA
          cn$data$foldChange[!is.na(cn$data$foldChange) & cn$data$foldChange < -input$scaleLimORA] <- -input$scaleLimORA
          # Add our own geom_points on top 
          cn <- cn + geom_point(data = cn$data, aes(x = x, y = y, size = size), shape = 21, stroke = input$nodeBorderORA)
          # Get two copies of the plot data table so we can get gene and node labels and modify them for plotting 
          cd1 <- cn$data
          cd1$name[!is.na(cd1$foldChange)] <- NA
          cd1$Label <- NA
          # cd1$Label[1:length(input$selectedTable_ORA_rows_selected)] <- er@result$Label[input$selectedTable_ORA_rows_selected]
          cd2 <- cn$data
          cd2$name[is.na(cd2$foldChange)] <- NA
          # Some reason, enrichplot changed the default such that up is blue and down is red. We prefer the opposite... 
          if (input$oraDirection == "upGenes") {
            cn <- cn + scale_color_gradientn(name = "fold change", colours = c("white", "firebrick"))
          } else if (input$oraDirection == "downGenes") {
            cn <- cn + scale_color_gradientn(name = "fold change", colours = c("deepskyblue4", "white"))
          } else if (input$oraDirection == "bothGenes") {
            cn <- cn + scale_colour_gradient2(name = "fold change", low = "deepskyblue4", mid = "white", high = "firebrick", midpoint = 0)
          }
          # Add gene and node labels in that order 
          if (input$showGeneLabelsORA) {
            cn <- cn +
              geom_text_repel(
                data = cd2, aes(label = name, x = x, y = y), size = input$textSizeORA/3, max.overlaps = input$oraMaxOver, nudge_x = input$oraDodgeX/10, nudge_y = input$oraDodgeY/10,
                fontface = "bold", color = "black", bg.color = "white", bg.r = .15, na.rm = T)
          }
          cn <- cn + geom_label_repel(
            data = cd1, aes(label = Label, x = x, y = y, fill = Label), size = input$textSizeORA/3, max.overlaps = input$oraMaxOver, nudge_x = input$oraDodgeX/10, nudge_y = input$oraDodgeY/10, fontface = "bold", show.legend = F, fill = alpha(c("white"),0.5), na.rm = T)
          
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
            theme_prism() +
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
            theme_prism() +
            theme(
              axis.text.y = element_text(
                vjust = -0.01, size = input$textSizeORA))
        }
        
        # make an upset plot of the genes in the selected pathways 
        geneSetsUpset <- setNames(lapply(input$selectedTable_ORA_rows_selected, function(i) {
          g <- er@result$geneName[i] %>% strsplit("/") %>% .[[1]]
          g <- g[nchar(g) > 1]
          g
        }), er@result$Label[input$selectedTable_ORA_rows_selected])
        m1 <- ComplexHeatmap::make_comb_mat(geneSetsUpset)
        up <- ComplexHeatmap::UpSet(
          m = m1, 
          comb_order = order(comb_size(m1), decreasing = T), 
          top_annotation = upset_top_annotation(m1, add_numbers = TRUE, annotation_name_rot = 90, annotation_name_offset = unit(10, "mm")), 
          left_annotation = upset_left_annotation(m1, add_numbers = TRUE)
          )
        
        return(list(
          "er" = er,
          "cnetplot" = cn,
          "barplot" = bp,
          "dotplot" = dp,
          "heatmap" = hp,
          "upset" = up
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
            pdf(file = file, width = (as.numeric(input$figWidthORA)/80), height = (as.numeric(input$figHeightORA)/80))
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
        
        # upset 
        output$upset_ORA <- renderPlot(
          {
            req(!is.null(oraPlotList()[["upset"]]))
            draw(oraPlotList()[["upset"]], padding = unit(c(5, 5, 5, 20), "mm"))
          },
          width = as.numeric(input$figWidthORA)*1.3,
          height = as.numeric(input$figHeightORA)
        )
        output$dlUpset_ORA <- downloadHandler(
          filename = function() {
            paste0(design, "ORA_", input$oraType, "_upset.pdf")
          },
          content = function(file) {
            pdf(file = file, width = (as.numeric(input$figWidthORA)/70), height = (as.numeric(input$figHeightORA)/90))
            print(oraPlotList()[["upset"]])
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