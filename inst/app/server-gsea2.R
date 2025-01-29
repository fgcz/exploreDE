if(inputDataReactive()$dataType == "RNASeq") {

  param <- inputDataReactive()$param
  
  observeEvent(input$tabs, {
    if (input$tabs == "gseaTab") {
      if (is_empty(allPathways)) {
        showModal(modalDialog(title = "No GSEA Results", "No pathway analysis was performed so this tab is empty."))
      }
    }
  })
  
  if(isTRUE(param$runGO & param$featureLevel != "isoform")) {
    gsea <- inputDataReactive()$gsea
    gseaHTML <- inputDataReactive()$gseaHTML
    se <- inputDataReactive()$se
    design <- inputDataReactive()$design
  
    req(!is.null(inputDataReactive()$gsea))
    
    output$gseaDesign <- renderText({design})
    
    # Some older datasets seem to be missing this parameter, so we will just assume it for now
    if (is.null(param$fdrThreshGSEA)) {
      param$fdrThreshGSEA <- 0.05
    }
    
    # Download full GSEA results file
    output$dlTableGSEA <- downloadHandler(
      filename = function() {
        basename(gseaHTML)
      },
      content = function(file) {
        file.copy(from = gseaHTML, to = file)
      }
    )
    # GSEA Overview Table
    gsea.overview <- ezFrame(
      "Number of Significant Pathways" = numeric()
    )
    for (i in 1:length(names(gsea))) {
      testName <- paste(names(gsea)[i])
      gsea.overview[testName, "Number of Significant Pathways"] <-
        nrow(gsea[[i]]@result)
    }
    output$gseaOverview <- function() {
      gsea.overview %>%
        kable(
          format = "html",
          caption =
            paste("GSEA:")
        ) %>%
        kable_styling(
          bootstrap_options = "striped",
          full_width = FALSE,
          position = "left"
        )
    }
    # ORA thresholds
    gseaThresholds <- ezFrame(
      "GSEA Cut-Offs" = character()
    )
    # gseaThresholds["Input DEG P-value", "GSEA Cut-Offs"] <- param$pValThreshGO
    # gseaThresholds["Input DEG Log2 Ratio", "GSEA Cut-Offs"] <- param$log2RatioThreshGO
    gseaThresholds["Output GSEA Term FDR", "GSEA Cut-Offs"] <- param$fdrThreshGSEA
    output$gseaThreshOverview <- function() {
      gseaThresholds %>%
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
    
    observeEvent(input$gseaType, {
      er <- gsea[[input$gseaType]]
      if (input$tabs == "gseaTab") {
        if (is.null(er) || nrow(er@result) == 0) {
          showModal(modalDialog(title = "No GSEA Results", "There are no results for this selection."))
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
        output$selectedTable_GSEA <- DT::renderDataTable({
          DT::datatable(
            data = as.data.frame(er)[, c("Label", "enrichmentScore", "NES", "geneName")],
            filter = "top", class = "cell-border stripe",
            rownames = FALSE, caption = "Click pathways in this table to add them to the figures on the right"
          ) %>%
            DT::formatStyle(
              columns = colnames(.$x$data), `font-size` = "12px"
            ) %>%
            DT::formatSignif(
              columns = c("enrichmentScore", "NES"),
              digits = 3
            ) %>%
            DT::formatStyle(
              columns = c("enrichmentScore", "NES"),
              color = styleInterval(cuts = 0, values = c("blue", "darkorange")),
              fontWeight = "bold"
            )
        })
      }
    })
    
    gseaPlotList <- eventReactive({
      input$textSizeGSEA
      input$nodeSizeGSEA
      input$selectedTable_GSEA_rows_selected
      input$scaleLimGSEA
      input$nodeSizeGSEA
      input$gseaMaxOver
      input$gseaDodgeY
      input$gseaDodgeX
      input$nodeBorderGSEA
      input$ridgePlotColourByGSEA
      input$gseaLabelAlpha
    }, {
      
      req(!is.null(input$selectedTable_GSEA_rows_selected))
      
      er <- gsea[[input$gseaType]]
      if (!is.null(er) && nrow(er@result) >= 2) {
        er@result$Label <- ""
        for (i in seq_along(er@result$Description)) {
          desc <- er@result$Description[i]
          label <- paste0(er@result$ID[i], "\n", desc)
          if (!is.na(label) && nchar(label) > 45) {
            label <- substr(label, 1, 45)
          }
          er@result$Label[i] <- label
        }
        
        erRP <- clusterProfiler::slice(er, input$selectedTable_GSEA_rows_selected)
        rp <- enrichplot::ridgeplot(x = erRP, fill = input$ridgePlotColourByGSEA)
        rp <- rp +
          scale_y_discrete(labels = (erRP@result$Label[match(levels(as.factor(rp$plot_env$gs2val.df$category)), erRP@result$Description)])) +
          theme(axis.text.y = element_text(vjust = -0.01, size = input$textSizeGSEA)) +
          geom_vline(xintercept = 0, linetype = "dashed") +
          scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4", mid = "white", midpoint = 0)
        ce_list <- NA
        for (i in 1:length(er@result$core_enrichment)) {
          core_enrich <- strsplit(
            er@result$core_enrichment[i],
            split = "/"
          ) %>% unlist()
          ce_list <- c(ce_list, core_enrich)
          core_enrich_symbol <- metadata(se)$enrichInput$seqAnno$gene_name[
            metadata(se)$enrichInput$seqAnno$gene_id %in% core_enrich
          ]
          core_enrich_symbol <- paste(core_enrich_symbol, collapse = "/")
          er@result$core_enrichment[i] <- core_enrich_symbol
        }
        sigGeneIds <- ce_list[!is.na(ce_list)]
        
        log2RatioGSEA <- metadata(se)$enrichInput$seqAnno$log2Ratio
        names(log2RatioGSEA) <- metadata(se)$enrichInput$seqAnno$gene_name
        log2RatioGSEA <- sort(log2RatioGSEA, decreasing = T)
        
        # Draw the main plot 
        cn <- enrichplot::cnetplot(
          x = er,
          node_label = "none",
          foldChange = log2RatioGSEA,
          showCategory = er@result$Description[input$selectedTable_GSEA_rows_selected],
          cex.params = list(category_label = (input$textSizeGSEA / 15), gene_label = (input$textSizeGSEA / 18), gene_node = input$nodeSizeGSEA/10, category_node = input$nodeSizeGSEA/8)
        )
        # Set the colour limits to the user-specified values 
        cn$data$foldChange[!is.na(cn$data$foldChange) & cn$data$color > input$scaleLimGSEA] <- input$scaleLimGSEA
        cn$data$foldChange[!is.na(cn$data$foldChange) & cn$data$foldChange < -input$scaleLimGSEA] <- -input$scaleLimGSEA
        # Add our own geom_points on top 
        cn <- cn + geom_point(data = cn$data, aes(x = x, y = y, size = size), shape = 21, stroke = input$nodeBorderGSEA)
        # Get two copies of the plot data table so we can get gene and node labels and modify them for plotting 
        cd1 <- cn$data
        cd1$name[!is.na(cd1$foldChange)] <- NA
        cd1$Label <- NA
        cd1$Label[1:length(input$selectedTable_GSEA_rows_selected)] <- er@result$Label[input$selectedTable_GSEA_rows_selected]
        cd2 <- cn$data
        cd2$name[is.na(cd2$foldChange)] <- NA
        # Some reason, enrichplot changed the default such that up is blue and down is red. We prefer the opposite... 
        cn <- cn + scale_colour_gradient2(name = "fold change", low = "deepskyblue4", mid = "white", high = "firebrick", midpoint = 0)
        # Add gene and node labels in that order 
        if (input$showGeneLabelsGSEA) {
          cn <- cn +
            geom_text_repel(
              data = cd2, aes(label = name, x = x, y = y), size = input$textSizeGSEA/3, max.overlaps = input$gseaMaxOver, nudge_x = input$gseaDodgeX/10, nudge_y = input$gseaDodgeY/10,
              fontface = "bold", color = "black", bg.color = "white", bg.r = .15, na.rm = T)
        }
        cn <- cn + geom_label_repel(
          data = cd1, aes(label = Label, x = x, y = y, fill = Label), size = input$textSizeGSEA/3, max.overlaps = input$gseaMaxOver, nudge_x = input$gseaDodgeX/10, nudge_y = input$gseaDodgeY/10, fontface = "bold", show.legend = F, fill = alpha(c("white"), input$gseaLabelAlpha), na.rm = T)
        
        hp <- enrichplot::heatplot(
          x = er,
          foldChange = log2RatioGSEA,
          showCategory = er@result$Description[input$selectedTable_GSEA_rows_selected],
        ) + scale_colour_distiller(palette = "RdBu")
        rs <- enrichplot::gseaplot2(
          x = er,
          geneSetID = input$selectedTable_GSEA_rows_selected,
          color = RColorBrewer::brewer.pal(length(input$selectedTable_GSEA_rows_selected), "Set2")
        )
        
        # make an upset plot of the genes in the selected pathways 
        geneSetsUpset <- setNames(lapply(input$selectedTable_GSEA_rows_selected, function(i) {
          g <- er@result$geneName[i] %>% strsplit("/") %>% .[[1]]
          g <- g[nchar(g) > 1]
          g
        }), er@result$Label[input$selectedTable_GSEA_rows_selected])
        m1 <- ComplexHeatmap::make_comb_mat(geneSetsUpset)
        geneLFCsUpset <- setNames(lapply(names(comb_size(m1)), function(x) {
          g <- extract_comb(m1, x)
          inputDataReactive()$seqAnno$log2Ratio[which(inputDataReactive()$seqAnno$gene_name %in% g)]
        }), names(comb_size(m1)))
        up <- ComplexHeatmap::UpSet(
          m = m1, 
          comb_order = order(comb_size(m1), decreasing = T), 
          top_annotation = upset_top_annotation(m1, add_numbers = TRUE, annotation_name_rot = 90, annotation_name_offset = unit(10, "mm")), 
          left_annotation = upset_left_annotation(m1, add_numbers = TRUE)
          # bottom_annotation = HeatmapAnnotation(
          #   "Log2 Ratios" = anno_boxplot(geneLFCsUpset),
          #   annotation_name_side = "right"
          # )
        )
        
        return(list(
          "cnetplot" = cn,
          "heatmap" = hp,
          "runningscore" = rs,
          "ridgeplot" = rp,
          "upset" = up
        ))
      }
    })
    
    # Custom GSEA results output
    observeEvent(
      {
        input$figHeightGSEA
        input$figWidthGSEA
        gseaPlotList()
      }, {
       
        # Network plot 
        output$cnetPlot_GSEA <- renderPlot(
          {
            gseaPlotList()[["cnetplot"]]
          },
          width = as.numeric(input$figWidthGSEA),
          height = as.numeric(input$figHeightGSEA)
        )
        output$dlCnetPlot_GSEA <- downloadHandler(
          filename = function() {
            paste0(design, "GSEA_", input$gseaType, "_network.pdf")
          },
          content = function(file) {
            pdf(file = file, width = (as.numeric(input$figWidthGSEA)/50), height = (as.numeric(input$figHeightGSEA)/50))
            print(gseaPlotList()[["cnetplot"]])
            dev.off()
          }
        )
        
        # Heatmap
        output$heatmap_GSEA <- renderPlot(
          {
            gseaPlotList()[["heatmap"]]
          },
          width = as.numeric(input$figWidthGSEA*2),
          height = as.numeric(input$figHeightGSEA)
        )
        output$dlHeatmap_GSEA <- downloadHandler(
          filename = function() {
            paste0(design, "GSEA_", input$gseaType, "_heatmap.pdf")
          },
          content = function(file) {
            pdf(file = file, width = (as.numeric(input$figWidthGSEA)/25), height = (as.numeric(input$figHeightGSEA)/50))
            print(gseaPlotList()[["heatmap"]])
            dev.off()
          }
        )
        
        # Running score
        output$runningScore_GSEA <- renderPlot(
          {
            gseaPlotList()[["runningscore"]]
          },
          width = as.numeric(input$figWidthGSEA),
          height = as.numeric(input$figHeightGSEA)
        )
        output$dlRunningScore_GSEA <- downloadHandler(
          filename = function() {
            paste0(design, "GSEA_", input$gseaType, "_runningscore.pdf")
          },
          content = function(file) {
            pdf(file = file, width = (as.numeric(input$figWidthGSEA)/50), height = (as.numeric(input$figHeightGSEA)/50))
            print(gseaPlotList()[["runningscore"]])
            dev.off()
          }
        )
        
        # Ridge plot 
        output$ridgePlot_GSEA <- renderPlot(
          {
            gseaPlotList()[["ridgeplot"]]
          },
          width = as.numeric(input$figWidthGSEA),
          height = as.numeric(input$figHeightGSEA)
        )
        output$dlRidgePlot_GSEA <- downloadHandler(
          filename = function() {
            paste0(design, "GSEA_", input$gseaType, "_ridgeplot.pdf")
          },
          content = function(file) {
            pdf(file = file, width = (as.numeric(input$figWidthGSEA)/50), height = (as.numeric(input$figHeightGSEA)/50))
            print(gseaPlotList()[["ridgeplot"]])
            dev.off()
          }
        )
        
        # Upset plot 
        output$upsetPlot_GSEA <- renderPlot(
          {
            draw(gseaPlotList()[["upset"]], padding = unit(c(5, 5, 5, 20), "mm"))
          },
          width = as.numeric(input$figWidthGSEA)*1.3,
          height = as.numeric(input$figHeightGSEA)
        )
        output$dlUpset_GSEA <- downloadHandler(
          filename = function() {
            paste0(design, "GSEA_", input$gseaType, "_upset.pdf")
          },
          content = function(file) {
            pdf(file = file, width = (as.numeric(input$figWidthGSEA)/70), height = (as.numeric(input$figHeightGSEA)/90))
            print(gseaPlotList()[["upset"]])
            dev.off()
          }
        )
        
      }
    )
  }
} else {
  observeEvent(input$tabs, {
    if (input$tabs == "gseaTab") {
      showModal(modalDialog(title = "No GSEA Results", "Unfortunately, our proteomics results don't include any pathway analyses! We hope to change this very soon."))
    }
  })
}