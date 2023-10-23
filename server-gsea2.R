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
  
    req(inputDataReactive()$gsea)
    
    output$gseaDesign <- renderText({design})
    
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
            data = as.data.frame(er)[, c("Label", "enrichmentScore", "geneName")],
            filter = "top", class = "cell-border stripe",
            rownames = FALSE, caption = "Click pathways in this table to add them to the figures on the right"
          ) %>%
            DT::formatStyle(
              columns = colnames(.$x$data), `font-size` = "12px"
            ) %>%
            DT::formatSignif(
              columns = c("enrichmentScore"),
              digits = 3
            ) %>%
            DT::formatStyle(
              columns = c("enrichmentScore"),
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
    }, {
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
        rp <- enrichplot::ridgeplot(x = erRP)
        rp <- rp +
          scale_y_discrete(labels = (erRP@result$Label[match(levels(as.factor(rp$plot_env$gs2val.df$category)), erRP@result$Description)])) +
          theme(axis.text.y = element_text(vjust = -0.01, size = input$textSizeGSEA)) +
          geom_vline(xintercept = 0, linetype = "dashed")
  
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
        log2RatioGSEA <- metadata(se)$enrichInput$log2Ratio[sigGeneIds]
        names(log2RatioGSEA) <- metadata(se)$enrichInput$seqAnno[
          sigGeneIds, "gene_name"
        ]
        log2RatioGSEA <- sort(log2RatioGSEA, decreasing = T)
        
        cn <- enrichplot::cnetplot(
          x = er,
          color.params = list(foldChange = log2RatioGSEA),
          showCategory = er@result$Description[input$selectedTable_GSEA_rows_selected],
          cex.params = list(category_label = (input$textSizeGSEA / 15), gene_label = (input$textSizeGSEA / 18), gene_node = (input$textSizeGSEA / 16), category_node = (input$textSizeGSEA / 13))
        )
        cn
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
        
        return(list(
          "cnetplot" = cn,
          "heatmap" = hp,
          "runningscore" = rs,
          "ridgeplot" = rp
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
            pdf(file = file, width = (as.numeric(input$figWidthGSEA)/50), height = (as.numeric(input$figHeightGSEA)/50), )
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
            pdf(file = file, width = (as.numeric(input$figWidthGSEA)/25), height = (as.numeric(input$figHeightGSEA)/50), )
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
            pdf(file = file, width = (as.numeric(input$figWidthGSEA)/50), height = (as.numeric(input$figHeightGSEA)/50), )
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
            pdf(file = file, width = (as.numeric(input$figWidthGSEA)/50), height = (as.numeric(input$figHeightGSEA)/50), )
            print(gseaPlotList()[["ridgeplot"]])
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