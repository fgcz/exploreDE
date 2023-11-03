req(!is.null(inputDataReactive()$dataType))

if (inputDataReactive()$dataType == "proteomics") {
  updateNumericInput(session = session, inputId = "xLimVolcano", value = 2)
}

volcanoColoursDefault <- c(
  "Highlight" = "goldenrod3",
  "NotSignificant" = "grey30",
  "SignificantUp" = "pink",
  "SignificantDown" = "lightblue",
  "FoldChangeUp" = "pink3",
  "FoldChangeDown" = "lightblue4",
  "FoldChange&SignificantUp" = "firebrick4",
  "FoldChange&SignificantDown" = "dodgerblue4"
)

output$volcanoColourPicker <- renderUI({
  lapply(seq_along(volcanoColoursDefault), function(i) {
    colourpicker::colourInput(
      inputId = paste0("volcanoColour", i),
      label = names(volcanoColoursDefault)[i],
      value = as.character(volcanoColoursDefault[i]),
      closeOnClick = TRUE,
      returnName = TRUE
    )
  })
})

updateSelectizeInput(
  session = session,
  inputId = "volcanoGenes",
  choices = inputDataReactive()$genes,
  selected = "",
  server = TRUE
)

observeEvent({
  genesReactive()
}, {
  output$geneBucket2 <- renderUI({
    bucket_list(
      header = "Drag and drop features in order to be plotted",
      group_name = "bucket_list_group",
      orientation = "horizontal",
      add_rank_list(
        text = "Include these features in this order",
        labels = genesReactive()$genes,
        input_id = "boxKeepBucketGenes"),
      add_rank_list(
        text = "Exclude these features",
        labels = NULL,
        input_id = "boxExcludeBucketGenes")
    )
  })
})

observeEvent(
  {
    input$pTypeVolcano
    input$contrastSelected
    inputDataReactive()
  }, ignoreInit = FALSE, ignoreNULL = TRUE, 
  {
    if (inputDataReactive()$dataType == "RNASeq") {
      seqAnno <- inputDataReactive()$seqAnno
      design <- inputDataReactive()$design
      output$volcanoDesign <- renderText({
        design
      })
      seqAnnoFilt <- seqAnno[which(seqAnno$usedInTest),]
    }
    
    if (inputDataReactive()$dataType == "proteomics") {
      req(!is.null(input$contrastSelected))
      seqAnnoFilt <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]
      seqAnnoFilt <- seqAnnoFilt %>% dplyr::select(gene_name, log2Ratio, pValue, fdr)
      
      output$volcanoDesign <- renderText({
        input$contrastSelected
      })
    }
    
    if (input$pTypeVolcano == "Raw") {
      pTypeVolcano <- "pValue"
    } else {
      pTypeVolcano <- "fdr"
    }
    
    updateNumericInput(
      session = session,
      inputId = "yLimVolcano",
      value = ceiling(max(-log10(seqAnnoFilt[[pTypeVolcano]]))*1.1)
    )
    updateNumericInput(
      session = session,
      inputId = "xLimVolcano",
      value = ceiling(max(seqAnnoFilt[["log2Ratio"]])*1.1)
    )
  }
)

# Volcano Plot Output ----
volcanoResultsList <- eventReactive(
  {
    input$lfcVolcano
    input$yLimVolcano
    input$xLimVolcano
    input$pTypeVolcano
    input$pThresholdVolcano
    input$volcanoShowGenes
    input$boxKeepBucketGenes
    input$contrastSelected
    # genesReactive()
    lapply(1:8, function(i) {
      input[[paste0("volcanoColour", i)]]
    })
  },
  ignoreInit = TRUE, ignoreNULL = FALSE, 
  {
    tryCatch(
      {
        if (inputDataReactive()$dataType == "RNASeq") {
          seqAnno <- inputDataReactive()$seqAnno
          design <- inputDataReactive()$design
          output$volcanoDesign <- renderText({
            design
          })
          seqAnnoFilt <- seqAnno[which(seqAnno$usedInTest),]
        }
        
        if (inputDataReactive()$dataType == "proteomics") {
          seqAnnoFilt <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]
          seqAnnoFilt <- seqAnnoFilt %>% dplyr::select(gene_name, description, log2Ratio, pValue, fdr)
          
          output$volcanoDesign <- renderText({
            input$contrastSelected
          })
        }
        
        if (input$pTypeVolcano == "Raw") {
          pTypeVolcano <- "pValue"
        } else {
          pTypeVolcano <- "fdr"
        }
        
        # Get the colours
        volcanoColours <- c(
          "Highlighted" = 1,
          "NotSignificant" = 1,
          "SignificantUp" = 1,
          "SignificantDown" = 1,
          "FoldChangeUp" = 1,
          "FoldChangeDown" = 1,
          "FoldChange&SignificantUp" = 1,
          "FoldChange&SignificantDown" = 1
        )
        for (i in 1:8) {
          volcanoColours[[i]] <- input[[paste0("volcanoColour", i)]]
        }
        
        # Create table of gene names, LFC, and p-value selected by user:
        volcanoTable <- seqAnnoFilt[, c("gene_name", "description", "log2Ratio", pTypeVolcano)]
        
        # Get the -log10 of the p-value:
        volcanoTable["log10p"] <- -log10(as.numeric(volcanoTable[[pTypeVolcano]]))
        
        # Get the gene labels for plotting
        volcanoTable$genelabels <- ""
        volcanoTable$genelabels[which(
          volcanoTable$gene_name %in% input$boxKeepBucketGenes
        )] <- volcanoTable$gene_name[which(
          volcanoTable$gene_name %in% input$boxKeepBucketGenes
        )]
        
        # Get the -log10p of p-value threshold as a single value:
        log10p <- -log10(as.numeric(input$pThresholdVolcano))
        
        # Get mean log2 count for MA plot
        if (inputDataReactive()$dataType == "RNASeq") {
          rMs <- rowMeans(inputDataReactive()$countList$`Normalised + Log2`)
        } else {
          rMs <- rowMeans(inputDataReactive()$countList[[2]])
          names(rMs) <- names(rMs) %>% gsub("\\~.*", "", .)
        }
        volcanoTable$Log2_Mean <- rMs[match(volcanoTable$gene_name, names(rMs))]
        
        # Determine what significance Status each gene is and label accordingly
        ## Not significant:
        volcanoTable["Status"] <- "NotSignificant"
        ## Significant up only:
        volcanoTable[which(volcanoTable["log10p"] >= log10p & volcanoTable["log2Ratio"] <= as.numeric(input$lfcVolcano) & volcanoTable["log2Ratio"] > 0), "Status"] <- "SignificantUp"
        ## Significant down only:
        volcanoTable[which(volcanoTable["log10p"] >= log10p & volcanoTable["log2Ratio"] < 0 & volcanoTable["log2Ratio"] >= -as.numeric(input$lfcVolcano)), "Status"] <- "SignificantDown"
        ## LFC up only:
        volcanoTable[which(volcanoTable["log10p"] < log10p & volcanoTable["log2Ratio"] >= as.numeric(input$lfcVolcano)), "Status"] <- "FoldChangeUp"
        ## LFC down only:
        volcanoTable[which(volcanoTable["log10p"] < log10p & volcanoTable["log2Ratio"] <= -as.numeric(input$lfcVolcano)), "Status"] <- "FoldChangeDown"
        # LFC & p-value up:
        volcanoTable[which(volcanoTable["log10p"] >= log10p & volcanoTable["log2Ratio"] >= as.numeric(input$lfcVolcano)), "Status"] <- "FoldChange&SignificantUp"
        # LFC & p-value down:
        volcanoTable[which(volcanoTable["log10p"] >= log10p & volcanoTable["log2Ratio"] <= -as.numeric(input$lfcVolcano)), "Status"] <- "FoldChange&SignificantDown"
        # Highlighted
        vc2 <- volcanoTable
        vc2 <- vc2[order(vc2$log10p, decreasing = TRUE), c("gene_name", "description", "log2Ratio", pTypeVolcano, "log10p", "Status")]
        
        # Set limits of data frame to x/y-axes limits
        volcanoTable$log10p_full <- volcanoTable$log10p
        volcanoTable$log2Ratio_full <- volcanoTable$log2Ratio
        volcanoTable$log10p[which(volcanoTable$log10p == "Inf")] <- max(volcanoTable$log10p[which(volcanoTable$log10p < Inf)])
        volcanoTable$log10p[volcanoTable$log10p > as.numeric(input$yLimVolcano)] <- as.numeric(input$yLimVolcano)
        volcanoTable$log2Ratio[volcanoTable$log2Ratio > as.numeric(input$xLimVolcano)] <- as.numeric(input$xLimVolcano)
        volcanoTable$log2Ratio[volcanoTable$log2Ratio < -as.numeric(input$xLimVolcano)] <- -as.numeric(input$xLimVolcano)
        
        # Summary table
        output$volcanoOverview <- function() {
          table(as.factor(volcanoTable$Status)) %>%
            kable(
              format = "html",
              caption = paste("DE Summary with Volcano Settings")
            ) %>%
            kable_styling(
              bootstrap_options = "striped",
              full_width = FALSE,
              position = "left"
            )
        }
        
        # Results table
        output$volcanoOverviewTable <- DT::renderDataTable({
          DT::datatable(
            data = vc2,
            filter = "top",
            class = "cell-border stripe",
            rownames = FALSE,
            colnames = c("Gene Symbol", "Description", "Log2 Ratio", pTypeVolcano, paste("-log10", pTypeVolcano), "Significance")
          ) %>%
            DT::formatSignif(columns = c("log10p", pTypeVolcano, "log2Ratio"), digits = 3) %>%
            DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
            DT::formatStyle(columns = "log2Ratio", color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold") %>%
            DT::formatStyle(columns = c("log10p", pTypeVolcano), color = styleInterval(cuts = log10p, values = c("black", "green")), fontWeight = "bold") %>%
            DT::formatStyle(columns = "Status", color = styleEqual(levels = names(volcanoColours), values = as.character(volcanoColours)), fontWeight = "bold")
        })
        
        return(list(
          "volcanoTable" = volcanoTable,
          "volcanoColours" = volcanoColours
        ))
      },
      error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
  }
)

observeEvent(
  {
    input$figWidthVolcano
    input$figHeightVolcano
    input$textSizeVolcano
    input$showLinesVolcano
    input$showAxesVolcano
    input$boldVolcano
    input$alphaVolcano
    input$dotSizeVolcano
    input$volcanoShowGenes
    input$boxKeepBucketGenes
    input$showBorderVolcano
    volcanoResultsList()
  },
  {
    volcanoTable <- volcanoResultsList()$volcanoTable
    volcanoColours <- volcanoResultsList()$volcanoColours

    # Get a df of the highlighted genes so we can make them bigger
    volcanoTableFull <- volcanoTable
    if (input$volcanoShowGenes == TRUE & length(input$boxKeepBucketGenes) >= 1) {
      volcanoHighlights <- volcanoTable[which(volcanoTable$gene_name %in% input$boxKeepBucketGenes), ]
      volcanoHighlights$Status <- "Highlighted"
      volcanoTable <- volcanoTable[which(!volcanoTable$gene_name %in% input$boxKeepBucketGenes), ]
    }

    # Create the MA plot
    if (input$volcanoShowGenes == TRUE & length(input$boxKeepBucketGenes) >= 1) {
      maplot <- ggplot(volcanoTable, aes(x = Log2_Mean, y = log2Ratio_full))
      if (input$showBorderVolcano) {
        maplot <- maplot + geom_point(data = volcanoHighlights, aes(y = log2Ratio_full, x = Log2_Mean, fill = Status), size = as.numeric(input$dotSizeVolcano), alpha = 0.9, pch = 21)
        maplot <- maplot + geom_point(pch = 21, size = as.numeric(input$dotSizeVolcano), alpha = as.numeric(input$alphaVolcano), aes(fill = Status))
        maplot <- maplot + scale_fill_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
      } else {
        maplot <- maplot + geom_point(data = volcanoHighlights, aes(y = log2Ratio_full, x = Log2_Mean, colour = Status), size = as.numeric(input$dotSizeVolcano), alpha = 0.9, pch = 16)
        maplot <- maplot + geom_point(pch = 16, size = as.numeric(input$dotSizeVolcano), alpha = as.numeric(input$alphaVolcano), aes(colour = Status))
        maplot <- maplot + scale_colour_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
      }
      maplot <- maplot +
        theme_prism(base_size = input$textSizeVolcano) +
        xlab("Log2 Normalised Mean") + ylab("Log2 Ratio") +
        geom_label_repel(
          data = volcanoHighlights,
          aes(label = genelabels),
          force = 2,
          nudge_y = 0.2,
          nudge_x = 0.2,
          color = "black",
          size = (input$textSizeVolcano / 3),
          max.overlaps = 500,
          fontface = "bold"
        )
    } else {
      maplot <- ggplot(volcanoTable, aes(x = Log2_Mean, y = log2Ratio_full))
      if (input$showBorderVolcano) {
        maplot <- maplot + geom_point(pch = 21, size = as.numeric(input$dotSizeVolcano), alpha = as.numeric(input$alphaVolcano), aes(fill = Status))
        maplot <- maplot + scale_fill_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
      } else {
        maplot <- maplot + geom_point(pch = 16, size = as.numeric(input$dotSizeVolcano), alpha = as.numeric(input$alphaVolcano), aes(colour = Status))
        maplot <- maplot + scale_colour_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
      }
      maplot <- maplot +
        theme_prism(base_size = input$textSizeVolcano) +
        xlab("Log2 Normalised Mean") + ylab("Log2 Ratio")

    }
    output$MAPlot <- renderPlot({
      maplot
    }, height = as.numeric(input$figHeightVolcano)/1.5, width = as.numeric(input$figWidthVolcano)*1.2)
    output$dlMAPlotButton <- downloadHandler(
      filename = function() {
        paste0(design, "_MA.pdf")
      },
      content = function(file) {
        ggsave(
          filename = file,
          plot = maplot,
          width = (as.numeric(input$figWidthVolcano) / 102),
          height = (as.numeric(input$figHeightVolcano)/ 115),
          limitsize = FALSE, units = "in"
        )
      }
    )
    output$MABrushTable <- renderDataTable({
      volcanoTableFullMA <- volcanoTableFull %>% dplyr::select(gene_name, description, log2Ratio_full, Log2_Mean, log10p_full, Status)
      DT::datatable(
        data = brushedPoints(volcanoTableFullMA, input$MABrush),
        colnames = c("Gene name", "Description", "Log2 Ratio", "Log2 Mean", "-Log10p", "Status"),
        filter = "top",
        caption = "Click and drag over dots on the volcano plot to see those genes in this table.",
        rownames = F) %>%
        DT::formatSignif(c("log2Ratio_full", "log10p_full", "Log2_Mean"), digits = 3) %>%
        DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
        DT::formatStyle(columns = c("log2Ratio_full"), color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold") %>%
        DT::formatStyle(columns = c("log10p_full"), color = styleInterval(cuts = -log10(as.numeric(input$pThresholdVolcano)), values = c("black", "green")), fontWeight = "bold") %>%
        DT::formatStyle(columns = "Status", color = styleEqual(levels = names(volcanoColours), values = as.character(volcanoColours)), fontWeight = "bold")
    })

    # Create the volcano plot
    volcanoStatic <- ggplot(
      data = volcanoTable, aes(x = log2Ratio, y = log10p)
    )
    if (input$showBorderVolcano) {
      volcanoStatic <- volcanoStatic + geom_point(aes(fill = Status), size = as.numeric(input$dotSizeVolcano), alpha = as.numeric(input$alphaVolcano), pch = 21)
      volcanoStatic <- volcanoStatic + scale_fill_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
    } else {
      volcanoStatic <- volcanoStatic + geom_point(aes(colour = Status), size = as.numeric(input$dotSizeVolcano), alpha = as.numeric(input$alphaVolcano), pch = 16)
      volcanoStatic <- volcanoStatic + scale_colour_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
    }
    volcanoStatic <- volcanoStatic + labs(
      x = "Log2 Fold Change",
      y = "-log10 P-value",
      fill = "Significance"
    ) +
      ylim(0, input$yLimVolcano) +
      xlim(-input$xLimVolcano, input$xLimVolcano) +
      theme_bw() +
      geom_vline(
        xintercept = c(
          -as.numeric(input$lfcVolcano),
          as.numeric(input$lfcVolcano)
        ),
        col = "black",
        linetype = "dashed"
      ) +
      geom_hline(
        yintercept = -log10(as.numeric(input$pThresholdVolcano)),
        col = "black",
        linetype = "dashed"
      ) +
      theme(
        axis.text.x = element_text(
          colour = "grey20", size = input$textSizeVolcano, angle = 0, hjust = .5,
          vjust = .5, face = "plain"
        ),
        axis.text.y = element_text(
          colour = "grey20", size = input$textSizeVolcano, angle = 0, hjust = 1,
          vjust = 0.5, face = "plain"
        ),
        axis.title.x = element_text(
          colour = "grey20", size = input$textSizeVolcano, angle = 0, hjust = .5,
          vjust = 0, face = "plain"
        ),
        axis.title.y = element_text(
          colour = "grey20", size = input$textSizeVolcano, angle = 90,
          hjust = .5, vjust = .5, face = "plain"
        ),
        legend.text = element_text(
          colour = "grey20", size = input$textSizeVolcano
        ),
        legend.title = element_text(
          colour = "grey20", size = input$textSizeVolcano
        ),
        title = element_text(colour = "grey20", size = input$textSizeVolcano),
        strip.text = element_text(size = input$textSizeVolcano),
        strip.text.x = element_text(size = input$textSizeVolcano),
        strip.text.y = element_text(size = input$textSizeVolcano)
      )

    if (input$volcanoShowGenes == TRUE & length(input$boxKeepBucketGenes) >= 1) {
      if (input$showBorderVolcano) {
        volcanoStatic <- volcanoStatic + geom_point(data = volcanoHighlights, aes(x = log2Ratio, y = log10p, fill = Status), size = as.numeric(input$dotSizeVolcano), alpha = 0.9, pch = 21)
        volcanoStatic <- volcanoStatic + scale_fill_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
      } else {
        volcanoStatic <- volcanoStatic + geom_point(data = volcanoHighlights, aes(x = log2Ratio, y = log10p, colour = Status), size = as.numeric(input$dotSizeVolcano), alpha = 0.9, pch = 16)
        volcanoStatic <- volcanoStatic + scale_colour_manual(breaks = names(volcanoColours), values = as.character(volcanoColours))
      }
      volcanoStatic <- volcanoStatic +
        geom_label_repel(
          data = volcanoHighlights,
          aes(label = genelabels),
          force = 3,
          nudge_y = 0.2,
          nudge_x = 0.2,
          color = "black",
          size = (input$textSizeVolcano / 3),
          max.overlaps = 500,
          fontface = "bold"
        )
    }
    # if (!input$showLinesVolcano) {
    #   volcanoStatic <- volcanoStatic +
    #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # }
    if (!input$showAxesVolcano) {
      volcanoStatic <- volcanoStatic +
        theme(
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank()
        )
    }
    if (input$boldVolcano) {
      volcanoStatic <- volcanoStatic +
        theme(
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.border = element_rect(linewidth = 1.5),
          axis.ticks = element_line(linewidth = 1.3)
        )
    }

    # Render the static plot
    output$volcanoStatic <- renderPlot(
      {
        volcanoStatic
      },
      width = as.numeric(input$figWidthVolcano),
      height = as.numeric(input$figHeightVolcano)
    )

    output$volcanoBrushTable <- renderDataTable({
      volcanoTableFull2 <- volcanoTableFull %>% dplyr::select(gene_name, description, log2Ratio, log2Ratio_full, log10p, log10p_full, Status)
      DT::datatable(
        data = brushedPoints(volcanoTableFull2, input$volcanoBrush),
        colnames = c("Gene name", "Description", "Log2 Ratio (plot)", "Log2 Ratio (full)", "-Log10p (plot)", "-Log10p (full)", "Status"),
        filter = "top",
        caption = "Click and drag over dots on the volcano plot to see those genes in this table.",
        rownames = F) %>%
        DT::formatSignif(c("log2Ratio", "log2Ratio_full", "log10p", "log10p_full"), digits = 3) %>%
        DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
        DT::formatStyle(columns = c("log2Ratio", "log2Ratio_full"), color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold") %>%
        DT::formatStyle(columns = c("log10p", "log10p_full"), color = styleInterval(cuts = -log10(as.numeric(input$pThresholdVolcano)), values = c("black", "green")), fontWeight = "bold") %>%
        DT::formatStyle(columns = "Status", color = styleEqual(levels = names(volcanoColours), values = as.character(volcanoColours)), fontWeight = "bold")
    })

    # Download button for the plot
    output$dlVolcanoPlotButton <- downloadHandler(
      filename = function() {
        paste0(design, "_volcano.pdf")
      },
      content = function(file) {
        ggsave(
          filename = file,
          plot = volcanoStatic,
          width = (as.numeric(input$figWidthVolcano) / 3.2),
          height = (as.numeric(input$figHeightVolcano) / 3.2),
          limitsize = FALSE, units = "mm"
        )
      }
    )

    output$dlVolcanoDFButton <- downloadHandler(
      filename = function() {
        paste0(design, "_volcano.xlsx")
      },
      content = function(file) {
        openxlsx::write.xlsx(volcanoTable, file = file)
      }
    )
  }
)