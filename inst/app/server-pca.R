if (inputDataReactive()$dataType == "RNASeq") {
  design <- inputDataReactive()$design
  output$pcaDesign <- renderText({design})
  updateSelectInput(session = session, inputId = "pcaCounts", choices = names(inputDataReactive()$countList), selected = "Normalised + Log2")
} else if (inputDataReactive()$dataType == "proteomics") {
  design <- input$contrastSelected
  output$pcaDesign <- renderText({design})
  updateSelectInput(session = session, inputId = "pcaCounts", choices = names(inputDataReactive()$countList), selected = "transformedData")
}

updateSelectInput(
  session = session,
  inputId = "pcaBatch",
  choices = c("None", inputDataReactive()$factors),
  selected = "None"
)

lapply(c(1:2), function(i) {
  selected <- switch(
    i,
    "1" = inputDataReactive()$factors[1],
    "2" = "None"
  )
  updateSelectInput(
    session = session,
    inputId = paste0("pcaFactor", i),
    choices = c("None", inputDataReactive()$factors),
    selected = selected
  )
})
observeEvent({
  input$pcaFactor1
  }, ignoreNULL = FALSE, ignoreInit = TRUE, {
    updateCheckboxGroupInput(
      session = session,
      inputId = paste0("pcaGroups"),
      choices = levels(as.factor(inputDataReactive()$dataset[, input$pcaFactor1])),
      selected = levels(as.factor(inputDataReactive()$dataset[, input$pcaFactor1]))
    )
})

observeEvent(input[["pcaFactor1"]], ignoreInit = T, {
  updateCheckboxGroupInput(
    session = session,
    inputId = "pcaGroups",
    choices = levels(as.factor(inputDataReactive()$dataset[[input[["pcaFactor1"]]]])),
    selected = levels(as.factor(inputDataReactive()$dataset[[input[["pcaFactor1"]]]]))
  )
})

if (inputDataReactive()$dataType == "proteomics") {
  req(!is.null(inputDataReactive()$seqAnnoList))
  if ("nrPeptides" %in% colnames(inputDataReactive()$seqAnnoList[[1]])) {
    output$nrPeptidesPCAUI <- renderUI({
      sliderInput(inputId = "nrPeptidePCA", label = "Minimum number of peptides", min = 0, max = 20, value = 2, step = 1, width = "85%")
    })
    outputOptions(output, "nrPeptidesPCAUI", suspendWhenHidden = FALSE)
  } else {
    output$nrPeptidesPCAUI <- renderUI({ NULL })
  }
}

observeEvent(
  {
    input$pcaFactor1
    input$pcaX
    input$pcaY
    input$pcaTopN
    input$pcaGroups
    input$pcaFactor2
    input$pcaCounts
    input$pcaAxesProp
    input$pcaBatch
    input$pcaCentre
    input$pcaScale
    input$textSizePCA
    input$showLinesPCA
    input$showAxesPCA
    input$boldPCA
    input$pcaShowNames
    input$pcaLog2
    input$dotSizePCA
    input$alphaPCA
    input$dotBorderPCA
    input$pcaLabelMaxOverlap
    input$geneLabelNudgePCAX
    input$geneLabelNudgePCAY
    input$geneLabelSizePCA
    input$downloadFormatPCA
    input$dpiPCA
    input$filnamePCA
    input$nrPeptidePCA
    input$pcaAddEllipses
    input$pcaEllipsesAlpha
    lapply(seq_along(inputDataReactive()$factorLevels), function (i) {
      input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]])]]
    })
  },
  ignoreInit = TRUE, ignoreNULL = FALSE,
  {
    
    req(!is.null(inputDataReactive()$dataset))
  
    # Get the colours for the condition levels
    coloursPCA <- setNames(lapply(input$pcaGroups, function(k) {
      paste0(col2hex(input[[paste0("GroupColour", input$pcaFactor1, k)]]), "FF")
    }), input$pcaGroups)
    
    # Keep only the groups selected:
    datasetPCA <- inputDataReactive()$dataset
    datasetPCA <- datasetPCA[datasetPCA[[input$pcaFactor1]] %in% input$pcaGroups, ]
    countsPCA <- inputDataReactive()$countList[[input$pcaCounts]]
    countsPCA <- countsPCA[, rownames(datasetPCA)]
    
    # If proteomics and minimum number of peptides selected, subset the counts
    if (inputDataReactive()$dataType == "proteomics") {
      if (!is.null(input$nrPeptidePCA)) {
        countsPCA <- countsPCA[which(rownames(countsPCA) %in% inputDataReactive()$seqAnnoList[[1]]$gene_name[which(inputDataReactive()$seqAnnoList[[1]]$nrPeptides >= input$nrPeptidePCA)]), ]
      }
    }
    
    if (input$pcaLog2) {
      if (input$pcaCounts %in% c("TPM", "FPKM", "Raw", "Normalised")) {
        countsPCA <- log2(countsPCA+1)
      }
    }
    
    if (input$pcaBatch != "None") {
      countsPCA <- limma::removeBatchEffect(countsPCA, datasetPCA[[input$pcaBatch]])
    }
    countsPCA <- as.data.frame(t(countsPCA))
    # Get the n genes with greatest standard deviation for the PCA plot
    if (ncol(countsPCA) > input$pcaTopN) {
      countsPCA <- countsPCA[, names(head(sort(sapply(countsPCA, sd, na.rm = FALSE), decreasing = TRUE), n = input$pcaTopN))]
    }
    req(length(dim(countsPCA)) == 2)
    req(ncol(countsPCA) > 2)
    
    # Get the PCA results
    pcaResults <- tryCatch({ prcomp(na.omit(countsPCA), center = input$pcaCentre, scale. = input$pcaScale) }, error = function(e){ cat("ERROR :", conditionMessage(e), "\n") })
    if (!is.null(pcaResults)) {
      
      pc_eigenvalues <- pcaResults$sdev^2
      pc_eigenvalues <- tibble(PC = paste0("PC", factor(1:length(pc_eigenvalues))), variance = pc_eigenvalues) %>%
        mutate(pct = variance / sum(variance) * 100) %>%
        mutate(pct_cum = cumsum(pct))
      pc_eigenvalues$PC <- factor(pc_eigenvalues$PC, levels = pc_eigenvalues$PC)
      
      # Get the gene loadings
      pc_loadings <- pcaResults$rotation
      pc_loadings <- pc_loadings %>%
        as_tibble(rownames = "gene")
      pc_scores <- pcaResults$x
      pc_scores <- pc_scores %>%
        as_tibble(rownames = "sample")
      pc_scores <- left_join(pc_scores, datasetPCA, by = c("sample" = "names"))
      
      updateSelectInput(session = session, inputId = "pcaX", choices = pc_eigenvalues$PC, selected = input$pcaX)
      updateSelectInput(session = session, inputId = "pcaY", choices = pc_eigenvalues$PC, selected = input$pcaY)
      
      # Output the scree plot
      output$pcaScree <- renderPlot(
        {
          pc_eigenvalues2 <- pc_eigenvalues
          pc_eigenvalues2$PC <- gsub("PC", "", pc_eigenvalues2$PC)
          pc_eigenvalues2$PC <- factor(pc_eigenvalues2$PC, levels = pc_eigenvalues2$PC)
          if (nrow(pc_eigenvalues2) > 15) {
            pc_eigenvalues2 <- pc_eigenvalues2[1:15, ]
          }
          pcaScree <- pc_eigenvalues2 %>%
            ggplot(aes(x = PC)) +
            geom_col(aes(y = pct)) +
            geom_line(aes(y = pct_cum, group = 1)) +
            geom_point(aes(y = pct_cum)) +
            labs(x = "Principal component", y = "Fraction variance explained") +
            theme_classic(base_size = as.numeric(input$textSizePCA))
          pcaScree
        },
        width = 400,
        height = 400
      )
      # Output the variance table
      output$pcaVars <- DT::renderDataTable({
        datatable(pc_eigenvalues[, c("PC", "pct", "pct_cum")], rownames = F, colnames = c("PC", "Variance explained", "Cumulative variance")) %>% formatRound(c("pct", "pct_cum"), digits = 3)
      })
      # Output the loading table
      top_genes <- pc_loadings %>%
        dplyr::select(all_of(c("gene", input$pcaX, input$pcaY))) %>%
        pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>%
        group_by(PC) %>%
        arrange(desc(abs(loading)))
      output$pcaLoadings <- DT::renderDataTable({
        datatable(top_genes, rownames = F) %>% formatRound("loading", digits = 3)
      })
      output$dlPCALoadDFButton <- downloadHandler(filename = function() {
        paste0(design, "_PCA_loadings.xlsx")
      },
      content = function(file) {
        openxlsx::write.xlsx(as.data.frame(pc_loadings), file = file)
      }
      )
      
      if (input$pcaFactor2 == "None") {
        pcaPlot <- ggplot(data = pc_scores, aes(x = .data[[input$pcaX]], y = .data[[input$pcaY]], fill = .data[[input$pcaFactor1]]))
        if (input$pcaAddEllipses) {
          pcaPlot <- pcaPlot + stat_ellipse(aes(colour = .data[[input$pcaFactor1]], fill = .data[[input$pcaFactor1]]), geom = "polygon", alpha = input$pcaEllipsesAlpha, show.legend = FALSE)
          pcaPlot <- pcaPlot + scale_colour_manual(breaks = names(coloursPCA), values = as.character(coloursPCA))
        }
        pcaPlot <- pcaPlot + geom_point(aes_string(fill = input$pcaFactor1), size = input$dotSizePCA, alpha = input$alphaPCA, shape = 21, stroke = input$dotBorderPCA)
        pcaPlot <- pcaPlot + guides(fill = guide_legend(override.aes = list(shape = 21, size = input$dotSizePCA, stroke = input$dotBorderPCA)))
      } else {
        pcaPlot <- ggplot(data = pc_scores, aes(x = .data[[input$pcaX]], y = .data[[input$pcaY]], shape = .data[[input$pcaFactor2]]))
        if (input$pcaAddEllipses) {
          pcaPlot <- pcaPlot + stat_ellipse(aes(colour = .data[[input$pcaFactor1]], fill = .data[[input$pcaFactor1]]), geom = "polygon", alpha = input$pcaEllipsesAlpha, show.legend = FALSE)
          pcaPlot <- pcaPlot + scale_colour_manual(breaks = names(coloursPCA), values = as.character(coloursPCA))
        }
        pcaPlot <- pcaPlot + geom_point(aes_string(fill = input$pcaFactor1), size = input$dotSizePCA, alpha = input$alphaPCA, stroke = input$dotBorderPCA)
        pcaPlot <- pcaPlot + scale_shape_manual(values = c(rep(c(21, 22, 23, 24, 25, 8, 3, 4), times = 10))[1:nlevels(as.factor(datasetPCA[[input$pcaFactor2]]))])
        pcaPlot <- pcaPlot + guides(fill = guide_legend(override.aes = list(shape = 21)))
      }
      if (input$pcaAxesProp) {
        pcaPlot <- pcaPlot + coord_fixed(ratio = pc_eigenvalues$pct[pc_eigenvalues$PC == input$pcaY] / pc_eigenvalues$pct[pc_eigenvalues$PC == input$pcaX])
      }
      pcaPlot <- pcaPlot + labs(
        x = paste0(input$pcaX, ": ", round(pc_eigenvalues$pct[pc_eigenvalues$PC == input$pcaX]), "% variance"),
        y = paste0(input$pcaY, ": ", round(pc_eigenvalues$pct[pc_eigenvalues$PC == input$pcaY]), "% variance"),
        fill = input$pcaFactor1, shape = input$pcaFactor2
      )
      if (input$pcaShowNames) {
        pcaPlot <- pcaPlot +
          geom_label_repel(
            aes(label = sample, fill = .data[[input$pcaFactor1]]),
            force = 2,
            nudge_y = as.numeric(input$geneLabelNudgePCAY)/10,
            nudge_x = as.numeric(input$geneLabelNudgePCAX)/10, 
            size = (input$textSizePCA / 3),
            max.overlaps = input$pcaLabelMaxOverlap,
            fontface = "bold", 
            show.legend = F,
            show_guide = FALSE
          ) 
      }
      pcaPlot <- pcaPlot + scale_fill_manual(breaks = names(coloursPCA), values = as.character(coloursPCA))
      pcaPlot <- pcaPlot + theme_bw()
      face <- switch(input$boldPCA, "TRUE" = "bold", "FALSE" = "plain")
      pcaPlot <- pcaPlot + theme(
        axis.text.x = element_text(size = input$textSizePCA, angle = 0, hjust = .5, vjust = .5, face = face, colour = "black"),
        axis.text.y = element_text(size = input$textSizePCA, angle = 0, hjust = 1, vjust = 0.5, face = face, colour = "black"),
        axis.title.x = element_text(size = input$textSizePCA, angle = 0, hjust = .5, vjust = 0, face = face),
        axis.title.y = element_text(size = input$textSizePCA, angle = 90, hjust = .5, vjust = .5, face = face),
        legend.text = element_text(size = input$textSizePCA),
        legend.title = element_text(size = input$textSizePCA),
        title = element_text(size = input$textSizePCA)
      )
      if (!input$showLinesPCA) {
        pcaPlot <- pcaPlot +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }
      if (!input$showAxesPCA) {
        pcaPlot <- pcaPlot +
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
      if (all(input$boldPCA, input$showAxesPCA)) {
        pcaPlot <- pcaPlot +
          theme(
            panel.border = element_rect(linewidth = 1.5),
            axis.ticks = element_line(linewidth = 1.3)
          )
      }
      figuresDataReactive$pcaStatic <- pcaPlot

      # Paired Plot 
      pairedPlot <- ggpairs(
        data = pc_scores, 
        columns = 2:6, 
        aes_string(fill = input$pcaFactor1, alpha = input$alphaPCA), 
        upper = list(continuous = wrap("points", alpha = input$alphaPCA, size = input$dotSizePCA, alpha = input$alphaPCA, pch = 21, stroke = input$dotBorderPCA)),
        lower = list(continuous = wrap("points", alpha = input$alphaPCA, size = input$dotSizePCA, alpha = input$alphaPCA, pch = 21, stroke = input$dotBorderPCA)),
        diag = list(continuous = wrap("densityDiag", size = 0.3), discrete = "barDiag", na = "naDiag")
      ) + 
        theme_bw(base_size = input$textSizePCA) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_colour_manual(values = as.character(coloursPCA), breaks = names(coloursPCA)) +
        scale_fill_manual(values = as.character(coloursPCA), breaks = names(coloursPCA))
      
      figuresDataReactive$pcaPaired <- pairedPlot
      
      output$dlPCAPlotDFButton <- downloadHandler(
        filename = function() {
          paste0(design, "_PCA_coords.xlsx")
        },
        content = function(file) {
          openxlsx::write.xlsx(pc_scores, file = file)
        }
      )
      
      output$pcaBrushTable <- renderDataTable({
        brushed_df <- brushedPoints(pc_scores, input$pcaBrush, xvar = input$pcaX, yvar = input$pcaY)
        DT::datatable(
          data = brushed_df[, c("sample", factorNames, input$pcaX, input$pcaY)],
          filter = "top",
          caption = "Click and drag over dots on the PCA plot to see those samples in this table.",
          rownames = F) %>% 
          DT::formatRound(columns = c(input$pcaX, input$pcaY), digits = 3)
      })
      
    } else {
      output$pcaNAWarning <- renderUI({
        print("I'm sorry! I think there are too many NAs in your count data for me to be able to calculate PCA results. If you think there is a mistake, please reach out to your bioinformatics coach for more help!")
      })
    }
})

output$pcaStatic <- renderPlot(
  {
    req(!is.null(figuresDataReactive$pcaStatic))
    figuresDataReactive$pcaStatic
  },
  width = function(){as.numeric(input$figWidthPCA)},
  height = function(){as.numeric(input$figHeightPCA)}
)
output$dlPCAPlotButton <- downloadHandler(
  filename = function() {
    paste(input$filnamePCA, design, tolower(input$downloadFormatPCA), sep = ".")
  },
  content = function(file) {
    if (input$downloadFormatPCA == "PDF") {
      pdf(file = file, width = as.numeric(input$figWidthPCA/80), height = as.numeric(input$figHeightPCA/80))
    } else if (input$downloadFormatPCA == "SVG") {
      svg(file = file, width = as.numeric(input$figWidthPCA/80), height = as.numeric(input$figHeightPCA/80))
    } else if (input$downloadFormatPCA == "PNG") {
      png(filename = file, width = as.numeric(input$figWidthPCA/80), height = as.numeric(input$figHeightPCA/80), units = "in", res = as.numeric(input$dpiPCA))
    }
    print(figuresDataReactive$pcaStatic)
    dev.off()
  }
)

output$pcaPaired <- renderPlot(
  {
    req(!is.null(figuresDataReactive$pcaPaired))
    figuresDataReactive$pcaPaired
  },
  width = function(){as.numeric(input$figWidthPCA)},
  height = function(){as.numeric(input$figHeightPCA)}
)
output$dlPCAPairButton <- downloadHandler(
  filename = function() {
    paste(input$filnamePCA, ".paired.", tolower(input$downloadFormatPCA), sep = ".")
  },
  content = function(file) {
    if (input$downloadFormatPCA == "PDF") {
      pdf(file = file, width = as.numeric(input$figWidthPCA/80), height = as.numeric(input$figHeightPCA/80))
    } else if (input$downloadFormatPCA == "SVG") {
      svg(file = file, width = as.numeric(input$figWidthPCA/80), height = as.numeric(input$figHeightPCA/80))
    } else if (input$downloadFormatPCA == "PNG") {
      png(filename = file, width = as.numeric(input$figWidthPCA/80), height = as.numeric(input$figHeightPCA/80), units = "in", res = as.numeric(input$dpiPCA))
    }
    print(figuresDataReactive$pcaPaired)
    dev.off()
  }
)