req(!is.null(inputDataReactive()$dataType))

if (inputDataReactive()$dataType == "RNASeq") {
  
  # Create the download link in XLSX format:
  output$dlDEGTable <- downloadHandler(
    filename = function() {
      basename(inputDataReactive()$degHTML)
    },
    content = function(file) {
      file.copy(from = inputDataReactive()$degHTML, to = file)
    }
  )
  #
  # Output the results summary table
  output$resultsTable <- DT::renderDataTable({
    DT::datatable(
      data = makeSignificantCounts(inputDataReactive()$se),
      rownames = TRUE
    )
  })
  output$DEGDesign <- renderText({
    inputDataReactive()$design
  })
  
  colsToPlotR <- c(
    "gene_id", "gene_name", "description", "type", "log2Ratio", "pValue", "fdr", "TPM_mean", paste0("TPM_", inputDataReactive()$param$sampleGroup, "_Mean"), paste0("TPM_", inputDataReactive()$param$refGroup, "_Mean")
  )
  colNamesToPlotR <- c(
    "Entrez ID", "Gene Symbol", "Type", "Description", "Log2 Ratio", "Raw p-value", "FDR", "Mean TPM", paste0(inputDataReactive()$param$sampleGroup, " Mean TPM"), paste0(inputDataReactive()$param$refGroup, " Mean TPM")
  )
  # Output the interactive data table
  output$degTable <- DT::renderDataTable({
    DT::datatable(
      data = inputDataReactive()$seqAnno %>% dplyr::select(all_of(colsToPlotR)),
      filter = "top",
      class = "cell-border stripe",
      rownames = FALSE,
      colnames = colNamesToPlotR
    ) %>%
      DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
      DT::formatSignif(columns = c("pValue", "fdr", "TPM_mean", paste0("TPM_", inputDataReactive()$param$sampleGroup, "_Mean"), paste0("TPM_", inputDataReactive()$param$refGroup, "_Mean")), digits = 3) %>%
      DT::formatSignif(columns = "log2Ratio", digits = 5) %>%
      DT::formatStyle(columns = "log2Ratio", color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold") %>%
      DT::formatStyle(columns = c("pValue", "fdr"), color = styleInterval(cuts = 0.05, values = c("green", "black")), fontWeight = "bold")
  })
  
  observeEvent(input$degTable_rows_selected, {
    seqAnnoSelected <- inputDataReactive()$seqAnno[input$degTable_rows_selected, ]
    output$dlDEGTable2 <- downloadHandler(
      filename = function() {
        paste0("selected_result_", inputDataReactive()$param$comparison, ".xlsx")
      },
      content = function(file) {
        openxlsx::write.xlsx(seqAnnoSelected, file = file)
      }
    )
  })
  observeEvent({
    input$featureCalcPType
    input$featureCalcPValue
    input$featureCalcLFC
  }, {
    if (input$featureCalcPType == "Raw") {
      featureCalcPType <- "pValue"
    } else {
      featureCalcPType <- "fdr"
    }
    seqAnnoCalc <- inputDataReactive()$seqAnno[inputDataReactive()$seqAnno$usedInTest == TRUE,] %>% dplyr::select(all_of(c(featureCalcPType, "log2Ratio")))
    seqAnnoCalc <- seqAnnoCalc[which(seqAnnoCalc[[featureCalcPType]] <= as.numeric(input$featureCalcPValue) & abs(seqAnnoCalc[["log2Ratio"]]) >= as.numeric(input$featureCalcLFC)),]
    seqAnnoCalcTable <- data.frame(
      "Up" = nrow(seqAnnoCalc[which(seqAnnoCalc[["log2Ratio"]] >= as.numeric(input$featureCalcLFC)),]),
      "Down" = nrow(seqAnnoCalc[which(seqAnnoCalc[["log2Ratio"]] <= -as.numeric(input$featureCalcLFC)),])
    )
    output$featureCalcTable <- function() {
      seqAnnoCalcTable %>%
        kable(
          format = "html"
        ) %>%
        kable_styling(
          bootstrap_options = "striped",
          full_width = FALSE,
          position = "left"
        )
    }
  })
}


if (inputDataReactive()$dataType == "proteomics") {
  observeEvent(
    {
      input$contrastSelected
    },
    {
      design <- input$contrastSelected
      designLevels <- input$contrastSelected %>% strsplit(split = "_vs_") %>% .[[1]]
      output$DEGDesign <- renderText({design})
      
      # Create the download link in XLSX format:
      output$dlDEGTable <- downloadHandler(
        filename = function() {
          paste0(design, "full_results_file.xlsx")
        },
        content = function(file) {
          tableToWrite <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]
          colnames(tableToWrite) <- c("protein_Id", "fasta.id", "diff", "p.value", "FDR", "description", "nrPeptides")
          openxlsx::write.xlsx(tableToWrite, file = file)
        }
      )
      
      seqAnnoDEG <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]
      
      # Output the results summary table
      signifCountsTable <- data.frame(
        check.names = F,
        row.names = c("p < 0.1", "p < 0.05", "p < 0.01", "p < 0.001", "p < 1e-04", "p < 1e-05")
      )
      significants <- unlist(setNames(lapply(c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), function(f) {
        isSig <- seqAnnoDEG$gene_name[seqAnnoDEG$pValue < f]
        length(isSig)
      }), paste0("p < ", c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001))))
      FDR <- unlist(lapply(c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), function(f) {
        isSig <- seqAnnoDEG$gene_name[seqAnnoDEG$pValue < f]
        if (length(isSig) >= 1) {
          max(seqAnnoDEG$fdr[seqAnnoDEG$pValue < f])
        } else {
          1
        }
      }))
      log2s <- setNames(lapply(c(0.5, 1, 1.5, 2, 3, 4, 8, 10), function(f) {
        unlist(lapply(c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), function(g) {
          isSig <- seqAnnoDEG$gene_name[which(seqAnnoDEG$pValue < g & abs(seqAnnoDEG$log2Ratio) >= f)]
          length(isSig)
        }))
      }), paste0("Log2 FC >= ", c(0.5, 1, 1.5, 2, 3, 4, 8, 10)))
      signifCountsTable <- cbind(signifCountsTable, significants, FDR, log2s)
      signifCountsTable$FDR <- round(signifCountsTable$FDR, digits = 3)
      output$resultsTable <- DT::renderDataTable({
        DT::datatable(
          data = signifCountsTable,
          rownames = TRUE
        )
      })
      
      colsToPlot1 <- c(
        "gene_name", "description", "nrPeptides", "log2Ratio", "pValue", "fdr"
      )
      colNamesToPlot1 <- c(
        "Protein ID", "Description", "Number of peptides", "Log2 Ratio", "Raw p-value", "FDR"
      )
      
      # Output the interactive data table
      output$degTable <- DT::renderDataTable({
        DT::datatable(
          data = seqAnnoDEG %>% dplyr::select(all_of(colsToPlot1)),
          filter = "top",
          class = "cell-border stripe",
          rownames = FALSE,
          colnames = colNamesToPlot1
        ) %>%
          DT::formatStyle(columns = colnames(.$x$data), `font-size` = "14px") %>%
          DT::formatSignif(columns = c("pValue", "fdr"), digits = 3) %>%
          DT::formatSignif(columns = "log2Ratio", digits = 5) %>%
          DT::formatStyle(columns = "log2Ratio", color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold") %>%
          DT::formatStyle(columns = c("pValue", "fdr"), color = styleInterval(cuts = 0.05, values = c("green", "black")), fontWeight = "bold")
      })
      
      observeEvent(input$degTable_rows_selected, {
        seqAnnoDEGSelected <- seqAnnoDEG[input$degTable_rows_selected, ]
        colnames(seqAnnoDEGSelected) <- c("protein_Id", "fasta.id", "diff", "p.value", "FDR", "description", "nrPeptides")
        output$dlDEGTable2 <- downloadHandler(
          filename = function() {
            paste0("selected_result_", design, ".xlsx")
          },
          content = function(file) {
            openxlsx::write.xlsx(seqAnnoDEGSelected, file = file)
          }
        )
      })
      
      observeEvent({
        input$featureCalcPType
        input$featureCalcPValue
        input$featureCalcLFC
      }, {
        if (input$featureCalcPType == "Raw") {
          featureCalcPType <- "pValue"
        } else {
          featureCalcPType <- "fdr"
        }
        seqAnnoCalc <- seqAnnoDEG %>% dplyr::select(all_of(c(featureCalcPType, "log2Ratio")))
        seqAnnoCalc <- seqAnnoCalc[which(seqAnnoCalc[[featureCalcPType]] <= as.numeric(input$featureCalcPValue) & abs(seqAnnoCalc[["log2Ratio"]]) >= as.numeric(input$featureCalcLFC)),]
        seqAnnoCalcTable <- data.frame(
          "Up" = nrow(seqAnnoCalc[which(seqAnnoCalc[["log2Ratio"]] >= as.numeric(input$featureCalcLFC)),]),
          "Down" = nrow(seqAnnoCalc[which(seqAnnoCalc[["log2Ratio"]] <= -as.numeric(input$featureCalcLFC)),])
        )
        output$featureCalcTable <- function() {
          seqAnnoCalcTable %>%
            kable(
              format = "html"
            ) %>%
            kable_styling(
              bootstrap_options = "striped",
              full_width = FALSE,
              position = "left"
            )
        }
      })
    }
  )
}
