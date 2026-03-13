req(inputDataReactive()$dataType)

# Bucket render removed — now in server-featureBucket.R

if (inputDataReactive()$dataType == "RNASeq") {
  # Create the download link in XLSX format:
  # Create the download link in XLSX format:
  output$dlTable_deg <- downloadHandler(
    filename = function() {
      generate_filename(
        base_name = paste0(input$filename_deg, "_Full"),
        design_name = inputDataReactive()$param$comparison,
        extension = "xlsx"
      )
    },
    content = function(file) {
      export_excel_standard(
        df = inputDataReactive()$seqAnno,
        filename = file,
        sheet_name = "DE Results",
        table_name = "DEResults"
      )
    }
  )

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
    "gene_id", "feature_id", "description", "type", "log2Ratio", "pValue", "fdr", "TPM_mean", paste0("TPM_", inputDataReactive()$param$sampleGroup, "_Mean"), paste0("TPM_", inputDataReactive()$param$refGroup, "_Mean")
  )
  colNamesToPlotR <- c(
    "Ensembl ID", "Feature", "Description", "Type", "Log2 Ratio", "Raw p-value", "FDR", "Mean TPM", paste0(inputDataReactive()$param$sampleGroup, " Mean TPM"), paste0(inputDataReactive()$param$refGroup, " Mean TPM")
  )
  # Output the interactive data table
  dt1 <- DT::datatable(
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
  output$degTable <- DT::renderDataTable(dt1, server = TRUE, selection = list(target = "row"))
  proxy <- dataTableProxy("degTable")
  observeEvent(input$resetDEGTableSelectRows, {
    proxy %>% selectRows(NULL)
  })

  observeEvent(input$degTable_rows_selected, {
    seqAnnoSelected <- inputDataReactive()$seqAnno[input$degTable_rows_selected, ]
    output$dlTableSelected_deg <- downloadHandler(
      filename = function() {
        generate_filename(
          base_name = paste0(input$filename_deg, "_Selected"),
          design_name = inputDataReactive()$param$comparison,
          extension = "xlsx"
        )
      },
      content = function(file) {
        export_excel_standard(
          df = seqAnnoSelected,
          filename = file,
          sheet_name = "Selected Results",
          table_name = "SelectedResults"
        )
      }
    )
  })
  observeEvent(
    {
      input$featureCalcPType
      input$featureCalcPValue
      input$featureCalcLFC
    },
    {
      if (input$featureCalcPType == "Raw") {
        featureCalcPType <- "pValue"
      } else {
        featureCalcPType <- "fdr"
      }

      seqAnnoCalc <- inputDataReactive()$seqAnno[inputDataReactive()$seqAnno$usedInTest == TRUE, ]
      seqAnnoCalc <- seqAnnoCalc[which(seqAnnoCalc[[featureCalcPType]] <= as.numeric(input$featureCalcPValue) & abs(seqAnnoCalc[["log2Ratio"]]) >= as.numeric(input$featureCalcLFC)), ]

      seqAnnoCalcTable <- data.frame(
        "Up" = nrow(seqAnnoCalc[which(seqAnnoCalc[["log2Ratio"]] >= as.numeric(input$featureCalcLFC)), ]),
        "Down" = nrow(seqAnnoCalc[which(seqAnnoCalc[["log2Ratio"]] <= -as.numeric(input$featureCalcLFC)), ])
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

      # Download the results
      # Download the results
      output$dlTableFiltered_deg <- downloadHandler(
        filename = function() {
          generate_filename(
            base_name = input$filename_deg_filtered,
            extension = "xlsx"
          )
        },
        content = function(file) {
          export_excel_standard(
            df = seqAnnoCalc,
            filename = file,
            sheet_name = "DE Features",
            table_name = "DEFeatures"
          )
        }
      )
    }
  )
}


if (inputDataReactive()$dataType == "proteomics") {
  observeEvent(
    {
      input$contrastSelected
    },
    {
      design <- input$contrastSelected
      designLevels <- input$contrastSelected %>%
        strsplit(split = "_vs_") %>%
        .[[1]]
      output$DEGDesign <- renderText({
        design
      })

      # Create the download link in XLSX format:
      # Create the download link in XLSX format:
      output$dlTable_deg <- downloadHandler(
        filename = function() {
          generate_filename(
            base_name = paste0(input$filename_deg, "_Full"),
            design_name = design,
            extension = "xlsx"
          )
        },
        content = function(file) {
          tableToWrite <- as.data.frame(rowData(se)[[grep(input$contrastSelected, names(rowData(se)))]])
          export_excel_standard(
            df = tableToWrite,
            filename = file,
            sheet_name = "Full Results",
            table_name = "FullResults"
          )
        }
      )

      seqAnnoDEG <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]

      # Output the results summary table
      signifCountsTable <- data.frame(
        check.names = F,
        row.names = c("p < 0.1", "p < 0.05", "p < 0.01", "p < 0.001", "p < 1e-04", "p < 1e-05")
      )
      significants <- unlist(setNames(lapply(c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), function(f) {
        isSig <- seqAnnoDEG$feature_id[seqAnnoDEG$pValue < f]
        length(isSig)
      }), paste0("p < ", c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001))))
      FDR <- unlist(lapply(c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), function(f) {
        isSig <- seqAnnoDEG$feature_id[seqAnnoDEG$pValue < f]
        if (length(isSig) >= 1) {
          max(seqAnnoDEG$fdr[seqAnnoDEG$pValue < f])
        } else {
          1
        }
      }))
      log2s <- setNames(lapply(c(0.5, 1, 1.5, 2, 3, 4, 8, 10), function(f) {
        unlist(lapply(c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), function(g) {
          isSig <- seqAnnoDEG$feature_id[which(seqAnnoDEG$pValue < g & abs(seqAnnoDEG$log2Ratio) >= f)]
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

      # Predefined columns and labels (canonical names from standardise_columns)
      colsToPlot1 <- colnames(seqAnnoDEG %>% dplyr::select(any_of(c(
        "feature_id", "protein_Id", "peptide_Id", "site", "gene_name",
        "description", "nrPeptides", "modelName", "log2Ratio", "pValue", "fdr"
      ))))
      col_labels <- c(
        "feature_id" = "Feature ID", "protein_Id" = "Protein ID",
        "peptide_Id" = "Peptide", "site" = "Site",
        "gene_name" = "Gene Symbol", "description" = "Description",
        "nrPeptides" = "Nr Peptides", "modelName" = "Imputed Model",
        "log2Ratio" = "Log2 Ratio", "pValue" = "Raw p-value", "fdr" = "FDR"
      )
      colNamesToPlot1 <- unname(col_labels[colsToPlot1])

      # Output the interactive data table
      dt1 <- DT::datatable(
        data = seqAnnoDEG %>% dplyr::select(all_of(colsToPlot1)),
        filter = "top",
        class = "cell-border stripe",
        rownames = FALSE,
        colnames = colNamesToPlot1,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(list(extend = "colvis", text = "Column Visibility"))
        )
      )
      # Apply formatting using column indices (colnames= renames break string lookups)
      dt1 <- dt1 %>% DT::formatStyle(columns = seq_along(colsToPlot1), `font-size` = "14px")
      # Wrap long Feature IDs (e.g. protein groups with many semicolon-separated accessions)
      fid_idx <- which(colsToPlot1 == "feature_id")
      if (length(fid_idx) > 0) {
        dt1 <- dt1 %>% DT::formatStyle(columns = fid_idx, `word-break` = "break-all", `max-width` = "250px")
      }
      pval_idx <- which(colsToPlot1 %in% c("pValue", "fdr"))
      lfc_idx <- which(colsToPlot1 == "log2Ratio")
      model_idx <- which(colsToPlot1 == "modelName")
      if (length(pval_idx) > 0) {
        dt1 <- dt1 %>%
          DT::formatSignif(columns = pval_idx, digits = 3) %>%
          DT::formatStyle(columns = pval_idx, color = styleInterval(cuts = 0.05, values = c("green", "black")), fontWeight = "bold")
      }
      if (length(lfc_idx) > 0) {
        dt1 <- dt1 %>%
          DT::formatSignif(columns = lfc_idx, digits = 5) %>%
          DT::formatStyle(columns = lfc_idx, color = styleInterval(cuts = 0, values = c("blue", "darkorange")), fontWeight = "bold")
      }
      if (length(model_idx) > 0) {
        dt1 <- dt1 %>% DT::formatStyle(columns = model_idx, color = styleEqual(levels = c("Linear_Model_moderated", "Imputed_Mean_moderated"), values = c("black", "darkblue")))
      }
      output$degTable <- DT::renderDataTable(dt1, server = TRUE, selection = list(target = "row"))
      proxy <- dataTableProxy("degTable")
      observeEvent(input$resetDEGTableSelectRows, {
        proxy %>% selectRows(NULL)
      })

      observeEvent(input$degTable_rows_selected, {
        seqAnnoDEGSelected <- seqAnnoDEG[input$degTable_rows_selected, ]
        output$dlTableSelected_deg <- downloadHandler(
          filename = function() {
            generate_filename(
              base_name = paste0(input$filename_deg, "_Selected"),
              design_name = design,
              extension = "xlsx"
            )
          },
          content = function(file) {
            export_excel_standard(
              df = seqAnnoDEGSelected,
              filename = file,
              sheet_name = "Selected Results",
              table_name = "SelectedResults"
            )
          }
        )
      })

      observeEvent(
        {
          input$featureCalcPType
          input$featureCalcPValue
          input$featureCalcLFC
        },
        {
          if (input$featureCalcPType == "Raw") {
            featureCalcPType <- "pValue"
          } else {
            featureCalcPType <- "fdr"
          }
          seqAnnoCalc <- seqAnnoDEG # %>% dplyr::select(all_of(c(featureCalcPType, "log2Ratio")))
          seqAnnoCalc <- seqAnnoCalc[which(seqAnnoCalc[[featureCalcPType]] <= as.numeric(input$featureCalcPValue) & abs(seqAnnoCalc[["log2Ratio"]]) >= as.numeric(input$featureCalcLFC)), ]

          seqAnnoCalcTable <- data.frame(
            "Up" = nrow(seqAnnoCalc[which(seqAnnoCalc[["log2Ratio"]] >= as.numeric(input$featureCalcLFC)), ]),
            "Down" = nrow(seqAnnoCalc[which(seqAnnoCalc[["log2Ratio"]] <= -as.numeric(input$featureCalcLFC)), ])
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
          # Download the results
          # Download the results
          output$dlTableFiltered_deg <- downloadHandler(
            filename = function() {
              generate_filename(
                base_name = input$filename_deg_filtered,
                extension = "xlsx"
              )
            },
            content = function(file) {
              export_excel_standard(
                df = seqAnnoCalc,
                filename = file,
                sheet_name = "DE Features",
                table_name = "DEFeatures"
              )
            }
          )
        }
      )
    }
  )
}
