# Server logic for Multi DE Table

# Reactive expression for the combined table
multiDETableReactive <- reactive({
    req(multiDeReactive$results)
    req(length(multiDeReactive$results) > 0)
    
    results <- multiDeReactive$results

    # Get all unique genes across all contrasts
    allGenes <- unique(unlist(lapply(results, function(r) {
        seqAnno <- r$seqAnno
        if ("feature_id" %in% colnames(seqAnno)) {
            seqAnno$feature_id
        } else if ("gene_id" %in% colnames(seqAnno)) {
            seqAnno$gene_id
        } else {
            character(0)
        }
    })))

    # Start with gene names and IDs from first contrast
    firstSeqAnno <- results[[1]]$seqAnno
    baseTable <- data.frame(
        feature_id = if ("feature_id" %in% colnames(firstSeqAnno)) {
            firstSeqAnno$feature_id
        } else if ("gene_id" %in% colnames(firstSeqAnno)) {
            firstSeqAnno$gene_id
        } else {
            allGenes
        },
        stringsAsFactors = FALSE
    )

    if ("gene_id" %in% colnames(firstSeqAnno)) {
        baseTable$gene_id <- firstSeqAnno$gene_id
    }

    # Proteomics: gene symbol / protein ID alongside the (UniProt) feature_id,
    # skipped when identical to it
    for (col in c("gene_name", "protein_Id")) {
        if (identical(results[[1]]$dataType, "proteomics") &&
            col %in% colnames(firstSeqAnno) &&
            !identical(as.character(firstSeqAnno[[col]]), as.character(baseTable$feature_id))) {
            baseTable[[col]] <- firstSeqAnno[[col]]
        }
    }

    if ("description" %in% colnames(firstSeqAnno)) {
        baseTable$description <- firstSeqAnno$description
    }

    # Add columns for each contrast
    for (contrastName in names(results)) {
        seqAnno <- results[[contrastName]]$seqAnno

        # Create a key for matching (feature_id)
        seqAnno$match_key <- if ("feature_id" %in% colnames(seqAnno)) {
            seqAnno$feature_id
        } else if ("gene_id" %in% colnames(seqAnno)) {
            seqAnno$gene_id
        } else {
            NA
        }

        # Match genes
        matchIdx <- match(baseTable$feature_id, seqAnno$match_key)

        # Log2 fold change
        lfcCol <- grep("^(log2?[._ ]?(fold)?[._ ]?(change|ratio|fc)|logfc|log2Ratio)$",
                      colnames(seqAnno), value = TRUE, ignore.case = TRUE)
        if (length(lfcCol) > 0) {
            baseTable[[paste0(contrastName, "_log2FC")]] <- seqAnno[[lfcCol[1]]][matchIdx]
        }

        # P-value
        if ("pvalue" %in% colnames(seqAnno)) {
            baseTable[[paste0(contrastName, "_pvalue")]] <- seqAnno$pvalue[matchIdx]
        }

        # Adjusted p-value
        if ("padj" %in% colnames(seqAnno)) {
            baseTable[[paste0(contrastName, "_padj")]] <- seqAnno$padj[matchIdx]
        } else if ("fdr" %in% colnames(seqAnno)) {
            baseTable[[paste0(contrastName, "_padj")]] <- seqAnno$fdr[matchIdx]
        }

        # Base mean
        if ("baseMean" %in% colnames(seqAnno)) {
            baseTable[[paste0(contrastName, "_baseMean")]] <- seqAnno$baseMean[matchIdx]
        } else if ("AveExpr" %in% colnames(seqAnno)) {
            baseTable[[paste0(contrastName, "_baseMean")]] <- seqAnno$AveExpr[matchIdx]
        }
    }

    baseTable
})

output$multiDETable <- DT::renderDataTable({
    req(multiDETableReactive())
    combinedData <- multiDETableReactive()

    # Find column names for formatting
    log2fc_cols <- grep("_log2FC$", colnames(combinedData), value = TRUE)
    pval_cols <- grep("_pvalue$", colnames(combinedData), value = TRUE)
    padj_cols <- grep("_padj$", colnames(combinedData), value = TRUE)
    basemean_cols <- grep("_baseMean$", colnames(combinedData), value = TRUE)

    # Render table
    dt <- DT::datatable(
        combinedData,
        filter = "top",
        class = "cell-border stripe",
        rownames = FALSE,
        options = list(
            scrollX = TRUE,
            autoWidth = TRUE,
            pageLength = 25,
            lengthMenu = c(25, 50, 100, 500)
        )
    )

    # Format numeric columns
    if (length(log2fc_cols) > 0) {
        dt <- dt %>% DT::formatSignif(columns = log2fc_cols, digits = 3)
    }
    if (length(pval_cols) > 0) {
        dt <- dt %>% DT::formatSignif(columns = pval_cols, digits = 4)
    }
    if (length(padj_cols) > 0) {
        dt <- dt %>% DT::formatSignif(columns = padj_cols, digits = 4)
    }
    if (length(basemean_cols) > 0) {
        dt <- dt %>% DT::formatSignif(columns = basemean_cols, digits = 2)
    }

    # Apply styling
    if (length(pval_cols) > 0) {
        dt <- dt %>% DT::formatStyle(
            columns = pval_cols,
            color = styleInterval(cuts = 0.05, values = c("green", "black")),
            fontWeight = "bold"
        )
    }
    if (length(padj_cols) > 0) {
        dt <- dt %>% DT::formatStyle(
            columns = padj_cols,
            color = styleInterval(cuts = 0.05, values = c("green", "black")),
            fontWeight = "bold"
        )
    }
    if (length(log2fc_cols) > 0) {
        dt <- dt %>% DT::formatStyle(
            columns = log2fc_cols,
            color = styleInterval(cuts = 0, values = c("blue", "orange")),
            fontWeight = "bold"
        )
    }

    dt
})

# Download handler
output$dlTable_multiDE <- downloadHandler(
    filename = function() {
        generate_filename(
            base_name = input$filename_multiTable,
            design_name = "Multi",
            extension = "xlsx"
        )
    },
    content = function(file) {
        req(multiDETableReactive())
        export_excel_standard(
            df = multiDETableReactive(),
            filename = file,
            sheet_name = "DE Results",
            table_name = "DEResults"
        )
    }
)
