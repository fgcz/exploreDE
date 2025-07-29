# Check if GSEA tab is accessed without results
observeEvent(input$tabs, {
  if (input$tabs == "gseaTab") {
    if (is.null(inputDataReactive()$goObjects$gsea) || 
        length(inputDataReactive()$goObjects$gsea) == 0) {
      showModal(modalDialog(
        title = "No GSEA Results", 
        "No pathway analysis was performed for this dataset, so this tab is empty.",
        footer = modalButton("OK")
      ))
      return()
    }
    
    # Check if current selection has no results
    if (!is.null(input$gseaType) && 
        !is.null(inputDataReactive()$goObjects$gsea[[input$gseaType]])) {
      er <- inputDataReactive()$goObjects$gsea[[input$gseaType]]
      if (is.null(er) || nrow(er@result) == 0) {
        showModal(modalDialog(
          title = "No GSEA Results", 
          "There are no results for this selection."
        ))
      }
    }
  }
})

# Initialize reactive values for GSEA plots
gseaPlotReactive <- reactiveValues(
  er = NULL,
  cnetplot = NULL,
  heatmap = NULL,
  runningscore = NULL,
  ridgeplot = NULL,
  upset = NULL
)


# Get data objects
gsea <- inputDataReactive()$goObjects$gsea
se <- inputDataReactive()$se
design <- inputDataReactive()$design
param <- inputDataReactive()$param

# Update GSEA type choices
updateSelectInput(
  session = session, 
  inputId = "gseaType", 
  choices = names(gsea), 
  selected = names(gsea)[[1]]
)

# Render design output
output$gseaDesign <- renderText({design})

# Set default FDR threshold if missing
if (is.null(param$fdrThreshGSEA)) {
  param$fdrThreshGSEA <- 0.05
}

# Create GSEA overview table
gsea.overview <- ezFrame("Number of Significant Pathways" = numeric())
for (i in seq_along(names(gsea))) {
  testName <- names(gsea)[i]
  gsea.overview[testName, "Number of Significant Pathways"] <- nrow(gsea[[i]]@result)
}

output$gseaOverview <- function() {
  gsea.overview %>%
    kable(format = "html", caption = "GSEA:") %>%
    kable_styling(
      bootstrap_options = "striped",
      full_width = FALSE,
      position = "left"
    )
}

# Create GSEA thresholds table
gseaThresholds <- ezFrame("GSEA Cut-Offs" = character())
gseaThresholds["Output GSEA Term FDR", "GSEA Cut-Offs"] <- param$fdrThreshGSEA

output$gseaThreshOverview <- function() {
  gseaThresholds %>%
    kable(format = "html", caption = "Cut-Offs:") %>%
    kable_styling(
      bootstrap_options = "striped",
      full_width = FALSE,
      position = "left"
    )
}

# Reactive values for current GSEA selection
gseaToUse <- reactiveValues(er = NULL)

# Main observer for GSEA type changes
observeEvent(input$gseaType, {
  
  gseaToUse$er = NULL
  gseaToUse$er_symbol = NULL
  er <- gsea[[input$gseaType]]
  
  # Check if results exist and show modal if empty
  if (input$tabs == "gseaTab" && (is.null(er) || nrow(er@result) == 0)) {
    showModal(modalDialog(
      title = "No GSEA Results", 
      "There are no results for this selection."
    ))
    return()
  }
  
  # Process results if they exist
  if (!is.null(er) && nrow(er@result) > 0) {
    
    # Handle GO terms (BP, MF, CC) differently
    if (input$gseaType %in% c("BP", "MF", "CC")) {
      er@result$Label <- ""
      
      # Create labels for GO terms
      for (i in seq_along(er@result$Description)) {
        desc <- er@result$Description[i]
        label <- paste0(er@result$ID[i], "\n", desc)
        if (!is.na(label) && nchar(label) > 45) {
          label <- substr(label, 1, 45)
        }
        er@result$Label[i] <- label
      }
      colsToTable <- c("Label", "enrichmentScore", "NES", "qvalue", "geneName")
      # Convert gene IDs to gene names for core enrichment
      # Split all core enrichment strings and collect unique gene IDs
      all_core_genes <- er@result$core_enrichment %>%
        strsplit("/") %>%
        unlist() %>%
        unique()
      
      # Create lookup table for gene ID to gene name conversion
      gene_lookup <- setNames(
        metadata(se)$enrichInput$seqAnno$gene_name,
        metadata(se)$enrichInput$seqAnno$gene_id
      )
      er_symbol <- er
      # Convert core enrichment from gene IDs to gene names
      er_symbol@result$core_enrichment <- vapply(er_symbol@result$core_enrichment, function(genes_string) {
        gene_ids <- strsplit(genes_string, "/")[[1]]
        gene_names <- gene_lookup[gene_ids]
        gene_names <- gene_names[!is.na(gene_names)]
        paste(gene_names, collapse = "/")
      }, character(1), USE.NAMES = FALSE)
    } else {
      # For other pathway databases (e.g., C2)
      er_symbol <- er
      colsToTable <- c("Description", "enrichmentScore", "NES", "qvalue", "core_enrichment")
    }
  }
  
  # Prepare log2 ratio data for plotting
  log2RatioGSEA <- metadata(se)$enrichInput$seqAnno$log2Ratio
  names(log2RatioGSEA) <- metadata(se)$enrichInput$seqAnno$gene_name
  log2RatioGSEA <- sort(log2RatioGSEA, decreasing = TRUE)
  
  # Update reactive values
  gseaToUse$er <- er
  gseaToUse$er_symbol <- er_symbol
  gseaToUse$log2RatioGSEA <- log2RatioGSEA
  gseaToUse$colsToTable <- colsToTable
  
  # Clear existing plots
  for (x in names(gseaPlotReactive)) {
    gseaPlotReactive[[x]] <- NULL
  }
  
  # Reset table selection AFTER updating reactive values
  shinyjs::delay(100, {
    proxy <- DT::dataTableProxy("selectedTable_GSEA")
    DT::selectRows(proxy, NULL)
  })
  
}, ignoreNULL = FALSE, ignoreInit = TRUE)

# Render the main GSEA table
output$selectedTable_GSEA <- DT::renderDataTable({
  req(gseaToUse$er, gseaToUse$colsToTable)
  
  DT::datatable(
    data = as.data.frame(gseaToUse$er)[, gseaToUse$colsToTable],
    filter = "top", 
    class = "cell-border stripe",
    rownames = FALSE, 
    caption = "Click pathways in this table to add them to the figures on the right"
  ) %>%
    DT::formatStyle(
      columns = colnames(.$x$data), 
      `font-size` = "12px"
    ) %>%
    DT::formatSignif(
      columns = c("enrichmentScore", "NES", "qvalue"),
      digits = 3
    ) %>%
    DT::formatStyle(
      columns = c("enrichmentScore", "NES"),
      color = styleInterval(cuts = 0, values = c("blue", "darkorange")),
      fontWeight = "bold"
    ) %>% 
    DT::formatStyle(
      columns = "qvalue",
      color = styleInterval(cuts = param$fdrThreshGSEA, values = c("green", "black")),
      fontWeight = "bold"
    )
})

# Download handler for GSEA results
output$dlTableGSEA <- downloadHandler(
  filename = function() {
    paste0(design, "_GSEA_", input$gseaType, ".xlsx")
  },
  content = function(file) {
    openxlsx::write.xlsx(
      x = as.data.frame(gsea[[input$gseaType]]@result), 
      file = file
    )
  }
)

# Observer for plot generation when rows are selected
observeEvent({
  input$textSizeGSEA
  input$nodeSizeGSEA
  input$selectedTable_GSEA_rows_selected
  input$scaleLimGSEA
  input$gseaMaxOver
  input$gseaDodgeY
  input$gseaDodgeX
  input$nodeBorderGSEA
  input$ridgePlotColourByGSEA
  input$gseaLabelAlpha
}, {
  
  req(gseaToUse$er_symbol)
  current_er <- isolate(gseaToUse$er)
  current_er_symbol <- isolate(gseaToUse$er_symbol)
  req(!is.null(input$selectedTable_GSEA_rows_selected))
  req(current_er)
  
  # Subset enrichment result for selected pathways
  erRP <- clusterProfiler::slice(current_er, input$selectedTable_GSEA_rows_selected)
  labelColumn <- ifelse(input$gseaType == "C2", "Description", "Label")
  
  # Network plot
  cn <- enrichplot::cnetplot(
    x = current_er_symbol,
    node_label = "none",
    foldChange = gseaToUse$log2RatioGSEA,
    showCategory = current_er_symbol@result$Description[input$selectedTable_GSEA_rows_selected],
    cex.params = list(
      category_label = input$textSizeGSEA / 15, 
      gene_label = input$textSizeGSEA / 18, 
      gene_node = input$nodeSizeGSEA / 10, 
      category_node = input$nodeSizeGSEA / 8
    )
  )
  
  # Apply colour scale limits
  cn$data$foldChange[!is.na(cn$data$foldChange) & cn$data$foldChange > input$scaleLimGSEA] <- input$scaleLimGSEA
  cn$data$foldChange[!is.na(cn$data$foldChange) & cn$data$foldChange < -input$scaleLimGSEA] <- -input$scaleLimGSEA
  
  # Add custom styling
  cn <- cn + 
    geom_point(data = cn$data, aes(x = x, y = y, size = size), shape = 21, stroke = input$nodeBorderGSEA) +
    scale_colour_gradient2(
      name = "fold change", 
      low = "deepskyblue4", 
      mid = "white", 
      high = "firebrick", 
      midpoint = 0
    )
  
  # Prepare data for labels
  cd1 <- cn$data
  cd1$name[!is.na(cd1$foldChange)] <- NA
  cd1$Label <- NA
  cd1$Label[seq_along(input$selectedTable_GSEA_rows_selected)] <- 
    current_er_symbol@result[[labelColumn]][input$selectedTable_GSEA_rows_selected]
  
  cd2 <- cn$data
  cd2$name[is.na(cd2$foldChange)] <- NA
  
  # Add gene labels if requested
  if (input$showGeneLabelsGSEA) {
    cn <- cn + geom_text_repel(
      data = cd2, 
      aes(label = name, x = x, y = y), 
      size = input$textSizeGSEA / 3, 
      max.overlaps = input$gseaMaxOver, 
      nudge_x = input$gseaDodgeX / 10, 
      nudge_y = input$gseaDodgeY / 10,
      fontface = "bold", 
      color = "black", 
      bg.color = "white", 
      bg.r = 0.15, 
      na.rm = TRUE
    )
  }
  
  # Add pathway labels
  cn <- cn + geom_label_repel(
    data = cd1, 
    aes(label = Label, x = x, y = y, fill = Label), 
    size = input$textSizeGSEA / 3, 
    max.overlaps = input$gseaMaxOver, 
    nudge_x = input$gseaDodgeX / 10, 
    nudge_y = input$gseaDodgeY / 10, 
    fontface = "bold", 
    show.legend = FALSE, 
    fill = alpha("white", input$gseaLabelAlpha), 
    na.rm = TRUE
  )
  
  gseaPlotReactive$cnetplot <- cn
  
  # Ridge plot
  rp <- enrichplot::ridgeplot(x = erRP, fill = input$ridgePlotColourByGSEA)
  rp <- rp +
    scale_y_discrete(labels = erRP@result$Label[match(
      levels(as.factor(rp$plot_env$gs2val.df$category)), 
      erRP@result$Description
    )]) +
    theme(axis.text.y = element_text(vjust = -0.01, size = input$textSizeGSEA)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_gradient2(
      low = "dodgerblue4", 
      high = "firebrick4", 
      mid = "white", 
      midpoint = 0
    )

  gseaPlotReactive$ridgeplot <- rp
  
  # Heatmap
  hp <- enrichplot::heatplot(
    x = current_er_symbol,
    foldChange = gseaToUse$log2RatioGSEA,
    showCategory = current_er_symbol@result$Description[input$selectedTable_GSEA_rows_selected]
  ) + 
    scale_colour_distiller(palette = "RdBu") + 
    coord_flip()
  
  gseaPlotReactive$heatmap <- hp
  
  # GSEA running score plot
  rs <- enrichplot::gseaplot2(
    x = current_er,
    geneSetID = input$selectedTable_GSEA_rows_selected,
    color = RColorBrewer::brewer.pal(length(input$selectedTable_GSEA_rows_selected), "Set2")
  )
  
  gseaPlotReactive$runningscore <- rs
  
  # UpSet plot
  geneSetsUpset <- setNames(
    lapply(input$selectedTable_GSEA_rows_selected, function(i) {
      colToUse <- if (input$gseaType %in% c("BP", "MF", "CC")) "geneName" else "core_enrichment"
      g <- current_er@result[[colToUse]][i] %>% 
        strsplit("/") %>% 
        .[[1]]
      g[nchar(g) > 1]
    }), 
    current_er@result[[labelColumn]][input$selectedTable_GSEA_rows_selected]
  )
  
  m1 <- ComplexHeatmap::make_comb_mat(geneSetsUpset)
  
  up <- ComplexHeatmap::UpSet(
    m = m1, 
    comb_order = order(comb_size(m1), decreasing = TRUE), 
    top_annotation = upset_top_annotation(
      m1, 
      add_numbers = TRUE, 
      annotation_name_rot = 90, 
      annotation_name_offset = unit(10, "mm")
    ), 
    left_annotation = upset_left_annotation(m1, add_numbers = TRUE)
  )
  
  gseaPlotReactive$upset <- up
})

# Render plot outputs ----

# Network plot 
output$cnetPlot_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["cnetplot"]]))
    gseaPlotReactive[["cnetplot"]]
  },
  width = as.numeric(input$figWidthGSEA),
  height = as.numeric(input$figHeightGSEA)
)

output$dlCnetPlot_GSEA <- downloadHandler(
  filename = function() {
    paste0(design, "_GSEA_", input$gseaType, "_network.pdf")
  },
  content = function(file) {
    pdf(
      file = file, 
      width = as.numeric(input$figWidthGSEA) / 50, 
      height = as.numeric(input$figHeightGSEA) / 50
    )
    print(gseaPlotReactive[["cnetplot"]])
    dev.off()
  }
)

# Heatmap
output$heatmap_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["heatmap"]]))
    gseaPlotReactive[["heatmap"]]
  },
  width = as.numeric(input$figWidthGSEA),
  height = as.numeric(input$figHeightGSEA * 4)
)

output$dlHeatmap_GSEA <- downloadHandler(
  filename = function() {
    paste0(design, "_GSEA_", input$gseaType, "_heatmap.pdf")
  },
  content = function(file) {
    pdf(
      file = file, 
      width = as.numeric(input$figWidthGSEA) / 25, 
      height = as.numeric(input$figHeightGSEA) / 50
    )
    print(gseaPlotReactive[["heatmap"]])
    dev.off()
  }
)

# Running score
output$runningScore_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["runningscore"]]))
    gseaPlotReactive[["runningscore"]]
  },
  width = as.numeric(input$figWidthGSEA),
  height = as.numeric(input$figHeightGSEA)
)

output$dlRunningScore_GSEA <- downloadHandler(
  filename = function() {
    paste0(design, "_GSEA_", input$gseaType, "_runningscore.pdf")
  },
  content = function(file) {
    pdf(
      file = file, 
      width = as.numeric(input$figWidthGSEA) / 50, 
      height = as.numeric(input$figHeightGSEA) / 50
    )
    print(gseaPlotReactive[["runningscore"]])
    dev.off()
  }
)

# Ridge plot 
output$ridgePlot_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["ridgeplot"]]))
    gseaPlotReactive[["ridgeplot"]]
  },
  width = as.numeric(input$figWidthGSEA),
  height = as.numeric(input$figHeightGSEA)
)

output$dlRidgePlot_GSEA <- downloadHandler(
  filename = function() {
    paste0(design, "_GSEA_", input$gseaType, "_ridgeplot.pdf")
  },
  content = function(file) {
    pdf(
      file = file, 
      width = as.numeric(input$figWidthGSEA) / 50, 
      height = as.numeric(input$figHeightGSEA) / 50
    )
    print(gseaPlotReactive[["ridgeplot"]])
    dev.off()
  }
)

# UpSet plot 
output$upsetPlot_GSEA <- renderPlot(
  {
    req(!is.null(gseaPlotReactive[["upset"]]))
    draw(gseaPlotReactive[["upset"]], padding = unit(c(5, 5, 5, 20), "mm"))
  },
  width = as.numeric(input$figWidthGSEA) * 1.3,
  height = as.numeric(input$figHeightGSEA)
)

output$dlUpset_GSEA <- downloadHandler(
  filename = function() {
    paste0(design, "_GSEA_", input$gseaType, "_upset.pdf")
  },
  content = function(file) {
    pdf(
      file = file, 
      width = as.numeric(input$figWidthGSEA) / 70, 
      height = as.numeric(input$figHeightGSEA) / 90
    )
    print(gseaPlotReactive[["upset"]])
    dev.off()
  }
)