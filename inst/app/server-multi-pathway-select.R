# Multi-Contrast Manual Pathway Selection
# Hand-pick GO terms (from the ORA-tested set, independent of significance) and compare
# them across contrasts with two views, both driven by the same selection:
#   1. lollipop plot   (x = -log10 p-value, one lollipop per contrast)
#   2. network         (enrichplot::cnetplot on a merged compareClusterResult)
# plus a shared stats table. Per-pathway stats and network membership come from the ORA
# enrichResults, matched by GO ID; a term only appears for a contrast where it was tested.

pwPlots <- reactiveValues(lollipop = NULL, cnet = NULL, cnetMsg = NULL, table = NULL)

# Distinct per-contrast colours: user's colour pickers if set, else the catalyst
# palette (mirrors server-multi-tile.R) so contrasts differ before the colour tab.
pwContrastColours <- function() {
  results <- multiDeReactive$results
  catalystCols <- c(
    "#DC050C",
    "#FB8072",
    "#1965B0",
    "#7BAFDE",
    "#882E72",
    "#B17BA6",
    "#FF7B00",
    "#FDC362",
    "#E7298A",
    "#E78AC3",
    "#33A02C",
    "#B2DF8A",
    "#55A1B1",
    "#8DD3C7",
    "#A6761D",
    "#E6AB02",
    "#7570B3",
    "#BEAED4",
    "#666666",
    "#999999",
    "#aa8282",
    "#d4b7b7",
    "#8600bf",
    "#ba5ce3",
    "#808000",
    "#aeae5c",
    "#1e90ff",
    "#00bfff",
    "#56ff0d",
    "#ffff00"
  )
  catalystCols <- c(catalystCols[c(seq(1, 30, by = 2), seq(2, 30, by = 2))])
  cols <- setNames(
    lapply(seq_along(names(results)), function(i) {
      k_clean <- gsub(" ", "_", names(results)[i])
      col <- input[[paste0("multiGroupColour", k_clean)]]
      if (is.null(col)) {
        return(paste0(rep(catalystCols, times = 5)[i], "FF"))
      }
      col <- tryCatch(gplots::col2hex(col), error = function(e) "#000000")
      paste0(col, "FF")
    }),
    names(results)
  )
  unlist(cols)
}

# Named-ID choices for the pathway selectize: every GO term tested in the ORA for the
# chosen direction (all of BP/MF/CC), independent of significance. This is the set that
# actually carries stats to plot; the selectize *value* is the bare GO ID so matching to
# the ORA results is exact, while the label shows the readable description.
pwChoicesReactive <- reactive({
  req(multiDeReactive$results)
  req(input$pwDirection)
  results <- multiDeReactive$results
  dir <- input$pwDirection
  onts <- c("BP", "MF", "CC")

  idDesc <- do.call(
    rbind,
    lapply(names(results), function(con) {
      ora <- results[[con]]$ora
      if (is.null(ora)) {
        return(NULL)
      }
      do.call(
        rbind,
        lapply(onts, function(ont) {
          er <- ora[[ont]][[dir]]
          if (is.null(er) || nrow(er@result) < 1) {
            return(NULL)
          }
          unique(er@result[, c("ID", "Description")])
        })
      )
    })
  )
  req(!is.null(idDesc) && nrow(idDesc) > 0)

  idDesc <- unique(idDesc)
  choices <- setNames(idDesc$ID, paste0(idDesc$ID, "  ", idDesc$Description))
  choices[order(names(choices))]
})

# Populate the pathway selectize. Preserve the current selection across refreshes
# (e.g. when contrasts are renamed, which mutates results but not the GO IDs); only
# terms that are no longer valid choices are dropped.
observe({
  choices <- pwChoicesReactive()
  keep <- isolate(input$pathwaySelectInput)
  keep <- keep[keep %in% choices]
  updateSelectizeInput(session, "pathwaySelectInput", choices = choices, selected = keep, server = TRUE)
})

# Populate the gene-search selectize with every feature detected in the data.
observe({
  req(multiDeReactive$results)
  genes <- sort(unique(unlist(lapply(multiDeReactive$results, function(r) r$seqAnno$feature_id))))
  req(length(genes) > 0)
  updateSelectizeInput(session, "pwGeneSearch", choices = genes, selected = character(0), server = TRUE)
})

# Find pathways by gene: add every tested pathway (chosen direction) whose gene set
# contains any of the selected genes to the current pathway selection.
observeEvent(input$pwGeneFindGo, {
  req(length(input$pwGeneSearch) >= 1)
  req(multiDeReactive$results)
  results <- multiDeReactive$results
  dir <- input$pwDirection
  onts <- c("BP", "MF", "CC")
  genes <- input$pwGeneSearch

  matchIds <- character(0)
  for (con in names(results)) {
    ora <- results[[con]]$ora
    if (is.null(ora)) {
      next
    }
    seqAnno <- results[[con]]$seqAnno
    id2name <- if (!is.null(seqAnno) && all(c("gene_id", "feature_id") %in% colnames(seqAnno))) {
      setNames(seqAnno$feature_id, seqAnno$gene_id)
    } else {
      NULL
    }
    for (ont in onts) {
      er <- ora[[ont]][[dir]]
      if (is.null(er) || nrow(er@result) < 1) {
        next
      }
      res <- er@result
      memberSymbols <- lapply(strsplit(res$geneID, "/"), function(ids) {
        if (!is.null(id2name)) {
          nm <- id2name[ids]
          ifelse(is.na(nm) | nchar(nm) == 0, ids, nm)
        } else {
          ids
        }
      })
      hit <- vapply(memberSymbols, function(ms) any(genes %in% ms), logical(1))
      matchIds <- c(matchIds, res$ID[hit])
    }
  }
  matchIds <- unique(matchIds)

  if (length(matchIds) == 0) {
    showNotification(
      "No tested pathways contain the selected gene(s) in this direction.",
      type = "warning",
      duration = 6
    )
    return()
  }

  newSel <- unique(c(input$pathwaySelectInput, matchIds))
  updateSelectizeInput(session, "pathwaySelectInput", choices = pwChoicesReactive(), selected = newSel, server = TRUE)
  showNotification(
    paste0("Added ", length(matchIds), " pathway(s) containing the selected gene(s)."),
    type = "message",
    duration = 5
  )
})

# Shared long-format table: one row per (contrast, selected term) that has an ORA
# result in the chosen direction. Matched by GO ID across BP/MF/CC. Consumed by the
# bubble, the lollipop, and the table.
pwDataReactive <- reactive({
  req(multiDeReactive$results)
  req(input$pwDirection)
  req(length(input$pathwaySelectInput) >= 1)
  results <- multiDeReactive$results
  dir <- input$pwDirection
  # The selectize value is already the bare GO ID.
  selIds <- input$pathwaySelectInput
  onts <- c("BP", "MF", "CC")

  dfList <- lapply(names(results), function(con) {
    ora <- results[[con]]$ora
    if (is.null(ora)) {
      return(NULL)
    }
    rows <- lapply(onts, function(ont) {
      er <- ora[[ont]][[dir]]
      if (is.null(er)) {
        return(NULL)
      }
      hit <- er@result[er@result$ID %in% selIds, , drop = FALSE]
      if (nrow(hit) == 0) {
        return(NULL)
      }
      data.frame(
        Contrast = con,
        ID = hit$ID,
        Description = hit$Description,
        GeneRatio = hit$GeneRatio,
        Count = hit$Count,
        pvalue = hit$pvalue,
        p.adjust = hit$p.adjust,
        stringsAsFactors = FALSE
      )
    })
    rows <- rows[!sapply(rows, is.null)]
    if (length(rows) == 0) {
      return(NULL)
    }
    do.call(rbind, rows)
  })
  dfList <- dfList[!sapply(dfList, is.null)]
  req(length(dfList) >= 1)

  df <- do.call(rbind, dfList)
  df$Contrast <- factor(df$Contrast, levels = names(results))
  df$negLogPadj <- -log10(df$p.adjust)
  df
})

# Keep the stats table in sync with the current selection.
observeEvent(pwDataReactive(), {
  pwPlots$table <- pwDataReactive()
})

# When the selection is emptied, clear the lollipop and table too (pwDataReactive
# req-fails on an empty selection, so the builders above won't fire to clear them).
observeEvent(
  input$pathwaySelectInput,
  {
    if (length(input$pathwaySelectInput) < 1) {
      pwPlots$lollipop <- NULL
      pwPlots$table <- NULL
    }
  },
  ignoreNULL = FALSE
)

# 2. Lollipop plot ----
# Straight horizontal stems (y == yend). Contrasts are separated with a manual y
# offset rather than position_dodge, which slants geom_segment.
observeEvent(
  {
    pwDataReactive()
    input$pwPType
    input$pwLollipopScale
    input$textSizePwSelect
    input$figWidthPwSelect
    input$figHeightPwSelect
    lapply(seq_along(multiDeReactive$results), function(i) {
      input[[(paste0("multiGroupColour", names(multiDeReactive$results)[[i]]) %>% gsub(" ", "_", .))]]
    })
  },
  ignoreNULL = TRUE,
  ignoreInit = TRUE,
  {
    df <- pwDataReactive()
    cols <- pwContrastColours()
    ts <- input$textSizePwSelect
    scl <- input$pwLollipopScale %||% 1
    pcol <- if (identical(input$pwPType, "Raw")) "pvalue" else "p.adjust"
    plab <- if (identical(input$pwPType, "Raw")) "-log10 raw p-value" else "-log10 adjusted p-value"
    df$negLogP <- -log10(df[[pcol]])

    descs <- sort(unique(as.character(df$Description)))
    df$ybase <- as.integer(factor(as.character(df$Description), levels = descs))
    nCon <- length(levels(df$Contrast))
    # Spread the contrasts within a tight band centred on each pathway row.
    df$ypos <- if (nCon > 1) {
      df$ybase + ((as.integer(df$Contrast) - (nCon + 1) / 2) / nCon) * 0.35
    } else {
      df$ybase
    }

    p <- ggplot(df, aes(colour = Contrast)) +
      geom_segment(aes(x = 0, xend = negLogP, y = ypos, yend = ypos), linewidth = scl) +
      geom_point(aes(x = negLogP, y = ypos, size = Count)) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
      scale_y_continuous(breaks = seq_along(descs), labels = descs, expand = expansion(add = 0.6)) +
      scale_colour_manual(values = cols) +
      scale_size(range = c(ts / 3, ts) * scl, name = "Gene Count") +
      xlab(plab) +
      ylab(NULL) +
      theme_classic() +
      theme(
        axis.text.x = element_text(size = ts, face = "bold"),
        axis.text.y = element_text(size = ts, face = "bold"),
        axis.title.x = element_text(size = ts, face = "bold"),
        legend.title = element_text(size = ts, face = "bold"),
        legend.text = element_text(size = ts, face = "bold")
      )

    pwPlots$lollipop <- p
  }
)

# 3. Network (cnetplot on a merged compareClusterResult) ----
observeEvent(
  {
    input$pathwaySelectInput
    input$pwDirection
    input$pwCnetGeneLabels
    input$pwCnetCatScale
    input$pwCnetGeneScale
    input$textSizePwSelect
    lapply(seq_along(multiDeReactive$results), function(i) {
      input[[(paste0("multiGroupColour", names(multiDeReactive$results)[[i]]) %>% gsub(" ", "_", .))]]
    })
  },
  ignoreNULL = FALSE,
  ignoreInit = TRUE,
  {
    req(multiDeReactive$results)
    results <- multiDeReactive$results
    dir <- input$pwDirection
    selIds <- input$pathwaySelectInput
    onts <- c("BP", "MF", "CC")

    if (length(selIds) < 1) {
      pwPlots$cnet <- NULL
      pwPlots$cnetMsg <- "Select one or more GO terms to build the network."
      return()
    }

    # One filtered enrichResult per contrast, restricted to the selected GO IDs.
    # We collapse BP/MF/CC into a single enrichResult per contrast so merge_result
    # yields one cluster (pie slice) per contrast.
    erList <- list()
    for (con in names(results)) {
      ora <- results[[con]]$ora
      if (is.null(ora)) {
        next
      }

      baseEr <- NULL
      hitRows <- list()
      for (ont in onts) {
        er <- ora[[ont]][[dir]]
        if (is.null(er)) {
          next
        }
        if (is.null(baseEr)) {
          baseEr <- er
        }
        hit <- er@result[er@result$ID %in% selIds, , drop = FALSE]
        if (nrow(hit) > 0) hitRows[[ont]] <- hit
      }
      if (is.null(baseEr) || length(hitRows) == 0) {
        next
      }

      merged <- do.call(rbind, hitRows)

      # Relabel geneID (gene_id) -> feature_id (display symbol) for readable nodes,
      # mirroring server-hypergeo2.R. No-op where feature_id already == geneID.
      seqAnno <- results[[con]]$seqAnno
      if (!is.null(seqAnno) && all(c("gene_id", "feature_id") %in% colnames(seqAnno))) {
        id2name <- setNames(seqAnno$feature_id, seqAnno$gene_id)
        merged$geneID <- sapply(strsplit(merged$geneID, "/"), function(ids) {
          nm <- id2name[ids]
          nm <- ifelse(is.na(nm) | nchar(nm) == 0, ids, nm)
          paste(nm, collapse = "/")
        })
      }

      erOut <- baseEr
      erOut@result <- merged
      erList[[con]] <- erOut
    }

    if (length(erList) == 0) {
      pwPlots$cnet <- NULL
      pwPlots$cnetMsg <- "None of the selected terms were tested in the ORA results for this direction."
      return()
    }

    catScale <- input$pwCnetCatScale %||% 1
    geneScale <- input$pwCnetGeneScale %||% 1
    cols <- pwContrastColours()
    cn <- tryCatch(
      {
        cc <- clusterProfiler::merge_result(erList)
        p <- enrichplot::cnetplot(
          cc,
          showCategory = length(selIds),
          node_label = if (isTRUE(input$pwCnetGeneLabels)) "all" else "category",
          size_category = catScale,
          size_item = geneScale
        )
        # Recolour the pathway-node pies with the app's contrast colours
        # (cnetplot otherwise uses ggplot's default hue palette for Cluster).
        p <- p + scale_fill_manual(values = cols, name = "Contrast")
        p +
          theme(
            legend.title = element_text(size = input$textSizePwSelect, face = "bold"),
            legend.text = element_text(size = input$textSizePwSelect)
          )
      },
      error = function(e) {
        showNotification(paste("Could not build the network plot:", e$message), type = "error", duration = 8)
        NULL
      }
    )
    pwPlots$cnet <- cn
    pwPlots$cnetMsg <- if (is.null(cn)) "Could not build the network plot from the current selection." else NULL
  }
)

# Outputs ----
output$pwLollipopPlot <- renderPlot(
  {
    validate(need(!is.null(pwPlots$lollipop), "Select one or more GO terms to plot."))
    pwPlots$lollipop
  },
  width = function() as.numeric(input$figWidthPwSelect),
  height = function() as.numeric(input$figHeightPwSelect)
)

output$pwCnetPlot <- renderPlot(
  {
    validate(need(!is.null(pwPlots$cnet), pwPlots$cnetMsg %||% "Select one or more GO terms to build the network."))
    pwPlots$cnet
  },
  width = function() as.numeric(input$figWidthPwSelect),
  height = function() as.numeric(input$figHeightPwSelect)
)

output$pwBubbleTable <- renderDataTable({
  req(pwPlots$table)
  DT::datatable(
    data = pwPlots$table,
    filter = "top",
    class = "cell-border stripe",
    rownames = FALSE,
    options = list(scrollX = TRUE, autoWidth = TRUE)
  ) %>%
    DT::formatSignif(columns = c("pvalue", "p.adjust", "negLogPadj"), digits = 4) %>%
    DT::formatStyle(
      columns = c("pvalue", "p.adjust"),
      color = styleInterval(cuts = 0.05, values = c("green", "black")),
      fontWeight = "bold"
    )
})

# Downloads ----
output$dlPlot_pwLollipop <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = input$filename_pwLollipop,
      design_name = "Multi",
      extension = tolower(input$downloadFormatPwSelect)
    )
  },
  content = function(file) {
    req(pwPlots$lollipop)
    export_plot_standard(
      plot_object = pwPlots$lollipop,
      filename = file,
      format = tolower(input$downloadFormatPwSelect),
      width = as.numeric(input$figWidthPwSelect) / 72,
      height = as.numeric(input$figHeightPwSelect) / 72,
      dpi = as.numeric(input$dpiPwSelect)
    )
  }
)

output$dlPlot_pwCnet <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = input$filename_pwCnet,
      design_name = "Multi",
      extension = tolower(input$downloadFormatPwSelect)
    )
  },
  content = function(file) {
    req(pwPlots$cnet)
    export_plot_standard(
      plot_object = pwPlots$cnet,
      filename = file,
      format = tolower(input$downloadFormatPwSelect),
      width = as.numeric(input$figWidthPwSelect) / 72,
      height = as.numeric(input$figHeightPwSelect) / 72,
      dpi = as.numeric(input$dpiPwSelect)
    )
  }
)

output$dlTable_pwBubble <- downloadHandler(
  filename = function() {
    generate_filename(base_name = "PathwaySelect_Table", design_name = "Multi", extension = "xlsx")
  },
  content = function(file) {
    req(pwPlots$table)
    export_excel_standard(
      df = pwPlots$table,
      filename = file,
      sheet_name = "Pathway Stats",
      table_name = "PathwayStats"
    )
  }
)
