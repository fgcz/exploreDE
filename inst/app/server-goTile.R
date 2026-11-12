eds <- inputDataReactive()$.eds

# The KEGG tile's state. Defined for every dataset, so saving a session can
# read it without checking that it exists.
goTileKeggReactive <- reactiveValues(
  keggCode = NULL,
  cache = list(),
  currentGenes = NULL,
  currentPathwayName = NULL,
  currentKeggId = NULL
)

if (eds$data_type() == "RNASeq" && eds$has_de()) {
  allPathways <- eds$all_pathways()
  # GO behaves like KEGG: The GO selector shows only when the features carry
  # GO annotation columns (any of GO BP / MF / CC).
  output$goTileGoSupported <- reactive({
    length(allPathways) > 0
  })
  outputOptions(output, "goTileGoSupported", suspendWhenHidden = FALSE)

  # Follows the app-wide selected contrast (constant for single-contrast data).
  output$goTilePlotDesign <- renderText({
    eds$contrast_label(selectedContrast())
  })

  # --- KEGG tile support (mirrors the heatmap KEGG option) -----------------
  # A fetched pathway is rendered as an additional tile alongside any GO terms.

  observe({
    # Detect species -> KEGG organism code (only hsa/mmu/rno supported)
    goTileKeggReactive$keggCode <- if (kegg_installed()) kegg_org_code(eds$ref_build())
    output$goTileKeggSupported <- reactive({
      !is.null(goTileKeggReactive$keggCode)
    })
    outputOptions(output, "goTileKeggSupported", suspendWhenHidden = FALSE)
  })

  observeEvent(input$tabs, {
    if (input$tabs == "goTileTab") {
      # Only warn if there is neither GO nor KEGG available for this dataset.
      if (length(allPathways) == 0 && is.null(goTileKeggReactive$keggCode)) {
        showModal(modalDialog(
          title = "No GO Annotations",
          "Your genes have no GO annotations, so this tab is not functional for this dataset."
        ))
      }
    }
  })

  if (length(allPathways) > 0) {
    updateSelectizeInput(
      session = session,
      inputId = "goTileInput",
      choices = allPathways,
      selected = "",
      server = TRUE
    )
  }

  # KEGG tile: fetch pathway genes on button click (mirrors server-heatmap2.R)
  observeEvent(input$goTileKeggGo, {
    req(input$goTileKeggInput)
    req(goTileKeggReactive$keggCode)

    keggCode <- goTileKeggReactive$keggCode
    pathwayInput <- trimws(input$goTileKeggInput)
    pathwayId <- if (grepl("^[0-9]", pathwayInput)) paste0(keggCode, pathwayInput) else pathwayInput

    # Serve from cache if we've fetched this pathway already this session
    if (pathwayId %in% names(goTileKeggReactive$cache)) {
      cached <- goTileKeggReactive$cache[[pathwayId]]
      goTileKeggReactive$currentGenes <- cached$genes
      goTileKeggReactive$currentPathwayName <- cached$name
      goTileKeggReactive$currentKeggId <- pathwayId
      return()
    }

    fetched <- withProgress(message = "Fetching KEGG pathway genes...", {
      tryCatch(
        kegg_fetch_pathway(pathwayId, keggCode),
        error = function(e) {
          showNotification(paste("Failed to fetch KEGG pathway:", e$message), type = "error", duration = 8)
          NULL
        }
      )
    })
    result <- if (!is.null(fetched)) {
      list(genes = kegg_map_genes(eds$de(selectedContrast()), fetched$entrez, keggCode), name = fetched$name)
    }

    if (!is.null(result)) {
      goTileKeggReactive$cache[[pathwayId]] <- result
      goTileKeggReactive$currentGenes <- result$genes
      goTileKeggReactive$currentPathwayName <- result$name
      goTileKeggReactive$currentKeggId <- pathwayId
    } else {
      goTileKeggReactive$currentGenes <- NULL
      goTileKeggReactive$currentPathwayName <- NULL
      goTileKeggReactive$currentKeggId <- NULL
    }
  })

  output$goTileKeggInfo <- renderUI({
    req(goTileKeggReactive$currentPathwayName)
    tagList(
      hr(),
      tags$b("Pathway: "),
      goTileKeggReactive$currentPathwayName,
      br(),
      tags$b("Genes found in dataset: "),
      length(goTileKeggReactive$currentGenes)
    )
  })

  # --- Tile plot builder (GO terms + optional KEGG pathway) ----------------
  observeEvent(
    {
      bananaMode() # rebuild with bananas / asterisks when banana mode toggles
      input$figHeightGoTile
      input$figWidthGoTile
      input$textSizeGoTile
      input$goTileInput
      input$pTypeGoTile
      input$clusterRowsGoTile
      input$colnamesGoTile
      input$geneNamesGoTile
      input$GoTileColourRed
      input$GoTileColourWhite
      input$GoTileColourBlue
      input$tileLimitColour
      input$degOnlyGoTile
      input$lfcGoTile
      input$pThresholdGoTile
      goTileKeggReactive$currentGenes
      selectedContrast()
    },
    # ignoreNULL = FALSE: goTileKeggReactive$currentGenes is NULL until a KEGG
    # pathway is fetched; if it ever ends the block again, the default
    # ignoreNULL = TRUE would stop GO terms alone from drawing a tile.
    ignoreNULL = FALSE,
    ignoreInit = TRUE,
    {
      figuresDataReactive$goTilePlot <- NULL

      pCol <- if (input$pTypeGoTile == "Raw") "pValue" else "fdr"
      seqAnno <- eds$de(selectedContrast())

      # Significance filter shared by GO and KEGG tiles (mirrors heatmap):
      # keep features used in the test, and — when requested — those passing
      # both the p-value threshold and the absolute log2FC threshold.
      applySig <- function(d1) {
        # A tile is coloured by log2FC, so rows without one can't be drawn (and
        # NA rows break the row clustering).
        d1 <- d1[d1$usedInTest %in% TRUE & !is.na(d1$log2Ratio), ]
        if (isTRUE(input$degOnlyGoTile)) {
          d1 <- d1[
            which(
              d1[[pCol]] <= as.numeric(input$pThresholdGoTile %||% 1) &
                abs(d1$log2Ratio) >= as.numeric(input$lfcGoTile %||% 0)
            ),
          ]
        }
        d1
      }

      # GO-term sources (one tile each), keyed by GO ID
      goToPlot <- gsub(" .*", "", input$goTileInput)
      goToPlot <- goToPlot[nzchar(goToPlot)]
      dfToPlotList <- list()
      # GO membership comes from features(), which always carries the GO
      # columns; .h5ad DE tables don't. Same feature_id space on every backend.
      goFeats <- eds$features()
      for (GO in goToPlot) {
        d1 <- seqAnno[seqAnno$feature_id %in% rownames(goFeats)[features_in_go_term(goFeats, GO)], ]
        dfToPlotList[[GO]] <- applySig(d1)
      }

      # KEGG source (if a pathway has been fetched), keyed by pathway name
      if (!is.null(goTileKeggReactive$currentGenes) && length(goTileKeggReactive$currentGenes) > 0) {
        keggTitle <- goTileKeggReactive$currentPathwayName %||% goTileKeggReactive$currentKeggId
        d1 <- seqAnno[seqAnno$feature_id %in% goTileKeggReactive$currentGenes, ]
        dfToPlotList[[keggTitle]] <- applySig(d1)
      }

      # Drop any source that has no features left after filtering
      dfToPlotList <- dfToPlotList[vapply(dfToPlotList, function(d) nrow(d) > 0, logical(1))]
      req(length(dfToPlotList) > 0)
      titles <- names(dfToPlotList)

      cols <- colorRamp2(
        breaks = c(input$tileLimitColour, 0, -input$tileLimitColour),
        colors = c(
          paste0(input$GoTileColourRed, "FF"),
          paste0(input$GoTileColourWhite, "FF"),
          paste0(input$GoTileColourBlue, "FF")
        )
      )

      bnGoTile <- banana_on(bananaMode())

      hmList <- lapply(titles, function(ttl) {
        d1 <- dfToPlotList[[ttl]]
        m1 <- as.matrix(d1$log2Ratio)
        colnames(m1) <- ttl
        rownames(m1) <- d1$feature_id
        m2 <- as.matrix(d1[[pCol]])
        colnames(m2) <- ttl
        rownames(m2) <- d1$feature_id

        Heatmap(
          matrix = m1,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if (bnGoTile) {
              banana_grid(banana_tier(m2[i, j]), x, y, h = grid::unit(input$textSizeGoTile * 3, "pt"))
              return(invisible(NULL))
            }
            # No stars for an untested p-value (e.g. DESeq2's independent
            # filtering leaves padj NA).
            grid.text(significance_stars(m2[i, j]), x, y)
          },
          name = "Log2FC",
          col = cols,
          cluster_rows = input$clusterRowsGoTile,
          show_column_names = input$colnamesGoTile,
          show_row_names = input$geneNamesGoTile,
          column_names_centered = TRUE,
          column_title = sub("\\s+", "\n", ttl),
          column_title_gp = gpar(fontsize = ceiling(input$textSizeGoTile * 1.1), fontface = "bold"),
          row_names_gp = gpar(fontsize = input$textSizeGoTile),
          column_names_gp = gpar(fontsize = input$textSizeGoTile),
          heatmap_legend_param = list(
            title = "Log2FC",
            at = c(-input$tileLimitColour, 0, input$tileLimitColour),
            labels = c(-input$tileLimitColour, "0", input$tileLimitColour)
          ),
          width = ncol(m1) * unit(15, "mm"),
          height = nrow(m1) * unit(5, "mm"),
          rect_gp = gpar(col = "white", lwd = 0.5)
        )
      }) %>%
        setNames(titles)

      figuresDataReactive$goTilePlot <- hmList
    }
  )
}

output$goTileOutput <- renderPlot(
  {
    req(!is.null(figuresDataReactive$goTilePlot))
    req(length(figuresDataReactive$goTilePlot) > 0)

    grobList <- lapply(figuresDataReactive$goTilePlot, function(hm) {
      grid.grabExpr(draw_heatmap(hm, legend_dx = input$goTileLegendDX, legend_dy = input$goTileLegendDY))
    })

    cowplot::plot_grid(plotlist = grobList, ncol = length(grobList))
  },
  width = function() {
    as.numeric(input$figWidthGoTile)
  },
  height = function() {
    as.numeric(input$figHeightGoTile)
  }
)

output$dlGoTilePlot <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = paste0(input$filename_goTile, "_", gsub(" .*", "", input$goTileInput)),
      design_name = if (eds$has_de()) eds$contrast_label(selectedContrast()) else "",
      extension = tolower(input$downloadFormatGoTile)
    )
  },
  content = function(file) {
    req(!is.null(figuresDataReactive$goTilePlot))
    req(length(figuresDataReactive$goTilePlot) > 0)

    export_plot_standard(
      plot_object = function() {
        grobList <- lapply(figuresDataReactive$goTilePlot, function(hm) {
          grid.grabExpr(draw_heatmap(hm, legend_dx = input$goTileLegendDX, legend_dy = input$goTileLegendDY))
        })
        print(cowplot::plot_grid(plotlist = grobList, ncol = length(grobList)))
      },
      filename = file,
      format = tolower(input$downloadFormatGoTile),
      width = as.numeric(input$figWidthGoTile / 80),
      height = as.numeric(input$figHeightGoTile / 60),
      dpi = as.numeric(input$dpiGoTile)
    )
  }
)
