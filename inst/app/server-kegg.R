eds <- inputDataReactive()$.eds

if (eds$data_type() == "RNASeq") {
  design <- eds$contrast_label(eds$primary_contrast())
  keggPlotReactive <- reactiveValues(outfile = NULL, gD = NULL, keggInput = NULL, spp = NULL)
  output$keggPlotDesign <- renderText({ design })

  ref <- eds$ref_build()
  spp <- ezRun::getSpecies(ref)
  if (grepl("^Rattus", ref)) spp <- "Rat"
  
  observeEvent(
    {
      input$keggInput
      input$lfcKEGG
      input$pTypeKEGG
      input$pThresholdKEGG
      input$keggLimitColour
    }, ignoreInit = TRUE,
    {
      # tryCatch({
      if (input$pTypeKEGG == "Raw") {
        pTypeKEGG <- "pValue"
      } else {
        pTypeKEGG <- "fdr"
      }
      
      switch(
        spp,
        "Human" = {spp = "hsa"},
        "Mouse" = {spp = "mmu"},
        "Rat" = {spp = "rno"}
      )
      
      if (grepl("^[0-9]", input$keggInput)) {
        keggInput <- paste0(spp, input$keggInput)
      } else {
        keggInput <- input$keggInput
      }
      
      ei     <- eds$enrichment_input(eds$primary_contrast())
      de_df  <- eds$de(eds$primary_contrast())
      sig_ids <- de_df$gene_id[
        which(de_df[[pTypeKEGG]] <= as.numeric(input$pThresholdKEGG) &
              abs(de_df$log2Ratio) >= as.numeric(input$lfcKEGG))
      ]
      gD <- ei$ranked_l2fc[names(ei$ranked_l2fc) %in% sig_ids]
      
      outfile <- paste0(keggInput, ".pathview.png")
      keggPlotReactive$outfile <- outfile
      keggPlotReactive$gD = gD
      keggPlotReactive$keggInput = keggInput
      keggPlotReactive$spp = spp
      
      # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  )
} else {
  observeEvent(input$tabs, {
    if (input$tabs == "keggTab") {
      showModal(modalDialog(title = "No KEGG Results", "Unfortunately, our proteomics results don't include any pathway analyses! We hope to change this very soon."))
    }
  })
}

output$keggOutput <- renderImage({
  req(!is.null(keggPlotReactive$outfile))
  
  # Use temporary directory to avoid Docker permissions issues
  temp_dir <- tempdir()
  old_wd <- getwd()
  
  tryCatch({
    # Change to temp directory where we have write permissions
    setwd(temp_dir)
    
    pathview(
      gene.data = keggPlotReactive$gD,
      gene.idtype = "ENSEMBL",
      pathway.id = keggPlotReactive$keggInput,
      key.pos = "bottomright",
      species = keggPlotReactive$spp, 
      low = list("gene" = "dodgerblue3"),
      mid = list("gene" = "lightgrey"),
      high = list("gene" = "indianred3"),
      kegg.native = TRUE,
      limit = list(gene = input$keggLimitColour))
    
    # Get the full path to the generated file
    output_file <- file.path(temp_dir, keggPlotReactive$outfile)
    
    list(src = output_file, 
         contentType = 'image/png')
    
  }, error = function(e) {
    cat("KEGG plot error:", conditionMessage(e), "\n")
    NULL
  }, finally = {
    # Always restore working directory
    setwd(old_wd)
  })
}, deleteFile = TRUE)

output$dlPlot_kegg <- downloadHandler(
  filename = function() {
    generate_filename(
      base_name = input$filename_kegg,
      suffix = keggPlotReactive$keggInput,
      design_name = design,
      extension = "png"
    )
  },
  content = function(file) {
    req(!is.null(keggPlotReactive$keggInput))
    
    # Use temporary directory
    temp_dir <- tempdir()
    old_wd <- getwd()
    
    tryCatch({
      setwd(temp_dir)
      
      # Run pathview to generate the file
      # Note: pathview writes the file to the current working directory
      pv.out <- pathview(
        gene.data = keggPlotReactive$gD,
        gene.idtype = "ENSEMBL",
        pathway.id = keggPlotReactive$keggInput,
        key.pos = "bottomright",
        species = keggPlotReactive$spp, 
        low = list("gene" = "dodgerblue3"),
        mid = list("gene" = "lightgrey"),
        high = list("gene" = "indianred3"),
        kegg.native = TRUE,
        limit = list(gene = input$keggLimitColour)
      )
      
      # Determine the generated filename
      # pathview defaults to [pathway.id].pathview.png
      generated_file <- paste0(keggPlotReactive$keggInput, ".pathview.png")
      
      if (file.exists(generated_file)) {
        file.copy(generated_file, file)
      } else {
        stop("Pathview failed to generate the output file.")
      }
      
    }, error = function(e) {
      cat("KEGG download error:", conditionMessage(e), "\n")
    }, finally = {
      setwd(old_wd)
    })
  }
)