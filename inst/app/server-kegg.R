if(inputDataReactive()$dataType == "RNASeq") {

  se <- inputDataReactive()$se
  param <- inputDataReactive()$param
  design <- inputDataReactive()$design
  keggPlotReactive <- reactiveValues(outfile = NULL, gD = NULL, keggInput = NULL, spp = NULL)
  output$keggPlotDesign <- renderText({design})
  
  spp <- ezRun::getSpecies(param$refBuild)
  if (grepl("^Rattus", param$refBuild)){
    spp <- "Rat"
  }
  
  observeEvent(input$tabs, {
    if (input$tabs == "keggTab") {
      cat(input$tabs)
      if (is_empty(allPathways)) {
        showModal(modalDialog(title = "No KEGG Data", "This dataset is not human or mouse, so this tab is not-functional."))
      }
    }
  })
  
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
        
        gD <- metadata(se)$enrichInput$log2Ratio[
          names(metadata(se)$enrichInput$log2Ratio) %in% inputDataReactive()$seqAnno$gene_id[
          which(inputDataReactive()$seqAnno[[pTypeKEGG]] <= as.numeric(input$pThresholdKEGG) & abs(inputDataReactive()$seqAnno$log2Ratio) >= as.numeric(input$lfcKEGG))
          ]
        ]
        
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
  list(src = keggPlotReactive$outfile, 
       contentType = 'image/png')
}, deleteFile = TRUE)