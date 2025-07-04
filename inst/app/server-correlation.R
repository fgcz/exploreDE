output$correlationStatic <- renderPlot({NULL})

# Update some inputs 
updateSelectizeInput(session = session, inputId = "correlationGene1", label = "Feature 1", choices = inputDataReactive()$genes, selected = "", server = TRUE)
updateSelectizeInput(session = session, inputId = "correlationGene2", label = "Feature 1", choices = inputDataReactive()$genes, selected = "", server = TRUE)
updateSelectInput(session = session, inputId = "correlationColourBy", label = "Select column to colour by", choices = inputDataReactive()$factors, selected =  inputDataReactive()$factors[[1]])
updateSelectInput(session = session, inputId = "correlationShapeBy", label = "Select column to shape by", choices = c("None", inputDataReactive()$factors), selected =  "None")
updateSelectInput(session = session, inputId = "correlationCounts", label = "Select count method to plot", choices = names(inputDataReactive()$countList), selected = switch(inputDataReactive()$dataType, "RNASeq" = "Normalised + Log2", "proteomics" = "transformedData"))
updateSelectInput(session = session, inputId = "correlationBatch", label = "Apply a batch-correction to the counts?", choices = c("None", inputDataReactive()$factors), selected = "None")
updateSliderInput(session = session, inputId = "correlationCorPosX", label = "Label position x", min = 1, max = 100, value = 10, step = 1)
updateSliderInput(session = session, inputId = "correlationCorPosY", label = "Label position y", min = 1, max = 100, value = 10, step = 1)

observeEvent({
  input$correlationColourBy
}, {
  output$correlationGroupBucket <- renderUI({
    req(input$tabs == "correlationsTab")
    bucket_list(
      header = "Drag and drop groups in order to be plotted",
      group_name = "correlation_bucket_list",
      orientation = "horizontal",
      add_rank_list(
        text = "Include these groups",
        labels = as.list(unique(levels(as.factor(inputDataReactive()$dataset[[input$correlationColourBy]])))),
        input_id = "correlationKeepBucket"),
      add_rank_list(
        text = "Exclude these groups",
        labels = NULL,
        input_id = "correlationExcludeBucket")
    )  
  })
  outputOptions(output, "correlationGroupBucket", suspendWhenHidden = FALSE)
})

correlationPlot <- reactiveValues(plot = NULL)
observeEvent({
  input$correlationMethod
  input$correlationGene1
  input$correlationGene2
  input$correlationColourBy
  input$correlationShapeBy
  input$correlationCounts
  input$correlationBatch
  input$correlationCorPosX
  input$correlationCorPosY
  input$correlationPointSize
  input$correlationPointStroke
  input$correlationPointAlpha
  input$textSizeCorrelation
  input$correlationKeepBucket
  input$correlationCountsLog
  input$correlationShowNames
  input$correlationMaxOverlaps
  input$correlationLabelPosX
  input$correlationLabelPosY
  lapply(seq_along(inputDataReactive()$factorLevels), function (i) {
    input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]])]]
  })
}, ignoreNULL = TRUE, ignoreInit = TRUE, {
  
  req(input$correlationGene1 != input$correlationGene2)
  req(nchar(input$correlationGene1) > 1 & nchar(input$correlationGene2) > 1)
  
  # get the counts 
  countsCorrelation <- inputDataReactive()$countList[[input$correlationCounts]][c(input$correlationGene1, input$correlationGene2),]
  
  # log2 if requested/necessary
  if (input$correlationCountsLog & !input$correlationCounts %in% c("Normalised + Log2", "VST", "transformedData")) {
    countsCorrelation <- log2(countsCorrelation+1)
  }
  
  # remove batch effect if requested 
  if (input$correlationBatch != "None") {
    countsCorrelation <- limma::removeBatchEffect(x = countsCorrelation, batch = inputDataReactive()$dataset[[input$correlationBatch]])
  }
  
  # prepare the counts 
  countsCorrelation <- countsCorrelation %>%
    t %>%
    data.frame(check.names = F) %>% 
    rownames_to_column("Name") %>% 
    left_join((inputDataReactive()$dataset %>% data.frame(check.names = F) %>% rownames_to_column("Name")), by = "Name")
  countsCorrelation <- countsCorrelation[which(countsCorrelation[[input$correlationColourBy]] %in% input$correlationKeepBucket),]
  
  # get the colours 
  coloursCorrelation <- setNames(lapply(input$correlationKeepBucket, function(k) {
    paste0(col2hex(input[[paste0("GroupColour", input$correlationColourBy, "__", k)]]), "FF")
  }), input$correlationKeepBucket)
  
  # get the correct coefficient 
  if (input$correlationMethod == "Spearman") {
    cor <- "rho"
  } else if (input$correlationMethod == "Kendall") {
    cor <- "tau"
  } else {
    cor <- "R"
  }
  
  # make the plot
  corr_plot1 <- ggplot(countsCorrelation, aes(x = .data[[input$correlationGene1]], y = .data[[input$correlationGene2]]))
  corr_plot1 <- corr_plot1 + geom_smooth(method = "lm", se = TRUE, colour = "black")
  corr_plot1 <- corr_plot1 + stat_cor(method = tolower(input$correlationMethod), cor.coef.name = cor, size = input$textSizeCorrelation/3, show.legend = FALSE, label.x.npc = tolower(input$correlationCorPosX), label.y.npc = tolower(input$correlationCorPosY), label.args = list(fontface = "bold"))
  if (input$correlationShapeBy == "None") {
    corr_plot1 <- corr_plot1 + geom_point(aes(fill = .data[[input$correlationColourBy]]), pch = 21, size = input$correlationPointSize, stroke = input$correlationPointStroke, alpha = input$correlationPointAlpha)
  } else {
    corr_plot1 <- corr_plot1 + geom_point(aes(fill = .data[[input$correlationColourBy]], shape = .data[[input$correlationShapeBy]]), size = input$correlationPointSize, stroke = input$correlationPointStroke, alpha = input$correlationPointAlpha)
    corr_plot1 <- corr_plot1 + scale_shape_manual(values = c(rep(c(21,22,23,24,25), times = 10))[1:nlevels(as.factor(countsCorrelation[[input$correlationShapeBy]]))]) 
    corr_plot1 <- corr_plot1 + guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)), shape = guide_legend(override.aes = list(fill = "black", size = 3)))
  }
  if (input$correlationShowNames) {
    corr_plot1 <- corr_plot1 + geom_label_repel(
      aes(label = Name, fill = .data[[input$correlationColourBy]]), 
      show.legend = FALSE, 
      alpha = input$correlationPointAlpha, 
      size = input$textSizeCorrelation/3, 
      nudge_x = input$correlationLabelPosX/100, 
      nudge_y = input$correlationLabelPosY/100, 
      segment.size = 0.5, 
      segment.color = "black", 
      max.overlaps = input$correlationMaxOverlaps
      )
  }
  corr_plot1 <- corr_plot1 + scale_fill_manual(breaks = names(coloursCorrelation), values = as.character(coloursCorrelation))
  corr_plot1 <- corr_plot1 + theme_prism(base_size = input$textSizeCorrelation)
  corr_plot1 <- corr_plot1 + theme(legend.text = element_text(size = input$textSizeCorrelation, color = "black", face = "bold"))
  corr_plot1
  correlationPlot$plot <- corr_plot1
})

# view the plot 
output$correlationStatic <- renderPlot({
  req(!is.null(correlationPlot$plot))
  print(correlationPlot$plot)
}, 
  width = function(){as.numeric(input$figWidthCorrelation)},
  height = function(){as.numeric(input$figHeightCorrelation)}
)
# download the plot 
output$dlCorrelationButton <- downloadHandler(
  filename = function() {paste(input$correlationGene1, input$correlationGene2, "correlation.pdf", sep = "_")},
  content = function(file) {
    req(!is.null(correlationPlot$plot))
    pdf(
      file = file, 
      width = (input$figWidthCorrelation/85), 
      height = (input$figHeightCorrelation/90),
      onefile = FALSE
    )
    print(correlationPlot$plot)
    dev.off()
  }
)