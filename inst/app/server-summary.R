req(!is.null(inputDataReactive()$se))

se <- inputDataReactive()$se
factorNames <- inputDataReactive()$factorNames
countList <- inputDataReactive()$countList
factorLevels <- inputDataReactive()$factorLevels

catalystCols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", 
  "#FF7B00", "#FDC362", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", 
  "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4", 
  "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", 
  "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
)
catalystCols <- c(catalystCols[c(seq(1,30, by = 2), seq(2, 30, by = 2))])
cc2 <- catalystCols
cc2 <- colorspace::darken(cc2, 0.4)
catalystCols <- c(catalystCols, cc2)
catalystCols <- paste0(catalystCols, "FF")

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colourPaletteList <- list(
  "Catalyst colours" = catalystCols,
  "Paired" = brewer.pal(12, "Paired"),
  "Set1" = brewer.pal(9, "Set1"),
  "Set2" = brewer.pal(8, "Set2"),
  "Set3" = brewer.pal(12, "Set3"),
  "Dark2" = brewer.pal(8, "Dark2"),
  "Dealer's choice" = sample(col_vector, 20)
)

output$colourPaletteUI <- renderUI({
  selectizeInput(
    inputId = "colourPalette", 
    label = "Colour Palette",
    choices = names(colourPaletteList), 
    selected = names(colourPaletteList)[1],
    multiple = TRUE, width = "85%",
    options = list(placeholder = 'Select features', plugins = list('remove_button', 'drag_drop', 'restore_on_backspace', 'clear_button'))
  )
})
outputOptions(output, "colourPaletteUI", suspendWhenHidden = FALSE)

# Colour picker for each of the groups in Condition:
observeEvent(input$colourPalette, {
  output$colourPickerUI <- renderUI({
    lapply(seq_along(factorLevels), function(i) {
      lapply(seq_along(factorLevels[[i]]), function(j) {
        colourpicker::colourInput(
          inputId = paste0("GroupColour", names(factorLevels)[[i]], factorLevels[[i]][[j]]),
          label = paste0(names(factorLevels)[i], ": ", factorLevels[[i]][j]),
          value =  rep(as.character(unlist(colourPaletteList[input$colourPalette])), times = 5)[j],
          palette = "square",
          closeOnClick = TRUE,
          returnName = TRUE, 
          width = "85%"
        )
      })
    })
  })
  outputOptions(output, "colourPickerUI", suspendWhenHidden = FALSE)
})

saved_inputs <- reactiveVal(data.frame(Input = character(), Value = character(), stringsAsFactors = FALSE))
# Observe all inputs and store their values
observe({
  current_inputs <- data.frame(
    Input = names(input),
    Value = sapply(names(input), function(x) paste(input[[x]], collapse = " ")),
    stringsAsFactors = FALSE
  )
  saved_inputs(current_inputs)
})
output$dlColourTemplate <- downloadHandler(
  filename = function() {paste0(design, "_colour_template.xlsx")},
  content = function(file) {
    colourTemplate <- saved_inputs()
    colnames(colourTemplate) <- c("Level", "Colour")
    colourTemplate <- colourTemplate[grep("^GroupColour", colourTemplate$Level), ]
    cat(colourTemplate$Level)
    writexl::write_xlsx(x = colourTemplate, path = file)
  }
)

observeEvent(input$importColourFile, {
  importedColours <- openxlsx::read.xlsx(input$importColourFile[1, 'datapath'], colNames = T)
  importedColours <- importedColours %>% data.frame(check.names = F)
  
  colourTemplate <- saved_inputs()
  colnames(colourTemplate) <- c("Level", "Colour")
  colourTemplate <- colourTemplate[grep("^GroupColour", colourTemplate$Level), ]
  
  if(any(!importedColours$Level %in% colourTemplate$Level)) {
    shinyalert::shinyalert(title = "Uh oh...", text = "You uploaded a file that has different factors and/or levels to the dataset. Are you sure it's from this study?", closeOnEsc = TRUE, closeOnClickOutside = TRUE, showCancelButton = TRUE)
  } else {
    lapply(importedColours$Level, function(l) {
      updateColourInput(session = session, inputId = l, value = importedColours$Colour[importedColours$Level == l])
    })
  }
})

observeEvent(inputDataReactive()$dataType, {
  if (inputDataReactive()$dataType == "RNASeq") {
    output$linkToMultiDEG <- renderUI({
      box(
        title = "Compare Multiple DE Tests!",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tags$p("Do you want to compare the results in this DE test with other DE tests you have run? Then, good news, because there is an app for just that! Head to https://fgcz-shiny.uzh.ch/multiDEG to find out more!"),
        h4("Try the MultiDEG app:", a("MultiDEG", href="https://fgcz-shiny.uzh.ch/multiDEG/", target = "_blank"))
      )
    })
    output$referencesText <- renderUI({
      box(
        title = "Cite us!",
        width = NULL, 
        solidHeader = TRUE,
        status = "primary",
        tags$p("If you plan to publish results generated in SUSHI, and/or figures from this app, please consider citing us!"),
        tags$b("SUSHI: Hatekeyama et al., 2016"),
        tags$p("Hatakeyama, M., Opitz, L., Russo, G. Qi, W., Schlapbach, R., Rehrauer, H. (2016) SUSHI: an exquisite recipe for fully documented, reproducible and reusable NGS data analysis. BMC Bioinformatics 17, 228. https://doi.org/10.1186/s12859-016-1104-8"),
        tags$b("exploreDE Shiny app: Leary and Rehrauer, 2023"),
        tags$p("Peter Leary, & Hubert Rehrauer. (2023). exploreDE Interactive Shiny App. Zenodo. https://doi.org/10.5281/zenodo.13927692"),
        tags$p("For the full list of citations required for the generation of these results, and/or for a written methods section, please email us at sequencing@fgcz.ethz.ch.")
      )
    })
  } else if (inputDataReactive()$dataType == "proteomics") {
    output$referencesText <- renderUI({
      box(
        title = "Cite us!",
        width = NULL, 
        solidHeader = TRUE,
        status = "primary",
        tags$p("If you plan to publish results generated in b-fabric and prolfqua, and/or figures from this app, please consider citing us!"),
        tags$b("prolfqua: A Comprehensive R‑Package for Proteomics Differential Expression Analysis, Wolski et al., 2023"),
        tags$p('Witold E. Wolski, Paolo Nanni, Jonas Grossmann, Maria d’Errico, Ralph Schlapbach, and Christian Panse. Journal of Proteome Research 2023 22 (4), 1092-1104. DOI: 10.1021/acs.jproteome.2c00441'),
        tags$b("b-fabric: Panse et al., 2022"),
        tags$p('Panse, Christian, Trachsel, Christian and Türker, Can. "Bridging data management platforms and visualization tools to enable ad-hoc and smart analytics in life sciences" Journal of Integrative Bioinformatics, vol. 19, no. 4, 2022, pp. 20220031. https://doi.org/10.1515/jib-2022-0031'),
        tags$b("exploreDE Shiny app: Leary and Rehrauer, 2023"),
        tags$p("Peter Leary, & Hubert Rehrauer. (2023). exploreDEG Interactive Shiny App. Zenodo. https://doi.org/10.5281/zenodo.10026461"),
        tags$p("For the full list of citations required for the generation of these results, and/or for a written methods section, please email us at proteomics@fgcz.ethz.ch.")
      )
    })
  }
})

updateSelectInput(session = session, inputId = "downloadCount", choices = names(countList), selected = names(countList)[1])

observe({
  
  # Create a table of the samples available in this app instance
  output$inputTable <- function() {
    cD <- inputDataReactive()$dataset
    cD <- cD %>% rownames_to_column("Sample ID")
    cD <- cD %>% dplyr::select(all_of(c("Sample ID", inputDataReactive()$factorNames)))
    kable(
      cD, 
      row.names = FALSE, format = "html") %>%
      kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE, position = "left")
  }
  
  # Create a link to download the counts as an Excel file
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste0(input$downloadCount, "_counts", ".xlsx")
    },
    content = function(file) { 
      if (inputDataReactive()$dataType == "RNASeq") {
        writexl::write_xlsx(as.data.frame(countList[[input$downloadCount]]) %>% rownames_to_column("gene_name") %>% left_join(inputDataReactive()$seqAnno[,c("gene_id", "transcript_id", "gene_name")]), path = file) 
      } else if (inputDataReactive()$dataType == "proteomics") {
        writexl::write_xlsx(as.data.frame(countList[[input$downloadCount]]) %>% rownames_to_column("gene_name"), path = file) 
      }
    }
  )
  
  # Download all currently selected inputs 
  output$downloadInputsE <- downloadHandler(
    filename = function() {
      "Input_Options.xlsx"
    },
    content = function(file) { 
      x <- saved_inputs()
      colnames(x) <- c("Input", "Value")
      writexl::write_xlsx(x, path = file)
    }
  )
  output$downloadInputsR <- downloadHandler(
    filename = function() {
      "Input_Options.qs"
    },
    content = function(file) { 
      qs::qsave(x = saved_inputs(), file = file, nthreads = 8)
    }
  )
  
  # Download current reactive data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("exploreDE_Data_", format(Sys.time(), "%Y-%m-%d_%H.%M.%S"), ".qs")
    },
    content = function(file) { 
      qs::qsave(x = inputDataReactive(), file = file, nthreads = 8)
    }
  )
  
  # Proteomics: Create metadata table
  if (inputDataReactive()$dataType == "proteomics") {
    output$settingsTable <- function() {
      if ("contrasts" %in% names(metadata(inputDataReactive()$se))) {
        contrastSummary <- as.data.frame(metadata(inputDataReactive()$se)$contrasts)
      } else {
        contrastSummary <- as.data.frame(metadata(inputDataReactive()$se))
      }
      kable(
        contrastSummary,
        row.names = FALSE,
        format = "html"
        ) %>%
        kable_styling(
          bootstrap_options = "striped", full_width = FALSE,
          position = "left"
        )
    }
  }
  
  # RNA Seq: Create settings and thresholds table
  if (inputDataReactive()$dataType == "RNASeq") {
    output$settingsTable <- function() {
      kable(makeCountResultSummary(inputDataReactive()$param, se),
            row.names = TRUE,
            col.names = "Setting", format = "html"
      ) %>%
        kable_styling(
          bootstrap_options = "striped", full_width = FALSE,
          position = "left"
        )
    }
    output$summaryTable <- function() {
      settings <- character()
      settings["Number of features:"] <- nrow(se)
      if (!is.null(rowData(se)$isPresentProbe)) {
        settings["Number of features with counts above threshold:"] <-
          sum(rowData(se)$isPresentProbe)
      }
      knitr::kable(as.data.frame(settings),
                   format = "html",
                   col.names = "Number", row.names = TRUE
      ) %>%
        kable_styling(
          bootstrap_options = "striped", full_width = FALSE,
          position = "left"
        )
    }
  }
})