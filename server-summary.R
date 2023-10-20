req(!is.null(inputDataReactive()$se))

se <- inputDataReactive()$se
factorNames <- inputDataReactive()$factorNames
countList <- inputDataReactive()$countList
factorLevels <- inputDataReactive()$factorLevels

lapply(seq_along(factorLevels), function(i) {
  updateSelectInput(
    session = session, 
    inputId = paste0("GroupColour", i),
    label = factorLevels[[i]]
  )
})

observeEvent(inputDataReactive()$dataType, {
  if (inputDataReactive()$dataType == "RNASeq") {
    output$linkToMultiDEG <- renderUI({
      box(
        title = "Compare Multiple DE Tests!",
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        tags$p("Do you want to compare the results in this DE test with other DE tests you have run? Then, good news, because there is an app for just that! Head to https://fgcz-shiny.uzh.ch/multiDEG to find out more!"),
        h4("Try the MultiDEG app:", a("MultiDEG", href="https://fgcz-shiny.uzh.ch/multiDEG/", target = "_blank"))
      )
    })
    output$referencesText <- renderUI({
      box(
        title = "Cite us!",
        width = 12, 
        solidHeader = TRUE,
        status = "primary",
        tags$p("If you plan to publish results generated in SUSHI, and/or figures from this app, please consider citing us!"),
        tags$b("SUSHI: Hatekeyama et al., 2016"),
        tags$p("Hatakeyama, M., Opitz, L., Russo, G. Qi, W., Schlapbach, R., Rehrauer, H. (2016) SUSHI: an exquisite recipe for fully documented, reproducible and reusable NGS data analysis. BMC Bioinformatics 17, 228. https://doi.org/10.1186/s12859-016-1104-8"),
        tags$b("exploreDEG Shiny app: Leary and Rehrauer, 2023"),
        tags$p("Leary, P. and Rehrauer, H. (2023) exploreDEG Interactive Shiny App. Zenodo, 20 July, available at: https://doi.org/10.5281/zenodo.8167437"),
        tags$p("For the full list of citations required for the generation of these results, and/or for a written methods section, please email us at sequencing@fgcz.ethz.ch.")
      )
    })
  } else if (inputDataReactive()$dataType == "proteomics") {
    output$referencesText <- renderUI({
      box(
        title = "Cite us!",
        width = 12, 
        solidHeader = TRUE,
        status = "primary",
        tags$p("If you plan to publish results generated in b-fabric and prolfqua, and/or figures from this app, please consider citing us!"),
        tags$b("B-Fabric: Panse et al., 2022"),
        tags$p('Panse, Christian, Trachsel, Christian and Türker, Can. "Bridging data management platforms and visualization tools to enable ad-hoc and smart analytics in life sciences" Journal of Integrative Bioinformatics, vol. 19, no. 4, 2022, pp. 20220031. https://doi.org/10.1515/jib-2022-0031'),
        tags$b("Prolfqua: Wolski et al., 2023"),
        tags$p('Witold E. Wolski, Paolo Nanni, Jonas Grossmann, Maria d’Errico, Ralph Schlapbach, and Christian Panse. Journal of Proteome Research 2023 22 (4), 1092-1104. DOI: 10.1021/acs.jproteome.2c00441'),
        tags$b("exploreDEG Shiny app: Leary and Rehrauer, 2023"),
        tags$p("Leary, P. and Rehrauer, H. (2023) exploreDEG Interactive Shiny App. Zenodo, 20 July, available at: https://doi.org/10.5281/zenodo.8167437"),
        tags$p("For the full list of citations required for the generation of these results, and/or for a written methods section, please email us at sequencing@fgcz.ethz.ch.")
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
  
  # Proteomics: Create metadata table
  if (inputDataReactive()$dataType == "proteomics") {
    output$settingsTable <- function() {
      kable(
        as.data.frame(metadata(inputDataReactive()$se)),
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