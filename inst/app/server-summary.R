eds <- inputDataReactive()$.eds
req(eds)
se <- eds$backend() # SE-only escape hatch for legacy helpers below

# The contrast selector (canonical `contrastSelected` + per-tab mirrors) is
# rendered and synced in server-initInputData.R; modules read selectedContrast().

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
cc2 <- catalystCols
cc2 <- colorspace::darken(cc2, 0.4)
catalystCols <- c(catalystCols, cc2)
catalystCols <- paste0(catalystCols, "FF")

n <- 60
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
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
    multiple = TRUE,
    width = "85%",
    options = list(
      placeholder = "Select features",
      plugins = list("remove_button", "drag_drop", "restore_on_backspace", "clear_button")
    )
  )
})
outputOptions(output, "colourPaletteUI", suspendWhenHidden = FALSE)

# Colour picker + optional display label for each level, one column per factor:
observeEvent(input$colourPalette, {
  output$colourPickerUI <- renderUI({
    factor_names <- eds$factor_names()
    fluidRow(lapply(factor_names, function(fname) {
      levs <- eds$factor_levels(fname)
      column(
        width = 3,
        h4(fname),
        lapply(seq_along(levs), function(j) {
          tagList(
            colourpicker::colourInput(
              inputId = paste0("GroupColour", fname, levs[[j]]),
              label = levs[[j]],
              # save-state restore: seed the saved colour (the picker is rebuilt on
              # colourPalette change, so this is where a restored colour lands).
              value = restore_default(
                paste0("GroupColour", fname, levs[[j]]),
                rep(as.character(unlist(colourPaletteList[input$colourPalette])), times = 5)[j]
              ),
              palette = "square",
              closeOnClick = TRUE,
              returnName = TRUE,
              width = "85%"
            ),
            groupLabelInput(
              inputId = paste0("GroupLabel", fname, levs[[j]]),
              value = restore_default(paste0("GroupLabel", fname, levs[[j]]), ""),
              placeholder = paste0("Plot label (default: ", levs[[j]], ")")
            )
          )
        })
      )
    }))
  })
  outputOptions(output, "colourPickerUI", suspendWhenHidden = FALSE)
})

# Label text boxes send on blur/Enter (not every keystroke), so plots don't
# redraw mid-typing.
groupLabelInput <- function(inputId, value, placeholder) {
  text_input_on_blur(inputId = inputId, label = NULL, value = value, placeholder = placeholder, width = "85%")
}

# Warn about labels that won't survive a PDF download (base pdf() is Latin-1;
# Greek is handled via the Symbol font), or that fell back to plain text
# because the markup couldn't be parsed.
output$groupLabelWarningsUI <- renderUI({
  ids <- grep("^GroupLabel", names(input), value = TRUE)
  vals <- vapply(ids, function(i) paste(input[[i]], collapse = ""), character(1))
  vals <- vals[nzchar(trimws(vals))]
  req(length(vals) > 0)
  parsed <- lapply(vals, parse_level_label)
  plain <- vapply(parsed, `[[`, character(1), "plain")
  msgs <- list()
  bad <- names(vals)[!vapply(parsed, `[[`, logical(1), "ok")]
  if (length(bad)) {
    msgs <- c(
      msgs,
      list(tags$li(paste0(
        "Couldn't parse the formatting in: ",
        paste(vals[bad], collapse = "; "),
        ". Shown as plain text instead - check for unclosed tags."
      )))
    )
  }
  if (label_has_greek(plain)) {
    msgs <- c(
      msgs,
      list(tags$li(
        "Greek letters are drawn with the standard PDF Symbol font in PDF downloads,",
        "so they look slightly different from the surrounding text."
      ))
    )
  }
  unsupported <- label_pdf_unsupported_chars(plain)
  if (length(unsupported)) {
    msgs <- c(
      msgs,
      list(tags$li(paste0(
        "These characters can't be drawn in PDF downloads and will appear as '.': ",
        paste(unsupported, collapse = " "),
        " (PNG and on-screen are fine)."
      )))
    )
  }
  req(length(msgs) > 0)
  div(class = "alert alert-warning", style = "margin-top: 10px;", tags$ul(msgs))
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
  filename = function() {
    generate_filename(
      base_name = "colour_template",
      design_name = designReactive(),
      extension = "xlsx"
    )
  },
  content = function(file) {
    colourTemplate <- saved_inputs()
    colnames(colourTemplate) <- c("Level", "Colour")
    colourTemplate <- colourTemplate[grep("^GroupColour", colourTemplate$Level), ]
    colourTemplate$Label <- vapply(
      sub("^GroupColour", "GroupLabel", colourTemplate$Level),
      function(id) {
        v <- input[[id]]
        if (is.null(v)) "" else v
      },
      character(1),
      USE.NAMES = FALSE
    )
    openxlsx::write.xlsx(colourTemplate, file)
  }
)

observeEvent(input$importColourFile, {
  importedColours <- openxlsx::read.xlsx(input$importColourFile[1, "datapath"], colNames = T)
  importedColours <- importedColours %>% data.frame(check.names = F)

  colourTemplate <- saved_inputs()
  colnames(colourTemplate) <- c("Level", "Colour")
  colourTemplate <- colourTemplate[grep("^GroupColour", colourTemplate$Level), ]

  if (any(!importedColours$Level %in% colourTemplate$Level)) {
    showModal(modalDialog(
      title = "Uh oh...",
      "You uploaded a file that has different factors and/or levels to the dataset. Are you sure it's from this study?",
      easyClose = TRUE
    ))
  } else {
    lapply(importedColours$Level, function(l) {
      updateColourInput(session = session, inputId = l, value = importedColours$Colour[importedColours$Level == l])
    })
    # Label column is optional: older templates only have Level + Colour
    if ("Label" %in% colnames(importedColours)) {
      lapply(importedColours$Level, function(l) {
        lab <- importedColours$Label[importedColours$Level == l]
        updateTextInput(
          session = session,
          inputId = sub("^GroupColour", "GroupLabel", l),
          value = if (is.na(lab)) "" else as.character(lab)
        )
      })
    }
  }
})

observeEvent(eds$data_type(), {
  if (eds$data_type() == "RNASeq" && eds$has_de()) {
    output$referencesText <- renderUI({
      box(
        title = "Cite us!",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tags$p(
          "If you plan to publish results generated in SUSHI, and/or figures from this app, please consider citing us!"
        ),
        tags$b("SUSHI: Hatakeyama et al., 2016"),
        tags$p(
          "Hatakeyama, M., Opitz, L., Russo, G. Qi, W., Schlapbach, R., Rehrauer, H. (2016) SUSHI: an exquisite recipe for fully documented, reproducible and reusable NGS data analysis. BMC Bioinformatics 17, 228. https://doi.org/10.1186/s12859-016-1104-8"
        ),
        tags$b("exploreDE Shiny app: Leary and Rehrauer, 2023"),
        tags$p(
          "Peter Leary, & Hubert Rehrauer. (2023). exploreDE Interactive Shiny App. Zenodo. https://doi.org/10.5281/zenodo.13927692"
        ),
        tags$p(
          "For the full list of citations required for the generation of these results, and/or for a written methods section, please email us at sequencing@fgcz.ethz.ch."
        )
      )
    })
    # Methods download UI - show for RNASeq (validation happens in download handler)
    output$methodsDownloadUI <- renderUI({
      tagList(
        # hr(style = "border-top: 1px solid #000000;"),
        h4("Download Methods Section"),
        tags$div(
          style = "background-color: #fff3cd; border: 2px solid #ffc107; padding: 10px; margin-bottom: 10px; border-radius: 5px;",
          tags$b(style = "color: #856404;", "⚠️ Beta feature warning!"),
          tags$p(
            style = "color: #856404; margin-top: 5px; margin-bottom: 0;",
            "This feature is very much a work-in-progress and limited. But, hopefully at least, it proves to be a useful starting point.
            Having said that, all text in this file must be very carefully inspected and verified prior to submission anywhere. If you have any questions or doubts, please reach out to your bioinformatics coach."
          )
        ),
        helpText(
          "Download a formatted methods section as a Word document (.docx). Only available for DESeq2 and edgeR RNA-seq analyses."
        ),
        downloadButton("dlMethods_summary", "Download methods (Word)")
      )
    })

    output$countsDownloadUI <- renderUI({
      tagList(
        # hr(style = "border-top: 1px solid #000000;"),
        h4("Download Count Files"),
        helpText("Select the type of count data you wish to download for all samples in the dataset."),
        selectizeInput(
          inputId = "downloadCount",
          "Select Count Table",
          choices = eds$assay_names(),
          selected = eds$assay_names()[1],
          width = "33%"
        ),
        textInput(
          inputId = "filename_summary",
          label = "File name",
          value = "Counts",
          placeholder = "Filename",
          width = "33%"
        ),
        downloadButton("dlTable_summary", "Download selected count file"),
        checkboxInput(inputId = "log2DownloadedCounts", label = "Log2 Counts?", value = FALSE),
        helpText("Log2 adds a psuedocount of 1, and only works on counts not already in log2 scale."),
        checkboxInput(
          inputId = "includeDEStats",
          label = "Include DE statistics (log2FC, p-value, FDR)?",
          value = FALSE
        ),
        helpText("Adds DE statistics from the current contrast as extra columns.")
      )
    })
  } else if (eds$omics_class() == "proteomics") {
    # No methods download for proteomics
    output$methodsDownloadUI <- renderUI({
      NULL
    })

    output$referencesText <- renderUI({
      box(
        title = "Cite us!",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tags$p(
          "If you plan to publish results generated in b-fabric and prolfqua/prolfquapp, and/or figures from this app, please consider citing us!"
        ),
        tags$b(
          "prolfqua: A Comprehensive R‑Package for Proteomics Differential Expression Analysis, Wolski et al., 2023"
        ),
        tags$p(
          "Witold E. Wolski, Paolo Nanni, Jonas Grossmann, Maria d’Errico, Ralph Schlapbach, and Christian Panse. Journal of Proteome Research 2023 22 (4), 1092-1104. DOI: 10.1021/acs.jproteome.2c00441"
        ),
        tags$b("prolfquapp: Shiny Application for Proteomics Differential Expression Analysis, Wolski et al."),
        tags$p(
          "Witold E. Wolski, et al. prolfquapp: Differential Expression Analysis for Proteomics Data. https://github.com/fgcz/prolfquapp"
        ),
        tags$b("b-fabric: Panse et al., 2022"),
        tags$p(
          "Panse, Christian, Trachsel, Christian and Türker, Can. 'Bridging data management platforms and visualization tools to enable ad-hoc and smart analytics in life sciences' Journal of Integrative Bioinformatics, vol. 19, no. 4, 2022, pp. 20220031. https://doi.org/10.1515/jib-2022-0031"
        ),
        tags$b("exploreDE Shiny app: Leary and Rehrauer, 2023"),
        tags$p(
          "Peter Leary, & Hubert Rehrauer. (2023). exploreDE Interactive Shiny App. Zenodo. https://doi.org/10.5281/zenodo.13927692"
        ),
        tags$p(
          "For the full list of citations required for the generation of these results, and/or for a written methods section, please email us at proteomics@fgcz.ethz.ch."
        )
      )
    })
    output$countsDownloadUI <- renderUI({
      tagList(
        hr(style = "border-top: 1px solid #000000;"),
        h4("Download Feature Value Files"),
        helpText(
          "Select the type of value data you wish to download for all samples in the dataset. For proteomics, this would be intensity values."
        ),
        selectizeInput(
          inputId = "downloadCount",
          "Select Value Table",
          choices = eds$assay_names(),
          selected = eds$assay_names()[1],
          width = "33%"
        ),

        textInput(inputId = "filename_summary", label = NULL, value = "Counts", placeholder = "Filename"),
        downloadButton("dlTable_summary", "Download selected file"),
        checkboxInput(
          inputId = "includeDEStats",
          label = "Include DE statistics (log2FC, p-value, FDR)?",
          value = FALSE
        ),
        helpText("Adds DE statistics from the current contrast as extra columns.")
      )
    })
  } else if (eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
    # No methods download for CountQC
    output$methodsDownloadUI <- renderUI({
      NULL
    })

    output$countsDownloadUI <- renderUI({
      tagList(
        # hr(style = "border-top: 1px solid #000000;"),
        h4("Download Count Files"),
        helpText("Select the type of count data you wish to download for all samples in the dataset."),
        selectizeInput(
          inputId = "downloadCount",
          "Select Count Table",
          choices = eds$assay_names(),
          selected = eds$assay_names()[1],
          width = "33%"
        ),

        textInput(inputId = "filename_summary", label = NULL, value = "Counts", placeholder = "Filename"),
        downloadButton("dlTable_summary", "Download selected count file"),
        checkboxInput(inputId = "log2DownloadedCounts", label = "Log2 Counts?", value = FALSE),
        helpText("Log2 adds a psuedocount of 1, and only works on counts not already in log2 scale.")
      )
    })
  }
})

# Workunit link — extract B-Fabric workunit ID from the dataset URL
output$workunitsLinkUI <- renderUI({
  wu_id <- NA_character_
  if (!is.null(dataUrl) && nzchar(dataUrl)) {
    # The workunit ID is only encoded in the URL for proteomics (…/workunit_<id>/…).
    # Genomics URLs carry only the order number (o<id>), which is NOT a workunit —
    # so we deliberately show no link there rather than link to the wrong entity.
    if (grepl("workunit_[0-9]+", dataUrl)) {
      wu_id <- sub(".*workunit_([0-9]+).*", "\\1", dataUrl)
    }
  }
  if (!is.na(wu_id)) {
    bf_url <- paste0("https://fgcz-bfabric.uzh.ch/bfabric/userlab/show-workunit.html?id=", wu_id)
    tags$p(
      tags$b("Input workunit: "),
      tags$a(href = bf_url, target = "_blank", paste0("Workunit ", wu_id, " on B-Fabric"))
    )
  }
})

# updateSelectInput(session = session, inputId = "downloadCount", choices = names(countList), selected = names(countList)[1])

observe({
  # Create a table of the samples available in this app instance
  output$inputTable <- function() {
    cD <- eds$samples()
    cD <- cD %>% rownames_to_column("Sample ID")
    cD <- cD %>% dplyr::select(all_of(c("Sample ID", eds$factor_names())))
    kable(
      cD,
      row.names = FALSE,
      format = "html"
    ) %>%
      kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE, position = "left")
  }

  # Create a link to download the counts as an Excel file
  # Create a link to download the counts as an Excel file
  output$dlTable_summary <- downloadHandler(
    filename = function() {
      generate_filename(
        base_name = input$filename_summary,
        design_name = designReactive(),
        extension = "xlsx"
      )
    },
    content = function(file) {
      if (eds$data_type() %in% c("RNASeq", "CountQC", "CrisprCountQC")) {
        ctdl <- eds$assay(input$downloadCount)
        if (input$log2DownloadedCounts & !input$downloadCount %in% c("normalised_log2", "vst")) {
          ctdl <- log2(ctdl + 1)
        }
        ctdl <- as.data.frame(ctdl) %>% rownames_to_column("feature_id")

        # Try to add annotation columns if they exist (RNA-seq has DE; CountQC does not)
        if (eds$data_type() == "RNASeq" && eds$has_de()) {
          seqAnno <- eds$de(selectedContrast())
          anno_cols <- c("gene_id", "transcript_id", "feature_id")
          valid_cols <- intersect(anno_cols, colnames(seqAnno))
          if (length(valid_cols) > 0) {
            ctdl <- ctdl %>% left_join(seqAnno[, valid_cols, drop = FALSE], by = "feature_id")
          }

          # Add DE statistics if requested (RNASeq only, not CountQC)
          if (!is.null(input$includeDEStats) && input$includeDEStats) {
            de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
            valid_de_cols <- intersect(de_cols, colnames(seqAnno))
            if (length(valid_de_cols) > 1) {
              de_data <- seqAnno[, valid_de_cols, drop = FALSE]
              ctdl <- ctdl %>% left_join(de_data, by = "feature_id")
            }
          }
        }

        export_excel_standard(
          df = ctdl,
          filename = file,
          sheet_name = "Count Data",
          table_name = "CountData"
        )
      } else if (eds$omics_class() == "proteomics") {
        ctdl <- as.data.frame(eds$assay(input$downloadCount)) %>% rownames_to_column("feature_id")

        # Add DE statistics if requested
        if (
          !is.null(input$includeDEStats) &&
            input$includeDEStats &&
            isTruthy(selectedContrast())
        ) {
          de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
          seqAnno_c <- eds$de(selectedContrast())
          valid_de_cols <- intersect(de_cols, colnames(seqAnno_c))
          if (length(valid_de_cols) > 1) {
            de_data <- seqAnno_c[, valid_de_cols, drop = FALSE]
            ctdl <- ctdl %>% left_join(de_data, by = "feature_id")
          }
        }

        export_excel_standard(
          df = ctdl,
          filename = file,
          sheet_name = "Feature Data",
          table_name = "FeatureData"
        )
      }
    }
  )

  # Download methods section as Word document
  output$dlMethods_summary <- downloadHandler(
    filename = function() {
      generate_filename(
        base_name = "methods",
        design_name = designReactive(),
        extension = "docx"
      )
    },
    content = function(file) {
      tryCatch(
        {
          # Only proceed if RNASeq data
          req(eds$data_type() == "RNASeq")

          # Methods generation is SUSHI/ezRun-only.
          if (eds$backend_kind() != "se") {
            showNotification(
              "Methods text is only available for SUSHI-processed datasets.",
              type = "message",
              duration = 8
            )
            return(NULL)
          }

          # Reach for the SE + param via the backend escape hatch.
          se <- eds$backend()
          param <- S4Vectors::metadata(se)$param

          # Extract dataDir from attribute (still emitted by the legacy shim)
          dataDir <- attr(inputDataReactive(), "dataDir")
          if (is.null(dataDir)) {
            dataDir <- ""
          }

          # Detect aligner and DE method
          aligner <- detectAligner(dataDir)
          de_method <- detectDEMethod(dataDir, param)

          # Additional safety check - only proceed for DESeq2 or edgeR (case insensitive)
          if (!tolower(de_method) %in% c("deseq2", "edger")) {
            showNotification(
              paste0("Methods download is only available for DESeq2 and edgeR analyses. Detected method: ", de_method),
              type = "error",
              duration = 10
            )
            return(NULL)
          }

          # Normalize DE method capitalization for consistency
          de_method <- if (tolower(de_method) == "deseq2") {
            "DESeq2"
          } else if (tolower(de_method) == "edger") {
            "edgeR"
          } else {
            de_method
          }

          fastp_version <- getFastpVersion(dataDir)

          # Get package versions
          pkg_names <- c("DESeq2", "edgeR", "clusterProfiler")
          pkg_versions <- getPackageVersions(pkg_names)

          # Check if GO was run
          run_go <- if (!is.null(param$runGO)) param$runGO else FALSE

          # Get reference build for GENCODE detection
          ref_build <- if (!is.null(param$refBuild)) param$refBuild else NULL

          # Build citations
          citations <- buildCitationsList(aligner, de_method, run_go, ref_build)

          # Prepare raw reads text
          raw_reads_text <- ""
          dataset <- eds$samples()
          if (!is.null(dataset)) {
            # Look for common read count column names (Reads, Total Reads, Read Count, etc)
            # Look for common read count column names (Reads, Total Reads, Read Count, etc)
            cols <- colnames(dataset)
            # Priority 1: Exact matches (with or without underscores/spaces)
            read_cols <- grep(
              "^Read_?Count$|^Reads$|^Raw_?Reads$|^Total_?Reads$",
              cols,
              value = TRUE,
              ignore.case = TRUE
            )

            # Priority 2: Contains "Read" and "Count" (e.g., "PF Read Count") if no exact match
            if (length(read_cols) == 0) {
              read_cols <- grep("Read", cols, value = TRUE, ignore.case = TRUE)
              read_cols <- grep("Count", read_cols, value = TRUE, ignore.case = TRUE)
              # Exclude filenames
              read_cols <- grep("File|Path|Name", read_cols, value = TRUE, invert = TRUE, ignore.case = TRUE)
            }

            # Priority 3: Just "Reads" but cautious
            if (length(read_cols) == 0) {
              read_cols <- grep("Read", cols, value = TRUE, ignore.case = TRUE)
              read_cols <- grep(
                "File|Path|Name|Mapped|Gene|Feature|Aligned",
                read_cols,
                value = TRUE,
                invert = TRUE,
                ignore.case = TRUE
              )
            }
            if (length(read_cols) > 0) {
              tryCatch(
                {
                  raw_val <- dataset[[read_cols[1]]]
                  raw_reads <- suppressWarnings(as.numeric(as.character(raw_val)))
                  if (!all(is.na(raw_reads))) {
                    min_raw <- round(min(raw_reads, na.rm = TRUE) / 1e6, 1)
                    max_raw <- round(max(raw_reads, na.rm = TRUE) / 1e6, 1)
                    # User requested format: "Individual library size ranged from X to Y million reads."
                    raw_reads_text <- paste0(
                      "Individual library size ranged from ",
                      min_raw,
                      " million to ",
                      max_raw,
                      " million reads. "
                    )
                  }
                },
                error = function(e) {}
              )
            }
          }

          # Generate methods text
          methods_sections <- generateMethodsText(
            param,
            aligner,
            de_method,
            fastp_version,
            pkg_versions,
            se,
            raw_reads_text = raw_reads_text
          )

          # Create Word document using officer package
          doc <- officer::read_docx()

          # Add title
          doc <- doc %>%
            officer::body_add_fpar(
              officer::fpar(officer::ftext("Methods", prop = officer::fp_text(bold = TRUE, font.size = 16)))
            ) %>%
            officer::body_add_par("")

          # Add disclaimer
          doc <- doc %>%
            officer::body_add_fpar(
              officer::fpar(officer::ftext("IMPORTANT NOTICE", prop = officer::fp_text(bold = TRUE, font.size = 14)))
            ) %>%
            officer::body_add_par(paste0(
              "This methods section was automatically generated by the exploreDE Shiny application (beta version). ",
              "All information contained herein MUST be carefully reviewed and verified for accuracy before inclusion in any publication or report. ",
              "Please verify: (1) all software versions are correct, (2) all analysis parameters match your study design, ",
              "(3) all citations are complete and properly formatted, and (4) the description accurately reflects your experimental procedures. ",
              "If you have any questions or concerns, please consult with your bioinformatics coach or contact the FGCZ team at sequencing@fgcz.ethz.ch. ",
              "Generated on: ",
              format(Sys.Date(), "%d %B %Y"),
              "."
            )) %>%
            officer::body_add_par("")

          # Add Library Preparation section (placeholder for user to fill in)
          doc <- doc %>%
            officer::body_add_fpar(
              officer::fpar(officer::ftext("Library Preparation", prop = officer::fp_text(bold = TRUE, font.size = 14)))
            ) %>%
            officer::body_add_par(
              "[User to complete: Include RNA quality assessment methods (Qubit, Fragment Analyzer, RIN values), library preparation kit details, and final library normalisation procedures.]"
            ) %>%
            officer::body_add_par("")

          # Add RNA-Seq Analysis section
          doc <- doc %>%
            officer::body_add_fpar(
              officer::fpar(officer::ftext(
                "RNA-Seq Data Analysis",
                prop = officer::fp_text(bold = TRUE, font.size = 14)
              ))
            ) %>%
            officer::body_add_par("")

          # Add Data Processing section
          doc <- doc %>%
            officer::body_add_fpar(
              officer::fpar(officer::ftext("Data Processing", prop = officer::fp_text(bold = TRUE, font.size = 12)))
            ) %>%
            officer::body_add_par(methods_sections$data_processing) %>%
            officer::body_add_par("")

          # Add Read Trimming section
          doc <- doc %>%
            officer::body_add_fpar(
              officer::fpar(officer::ftext(
                "Read Trimming and Quality Control",
                prop = officer::fp_text(bold = TRUE, font.size = 12)
              ))
            ) %>%
            officer::body_add_par(methods_sections$trimming) %>%
            officer::body_add_par("")

          # Add Mapping section
          doc <- doc %>%
            officer::body_add_fpar(
              officer::fpar(officer::ftext(
                "Read Mapping and Quantification",
                prop = officer::fp_text(bold = TRUE, font.size = 12)
              ))
            ) %>%
            officer::body_add_par(methods_sections$mapping) %>%
            officer::body_add_par("")

          # Add DE Analysis section
          doc <- doc %>%
            officer::body_add_fpar(
              officer::fpar(officer::ftext(
                "Differential Expression Analysis",
                prop = officer::fp_text(bold = TRUE, font.size = 12)
              ))
            ) %>%
            officer::body_add_par(methods_sections$de_analysis) %>%
            officer::body_add_par("")

          # Add GO Analysis section if applicable
          if (!is.null(methods_sections$go_analysis)) {
            doc <- doc %>%
              officer::body_add_fpar(
                officer::fpar(officer::ftext(
                  "Gene Ontology Enrichment Analysis",
                  prop = officer::fp_text(bold = TRUE, font.size = 12)
                ))
              ) %>%
              officer::body_add_par(methods_sections$go_analysis) %>%
              officer::body_add_par("")
          }

          # Add References section
          doc <- doc %>%
            officer::body_add_fpar(
              officer::fpar(officer::ftext("References", prop = officer::fp_text(bold = TRUE, font.size = 14)))
            ) %>%
            officer::body_add_par("")

          for (citation_name in names(citations)) {
            doc <- doc %>%
              officer::body_add_par(citations[[citation_name]]) %>%
              officer::body_add_par("")
          }

          # Save the document
          print(doc, target = file)
        },
        error = function(e) {
          showNotification(
            paste0("Error generating methods document: ", e$message),
            type = "error",
            duration = 10
          )
          cat("Methods download error:", e$message, "\n")
        }
      )
    }
  )

  # Proteomics: Create metadata table
  if (eds$omics_class() == "proteomics") {
    output$settingsTable <- function() {
      if (eds$backend_kind() != "se") {
        settings <- character()
        settings["Data type:"] <- eds$data_type()
        settings["Number of features:"] <- eds$n_features()
        settings["Number of samples:"] <- eds$n_samples()
        return(
          kable(as.data.frame(settings), row.names = TRUE, col.names = "Setting", format = "html") %>%
            kable_styling(bootstrap_options = "striped", full_width = FALSE, position = "left")
        )
      }
      md <- S4Vectors::metadata(eds$backend())
      if ("contrasts" %in% names(md)) {
        contrastSummary <- as.data.frame(md$contrasts)
      } else {
        contrastSummary <- as.data.frame(md)
      }
      kable(
        contrastSummary,
        row.names = FALSE,
        format = "html"
      ) %>%
        kable_styling(
          bootstrap_options = "striped",
          full_width = FALSE,
          position = "left"
        )
    }
  }

  # RNA Seq: Create settings and thresholds table
  if (eds$data_type() == "RNASeq") {
    output$settingsTable <- function() {
      # Uploaded (non-SUSHI) dataset has no ezRun param object — show an
      # accessor-derived summary instead of the SUSHI settings block.
      if (eds$backend_kind() != "se") {
        settings <- character()
        settings["Data type:"] <- eds$data_type()
        settings["Organism:"] <- as.character(eds$organism())
        settings["Reference build:"] <- as.character(eds$ref_build())
        settings["Number of features:"] <- eds$n_features()
        settings["Number of samples:"] <- eds$n_samples()
        return(
          kable(as.data.frame(settings), row.names = TRUE, col.names = "Setting", format = "html") %>%
            kable_styling(bootstrap_options = "striped", full_width = FALSE, position = "left")
        )
      }
      kable(
        count_result_summary(S4Vectors::metadata(se)$param, se),
        row.names = TRUE,
        col.names = "Setting",
        format = "html"
      ) %>%
        kable_styling(
          bootstrap_options = "striped",
          full_width = FALSE,
          position = "left"
        )
    }
    output$summaryTable <- function() {
      settings <- character()
      settings["Number of features:"] <- eds$n_features()
      seqAnno <- if (eds$has_de()) eds$de(selectedContrast()) else NULL
      if (!is.null(seqAnno$isPresentProbe)) {
        settings["Number of features with counts above threshold:"] <-
          sum(seqAnno$isPresentProbe)
      }
      knitr::kable(as.data.frame(settings), format = "html", col.names = "Number", row.names = TRUE) %>%
        kable_styling(
          bootstrap_options = "striped",
          full_width = FALSE,
          position = "left"
        )
    }
  } else if (eds$data_type() %in% c("CountQC", "CrisprCountQC")) {
    output$settingsTable <- function() {
      if (eds$backend_kind() != "se") {
        settings <- character()
        settings["Data type:"] <- eds$data_type()
        settings["Reference:"] <- as.character(eds$ref_build())
        settings["Number of features:"] <- eds$n_features()
        return(
          kable(as.data.frame(settings), row.names = TRUE, col.names = "Setting", format = "html") %>%
            kable_styling(bootstrap_options = "striped", full_width = FALSE, position = "left")
        )
      }
      param <- S4Vectors::metadata(se)$param
      settings <- character()

      if (!is.null(param$refBuild)) {
        settings["Reference:"] <- param$refBuild
      }
      if (!is.null(param$normMethod)) {
        settings["Normalization method:"] <- param$normMethod
      }
      if (!is.null(param$signalThreshold)) {
        settings["Log2 signal threshold:"] <- param$signalThreshold
      }
      if (!is.null(param$signalThreshold)) {
        settings["Linear signal threshold:"] <- round(2^param$signalThreshold, 2)
      }
      if (!is.null(param$featureLevel)) {
        settings["Feature level:"] <- param$featureLevel
      }
      settings["Number of features:"] <- eds$n_features()
      if (!is.null(param$countColumn)) {
        settings["Data Column Used:"] <- param$countColumn
      }

      kable(as.data.frame(settings), row.names = TRUE, col.names = "Setting", format = "html") %>%
        kable_styling(
          bootstrap_options = "striped",
          full_width = FALSE,
          position = "left"
        )
    }
    output$summaryTable <- function() {
      NULL
    }
  }
})
