# Create reactive for design to use across all outputs and download handlers
designReactive <- reactive({
  req(inputDataReactive())
  if (inputDataReactive()$dataType == "RNASeq") {
    return(inputDataReactive()$design)
  } else if (inputDataReactive()$dataType == "proteomics") {
    if (isTruthy(input$contrastSelected)) {
      return(input$contrastSelected)
    } else {
      return("No contrast selected")
    }
  } else if (inputDataReactive()$dataType %in% c("CountQC", "CrisprCountQC")) {
    return(inputDataReactive()$design)
  }
  return("Unknown")
})

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
    multiple = TRUE, width = "85%",
    options = list(placeholder = "Select features", plugins = list("remove_button", "drag_drop", "restore_on_backspace", "clear_button"))
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
          value = rep(as.character(unlist(colourPaletteList[input$colourPalette])), times = 5)[j],
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
  filename = function() {
    generate_filename(
      base_name = "colour_template",
      design_name = design,
      extension = "xlsx"
    )
  },
  content = function(file) {
    colourTemplate <- saved_inputs()
    colnames(colourTemplate) <- c("Level", "Colour")
    colourTemplate <- colourTemplate[grep("^GroupColour", colourTemplate$Level), ]
    writexl::write_xlsx(x = colourTemplate, path = file)
  }
)

observeEvent(input$importColourFile, {
  importedColours <- openxlsx::read.xlsx(input$importColourFile[1, "datapath"], colNames = T)
  importedColours <- importedColours %>% data.frame(check.names = F)

  colourTemplate <- saved_inputs()
  colnames(colourTemplate) <- c("Level", "Colour")
  colourTemplate <- colourTemplate[grep("^GroupColour", colourTemplate$Level), ]

  if (any(!importedColours$Level %in% colourTemplate$Level)) {
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
        h4("Try the MultiDEG app:", a("MultiDEG", href = "https://fgcz-shiny.uzh.ch/multiDEG/", target = "_blank"))
      )
    })
    output$referencesText <- renderUI({
      box(
        title = "Cite us!",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tags$p("If you plan to publish results generated in SUSHI, and/or figures from this app, please consider citing us!"),
        tags$b("SUSHI: Hatakeyama et al., 2016"),
        tags$p("Hatakeyama, M., Opitz, L., Russo, G. Qi, W., Schlapbach, R., Rehrauer, H. (2016) SUSHI: an exquisite recipe for fully documented, reproducible and reusable NGS data analysis. BMC Bioinformatics 17, 228. https://doi.org/10.1186/s12859-016-1104-8"),
        tags$b("exploreDE Shiny app: Leary and Rehrauer, 2023"),
        tags$p("Peter Leary, & Hubert Rehrauer. (2023). exploreDE Interactive Shiny App. Zenodo. https://doi.org/10.5281/zenodo.13927692"),
        tags$p("For the full list of citations required for the generation of these results, and/or for a written methods section, please email us at sequencing@fgcz.ethz.ch.")
      )
    })
    # Methods download UI - show for RNASeq (validation happens in download handler)
    output$methodsDownloadUI <- renderUI({
      tagList(
        hr(style = "border-top: 1px solid #000000;"),
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
        helpText("Download a formatted methods section as a Word document (.docx). Only available for DESeq2 and edgeR RNA-seq analyses."),
        downloadButton("dlMethods_summary", "Download methods (Word)")
      )
    })

    output$countsDownloadUI <- renderUI({
      tagList(
        hr(style = "border-top: 1px solid #000000;"),
        h4("Download Count Files"),
        helpText("Select the type of count data you wish to download for all samples in the dataset."),
        selectizeInput(inputId = "downloadCount", "Select Count Table", choices = names(inputDataReactive()$countList), selected = names(inputDataReactive()$countList)[1], width = "33%"),

        textInput(inputId = "filename_summary", label = NULL, value = "Counts", placeholder = "Filename"),
        downloadButton("dlTable_summary", "Download selected count file"),
        checkboxInput(inputId = "log2DownloadedCounts", label = "Log2 Counts?", value = FALSE),
        helpText("Log2 adds a psuedocount of 1, and only works on counts not already in log2 scale."),
        checkboxInput(inputId = "includeDEStats", label = "Include DE statistics (log2FC, p-value, FDR)?", value = FALSE),
        helpText("Adds DE statistics from the current contrast as extra columns.")
      )
    })
  } else if (inputDataReactive()$dataType == "proteomics") {
    # No methods download for proteomics
    output$methodsDownloadUI <- renderUI({ NULL })

    output$referencesText <- renderUI({
      box(
        title = "Cite us!",
        width = NULL,
        solidHeader = TRUE,
        status = "primary",
        tags$p("If you plan to publish results generated in b-fabric and prolfqua, and/or figures from this app, please consider citing us!"),
        tags$b("prolfqua: A Comprehensive R‑Package for Proteomics Differential Expression Analysis, Wolski et al., 2023"),
        tags$p("Witold E. Wolski, Paolo Nanni, Jonas Grossmann, Maria d’Errico, Ralph Schlapbach, and Christian Panse. Journal of Proteome Research 2023 22 (4), 1092-1104. DOI: 10.1021/acs.jproteome.2c00441"),
        tags$b("b-fabric: Panse et al., 2022"),
        tags$p('Panse, Christian, Trachsel, Christian and Türker, Can. "Bridging data management platforms and visualization tools to enable ad-hoc and smart analytics in life sciences" Journal of Integrative Bioinformatics, vol. 19, no. 4, 2022, pp. 20220031. https://doi.org/10.1515/jib-2022-0031'),
        tags$b("exploreDE Shiny app: Leary and Rehrauer, 2023"),
        tags$p("Peter Leary, & Hubert Rehrauer. (2023). exploreDEG Interactive Shiny App. Zenodo. https://doi.org/10.5281/zenodo.10026461"),
        tags$p("For the full list of citations required for the generation of these results, and/or for a written methods section, please email us at proteomics@fgcz.ethz.ch.")
      )
    })
    output$countsDownloadUI <- renderUI({
      tagList(
        hr(style = "border-top: 1px solid #000000;"),
        h4("Download Feature Value Files"),
        helpText("Select the type of value data you wish to download for all samples in the dataset. For proteomics, this would be intensity values."),
        selectizeInput(inputId = "downloadCount", "Select Value Table", choices = names(inputDataReactive()$countList), selected = names(inputDataReactive()$countList)[1], width = "33%"),

        textInput(inputId = "filename_summary", label = NULL, value = "Counts", placeholder = "Filename"),
        downloadButton("dlTable_summary", "Download selected file"),
        checkboxInput(inputId = "includeDEStats", label = "Include DE statistics (log2FC, p-value, FDR)?", value = FALSE),
        helpText("Adds DE statistics from the current contrast as extra columns.")
      )
    })
  } else if (inputDataReactive()$dataType %in% c("CountQC", "CrisprCountQC")) {
    # No methods download for CountQC
    output$methodsDownloadUI <- renderUI({ NULL })

    output$countsDownloadUI <- renderUI({
      tagList(
        hr(style = "border-top: 1px solid #000000;"),
        h4("Download Count Files"),
        helpText("Select the type of count data you wish to download for all samples in the dataset."),
        selectizeInput(inputId = "downloadCount", "Select Count Table", choices = names(inputDataReactive()$countList), selected = names(inputDataReactive()$countList)[1], width = "33%"),

        textInput(inputId = "filename_summary", label = NULL, value = "Counts", placeholder = "Filename"),
        downloadButton("dlTable_summary", "Download selected count file"),
        checkboxInput(inputId = "log2DownloadedCounts", label = "Log2 Counts?", value = FALSE),
        helpText("Log2 adds a psuedocount of 1, and only works on counts not already in log2 scale.")
      )
    })
  }
})

# Helper functions for methods text generation ----

# Detect aligner from input_dataset.tsv
detectAligner <- function(dataDir) {
  tryCatch({
    dataset_path <- file.path(dataDir, "../input_dataset.tsv")
    if (file.exists(dataset_path)) {
      dataset <- read_tsv(dataset_path, show_col_types = FALSE)

      # Check Count [File] column for patterns
      count_col_name <- grep("Count.*File", names(dataset), value = TRUE, ignore.case = TRUE)
      if (length(count_col_name) > 0) {
        count_col <- dataset[[count_col_name[1]]][1]

        if (grepl("Kallisto", count_col, ignore.case = TRUE)) {
          return("Kallisto")
        } else if (grepl("FeatureCounts", count_col, ignore.case = TRUE)) {
          return("FeatureCounts")
        }
      }

      # Check if BAM column exists (indicates STAR was used)
      if (any(grepl("BAM.*File", names(dataset), ignore.case = TRUE))) {
        return("STAR")
      }
    }
    return("Unknown")
  }, error = function(e) {
    return("Unknown")
  })
}

# Detect DE method
detectDEMethod <- function(dataDir, param) {
  # First check directory path - most reliable
  if (grepl("DESeq2", dataDir, ignore.case = TRUE)) {
    return("DESeq2")
  } else if (grepl("EdgeR|edgeR", dataDir, ignore.case = FALSE)) {
    return("edgeR")
  }

  # Check normalization method
  if (!is.null(param$normMethod)) {
    norm_method <- tolower(param$normMethod)
    if (grepl("tmm|rle", norm_method)) {
      return("edgeR")
    }
  }

  # Check param$testMethod if it exists
  if (!is.null(param$testMethod)) {
    # Extract first word before "using" if present
    test_method <- tolower(param$testMethod)
    test_method <- strsplit(test_method, " using ")[[1]][1]
    test_method <- trimws(test_method)

    # edgeR test methods
    if (test_method %in% c("glm", "lrt", "exacttest", "ql", "quasi-likelihood")) {
      return("edgeR")
    }
    # DESeq2 test methods
    if (test_method %in% c("wald")) {
      return("DESeq2")
    }
  }

  return("DESeq2")  # Default assumption
}

# Format reference build nicely
formatRefBuild <- function(ref_build) {
  # If ref_build is already formatted nicely, return as-is
  if (!grepl("/", ref_build)) {
    return(ref_build)
  }

  # Parse path like "Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30"
  parts <- strsplit(ref_build, "/")[[1]]

  # Extract components
  species <- parts[1]
  database <- parts[2]
  build <- parts[3]

  # Format species
  species_formatted <- if (grepl("Homo_sapiens", species, ignore.case = TRUE)) {
    "human"
  } else if (grepl("Mus_musculus", species, ignore.case = TRUE)) {
    "mouse"
  } else {
    gsub("_", " ", species)
  }

  # Extract release from the last part if it exists
  release <- ""
  if (length(parts) >= 5) {
    release_part <- parts[5]
    # Extract release number from "Release_42-2023-01-30" or "Release_M31-2023-01-30"
    if (grepl("Release_", release_part)) {
      release_num <- sub("Release_([^-]+).*", "\\1", release_part)
      release <- paste0(" (Release ", release_num, ")")
    }
  }

  # Format final string
  formatted <- paste0(database, " ", species_formatted, " genome build ", build, release)
  return(formatted)
}

# Get fastp version from preprocessing logs
getFastpVersion <- function(dataDir) {
  tryCatch({
    parent_dir <- dirname(dataDir)
    log_files <- list.files(parent_dir, pattern = "_preprocessing\\.log$", full.names = TRUE)

    if (length(log_files) > 0) {
      log_content <- readLines(log_files[1], n = 100, warn = FALSE)
      version_line <- grep("fastp v", log_content, value = TRUE)
      if (length(version_line) > 0) {
        version <- sub(".*fastp v([0-9.]+).*", "\\1", version_line[1])
        return(version)
      }
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

# Get package versions
getPackageVersions <- function(packages) {
  versions <- sapply(packages, function(pkg) {
    tryCatch({
      as.character(packageVersion(pkg))
    }, error = function(e) {
      "version unavailable"
    })
  })
  return(versions)
}

# Build citations list
buildCitationsList <- function(aligner, de_method, run_go = FALSE, ref_build = NULL) {
  citations <- list()

  # Always cite SUSHI
  citations$SUSHI <- "Hatakeyama M, Opitz L, Russo G, Qi W, Schlapbach R, Rehrauer H (2016) SUSHI: an exquisite recipe for fully documented, reproducible and reusable NGS data analysis. BMC Bioinformatics 17, 228. https://doi.org/10.1186/s12859-016-1104-8"

  # FastQC
  citations$FastQC <- "Andrews S (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data. Available online at: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
  
  # MultiQC
  citations$MultiQC <- "Ewels P, Magnusson M, Lundin S, Käller M (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354"

  # fastp
  citations$fastp <- "Chen S, Zhou Y, Chen Y, Gu J (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884–i890. https://doi.org/10.1093/bioinformatics/bty560"

  # Aligner-specific citations
  if (aligner == "STAR") {
    citations$STAR <- "Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013) STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21. https://doi.org/10.1093/bioinformatics/bts635"
  } else if (aligner == "Kallisto") {
    citations$Kallisto <- "Bray NL, Pimentel H, Melsted P, Pachter L (2016) Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology 34, 525–527. https://doi.org/10.1038/nbt.3519"
  } else if (aligner == "FeatureCounts") {
    citations$Rsubread <- "Liao Y, Smyth GK, Shi W (2013) The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research 41, e108. https://doi.org/10.1093/nar/gkt214"
    citations$STAR <- "Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013) STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21. https://doi.org/10.1093/bioinformatics/bts635"
  }

  # GENCODE citation if detected in reference build
  if (!is.null(ref_build) && grepl("GENCODE", ref_build, ignore.case = TRUE)) {
    citations$GENCODE <- "Frankish A, Diekhans M, Jungreis I, Lagarde J, Loveland JE, Mudge JM, Sisu C, Wright JC, Armstrong J, Barnes I, et al. (2021) GENCODE: reference annotation for the human and mouse genomes. Nucleic Acids Research 49, D916–D923. https://doi.org/10.1093/nar/gkaa1087"
  }

  # DE method citations
  if (de_method == "DESeq2") {
    citations$DESeq2 <- "Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550. https://doi.org/10.1186/s13059-014-0550-8"
  } else if (de_method == "edgeR") {
    citations$edgeR <- "Robinson MD, McCarthy DJ, Smyth GK (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139–140. https://doi.org/10.1093/bioinformatics/btp616"
  }

  # GO analysis citation if applicable
  if (run_go) {
    citations$clusterProfiler <- "Yu G, Wang LG, Han Y, He QY (2012) clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16, 284–287. https://doi.org/10.1089/omi.2011.0118"
  }

  # R and exploreDE citations
  citations$R <- "R Core Team (2020) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org"
  citations$exploreDE <- "Leary P, Rehrauer H (2024) exploreDE Interactive Shiny App. Zenodo. https://doi.org/10.5281/zenodo.13927692"

  return(citations)
}

# Generate methods text sections
generateMethodsText <- function(param, aligner, de_method, fastp_version, pkg_versions, se, raw_reads_text = "") {
  sections <- list()

  # Final Mapped/Counted Reads (from SE)
  final_reads_info <- ""
  if (!is.null(se)) {
    tryCatch({
      lib_sizes <- colSums(assay(se, "counts"))
      min_size <- round(min(lib_sizes) / 1e6, 1)
      max_size <- round(max(lib_sizes) / 1e6, 1)
      final_reads_info <- paste0("Final library sizes ranged from ", min_size, " million to ", max_size, " million reads. ")
    }, error = function(e) {
      final_reads_info <- ""
    })
  }

  # Read trimming section
  fastp_ver_text <- if (!is.null(fastp_version)) paste0(" (version ", fastp_version, ")") else ""
  sections$trimming <- paste0(
    raw_reads_text,
    "Read quality was inspected using FastQC (Andrews, 2010) and MultiQC (Ewels et al., 2016), and sequencing adaptors removed using fastp",
    fastp_ver_text,
    " (Chen et al., 2018). ",
    "Quality trimming was performed to remove low-quality bases, with reads shorter than 18 bases discarded. ",
    "Poly-X tails were trimmed with a minimum length threshold of 10 bases."
  )

  # Mapping and quantification section
  ref_build_raw <- if (!is.null(param$refBuild)) param$refBuild else "the reference genome"
  ref_build <- formatRefBuild(ref_build_raw)
  feature_level <- if (!is.null(param$featureLevel)) param$featureLevel else "gene"
  transcript_types <- if (!is.null(param$transcriptTypes)) param$transcriptTypes else "protein-coding"

  if (aligner == "Kallisto") {
    sections$mapping <- paste0(
      "Transcript-level abundance was estimated using Kallisto ",
      "(Bray et al., 2016), which performs pseudo-alignment of reads to the ", ref_build, " reference transcriptome. ",
      "Kallisto uses an expectation-maximisation algorithm to quantify transcript abundances. ",
      "Only ", transcript_types, " transcripts were retained for downstream analysis."
    )
  } else if (aligner == "STAR") {
    sections$mapping <- paste0(
      "Alignment of the RNA-seq reads was performed using the STAR aligner ",
      "(Dobin et al., 2013) with the ", ref_build, " reference genome as the reference. ",
      "STAR employs a seed-and-extend approach optimised for RNA-seq data, ",
      "enabling accurate detection of splice junctions. ",
      "Only ", transcript_types, " genes were retained for downstream analysis."
    )
  } else if (aligner == "FeatureCounts") {
    sections$mapping <- paste0(
      "Alignment of the RNA-seq reads was performed using the STAR aligner ",
      "(Dobin et al., 2013) with the ", ref_build, " reference genome as the reference. ",
      "The counting of ", feature_level, "-level expression values was performed using the 'featureCounts' function ",
      "of the R package Rsubread (Liao et al., 2013). ",
      "Only ", transcript_types, " genes were retained for downstream analysis."
    )
  } else {
    sections$mapping <- paste0(
      "Reads were aligned to the ", ref_build, " reference genome and ",
      "quantified at the ", feature_level, " level using standard RNA-seq ",
      "alignment and quantification methods, ",
      "with only ", transcript_types, " genes retained for downstream analysis."
    )
  }
  
  # Append final library size info
  sections$mapping <- paste0(sections$mapping, " ", final_reads_info)

  # Filter, Normalisation and Differential Expression section
  
  deseq2_ver <- if (pkg_versions["DESeq2"] != "version unavailable") paste0(" (version ", pkg_versions["DESeq2"], ")") else ""
  edger_ver <- if (pkg_versions["edgeR"] != "version unavailable") paste0(" (version ", pkg_versions["edgeR"], ")") else ""

  if (de_method == "DESeq2") {
    sections$de_analysis <- paste0(
      "Count data were normalised using DESeq2's median-of-ratios method, which accounts for differences in library size ",
      "and RNA composition bias between samples. ",
      "Differential expression was performed using the DESeq2 Bioconductor R package", deseq2_ver, " ",
      "(Love et al., 2014), which implements a generalised linear model with a negative binomial distribution ",
      "to account for technical and biological variability. ",
      "Dispersion and log2 fold change estimates were moderated using shrinkage methods to improve stability. ",
      "Statistical significance was assessed using the Wald test, with p-values adjusted for multiple testing ",
      "using the Benjamini-Hochberg method to control the false discovery rate (FDR)."
    )
  } else if (de_method == "edgeR") {
    sections$de_analysis <- paste0(
      "Count data were normalised using the trimmed mean of M-values (TMM) method implemented in edgeR, ",
      "which corrects for differences in library size and RNA composition between samples. ",
      "Differential expression was performed using the edgeR Bioconductor R package", edger_ver, " ",
      "(Robinson et al., 2010), which implements a generalised linear model with a negative binomial distribution ",
      "to account for technical and biological variability. ",
      "Dispersion estimates were moderated using empirical Bayes shrinkage to improve reliability for genes with low counts. ",
      "Statistical significance was assessed using the likelihood ratio test (LRT), ",
      "with p-values adjusted for multiple testing using the Benjamini-Hochberg method to control the false discovery rate (FDR)."
    )
  }

  # Gene ontology section (if applicable)
  if (!is.null(param$runGO) && param$runGO) {
    cluster_ver <- if (pkg_versions["clusterProfiler"] != "version unavailable") paste0(" (version ", pkg_versions["clusterProfiler"], ")") else ""
    fdr_ora <- if (!is.null(param$fdrThreshORA)) param$fdrThreshORA else 0.05
    fdr_gsea <- if (!is.null(param$fdrThreshGSEA)) param$fdrThreshGSEA else 0.05

    sections$go_analysis <- paste0(
      "Gene Ontology (GO) term pathway analysis was performed using both the hypergeometric over-representation test ",
      "via the 'enricher' function, and gene-set enrichment analysis via the 'GSEA' function, ",
      "of the clusterProfiler Bioconductor R package", cluster_ver, " ",
      "(Yu et al., 2012), ",
      "with FDR thresholds of ", fdr_ora, " and ", fdr_gsea, ", respectively."
    )
  }

  # Data processing section
  # Get R version info if available
  r_version <- paste0("R version ", R.version$major, ".", R.version$minor)

  # Get Bioconductor version if available
  bioc_version <- ""
  tryCatch({
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      bv <- as.character(BiocManager::version())
      bioc_version <- paste0(" and Bioconductor version ", bv)
    }
  }, error = function(e) {
    bioc_version <- ""
  })

  sections$data_processing <- paste0(
    "All analyses were conducted using the SUSHI framework (Hatakeyama et al., 2016), ",
    "which provides a reproducible and documented workflow for NGS data analysis. ",
    "All R functions were executed on ", r_version, bioc_version, " (R Core Team, 2020). ",
    "Figures were generated using the exploreDE Shiny application using the data generated from the above steps (Leary and Rehrauer, 2024)."
  )

  return(sections)
}

# updateSelectInput(session = session, inputId = "downloadCount", choices = names(countList), selected = names(countList)[1])

observe({
  # Create a table of the samples available in this app instance
  output$inputTable <- function() {
    cD <- inputDataReactive()$dataset
    cD <- cD %>% rownames_to_column("Sample ID")
    cD <- cD %>% dplyr::select(all_of(c("Sample ID", inputDataReactive()$factors)))
    kable(
      cD,
      row.names = FALSE, format = "html"
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
      if (inputDataReactive()$dataType %in% c("RNASeq", "CountQC", "CrisprCountQC")) {
        ctdl <- countList[[input$downloadCount]]
        if (input$log2DownloadedCounts & !input$downloadCount %in% c("Normalised + Log2", "VST")) {
          ctdl <- log2(ctdl + 1)
        }
        ctdl <- as.data.frame(ctdl) %>% rownames_to_column("feature_id")

        # Try to add annotation columns if they exist
        anno_cols <- c("gene_id", "transcript_id", "feature_id")
        valid_cols <- intersect(anno_cols, colnames(inputDataReactive()$seqAnno))
        if (length(valid_cols) > 0) {
          ctdl <- ctdl %>% left_join(inputDataReactive()$seqAnno[, valid_cols, drop = FALSE], by = "feature_id")
        }

        # Add DE statistics if requested (RNASeq only, not CountQC)
        if (inputDataReactive()$dataType == "RNASeq" &&
            !is.null(input$includeDEStats) &&
            input$includeDEStats) {
          de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
          valid_de_cols <- intersect(de_cols, colnames(inputDataReactive()$seqAnno))
          if (length(valid_de_cols) > 1) {  # More than just feature_id
            de_data <- inputDataReactive()$seqAnno[, valid_de_cols, drop = FALSE]
            ctdl <- ctdl %>% left_join(de_data, by = "feature_id")
          }
        }

        export_excel_standard(
          df = ctdl,
          filename = file,
          sheet_name = "Count Data",
          table_name = "CountData"
        )
      } else if (inputDataReactive()$dataType == "proteomics") {
        ctdl <- as.data.frame(countList[[input$downloadCount]]) %>% rownames_to_column("feature_id")

        # Add DE statistics if requested
        if (!is.null(input$includeDEStats) &&
            input$includeDEStats &&
            !is.null(input$contrastSelected)) {
          de_cols <- c("feature_id", "log2Ratio", "pValue", "fdr")
          valid_de_cols <- intersect(de_cols, colnames(inputDataReactive()$seqAnnoList[[input$contrastSelected]]))
          if (length(valid_de_cols) > 1) {
            de_data <- inputDataReactive()$seqAnnoList[[input$contrastSelected]][, valid_de_cols, drop = FALSE]
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

  # Download all currently selected inputs
  output$dlInputsExcel_summary <- downloadHandler(
    filename = function() {
      generate_filename(
        base_name = "Input_Options",
        design_name = designReactive(),
        extension = "xlsx"
      )
    },
    content = function(file) {
      x <- saved_inputs()
      colnames(x) <- c("Input", "Value")
      writexl::write_xlsx(x, path = file)
    }
  )
  output$dlInputsQS_summary <- downloadHandler(
    filename = function() {
      generate_filename(
        base_name = "Input_Options",
        design_name = designReactive(),
        extension = "qs"
      )
    },
    content = function(file) {
      qs::qsave(x = saved_inputs(), file = file, nthreads = 8)
    }
  )

  # Download current reactive data
  output$dlData_summary <- downloadHandler(
    filename = function() {
      generate_filename(
        base_name = "exploreDE_Data",
        design_name = designReactive(),
        extension = "qs"
      )
    },
    content = function(file) {
      qs::qsave(x = inputDataReactive(), file = file, nthreads = 8)
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
      tryCatch({
        # Only proceed if RNASeq data
        req(inputDataReactive()$dataType == "RNASeq")

        param <- inputDataReactive()$param
        se <- inputDataReactive()$se

        # Extract dataDir from attribute
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
        de_method <- if (tolower(de_method) == "deseq2") "DESeq2" else if (tolower(de_method) == "edger") "edgeR" else de_method

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
      dataset <- inputDataReactive()$dataset
      if (!is.null(dataset)) {
        # Look for common read count column names (Reads, Total Reads, Read Count, etc)
        # Look for common read count column names (Reads, Total Reads, Read Count, etc)
        cols <- colnames(dataset)
        # Priority 1: Exact matches (with or without underscores/spaces)
        read_cols <- grep("^Read_?Count$|^Reads$|^Raw_?Reads$|^Total_?Reads$", cols, value = TRUE, ignore.case = TRUE)
        
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
           read_cols <- grep("File|Path|Name|Mapped|Gene|Feature|Aligned", read_cols, value = TRUE, invert = TRUE, ignore.case = TRUE)
        }
        if (length(read_cols) > 0) {
          tryCatch({
            raw_val <- dataset[[read_cols[1]]]
            raw_reads <- suppressWarnings(as.numeric(as.character(raw_val)))
            if (!all(is.na(raw_reads))) {
               min_raw <- round(min(raw_reads, na.rm=TRUE) / 1e6, 1)
               max_raw <- round(max(raw_reads, na.rm=TRUE) / 1e6, 1)
               # User requested format: "Individual library size ranged from X to Y million reads."
               raw_reads_text <- paste0("Individual library size ranged from ", min_raw, " million to ", max_raw, " million reads. ")
            }
          }, error = function(e) {})
        }
      }

      # Generate methods text
      methods_sections <- generateMethodsText(
        param, aligner, de_method, fastp_version, pkg_versions, se,
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
          "Generated on: ", format(Sys.Date(), "%d %B %Y"), "."
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
          officer::fpar(officer::ftext("RNA-Seq Data Analysis", prop = officer::fp_text(bold = TRUE, font.size = 14)))
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
          officer::fpar(officer::ftext("Read Trimming and Quality Control", prop = officer::fp_text(bold = TRUE, font.size = 12)))
        ) %>%
        officer::body_add_par(methods_sections$trimming) %>%
        officer::body_add_par("")

      # Add Mapping section
      doc <- doc %>%
        officer::body_add_fpar(
          officer::fpar(officer::ftext("Read Mapping and Quantification", prop = officer::fp_text(bold = TRUE, font.size = 12)))
        ) %>%
        officer::body_add_par(methods_sections$mapping) %>%
        officer::body_add_par("")



      # Add DE Analysis section
      doc <- doc %>%
        officer::body_add_fpar(
          officer::fpar(officer::ftext("Differential Expression Analysis", prop = officer::fp_text(bold = TRUE, font.size = 12)))
        ) %>%
        officer::body_add_par(methods_sections$de_analysis) %>%
        officer::body_add_par("")

      # Add GO Analysis section if applicable
      if (!is.null(methods_sections$go_analysis)) {
        doc <- doc %>%
          officer::body_add_fpar(
            officer::fpar(officer::ftext("Gene Ontology Enrichment Analysis", prop = officer::fp_text(bold = TRUE, font.size = 12)))
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
      }, error = function(e) {
        showNotification(
          paste0("Error generating methods document: ", e$message),
          type = "error",
          duration = 10
        )
        cat("Methods download error:", e$message, "\n")
      })
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
  } else if (inputDataReactive()$dataType %in% c("CountQC", "CrisprCountQC")) {
    output$settingsTable <- function() {
      param <- inputDataReactive()$param
      settings <- character()

      if (!is.null(param$refBuild)) settings["Reference:"] <- param$refBuild
      if (!is.null(param$normMethod)) settings["Normalization method:"] <- param$normMethod
      if (!is.null(param$signalThreshold)) settings["Log2 signal threshold:"] <- param$signalThreshold
      if (!is.null(param$signalThreshold)) settings["Linear signal threshold:"] <- round(2^param$signalThreshold, 2)
      if (!is.null(param$featureLevel)) settings["Feature level:"] <- param$featureLevel
      settings["Number of features:"] <- nrow(se)
      if (!is.null(param$countColumn)) settings["Data Column Used:"] <- param$countColumn

      kable(as.data.frame(settings),
        row.names = TRUE,
        col.names = "Setting", format = "html"
      ) %>%
        kable_styling(
          bootstrap_options = "striped", full_width = FALSE,
          position = "left"
        )
    }
    output$summaryTable <- function() {
      NULL
    }
  }
})
