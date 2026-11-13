# Summary tab: SUSHI settings table, tool detection, citations and the Methods text.

#' count_result_summary
#' SUSHI RNA-seq settings table for the Summary tab. A port of
#' `ezRun::makeCountResultSummary()`, so exploreDE does not need ezRun for it.
#' @param param SUSHI parameter list (`metadata(se)$param`).
#' @param se The SUSHI deResult `SummarizedExperiment`.
#' @return A one-column data.frame of settings.
#' @noRd
count_result_summary <- function(param, se) {
  specified <- function(x) !is.null(x) && length(x) > 0 && !is.na(x[1]) && nzchar(x[1])
  md <- S4Vectors::metadata(se)
  settings <- character()
  # A setting the result does not record (NULL) gets no row
  put <- function(settings, label, value) {
    if (length(value) > 0) {
      settings[label] <- value
    }
    settings
  }
  settings <- put(settings, "Analysis:", md$analysis)
  settings <- put(settings, "Reference:", param$refBuild)
  settings <- put(settings, "Feature level:", md$featureLevel)
  settings <- put(settings, "Data Column Used:", md$countName)
  settings <- put(settings, "Method:", md$method)
  if (specified(param$sampleGroupBaseline) && specified(param$refGroupBaseline)) {
    settings["Baseline correction:"] <- paste(param$sampleGroupBaseline, "and", param$refGroupBaseline)
  }
  settings <- put(settings, "Comparison:", param$comparison)
  if (!is.null(param$normMethod)) {
    settings["Normalization:"] <- param$normMethod
  }
  if (!is.null(param$deTest)) {
    settings["Differential expression test:"] <- param$deTest
  }
  if (isTRUE(param$useSigThresh)) {
    settings["Log2 signal threshold:"] <- signif(log2(param$sigThresh), digits = 4)
    settings["Linear signal threshold:"] <- signif(param$sigThresh, digits = 4)
  }
  as.data.frame(settings)
}

#' detectAligner
#' Inspect input_dataset.tsv to guess whether STAR / Kallisto / FeatureCounts was used.
#' @param dataDir Analysis directory containing ../input_dataset.tsv
#' @return Character: one of "STAR", "Kallisto", "FeatureCounts", "Unknown"
#' @noRd
detectAligner <- function(dataDir) {
  dataset_path <- file.path(dataDir, "../input_dataset.tsv")
  if (!file.exists(dataset_path)) {
    return("Unknown")
  }
  dataset <- readr::read_tsv(dataset_path, show_col_types = FALSE)

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
  "Unknown"
}

#' detectDEMethod
#' Guess DE method (DESeq2 / edgeR) from dataDir string and param fields.
#' @param dataDir Analysis directory string
#' @param param Parameter list
#' @return Character: "DESeq2" or "edgeR"
#' @noRd
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

  return("DESeq2") # Default assumption
}

#' formatRefBuild
#' Pretty-print a reference build path like "Homo_sapiens/GENCODE/GRCh38.p13/.../Release_42-..."
#' @param ref_build Character reference build path
#' @return Human-readable string (e.g. "GENCODE human genome build GRCh38.p13 (Release 42)")
#' @noRd
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

#' getFastpVersion
#' Scan preprocessing logs in `dirname(dataDir)` for the fastp version string.
#' @param dataDir Analysis directory
#' @return Version string or NULL
#' @noRd
getFastpVersion <- function(dataDir) {
  log_files <- list.files(dirname(dataDir), pattern = "_preprocessing\\.log$", full.names = TRUE)
  if (length(log_files) == 0) {
    return(NULL)
  }
  version_line <- grep("fastp v", readLines(log_files[1], n = 100, warn = FALSE), value = TRUE)
  if (length(version_line) == 0) {
    return(NULL)
  }
  sub(".*fastp v([0-9.]+).*", "\\1", version_line[1])
}

#' getPackageVersions
#' Look up installed versions of the given R packages (used for methods-text).
#' @param packages Character vector of package names
#' @return Named character vector of versions (or "version unavailable")
#' @noRd
getPackageVersions <- function(packages) {
  versions <- vapply(
    packages,
    function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else "version unavailable"
    },
    character(1)
  )
  return(versions)
}

#' buildCitationsList
#' Build a named list of citation strings for methods section of the report.
#' @param aligner Aligner name from detectAligner()
#' @param de_method DE method name from detectDEMethod()
#' @param run_go Logical; include clusterProfiler citation
#' @param ref_build Reference build path (optional; drives GENCODE citation)
#' @return Named list of citation strings
#' @noRd
buildCitationsList <- function(aligner, de_method, run_go = FALSE, ref_build = NULL) {
  citations <- list()

  # Always cite SUSHI
  citations$SUSHI <- paste0(
    "Hatakeyama M, Opitz L, Russo G, Qi W, Schlapbach R, Rehrauer H (2016) SUSHI: an exquisite recipe for fully ",
    "documented, reproducible and reusable NGS data analysis. BMC Bioinformatics 17, 228. ",
    "https://doi.org/10.1186/s12859-016-1104-8"
  )

  # FastQC
  citations$FastQC <- paste0(
    "Andrews S (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data. Available online at: ",
    "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
  )

  # MultiQC
  citations$MultiQC <- paste0(
    "Ewels P, Magnusson M, Lundin S, K\u00e4ller M (2016) MultiQC: summarize analysis results for multiple tools and ",
    "samples in a single report. Bioinformatics 32(19), 3047\u20133048. https://doi.org/10.1093/bioinformatics/btw354"
  )

  # fastp
  citations$fastp <- paste0(
    "Chen S, Zhou Y, Chen Y, Gu J (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, ",
    "i884\u2013i890. https://doi.org/10.1093/bioinformatics/bty560"
  )

  # Aligner-specific citations
  if (aligner == "STAR") {
    citations$STAR <- paste0(
      "Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013) STAR: ",
      "ultrafast universal RNA-seq aligner. Bioinformatics 29, 15\u201321. ",
      "https://doi.org/10.1093/bioinformatics/bts635"
    )
  } else if (aligner == "Kallisto") {
    citations$Kallisto <- paste0(
      "Bray NL, Pimentel H, Melsted P, Pachter L (2016) Near-optimal probabilistic RNA-seq quantification. Nature ",
      "Biotechnology 34, 525\u2013527. https://doi.org/10.1038/nbt.3519"
    )
  } else if (aligner == "FeatureCounts") {
    citations$Rsubread <- paste0(
      "Liao Y, Smyth GK, Shi W (2013) The Subread aligner: fast, accurate and scalable read mapping by ",
      "seed-and-vote. Nucleic Acids Research 41, e108. https://doi.org/10.1093/nar/gkt214"
    )
    citations$STAR <- paste0(
      "Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013) STAR: ",
      "ultrafast universal RNA-seq aligner. Bioinformatics 29, 15\u201321. ",
      "https://doi.org/10.1093/bioinformatics/bts635"
    )
  }

  # GENCODE citation if detected in reference build
  if (!is.null(ref_build) && grepl("GENCODE", ref_build, ignore.case = TRUE)) {
    citations$GENCODE <- paste0(
      "Frankish A, Diekhans M, Jungreis I, Lagarde J, Loveland JE, Mudge JM, Sisu C, Wright JC, Armstrong J, Barnes ",
      "I, et al. (2021) GENCODE: reference annotation for the human and mouse genomes. Nucleic Acids Research 49, ",
      "D916\u2013D923. https://doi.org/10.1093/nar/gkaa1087"
    )
  }

  # DE method citations
  if (de_method == "DESeq2") {
    citations$DESeq2 <- paste0(
      "Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with ",
      "DESeq2. Genome Biology 15, 550. https://doi.org/10.1186/s13059-014-0550-8"
    )
  } else if (de_method == "edgeR") {
    citations$edgeR <- paste0(
      "Robinson MD, McCarthy DJ, Smyth GK (2010) edgeR: a Bioconductor package for differential expression analysis ",
      "of digital gene expression data. Bioinformatics 26, 139\u2013140. ",
      "https://doi.org/10.1093/bioinformatics/btp616"
    )
  }

  # GO analysis citation if applicable
  if (run_go) {
    citations$clusterProfiler <- paste0(
      "Yu G, Wang LG, Han Y, He QY (2012) clusterProfiler: an R package for comparing biological themes among gene ",
      "clusters. OMICS 16, 284\u2013287. https://doi.org/10.1089/omi.2011.0118"
    )
  }

  # R and exploreDE citations
  citations$R <- paste0(
    "R Core Team (2020) R: A language and environment for statistical computing. R Foundation for Statistical ",
    "Computing, Vienna, Austria. https://www.R-project.org"
  )
  citations$exploreDE <- paste0(
    "Leary P, Rehrauer H (2024) exploreDE Interactive Shiny App. Zenodo. https://doi.org/10.5281/zenodo.13927692"
  )

  return(citations)
}

#' generateMethodsText
#' Build the named list of methods-section paragraphs (trimming, mapping, DE analysis, GO, processing).
#' @param param Parameter list from the SE
#' @param aligner Aligner name from detectAligner()
#' @param de_method DE method name from detectDEMethod()
#' @param fastp_version Version string from getFastpVersion() (or NULL)
#' @param pkg_versions Named character vector from getPackageVersions()
#' @param se SummarizedExperiment (for library-size info)
#' @param raw_reads_text Optional prefix describing raw-read counts
#' @return Named list of methods paragraphs
#' @noRd
generateMethodsText <- function(param, aligner, de_method, fastp_version, pkg_versions, se, raw_reads_text = "") {
  sections <- list()
  sections$trimming <- methods_trimming_text(raw_reads_text, fastp_version)
  sections$mapping <- paste0(methods_mapping_text(param, aligner), " ", final_library_size_text(se))
  sections$de_analysis <- methods_de_text(de_method, pkg_versions)
  if (!is.null(param$runGO) && param$runGO) {
    sections$go_analysis <- methods_go_text(param, pkg_versions)
  }
  sections$data_processing <- methods_processing_text()
  sections
}

# " (version x.y.z)" for a package of getPackageVersions(), or "".
methods_version_text <- function(pkg_versions, pkg) {
  if (pkg_versions[pkg] != "version unavailable") paste0(" (version ", pkg_versions[pkg], ")") else ""
}

# Final mapped / counted reads, from the SE's counts.
final_library_size_text <- function(se) {
  if (is.null(se) || !"counts" %in% SummarizedExperiment::assayNames(se)) {
    return("")
  }
  lib_sizes <- colSums(SummarizedExperiment::assay(se, "counts"))
  paste0(
    "Final library sizes ranged from ",
    round(min(lib_sizes) / 1e6, 1),
    " million to ",
    round(max(lib_sizes) / 1e6, 1),
    " million reads. "
  )
}

methods_trimming_text <- function(raw_reads_text, fastp_version) {
  fastp_ver_text <- if (!is.null(fastp_version)) paste0(" (version ", fastp_version, ")") else ""
  paste0(
    raw_reads_text,
    "Read quality was inspected using FastQC (Andrews, 2010) and MultiQC (Ewels et al., 2016), and sequencing ",
    "adaptors removed using fastp",
    fastp_ver_text,
    " (Chen et al., 2018). ",
    "Quality trimming was performed to remove low-quality bases, with reads shorter than 18 bases discarded. ",
    "Poly-X tails were trimmed with a minimum length threshold of 10 bases."
  )
}

methods_mapping_text <- function(param, aligner) {
  ref_build <- formatRefBuild(if (!is.null(param$refBuild)) param$refBuild else "the reference genome")
  feature_level <- if (!is.null(param$featureLevel)) param$featureLevel else "gene"
  transcript_types <- if (!is.null(param$transcriptTypes)) param$transcriptTypes else "protein-coding"
  # STAR and featureCounts share the alignment sentence and what was kept
  star <- paste0(
    "Alignment of the RNA-seq reads was performed using the STAR aligner (Dobin et al., 2013) with the ",
    ref_build,
    " reference genome as the reference. "
  )
  genes_kept <- paste0("Only ", transcript_types, " genes were retained for downstream analysis.")
  switch(
    aligner,
    Kallisto = paste0(
      "Transcript-level abundance was estimated using Kallisto (Bray et al., 2016), which performs ",
      "pseudo-alignment of reads to the ",
      ref_build,
      " reference transcriptome. Kallisto uses an expectation-maximisation algorithm to quantify transcript ",
      "abundances. Only ",
      transcript_types,
      " transcripts were retained for downstream analysis."
    ),
    STAR = paste0(
      star,
      "STAR employs a seed-and-extend approach optimised for RNA-seq data, enabling accurate detection of splice ",
      "junctions. ",
      genes_kept
    ),
    FeatureCounts = paste0(
      star,
      "The counting of ",
      feature_level,
      "-level expression values was performed using the 'featureCounts' function of the R package Rsubread ",
      "(Liao et al., 2013). ",
      genes_kept
    ),
    paste0(
      "Reads were aligned to the ",
      ref_build,
      " reference genome and quantified at the ",
      feature_level,
      " level using standard RNA-seq alignment and quantification methods, with only ",
      transcript_types,
      " genes retained for downstream analysis."
    )
  )
}

# Filter, normalisation and DE; NULL for a method other than DESeq2 or edgeR.
methods_de_text <- function(de_method, pkg_versions) {
  switch(
    de_method,
    DESeq2 = paste0(
      "Count data were normalised using DESeq2's median-of-ratios method, which accounts for differences in library ",
      "size and RNA composition bias between samples. ",
      "Differential expression was performed using the DESeq2 Bioconductor R package",
      methods_version_text(pkg_versions, "DESeq2"),
      " (Love et al., 2014), which implements a generalised linear model with a negative binomial distribution ",
      "to account for technical and biological variability. ",
      "Dispersion and log2 fold change estimates were moderated using shrinkage methods to improve stability. ",
      "Statistical significance was assessed using the Wald test, with p-values adjusted for multiple testing ",
      "using the Benjamini-Hochberg method to control the false discovery rate (FDR)."
    ),
    edgeR = paste0(
      "Count data were normalised using the trimmed mean of M-values (TMM) method implemented in edgeR, ",
      "which corrects for differences in library size and RNA composition between samples. ",
      "Differential expression was performed using the edgeR Bioconductor R package",
      methods_version_text(pkg_versions, "edgeR"),
      " (Robinson et al., 2010), which implements a generalised linear model with a negative binomial distribution ",
      "to account for technical and biological variability. ",
      "Dispersion estimates were moderated using empirical Bayes shrinkage to improve reliability for genes with low ",
      "counts. Statistical significance was assessed using the likelihood ratio test (LRT), ",
      "with p-values adjusted for multiple testing using the Benjamini-Hochberg method to control the false ",
      "discovery rate (FDR)."
    )
  )
}

methods_go_text <- function(param, pkg_versions) {
  paste0(
    "Gene Ontology (GO) term pathway analysis was performed using both the hypergeometric over-representation test ",
    "via the 'enricher' function, and gene-set enrichment analysis via the 'GSEA' function, ",
    "of the clusterProfiler Bioconductor R package",
    methods_version_text(pkg_versions, "clusterProfiler"),
    " (Yu et al., 2012), with FDR thresholds of ",
    if (!is.null(param$fdrThreshORA)) param$fdrThreshORA else 0.05,
    " and ",
    if (!is.null(param$fdrThreshGSEA)) param$fdrThreshGSEA else 0.05,
    ", respectively."
  )
}

methods_processing_text <- function() {
  bioc_version <- if (requireNamespace("BiocManager", quietly = TRUE)) {
    paste0(" and Bioconductor version ", as.character(BiocManager::version()))
  } else {
    ""
  }
  paste0(
    "All analyses were conducted using the SUSHI framework (Hatakeyama et al., 2016), ",
    "which provides a reproducible and documented workflow for NGS data analysis. ",
    "All R functions were executed on R version ",
    R.version$major,
    ".",
    R.version$minor,
    bioc_version,
    " (R Core Team, 2020). ",
    "Figures were generated using the exploreDE Shiny application using the data generated from the above steps ",
    "(Leary and Rehrauer, 2024)."
  )
}
