# Helper functions

#' batch_correction_design
#' Build a design matrix for limma::removeBatchEffect that preserves biological
#' signal. Never throws — returns a tagged result so callers can notify the
#' user and decide whether to skip batch correction or proceed without a
#' design matrix.
#' @param dataset Sample metadata data.frame
#' @param designColumn Column name to build the design from. "None" / empty
#'   means "no design" (batch correction still runs, without preserving any
#'   biology).
#' @param batchColumn Column name of the batch factor — required for the
#'   collinearity check.
#' @return list(design = <matrix or NULL>, skip = <logical>, msg = <chr or NULL>)
#' @export
batch_correction_design <- function(dataset, designColumn, batchColumn) {
  null_design <- list(design = NULL, skip = FALSE, msg = NULL)
  skipped     <- function(msg) list(design = NULL, skip = TRUE,  msg = msg)
  warned      <- function(msg) list(design = NULL, skip = FALSE, msg = msg)

  if (is.null(designColumn) || !nzchar(designColumn) || designColumn == "None") {
    return(null_design)
  }
  if (!designColumn %in% colnames(dataset)) {
    return(skipped(sprintf(
      "Design column '%s' missing from sample metadata — batch correction skipped.",
      designColumn
    )))
  }

  grp <- factor(dataset[[designColumn]])
  if (nlevels(grp) < 2) {
    return(warned(sprintf(
      "Design factor '%s' has only one level — batch correction applied without a design matrix.",
      designColumn
    )))
  }

  design <- tryCatch(stats::model.matrix(~grp), error = function(e) e)
  if (inherits(design, "error")) {
    return(skipped(sprintf(
      "Could not build design matrix from '%s': %s — batch correction skipped.",
      designColumn, conditionMessage(design)
    )))
  }

  if (!is.null(batchColumn) && nzchar(batchColumn) && batchColumn != "None" &&
      batchColumn %in% colnames(dataset)) {
    batch <- factor(dataset[[batchColumn]])
    if (nlevels(batch) >= 2) {
      batch_mm <- tryCatch(stats::model.matrix(~batch), error = function(e) NULL)
      if (!is.null(batch_mm)) {
        combined <- cbind(batch_mm, design[, -1, drop = FALSE])
        if (qr(combined)$rank < ncol(combined)) {
          return(warned(sprintf(
            "Batch ('%s') and design ('%s') are confounded — batch correction applied without a design matrix.",
            batchColumn, designColumn
          )))
        }
      }
    }
  }

  list(design = design, skip = FALSE, msg = NULL)
}

#' apply_batch_correction
#' Graceful-failure wrapper around limma::removeBatchEffect. Never throws:
#' on any validation or limma error, returns the original counts unchanged
#' and surfaces a Shiny notification.
#' @export
apply_batch_correction <- function(counts, dataset, batchCol, designCol) {
  if (is.null(batchCol) || !nzchar(batchCol) || batchCol == "None") return(counts)

  notify <- function(msg, type = "warning") {
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::showNotification(msg, type = type, duration = 8)
    } else {
      warning(msg, call. = FALSE)
    }
  }

  if (!batchCol %in% colnames(dataset) ||
      nlevels(factor(dataset[[batchCol]])) < 2) {
    notify(sprintf("Batch factor '%s' has <2 levels — batch correction skipped.", batchCol))
    return(counts)
  }

  info <- batch_correction_design(dataset, designCol, batchCol)
  if (!is.null(info$msg)) {
    notify(info$msg, type = if (isTRUE(info$skip)) "error" else "warning")
  }
  if (isTRUE(info$skip)) return(counts)

  args <- list(x = counts, batch = dataset[[batchCol]])
  if (!is.null(info$design)) args$design <- info$design

  tryCatch(
    do.call(limma::removeBatchEffect, args),
    error = function(e) {
      notify(sprintf("Batch correction failed: %s — showing uncorrected counts.",
                     conditionMessage(e)),
             type = "error")
      counts
    }
  )
}

#' standardise_columns
#' Fuzzy-match incoming column names to canonical names.
#' Upstream data may call columns anything; this maps them to a known schema
#' so downstream code can rely on consistent names.
#' @param df A data.frame (typically from rowData contrast)
#' @return The same data.frame with standardised column names
standardise_columns <- function(df) {
  # Canonical name → regex patterns (case-insensitive) that map to it
  patterns <- c(
    "protein_Id"  = "^protein[._]?id$",
    "peptide_Id"  = "^peptide[._]?id$|^peptide$",
    "site"        = "^site$",
    "gene_name"   = "^gene[._]?name$|^gene[._]?symbol$",
    "description" = "^description$|^desc$",
    "fasta.id"    = "^fasta[._]?id$",
    "modelName"   = "^model[._]?name$",
    "nrPeptides"  = "^nr[._]?peptides$|^n[._]?peptides$|^num[._]?peptides$",
    "log2Ratio"   = "^diff$|^log2[._]?(fold)?[._]?(change|ratio|fc)?$|^logfc$",
    "pValue"      = "^p[._]?value$",
    "fdr"         = "^fdr$|^adj[._]?p[._]?value$|^padj$|^q[._]?value$"
  )

  current <- colnames(df)
  for (canonical in names(patterns)) {
    idx <- grep(patterns[[canonical]], current, ignore.case = TRUE)
    if (length(idx) == 1) {
      current[idx] <- canonical
    }
  }
  colnames(df) <- current
  df
}

#' build_feature_id
#' Build a unique, human-readable feature identifier for proteomics data.
#' Uses canonical column names (after standardisation).
#' @param df A standardised data.frame with canonical column names
#' @return Character vector of unique feature IDs
build_feature_id <- function(df) {
  pid <- df$protein_Id

  # If protein_Id is already unique, use it directly
  if (!anyDuplicated(pid)) return(pid)

  # Try known sub-feature columns (canonical names) in priority order
  # Always prepend protein_Id so users see which protein a sub-feature belongs to
  for (col in intersect(c("site", "peptide_Id"), colnames(df))) {
    combined <- paste(pid, gsub("~", "_", df[[col]]), sep = "_")
    if (!anyDuplicated(combined)) return(combined)
  }

  # Fallback: cleaned raw rownames
  gsub("~[^~]+~", "_", rownames(df))
}

#' collapse_guides_to_genes
#' Collapse sgRNA-level count matrix to gene-level using row sums.
#' Filters out negative control genes (matching ^negCtrl).
#' @param counts A numeric matrix with sgRNA rows and sample columns
#' @param sgRNA_to_gene Named character vector: names are sgRNA feature IDs, values are gene names
#' @return A numeric matrix with gene-level rows
collapse_guides_to_genes <- function(counts, sgRNA_to_gene, exclude = NULL) {
  common <- intersect(rownames(counts), names(sgRNA_to_gene))
  counts <- counts[common, , drop = FALSE]
  collapsed <- rowsum(counts, group = sgRNA_to_gene[common], reorder = FALSE)
  if (length(exclude) > 0) {
    collapsed <- collapsed[!rownames(collapsed) %in% exclude, , drop = FALSE]
  }
  collapsed
}

#' convert_genomics_se
#' @export
convert_genomics_se <- function(se, dataDir, param = NULL) {
  # Load the parameters
  if (is.null(param)) {
    param <- metadata(se)$param
  }

  # Generic for RNA seq AND CountQC ----
  # Prepare the seqAnno results table
  # Suppress S4Vectors deprecation warnings when accessing genomic ranges
  seqAnno <- suppressWarnings(data.frame(rowData(se), row.names = rownames(se), check.names = FALSE))
  # Tidy the description column to remove extraneous accessors
  if ("description" %in% colnames(seqAnno)) {
    seqAnno$description <- gsub(" \\[..*", "", seqAnno$description)
  } else {
    seqAnno$description <- seqAnno$gene_name
  }
  # Uniquify the gene symbols to create feature_id
  if (param$featureLevel == "gene") {
    seqAnno$feature_id <- scuttle::uniquifyFeatureNames(
      ID = seqAnno$gene_id,
      names = seqAnno$gene_name
    )
  } else if (param$featureLevel == "isoform") {
    seqAnno$feature_id <- scuttle::uniquifyFeatureNames(
      ID = seqAnno$transcript_id,
      names = seqAnno$gene_name
    )
  }
  # If uniquifying broke some names, fix them
  seqAnno$feature_id[grep("^_", seqAnno$feature_id)] <- gsub("_", "", seqAnno$feature_id[grep("^_", seqAnno$feature_id)])
  genes <- seqAnno$feature_id

  # Prepare the metadata table ----
  if (file.exists(file.path(dataDir, "../input_dataset.tsv"))) {
    dataset <- read_tsv(file.path(dataDir, "../input_dataset.tsv"), show_col_types = F) %>% data.frame(check.names = F)
  } else {
    dataset <- data.frame(colData(se), check.names = FALSE)
  }
  if (length(param$samples) > 0) {
    dataset <- dataset[which(dataset$Name %in% param$samples), ]
  }
  dataset[["Name"]] <- factor(dataset[["Name"]], levels = dataset[["Name"]])
  rownames(dataset) <- dataset$Name
  # Get the factors that should be available for plotting
  factors <- colnames(dataset)[grep("\\[Factor]", colnames(dataset))]
  # Remove [Factor] from colnames and factors
  colnames(dataset) <- gsub(" \\[.*", "", colnames(dataset)) %>% gsub(" ", "_", .)
  factors <- gsub(" \\[.*", "", factors) %>% gsub(" ", "_", .)
  # Get the individual levels for each factor
  factorLevels <- lapply(factors, function(f) {
    levels(as.factor(dataset[[f]]))
  }) %>% setNames(factors)
  # If someone manually enters a duplicate column (but sans [Factor]), things will break. So, remove!
  # Identify duplicated column names
  dups <- which(duplicated(names(dataset)) | duplicated(names(dataset), fromLast = TRUE))
  dup_names <- unique(names(dataset)[dups])
  # Go through each duplicate group
  for (dup in dup_names) {
    cols <- which(names(dataset) == dup)
    # Compare all pairwise combinations (in case there are >2)
    keep <- cols[1]
    drop_these <- integer(0)
    for (i in 2:length(cols)) {
      if (all(dataset[[cols[i]]] == dataset[[keep]], na.rm = TRUE)) {
        # If identical, mark for dropping
        drop_these <- c(drop_these, cols[i])
      } else {
        # If not identical, make unique
        names(dataset)[cols[i]] <- paste0(dup, ".", i - 1)
      }
    }
    # Drop identical duplicates
    if (length(drop_these) > 0) {
      dataset <- dataset[, -drop_these, drop = FALSE]
    }
  }
  # The app doesn't like empty strings in the factor columns, so for now, add "None" to those entries
  # dataset <- dataset %>% mutate_all(~ifelse(. == "", "None", .))
  if (is.null(param$groupingName) || !nzchar(param$groupingName)) {
    param$groupingName <- factors[[1]]
  }
  # Prepare the count files
  countList <- list()
  countList[["Raw"]] <- assay(se, "counts")
  normCounts <- grep("^signal$|^xNorm$", names(assays(se)), value = TRUE)
  countList[["Normalised"]] <- assay(se, normCounts)
  countList[["Normalised + Log2"]] <- log2(assay(se, normCounts) + param$backgroundExpression)
  if (!is.null(param$refBuild)) {
    countList[["FPKM"]] <- getRpkm(se)
    countList[["TPM"]] <- getTpm(se)
  }
  forVST <- as.matrix(assay(se, "counts"))
  mode(forVST) <- "integer"
  vstCountMatrix <- tryCatch({
    if (nrow(forVST) < 1000) {
      DESeq2::varianceStabilizingTransformation(forVST, blind = FALSE)
    } else {
      DESeq2::vst(forVST, blind = FALSE)
    }
  }, error = function(e) {
    message("VST failed (", e$message, "), falling back to varianceStabilizingTransformation")
    DESeq2::varianceStabilizingTransformation(forVST, blind = FALSE)
  })
  countList[["VST"]] <- vstCountMatrix
  # For all the counts, match the assays to the seqAnno order
  for (i in seq_along(countList)) {
    if (param$featureLevel == "gene") {
      geneOrder <- intersect(seqAnno$gene_id, rownames(countList[[i]]))
      countList[[i]] <- countList[[i]][geneOrder, ]
      rownames(countList[[i]]) <- seqAnno$feature_id
    } else if (param$featureLevel == "isoform") {
      geneOrder <- intersect(seqAnno$transcript_id, rownames(countList[[i]]))
      countList[[i]] <- countList[[i]][geneOrder, ]
      rownames(countList[[i]]) <- seqAnno$feature_id %>%
        gsub("[:/()-]", ".", .) %>%
        gsub(" ", ".", .)
    }
  }
  # Prepare the GO info ----
  allPathwaysDEG <- c(
    unique(unlist(strsplit(x = c(seqAnno$`GO BP`), "; "))),
    unique(unlist(strsplit(x = c(seqAnno$`GO MF`), "; "))),
    unique(unlist(strsplit(x = c(seqAnno$`GO CC`), "; ")))
  )
  goTerms <- clusterProfiler::go2term(allPathwaysDEG)
  allPathways <- paste(goTerms[, 1], goTerms[, 2])

  dataTypeGStore <- if (grepl("ExploreMageckCounts", dataDir)) {
    "CrisprCountQC"
  } else if (grepl("CountQC", dataDir)) {
    "CountQC"
  } else {
    "RNASeq"
  }

  # Create the first part of the import list
  gStoreList <- list(
    "dataType" = dataTypeGStore,
    "se" = se,
    "param" = param,
    "seqAnno" = seqAnno,
    "dataset" = dataset,
    "factors" = factors,
    "factorLevels" = factorLevels,
    "groupingName" = param$groupingName,
    "countList" = countList,
    "genes" = genes,
    "allPathways" = allPathways
  )

  # DE results ----
  if (dataTypeGStore %in% c("CountQC", "CrisprCountQC")) {
    # CountQC/CrisprCountQC don't have DE results, but some tabs expect these fields
    gStoreList$design <- dataTypeGStore
    gStoreList$designLevels <- NULL
    if (dataTypeGStore == "CrisprCountQC") {
      # Extract gene names from sgRNA IDs
      # Handles both GENE_N_SEQUENCE and GENE_sgrnaN_SEQUENCE patterns
      gStoreList$sgRNA_to_gene <- setNames(
        gsub("_(sgrna)?[0-9]+_[ACGT]+$", "", seqAnno$feature_id, ignore.case = TRUE),
        seqAnno$feature_id
      )
    }
  } else if (dataTypeGStore == "RNASeq") {
    # Get the p-value and log2 columns
    pValueColumn <- grep("^p[-_]?value$", colnames(seqAnno), value = TRUE, ignore.case = TRUE)[[1]]
    pAdjColumn <- grep("^(padj|adj[._]?p[-_]?value|adjusted[._]?p[-_]?value|fdr|q[-_]?value)$", colnames(seqAnno), value = TRUE, ignore.case = TRUE)[[1]]
    log2FCColumn <- grep("^(log2?[._ ]?(fold)?[._ ]?(change|ratio|fc)|logfc)$", colnames(seqAnno), value = TRUE, ignore.case = TRUE)[[1]]
    # Sort the dataframe by raw pvalue
    seqAnno <- seqAnno[order(seqAnno[[pValueColumn]]), ]
    # get the list of genes to go in as input options:
    genes <- seqAnno$feature_id[order(seqAnno$pValue)]
    # Set the comparison column levels to be in order of test, e.g., diff, undiff
    dataset[[param$groupingName]] <- factor(
      dataset[[param$groupingName]],
      levels = c(
        param$sampleGroup,
        param$refGroup,
        unique(dataset[[param$groupingName]])[!unique(dataset[[param$groupingName]]) %in% c(param$sampleGroup, param$refGroup)]
      )
    )
    # Add global and per group means to the seqAnno table
    for (count in names(countList)) {
      cts <- countList[[count]]
      seqAnno[[paste0(count, "_mean")]] <- rowMeans(cts)[seqAnno$feature_id]
      seqAnno[[paste0(count, "_SD")]] <- matrixStats::rowSds(cts)[seqAnno$feature_id]
      isSample <- dataset[[param$groupingName]] == param$sampleGroup
      seqAnno[[paste0(count, "_", param$sampleGroup, "_Mean")]] <- rowMeans(cts[, isSample, drop = FALSE])[seqAnno$feature_id]
      isRef <- dataset[[param$groupingName]] == param$refGroup
      seqAnno[[paste0(count, "_", param$refGroup, "_Mean")]] <- rowMeans(cts[, isRef, drop = FALSE])[seqAnno$feature_id]
    }
    # Prepare the contrast info ----
    design <- param$comparison
    designLevels <- design %>%
      strsplit(split = "--over--") %>%
      .[[1]]

    goObjects <- list(runGO = param$runGO)
    if (param$runGO) {
      goObjects$ora <- metadata(se)$enrichResult
      goObjects$gsea <- metadata(se)$enrichResultGSEA
    }
    groupingName <- param$groupingName
    for (f in c("goObjects", "groupingName", "design", "designLevels", "seqAnno", "dataset", "genes")) {
      gStoreList[[f]] <- get(f)
    }
  }

  # Store dataDir as attribute for methods text generation
  attr(gStoreList, "dataDir") <- dataDir

  return(
    gStoreList
  )
}

#' convert_proteomics_se
#' @export
convert_proteomics_se <- function(se) {
  # Get the contrasts available
  contrasts <- names(rowData(se))[grep("^constrast_", names(rowData(se)))] %>% gsub("constrast_", "", .)

  # Convert contrast DataFrames once, reuse for seqAnnoList and countList lookup
  contrast_dfs <- setNames(lapply(contrasts, function(con) {
    as.data.frame(rowData(se)[[paste0("constrast_", con)]])
  }), contrasts)

  # Make a list of seqAnno tables for each contrast
  seqAnnoList <- setNames(lapply(contrasts, function(con) {
    sa <- standardise_columns(contrast_dfs[[con]])
    sa <- sa[order(sa$pValue), ]

    # Build feature_id from canonical column names
    sa$feature_id <- build_feature_id(sa)

    # If feature_id == protein_Id (protein-level data), drop redundant protein_Id
    if ("protein_Id" %in% colnames(sa) && identical(sa$feature_id, sa$protein_Id)) {
      sa$protein_Id <- NULL
    }

    # Keep feature_id, all non-numeric columns (domain/annotation), and selected stats
    stats_to_keep <- c("log2Ratio", "pValue", "fdr", "nrPeptides")
    non_numeric <- colnames(sa)[!vapply(sa, is.numeric, logical(1))]
    keep <- unique(c("feature_id", non_numeric, intersect(stats_to_keep, colnames(sa))))
    sa <- sa[, keep, drop = FALSE]

    # Guarantee description exists (volcano, brush tables, etc. expect it)
    if (!"description" %in% colnames(sa)) {
      sa$description <- if ("gene_name" %in% colnames(sa)) sa$gene_name else sa$feature_id
    }

    sa
  }), contrasts)

  # Prepare the metadata
  dataset <- as.data.frame(colData(se))
  names <- colnames(dataset)[grep("^Sample$|^Name$|^sampleName$", colnames(dataset), ignore.case = T)]
  dataset$Name <- dataset[[names]]
  rownames(dataset) <- dataset[[names]]

  # Get the main columns from the metadata that are plottable
  factors <- colnames(dataset)[grep("Sample|Name|names|channel|raw\\.file", colnames(dataset), invert = T, ignore.case = T)]
  factorLevels <- lapply(factors, function(f) {
    levels(as.factor(dataset[[f]]))
  }) %>% setNames(factors)

  # Sort numeric factor levels
  for (x in factors) {
    dataset[[x]] <- factor(dataset[[x]], levels = unique(gtools::mixedsort(dataset[[x]])))
  }

  # Map raw SE rownames → clean feature_ids via lookup (no quant-type regex)
  # NB: seqAnnoList is sorted by pValue, so we must build the lookup from the
  # UNSORTED contrast_df to keep rowname→feature_id alignment correct.
  unsorted_sa <- standardise_columns(contrast_dfs[[contrasts[1]]])
  unsorted_sa$feature_id <- build_feature_id(unsorted_sa)
  raw_rn <- rownames(contrast_dfs[[contrasts[1]]])
  fid_lookup <- setNames(unsorted_sa$feature_id, raw_rn)

  countList <- setNames(lapply(names(assays(se)), function(ass) {
    x <- assay(se, ass)
    mapped <- fid_lookup[rownames(x)]
    valid <- !is.na(mapped)
    x <- x[valid, , drop = FALSE]
    rownames(x) <- unname(mapped[valid])
    x
  }), names(assays(se)))

  # Get the feature IDs
  genes <- unique(unlist(lapply(seqAnnoList, `[[`, "feature_id")))

  # Build lookup maps for feature selection modes
  # Use first contrast's seqAnno as representative (all contrasts share the same features)
  sa1 <- seqAnnoList[[1]]

  # gene_lookup: gene_name → feature_ids (only if gene_name column exists and differs from feature_id)
  gene_lookup <- NULL
  gene_choices <- NULL
  if ("gene_name" %in% colnames(sa1) && !identical(sa1$feature_id, sa1$gene_name)) {
    gene_lookup <- split(sa1$feature_id, sa1$gene_name)
    gene_choices <- unique(sa1$gene_name)
  }

  # protein_lookup: protein_Id → feature_ids (only for sub-protein data where protein_Id still exists)
  protein_lookup <- NULL
  protein_choices <- NULL
  if ("protein_Id" %in% colnames(sa1)) {
    protein_lookup <- split(sa1$feature_id, sa1$protein_Id)
    protein_choices <- unique(sa1$protein_Id)
  }

  # Return the final list
  return(list(
    "dataType" = "proteomics",
    "se" = se,
    "seqAnnoList" = seqAnnoList,
    "contrasts" = contrasts,
    "dataset" = dataset,
    "factors" = factors,
    "countList" = countList,
    "genes" = genes,
    "factorLevels" = factorLevels,
    "groupingName" = factors[[1]],
    "gene_lookup" = gene_lookup,
    "protein_lookup" = protein_lookup,
    "gene_choices" = gene_choices,
    "protein_choices" = protein_choices
  ))
}

#' generateParams
#' @export
generateParams <- function(
  sampleGroup = "diff", refGroup = "undiff", groupingName = "Condition",
  comparison = "diff--over--undiff", featureLevel = "gene",
  runGO = TRUE, fdrThreshORA = 0.05, fdrThreshGSEA = 0.05, pValThreshGO = 0.01,
  log2RatioThreshGO = 0, backgroundExpression = 10, refBuild = "GRCh38"
) {
  return(list(
    sampleGroup = sampleGroup,
    refGroup = refGroup,
    groupingName = groupingName,
    comparison = comparison,
    featureLevel = featureLevel,
    runGO = runGO,
    fdrThreshORA = fdrThreshORA,
    fdrThreshGSEA = fdrThreshGSEA,
    pValThreshGO = pValThreshGO,
    log2RatioThreshGO = log2RatioThreshGO,
    backgroundExpression = backgroundExpression,
    refBuild = refBuild
  ))
}

# Extract the dataset path from the URL string
queryList <- parseQueryString(session$clientData$url_search)
if (is.list(queryList)) {
  dataUrl <- queryList$data
} else {
  dataUrl <- NULL
}

# whoami?
message(ezTime(), " app:exploreDE; ", "username:", username, "; ", "dataUrl:", dataUrl)

# complete the full file path for both proteomics and genomics servers
if (!is.null(dataUrl)) {
  if (grepl("Proteomics|prolfqua", dataUrl, ignore.case = TRUE)) {
    dataDir <- paste0("https://fgcz-ms.uzh.ch/public/pStore/", dataUrl)
    dataType <- "proteomics"
  } else if (grepl("CountQC|DESeq2|EdgeR|Limma|ExploreMageckCounts", dataUrl, ignore.case = TRUE)) {
    urlDataRoot <- c("/srv/gstore/projects", "/srv/GT/analysis/course_sushi/public/gstore/projects")
    dataDir <- file.path(urlDataRoot, dataUrl)
    dataDir <- dataDir[file.exists(dataDir)][1]
    if (!file.exists(dataDir)) {
      showModal(modalDialog(
        title = "Something went wrong",
        "It looks like either the dataset you're looking for doesn't exist, or has not finished being processed in SUSHI yet."
      ))
      Sys.sleep(30)
      stopApp(returnValue = invisible())
    }
    dataType <- "gStore"
  }
  projectFromUrl <- regmatches(dataUrl, regexec("p[0-9]{4,}", dataUrl))[[1]][1] # JLR 2025
} else if (is.null(dataUrl)) {
  # dataDir <- "/srv/gstore/projects/p3009/o5638_DESeq2_diff--over--undiff_2024-11-20--12-29-40/diff--over--undiff/"
  # dataDir <- "/srv/gstore/projects/p36911/o39664_o37783_o39124_o38710_o41353_ExploreMageckCounts_2026-04-09--15-07-24/ExploreMageckCounts/counts-yuunbieawvpw-EzResult.RData"
  dataDir <- "/srv/gstore/projects/p29779/o41056_ExploreMageckCounts_2026-04-10--09-29-57/ExploreMageckCounts/counts-fzvwfvxsycwz-EzResult.RData"
  dataType <- "gStore"
  # dataDir <- "https://fgcz-ms.uzh.ch/public/pStore/p39819/bfabric/Proteomics/exploreDE/2025/2025-10/2025-10-22/workunit_335485/2997996.rds"
  # dataType <- "proteomics"
  projectFromUrl <- dataDir %>%
    gsub(".*projects\\/", "", .) %>%
    gsub("\\/.*", "", .)
  showNotification("Since you did not specify a dataset in the URL, you are seeing a demo dataset.", type = "message", duration = 600, closeButton = TRUE)
}

# Store project directory and current dataset path for multi-analysis browsing
if (!is.null(projectFromUrl) && projectFromUrl != "") {
  # Determine which root path contains this project
  projectPath <- file.path("/srv/gstore/projects", projectFromUrl)
  if (dir.exists(projectPath)) {
    multiProjectInfo$projectDir <- projectPath
  } else {
    # Try course directory
    projectPath <- file.path("/srv/GT/analysis/course_sushi/public/gstore/projects", projectFromUrl)
    if (dir.exists(projectPath)) {
      multiProjectInfo$projectDir <- projectPath
    }
  }
}

# Store current dataset path (will be set after dataDir is determined)
multiProjectInfo$currentDataDir <- NULL

if (!exists("dataDir")) {
  showModal(modalDialog(
    title = "Something went wrong",
    "It looks like either the dataset you're looking for doesn't exist, or has not finished being processed in SUSHI yet."
  ))
  stop()
}
is_url <- function(dataDir) {
  return(grepl("^https?://", dataDir))
}

# Import proteomics data from pStore ----
if (dataType == "proteomics") {
  if (is_url(dataDir)) {
    se <- readRDS(url(dataDir)) # If it's a URL
  } else {
    se <- readRDS(dataDir) # If it's a local file path
  }
  # Update the available contrasts from the SE
  contrasts <- names(rowData(se))[grep("^constrast_", names(rowData(se)))] %>% gsub("constrast_", "", .)
  output$proteomicsContrastSelectorUI <- renderUI({
    selectInput(inputId = "contrastSelected", label = "Select contrast to view", choices = contrasts, selected = contrasts[1], multiple = F, selectize = T)
  })
  outputOptions(output, "proteomicsContrastSelectorUI", suspendWhenHidden = FALSE)
} else if (grepl("zip", dataDir)) {
  myTempFile <- tempfile(fileext = ".rds", tmpdir = ".") #
  res <- system(
    paste("ssh fgcz-ms.uzh.ch unzip -c -qq /srv/www/htdocs/p30385/bfabric/Proteomics/DEA_FragPipe-DiaNN/2023/2023-09/2023-09-22/workunit_293901/2360252.zip '\\*rds' > ", myTempFile)
  )
  if (res != 0) {
    message("download failed")
  }
  se <- readRDS(myTempFile)

  # Update the available contrasts from the SE
  contrasts <- names(rowData(se))[grep("^constrast_", names(rowData(se)))] %>% gsub("constrast_", "", .)
  output$proteomicsContrastSelectorUI <- renderUI({
    selectInput(inputId = "contrastSelected", label = "Select contrast to view", choices = contrasts, selected = contrasts[1], multiple = F, selectize = T)
  })
  outputOptions(output, "proteomicsContrastSelectorUI", suspendWhenHidden = FALSE)
}


# Import RNA-seq/CountQC data from SUSHI ----
if (dataType %in% c("gStore", "demo")) {
  # Run LDAP queries ----
  # Check for demo data first (no auth required)
  if (dataType == "demo") {
    message("Loading demo data, so no LDAP checks required.")
    allowed <- TRUE
    ldap_role <- "demo"
    # Check for course data (no auth required)
  } else if (grepl("/srv/GT/analysis/course_sushi/public/gstore/projects", dataDir)) {
    message("Loading SUSHI course data, so no LDAP checks required.")
    allowed <- TRUE
    ldap_role <- "course"
  } else {
    # Only run LDAP checks for gStore data that's not demo or course
    roles <- system(paste0("ldapsearch -x -H ldaps://fgcz-bfabric-ldap:636 -b 'dc=bfabric,dc=org' '(cn=", username, ")' memberof | grep Roles | sed 's/,ou=.*//g;s,.*cn=,,g'"), intern = T)
    allowedProjects <- system(paste0("ldapsearch -x -H ldaps://fgcz-bfabric-ldap:636 -b 'dc=bfabric,dc=org' '(cn=", username, ")' memberof | grep Projec | sed 's/,ou=.*//g;s,.*cn=P_,p,g' | sort | uniq"), intern = T)

    # Determine ldap_role and allowed status
    if ("R_2" %in% roles) {
      ldap_role <- "employee"
      allowed <- TRUE
    } else if ("R_3" %in% roles) {
      ldap_role <- "user"
      allowed <- ifelse(projectFromUrl %in% allowedProjects, TRUE, FALSE)
    } else {
      ldap_role <- "unauthorized"
      allowed <- FALSE
    }
    message("Role: ", roles, "; LDAP: ", ldap_role, "; Allowed to all projects?: ", allowed, "; Allowed projects: ", allowedProjects)
  }

  if (allowed) {
    # If RNA Seq, load that
    if (grepl("DESeq2|EdgeR|Limma", dataDir, ignore.case = TRUE)) {
      # Check if it's an old format URL, if so, reformat
      if (grepl("EzResult.RData", dataDir)) {
        dataDir <- gsub("\\/result-.*.-EzResult.RData", "", dataDir)
      }
      # Then load the SE
      if (file.exists(file.path(dataDir, "deResult.rds"))) {
        se <- readRDS(file.path(dataDir, "deResult.rds"))
      } else if (file.exists(file.path(dataDir, "deResult.qs2"))) {
        se <- qs2::qs_read(file.path(dataDir, "deResult.qs2"), nthreads = 4)
      }
      # Store current dataset path for multi-analysis
      multiProjectInfo$currentDataDir <- dataDir
      # If CountQC, load that (we use only the 'present' genes because it has the normalised counts in the assays)
    } else if (grepl("CountQC|ExploreMageckCounts", dataDir)) {
      dataDir <- dirname(dataDir)
      if (file.exists(file.path(dataDir, "rawData_final.rds"))) {
        se <- readRDS(file.path(dataDir, "rawData_final.rds"))
      } else if (file.exists(file.path(dataDir, "rawData_final.qs2"))) {
        se <- qs2::qs_read(file.path(dataDir, "rawData_final.qs2"), nthreads = 4)
      }
    } else {
      # Or tell the user their URL is wrong
      showModal(modalDialog(
        title = "The file does not exist",
        "Either the analysis has not yet finished running, you have made a mistake in the URL, or you have not pointed to any dataset. Please try again! If the issue persists, email peter.leary@uzh.ch",
        easyClose = TRUE,
        footer = NULL
      ))
    }
  }
}


# Generate inputDataReactive ----
inputDataReactive <- reactive({
  waiter <- waiter::Waiter$new(fadeout = TRUE, color = "#527085")
  waiter$show()
  on.exit(waiter$hide())

  # Load proteomics data ----
  if (dataType == "proteomics") {
    return(convert_proteomics_se(se))
  }

  # Load RNA Seq data ----
  if (dataType %in% c("demo", "gStore")) {
    return(convert_genomics_se(se, dataDir))
  }
})
inputDataReactive()$dataType

figuresDataReactive <- reactiveValues(
  "volcanoStatic" = NULL, "volcanoMA" = NULL, "volcanoPlotly" = NULL,
  "heatmapDEFeatures" = NULL, "heatmapCustom" = NULL, "heatmapGO" = NULL,
  "pcaStatic" = NULL, "pcaPaired" = NULL,
  "boxplotStatic" = NULL, "barplotStatic" = NULL,
  "goTilePlot" = NULL
)
# genesReactive removed — replaced by featureBucket in server-featureBucket.R

output$test <- renderText({
  inputDataReactive()$dataType
})
