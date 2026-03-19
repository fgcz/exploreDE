#' standardise_columns
#' Fuzzy-match incoming column names to canonical names.
#' Upstream data may call columns anything; this maps them to a known schema
#' so downstream code can rely on consistent names.
#' @param df A data.frame (typically from rowData contrast)
#' @return The same data.frame with standardised column names
#' @export
standardise_columns <- function(df) {
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
    if (length(idx) == 1) current[idx] <- canonical
  }
  colnames(df) <- current
  df
}

#' build_feature_id
#' Build a unique, human-readable feature identifier for proteomics data.
#' Uses canonical column names (after standardisation).
#' @param df A standardised data.frame with canonical column names
#' @return Character vector of unique feature IDs
#' @export
build_feature_id <- function(df) {
  pid <- df$protein_Id
  if (!anyDuplicated(pid)) return(pid)

  for (col in intersect(c("site", "peptide_Id"), colnames(df))) {
    combined <- paste(pid, gsub("~", "_", df[[col]]), sep = "_")
    if (!anyDuplicated(combined)) return(combined)
  }

  gsub("~[^~]+~", "_", rownames(df))
}

#' convert_genomics_se
#' @export
convert_genomics_se <- function(se, dataDir, param = NULL) {

  if (is.null(param)) {
    param <- metadata(se)$param
  }

  # Prepare the seqAnno results table ----
  seqAnno <- data.frame(rowData(se), row.names = rownames(se), check.names = FALSE)
  # Get the p-value and log2 columns
  pValueColumn <- grep("^p[-_]?value$", colnames(seqAnno), value = TRUE, ignore.case = TRUE)[[1]]
  pAdjColumn <- grep("^(padj|adj[._]?p[-_]?value|adjusted[._]?p[-_]?value|fdr|q[-_]?value)$", colnames(seqAnno), value = TRUE, ignore.case = TRUE)[[1]]
  log2FCColumn <- grep("^(log2?[._ ]?(fold)?[._ ]?(change|ratio|fc)|logfc)$", colnames(seqAnno), value = TRUE, ignore.case = TRUE)[[1]]
  # Sort the dataframe by raw pvalue
  seqAnno <- seqAnno[order(seqAnno[[pValueColumn]]),]
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
      names = seqAnno$gene_name)
  } else if (param$featureLevel == "isoform") {
    seqAnno$feature_id <- scuttle::uniquifyFeatureNames(
      ID = seqAnno$transcript_id,
      names = seqAnno$gene_name)
  }
  # If uniquifying broke some names, fix them
  seqAnno$feature_id[grep("^_", seqAnno$feature_id)] <- gsub("_", "", seqAnno$feature_id[grep("^_", seqAnno$feature_id)])
  # get the list of genes to go in as input options:
  genes <- seqAnno$feature_id[order(seqAnno$pValue)]
  
  
  
  # Prepare the metadata table ----
  dataset <- data.frame(colData(se), check.names = FALSE)
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
  dataset <- dataset %>% mutate_all(~ifelse(. == "", "None", .))
  if (is.null(param$groupingName)) {
    param$groupingName <- factors[[1]]
  }
  # Set the comparison column levels to be in order of test, e.g., diff, undiff
  dataset[[param$groupingName]] <- factor(
    dataset[[param$groupingName]], 
    levels = c(
      param$sampleGroup, 
      param$refGroup, 
      unique(dataset[[param$groupingName]])[!unique(dataset[[param$groupingName]]) %in% c(param$sampleGroup, param$refGroup)])
  )
  # make the sample names an actual column to make it easier to bind to things:
  dataset$names <- rownames(dataset)
  
  
  
  # Prepare the count files ----
  countList <- list()
  countList[["Raw"]] <- assay(se, "counts")
  countList[["Normalised"]] <- assay(se, "xNorm")
  countList[["Normalised + Log2"]] <- log2(assay(se, "xNorm")+param$backgroundExpression)
  if (!is.null(param$refBuild)) {
    countList[["FPKM"]] <- getRpkm(se)
    countList[["TPM"]] <- getTpm(se)
  }
  forVST <- as.matrix(assay(se, "counts"))
  mode(forVST) <- "integer"
  vstCountMatrix <- varianceStabilizingTransformation(forVST)
  countList[["VST"]] <- vstCountMatrix
  # For all the counts, match the matrices to the seqAnno order
  for (i in seq_along(countList)) {
    if(param$featureLevel == "gene") {
      geneOrder <- intersect(seqAnno$gene_id, rownames(countList[[i]]))
      countList[[i]] <- countList[[i]][geneOrder, ]
      rownames(countList[[i]]) <- seqAnno$feature_id
    } else if(param$featureLevel == "isoform") {
      geneOrder <- intersect(seqAnno$transcript_id, rownames(countList[[i]]))
      countList[[i]] <- countList[[i]][geneOrder, ]
      rownames(countList[[i]]) <- seqAnno$feature_id %>% gsub("[:/()-]", ".", .) %>% gsub(" ", ".", .)
    }
  }
  # Add global and per group means to the seqAnno table 
  for (count in names(countList)) {
    cts <- countList[[count]]
    seqAnno[[paste0(count, "_mean")]] <- rowMeans(cts)[seqAnno$feature_id]
    seqAnno[[paste0(count, "_SD")]] <- matrixStats::rowSds(cts)[seqAnno$feature_id]
    isSample <- dataset[[param$groupingName]] == param$sampleGroup
    seqAnno[[paste0(count, "_", param$sampleGroup, "_Mean")]] <- rowMeans(cts[ , isSample, drop=FALSE])[seqAnno$feature_id]
    isRef <- dataset[[param$groupingName]] == param$refGroup
    seqAnno[[paste0(count, "_", param$refGroup, "_Mean")]] <- rowMeans(cts[ , isRef, drop=FALSE])[seqAnno$feature_id]
  }
  
  
  
  # Prepare the contrast info ----
  design <- param$comparison
  designLevels <- design %>% strsplit(split = "--over--") %>% .[[1]]
  
  
  
  # Prepare the GO info ----
  allPathwaysDEG <- c(
    unique(unlist(strsplit(x = c(seqAnno$`GO BP`), "; "))),
    unique(unlist(strsplit(x = c(seqAnno$`GO MF`), "; "))),
    unique(unlist(strsplit(x = c(seqAnno$`GO CC`), "; ")))
  )
  goTerms <- clusterProfiler::go2term(allPathwaysDEG)
  allPathways <- paste(goTerms[,1], goTerms[,2])
  goObjects <- list(runGO = param$runGO)
  if (param$runGO) {
    goObjects$ora <- metadata(se)$enrichResult
    goObjects$gsea <- metadata(se)$enrichResultGSEA
  }
  
  if (param$runGO) {
    return(list(
      "dataType" = "RNASeq",
      "se" = se,
      "param" = param,
      "seqAnno" = seqAnno,
      "design" = param$comparison,
      "designLevels" = designLevels,
      "dataset" = dataset,
      "factors" = factors,
      "factorLevels" = factorLevels,
      "countList" = countList,
      "genes" = genes,
      "allPathways" = allPathways,
      "goObjects" = goObjects,
      "groupingName" = param$groupingName
    )
    )
  }
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
    sa
  }), contrasts)

  # Prepare the metadata
  dataset <- as.data.frame(colData(se))
  names <- colnames(dataset)[grep("^Sample$|^Name$|^sampleName$", colnames(dataset), ignore.case = T)]
  dataset$names <- dataset[[names]]
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

  # Build gene → feature_id lookup (only useful for sub-protein data)
  gene_lookup <- NULL
  if ("gene_name" %in% colnames(seqAnnoList[[1]])) {
    all_sa <- do.call(rbind, lapply(seqAnnoList, function(sa) {
      sa[, c("feature_id", "gene_name"), drop = FALSE]
    }))
    all_sa <- unique(all_sa)
    if (!identical(all_sa$feature_id, all_sa$gene_name)) {
      gene_lookup <- split(all_sa$feature_id, all_sa$gene_name)
    }
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
    "gene_lookup" = gene_lookup,
    "factorLevels" = factorLevels,
    "groupingName" = factors[[1]]
    )
  )
}

#' zscore
#' @export
zscore <- function(x) {
  # Remove NAs for mean and sd calculation
  if (all(is.na(x))) {
    # If all values are NA, return all NAs
    return(rep(NA, length(x)))
  }
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    # If sd is NA (all NAs) or zero (no variance), return zeros (or NAs)
    return(rep(0, length(x)))
  }
  (x - m) / s
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
