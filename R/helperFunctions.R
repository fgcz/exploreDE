#' convert_genomics_se
#' @export
convert_genomics_se <- function(se, dataDir) {
  param <- metadata(se)$param
  seqAnno <- data.frame(
    rowData(se),
    row.names = rownames(se),
    check.names = FALSE)
  
  # Uniquify the gene symbols 
  seqAnno <- seqAnno[order(seqAnno$pValue),]
  seqAnno$description <- gsub(" \\[..*", "", seqAnno$description)
  
  # change any NA FDRs to 1:
  if (param$featureLevel == "gene") {
    seqAnno$gene_name <- scuttle::uniquifyFeatureNames(
      ID = seqAnno$gene_id,
      names = seqAnno$gene_name)
  } else if (param$featureLevel == "isoform") {
    seqAnno$gene_name <- scuttle::uniquifyFeatureNames(
      ID = seqAnno$transcript_id, 
      names = seqAnno$gene_name)
  }
  seqAnno$gene_name[grep("^_", seqAnno$gene_name)] <- gsub("_", "", seqAnno$gene_name[grep("^_", seqAnno$gene_name)])
  dataset <- data.frame(colData(se), check.names = FALSE)
  
  # remove things like [] & whitespace from colnames to make subsetting easier:
  # I appear to have got confused somewhere with my need for some of these factor variables... will fix 
  factorNs <- grep("\\[Factor]", colnames(dataset))
  factors <- colnames(dataset)[grep("\\[Factor]", colnames(dataset))]
  factorNames <- colnames(dataset)[grep("\\[Factor]", colnames(dataset))]
  colnames(dataset) <- gsub(" \\[.*", "", colnames(dataset)) %>% gsub(" ", "_", .)
  factorNames <- gsub(" \\[.*", "", factorNames) %>% gsub(" ", "_", .)
  factors <- gsub(" \\[.*", "", factors) %>% gsub(" ", "_", .)
  factorLevels <- NULL
  for (i in seq_along(dataset[factorNs])){
    factorLevels[[i]] <- paste0(
      colnames(dataset[factorNs[[i]]]),
      ": ",
      levels(as.factor(dataset[, factorNs[i]])))
  }
  factorLevels <- unlist(factorLevels)
  names(factorLevels) <- factorLevels %>% gsub("\\ |\\:", "_", .)
  
  # The app doesn't like empty strings in the factor columns, so for now, add "None" to those entries 
  dataset <- dataset %>% mutate_all(~ifelse(. == "", "None", .))
  if (is.null(param$groupingName)) {
    param$groupingName <- "Condition"
  }
  dataset[[param$groupingName]] <- factor(
    dataset[[param$groupingName]], 
    levels = c(
      param$sampleGroup, 
      param$refGroup, 
      unique(dataset[[param$groupingName]])[!unique(dataset[[param$groupingName]]) %in% c(param$sampleGroup, param$refGroup)]))
  
  # Make a list of counts:
  countList <- list()
  countList[["Raw"]] <- assay(se, "counts")
  countList[["Normalised"]] <- assay(se, "xNorm")
  countList[["Normalised + Log2"]] <- log2(assay(se, "xNorm")+param$backgroundExpression)
  countList[["FPKM"]] <- getRpkm(se)
  countList[["TPM"]] <- getTpm(se)
  for (i in seq_along(countList)) {
    if(param$featureLevel == "gene") {
      countList[[i]] <- countList[[i]][which(rownames(countList[[i]]) %in% seqAnno$gene_id), ]
      rownames(countList[[i]]) <- seqAnno$gene_name[match(rownames(countList[[i]]), seqAnno$gene_id)]
    } else if(param$featureLevel == "isoform") {
      countList[[i]] <- countList[[i]][rownames(countList[[i]]) %in% seqAnno$transcript_id, ]
      rownames(countList[[i]]) <- seqAnno$gene_name[match(rownames(countList[[i]]), seqAnno$transcript_id)] %>% gsub("[:/()-]", ".", .) %>% gsub(" ", ".", .)
    }
  }
  forVST <- as.matrix(assay(se, "counts"))
  mode(forVST) <- "integer"
  vstCountMatrix <- varianceStabilizingTransformation(forVST)
  if(param$featureLevel == "gene") {
    vstCountMatrix <- vstCountMatrix[intersect(rownames(vstCountMatrix), seqAnno$gene_id), ]
    vstCountMatrix <- vstCountMatrix[match(seqAnno$gene_id, rownames(vstCountMatrix)), ]
    rownames(vstCountMatrix) <- seqAnno$gene_name
  } else if (param$featureLevel == "isoform") {
    vstCountMatrix <- vstCountMatrix[intersect(rownames(vstCountMatrix), seqAnno$transcript_id), ]
    vstCountMatrix <- vstCountMatrix[match(seqAnno$transcript_id, rownames(vstCountMatrix)), ]
    rownames(vstCountMatrix) <- seqAnno$gene_name %>% gsub("[:/()-]", ".", .) %>% gsub(" ", ".", .)
  }
  countList[["VST"]] <- vstCountMatrix
  for (count in names(countList)) {
    cts <- countList[[count]]
    seqAnno[[paste0(count, "_mean")]] <- rowMeans(cts)[seqAnno$gene_name]
    seqAnno[[paste0(count, "_SD")]] <- matrixStats::rowSds(cts)[seqAnno$gene_name]
    isSample <- dataset[[param$groupingName]] == param$sampleGroup
    seqAnno[[paste0(count, "_", param$sampleGroup, "_Mean")]] <- rowMeans(cts[ , isSample, drop=FALSE])[seqAnno$gene_name]
    isRef <- dataset[[param$groupingName]] == param$refGroup
    seqAnno[[paste0(count, "_", param$refGroup, "_Mean")]] <- rowMeans(cts[ , isRef, drop=FALSE])[seqAnno$gene_name]
  }
  
  # make the sample names an actual column to make it easier to bind to things:
  dataset$names <- rownames(dataset)
  
  # get the list of genes to go in as input options:
  genes <- seqAnno$gene_name[order(seqAnno$pValue)]
  
  # get contrast 
  design <- param$comparison
  designLevels <- design %>% strsplit(split = "--over--") %>% .[[1]]
  
  allPathwaysDEG <- c(
    unique(unlist(strsplit(x = c(seqAnno$`GO BP`), "; "))),
    unique(unlist(strsplit(x = c(seqAnno$`GO MF`), "; "))),
    unique(unlist(strsplit(x = c(seqAnno$`GO CC`), "; ")))
  )
  goTerms <- clusterProfiler::go2term(allPathwaysDEG)
  allPathways <- paste(goTerms[,1], goTerms[,2])
  
  if (param$runGO) {
    # Get the pathway objects 
    ora <- metadata(se)$enrichResult
    gsea <- metadata(se)$enrichResultGSEA
    oraHTML <- file.path(dataDir, list.files(dataDir)[grep("ORA_results.xlsx", list.files(dataDir))])
    gseaHTML <- file.path(dataDir, list.files(dataDir)[grep("GSEA_results.xlsx", list.files(dataDir))])
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
      "factorNames" = factorNames,
      "degHTML" = file.path(dataDir, list.files(dataDir)[grep("xlsx", list.files(dataDir))]) %>% .[!grepl("ORA|GSEA", .)],
      "countList" = countList,
      "genes" = genes,
      "factorLevels" = factorLevels,
      "allPathways" = allPathways,
      "ora" = ora,
      "oraHTML" = oraHTML,
      "gsea" = gsea,
      "gseaHTML" = gseaHTML
    )
    )
  } else {
    return(list(
      "dataType" = "RNASeq",
      "se" = se,
      "param" = param,
      "seqAnno" = seqAnno,
      "design" = design,
      "designLevels" = designLevels,
      "dataset" = dataset,
      "factors" = factors,
      "factorNames" = factorNames,
      "degHTML" = file.path(dataDir, list.files(dataDir)[grep("xlsx", list.files(dataDir))]) %>% .[!grepl("ORA|GSEA", .)],
      "countList" = countList,
      "genes" = genes,
      "factorLevels" = factorLevels,
      "allPathways" = allPathways
    )
    )
  }
}

#' convert_proteomics_se
#' @export
convert_proteomics_se <- function(se) {
  # Get the contrasts available 
  contrasts <- names(rowData(se))[grep("^constrast_", names(rowData(se)))] %>% gsub("constrast_", "", .)
  
  # Make a list of seqAnno tables for each contrast
  seqAnnoList <- setNames(lapply(contrasts, function(con) {
    sa <- rowData(se)[[paste0("constrast_", con)]]
    sa <- sa[order(sa$p.value),]
    sa <- sa %>% dplyr::select(any_of(c("protein_Id", "site", "fasta.id", "IDcolumn", "diff", "p.value", "FDR", "description", "nrPeptides", "modelName")))
    if (any(grepl("site", colnames(sa)))) {
      sa$site <- gsub("\\~", ".", sa$site)
      sa <- sa %>% dplyr::select(-protein_Id)
      colnames(sa)[grep("site", colnames(sa))] <- "protein_Id"
    }
    colnames(sa)[grep(paste(c("protein_Id", "diff", "p.value", "FDR"), collapse = "|"), colnames(sa))] <- c("gene_name", "log2Ratio", "pValue", "fdr")
    sa
  }), contrasts)
  
  # Prepare the metadata
  dataset <- as.data.frame(colData(se))
  names <- colnames(dataset)[grep("^Sample$|^Name$|^sampleName$", colnames(dataset), ignore.case = T)]
  dataset$names <- dataset[[names]]
  rownames(dataset) <- dataset[[names]]
  
  # Get the main columns from the metadata that are plottable 
  factors <- colnames(dataset)[grep("Sample|Name|names|channel|raw\\.file", colnames(dataset), invert = T, ignore.case = T)]
  factorNames <- colnames(dataset)[grep("Sample|Name|names|channel|raw\\.file", colnames(dataset), invert = T, ignore.case = T)]
  
  # Sort numeric factor levels
  for (x in factors) {
    dataset[[x]] <- factor(dataset[[x]], levels = unique(gtools::mixedsort(dataset[[x]])))
  }
  
  # Currently only transformed counts available
  countList <- setNames(lapply(names(assays(se)), function(ass) {
    x <- assay(se, ass)
    if (any(grepl("site", colnames(rowData(se)[[1]])))) {
      rownames(x) <- gsub("\\~", ".", rownames(x))
      rownames(x) <- gsub("\\.lfq\\.light", "", rownames(x))
      rownames(x) <- gsub(".*lfq\\.", "", rownames(x))
    } else {
      rownames(x) <- gsub("\\~.*", "", rownames(x))
    }
    x
  }), names(assays(se)))
  
  # Get the protein IDs
  genes <- unique(unlist(lapply(seq_along(seqAnnoList), function(x) {seqAnnoList[[x]][["gene_name"]]})))
  
  # Get a vector of all the levels of each plottable column
  factorLevels <- NULL
  for (i in seq_along(dataset[factors])){
    factorLevels[[i]] <- paste0(
      colnames(dataset[factors[[i]]]),
      ": ",
      levels(as.factor(dataset[, factors[i]])))
  }
  factorLevels <- unlist(factorLevels)
  names(factorLevels) <- factorLevels %>% gsub("\\ |\\:", "_", .)
  
  # Return the final list
  return(list(
    "dataType" = "proteomics",
    "se" = se,
    "seqAnnoList" = seqAnnoList,
    "contrasts" = contrasts,
    "dataset" = dataset,
    "factors" = factors,
    "factorNames" = factorNames,
    "countList" = countList,
    "genes" = genes,
    "factorLevels" = factorLevels
  )
  )
}
