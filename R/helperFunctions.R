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
  names <- colnames(dataset)[grep("Sample|Name", colnames(dataset))]
  dataset$names <- dataset[[names]]
  rownames(dataset) <- dataset[[names]]
  
  # Get the main columns from the metadata that are plottable 
  factors <- colnames(dataset)[grep("Sample|Name|names|channel|raw.file", colnames(dataset), invert = T)]
  factorNames <- colnames(dataset)[grep("Sample|Name|names|channel|raw.file", colnames(dataset), invert = T)]
  
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