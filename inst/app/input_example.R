# Helper functions

#' standardise_columns
#' Fuzzy-match incoming column names to canonical names.
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
build_feature_id <- function(df) {
  pid <- df$protein_Id
  if (!anyDuplicated(pid)) return(pid)

  # Try known sub-feature columns (canonical names) in priority order
  for (col in intersect(c("site", "peptide_Id"), colnames(df))) {
    combined <- paste(pid, gsub("~", "_", df[[col]]), sep = "_")
    if (!anyDuplicated(combined)) return(combined)
  }

  # Fallback: cleaned raw rownames
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
  # dataset <- data.frame(colData(se), check.names = FALSE)
  dataset <- read_tsv(file.path(dataDir, "../input_dataset.tsv"), show_col_types = F) %>% data.frame(check.names = F)
  dataset[["Name"]] <- factor(dataset[["Name"]], levels = dataset[["Name"]])
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
  # forVST <- metadata(se)$nativeResult$dds
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
  
  return(
    list(
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


#' convert_proteomics_se
#' @export
convert_proteomics_se <- function(se) {
  # Get the contrasts available 
  contrasts <- names(rowData(se))[grep("^constrast_", names(rowData(se)))] %>% gsub("constrast_", "", .)
  
  # Make a list of seqAnno tables for each contrast
  seqAnnoList <- setNames(lapply(contrasts, function(con) {
    sa <- as.data.frame(rowData(se)[[paste0("constrast_", con)]])
    sa <- standardise_columns(sa)
    sa <- sa[order(sa$pValue), ]

    # Build feature_id from canonical column names
    sa$feature_id <- build_feature_id(sa)

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
  raw_rn <- rownames(as.data.frame(rowData(se)[[paste0("constrast_", contrasts[1])]]))
  fid_lookup <- setNames(seqAnnoList[[1]]$feature_id, raw_rn)

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
    "groupingName" = factors[[1]]
  )
  )
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



is_url <- function(dataDir) {
  return(grepl("^https?://", dataDir)) 
}
# dataDir <- "/srv/gstore/projects/p3009/o5638_DESeq2_diff--over--undiff_2024-11-20--12-29-40/diff--over--undiff/"
# dataDir <- "/srv/gstore/projects/p28466/eq2_dnaj1dnaj2dnaj3pepe--over--WT_2025-07-07--12-43-22/dnaj1dnaj2dnaj3pepe--over--WT"
# dataDir <- "/srv/gstore/projects/p35043/o39540_DESeq2_AGD_Tumor--over--Control_Tumor_2025-10-01--12-34-19/AGD_Tumor--over--Control_Tumor"
# dataDir <- "/srv/gstore/projects/p26523/o36354_DESeq2_Variant_CD--over--WT_CD_2025-07-23--16-53-53/Variant_CD--over--WT_CD"
# dataDir <- "/srv/gstore/projects/p33393/o37634_EdgeR_P3NS-shKCNQ2--over--P3NS-shSCR_2025-05-16--13-06-04/P3NS-shKCNQ2--over--P3NS-shSCR"
# dataDir <- "/srv/gstore/projects/p33393/o37976_DESeq2_BLZ--over--Vehicle_2025-05-16--12-57-24/BLZ--over--Vehicle"
# dataDir <- "https://fgcz-ms.uzh.ch/public/pStore/../jStore/p2699/bfabric/Proteomics/SummarizedExperiment/2023/2023-11/2023-11-23/workunit_o31423_p2699_GeneSymbolRollUp_unpaired/SummarizedExperiment.rds"
# dataDir <- "/srv/gstore/projects/p26523/o36354_DESeq2_Variant_Inflamed--over--Variant_Non_Inflamed_2025-07-23--16-57-19/Variant_Inflamed--over--Variant_Non_Inflamed"
# dataDir <- "/srv/gstore/projects/p38300/o38603_DESeq2_Sorbitol_132mg--over--NaCl_corrected_2025-06-30--12-37-44/Sorbitol_132mg--over--NaCl"
# dataDir <- "/srv/gstore/projects/p2753/EdgeR_Shmt1_KO--over--Shmt1_Ung_dKO_2024-08-09--05-44-48/Shmt1_KO--over--Shmt1_Ung_dKO"
# dataDir <- "/srv/gstore/projects/p39127/o39272_DESeq2_J1_withDox_withIFN--over--J1_withDox_withoutIFN_2025-09-05--12-17-56/J1_withDox_withIFN--over--J1_withDox_withoutIFN"
# dataDir <- "https://fgcz-ms.uzh.ch/public/pStore/p38633/bfabric/Proteomics/exploreDE/2025/2025-07/2025-07-01//workunit_328724//2913198.rds"
# dataDir <- "https://fgcz-ms.uzh.ch/public/pStore/p39819/bfabric/Proteomics/exploreDE/2025/2025-10/2025-10-22/workunit_335484/2997993.rds"
# dataDir <- "https://fgcz-ms.uzh.ch/public/pStore/p39819/bfabric/Proteomics/exploreDE/2025/2025-10/2025-10-22/workunit_335485/2997996.rds"
# dataDir <- "https://fgcz-ms.uzh.ch/public/pStore/p39763/bfabric/Proteomics/exploreDE/2025/2025-11/2025-11-03/workunit_335997/3012890.rds"
# dataDir <- "/srv/gstore/projects/p39374/o39634_DESeq2_CEPC_COMP--over--CEPC_unstimulated_2025-11-13--09-14-39/CEPC_COMP--over--CEPC_unstimulated"
# dataDir <- "/srv/gstore/projects/p23793/o32655_CountQC_2023-10-05--18-33-15/Count_QC/counts-sswoowrqapiv-EzResult.RData"
# dataDir <- "/srv/gstore/projects/p23793/o32655_CountQC_2023-10-05--18-33-15/Count_QC/counts-sswoowrqapiv-EzResult.RData" # CountQC
# dataDir <- "https://fgcz-ms.uzh.ch/public/pStore/p39763/bfabric/Proteomics/exploreDE/2025/2025-11/2025-11-03/workunit_335997/3012890.rds" # Proteomics
# dataDir <- "/srv/gstore/projects/p39374/o39634_DESeq2_IVD_Disc_Matrix--over--IVD_unstimulated_2025-11-13--09-20-35/IVD_Disc_Matrix--over--IVD_unstimulated" # twoFactor
# dataDir <- "/srv/gstore/projects/p38738/o40056_EdgeR_HDM_F2_STI--over--HDM_F2_2025-11-30--08-06-19/HDM_F2_STI--over--HDM_F2"
# dataDir <- "/srv/gstore/projects/p38738/o40056_EdgeR_HDM_F2--over--HDM_FR2_2025-12-03--11-48-47/HDM_F2--over--HDM_FR2"
dataDir <- "/srv/gstore/projects/p29284/o36981_EdgeR_HR2_B_cells--over--HR1_B_cells_2025-06-27--09-16-12/HR2_B_cells--over--HR1_B_cells"
dataUrl <- "p29284/o36981_EdgeR_HR2_B_cells--over--HR1_B_cells_2025-06-27--09-16-12/HR2_B_cells--over--HR1_B_cells"
dataDir <- "/srv/gstore/projects/p29779/o38737_DESeq2_AP_C_shp65_w_HP--over--AP_C_shp65_wo_HP_2025-07-29--10-21-00/AP_C_shp65_w_HP--over--AP_C_shp65_wo_HP"
username <- Sys.getenv("SHINYPROXY_USERNAME")
if (is_url(dataDir)) {
  se <- readRDS(url(dataDir))  # If it's a URL
} else {
  se <- readRDS(file.path(dataDir, "deResult.rds"))
}
inputDataReactive <- reactive({
  # Load proteomics data ----
  if (grepl("rds|zip", dataDir)) {
    return(convert_proteomics_se(se))
  }
  
  # Load RNA Seq data ----
  if (file.exists(file.path(dataDir, "deResult.rds"))) {
    return(convert_genomics_se(se, dataDir))
  }
})
inputDataReactive()$dataType
input <- list()
names(inputDataReactive())

colourPaletteList <- list(
  "House colours" = c(
    "indianred3", "steelblue4", "chartreuse4", "grey30", 
    "goldenrod3", "firebrick4", "royalblue4", "mediumorchid3",
    "turquoise4", "darkolivegreen", "thistle4", "darkorange3", 
    "hotpink2", "burlywood3", "cadetblue4", "chocolate4", "firebrick"
  ),
  "Paired" = brewer.pal(12, "Paired"),
  "Set1" = brewer.pal(9, "Set1"),
  "Set2" = brewer.pal(8, "Set2"),
  "Set3" = brewer.pal(12, "Set3"),
  "Dark2" = brewer.pal(8, "Dark2")
)

input$colourPalette <- c("Dark2", "Paired")
i <- 1
j <- 13
for (i in seq_along(inputDataReactive()$factorLevels)) {
  for (j in seq_along(inputDataReactive()$factorLevels[[i]])) {
    input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[i], inputDataReactive()$factorLevels[[i]][j])]] <- rep(as.character(unlist(colourPaletteList[input$colourPalette])), times = 5)[j]
  }
}
input$selectedTable_ORA_rows_selected <- 1:5
input$textSizeORA <- 12
input$degTable_rows_selected <- c(12,68)
input$featureCalcPType <- "Raw"
input$featureCalcPValue <- 0.05
input$featureCalcLFC <- 0

input$pcaFactor1 <- inputDataReactive()$factors[1]
input$pcaGroups <- inputDataReactive()$factorLevels[[1]]
# input$pcaGroups <- levels(inputDataReactive()$dataset$G_)
input$pcaCounts <- "Normalised + Log2"
input$pcaAddEllipses <- TRUE
input$pcaEllipsesAlpha <- 0.2
# input$pcaCounts <- "transformedData"
input$pcaBatch <- "None"
input$pcaTopN <- 2000
input$pcaCentre = TRUE 
input$pcaScale = FALSE
input$textSizePCA <- 12
input$pcaX = "PC1"
input$pcaY = "PC2"
input$pcaFactor2 <- "None"
input$pcaAxesProp <- TRUE
input$pcaShowNames <- TRUE
input$showLinesPCA <- TRUE
input$showAxesPCA <- TRUE
input$boldPCA <- TRUE
input$pcaLog2 <- TRUE
input$figWidthPCA <- 600
input$figHeightPCA <- 600
input$textSizePCA <- 12
input$dotSizePCA <- 3
input$alphaPCA <- 0.9
input$dotBorderPCA <- 0.2
input$showAxesPCA <- TRUE
input$boldPCA <- TRUE
input$showLinesPCA <- TRUE
input$pcaLabelMaxOverlap <- 10
input$geneLabelNudgePCAX <- 0
input$geneLabelNudgePCAY <- 0
input$geneLabelSizePCA <- 12
input$downloadFormatPCA <- "PDF"
input$dpiPCA <- 600
input$filnamePCA <- "PCA"

# Ordination method options (PCA, PCoA, NMDS, GLM-PCA)
input$ordinationMethod <- "pca"
input$ordinationPlotType <- "samples"
input$pcoaDistance <- "euclidean"
input$glmpcaFamily <- "poi"
input$glmpcaMinibatch <- 100
input$biplotNGenes <- 10
input$biplotArrowScale <- 1
input$biplotArrowColour <- "#E41A1C"
input$nrPeptidePCA <- 2

input$boxplotFactor1 <- inputDataReactive()$factors[1]
input$fb_keep_boxplot <- inputDataReactive()$genes[1:5]
# input$keepBucketBoxplot <- grep("Hrh1", inputDataReactive()$genes, value = T)
input$boxplotCounts <- names(inputDataReactive()$countList)[2]
input$boxplotBatch <- "None"
input$boxplotFactor2 <- "None"
input$boxKeepBucket <- unique(levels(as.factor(inputDataReactive()$dataset[[input$boxplotFactor1]])))
input$boxplotShowP <- FALSE
input$boxplotShowBox <- TRUE
input$boxplotShowViolin <- FALSE
input$boxplotShowMeanBar <- TRUE
input$boxplotShowPoint <- TRUE
input$boxplotPointDodge <- 2
input$boxplotBoxAlpha <- 0.2
input$boxplotPoints <- "Fancy"
input$textSizeBoxplot <- 12
input$boxplotGrey <- FALSE
input$boxplotThemeChoice <- "Prism"
input$boxplotNCol <- 3
input$boxplotVertLines <- FALSE 
input$figHeightBoxplot <- 800
input$figWidthBoxplot <- 800
input$showDotsBarplot <- TRUE
input$boxplotPointSize <- 5
input$boxplotPointAlpha <- 0.9
input$boxplotMeanLine <- 0.5
input$boxplotShowPHeight <- 1
input$boxplotShowPLabelSize <- 3
input$boxplotShowPBracketSize <- 3
input$boxplotShowPTipSizeA <- 0
input$boxplotShowPTipSizeB <- 0
input$boxplotShowPDodge <- 1
input$boxplotPointBorder <- 0.5
input$boxplotCountsLog <- FALSE
input$boxplotZScore <- FALSE

input$fb_select_boxplot <- inputDataReactive()$genes[1:5]
input$fb_text_boxplot <- inputDataReactive()$genes[6:10]
input$boxplotFactor1 <- inputDataReactive()$factors[1]
input$boxplotFactor2 <- "None"
input$boxplotCounts <- names(inputDataReactive()$countList)[2]
input$boxplotBatch <- "None"
input$boxplotShowP <- TRUE
input$boxplotShowPoint <- TRUE
input$boxplotShowMeanBar <- TRUE
input$boxplotShowBox <- FALSE
input$boxplotShowViolin <- TRUE
input$boxplotPointSize <- 3
input$boxplotPointAlpha <- 0.9
input$boxplotBoxAlpha <- 0.2
input$boxplotMeanLine <- 0.5
# input$boxplotPointJitter
# input$boxplotMeanLineTop
# input$boxplotPoints
input$boxplotThemeChoice <- "Prism"
input$boxplotVertLines <- FALSE
input$boxplotGrey <- FALSE
input$boxplotNCol <- 3
input$textSizeBoxplot <- 12
input$figWidthBoxplot <- 600
input$figHeightBoxplot <- 900
input$showDotsBarplot <- TRUE


input$heatmapColumnSelect <- inputDataReactive()$factors[1:2]
if (inputDataReactive()$dataType == "proteomics") {
  input$contrastSelected <- inputDataReactive()$contrasts[1]
  input$downloadCount <- "transformed"
}

# input$heatmapScale <- "Continuous"
input$pTypeHeatmap <- "Raw"
# input$hmKeepBucket <- unique(levels(as.factor(inputDataReactive()$dataset[,inputDataReactive()$factors[3]])))
input$heatmapCounts <- switch(inputDataReactive()$dataType, "RNASeq" = "TPM", "proteomics" = "transformedData")
input$xCoord <- 10
input$yCoord <- 10

input$goTileInput <- "GO:0042254 ribosome biogenesis"
input$pTypeGoTile <- "FDR"
input$tileLimitColour <- 4
input$degOnlyGoTile <- FALSE
input$clusterRowsGoTile <- TRUE
input$colnamesGoTile <- TRUE
input$geneNamesGoTile <- TRUE
input$GoTileColourRed <- "#801717"
input$GoTileColourWhite <- "#FFFFFF"
input$GoTileColourBlue <- "#113D69"
input$textSizeGoTile <- 12
input$figWidthGoTile <- 400
input$figHeightGoTile <- 800

input$heatmapDPI <- 600
input$heatmapFactors <- inputDataReactive()$factors
input$heatmapDownloadFormat <- "PDF"
input$heatmapFilename <- "both_regulated_heatmap"
input$heatmapLimitD <- 2
input$heatmapColourRed <- "#801717"
input$heatmapColourWhite <- "#FFFFFF"
input$heatmapColourBlue <- "#113D69"
input$pThresholdHeatmap <- 0.05
input$lfcHeatmap <- 0
input$heatmapGeneNumber <- 50
input$clusterColsHeatmap <- TRUE
input$clusterRowsHeatmap <- TRUE
input$colnamesHeatmap <- TRUE
input$geneNamesHeatmap <- TRUE
input$figHeightHeatmap <- 800
input$figWidthHeatmap <- 1000
input$textSizeHeatmap <- 12
input$heatmapBatch <- "None"
input$heatmapLimitCPalette <- "PRGn"
input$heatmapLimitCPaletteRev <- FALSE
input$heatmapLimitCHigh <- 1e5
input$heatmapLimitCMid <- 1e4
input$heatmapLimitCLow <- 1e3
input$heatmapFactor2 <- ""
input$showClusterColDend <- TRUE
input$showClusterRowDend <- TRUE
input$heatmapGoInput <- "GO:0042254 ribosome biogenesis"
# input$heatmapGoInput <- inputDataReactive()$seqAnno$`GO BP`[[2]] %>% strsplit("; ") %>% .[[1]] %>% .[1]
input$heatmapGoExpr <- 0.05
input$heatmapLog2 <- TRUE
input$heatmapZScore <- "Centred"
input$heatmapFeatureDirection <- "Both"
input$heatmapAtLow <- -4
input$heatmapAtMid <- 0
input$heatmapAtHigh <- 4

input$pTypeVolcano <- "FDR"
input$lfcVolcano <- 0.25
input$pThresholdVolcano <- 0.05
input$textSizeVolcano <- 12
input$alphaVolcano <- 0.8
input$fb_select_volcano <- inputDataReactive()$genes[1:5]
input$fb_text_volcano <- inputDataReactive()$genes[8:9]
input$volcanoShowGenes <- TRUE
input$fb_keep_volcano <- inputDataReactive()$genes[1:5]
input$dotSizeVolcano <- 3
input$yLimVolcano <- 25
input$xLimVolcano <- 30
input$volcanoLabelAllUp <- TRUE
input$volcanoLabelAllDown <- TRUE
input$volcanoAnnotationHighlightColour <- FALSE
input$volcanoLabelMaxOverlap <- 50
input$volcanoPointBorder <- 0.5
input$downloadFormatVolcano <- "PDF"
input$dpiVolcano <- 600
input$filnameVolcano <- "Volcano"

input$oraType <- "BP"
input$oraDirection <- "upGenes"
input$gseaType <- "BP"
# input$gseaType <- "C2"
input$selectedTable_GSEA_rows_selected <- c(1,2,3,4)
input$textSizeGSEA <- 12

input$nrPeptideVolcano <- 1
input$showImputedVolcano <- TRUE
input$highlightImputedVolcano <- TRUE


volcanoColoursDefault <- c(
  "Highlight" = "goldenrod3",
  "NotSignificant" = "grey30",
  "SignificantUp" = "pink",
  "SignificantDown" = "lightblue",
  "FoldChangeUp" = "pink3",
  "FoldChangeDown" = "lightblue4",
  "FoldChangeSignificantUp" = "firebrick4",
  "FoldChangeSignificantDown" = "dodgerblue4"
)

if (inputDataReactive()$dataType == "proteomics") {
  req(!is.null(inputDataReactive()$seqAnnoList))
  if("modelName" %in% colnames(inputDataReactive()$seqAnnoList[[1]])) {
    volcanoColoursDefault <- c(volcanoColoursDefault, "Imputed" = "darkturquoise")
  }
}
for (x in names(volcanoColoursDefault)) {
  input[[paste0("volcanoColour", x) ]] <- as.character(volcanoColoursDefault[x])
}

input$pTypeKEGG <- "FDR"
input$pThresholdKEGG <- 0.05
input$lfcKEGG <- 0.5
input$keggLimitColour <- 1


input$boxplotShowConditions <- TRUE
input$boxplotConditionAngle <- 45
input$boxplotConditionFormat <- TRUE
input$boxplotMeanBarFront <- TRUE

input$scaleLimORA <- 4
input$nodeSizeORA <- 10
input$oraMaxOver <- 10
input$oraDodgeY <- 1
input$oraDodgeX <- 1
input$nodeBorderORA <- 1

input$showGeneLabelsGSEA <- TRUE
input$gseaLabelAlpha <- 0.2
input$scaleLimGSEA <- 4
input$nodeSizeGSEA <- 10
input$gseaMaxOver <- 10
input$gseaDodgeY <- 1
input$gseaDodgeX <- 1
input$nodeBorderGSEA <- 1
input$ridgePlotColourByGSEA <- "qvalue"

input$showGeneLabelsORA <- TRUE

input$contrastSelected

input$correlationGene1 <- inputDataReactive()$genes[1]
input$correlationGene2 <- inputDataReactive()$genes[2]
input$correlationColourBy <- inputDataReactive()$factors[1]
input$correlationShapeBy <- "None"
input$correlationCounts <- switch(inputDataReactive()$dataType, "RNASeq" = "Normalised + Log2", "proteomics" = "transformedData")
input$correlationBatch <- "None"
input$correlationCorPosX <- "left"
input$correlationCorPosY <- "top"
input$correlationPointSize <- 2
input$correlationPointStroke <- 0.8
input$correlationPointAlpha <- 1
input$textSizeCorrelation <- 12
input$correlationMethod <- "Pearson"
input$correlationCountsLog <- TRUE
input$correlationShowNames <- TRUE
input$correlationMaxOverlaps <- 10
input$correlationLabelPosX <- 0
input$correlationLabelPosY <- 0
input$correlationKeepBucket <- levels(inputDataReactive()$dataset[[inputDataReactive()$factors[1]]])

input$fb_keep_heatmap <- inputDataReactive()$genes[c(1, 4, 5, 6, 9)]
input$boxplotBatch <- "None"
dataSummary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sqrt(sd(x[[col]]) / length(x[[col]])))
  }
  data_sum <- ddply(
    data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
if (inputDataReactive()$dataType == "RNASeq") {
  seqAnnoReactive <- reactiveValues("sa" = inputDataReactive()$seqAnno, "saF" = NULL)
  design <- inputDataReactive()$design
  designLevels <- inputDataReactive()$designLevels
} else if (inputDataReactive()$dataType == "proteomics") {
  seqAnnoReactive <- reactiveValues("sa" = NULL, "saF" = NULL)
  seqAnnoReactive$sa = inputDataReactive()$seqAnnoList[[input$contrastSelected]]
  seqAnnoReactive$sa <- seqAnnoReactive$sa %>% dplyr::select("feature_id", "log2Ratio", "pValue", "fdr")
  seqAnnoReactive$saF <- seqAnnoReactive$sa
  design <- input$contrastSelected
}
input$barplotSDorSE <- "SE"
input$boxplotShowPJust <- 1
input$boxplotYScaleFree <- TRUE
input$boxplotPlotBorder <- TRUE
input$borderVolcano <- TRUE
input$oraLabelAlpha <- 0.9
input$geneLabelSizeVolcano <- 2
input$showAxesVolcano <- TRUE
input$showLinesVolcano <- TRUE
input$pcaBatch <- "None"
genesToPlot <- inputDataReactive()$seqAnno$feature_id[1:5]
input$cellBorderSizeHeatmap1 <- 5
input$cellBorderSizeHeatmap2 <- 5
input$cellBorderColourHeatmap1 <- "black"
input$cellBorderColourHeatmap2 <- "black"

i <- 1
for (i in seq_along(input$heatmapColumnSelect)) {
  input[[paste0("hmKeepBucket", i)]] <- unique(levels(as.factor(inputDataReactive()$dataset[,input$heatmapColumnSelect[i]])))
}
input$pcaCounts <- names(inputDataReactive()$countList)[[2]]

featureBucket$features <- inputDataReactive()$genes[c(1, 4, 5, 6, 9)]
featureBucket$version <- 1L
