# Helper functions
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
  # Uniquify the gene symbols 
  if (param$featureLevel == "gene") {
    seqAnno$gene_name <- scuttle::uniquifyFeatureNames(
      ID = seqAnno$gene_id,
      names = seqAnno$gene_name)
  } else if (param$featureLevel == "isoform") {
    seqAnno$gene_name <- scuttle::uniquifyFeatureNames(
      ID = seqAnno$transcript_id, 
      names = seqAnno$gene_name)
  }
  # If uniquifying broke some names, fix them
  seqAnno$gene_name[grep("^_", seqAnno$gene_name)] <- gsub("_", "", seqAnno$gene_name[grep("^_", seqAnno$gene_name)])
  # get the list of genes to go in as input options:
  genes <- seqAnno$gene_name[order(seqAnno$pValue)]
  
  
  
  # Prepare the metadata table ----
  if (file.exists(file.path(dataDir, "../input_dataset.tsv"))) {
    dataset <- read_tsv(file.path(dataDir, "../input_dataset.tsv"), show_col_types = F) %>% data.frame(check.names = F)
  } else {
    dataset <- data.frame(colData(se), check.names = FALSE)
  }
  dataset <- dataset[which(dataset$Name %in% param$samples), ]
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
  vstCountMatrix <- varianceStabilizingTransformation(forVST)
  countList[["VST"]] <- vstCountMatrix
  # For all the counts, match the matrices to the seqAnno order
  for (i in seq_along(countList)) {
    if(param$featureLevel == "gene") {
      geneOrder <- intersect(seqAnno$gene_id, rownames(countList[[i]]))
      countList[[i]] <- countList[[i]][geneOrder, ]
      rownames(countList[[i]]) <- seqAnno$gene_name
    } else if(param$featureLevel == "isoform") {
      geneOrder <- intersect(seqAnno$transcript_id, rownames(countList[[i]]))
      countList[[i]] <- countList[[i]][geneOrder, ]
      rownames(countList[[i]]) <- seqAnno$gene_name %>% gsub("[:/()-]", ".", .) %>% gsub(" ", ".", .)
    }
  }
  # Add global and per group means to the seqAnno table 
  for (count in names(countList)) {
    cts <- countList[[count]]
    seqAnno[[paste0(count, "_mean")]] <- rowMeans(cts)[seqAnno$gene_name]
    seqAnno[[paste0(count, "_SD")]] <- matrixStats::rowSds(cts)[seqAnno$gene_name]
    isSample <- dataset[[param$groupingName]] == param$sampleGroup
    seqAnno[[paste0(count, "_", param$sampleGroup, "_Mean")]] <- rowMeans(cts[ , isSample, drop=FALSE])[seqAnno$gene_name]
    isRef <- dataset[[param$groupingName]] == param$refGroup
    seqAnno[[paste0(count, "_", param$refGroup, "_Mean")]] <- rowMeans(cts[ , isRef, drop=FALSE])[seqAnno$gene_name]
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
  factorLevels <- lapply(factors, function(f) {
    levels(as.factor(dataset[[f]]))
  }) %>% setNames(factors)
  
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


# Extract the dataset path from the URL string 
queryList = parseQueryString(session$clientData$url_search) 
if (is.list(queryList)){
  dataUrl <- queryList$data
} else {
  dataUrl <- NULL
}

# whoami?
message(ezTime(), " app:exploreDE; ", "username:", username, "; ", "dataUrl:", dataUrl)

# complete the full file path for both proteomics and genomics servers 
if (!is.null(dataUrl)) {
  if (!grepl("Proteomics|prolfqua", dataUrl, ignore.case = TRUE)) {
    urlDataRoot = c("/srv/gstore/projects", "/srv/GT/analysis/course_sushi/public/gstore/projects")
    dataDir <- file.path(urlDataRoot, dataUrl)
    dataDir <- dataDir[file.exists(dataDir)][1]
    if (!file.exists(dataDir)){
      showModal(modalDialog(
        title = "Something went wrong",
        "It looks like either the dataset you're looking for doesn't exist, or has not finished being processed in SUSHI yet."
      ))
      Sys.sleep(30)
      stopApp(returnValue = invisible())
    }
  } else {
    dataDir <- paste0("https://fgcz-ms.uzh.ch/public/pStore/", dataUrl)
  }
  projectFromUrl <- regmatches(dataUrl, regexec("p[0-9]{4,}", dataUrl))[[1]][1] # JLR 2025
} else if (is.null(dataUrl) & !exists("fileSE")) {
  dataDir <- "/srv/gstore/projects/p3009/o5638_DESeq2_diff--over--undiff_2024-11-20--12-29-40/diff--over--undiff/"
  projectFromUrl <- "p3009" # JLR 2025
  showNotification("Since you did not specify a dataset in the URL, you are seeing a demo dataset.", type = "message", duration = NULL, closeButton = TRUE)
}

if(exists("fileSE")) {
  dataDir <- fileSE
}
if(!exists("dataDir")) {
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
if (grepl("rds", dataDir) & grepl("Proteomics|prolfqua", dataDir) || grepl("rds", dataDir) & exists("fileSE")) {
  if (is_url(dataDir)) {
    se <- readRDS(url(dataDir))  # If it's a URL
  } else {
    se <- readRDS(dataDir)       # If it's a local file path
  }
} else if (grepl("zip", dataDir)) {
  myTempFile <- tempfile(fileext=".rds", tmpdir = ".") #
  res = system(
    paste("ssh fgcz-ms.uzh.ch unzip -c -qq /srv/www/htdocs/p30385/bfabric/Proteomics/DEA_FragPipe-DiaNN/2023/2023-09/2023-09-22/workunit_293901/2360252.zip '\\*rds' > ", myTempFile))
  if (res != 0){
    message("download failed")
  }
  se <- readRDS(myTempFile)
}
if (grepl("rds|zip", dataDir) & grepl("Proteomics|prolfqua", dataDir) || grepl("rds|zip", dataDir) & exists("fileSE")) {
  contrasts <- names(rowData(se))[grep("^constrast_", names(rowData(se)))] %>% gsub("constrast_", "", .)
  output$proteomicsContrastSelectorUI <- renderUI({
    selectInput(inputId = "contrastSelected", label = "Select contrast to view", choices = contrasts, selected = contrasts[1], multiple = F, selectize = T)
  })
}


# Import RNA seq data from SUSHI ----
if (grepl("gstore", dataDir) | exists("fileSE")) {
  roles <- system(paste0("ldapsearch -x -H ldaps://fgcz-bfabric-ldap:636 -b 'dc=bfabric,dc=org' '(cn=",username,")' memberof | grep Roles | sed 's/,ou=.*//g;s,.*cn=,,g'"),intern=T) 
  allowedProjects <- system(paste0("ldapsearch -x -H ldaps://fgcz-bfabric-ldap:636 -b 'dc=bfabric,dc=org' '(cn=",username,")' memberof | grep Projec | sed 's/,ou=.*//g;s,.*cn=P_,p,g' | sort | uniq"),intern=T)
  if ("R_2" %in% roles){
    ldap_role <- "employee"
    allowed <- TRUE
  } else if ("R_3" %in% roles){
    ldap_role <- "user"
    allowed <- ifelse(projectFromUrl %in% allowedProjects, TRUE, FALSE)
  } else {
    allowed <- FALSE
  }
  message("Role: ", roles, "; LDAP: ", ldap_role, "; Allowed to all projects?: ", allowed, "; Allowed projects: ", allowedProjects)
  if (is.null(dataUrl)) {
    allowed = TRUE
  }
  if (allowed) {
    if (grepl("EzResult.RData", dataDir)) {
      dataDir <- gsub("\\/result-.*.-EzResult.RData", "", dataDir)
    }
    if (file.exists(file.path(dataDir, "deResult.rds"))) {
      se <- readRDS(file.path(dataDir, "deResult.rds"))
    } else if (file.exists(file.path(dataDir, "deResult.qs2"))) {
      se <- qs2::qs_read(file.path(dataDir, "deResult.qs2"), nthreads = 4)
    } else {
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
  waiter <- waiter::Waiter$new(fadeout = TRUE, color = "#86a5bf")
  waiter$show()
  on.exit(waiter$hide())
  
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
figuresDataReactive <- reactiveValues(
  "volcanoStatic" = NULL, "volcanoMA" = NULL, "volcanoPlotly" = NULL,
  "heatmapDEFeatures" = NULL, "heatmapCustom" = NULL, "heatmapGO" = NULL,
  "pcaStatic" = NULL, "pcaPaired" = NULL,
  "boxplotStatic" = NULL, "barplotStatic" = NULL,
  "goTilePlot" = NULL
)
genesReactive <- reactiveValues(genes = NULL)

output$test <- renderText({
  inputDataReactive()$dataType
})
