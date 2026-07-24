# What's going on here?
# The app should first check if there's a URL query leary string, and if so, load that data 
# If not, users should be able to upload their own version of summarizedExperiment/AnnData 
username <- Sys.getenv("SHINYPROXY_USERNAME")
se1 <- readRDS("/srv/gstore/projects/p3009/o5638_DESeq2_diff--over--undiff_2024-11-20--12-29-40/diff--over--undiff/deResult.rds")
se2 <- readRDS("/srv/gstore/projects/p3009/o5638_CountQC_2025-11-14--15-14-04/Count_QC/rawData.rds")
se3 <- readRDS(url("https://fgcz-ms.uzh.ch/public/pStore/p39763/bfabric/Proteomics/exploreDE/2025/2025-11/2025-11-03/workunit_335997/3012890.rds"))

# 1. Get URL Data ----

# Fake session with the query string only
session <- list(
  clientData = list(
    url_search = "?data=p3009/o5638_DESeq2_diff--over--undiff_2024-11-20--12-29-40/diff--over--undiff/" # SUSHI RNA Seq
    # url_search = "?data=p3009/o5638_CountQC_2025-11-14--15-14-04/Count_QC" # SUSHI CountQC
    # url_search = "?data=p39763/bfabric/Proteomics/exploreDE/2025/2025-11/2025-11-03/workunit_335997/3012890.rds" # B-fabric proteomics
  )
)

# Extract the dataset path from the URL string 
queryList <- parseQueryString(session$clientData$url_search)
if (is.list(queryList)) {
  dataUrl <- queryList$data
} else {
  dataUrl <- NULL
}

# Prepare the dataDir variable, from which we'll load the SE object 
if (!is.null(dataUrl)) {
  # If the URL has pX/oX, assume it's from genomics...
  if (grepl("^p\\d{4,5}/o\\d{4,5}", dataUrl)) {
    urlDataRoot = c("/srv/gstore/projects", "/srv/GT/analysis/course_sushi/public/gstore/projects")
    dataDir <- file.path(urlDataRoot, dataUrl)
    dataDir <- dataDir[file.exists(dataDir)][1]
  # Otherwise, if it's clearly labelled proteomics...
  } else if (grepl("Proteomics|prolfqua", dataUrl, ignore.case = TRUE)) {
    dataDir <- paste0("https://fgcz-ms.uzh.ch/public/pStore/", dataUrl)
  }

  # Check the user has permission to access the data
  projectFromUrl <- regmatches(dataUrl, regexec("p[0-9]{4,}", dataUrl))[[1]][1] # JLR 2025
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
  # If loading course data, always allow 
  if (grepl("/srv/GT/analysis/course_sushi/public/gstore/projects", dataDir)) {
    message("Loading SUSHI course data, so no LDAP checks required.")    
    return(list(allowed = TRUE, ldap_role = "course"))  
  }
  
  # For proteomics, load the SE right away...
  is_url <- function(dataDir) {
    return(grepl("^https?://", dataDir)) 
  }
  if (is_url(dataDir)) {
    se <- readRDS(url(dataDir))
  }
  
  # For genomics, find the correct file 
  if (grepl("^p\\d{4,5}/o\\d{4,5}", dataUrl)) {
    # old links point to a specific file
    if (grepl("EzResult.RData", dataDir)) {
      dataDir <- gsub("\\/result-.*.-EzResult.RData", "", dataDir)
    }
    # Or find the correct se for CountQC or newer DE 
    if (grepl("CountQC", dataUrl)) {
      fileName <- "rawData"
    } else {
      fileName <- "deResult"
    }
    # construct the final file path 
    dataDirFinal <- list.files(path = dataDir, pattern = fileName, full.names = TRUE)
    # and load, whether rds or qs2 
    if (file.exists(dataDirFinal)) {
      if (grepl(".rds", dataDirFinal)) {
        se <- readRDS(dataDirFinal)
      } else if (grepl(".qs2", dataDirFinal)) {
        se <- qs2::qs_read(dataDirFinal)
      }
    } else {
      stop("no likey")
    }
  }

}

# 2. Prepare the SE ---- 
# There are four data types: 1. RNA seq, 2. proteomics, 3. CountQC, 4. AnnData
detect_se_type <- function(se) {
  assays <- assayNames(se)
  cdn    <- colnames(colData(se))
  
  # 1. Proteomics: SummarizedExperiment (non-Ranged) with rawData/transformedData
  if (inherits(se, "SummarizedExperiment") && 
      !inherits(se, "RangedSummarizedExperiment") &&
      any(c("rawData", "transformedData") %in% assays)) {
    return("proteomics")
  }
  
  # 2. RNA: RangedSE with xNorm / isPresent / sf
  if (inherits(se, "RangedSummarizedExperiment") &&
      (any(c("xNorm", "isPresent") %in% assays) || "sf" %in% cdn)) {
    return("rna")
  }
  
  # 3. CountQC: RangedSE with counts + presentFlag and none of the RNA extras
  if (inherits(se, "RangedSummarizedExperiment") &&
      all(c("counts", "presentFlag") %in% assays) &&
      !any(c("xNorm", "isPresent") %in% assays) &&
      !"sf" %in% cdn) {
    return("countqc")
  }
  
  # 4. Fallback
  "unknown"
}

prepare_se <- function(se) {
  # What are we dealing with here?
  se_type <- detect_se_type(se)
  message("Detected SE type: ", se_type)
  
  # get the SUSHI params, if they're there
  if (!is.null(metadata(se)$param)) {
    param <- metadata(se)$param
  }
  
  seqAnno <- data.frame(rowData(se), row.names = rownames(se), check.names = FALSE)
  
}

