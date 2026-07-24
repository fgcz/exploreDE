# ---------------------------------------------------------------------------
# Phase 0/1 bootstrap (added 2026-05-15)
#
# The ExploreDataset R6 hierarchy + legacy_gStoreList_from_eds() shim live in
# the package proper (exploreDE/R/), but the deployed app runs against the
# package version baked into fgcz-base-image — which won't include these until
# the image is rebuilt. Until then we source the source-tree R/ files directly
# so the running app picks up the new infrastructure without an image rebuild.
#
# Once the rebuilt image deploys, `exists()` finds the package-provided
# new_explore_dataset() and the loop becomes a no-op.
# ---------------------------------------------------------------------------
if (!exists("new_explore_dataset", mode = "function")) {
  .eds_R_dir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
  for (.f in c(
    "utils_explore_dataset.R",
    "ExploreDataset.R",
    "ExploreDatasetSE.R",
    "ExploreDatasetSE_RNASeq.R",
    "ExploreDatasetSE_CountQC.R",
    "ExploreDatasetSE_CrisprCountQC.R",
    "ExploreDatasetSE_Proteomics.R",
    "new_explore_dataset.R",
    "validate_explore_dataset.R",
    "legacy_glue.R"   # multi-* proteomics auto-load still needs the shim
  )) {
    .fp <- file.path(.eds_R_dir, .f)
    if (file.exists(.fp)) source(.fp)
  }
  rm(.eds_R_dir, .f, .fp)
}


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
  dataDir <- "/srv/gstore/projects/p3009/o5638_DESeq2_diff--over--undiff_2024-11-20--12-29-40/diff--over--undiff/" # teaching dataset
  # Infer dataType from whichever demo dataDir was last (un)commented above —
  # avoids the "I uncommented mageck but forgot to flip dataType" footgun.
  dataType <- if (grepl("Proteomics|prolfqua", dataDir, ignore.case = TRUE)) {
    "proteomics"
  } else {
    "gStore"
  }
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
      # Also load the unfiltered SE (pre-isValid) so the Features Above
      # Threshold bar plot can use the full transcriptome as denominator
      # (matching the ezRun CountQC.Rmd report).
      seUnfiltered <- NULL
      if (file.exists(file.path(dataDir, "rawData.rds"))) {
        seUnfiltered <- readRDS(file.path(dataDir, "rawData.rds"))
      } else if (file.exists(file.path(dataDir, "rawData.qs2"))) {
        seUnfiltered <- qs2::qs_read(file.path(dataDir, "rawData.qs2"), nthreads = 4)
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
# The reactive now exposes the new ExploreDataset hierarchy directly via $.eds.
# Two legacy carry-overs remain:
#   - $seUnfiltered: unfiltered CountQC SE for the "Features Above Threshold"
#       barplot in server-countQc.R (the new eds API has no equivalent assay yet).
#   - attr(., "dataDir"): consumed by server-summary.R for methods text generation
#       (aligner / DE method detection from SUSHI directory names).
inputDataReactive <- reactive({
  waiter <- waiter::Waiter$new(fadeout = TRUE, color = "#527085")
  waiter$show()
  on.exit(waiter$hide())

  res <- list(.eds = new_explore_dataset(se))
  if (exists("seUnfiltered", inherits = TRUE) && !is.null(seUnfiltered)) {
    res$seUnfiltered <- seUnfiltered
  }
  attr(res, "dataDir") <- if (exists("dataDir")) dataDir else NULL
  res
})
# Force eager evaluation at session start so the eds is built before any tab renders
invisible(inputDataReactive())

figuresDataReactive <- reactiveValues(
  "volcanoStatic" = NULL, "volcanoMA" = NULL, "volcanoPlotly" = NULL,
  "heatmapDEFeatures" = NULL, "heatmapCustom" = NULL, "heatmapGO" = NULL,
  "pcaStatic" = NULL, "pcaPaired" = NULL,
  "boxplotStatic" = NULL, "barplotStatic" = NULL,
  "goTilePlot" = NULL
)
# genesReactive removed — replaced by featureBucket in server-featureBucket.R

output$test <- renderText({
  inputDataReactive()$.eds$data_type()
})
