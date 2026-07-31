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


# ---------------------------------------------------------------------------
# Interactive data loading (M2, 2026-07-30)
#
# Data can arrive two ways: (1) a ?data= SUSHI URL, resolved + loaded eagerly at
# session start exactly as before; or (2) a user upload / demo button, mid-
# session (server-landing.R). Both funnel into ONE nullable reactive `edsState`,
# so every downstream `inputDataReactive()$.eds` read is unchanged — it just
# req()-suspends until data is present, and the landing screen shows meanwhile.
#
# Downstream contract preserved: $.eds (always), $seUnfiltered (CountQC only),
# attr(., "dataDir") (server-summary methods text). Uploads supply .eds only.
# ---------------------------------------------------------------------------

edsState <- reactiveVal(NULL)

inputDataReactive <- reactive({
  req(edsState())
  edsState()
})

# The contrast currently selected in the UI: the shared `contrastSelected`
# selector (rendered for proteomics in load_sushi_se, and for multi-contrast
# non-proteomics data in server-summary.R), falling back to the dataset's
# primary contrast. DE + enrichment modules read this so their tabs follow the
# selector whenever a dataset carries more than one contrast; for single-contrast
# data it is a constant (== primary_contrast), so behaviour is unchanged.
selectedContrast <- reactive({
  eds <- inputDataReactive()$.eds
  # Counts-only datasets carry no contrasts; return NULL so DE-family
  # observeEvent(selectedContrast()) blocks (ignoreNULL = TRUE) never fire.
  if (!eds$has_de()) return(NULL)
  cs  <- eds$contrasts()
  sel <- input$contrastSelected
  if (!is.null(sel) && nzchar(sel) && sel %in% cs) sel else eds$primary_contrast()
})

is_url <- function(dataDir) {
  return(grepl("^https?://", dataDir))
}

# Publish an already-built eds into the shared reactive.
set_eds <- function(eds, seUnfiltered = NULL, dataDir = NULL) {
  res <- list(.eds = eds)
  if (!is.null(seUnfiltered)) res$seUnfiltered <- seUnfiltered
  attr(res, "dataDir") <- dataDir
  edsState(res)

  # Proteomics carries no primary contrast — DE + heatmap modules follow the
  # `contrastSelected` picker. load_sushi_se() renders this for the SUSHI URL
  # path; render it here from the eds too so ANY load path (h5ad upload
  # included) populates the selector. Without it input$contrastSelected stays
  # NULL and every DE tab / heatmap feature picker starves.
  if (identical(eds$omics_class(), "proteomics") && length(eds$contrasts())) {
    contrasts <- eds$contrasts()
    output$proteomicsContrastSelectorUI <- renderUI({
      selectInput(inputId = "contrastSelected", label = "Select contrast to view",
                  choices = contrasts, selected = contrasts[1], multiple = FALSE, selectize = TRUE)
    })
    outputOptions(output, "proteomicsContrastSelectorUI", suspendWhenHidden = FALSE)
  }
}

# Wrap a loaded SE into an eds and publish it. The full-screen spinner is managed
# by the landing observers + the deferred module-source observer (functional
# waiter_show/hide), so no separate waiter here (it would layer a 2nd overlay).
set_eds_from_se <- function(se, seUnfiltered = NULL, dataDir = NULL) {
  set_eds(new_explore_dataset(se), seUnfiltered, dataDir)
}

# Resolve a ?data= value into list(dataDir, dataType, projectFromUrl); NULL
# (after a modal) when it can't be resolved.
resolve_from_url <- function(dataUrl) {
  proj <- regmatches(dataUrl, regexec("p[0-9]{4,}", dataUrl))[[1]][1] # JLR 2025
  if (grepl("Proteomics|prolfqua", dataUrl, ignore.case = TRUE)) {
    return(list(dataDir = paste0("https://fgcz-ms.uzh.ch/public/pStore/", dataUrl),
                dataType = "proteomics", projectFromUrl = proj))
  }
  if (grepl("CountQC|DESeq2|EdgeR|Limma|ExploreMageckCounts", dataUrl, ignore.case = TRUE)) {
    urlDataRoot <- c("/srv/gstore/projects", "/srv/GT/analysis/course_sushi/public/gstore/projects")
    dataDir <- file.path(urlDataRoot, dataUrl)
    dataDir <- dataDir[file.exists(dataDir)][1]
    if (is.na(dataDir) || !file.exists(dataDir)) {
      showModal(modalDialog(
        title = "Something went wrong",
        "It looks like either the dataset you're looking for doesn't exist, or has not finished being processed in SUSHI yet."
      ))
      return(NULL)
    }
    return(list(dataDir = dataDir, dataType = "gStore", projectFromUrl = proj))
  }
  showModal(modalDialog(
    title = "The file does not exist",
    "Either the analysis has not yet finished running, you have made a mistake in the URL, or you have not pointed to any dataset. Please try again! If the issue persists, email peter.leary@uzh.ch",
    easyClose = TRUE, footer = NULL
  ))
  NULL
}

# Record project dir for multi-analysis browsing.
record_project <- function(projectFromUrl) {
  multiProjectInfo$currentDataDir <- NULL
  if (is.null(projectFromUrl) || is.na(projectFromUrl) || projectFromUrl == "") return(invisible())
  projectPath <- file.path("/srv/gstore/projects", projectFromUrl)
  if (dir.exists(projectPath)) {
    multiProjectInfo$projectDir <- projectPath
  } else {
    projectPath <- file.path("/srv/GT/analysis/course_sushi/public/gstore/projects", projectFromUrl)
    if (dir.exists(projectPath)) multiProjectInfo$projectDir <- projectPath
  }
}

# Load a SUSHI SE (proteomics or genomics) from a resolved dataDir.
# Returns list(se, seUnfiltered, dataDir) or NULL (after a modal). LDAP checks
# and every load branch are preserved verbatim from the original eager block.
load_sushi_se <- function(dataDir, dataType, projectFromUrl) {

  # ---- Import proteomics data from pStore ----
  if (dataType == "proteomics") {
    se <- if (is_url(dataDir)) readRDS(url(dataDir)) else readRDS(dataDir)
    contrasts <- names(rowData(se))[grep("^constrast_", names(rowData(se)))] %>% gsub("constrast_", "", .)
    output$proteomicsContrastSelectorUI <- renderUI({
      selectInput(inputId = "contrastSelected", label = "Select contrast to view", choices = contrasts, selected = contrasts[1], multiple = F, selectize = T)
    })
    outputOptions(output, "proteomicsContrastSelectorUI", suspendWhenHidden = FALSE)
    return(list(se = se, seUnfiltered = NULL, dataDir = dataDir))
  }
  if (grepl("zip", dataDir)) {
    myTempFile <- tempfile(fileext = ".rds", tmpdir = ".") #
    res <- system(
      paste("ssh fgcz-ms.uzh.ch unzip -c -qq /srv/www/htdocs/p30385/bfabric/Proteomics/DEA_FragPipe-DiaNN/2023/2023-09/2023-09-22/workunit_293901/2360252.zip '\\*rds' > ", myTempFile)
    )
    if (res != 0) {
      message("download failed")
    }
    se <- readRDS(myTempFile)
    contrasts <- names(rowData(se))[grep("^constrast_", names(rowData(se)))] %>% gsub("constrast_", "", .)
    output$proteomicsContrastSelectorUI <- renderUI({
      selectInput(inputId = "contrastSelected", label = "Select contrast to view", choices = contrasts, selected = contrasts[1], multiple = F, selectize = T)
    })
    outputOptions(output, "proteomicsContrastSelectorUI", suspendWhenHidden = FALSE)
    return(list(se = se, seUnfiltered = NULL, dataDir = dataDir))
  }

  # ---- Import RNA-seq / CountQC / Mageck data from SUSHI ----
  # Run LDAP queries — demo + course skip checks.
  if (dataType == "demo") {
    message("Loading demo data, so no LDAP checks required.")
    allowed <- TRUE
    ldap_role <- "demo"
  } else if (grepl("/srv/GT/analysis/course_sushi/public/gstore/projects", dataDir)) {
    message("Loading SUSHI course data, so no LDAP checks required.")
    allowed <- TRUE
    ldap_role <- "course"
  } else {
    roles <- system(paste0("ldapsearch -x -H ldaps://fgcz-bfabric-ldap:636 -b 'dc=bfabric,dc=org' '(cn=", username, ")' memberof | grep Roles | sed 's/,ou=.*//g;s,.*cn=,,g'"), intern = T)
    allowedProjects <- system(paste0("ldapsearch -x -H ldaps://fgcz-bfabric-ldap:636 -b 'dc=bfabric,dc=org' '(cn=", username, ")' memberof | grep Projec | sed 's/,ou=.*//g;s,.*cn=P_,p,g' | sort | uniq"), intern = T)

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

  if (!isTRUE(allowed)) {
    showModal(modalDialog(
      title = "Access denied",
      "You do not have permission to view this project's data. If you believe this is a mistake, email peter.leary@uzh.ch"
    ))
    return(NULL)
  }

  se <- NULL
  seUnfiltered <- NULL
  # If RNA Seq, load that
  if (grepl("DESeq2|EdgeR|Limma", dataDir, ignore.case = TRUE)) {
    if (grepl("EzResult.RData", dataDir)) {
      dataDir <- gsub("\\/result-.*.-EzResult.RData", "", dataDir)
    }
    if (file.exists(file.path(dataDir, "deResult.rds"))) {
      se <- readRDS(file.path(dataDir, "deResult.rds"))
    } else if (file.exists(file.path(dataDir, "deResult.qs2"))) {
      se <- qs2::qs_read(file.path(dataDir, "deResult.qs2"), nthreads = 4)
    }
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
    if (file.exists(file.path(dataDir, "rawData.rds"))) {
      seUnfiltered <- readRDS(file.path(dataDir, "rawData.rds"))
    } else if (file.exists(file.path(dataDir, "rawData.qs2"))) {
      seUnfiltered <- qs2::qs_read(file.path(dataDir, "rawData.qs2"), nthreads = 4)
    }
  }

  if (is.null(se)) {
    showModal(modalDialog(
      title = "The file does not exist",
      "Either the analysis has not yet finished running, you have made a mistake in the URL, or you have not pointed to any dataset. Please try again! If the issue persists, email peter.leary@uzh.ch",
      easyClose = TRUE,
      footer = NULL
    ))
    return(NULL)
  }
  list(se = se, seUnfiltered = seUnfiltered, dataDir = dataDir)
}

# ---- Entry points (shared by URL boot, the demo button, and uploads) --------

# Resolve + load a ?data= URL, then publish the eds.
load_from_url <- function(dataUrl) {
  resolved <- resolve_from_url(dataUrl)
  if (is.null(resolved)) return(invisible())
  record_project(resolved$projectFromUrl)
  loaded <- load_sushi_se(resolved$dataDir, resolved$dataType, resolved$projectFromUrl)
  if (is.null(loaded)) return(invisible())
  set_eds_from_se(loaded$se, loaded$seUnfiltered, loaded$dataDir)
}

# The demo dataset (no ?data=, triggered by the landing-screen button).
# A directory (as a ?data= URL resolves to), NOT the deResult file — load_sushi_se
# appends deResult.rds/.qs2 itself for the DESeq2/EdgeR/Limma branch.
DEMO_DATADIR <- "/srv/gstore/projects/p3009/o5638_EdgeR_diff--over--undiff_2026-07-30--12-08-00/diff--over--undiff"
load_demo <- function() {
  record_project("p3009")
  loaded <- load_sushi_se(DEMO_DATADIR, "demo", "p3009")
  if (is.null(loaded)) return(invisible())
  set_eds_from_se(loaded$se, loaded$seUnfiltered, loaded$dataDir)
}

# A user-uploaded canonical .h5ad. Returns list(ok, msg) so the landing screen
# can surface a friendly validation message; publishes the eds only on success.
load_h5ad_upload <- function(path) {
  eds <- tryCatch(new_explore_dataset(path), error = function(e) e)
  if (inherits(eds, "error")) {
    m <- conditionMessage(eds)
    hint <- if (grepl("HDF5|link|iteration|truncat|superblock|addr", m, ignore.case = TRUE)) {
      paste0("\n\nThe file could not be opened as HDF5. It may be corrupt or incompletely ",
             "written — e.g. an interrupted or out-of-memory build. Try rebuilding it in a ",
             "fresh R session (plain Rscript, rm() large objects + gc()).")
    } else ""
    return(list(ok = FALSE, msg = paste0("Could not read the file: ", m, hint)))
  }
  v <- tryCatch(validate_explore_dataset(eds, strict = FALSE),
                error = function(e) list(ok = FALSE, errors = conditionMessage(e)))
  if (!isTRUE(v$ok)) {
    return(list(ok = FALSE, msg = paste(v$errors, collapse = "\n")))
  }
  set_eds(eds)
  list(ok = TRUE, msg = NULL)
}

# ---- Boot: resolve + load a ?data= URL eagerly (unchanged UX for SUSHI) ------
queryList <- parseQueryString(session$clientData$url_search)
dataUrl   <- if (is.list(queryList)) queryList$data else NULL
message(ezTime(), " app:exploreDE; ", "username:", username, "; ", "dataUrl:", dataUrl)

if (!is.null(dataUrl)) {
  # `dataType` at module scope for the telemetry-app-name check in app.R.
  dataType <- if (grepl("Proteomics|prolfqua", dataUrl, ignore.case = TRUE)) "proteomics" else "gStore"
  load_from_url(dataUrl)
}
# No ?data= → edsState stays NULL → the landing screen (server-landing.R) shows.

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
