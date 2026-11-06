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

# ---- App-wide contrast (linked per-tab selectors) ---------------------------
# ONE server-side value, `activeContrast`, is the single source of truth. Every
# contrast-dependent tab carries its own selector (Shiny input ids must be
# unique): the canonical `contrastSelected` on the Summary tab, plus mirrors
# `contrastSelected_<tab>`. Changing any of them writes activeContrast(); one
# observer pushes the value back out to all the others, so a change made on
# Volcano is already in place when Heatmap is opened. The loop terminates: a
# reactiveVal set to its current value doesn't invalidate, and a selectize
# updated to its current value doesn't fire. Save/Load state only snapshots the
# canonical id (mirrors are dropped in save_state_filter_inputs()).
# Selectors render only for datasets with more than one contrast.
activeContrast <- reactiveVal(NULL)

contrast_mirror_tabs <- c(
  "degTable",
  "volcano",
  "deQc",
  "heatmap",
  "boxplot",
  "ora",
  "gsea",
  "goTile",
  "kegg",
  "correlation",
  "replicability"
)
contrast_selector_ids <- c(
  summary = "contrastSelected",
  stats::setNames(paste0("contrastSelected_", contrast_mirror_tabs), contrast_mirror_tabs)
)

# The contrast currently selected: activeContrast(), validated against the
# dataset and falling back to its primary contrast. Every contrast-dependent
# module reads this; for single-contrast data it is a constant
# (== primary_contrast), so behaviour is unchanged there.
selectedContrast <- reactive({
  eds <- inputDataReactive()$.eds
  # Counts-only datasets carry no contrasts; return NULL so DE-family
  # observeEvent(selectedContrast()) blocks (ignoreNULL = TRUE) never fire.
  if (!eds$has_de()) {
    return(NULL)
  }
  cs <- eds$contrasts()
  sel <- activeContrast()
  if (!is.null(sel) && nzchar(sel) && sel %in% cs) sel else eds$primary_contrast()
})

base::lapply(names(contrast_selector_ids), function(tab) {
  id <- contrast_selector_ids[[tab]]
  output[[paste0("contrastSelectorUI_", tab)]] <- renderUI({
    eds <- inputDataReactive()$.eds
    cs <- eds$contrasts()
    if (length(cs) < 2) {
      return(NULL)
    }
    selectizeInput(
      inputId = id,
      label = "Select contrast to view",
      choices = cs,
      selected = isolate(selectedContrast()),
      multiple = FALSE
    )
  })
  outputOptions(output, paste0("contrastSelectorUI_", tab), suspendWhenHidden = FALSE)
  observeEvent(
    input[[id]],
    {
      v <- input[[id]]
      if (nzchar(v) && !identical(v, isolate(activeContrast()))) activeContrast(v)
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
  )
})

observeEvent(
  activeContrast(),
  {
    v <- activeContrast()
    for (id in contrast_selector_ids) {
      if (!identical(isolate(input[[id]]), v)) {
        updateSelectizeInput(session = session, inputId = id, selected = v)
      }
    }
  },
  ignoreNULL = TRUE
)

# Publish an already-built eds into the shared reactive.
set_eds <- function(eds, seUnfiltered = NULL, dataDir = NULL) {
  res <- list(.eds = eds)
  if (!is.null(seUnfiltered)) {
    res$seUnfiltered <- seUnfiltered
  }
  attr(res, "dataDir") <- dataDir
  edsState(res)

  # Reset the app-wide contrast for the new dataset. Deliberately NOT seeded to
  # the primary contrast: selectedContrast() already falls back to it, and the
  # first selector to report its value sets activeContrast(), which invalidates
  # selectedContrast() *after* the modules are sourced. Several modules rely on
  # that late invalidation for their first render (observeEvent(..., ignoreInit
  # = TRUE)), exactly as they did when the old selector's input first arrived.
  activeContrast(NULL)
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
  # The LDAP check below uses the project named in the link, so a ".." could
  # name an allowed project and load another one. Real links never contain it.
  if (any(strsplit(dataUrl, "/", fixed = TRUE)[[1]] == "..")) {
    showModal(modalDialog(
      title = "The file does not exist",
      "This link is not a valid dataset path. If the issue persists, email peter.leary@uzh.ch",
      easyClose = TRUE,
      footer = NULL
    ))
    return(NULL)
  }
  proj <- regmatches(dataUrl, regexec("p[0-9]{4,}", dataUrl))[[1]][1] # JLR 2025
  if (grepl("Proteomics|prolfqua", dataUrl, ignore.case = TRUE)) {
    return(list(
      dataDir = paste0("https://fgcz-ms.uzh.ch/public/pStore/", dataUrl),
      dataType = "proteomics",
      projectFromUrl = proj
    ))
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
    easyClose = TRUE,
    footer = NULL
  ))
  NULL
}

# Record the project dir (and the loaded dataset's dir) for multi-analysis
# browsing; NULL forgets both. Call only after the dataset has passed the
# access check, so a denied project never becomes browsable.
record_project <- function(projectFromUrl, dataDir = NULL) {
  multiProjectInfo$projectDir <- NULL
  multiProjectInfo$currentDataDir <- dataDir
  if (is.null(projectFromUrl) || is.na(projectFromUrl) || projectFromUrl == "") {
    return(invisible())
  }
  projectPath <- file.path("/srv/gstore/projects", projectFromUrl)
  if (dir.exists(projectPath)) {
    multiProjectInfo$projectDir <- projectPath
  } else {
    projectPath <- file.path("/srv/GT/analysis/course_sushi/public/gstore/projects", projectFromUrl)
    if (dir.exists(projectPath)) multiProjectInfo$projectDir <- projectPath
  }
}

# The user's B-Fabric roles (R_2 = employee, R_3 = user) and projects, from
# LDAP. Looked up once per session and shared by the ?data= check and every
# multi-DE dataset (ldap_allows_project() in R/utils_multi_de.R).
ldap_lookup <- function(user) {
  query <- paste0(
    "ldapsearch -x -H ldaps://fgcz-bfabric-ldap:636 -b 'dc=bfabric,dc=org' ",
    shQuote(paste0("(cn=", user, ")")),
    " memberof"
  )
  list(
    roles = system(paste0(query, " | grep Roles | sed 's/,ou=.*//g;s,.*cn=,,g'"), intern = TRUE),
    projects = system(paste0(query, " | grep Projec | sed 's/,ou=.*//g;s,.*cn=P_,p,g' | sort | uniq"), intern = TRUE)
  )
}
.ldap_cache <- NULL
session_ldap <- function() {
  if (is.null(.ldap_cache)) {
    .ldap_cache <<- ldap_lookup(username)
  }
  .ldap_cache
}

# Load a SUSHI SE (proteomics or genomics) from a resolved dataDir.
# Returns list(se, seUnfiltered, dataDir) or NULL (after a modal). LDAP checks
# and every load branch are preserved verbatim from the original eager block.
load_sushi_se <- function(dataDir, dataType, projectFromUrl) {
  # ---- Import proteomics data from pStore ----
  if (dataType == "proteomics") {
    se <- if (grepl("\\.h5ad$", dataDir, ignore.case = TRUE)) {
      # prolfquapp DEA AnnData; the reader needs a local file.
      path <- dataDir
      if (is_url(dataDir)) {
        path <- tempfile(fileext = ".h5ad")
        on.exit(unlink(path), add = TRUE)
        utils::download.file(dataDir, path, mode = "wb", quiet = TRUE)
      }
      se_from_prolfquapp_anndata(.h5ad_read(path))
    } else if (is_url(dataDir)) {
      readRDS(url(dataDir))
    } else {
      readRDS(dataDir)
    }
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
    lk <- session_ldap()
    roles <- lk$roles
    allowedProjects <- lk$projects

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
    message(
      "Role: ",
      roles,
      "; LDAP: ",
      ldap_role,
      "; Allowed to all projects?: ",
      allowed,
      "; Allowed projects: ",
      allowedProjects
    )
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
  # A full link (?data=https://...) to an .h5ad on an allowed site.
  if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://", dataUrl)) {
    return(load_h5ad_from_link(dataUrl))
  }
  resolved <- resolve_from_url(dataUrl)
  if (is.null(resolved)) {
    return(invisible())
  }
  loaded <- load_sushi_se(resolved$dataDir, resolved$dataType, resolved$projectFromUrl)
  if (is.null(loaded)) {
    return(invisible())
  }
  record_project(resolved$projectFromUrl, loaded$dataDir)
  set_eds_from_se(loaded$se, loaded$seUnfiltered, loaded$dataDir)
}

# ?data=https://... : Download the .h5ad (allowed sites only, no redirects,
# 1 GB limit; R/utils_url_load.R), then the same checks as an upload. No LDAP
# check: Whoever has the link can already download the file.
load_h5ad_from_link <- function(url) {
  record_project(NULL)
  fail <- function(title, msg) {
    showModal(modalDialog(title = title, msg, easyClose = TRUE))
    invisible()
  }
  ok <- url_load_allowed(url)
  if (!isTRUE(ok)) {
    return(fail("This link can't be opened", paste0("exploreDE can't load it: ", attr(ok, "reason"), ".")))
  }
  path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(path), add = TRUE)
  fetched <- tryCatch(url_load_fetch(url, path), error = function(e) e)
  if (inherits(fetched, "error")) {
    return(fail(
      "Could not download this file",
      paste0(conditionMessage(fetched), ". Links to finished jobs expire; check that the job still exists.")
    ))
  }
  res <- load_h5ad_upload(path)
  if (!isTRUE(res$ok)) {
    fail("Could not load this file", res$msg)
  }
  invisible()
}

# The demo dataset (no ?data=, triggered by the landing-screen button): the
# SUSHI RNA-seq result shipped with the package, so the demo works anywhere, not
# only where /srv/gstore is mounted. The app always runs with inst/app (or an
# installed package's app/) as its working directory, so the examples are at
# ../examples in the source tree, the container and an installed package alike.
DEMO_DATADIR <- normalizePath("../examples/genomics", mustWork = FALSE)
load_demo <- function() {
  record_project(NULL)
  f <- file.path(DEMO_DATADIR, "deResult.rds")
  if (!file.exists(f)) {
    showModal(modalDialog(title = "Demo data not found", paste("Expected the demo dataset at", f), easyClose = TRUE))
    return(invisible())
  }
  set_eds_from_se(readRDS(f), dataDir = DEMO_DATADIR)
}

# A file handed to exploreDE(file = ...): an .h5ad goes through the same checks
# as an upload; anything else (a SummarizedExperiment .rds/.qs2, or a SUSHI
# result directory) through new_explore_dataset().
load_local_file <- function(path) {
  if (grepl("\\.h5ad$", path, ignore.case = TRUE)) {
    res <- load_h5ad_upload(path)
    if (!isTRUE(res$ok)) {
      showModal(modalDialog(title = "Could not load this file", res$msg, easyClose = TRUE))
    }
    return(invisible())
  }
  eds <- tryCatch(new_explore_dataset(path), error = function(e) e)
  if (inherits(eds, "error")) {
    showModal(modalDialog(title = "Could not load this file", conditionMessage(eds), easyClose = TRUE))
    return(invisible())
  }
  set_eds(eds, dataDir = if (dir.exists(path)) path else dirname(path))
}

# A user-uploaded canonical .h5ad. Returns list(ok, msg) so the landing screen
# can surface a friendly validation message; publishes the eds only on success.
load_h5ad_upload <- function(path) {
  eds <- tryCatch(new_explore_dataset(path), error = function(e) e)
  if (inherits(eds, "error")) {
    m <- conditionMessage(eds)
    hint <- if (grepl("HDF5|link|iteration|truncat|superblock|addr", m, ignore.case = TRUE)) {
      paste0(
        "\n\nThe file could not be opened as HDF5. It may be corrupt or incompletely ",
        "written — e.g. an interrupted or out-of-memory build. Try rebuilding it in a ",
        "fresh R session (plain Rscript, rm() large objects + gc())."
      )
    } else {
      ""
    }
    return(list(ok = FALSE, msg = paste0("Could not read the file: ", m, hint)))
  }
  v <- tryCatch(validate_explore_dataset(eds, strict = FALSE), error = function(e) {
    list(ok = FALSE, errors = conditionMessage(e))
  })
  if (!isTRUE(v$ok)) {
    return(list(ok = FALSE, msg = paste(v$errors, collapse = "\n")))
  }
  set_eds(eds)
  list(ok = TRUE, msg = NULL)
}

# ---- Boot: resolve + load a ?data= URL eagerly (unchanged UX for SUSHI) ------
queryList <- parseQueryString(session$clientData$url_search)
dataUrl <- if (is.list(queryList)) queryList$data else NULL
message(format(Sys.time(), "%Y-%m-%d--%H-%M-%S"), " app:exploreDE; ", "username:", username, "; ", "dataUrl:", dataUrl)

if (!is.null(dataUrl)) {
  # `dataType` at module scope for the telemetry-app-name check in app.R.
  dataType <- if (grepl("Proteomics|prolfqua", dataUrl, ignore.case = TRUE)) "proteomics" else "gStore"
  load_from_url(dataUrl)
} else if (!is.null(getOption("exploreDE.file"))) {
  # Launched locally as exploreDE(file = ...).
  load_local_file(getOption("exploreDE.file"))
}
# Nothing to load → edsState stays NULL → the landing screen (server-landing.R) shows.

figuresDataReactive <- reactiveValues(
  "volcanoStatic" = NULL,
  "volcanoMA" = NULL,
  "volcanoPlotly" = NULL,
  "heatmapDEFeatures" = NULL,
  "heatmapCustom" = NULL,
  "heatmapGO" = NULL,
  "pcaStatic" = NULL,
  "pcaPaired" = NULL,
  "boxplotStatic" = NULL,
  "barplotStatic" = NULL,
  "goTilePlot" = NULL
)
# genesReactive removed — replaced by featureBucket in server-featureBucket.R

output$test <- renderText({
  inputDataReactive()$.eds$data_type()
})
