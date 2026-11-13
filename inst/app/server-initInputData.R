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

# The label of the selected contrast for titles and download file names, used
# by every tab. Defined once here: the data modules share one environment, so
# a copy in a module would replace this one for all of them.
designReactive <- reactive({
  req(inputDataReactive())
  eds <- inputDataReactive()$.eds
  if (eds$has_de()) {
    contrast <- selectedContrast()
    eds$contrast_label(contrast)
  } else if (eds$omics_class() == "proteomics") {
    "No contrast selected"
  } else {
    eds$data_type() # QC fallback — no DE, use the data type name
  }
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
  # No user (SHINYPROXY_USERNAME unset): no roles and no projects, without
  # asking LDAP for "(cn=)".
  if (!nzchar(user)) {
    return(list(roles = character(), projects = character()))
  }
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
# The LDAP role decided by the last ?data= load ("employee", "user", "course",
# ...); NA until then. The Ordination tab shows its staff-only panel for
# "employee".
ldapRole <- reactiveVal(NA_character_)
session_ldap <- function() {
  if (is.null(.ldap_cache)) {
    .ldap_cache <<- ldap_lookup(username)
  }
  .ldap_cache
}

# A proteomics result on pStore: an .rds, or a prolfquapp DEA AnnData.
read_pstore_se <- function(dataDir) {
  if (grepl("\\.h5ad$", dataDir, ignore.case = TRUE)) {
    # The reader needs a local file.
    path <- dataDir
    if (is_url(dataDir)) {
      path <- tempfile(fileext = ".h5ad")
      on.exit(unlink(path), add = TRUE)
      utils::download.file(dataDir, path, mode = "wb", quiet = TRUE)
    }
    return(se_from_prolfquapp_anndata(.h5ad_read(path)))
  }
  if (is_url(dataDir)) readRDS(url(dataDir)) else readRDS(dataDir)
}

# Load a SUSHI SE (proteomics or genomics) from a resolved dataDir.
# Returns list(se, seUnfiltered, dataDir) or NULL (after a modal). LDAP checks
# and every load branch are preserved verbatim from the original eager block.
load_sushi_se <- function(dataDir, dataType, projectFromUrl) {
  # ---- Import proteomics data from pStore ----
  if (dataType == "proteomics") {
    # A link to a missing file, or pStore out of reach: a message, not the end
    # of the session (this runs while the server function starts).
    se <- tryCatch(read_pstore_se(dataDir), error = function(e) e)
    if (inherits(se, "error")) {
      showModal(modalDialog(
        title = "This link can't be opened",
        paste0("exploreDE could not read ", dataDir, ": ", conditionMessage(se)),
        easyClose = TRUE
      ))
      return(NULL)
    }
    return(list(se = se, seUnfiltered = NULL, dataDir = dataDir))
  }

  # ---- Import RNA-seq / CountQC / Mageck data from SUSHI ----
  # Run LDAP queries — demo + course skip checks.
  if (grepl("/srv/GT/analysis/course_sushi/public/gstore/projects", dataDir)) {
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
  ldapRole(ldap_role)

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
    se <- read_se_any(dataDir, "deResult")
    # If CountQC, load that (we use only the 'present' genes because it has the normalised counts in the assays)
  } else if (grepl("CountQC|ExploreMageckCounts", dataDir)) {
    # Current SUSHI links point at the result directory; old CountQC links
    # point at a file inside it (counts-...-EzResult.RData).
    if (!dir.exists(dataDir)) {
      dataDir <- dirname(dataDir)
    }
    se <- read_se_any(dataDir, "rawData_final")
    # Also load the unfiltered SE (pre-isValid) so the Features Above
    # Threshold bar plot can use the full transcriptome as denominator
    # (matching the ezRun CountQC.Rmd report).
    seUnfiltered <- read_se_any(dataDir, "rawData")
    # Apps without a filtered object (ExploreMageckCounts) only write rawData.
    if (is.null(se) && !is.null(seUnfiltered)) {
      se <- seUnfiltered
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
# Returns NULL when the file opened, otherwise why it did not.
load_local_file <- function(path) {
  if (grepl("\\.h5ad$", path, ignore.case = TRUE)) {
    res <- load_h5ad_upload(path)
    return(if (isTRUE(res$ok)) NULL else res$msg)
  }
  eds <- tryCatch(new_explore_dataset(path), error = function(e) e)
  if (inherits(eds, "error")) {
    return(conditionMessage(eds))
  }
  set_eds(eds, dataDir = if (dir.exists(path)) path else dirname(path))
  NULL
}

# A user-uploaded canonical .h5ad, read and validated but not published.
# Returns list(ok, msg, eds); msg is a friendly validation message.
read_h5ad_upload <- function(path) {
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
  list(ok = TRUE, msg = NULL, eds = eds)
}

# A user-uploaded canonical .h5ad. Returns list(ok, msg) so the landing screen
# can surface a friendly validation message; publishes the eds only on success.
# The rows of a fileInput value whose names end in .h5ad; the others are
# refused with a message. readRDS() on an uploaded .rds could run code, so the
# upload boxes take .h5ad only (exploreDE(file) still opens a local .rds).
h5ad_uploads_only <- function(files) {
  ok <- grepl("\\.h5ad$", files$name, ignore.case = TRUE)
  for (name in files$name[!ok]) {
    showNotification(
      paste0(name, " is not an .h5ad file; only .h5ad files can be uploaded."),
      type = "error",
      duration = 10
    )
  }
  files[ok, , drop = FALSE]
}

load_h5ad_upload <- function(path) {
  res <- read_h5ad_upload(path)
  if (!isTRUE(res$ok)) {
    return(res[c("ok", "msg")])
  }
  set_eds(res$eds)
  list(ok = TRUE, msg = NULL)
}

# Report the contrasts multi_records_from_dataset() could not build.
multi_notify_failed <- function(records, what) {
  failed <- attr(records, "failed")
  for (con in names(failed)) {
    showNotification(
      paste0("Could not load contrast '", con, "' of ", what, ": ", failed[[con]]),
      type = "warning",
      duration = 10
    )
  }
}

# Files added for Multi DE: the landing page's extra files, the Data Selection
# tab's upload, or the extra paths of exploreDE(c(...)). The contrasts of each
# become Multi DE records named "<file>_<contrast>", the file's label being
# unique among the files still in Multi DE. When Multi DE is empty, the loaded
# dataset's own contrasts join first, under their own names. `files` is a
# fileInput value (columns name, datapath). Returns list(contrasts, files):
# how many were added; every file or contrast that could not be used is
# reported.
multi_add_uploads <- function(files) {
  none <- list(contrasts = 0L, files = 0L)
  main <- isolate(edsState())$.eds
  if (is.null(main)) {
    return(none)
  }
  current <- isolate(multiDeReactive$results)
  labels <- unlist(lapply(current, function(x) x[["source"]][["file"]]), use.names = FALSE)
  records <- list()
  n_files <- 0L
  for (i in seq_len(nrow(files))) {
    name <- files$name[i]
    res <- read_h5ad_upload(files$datapath[i])
    if (!isTRUE(res$ok)) {
      showNotification(paste0("Could not load ", name, ": ", res$msg), type = "error", duration = 10)
      next
    }
    eds <- res$eds
    if (!identical(eds$omics_class(), main$omics_class())) {
      showNotification(
        paste0(
          name,
          " holds ",
          eds$omics_class(),
          " data, so it cannot be compared with the loaded ",
          main$omics_class(),
          " dataset."
        ),
        type = "error",
        duration = 10
      )
      next
    }
    if (!eds$has_de()) {
      showNotification(
        paste0(name, " has no DE contrasts, so there is nothing to compare."),
        type = "warning",
        duration = 10
      )
      next
    }
    label <- multi_upload_label(name, labels)
    recs <- multi_records_from_dataset(eds, label)
    multi_notify_failed(recs, name)
    if (!length(recs)) {
      next
    }
    labels <- c(labels, label)
    records <- c(records, recs)
    n_files <- n_files + 1L
  }
  if (!length(records)) {
    return(none)
  }
  if (!length(current) && main$has_de()) {
    own <- multi_records_from_dataset(main)
    multi_notify_failed(own, "the loaded dataset")
    current <- c(current, own)
  }
  # A record name is an input id too, so it must be unique among the datasets.
  results <- c(current, records)
  names(results) <- make.unique(names(results), sep = "_")
  multiDeReactive$results <- results
  multiDeReactive$menuUpdateTrigger <- isolate(multiDeReactive$menuUpdateTrigger) + 1
  list(contrasts = length(records), files = n_files)
}

# ---- Boot: resolve + load a ?data= URL eagerly (unchanged UX for SUSHI) ------
queryList <- parseQueryString(isolate(session$clientData$url_search))
dataUrl <- if (is.list(queryList)) queryList$data else NULL
message(format(Sys.time(), "%Y-%m-%d--%H-%M-%S"), " app:exploreDE; ", "username:", username, "; ", "dataUrl:", dataUrl)

if (!is.null(dataUrl)) {
  # `dataType` at module scope for the telemetry-app-name check in app.R.
  dataType <- if (grepl("Proteomics|prolfqua", dataUrl, ignore.case = TRUE)) "proteomics" else "gStore"
  load_from_url(dataUrl)
} else if (!is.null(getOption("exploreDE.file"))) {
  # Launched locally as exploreDE(file = ...): the first file that loads opens,
  # the ones after it go to Multi DE (as on the landing page).
  files <- getOption("exploreDE.file")
  failed <- character()
  opened <- 0L
  for (i in seq_along(files)) {
    why <- load_local_file(files[i])
    if (is.null(why)) {
      opened <- i
      break
    }
    failed <- c(failed, paste0(basename(files[i]), ": ", why))
  }
  if (length(failed)) {
    showModal(modalDialog(
      title = if (length(failed) == 1) "Could not load this file" else "Could not load these files",
      tags$ul(lapply(failed, tags$li)),
      easyClose = TRUE
    ))
  }
  rest <- if (opened > 0) files[-seq_len(opened)] else character()
  if (length(rest)) {
    tryCatch(
      multi_add_uploads(data.frame(name = basename(rest), datapath = rest, stringsAsFactors = FALSE)),
      error = function(e) {
        showNotification(paste("Could not add the other files to Multi DE:", conditionMessage(e)), type = "error")
      }
    )
  }
}
# Nothing to load → edsState stays NULL → the landing screen (server-landing.R) shows.

figuresDataReactive <- reactiveValues(
  "volcanoStatic" = NULL,
  "volcanoMA" = NULL,
  "volcanoPlotly" = NULL,
  "heatmapDEFeatures" = NULL,
  "heatmapCustom" = NULL,
  "heatmapGO" = NULL,
  "heatmapRatios" = NULL,
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
