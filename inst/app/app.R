## Sys.setenv(SHINYPROXY_USERNAME=Sys.getenv("USER"))

# JLR 2025 employee and project restrictions, username coming from app.R, projectFromUrl from above, be careful that ldap-utils in installed in bash and you can access LDAP (certificate installed?)
username <- Sys.getenv("SHINYPROXY_USERNAME")

cat("loading packages...\n\n")
packagesToLoad <- c(
  "shiny", "shinydashboard", "tidyverse", "ggpubr", "plotly", "DESeq2", "RColorBrewer",
  "ComplexHeatmap", "clusterProfiler", "DT", "colourpicker", "writexl", "circlize",
  "kableExtra", "ggrepel", "gplots", "sortable", "waiter", "ggprism", "ggbeeswarm",
  "rstatix", "gridExtra", "shinylogs", "parallel", "plyr", "reshape2", "shinycssloaders", "shinyjs", "GGally", "patchwork",
  "Matrix", "matrixStats", "SingleCellExperiment", "fresh", "pathview", "openxlsx", "dendextend", "officer"
)
invisible(lapply(packagesToLoad, function(pkg) {
  suppressPackageStartupMessages(suppressWarnings(library(pkg, character.only = TRUE, quietly = TRUE)))
}))
cat("... packages loaded!")
# Ensure shiny::validate is not masked by jsonlite::validate
validate <- shiny::validate
reactiveConsole(TRUE)
# Allow larger uploads for user-supplied .h5ad datasets (BYO data, M2).
options(shiny.maxRequestSize = 1000 * 1024^2)
if (file.exists("/srv/shiny-server/fgcz_apps/telemetry_module.R")) source("/srv/shiny-server/fgcz_apps/telemetry_module.R")

# Source package R/ helpers at startup (dev mode — no install required).
# Provides global definitions of the conversion / plotting helpers.
# server-*.R files source() with local = TRUE to redefine reactive-coupled
# helpers (bucket_*, resolve_features, etc.) inside the server scope.
r_dir <- if (dir.exists("../../R")) "../../R" else if (dir.exists("R")) "R" else NULL
if (!is.null(r_dir)) {
  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  invisible(lapply(r_files, source))
  cat("sourced", length(r_files), "R/ helper file(s) from", r_dir, "\n")
}

my_theme <- create_theme(
  adminlte_color(
    light_blue = "#527085"
  )
)

ui <- dashboardPage(
  dashboardHeader(
    title = "Explore DE",
    tags$li(
      title = "Please include the URL to the dataset in your email.",
      a(
        href = "mailto:sequencing@fgcz.ethz.ch?subject=exploreDE-shiny-app-feedback&body=Please%20include%20the%20URL%20you%20are%20having%20issues%20with.%20Thanks!%0A%0A",
        "Request Features/Report Bugs"
      ),
      class = "dropdown"
    ),
    tags$li(
      a(
        href = "http://www.fgcz.ch",
        target = "_blank",
        img(src = "fgcz_logo.png", title = "FGCZ", height = "30px"),
        style = "padding-top:10px; padding-bottom:5px;"
      ),
      class = "dropdown"
    ),
    tags$li(
      a(
        href = "http://www.ethz.ch/en.html",
        target = "_blank",
        img(src = "eth_logo.png", title = "FGCZ", height = "22px"),
        style = "padding-top:13px; padding-bottom:10px;"
      ),
      class = "dropdown"
    ),
    tags$li(
      a(
        href = "http://www.uzh.ch/en.html",
        target = "_blank",
        img(src = "University_of_Zurich_Logo.png", title = "FGCZ", height = "30px"),
        style = "padding-top:10px; padding-bottom:5px;"
      ),
      class = "dropdown"
    )
  ),
  dashboardSidebar(
    shinyjs::useShinyjs(),
    sidebarMenuOutput("sidebarMenu")
  ),
  dashboardBody(
    if (exists("use_telemetry_js")) use_telemetry_js(),
    use_theme(my_theme),
    tags$head(
      tags$link(rel = "shortcut icon", href = "sushi.png"),
      tags$head(
        tags$style(
          HTML("
            .shiny-split-layout > div {
              overflow: visible;
            }
            .box {
              box-shadow: 0 4px 8px rgba(0,0,0,0.1);
              transition: all 0.3s ease;
            }
            .box:hover {
              box-shadow: 0 8px 16px rgba(0,0,0,0.2);
            }
            .box.box-solid.box-primary>.box-header {
              color:#fff;
              background: linear-gradient(135deg, #527085 0%, #8FA9B8 100%);
            }
            .box.box-solid.box-primary{
              border-bottom-color:#527085;
              border-left-color:#527085;
              border-right-color:#527085;
              border-top-color:#527085;
            }
            .skin-blue .box.box-solid.box-success>.box-header {
              color:#fff !important;
              background: #2d7a3e !important;
              background: linear-gradient(135deg, #2d7a3e 0%, #3a9950 100%) !important;
            }
            .skin-blue .box.box-solid.box-success{
              border-bottom-color:#2d7a3e !important;
              border-left-color:#2d7a3e !important;
              border-right-color:#2d7a3e !important;
              border-top-color:#2d7a3e !important;
            }
            .content-wrapper .box.box-solid.box-success>.box-header {
              color:#fff !important;
              background: #2d7a3e !important;
              background: linear-gradient(135deg, #2d7a3e 0%, #3a9950 100%) !important;
            }
            .content-wrapper .box.box-solid.box-success{
              border-bottom-color:#2d7a3e !important;
              border-left-color:#2d7a3e !important;
              border-right-color:#2d7a3e !important;
              border-top-color:#2d7a3e !important;
            }
            .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
              background-color: #6E8D9E;
              transition: background-color 0.3s ease;
            }
            .skin-blue .sidebar-menu > li:hover > a {
              border-left-color: #3A5060;
            }
            /* Main Header Gradient */
            .skin-blue .main-header .logo {
              background: linear-gradient(135deg, #527085 0%, #8FA9B8 100%) !important;
            }
            .skin-blue .main-header .navbar {
              background: linear-gradient(135deg, #527085 0%, #8FA9B8 100%) !important;
            }
            /* body */
            .content-wrapper, .right-side {
              background-color: #FFFFFF;
            }
            .sidebar-menu > li.header {
              font-size: 16px;
              font-weight: bold;
              color: #b8c7ce;
              padding: 10px 25px 10px 15px;
              text-shadow: 1px 1px 2px rgba(0,0,0,0.1);
            }
            /* FORCE GREEN BOXES - Must be last! */
            .box.box-success > .box-header,
            .box.box-solid.box-success > .box-header {
              background: linear-gradient(135deg, #2d7a3e 0%, #3a9950 100%) !important;
            }
            .box.box-success,
            .box.box-solid.box-success {
              border-color: #2d7a3e !important;
              border-top-color: #2d7a3e !important;
            }
            /* Green hover for multi-menu items only */
            .skin-blue .sidebar-menu li.multi-menu-item:hover > a {
              background-color: #2d7a3e !important;
              border-left-color: #2d7a3e !important;
            }
            .skin-blue .sidebar-menu li.multi-menu-item.active > a {
              background-color: #2d7a3e !important;
              border-left-color: #2d7a3e !important;
            }
            /* Green header when on multi tabs */
            body.multi-active.skin-blue .main-header .logo {
              background: linear-gradient(135deg, #2d7a3e 0%, #3a9950 100%) !important;
            }
            body.multi-active.skin-blue .main-header .navbar {
              background: linear-gradient(135deg, #2d7a3e 0%, #3a9950 100%) !important;
            }
            ")
        ),
        tags$script(HTML("
          $(document).on('shiny:inputchanged', function(event) {
            if (event.name === 'tabs') {
              var multiTabs = ['multiDataSelectionTab', 'multiDETableTab', 'multiBoxplotTab',
                               'multiUpsetPlotTab', 'multiTilePlotTab', 'multiBubbleTab',
                               'multiPathwaySelectTab', 'multiHeatmapTab', 'multiKeggTab'];
              if (multiTabs.indexOf(event.value) !== -1) {
                $('body').addClass('multi-active');
              } else {
                $('body').removeClass('multi-active');
              }
            }
          });

          // Hidden easter egg: triple-click the footer heart (#bananaHeart)
          // within 1.2s to toggle banana mode (server: input$bananaToggle).
          (function () {
            var n = 0, t = null;
            $(document).on('click', '#bananaHeart', function () {
              clearTimeout(t);
              if (++n >= 3) { n = 0; Shiny.setInputValue('bananaToggle', Date.now()); }
              t = setTimeout(function () { n = 0; }, 1200);
            });
            console.log('%c🍌  …there is a banana in here somewhere', 'font-size:13px');
          })();
        "))
      )
    ),
    use_waiter(),
    # Hide the full-page loading spinner the instant the summary page's primary
    # table (settingsTable) actually renders on the client. This is a UI helper:
    # it injects a document-level `shiny:value` listener that calls
    # waiter.hide() when that output paints — event-driven, no fixed delay.
    # (Must live in the UI; called server-side it only builds an unused tag.)
    waiter::waiter_hide_on_render("settingsTable"),
    tabItems(
      source("ui-tab-landing.R", local = TRUE)$value,
      source("ui-tab-summary.R", local = TRUE)$value,
      source("ui-tab-boxplot.R", local = TRUE)$value,
      source("ui-tab-heatmap.R", local = TRUE)$value,
      source("ui-tab-ordination.R", local = TRUE)$value,
      source("ui-tab-correlation.R", local = TRUE)$value,
      source("ui-tab-colour.R", local = TRUE)$value,
      source("ui-tab-DEGTable.R", local = TRUE)$value,
      source("ui-tab-volcano.R", local = TRUE)$value,
      source("ui-tab-deQc.R", local = TRUE)$value,
      source("ui-tab-hypergeo2.R", local = TRUE)$value,
      source("ui-tab-gsea2.R", local = TRUE)$value,
      source("ui-tab-goTile.R", local = TRUE)$value,
      source("ui-tab-kegg.R", local = TRUE)$value,
      source("ui-tab-replicability.R", local = TRUE)$value,
      source("ui-tab-countQc.R", local = TRUE)$value,
      source("ui-tab-multi-data-selection.R", local = TRUE)$value,
      source("ui-tab-multi-de-table.R", local = TRUE)$value,
      source("ui-tab-multi-boxplot.R", local = TRUE)$value,
      source("ui-tab-multi-upset.R", local = TRUE)$value,
      source("ui-tab-multi-tile.R", local = TRUE)$value,
      source("ui-tab-multi-bubble.R", local = TRUE)$value,
      source("ui-tab-multi-pathway-select.R", local = TRUE)$value,
      source("ui-tab-multi-heatmap.R", local = TRUE)$value,
      source("ui-tab-multi-kegg.R", local = TRUE)$value,
      source("ui-tab-saveState.R", local = TRUE)$value
    ),
    tags$footer(
      style = "text-align: center; padding: 20px; color: #999; font-size: 16px; border-top: 1px solid #eee; margin-top: 20px;",
      HTML("Made with "),
      # Bigger click target via padding cancelled by an equal negative margin,
      # so the footer text doesn't shift; no :hover styling, so nothing moves.
      tags$span(
        id = "bananaHeart",
        style = paste(
          "display:inline-block;",
          "padding:.35em .7em;",
          "margin:-.35em -.7em;",
          "cursor:pointer;",
          "user-select:none; -webkit-user-select:none;"
        ),
        HTML("&hearts;")
      ),
      HTML(" in Switzerland 🇨🇭")
    )
  )
)

server <- function(input, output, session) {
  # ── Session save/restore primitives ────────────────────────────────────────
  # `restoring` is TRUE only while a saved-state bundle is being applied (set in
  # server-saveState.R). `restore_seed` holds the saved input snapshot so widgets
  # rebuilt mid-restore can seed their saved value. `restore_default()` is the
  # single convention: during a restore it returns the saved value for `id`
  # (optionally intersected with `valid` = the widget's current choices), else the
  # computed default — a no-op outside restore, so live behaviour is unchanged.
  # isolate() so reading the store never adds a reactive dependency to the caller.
  restoring    <- reactiveVal(FALSE)
  restore_seed <- reactiveValues()
  # Thin closure over this session's reactives; logic lives in the exported,
  # unit-tested restore_default_value() (R/utils_savestate.R).
  restore_default <- function(id, computed_default, valid = NULL) {
    restore_default_value(id, computed_default, valid, restoring, restore_seed)
  }

  # ── Banana mode (hidden easter egg) ────────────────────────────────────────
  # Swaps p-value significance markers for a pile of banana images across
  # boxplots and tile/GO-tile heatmaps. Triple-click the footer heart (JS in
  # dashboardBody) or append ?banana=1. Defined here so it's in scope for every
  # lazily-source()d server-*.R module. Rendering logic: R/banana_mode.R.
  bananaMode <- reactiveVal(FALSE)
  observeEvent(input$bananaToggle, {
    on <- !bananaMode(); bananaMode(on)
    showNotification(
      if (on) "\U0001F34C Banana mode engaged" else "Banana mode off",
      duration = 2, type = "message"
    )
  }, ignoreInit = TRUE)
  observe({
    q <- parseQueryString(session$clientData$url_search)
    if (identical(q[["banana"]], "1")) bananaMode(TRUE)
  })

  # Multi DE Reactive (defined early so it's available to server-initInputData.R)
  multiProjectInfo <- reactiveValues(projectDir = NULL)

  source("server-initInputData.R", local = TRUE)
  source("server-landing.R", local = TRUE)

  # Telemetry (skipped when DB directory is unavailable, e.g. local dev)
  telemetry_app_name <- "exploreDE_genomics"
  if (exists("dataType") && dataType == "proteomics") {
    telemetry_app_name <- "exploreDE_proteomics"
  }
  telemetry_db_dir <- "/srv/shiny-server/fgcz_apps/telemetry_reader/apps_data"
  telemetry_ready <- dir.exists(telemetry_db_dir) ||
    dir.create(telemetry_db_dir, recursive = TRUE, showWarnings = FALSE)
  if (exists("setup_telemetry_global") && telemetry_ready) {
    telemetry <- setup_telemetry_global(telemetry_app_name, data_dir = "/srv/shiny-server/fgcz_apps/telemetry_reader/apps_data")
    start_telemetry_server(input, session, telemetry)
  }
  # Source download helpers (adjust path as needed for deployment structure)
  if (file.exists("../../R/downloadHelpers.R")) {
    source("../../R/downloadHelpers.R")
  } else if (file.exists("R/downloadHelpers.R")) {
     source("R/downloadHelpers.R")
  }
  # Data-dependent server modules run `eds <- inputDataReactive()$.eds; req(eds)`
  # at source time, so they can only be sourced once a dataset exists (otherwise
  # req() aborts the whole server function and the app never starts). Defer them
  # to the first moment edsState is populated — a ?data= URL fires this on the
  # first flush; an upload/demo fires it later. Sourced ONCE, all into the same
  # handler environment so cross-module references resolve, and into a scope that
  # still sees the server's output/input/session and the reactives below.
  .dataModulesSourced <- FALSE
  observeEvent(edsState(), {
    if (.dataModulesSourced || is.null(edsState())) return()
    .dataModulesSourced <<- TRUE
    # Keep the loading spinner up across module sourcing (the landing observers
    # show it on click; a ?data= URL has no click, so ensure it's shown here).
    waiter::waiter_show(
      html  = tagList(waiter::spin_fading_circles(),
                      tags$h4("Preparing your dataset…", style = "color:#fff; margin-top:1em;")),
      color = "#527085"
    )
    # NB: the spinner is hidden client-side by the waiter_hide_on_render() helper
    # wired into the UI (see dashboardBody) — it fires when settingsTable paints,
    # which is exactly "summary is ready". Nothing to hide here: hiding on exit
    # would drop the overlay the instant module-sourcing finishes (several
    # flushes before any output renders), flashing the landing screen through.
    source("server-featureBucket.R", local = TRUE)
    source("server-summary.R", local = TRUE)
    source("server-DEGTable.R", local = TRUE)
    source("server-volcano2.R", local = TRUE)
    source("server-deQc.R", local = TRUE)
    source("server-heatmap2.R", local = TRUE)
    source("server-ordination.R", local = TRUE)
    source("server-boxplot.R", local = TRUE)
    source("server-correlation.R", local = TRUE)
    source("server-hypergeo2.R", local = TRUE)
    source("server-gsea2.R", local = TRUE)
    source("server-goTile.R", local = TRUE)
    source("server-kegg.R", local = TRUE)
    source("server-replicability.R", local = TRUE)
    source("server-countQc.R", local = TRUE)
    source("server-multi-data-selection.R", local = TRUE)
    source("server-multi-de-table.R", local = TRUE)
    source("server-multi-upset.R", local = TRUE)
    source("server-multi-bubble.R", local = TRUE)
    source("server-multi-pathway-select.R", local = TRUE)
    source("server-multi-heatmap.R", local = TRUE)
    source("server-multi-boxplot.R", local = TRUE)
    source("server-multi-tile.R", local = TRUE)
    source("server-multi-kegg.R", local = TRUE)
    source("server-saveState.R", local = TRUE)
  }, ignoreNULL = TRUE)

  # Multi DE Reactive
  multiDeReactive <- reactiveValues(
    results = list(),
    available = NULL,
    menuUpdateTrigger = 0  # Add trigger to reactive values so it's accessible everywhere
  )

  # Reactive to track if we have multi-DE results
  hasMultiData <- reactive({
    length(multiDeReactive$results) > 0
  })

  output$sidebarMenu <- renderMenu({
    # No data yet (opened without ?data=): show only the landing/upload entry.
    if (is.null(edsState())) {
      return(sidebarMenu(
        id = "tabs",
        menuItem(text = "Load data", tabName = "landingTab", icon = icon("upload"))
      ))
    }
    eds      <- inputDataReactive()$.eds
    dataType <- eds$data_type()
    # Column-presence check: cheap, no DB lookup. goTile tab is only shown when
    # SUSHI has annotated features with GO terms (GO BP / MF / CC columns).
    has_go   <- any(c("GO BP", "GO MF", "GO CC") %in% colnames(eds$features()))
    # Create dependency on menuUpdateTrigger for manual updates
    multiDeReactive$menuUpdateTrigger
    # Use isolate to check hasMultiData without automatic re-rendering
    hasData <- isolate(hasMultiData())

    # 0. Tabs
    menuItems <- list(

      # 0. Session (always shown, all data types)
      tags$li(class = "header", "Session"),
      menuItem(text = "Results Summary", tabName = "summaryTab", icon = icon("list")),
      menuItem(text = "Colours / Labels", tabName = "colourTab", icon = icon("palette")),
      menuItem(text = "Save / Load Session (beta)", tabName = "saveStateTab", icon = icon("floppy-disk")),

      # 1. Count tabs (Always shown)
      tags$li(class = "header", "Counts"),
      menuItem(text = "Ordination", tabName = "ordinationTab", icon = icon("meteor")),
      menuItem(text = "Heatmaps", tabName = "heatmapTab", icon = icon("map")),
      menuItem(text = "Boxplots", tabName = "boxplotTab", icon = icon("box-open")),
      menuItem(text = "Correlations", tabName = "correlationsTab", icon = icon("chart-line")),
      menuItem(text = "Replicability", tabName = "replicabilityTab", icon = icon("copy"))
    )
    # Count QC — RNASeq / CountQC / CrisprCountQC only (hidden for proteomics)
    if (dataType %in% c("RNASeq", "CountQC", "CrisprCountQC")) {
      menuItems <- c(menuItems, list(
        menuItem(text = "Count QC", tabName = "countQcTab", icon = icon("chart-bar"))
      ))
    }

    # 2. DE tabs (Shown if NOT CountQC, and only when DE results exist —
    #    a counts-only BYO upload carries no contrasts)
    if (!dataType %in% c("CountQC", "CrisprCountQC") && eds$has_de()) {
      deMenuItems <- list(
        tags$li(class = "header", "DE Results"),
        menuItem(text = "DE Table", tabName = "degTableTab", icon = icon("table")),
        menuItem(text = "Volcano Plots", tabName = "volcanoTab", icon = icon("wifi")),
        menuItem(text = "DE QC", tabName = "deQcTab", icon = icon("chart-line"))
      )
      menuItems <- c(menuItems, deMenuItems)
    }

    # 3. Pathway tabs (Shown ONLY if RNASeq with DE results)
    if (dataType == "RNASeq" && eds$has_de()) {
      # Check if species is supported for KEGG (human, mouse, rat)
      keggSpeciesSupported <- tryCatch({
        !is.null(kegg_species_label(eds$ref_build()))
      }, error = function(e) FALSE)

      pathwayItems <- list(
        tags$li(class = "header", "Pathway Analyses"),
        menuItem(text = "Over-Representation Analysis", tabName = "oraTab", icon = icon("hubspot")),
        menuItem(text = "Gene Set Enrichment Analysis", tabName = "gseaTab", icon = icon("hubspot"))
      )
      if (has_go) {
        pathwayItems <- c(pathwayItems, list(
          menuItem(text = "Pathway Tile Plots", tabName = "goTileTab", icon = icon("table"))
        ))
      }
      if (keggSpeciesSupported) {
        pathwayItems <- c(pathwayItems, list(
          menuItem(text = "KEGG Pathway Plots", tabName = "keggTab", icon = icon("egg"))
        ))
      }
      menuItems <- c(menuItems, pathwayItems)
    }

    # 4. Multi DE Analysis (NOT shown for CountQC)
    if (!dataType %in% c("CountQC", "CrisprCountQC")) {
      menuItems <- c(menuItems, list(
        tags$li(class = "header", "Multi DE Analysis")
      ))

      # Add multi menu items with custom class attribute
      multiDataItem <- menuItem(text = "Data Selection", tabName = "multiDataSelectionTab", icon = icon("database"))
      multiDataItem$attribs$class <- paste(multiDataItem$attribs$class, "multi-menu-item")
      menuItems <- c(menuItems, list(multiDataItem))

      if (hasData) {
        # Generic tabs (always shown if data exists)
        multiDETableItem <- menuItem(text = "Integrated DE Table", tabName = "multiDETableTab", icon = icon("table"))
        multiDETableItem$attribs$class <- paste(multiDETableItem$attribs$class, "multi-menu-item")

        multiBoxplotItem <- menuItem(text = "Multi-Contrast Boxplots", tabName = "multiBoxplotTab", icon = icon("chart-bar"))
        multiBoxplotItem$attribs$class <- paste(multiBoxplotItem$attribs$class, "multi-menu-item")

        multiUpsetItem <- menuItem(text = "Upset/Venn Plots", tabName = "multiUpsetPlotTab", icon = icon("ellipsis-h"))
        multiUpsetItem$attribs$class <- paste(multiUpsetItem$attribs$class, "multi-menu-item")

        multiTileItem <- menuItem(text = "DE Tile Plots", tabName = "multiTilePlotTab", icon = icon("th"))
        multiTileItem$attribs$class <- paste(multiTileItem$attribs$class, "multi-menu-item")

        menuItems <- c(menuItems, list(multiDETableItem, multiBoxplotItem, multiUpsetItem, multiTileItem))

        # Pathway tabs (only for RNASeq)
        if (dataType == "RNASeq") {
          multiBubbleItem <- menuItem(text = "Pathway Bubble Plots", tabName = "multiBubbleTab", icon = icon("circle"))
          multiBubbleItem$attribs$class <- paste(multiBubbleItem$attribs$class, "multi-menu-item")

          multiPathwaySelectItem <- menuItem(text = "Pathway Selection", tabName = "multiPathwaySelectTab", icon = icon("diagram-project"))
          multiPathwaySelectItem$attribs$class <- paste(multiPathwaySelectItem$attribs$class, "multi-menu-item")

          multiHeatmapItem <- menuItem(text = "Pathway Heatmaps", tabName = "multiHeatmapTab", icon = icon("fire"))
          multiHeatmapItem$attribs$class <- paste(multiHeatmapItem$attribs$class, "multi-menu-item")

          multiPathwayItems <- list(multiBubbleItem, multiPathwaySelectItem, multiHeatmapItem)

          if (keggSpeciesSupported) {
            multiKeggItem <- menuItem(text = "KEGG Pathway Plots", tabName = "multiKeggTab", icon = icon("egg"))
            multiKeggItem$attribs$class <- paste(multiKeggItem$attribs$class, "multi-menu-item")
            multiPathwayItems <- c(multiPathwayItems, list(multiKeggItem))
          }

          menuItems <- c(menuItems, multiPathwayItems)
        }
      }
    }

    # Use isolate(input$tabs) to preserve selection without creating a dependency loop.
    # Arriving from the landing screen (or first load): land on the summary tab.
    currentTab <- isolate(input$tabs)
    if (is.null(currentTab) || currentTab == "landingTab") currentTab <- "summaryTab"
    sidebarMenu(id = "tabs", .list = menuItems, selected = currentTab)
  })
}

shinyApp(ui = ui, server = server)
