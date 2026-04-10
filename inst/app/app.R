## Sys.setenv(SHINYPROXY_USERNAME=Sys.getenv("USER"))

# JLR 2025 employee and project restrictions, username coming from app.R, projectFromUrl from above, be careful that ldap-utils in installed in bash and you can access LDAP (certificate installed?)
username <- Sys.getenv("SHINYPROXY_USERNAME")

cat("loading packages...\n\n")
packagesToLoad <- c(
  "shiny", "shinydashboard", "tidyverse", "ggpubr", "plotly", "DESeq2", "RColorBrewer",
  "ComplexHeatmap", "clusterProfiler", "DT", "colourpicker", "writexl", "circlize",
  "ezRun", "kableExtra", "ggrepel", "gplots", "sortable", "waiter", "ggprism", "ggbeeswarm",
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
if (file.exists("/srv/shiny-server/fgcz_apps/telemetry_module.R")) source("/srv/shiny-server/fgcz_apps/telemetry_module.R")

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
                               'multiHeatmapTab', 'multiKeggTab'];
              if (multiTabs.indexOf(event.value) !== -1) {
                $('body').addClass('multi-active');
              } else {
                $('body').removeClass('multi-active');
              }
            }
          });
        "))
      )
    ),
    use_waiter(),
    tabItems(
      source("ui-tab-summary.R", local = TRUE)$value,
      source("ui-tab-boxplot.R", local = TRUE)$value,
      source("ui-tab-heatmap.R", local = TRUE)$value,
      source("ui-tab-ordination.R", local = TRUE)$value,
      source("ui-tab-correlation.R", local = TRUE)$value,
      source("ui-tab-colour.R", local = TRUE)$value,
      source("ui-tab-DEGTable.R", local = TRUE)$value,
      source("ui-tab-volcano.R", local = TRUE)$value,
      source("ui-tab-hypergeo2.R", local = TRUE)$value,
      source("ui-tab-gsea2.R", local = TRUE)$value,
      source("ui-tab-goTile.R", local = TRUE)$value,
      source("ui-tab-kegg.R", local = TRUE)$value,
      source("ui-tab-replicability.R", local = TRUE)$value,
      source("ui-tab-multi-data-selection.R", local = TRUE)$value,
      source("ui-tab-multi-de-table.R", local = TRUE)$value,
      source("ui-tab-multi-boxplot.R", local = TRUE)$value,
      source("ui-tab-multi-upset.R", local = TRUE)$value,
      source("ui-tab-multi-tile.R", local = TRUE)$value,
      source("ui-tab-multi-bubble.R", local = TRUE)$value,
      source("ui-tab-multi-heatmap.R", local = TRUE)$value,
      source("ui-tab-multi-kegg.R", local = TRUE)$value
    ),
    tags$footer(
      style = "text-align: center; padding: 20px; color: #999; font-size: 16px; border-top: 1px solid #eee; margin-top: 20px;",
      HTML("Made with &hearts; in Switzerland 🇨🇭")
    )
  )
)

server <- function(input, output, session) {
  # Multi DE Reactive (defined early so it's available to server-initInputData.R)
  multiProjectInfo <- reactiveValues(projectDir = NULL)

  source("server-initInputData.R", local = TRUE)

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
  source("server-featureBucket.R", local = TRUE)
  source("server-summary.R", local = TRUE)
  source("server-DEGTable.R", local = TRUE)
  source("server-volcano2.R", local = TRUE)
  source("server-heatmap2.R", local = TRUE)
  source("server-ordination.R", local = TRUE)
  source("server-boxplot.R", local = TRUE)
  source("server-correlation.R", local = TRUE)
  source("server-hypergeo2.R", local = TRUE)
  source("server-gsea2.R", local = TRUE)
  source("server-goTile.R", local = TRUE)
  source("server-kegg.R", local = TRUE)
  source("server-replicability.R", local = TRUE)
  source("server-multi-data-selection.R", local = TRUE)
  source("server-multi-de-table.R", local = TRUE)
  source("server-multi-upset.R", local = TRUE)
  source("server-multi-bubble.R", local = TRUE)
  source("server-multi-heatmap.R", local = TRUE)
  source("server-multi-boxplot.R", local = TRUE)
  source("server-multi-tile.R", local = TRUE)
  source("server-multi-kegg.R", local = TRUE)

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
    req(inputDataReactive())
    dataType <- inputDataReactive()$dataType
    # Create dependency on menuUpdateTrigger for manual updates
    multiDeReactive$menuUpdateTrigger
    # Use isolate to check hasMultiData without automatic re-rendering
    hasData <- isolate(hasMultiData())

    # 0. Tabs
    menuItems <- list(
      menuItem(text = "Results Summary", tabName = "summaryTab", icon = icon("list")),
      menuItem(text = "Colours", tabName = "colourTab", icon = icon("palette")),

      # 1. Count tabs (Always shown)
      tags$li(class = "header", "Counts"),
      menuItem(text = "Ordination", tabName = "ordinationTab", icon = icon("meteor")),
      menuItem(text = "Heatmaps", tabName = "heatmapTab", icon = icon("map")),
      menuItem(text = "Boxplots", tabName = "boxplotTab", icon = icon("box-open")),
      menuItem(text = "Correlations", tabName = "correlationsTab", icon = icon("chart-line")),
      menuItem(text = "Replicability", tabName = "replicabilityTab", icon = icon("copy"))
    )

    # 2. DE tabs (Shown if NOT CountQC)
    if (!dataType %in% c("CountQC", "CrisprCountQC")) {
      menuItems <- c(menuItems, list(
        tags$li(class = "header", "DE Results"),
        menuItem(text = "DE Table", tabName = "degTableTab", icon = icon("table")),
        menuItem(text = "Volcano Plots", tabName = "volcanoTab", icon = icon("wifi"))
      ))
    }

    # 3. Pathway tabs (Shown ONLY if RNASeq)
    if (dataType == "RNASeq") {
      # Check if species is supported for KEGG (human, mouse, rat)
      keggSpeciesSupported <- tryCatch({
        spp <- ezRun::getSpecies(inputDataReactive()$param$refBuild)
        if (grepl("^Rattus", inputDataReactive()$param$refBuild)) spp <- "Rat"
        spp %in% c("Human", "Mouse", "Rat")
      }, error = function(e) FALSE)

      pathwayItems <- list(
        tags$li(class = "header", "Pathway Analyses"),
        menuItem(text = "Over-Representation Analysis", tabName = "oraTab", icon = icon("hubspot")),
        menuItem(text = "Gene Set Enrichment Analysis", tabName = "gseaTab", icon = icon("hubspot")),
        menuItem(text = "GO Tile Plots", tabName = "goTileTab", icon = icon("table"))
      )
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

          multiHeatmapItem <- menuItem(text = "Pathway Heatmaps", tabName = "multiHeatmapTab", icon = icon("fire"))
          multiHeatmapItem$attribs$class <- paste(multiHeatmapItem$attribs$class, "multi-menu-item")

          multiPathwayItems <- list(multiBubbleItem, multiHeatmapItem)

          if (keggSpeciesSupported) {
            multiKeggItem <- menuItem(text = "KEGG Pathway Plots", tabName = "multiKeggTab", icon = icon("egg"))
            multiKeggItem$attribs$class <- paste(multiKeggItem$attribs$class, "multi-menu-item")
            multiPathwayItems <- c(multiPathwayItems, list(multiKeggItem))
          }

          menuItems <- c(menuItems, multiPathwayItems)
        }
      }
    }

    # Use isolate(input$tabs) to preserve selection without creating a dependency loop
    currentTab <- isolate(input$tabs)
    sidebarMenu(id = "tabs", .list = menuItems, selected = currentTab)
  })
}

shinyApp(ui = ui, server = server)
