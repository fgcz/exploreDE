## Sys.setenv(SHINYPROXY_USERNAME=Sys.getenv("USER"))

# JLR 2025 employee and project restrictions, username coming from app.R, projectFromUrl from above, be careful that ldap-utils in installed in bash and you can access LDAP (certificate installed?)
username <- Sys.getenv("SHINYPROXY_USERNAME")

cat("loading packages...\n\n")
packagesToLoad <- c(
  "shiny", "shinydashboard", "tidyverse", "ggpubr", "plotly", "DESeq2", "RColorBrewer",
  "ComplexHeatmap", "clusterProfiler", "DT", "colourpicker", "writexl", "circlize",
  "ezRun", "kableExtra", "ggrepel", "gplots", "sortable", "waiter", "ggprism", "ggbeeswarm",
  "rstatix", "gridExtra", "shinylogs", "parallel", "plyr", "shinycssloaders", "GGally", "patchwork",
  "Matrix", "SingleCellExperiment", "fresh", "pathview", "exploreDE"
)
invisible(lapply(packagesToLoad, function(pkg) {
  suppressPackageStartupMessages(suppressWarnings(library(pkg, character.only = TRUE, quietly = TRUE)))
}))
cat("... packages loaded!")
reactiveConsole(TRUE)

my_theme = create_theme(
  adminlte_color(
    light_blue = "#86a5bf"
  )
)

ui = dashboardPage(
  dashboardHeader(
    title = "Explore DE",
    tags$li(
      title = "Please include the URL to the dataset in your email.",
      a(
        href = 'mailto:sequencing@fgcz.ethz.ch?subject=exploreDE-shiny-app-feedback&body=Please%20include%20the%20URL%20you%20are%20having%20issues%20with.%20Thanks!%0A%0A', 
        "Request Features/Report Bugs"), 
      class = "dropdown"
    ),
    tags$li(
      a(href = 'http://www.fgcz.ch', 
        target = "_blank",
        img(src = 'fgcz_logo.png', title = "FGCZ", height = "30px"),
        style = "padding-top:10px; padding-bottom:5px;"),
      class = "dropdown"),
    tags$li(
      a(href = 'http://www.ethz.ch/en.html',
        target = "_blank",
        img(src = 'eth_logo.png', title = "FGCZ", height = "22px"),
        style = "padding-top:13px; padding-bottom:10px;"),
      class = "dropdown"),
    tags$li(
      a(href = 'http://www.uzh.ch/en.html',
        target = "_blank",
        img(src = 'University_of_Zurich_Logo.png', title = "FGCZ", height = "30px"),
        style = "padding-top:10px; padding-bottom:5px;"),
      class = "dropdown")
  ),
  dashboardSidebar(
    shinyjs::useShinyjs(),
    sidebarMenu(
      id = "tabs",
      menuItem(text = "Results Summary", tabName = "summaryTab",  icon = icon("list")),
      menuItem(text = "DE Table", tabName = "degTableTab", icon = icon("table")),
      menuItem(text = "Volcano Plots", tabName = "volcanoTab", icon = icon("wifi")),
      menuItem(text = "Heatmaps", tabName = "heatmapTab", icon = icon("map")),
      menuItem(text = "PCA Plots", tabName = "pcaPlotTab", icon = icon("meteor")),
      menuItem(text = "Boxplots", tabName = "boxplotTab", icon = icon("box-open")),
      menuItem(text = "Correlations", tabName = "correlationsTab", icon = icon("chart-line")),
      conditionalPanel(
        condition = "output.test=='RNASeq'",
        sidebarMenu(
          menuItem(text = "Over-Representation Analysis", tabName = "oraTab", icon = icon("hubspot")),
          menuItem(text = "Gene Set Enrichment Analysis", tabName = "gseaTab", icon = icon("hubspot")),
          menuItem(text = "GO Tile Plots", tabName = "goTileTab", icon = icon("table")),
          menuItem(text = "KEGG Pathway Plots", tabName = "keggTab", icon = icon("egg"))
        )
      ),
      menuItem(text= "Colours", tabName = "colourTab", icon = icon("palette"))
    )
  ), 
  dashboardBody(
    use_theme(my_theme),
    tags$head(
      tags$link(rel = "shortcut icon", href = "sushi.png"),
      tags$head(
        tags$style(
          HTML("
            .shiny-split-layout > div {
              overflow: visible;
            }
            .box.box-solid.box-primary>.box-header {
            color:#fff; background:#86a5bf}
            .box.box-solid.box-primary{
            border-bottom-color:#86a5bf;
            border-left-color:#86a5bf;
            border-right-color:#86a5bf;
            border-top-color:#86a5bf;
            }
            .box.box-solid.box-success>.box-header {
            color:#fff; background:#99adad}
            .box.box-solid.box-success{
            border-bottom-color:#99adad;
            border-left-color:#99adad;
            border-right-color:#99adad;
            border-top-color:#99adad;
            }
            .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
            background-color: #93bac2;
            }
            .skin-blue .sidebar-menu > li:hover > a {
              border-left-color: #455b73;
            }
            /* body */
            .content-wrapper, .right-side {
            background-color: #FFFFFF;
            }
            "
            )
          )
        )
      ),
    use_waiter(),
    tabItems(
      source("ui-tab-summary.R", local = TRUE)$value,
      source("ui-tab-DEGTable.R", local = TRUE)$value,
      source("ui-tab-volcano.R", local = TRUE)$value,
      source("ui-tab-heatmap.R", local = TRUE)$value,
      source("ui-tab-pca.R", local = TRUE)$value,
      source("ui-tab-boxplot.R", local = TRUE)$value,
      source("ui-tab-correlation.R", local = TRUE)$value,
      source("ui-tab-hypergeo2.R", local = TRUE)$value,
      source("ui-tab-gsea2.R", local = TRUE)$value,
      source("ui-tab-goTile.R", local = TRUE)$value,
      source("ui-tab-kegg.R", local = TRUE)$value,
      source("ui-tab-colour.R", local = TRUE)$value
    )
  )
)

server = function(input, output, session) {
  source("server-initInputData.R", local = TRUE)
  source("server-summary.R", local = TRUE)
  source("server-DEGTable.R", local = TRUE)
  source("server-volcano2.R", local = TRUE)
  source("server-heatmap2.R", local = TRUE)
  source("server-pca.R", local = TRUE)
  source("server-boxplot.R", local = TRUE)
  source("server-correlation.R", local = TRUE)
  source("server-hypergeo2.R", local = TRUE)
  source("server-gsea2.R", local = TRUE)
  source("server-goTile.R", local = TRUE)
  source("server-kegg.R", local = TRUE)
}

shinyApp(ui = ui, server = server)
