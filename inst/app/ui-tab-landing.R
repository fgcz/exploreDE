tabItem(
  tabName = "landingTab",
  fluidPage(
    fluidRow(
      column(
        width = 9,
        box(
          title = tagList(
            icon("upload"), "Load your own data",
            tags$span(
              "BETA",
              style = paste0(
                "margin-left:8px; padding:1px 7px; font-size:11px; font-weight:700; ",
                "letter-spacing:0.5px; color:#fff; background:#f0ad4e; border-radius:3px; ",
                "vertical-align:middle;"))),
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          helpText(
            "Welcome to exploreDE. You have opened the app without pointing it at a SUSHI
            dataset, so you can load your own data here instead.",
            br(), br(),
            "Upload a single ", tags$code(".h5ad"), " file built with exploreDE's
            ", tags$code("build_explore_h5ad()"), " constructor. It bundles your counts (and,
            optionally, sample metadata, DE results, and pathway/enrichment results) into the
            format the app reads. The file is parsed into memory for this session only — nothing
            is written to the server, and it is discarded when you disconnect.",
            br(), br(),
            "Only the tabs your data can support will appear; the rest stay hidden until the
            required pieces (counts, DE, pathways) are present."
          ),
          fileInput(
            inputId = "uploadH5ad",
            label = "Choose an exploreDE .h5ad file",
            accept = c(".h5ad", "application/octet-stream"),
            width = "85%"
          ),
          verbatimTextOutput("upload_h5ad_status")
        )
      )
    ),
    fluidRow(
      column(
        width = 9,
        box(
          title = tagList(icon("flask"), "Just exploring?"),
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          helpText(
            "No data of your own to hand? Load a demo dataset to see what exploreDE can do."
          ),
          actionButton(
            inputId = "loadDemo",
            label = "Load a demo dataset instead",
            icon = icon("play"),
            class = "btn-primary"
          )
        )
      )
    )
  )
)
