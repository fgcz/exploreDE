tabItem(
  tabName = "landingTab",
  fluidPage(
    fluidRow(
      column(
        width = 9,
        box(
          title = tagList(
            icon("upload"),
            "Load your own data",
            tags$span(
              "BETA",
              style = paste0(
                "margin-left:8px; padding:1px 7px; font-size:11px; font-weight:700; ",
                "letter-spacing:0.5px; color:#fff; background:#f0ad4e; border-radius:3px; ",
                "vertical-align:middle;"
              )
            )
          ),
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          helpText(
            "Welcome to exploreDE. You have opened the app without pointing it at a SUSHI
            dataset, so you can load your own data here instead.",
            br(),
            br(),
            "Upload one or more ",
            tags$code(".h5ad"),
            " files built with exploreDE's
            ",
            tags$code("build_explore_h5ad()"),
            " constructor. A file bundles your counts (and,
            optionally, sample metadata, DE results, and pathway/enrichment results) into the
            format the app reads. Files are parsed into memory for this session only — nothing
            is written to the server, and they are discarded when you disconnect. Up to 1 GB per
            upload.",
            br(),
            br(),
            "Only the tabs your data can support will appear; the rest stay hidden until the
            required pieces (counts, DE, pathways) are present.",
            br(),
            br(),
            "With several files, the first opens in the tabs and the contrasts of all of them
            are loaded into Multi DE Analysis, where the Data Selection tab can add more files
            later."
          ),
          fileInput(
            inputId = "uploadH5ad",
            label = "Choose exploreDE .h5ad file(s)",
            accept = c(".h5ad", "application/octet-stream"),
            multiple = TRUE,
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
