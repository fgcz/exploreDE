tabItem(
  tabName = "saveStateTab",
  fluidPage(
    fluidRow(
      column(
        width = 9,
        box(
          title = tagList(
            icon("floppy-disk"), "Save current session",
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
            "Download a single .rds bundle that captures your current settings across every tab —
            selectors, thresholds, colours, the feature bucket, and the GO/KEGG and pathway
            selections. The file stays on your machine; nothing is uploaded to the server.
            Re-upload it later (below) to pick up exactly where you left off.",
            br(), br(),
            tags$b("Note: "),
            "a saved session can only be restored on the exact same dataset it was made from —
            open the app with the original data link before restoring, or the upload will be
            refused. Multi-contrast (Multi DE) selections are not included in this version."
          ),
          downloadButton(outputId = "download_state", label = "Save session (.rds)", class = "btn-primary"),
          br(), br(),
          verbatimTextOutput("save_state_status")
        )
      )
    ),
    fluidRow(
      column(
        width = 9,
        box(
          title = tagList(icon("upload"), "Load a saved session"),
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          helpText(
            "Upload an .rds bundle exported from a previous session to restore the same settings
            across every tab. The bundle must have been saved from this exact dataset — one from a
            different dataset is refused. Restore is applied progressively as tabs mount — open a
            tab to see its settings snap into place."
          ),
          fileInput(
            inputId = "load_state_file",
            label = "Choose a saved session file (.rds)",
            accept = c(".rds", "application/octet-stream"),
            width = "85%"
          ),
          verbatimTextOutput("load_state_status")
        )
      )
    )
  )
)
