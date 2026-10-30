tabItem(
  tabName = "saveStateTab",
  fluidPage(
    fluidRow(
      column(
        width = 9,
        box(
          title = tagList(icon("floppy-disk"), "Save current session"),
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          helpText(
            "Download one file that captures your current settings across every tab —
            selectors, thresholds, colours, the feature bucket, the GO/KEGG and pathway
            selections, and the Multi DE datasets and their settings. The file stays on your
            machine; nothing is uploaded to the server. It is a .json file, so it can be opened
            in any text editor. Re-upload it later (below) to pick up exactly where you left off.",
            br(),
            br(),
            tags$b("Note: "),
            "a saved session can only be restored on the exact same dataset it was made from —
            open the app with the original data link before restoring, or the upload will be
            refused. Multi DE datasets are re-read on restore, so ones that have since moved,
            or that you can no longer access, are skipped."
          ),
          downloadButton(outputId = "download_state_json", label = "Save session (.json)", class = "btn-primary"),
          br(),
          br(),
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
            "Upload a .json file saved from a previous session to restore the same settings
            across every tab. The bundle must have been saved from this exact dataset — one from a
            different dataset is refused. Restore is applied progressively as tabs mount — open a
            tab to see its settings snap into place."
          ),
          fileInput(
            inputId = "load_state_file",
            label = "Choose a saved session file (.json)",
            accept = c(".json", "application/json"),
            width = "85%"
          ),
          verbatimTextOutput("load_state_status")
        )
      )
    )
  )
)
