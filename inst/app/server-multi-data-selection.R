# Server logic for Multi DE Data Selection
# Note: multiProjectInfo is defined in app.R

# Where multi-DE datasets may come from (as for ?data=). Every link, typed or
# ticked, goes through multi_resolve_dataset() (R/utils_multi_de.R) and then
# the same project check as ?data=: Course data needs none, gStore projects
# need the user's LDAP role or project membership.
MULTI_COURSE_ROOT <- "/srv/GT/analysis/course_sushi/public/gstore/projects"
MULTI_DATA_ROOTS  <- c("/srv/gstore/projects", MULTI_COURSE_ROOT)
multi_access_ok <- function(resolved) {
    if (identical(resolved$root, MULTI_COURSE_ROOT)) return(TRUE)
    ldap_allows_project(session_ldap(), resolved$project)
}

# Conditional UI - Only show data selection box for genomics (not proteomics)
output$multiDataSelectionBoxUI <- renderUI({
    # Check if we have proteomics data loaded
    .eds <- inputDataReactive()$.eds
    # Proteomics fans a single SE out into per-contrast objects; an uploaded
    # .h5ad is a single self-contained dataset whose contrasts are auto-loaded
    # below (section 5). In both cases the SUSHI project-browser / URL-import
    # box is meaningless, so hide it and let the contrasts populate directly.
    if (!is.null(.eds) && .eds$omics_class() == "proteomics") {
        return(NULL)
    }
    if (!is.null(.eds) && .eds$backend_kind() == "h5ad") {
        # Multi-contrast DE needs at least two contrasts to compare. An upload
        # with fewer can't populate the multi-analysis tabs, so say so plainly
        # instead of showing an empty (and useless) data-selection box.
        nCon <- length(.eds$contrasts())
        if (nCon < 2) {
            return(box(
                title = "Data Selection",
                width = NULL, solidHeader = TRUE, status = "warning",
                helpText(
                    if (nCon == 0)
                        paste0("The uploaded dataset contains no differential-expression ",
                               "contrasts, so multi-contrast DE analysis is unavailable.")
                    else
                        paste0("The uploaded dataset contains only one contrast",
                               " (\"", .eds$contrasts()[[1]], "\"). Multi-contrast DE ",
                               "analysis compares two or more contrasts, so these tabs ",
                               "are unavailable for this dataset.")
                )
            ))
        }
        # >= 2 contrasts: auto-loaded below; the Selected Datasets list shows them.
        return(NULL)
    }

    # Show the data selection box for genomics or when no data is loaded yet
    box(
        title = "Data Selection",
        width = NULL,
        solidHeader = TRUE,
        status = "success",
        tabBox(
            width = 12,
            tabPanel(
                title = "Browse Projects",
                icon = icon("folder-open"),
                fluidRow(
                    column(
                        width = 12,
                        helpText("Click Refresh to find available datasets from projects you have access to. Select up to 10 datasets to load."),
                        actionButton(
                            inputId = "multiRefreshProjects",
                            label = "Refresh Available Datasets",
                            icon = icon("sync")
                        ),
                        br(), br(),
                        uiOutput("multiAvailableDatasetsCheckboxUI"),
                        br(),
                        actionButton(
                            inputId = "multiLoadSelected",
                            label = "Load Selected Datasets",
                            icon = icon("download"),
                            class = "btn-success"
                        )
                    )
                )
            ),
            tabPanel(
                title = "Manual Import",
                icon = icon("edit"),
                fluidRow(
                    column(
                        width = 6,
                        helpText("You can upload an Excel file with two columns here. Column 1 must contain the links to the live reports, and the second column the name you'd like to give to each contrast."),
                        splitLayout(
                            fileInput(inputId = "multiImportFile", accept = ".xlsx", width = "85%", label = NULL),
                            downloadButton(outputId = "multiDlExampleFile", label = "See example import file", inline = T)
                        ),
                        actionButton(inputId = "multiImportButton", label = "Click to load the data"),
                        br(), br(),
                        lapply(1:12, function(i) {
                            splitLayout(
                                textInput(
                                    inputId = paste0("multiContrastLink", i),
                                    label = paste0("Contrast Link: ", i),
                                    width = "100%",
                                    placeholder = "https://fgcz-shiny.uzh.ch/exploreDEG/?data=..."
                                ),
                                textInput(
                                    inputId = paste0("multiContrastName", i),
                                    label = paste0("Contrast Name: ", i),
                                    width = "70%",
                                    placeholder = "Treated over Untreated"
                                )
                            )
                        })
                    )
                )
            )
        )
    )
})

# 0. Colour Picker Logic (Ported from multiDEG)
# User requested nice + colourblind friendly with >10 colours.
# Mixing Okabe-Ito (8) with Paul Tol's Muted/Bright (4) to get 12 distinct, accessible colours.
houseCols <- c(
    "#E69F00", # Orange (OI)
    "#56B4E9", # Sky Blue (OI)
    "#009E73", # Bluish Green (OI)
    "#F0E442", # Yellow (OI)
    "#0072B2", # Blue (OI)
    "#D55E00", # Vermilion (OI)
    "#CC79A7", # Reddish Purple (OI)
    "#332288", # Indigo (Tol)
    "#88CCEE", # Cyan (Tol)
    "#44AA99", # Teal (Tol)
    "#117733", # Green (Tol)
    "#999933"  # Olive (Tol)
)
houseCols <- paste0(houseCols, "FF")
n <- 60
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colourPaletteList2 <- list(
    "House colours" = houseCols,
    "Paired" = brewer.pal(12, "Paired"),
    "Set1" = brewer.pal(9, "Set1"),
    "Set2" = brewer.pal(8, "Set2"),
    "Set3" = brewer.pal(12, "Set3"),
    "Dark2" = brewer.pal(8, "Dark2"),
    "Dealer's choice" = sample(col_vector, 20)
)

output$multiColourPaletteUI <- renderUI({
    selectizeInput(
        inputId = "multiColourPalette",
        label = "Colour Palette",
        choices = names(colourPaletteList2),
        selected = names(colourPaletteList2)[1],
        multiple = TRUE,
        options = list(plugins = list("remove_button", "drag_drop"))
    )
})
outputOptions(output, "multiColourPaletteUI", suspendWhenHidden = FALSE)

observeEvent(
    {
        multiDeReactive$results
        input$multiColourPalette
    },
    {
        req(multiDeReactive$results)
        output$multiColourPickerUI <- renderUI({
            # Fallback to first palette if multiColourPalette not yet available
            selectedPalette <- input$multiColourPalette
            if (is.null(selectedPalette) || length(selectedPalette) == 0) {
                selectedPalette <- names(colourPaletteList2)[1]
            }
            lapply(seq_along(names(multiDeReactive$results)), function(i) {
                colourpicker::colourInput(
                    inputId = (paste0("multiGroupColour", names(multiDeReactive$results)[[i]]) %>% gsub(" ", "_", .)),
                    label = names(multiDeReactive$results)[[i]],
                    value = rep(as.character(unlist(colourPaletteList2[selectedPalette])), times = 5)[i],
                    palette = "square",
                    closeOnClick = TRUE,
                    returnName = TRUE
                )
            })
        })
        outputOptions(output, "multiColourPickerUI", suspendWhenHidden = FALSE)
    }
)

# 1. Browse Projects Logic - Simplified to use current project only
observeEvent(input$multiRefreshProjects, {
    # The demo, uploads and local files belong to no project, so there is
    # nothing to browse: Say so rather than doing nothing.
    if (is.null(multiProjectInfo$projectDir)) {
        showNotification(paste("There is no project to browse: This dataset wasn't opened from a SUSHI link.",
                               "Use the Manual Import tab to load result links instead."),
                         type = "warning", duration = 10)
        return(NULL)
    }

    projectDir <- multiProjectInfo$projectDir

    if (!dir.exists(projectDir)) {
        showNotification(paste("Could not find project directory:", projectDir), type = "warning")
        return(NULL)
    }

    # List all subdirectories in the project
    subdirs <- list.dirs(projectDir, full.names = TRUE, recursive = FALSE)

    # Filter for folders with DESeq2, EdgeR, or Limma in their name
    analysis_folders <- subdirs[grepl("DESeq2|EdgeR|Limma", basename(subdirs), ignore.case = TRUE)]

    # For each analysis folder, look for subfolders containing deResult.rds
    available <- lapply(analysis_folders, function(d) {
        name <- basename(d)

        # Determine analysis type from folder name
        type <- "Unknown"
        if (grepl("DESeq2", name, ignore.case = TRUE)) type <- "DESeq2"
        if (grepl("EdgeR", name, ignore.case = TRUE)) type <- "EdgeR"
        if (grepl("Limma", name, ignore.case = TRUE)) type <- "Limma"

        # Look for subfolders containing deResult.rds or deResult.qs2
        subfolders <- list.dirs(d, full.names = TRUE, recursive = FALSE)

        for (subfolder in subfolders) {
            has_deresult <- file.exists(file.path(subfolder, "deResult.rds"))
            has_qs2 <- file.exists(file.path(subfolder, "deResult.qs2"))

            if (has_deresult || has_qs2) {
                return(data.frame(
                    Name = name,
                    Path = subfolder,  # Use the subfolder path where deResult is located
                    Type = type,
                    stringsAsFactors = FALSE
                ))
            }
        }

        return(NULL)
    }) %>% bind_rows()

    if (is.null(available) || nrow(available) == 0) {
        showNotification("No datasets found in current project.", type = "warning")
        multiDeReactive$available <- NULL
        return(NULL)
    }

    multiDeReactive$available <- available
    showNotification(paste("Found", nrow(available), "datasets in project", basename(projectDir)),
                     type = "message")
})

# Render checkbox UI for available datasets
output$multiAvailableDatasetsCheckboxUI <- renderUI({
    req(multiDeReactive$available)
    req(nrow(multiDeReactive$available) > 0)

    available <- multiDeReactive$available

    # Create checkbox labels with dataset name and type
    choices <- setNames(
        available$Path,
        paste0(available$Name, " [", available$Type, "]")
    )

    # Pre-select the currently loaded dataset if it's in the list
    preselected <- NULL
    if (!is.null(multiProjectInfo$currentDataDir)) {
        # Normalize paths for comparison (remove trailing slashes)
        currentPath <- gsub("/$", "", multiProjectInfo$currentDataDir)
        availablePaths <- gsub("/$", "", available$Path)

        # Check if current path matches any available path
        matchIndex <- which(availablePaths == currentPath)
        if (length(matchIndex) > 0) {
            # Use the original path from available$Path (with or without trailing slash)
            preselected <- available$Path[matchIndex[1]]
        }
    }

    box(
        title = "Available Datasets",
        width = NULL,
        solidHeader = TRUE,
        status = "info",
        helpText(paste0("Found ", nrow(available), " datasets. Select up to 10 to load.")),
        checkboxGroupInput(
            inputId = "multiDatasetsCheckbox",
            label = NULL,
            choices = choices,
            selected = preselected
        )
    )
})

observeEvent(input$multiLoadSelected, {
    req(input$multiDatasetsCheckbox)
    req(length(input$multiDatasetsCheckbox) > 0)

    # Validate at least 2 datasets selected
    if (length(input$multiDatasetsCheckbox) < 2) {
        showNotification("Please select at least 2 datasets for multi-contrast analysis.", type = "warning")
        return(NULL)
    }

    # Validate selection limit
    if (length(input$multiDatasetsCheckbox) > 12) {
        showNotification("Please select 12 or fewer datasets.", type = "warning")
        return(NULL)
    }

    # Show waiter on the data selection box only
    waiter <- Waiter$new(
        id = "multiDataSelectionBoxUI",
        html = spin_fading_circles(),
        color = "#2d7a3e", 
        fadeout = TRUE
    )
    waiter$show()

    # Only paths the project browser offered: Shiny does not check checkbox
    # values against the choices, so a modified client could send any path.
    selectedPaths <- intersect(input$multiDatasetsCheckbox, multiDeReactive$available$Path)

    tryCatch({
        for (i in seq_along(selectedPaths)) {
            # Update waiter message
            waiter$update(
                html = tagList(
                    spin_fading_circles(),
                    h4(paste("Loading dataset", i, "of", length(selectedPaths)))
                )
            )

            fullPath <- selectedPaths[i]
            # Extract name from parent directory (contains DESeq2/EdgeR/Limma and date)
            # This ensures different methods of same contrast get unique names
            name <- basename(dirname(fullPath))

            resolved <- multi_resolve_dataset(fullPath, MULTI_DATA_ROOTS)
            if (is.null(resolved)) {
                showNotification(paste("No deResult found in", fullPath), type = "error")
                next
            }
            if (!multi_access_ok(resolved)) {
                showNotification(paste("You do not have access to the project of", name), type = "error")
                next
            }
            formatted <- multi_load_de_result(resolved$path)
            if (is.null(formatted)) {
                showNotification(paste("No deResult found in", fullPath), type = "error")
                next
            }
            formatted$source <- list(kind = "link", link = resolved$link)  # for saved sessions
            multiDeReactive$results[[name]] <- formatted
        }

        # Hide waiter on success
        waiter$hide()
        showNotification("Datasets loaded successfully!", type = "message")

        # Store current tab before menu update
        currentTab <- isolate(input$tabs)

        # Trigger menu update to show new multi-analysis tabs
        multiDeReactive$menuUpdateTrigger <- multiDeReactive$menuUpdateTrigger + 1

        # After a brief delay, restore the tab selection
        shinyjs::delay(100, {
            if (!is.null(currentTab)) {
                updateTabItems(session, "tabs", selected = currentTab)
            }
        })
    }, error = function(e) {
        # Hide waiter on error
        waiter$hide()
        showNotification(paste("Error loading datasets:", e$message), type = "error")
    })
})


# 2. Manual Import Logic (Ported from multiDEG)
observeEvent(input$multiImportFile, {
    req(input$multiImportFile)
    userTable <- openxlsx::read.xlsx(input$multiImportFile[1, "datapath"], colNames = F)
    if (userTable[1, 1] == "Link") {
        userTable <- userTable[-1, ]
    }
    lapply(1:nrow(userTable), function(i) {
        updateTextInput(
            session = session,
            inputId = paste0("multiContrastLink", i),
            value = userTable[i, 1]
        )
        updateTextInput(
            session = session,
            inputId = paste0("multiContrastName", i),
            value = userTable[i, 2]
        )
    })
})

output$multiDlExampleFile <- downloadHandler(
    filename = function() {
        paste0("example.xlsx")
    },
    content = function(file) {
        eg <- data.frame(
            "Link" = c(
                "https://fgcz-shiny.uzh.ch/exploreDEG/?data=p1234/o20000_DESeq2_KO--over--WT_2023-04-25--13-42-40/KO--over--WT/",
                "https://fgcz-shiny.uzh.ch/exploreDEG/?data=p1234/o20000_DESeq2_Treated--over--Untreated_2023-04-25--13-44-40/Treated--over--Untreated/"
            ),
            "Name" = c(
                "KO over WT",
                "Treated vs Untreated"
            )
        )
        openxlsx::write.xlsx(eg, file, rowNames = FALSE)
    }
)

observeEvent(input$multiImportButton, {
    withProgress(message = "Loading manually imported data...", value = 0, {
        inputTable <- data.frame("Link" = character(length = 10L), "Name" = character(length = 10L))
        for (i in 1:12) {
            inputTable[i, 1] <- input[[paste0("multiContrastLink", i)]]
            inputTable[i, 2] <- input[[paste0("multiContrastName", i)]]
        }
        inputTable <- inputTable[-which(inputTable$Link == ""), ]

        if (nrow(inputTable) == 0) {
            return(NULL)
        }

        # Loop through and load
        for (i in 1:nrow(inputTable)) {
            incProgress(1 / nrow(inputTable), detail = paste("Loading", inputTable[i, "Name"]))

            link <- inputTable[i, "Link"]
            name <- inputTable[i, "Name"]

            # Link -> result directory under a data root (strips spaces and an
            # exploreDE URL down to its ?data= path; refuses "..").
            resolved <- multi_resolve_dataset(link, MULTI_DATA_ROOTS)
            if (is.null(resolved)) {
                showNotification(paste("Could not find data for:", name), type = "warning")
                next
            }
            if (!multi_access_ok(resolved)) {
                showNotification(paste("You do not have access to the project of", name), type = "error")
                next
            }
            formatted <- multi_load_de_result(resolved$path)
            if (is.null(formatted)) {
                showNotification(paste("Could not find data for:", name), type = "warning")
                next
            }
            formatted$source <- list(kind = "link", link = resolved$link)  # for saved sessions
            multiDeReactive$results[[name]] <- formatted
        }
    })

    # Trigger menu update to show new multi-analysis tabs
    currentTab <- isolate(input$tabs)
    multiDeReactive$menuUpdateTrigger <- multiDeReactive$menuUpdateTrigger + 1

    # After a brief delay, restore the tab selection
    shinyjs::delay(100, {
        if (!is.null(currentTab)) {
            updateTabItems(session, "tabs", selected = currentTab)
        }
    })
})

# 3. Selected Datasets Management
output$multiSelectedDatasetsTable <- DT::renderDataTable({
    req(length(multiDeReactive$results) > 0)
    data.frame(
        Name = names(multiDeReactive$results),
        Type = sapply(multiDeReactive$results, function(x) x$dataType),
        stringsAsFactors = FALSE
    ) %>%
        DT::datatable(
            selection = "multiple",
            editable = FALSE,  # Disable inline editing - use text inputs below instead
            options = list(dom = "t")
        )
})

# Render rename text inputs dynamically
observeEvent(
    multiDeReactive$results,
    {
        req(multiDeReactive$results)
        req(length(multiDeReactive$results) > 0)

        output$multiRenameInputsUI <- renderUI({
            tagList(
                h5("Edit Dataset Names:"),
                lapply(seq_along(names(multiDeReactive$results)), function(i) {
                    currentName <- names(multiDeReactive$results)[[i]]
                    textInput(
                        inputId = paste0("multiRename_", i),
                        label = paste0("Dataset ", i, ":"),
                        value = currentName,
                        placeholder = "Enter dataset name...",
                        width = "100%"
                    )
                })
            )
        })
    },
    ignoreNULL = TRUE,
    ignoreInit = FALSE
)

# Handle "Update Names" button click
observeEvent(input$multiUpdateNames, {
    req(multiDeReactive$results)
    req(length(multiDeReactive$results) > 0)

    # Store current tab before making changes
    currentTab <- isolate(input$tabs)

    # Collect all new names from the text inputs
    newNames <- character(length(multiDeReactive$results))
    hasChanges <- FALSE

    for (i in seq_along(names(multiDeReactive$results))) {
        oldName <- names(multiDeReactive$results)[i]
        inputId <- paste0("multiRename_", i)
        newName <- trimws(input[[inputId]])

        # Validate the new name
        if (is.null(newName) || newName == "") {
            showNotification("Dataset names cannot be empty!", type = "warning")
            return()
        }

        newNames[i] <- newName

        if (newName != oldName) {
            hasChanges <- TRUE
        }
    }

    # Check for duplicate names
    if (any(duplicated(newNames))) {
        showNotification("Dataset names must be unique!", type = "warning")
        return()
    }

    # If no changes, just notify and return
    if (!hasChanges) {
        showNotification("No changes to apply.", type = "message")
        return()
    }

    # Clear figure reactives before updating to prevent stale observers crashing
    multiBoxplotFiguresReactive$boxplot <- NULL
    multiBoxplotFiguresReactive$barplot <- NULL
    tilePlots$top <- NULL
    tilePlots$custom <- NULL
    tilePlots$go <- NULL
    multiBubbleReactive$plot <- NULL
    multiUpsetReactive$plot <- NULL

    # Apply the changes by rebuilding the list
    temp_results <- list()
    for (i in seq_along(multiDeReactive$results)) {
        temp_results[[newNames[i]]] <- multiDeReactive$results[[i]]
    }

    # Update reactive values
    multiDeReactive$results <- temp_results
    multiDeReactive$menuUpdateTrigger <- multiDeReactive$menuUpdateTrigger + 1

    showNotification("Dataset names updated successfully!", type = "message")

    # Restore the current tab after a brief delay
    shinyjs::delay(100, {
        if (!is.null(currentTab)) {
            updateTabItems(session, "tabs", selected = currentTab)
        }
    })
})

observeEvent(input$multiRemoveSelected, {
    req(input$multiSelectedDatasetsTable_rows_selected)
    selected_indices <- input$multiSelectedDatasetsTable_rows_selected
    names_to_remove <- names(multiDeReactive$results)[selected_indices]

    # Create a temporary copy to modify
    temp_results <- multiDeReactive$results
    for (n in names_to_remove) {
        temp_results[[n]] <- NULL
    }
    # Assign back
    multiDeReactive$results <- temp_results

    showNotification("Selected datasets removed.", type = "message")
})

observeEvent(input$multiClearAll, {
    multiDeReactive$results <- list()
})

# 3b. Contrast Order Bucket (drag-and-drop reordering)
observeEvent(
    multiDeReactive$results,
    {
        output$multiContrastOrderBucket <- renderUI({
            req(multiDeReactive$results)
            req(length(multiDeReactive$results) > 0)
            sortable::rank_list(
                text = "Drag to reorder contrasts:",
                labels = as.list(names(multiDeReactive$results)),
                input_id = "multiContrastOrder",
                options = sortable::sortable_options(multiDrag = TRUE)
            )
        })
    },
    ignoreNULL = TRUE,
    ignoreInit = FALSE
)

observeEvent(input$multiApplyOrder, {
    req(multiDeReactive$results)
    req(input$multiContrastOrder)
    newOrder <- input$multiContrastOrder
    req(length(newOrder) == length(multiDeReactive$results))
    req(all(newOrder %in% names(multiDeReactive$results)))

    # Store current tab before making changes
    currentTab <- isolate(input$tabs)

    multiDeReactive$results <- multiDeReactive$results[newOrder]
    multiDeReactive$menuUpdateTrigger <- multiDeReactive$menuUpdateTrigger + 1

    showNotification("Contrast order updated.", type = "message")

    shinyjs::delay(100, {
        if (!is.null(currentTab)) {
            updateTabItems(session, "tabs", selected = currentTab)
        }
    })
})

# TRUE when at least one loaded contrast carries a usable raw p-value. FDR-only
# producers (SAINTexpress: pValue all NA, fdr = BFDR) fall back to FDR in the
# multi Tile / Upset tabs, mirroring eds$has_pvalue() in the single-contrast tabs.
multiHasPvalue <- reactive({
    any(vapply(multiDeReactive$results, function(x) {
        p <- x$seqAnno$pValue
        !is.null(p) && any(!is.na(p))
    }, logical(1)))
})

observeEvent(multiHasPvalue(), {
    req(length(multiDeReactive$results) > 0)
    if (multiHasPvalue()) {
        updateSelectInput(session, "pTypeTile", choices = c("FDR", "Raw"),
                          selected = isolate(input$pTypeTile) %||% "Raw")
        updateRadioButtons(session, "pTypeUpset", choices = c("FDR", "Raw"),
                           selected = isolate(input$pTypeUpset) %||% "FDR", inline = TRUE)
    } else {
        updateSelectInput(session, "pTypeTile", choices = c("FDR"), selected = "FDR")
        updateRadioButtons(session, "pTypeUpset", choices = c("FDR"), selected = "FDR", inline = TRUE)
    }
})

# 4. Proteomics Auto-Load (Triggered once)
observeEvent(inputDataReactive(), {
    eds <- inputDataReactive()$.eds
    req(eds)

    if (eds$omics_class() == "proteomics" && eds$backend_kind() == "se") {
        # Multi-* downstream modules still consume the legacy gStoreList shape
        # (seqAnnoList, countList, factors, factorLevels, gene_lookup, etc.).
        # Rebuild it on-demand via the package-level shim. Slated for proper
        # eds-native migration in the multi-* subsystem rewrite.
        data <- legacy_gStoreList_from_eds(eds, se = eds$backend())

        # Auto-populate multiDeReactive
        # The 'convert_proteomics_se' function returns 'seqAnnoList' which has per-contrast data.
        # We need to structure this as if they were separate datasets.

        # Actually, multiDEG expects a list of objects where each object has 'se', 'seqAnno', etc.
        # For proteomics, we have one SE but multiple contrasts in rowData.
        # We should create a 'dataset' object for each contrast.

        contrasts <- data$contrasts

        for (con in contrasts) {
            # SeqAnno for this contrast
            sa <- data$seqAnnoList[[con]]

            # Determine grouping info
            # Proteomics contrasts are usually "ConditionA vs ConditionB" or "ConditionA-ConditionB"
            # We need to extract this to set sampleGroup and refGroup
            # And we need a groupingName.

            # Try to find groupingName from metadata or default to first factor
            groupingName <- NULL
            if (!is.null(metadata(data$se)$param$groupingName)) {
                groupingName <- metadata(data$se)$param$groupingName
            } else if (length(data$factors) > 0) {
                groupingName <- data$factors[1]
            } else {
                # Fallback: Try to find a likely grouping column
                potential_groups <- colnames(data$dataset)[grep("condition|group|treatment|genotype|batch", colnames(data$dataset), ignore.case = TRUE)]
                if (length(potential_groups) > 0) {
                    groupingName <- potential_groups[1]
                } else {
                    # Last resort: use Name (one sample per group)
                    groupingName <- "Name"
                }
            }

            # Try to parse contrast for sample/ref groups
            # Common patterns: "A_vs_B", "A vs B", "A--over--B"
            sampleGroup <- NULL
            refGroup <- NULL

            if (grepl("_vs_", con)) {
                parts <- strsplit(con, "_vs_")[[1]]
                sampleGroup <- parts[1]
                refGroup <- parts[2]
            } else if (grepl(" vs ", con)) {
                parts <- strsplit(con, " vs ")[[1]]
                sampleGroup <- parts[1]
                refGroup <- parts[2]
            } else if (grepl("--over--", con)) {
                parts <- strsplit(con, "--over--")[[1]]
                sampleGroup <- parts[1]
                refGroup <- parts[2]
            } else {
                # Fallback: assume the contrast name is the sample group and we don't know ref?
                # Or just leave them NULL and let boxplot handle it (it might skip brackets)
                sampleGroup <- con
            }

            # Let's create a lightweight object
            pseudoObj <- list(
                dataType = "proteomics",
                se = data$se, # Shared SE
                seqAnno = sa,
                param = list(
                    runGO = FALSE,
                    featureLevel = "gene", # Default to gene/protein level
                    groupingName = groupingName,
                    sampleGroup = sampleGroup,
                    refGroup = refGroup,
                    grouping = if (!is.null(groupingName)) setNames(as.character(data$dataset[[groupingName]]), rownames(data$dataset)) else NULL
                ),
                dataset = data$dataset,
                countList = data$countList, # Pass pre-computed countList
                design = con,
                genes = data$genes,
                factors = data$factors,
                factorLevels = data$factorLevels,
                gene_lookup = data$gene_lookup,
                protein_lookup = data$protein_lookup,
                gene_choices = data$gene_choices,
                protein_choices = data$protein_choices,
                groupingName = groupingName
            )

            pseudoObj$source <- list(kind = "contrast", contrast = con)  # for saved sessions
            multiDeReactive$results[[con]] <- pseudoObj
        }

        # Reveal the multi-analysis tabs, preserving the current tab selection.
        # The sidebar reads hasMultiData() under isolate(), so without this bump
        # it stays on "Data Selection" only (mirrors the h5ad auto-load below).
        if (length(multiDeReactive$results) > 0) {
            currentTab <- isolate(input$tabs)
            multiDeReactive$menuUpdateTrigger <- multiDeReactive$menuUpdateTrigger + 1
            shinyjs::delay(100, {
                if (!is.null(currentTab)) updateTabItems(session, "tabs", selected = currentTab)
            })
        }

        # Also populate the Manual Import text boxes so user can see them?
        # Or just leave them in the "Selected Datasets" list.
        # The user requirement was: "Auto-populate the manual import inputs... allowing you to remove unwanted contrasts."
        # So we should update the text inputs.

        # Update inputs
        lapply(seq_along(contrasts), function(i) {
            if (i <= 10) {
                updateTextInput(session, paste0("multiContrastName", i), value = contrasts[i])
                updateTextInput(session, paste0("multiContrastLink", i), value = "Internal (Proteomics)")
            }
        })
    }
})

# 5. Uploaded .h5ad Auto-Load (genomics, multi-contrast)
# An uploaded canonical .h5ad is a single ExploreDataset carrying all of its
# contrasts, so — unlike the SUSHI genomics path — there is no per-contrast SE
# to browse for. Fan the dataset out into one Multi-DE object per contrast (the
# genomics analogue of the proteomics auto-load above), so the Data Selection
# tab offers exactly the contrasts stored in the uploaded file.
observeEvent(inputDataReactive(), {
    eds <- inputDataReactive()$.eds
    req(eds)
    req(eds$backend_kind() == "h5ad")
    req(eds$omics_class() == "genomics")
    req(length(eds$contrasts()) >= 2)            # multi-DE needs >=2 contrasts
    req(length(multiDeReactive$results) == 0)   # one-shot; don't clobber user edits

    contrasts <- eds$contrasts()

    for (con in contrasts) {
        obj <- tryCatch(
            legacy_multi_contrast_from_eds(eds, con),
            error = function(e) {
                showNotification(paste0("Could not load contrast '", con, "': ", e$message),
                                 type = "warning")
                NULL
            }
        )
        if (!is.null(obj)) {
            obj$source <- list(kind = "contrast", contrast = con)  # for saved sessions
            multiDeReactive$results[[con]] <- obj
        }
    }

    if (length(multiDeReactive$results) > 0) {
        # Reveal the multi-analysis tabs, preserving the current tab selection.
        currentTab <- isolate(input$tabs)
        multiDeReactive$menuUpdateTrigger <- multiDeReactive$menuUpdateTrigger + 1
        shinyjs::delay(100, {
            if (!is.null(currentTab)) updateTabItems(session, "tabs", selected = currentTab)
        })
    }
})
