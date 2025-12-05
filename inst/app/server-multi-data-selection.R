# Server logic for Multi DE Data Selection
# Note: multiProjectInfo is defined in app.R

# Conditional UI - Only show data selection box for genomics (not proteomics)
output$multiDataSelectionBoxUI <- renderUI({
    # Check if we have proteomics data loaded
    if (!is.null(inputDataReactive()) && inputDataReactive()$dataType == "proteomics") {
        # Don't show the data selection box for proteomics
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
                        lapply(1:10, function(i) {
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
catalystCols <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
    "#FF7B00", "#FDC362", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
    "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
    "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
    "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
)
catalystCols <- c(catalystCols[c(seq(1, 30, by = 2), seq(2, 30, by = 2))])
cc2 <- catalystCols
cc2 <- colorspace::darken(cc2, 0.4)
catalystCols <- c(catalystCols, cc2)
catalystCols <- paste0(catalystCols, "FF")
n <- 60
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colourPaletteList <- list(
    "Catalyst colours" = catalystCols,
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
        choices = names(colourPaletteList),
        selected = names(colourPaletteList)[1],
        multiple = TRUE,
        options = list(plugins = list("remove_button", "drag_drop"))
    )
})

observeEvent(
    {
        multiDeReactive$results
        input$multiColourPalette
    },
    {
        req(multiDeReactive$results)
        output$multiColourPickerUI <- renderUI({
            lapply(seq_along(names(multiDeReactive$results)), function(i) {
                colourpicker::colourInput(
                    inputId = (paste0("multiGroupColour", names(multiDeReactive$results)[[i]]) %>% gsub(" ", "_", .)),
                    label = names(multiDeReactive$results)[[i]],
                    value = rep(as.character(unlist(colourPaletteList[input$multiColourPalette])), times = 5)[i],
                    palette = "square",
                    closeOnClick = TRUE,
                    returnName = TRUE
                )
            })
        })
    }
)

# 1. Browse Projects Logic - Simplified to use current project only
observeEvent(input$multiRefreshProjects, {
    req(multiProjectInfo$projectDir)

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
    if (length(input$multiDatasetsCheckbox) > 10) {
        showNotification("Please select 10 or fewer datasets.", type = "warning")
        return(NULL)
    }

    # Show waiter on the data selection box only
    waiter <- Waiter$new(
        id = "multiDataSelectionBoxUI",
        html = spin_fading_circles(),
        color = "#527085"
    )
    waiter$show()

    selectedPaths <- input$multiDatasetsCheckbox

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

            # Load the dataset
            if (file.exists(file.path(fullPath, "deResult.rds"))) {
                se <- readRDS(file.path(fullPath, "deResult.rds"))
                formatted <- convert_genomics_se(se, fullPath)

                # Compatibility Fix: Ensure log2Ratio exists
                if (!"log2Ratio" %in% colnames(formatted$seqAnno)) {
                    lfcCol <- grep("^(log2?[._ ]?(fold)?[._ ]?(change|ratio|fc)|logfc)$",
                                   colnames(formatted$seqAnno), value = TRUE, ignore.case = TRUE)
                    if (length(lfcCol) > 0) {
                        formatted$seqAnno$log2Ratio <- formatted$seqAnno[[lfcCol[1]]]
                    }
                }

                # Extract ORA and GSEA data from metadata (like multiDEG does)
                if (!is.null(formatted$param$runGO) && formatted$param$runGO) {
                    # ORA data
                    if (!is.null(metadata(se)$enrichResult)) {
                        ora <- metadata(se)$enrichResult

                        # Add Label column to each ORA result
                        for (ont in names(ora)) {
                            for (sig in names(ora[[ont]])) {
                                if (nrow(ora[[ont]][[sig]]) >= 2) {
                                    er <- ora[[ont]][[sig]]
                                    er@result$Label <- ""
                                    for (j in seq_along(er@result$Description)) {
                                        desc <- er@result$Description[j]
                                        if (nchar(desc) > 35) {
                                            desc <- substr(desc, 1, 45)
                                        }
                                        er@result$Label[j] <- paste0(er@result$ID[j], "\n", desc)
                                    }
                                    ora[[ont]][[sig]] <- er
                                }
                            }
                        }
                        formatted$ora <- ora
                    }

                    # GSEA data
                    if (!is.null(metadata(se)$enrichResultGSEA)) {
                        gsea <- metadata(se)$enrichResultGSEA

                        # Add Label column to each GSEA result
                        for (ont in names(gsea)) {
                            if (nrow(gsea[[ont]]) >= 2) {
                                er <- gsea[[ont]]
                                er@result$Label <- ""
                                for (j in seq_along(er@result$Description)) {
                                    desc <- er@result$Description[j]
                                    if (nchar(desc) > 35) {
                                        desc <- substr(desc, 1, 45)
                                    }
                                    er@result$Label[j] <- paste0(er@result$ID[j], "\n", desc)
                                }
                                gsea[[ont]] <- er
                            }
                        }
                        formatted$gsea <- gsea
                    }
                }

                multiDeReactive$results[[name]] <- formatted
            }
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
        for (i in 1:10) {
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

            # Clean link logic (from multiDEG)
            if (grepl(" ", link)) {
                link <- gsub(" ", "", link)
            }
            if (isTRUE(as.logical(grep("^http", link)))) {
                link <- gsub(".*\\?data\\=", "", link)
            }

            # Construct path
            # We need to handle the path construction carefully.
            # multiDEG had logic to check multiple roots.
            urlDataRoot <- c("/srv/gstore/projects", "/srv/GT/analysis/course_sushi/public/gstore/projects")

            # If link is already a full path?
            # Usually it's p1234/o5678...

            found <- FALSE
            for (root in urlDataRoot) {
                fullPath <- file.path(root, link)
                # Check for deResult.rds
                if (file.exists(file.path(fullPath, "deResult.rds"))) {
                    se <- readRDS(file.path(fullPath, "deResult.rds"))
                    formatted <- convert_genomics_se(se, fullPath)

                    # Compatibility Fix: Ensure log2Ratio exists (for multiDEG tabs)
                    if (!"log2Ratio" %in% colnames(formatted$seqAnno)) {
                        lfcCol <- grep("^(log2?[._ ]?(fold)?[._ ]?(change|ratio|fc)|logfc)$", colnames(formatted$seqAnno), value = TRUE, ignore.case = TRUE)
                        if (length(lfcCol) > 0) {
                            formatted$seqAnno$log2Ratio <- formatted$seqAnno[[lfcCol[1]]]
                        }
                    }

                    # Extract ORA and GSEA data from metadata (like multiDEG does)
                    if (!is.null(formatted$param$runGO) && formatted$param$runGO) {
                        # ORA data
                        if (!is.null(metadata(se)$enrichResult)) {
                            ora <- metadata(se)$enrichResult

                            # Add Label column to each ORA result (copied from multiDEG)
                            for (ont in names(ora)) {
                                for (sig in names(ora[[ont]])) {
                                    if (nrow(ora[[ont]][[sig]]) >= 2) {
                                        er <- ora[[ont]][[sig]]
                                        er@result$Label <- ""
                                        for (j in seq_along(er@result$Description)) {
                                            desc <- er@result$Description[j]
                                            if (nchar(desc) > 35) {
                                                desc <- substr(desc, 1, 45)
                                            }
                                            er@result$Label[j] <- paste0(er@result$ID[j], "\n", desc)
                                        }
                                        ora[[ont]][[sig]] <- er
                                    }
                                }
                            }
                            formatted$ora <- ora
                        }

                        # GSEA data
                        if (!is.null(metadata(se)$enrichResultGSEA)) {
                            gsea <- metadata(se)$enrichResultGSEA

                            # Add Label column to each GSEA result (copied from multiDEG)
                            for (ont in names(gsea)) {
                                if (nrow(gsea[[ont]]) >= 2) {
                                    er <- gsea[[ont]]
                                    er@result$Label <- ""
                                    for (j in seq_along(er@result$Description)) {
                                        desc <- er@result$Description[j]
                                        if (nchar(desc) > 35) {
                                            desc <- substr(desc, 1, 45)
                                        }
                                        er@result$Label[j] <- paste0(er@result$ID[j], "\n", desc)
                                    }
                                    gsea[[ont]] <- er
                                }
                            }
                            formatted$gsea <- gsea
                        }
                    }

                    multiDeReactive$results[[name]] <- formatted
                    found <- TRUE
                    break
                }
            }

            if (!found) {
                showNotification(paste("Could not find data for:", name), type = "warning")
            }
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
            editable = TRUE,  # Make all cells editable by default
            options = list(dom = "t")
        )
})

# Handle edits to the Selected Datasets table
observeEvent(input$multiSelectedDatasetsTable_cell_edit, {
    info <- input$multiSelectedDatasetsTable_cell_edit
    row <- info$row
    col <- info$col
    newValue <- info$value

    # Only allow editing the Name column (column 0)
    if (col == 0) {
        oldName <- names(multiDeReactive$results)[row]
        newName <- newValue

        # Don't allow empty names
        if (is.null(newName) || newName == "" || trimws(newName) == "") {
            showNotification("Dataset name cannot be empty!", type = "warning")
            return(NULL)
        }

        # Check if new name already exists
        if (newName %in% names(multiDeReactive$results) && newName != oldName) {
            showNotification("A dataset with this name already exists!", type = "warning")
            return(NULL)
        }

        # Rename by changing the name in place to preserve order
        temp_results <- multiDeReactive$results
        result_names <- names(temp_results)
        result_names[row] <- newName
        names(temp_results) <- result_names

        # Reassign to trigger reactive updates
        multiDeReactive$results <- temp_results
        showNotification(paste("Renamed to:", newName), type = "message")
    } else {
        # Ignore edits to Type column
        showNotification("Type column cannot be edited.", type = "warning")
    }
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

# 4. Proteomics Auto-Load (Triggered once)
observeEvent(inputDataReactive(), {
    req(inputDataReactive())
    data <- inputDataReactive()

    if (data$dataType == "proteomics") {
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
                design = con
            )

            multiDeReactive$results[[con]] <- pseudoObj
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
