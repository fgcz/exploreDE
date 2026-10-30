# Helper functions

#' kegg_species_label
#' Maps a `ref_build`/`refBuild` string's leading "Genus_species" segment
#' (e.g. `"Homo_sapiens/GENCODE/..."`) to the KEGG-organism display name the
#' pathway tabs use to decide whether KEGG lookups are supported. Self-
#' contained so exploreDE does not need `ezRun::getSpecies()` just for this —
#' `ezRun` is otherwise optional (`Suggests`), and a fresh install without it
#' would hard-error on the four call sites this replaces (app.R, server-kegg.R,
#' server-goTile.R, server-heatmap2.R). Only the three organisms those call
#' sites ever recognised are mapped; anything else returns `NULL`, matching
#' their existing `switch(spp, "Human"=, "Mouse"=, "Rat"=, NULL)` behaviour.
#'
#' @param ref_build A `ref_build`/`refBuild` string, or `NA`/`NULL`.
#' @return `"Human"`, `"Mouse"`, `"Rat"`, or `NULL`.
#' @export
kegg_species_label <- function(ref_build) {
  if (is.null(ref_build) || is.na(ref_build) || !nzchar(ref_build)) return(NULL)
  genus_species <- strsplit(ref_build, "/", fixed = TRUE)[[1]][1]
  switch(genus_species,
    Homo_sapiens      = "Human",
    Mus_musculus      = "Mouse",
    Rattus_norvegicus = "Rat",
    NULL)
}

#' go_annotation_columns
#' The GO annotation columns a feature table carries, out of `GO BP`, `GO MF`,
#' and `GO CC`. Any subset is valid: Non-reference organisms often have only
#' one ontology. An empty result means the GO tabs are hidden, as KEGG is for
#' unsupported species.
#'
#' @param df A feature or DE table (`data.frame`).
#' @return Character vector of the GO columns present, possibly empty.
#' @export
go_annotation_columns <- function(df) {
  intersect(c("GO BP", "GO MF", "GO CC"), colnames(df))
}

#' features_in_go_term
#' Which rows of a feature table are annotated to one GO term, looking in
#' whichever GO columns the table has (see [go_annotation_columns()]). Each
#' column holds `"; "`-separated GO ids; `NA` means no annotation.
#'
#' @param df A feature or DE table (`data.frame`).
#' @param go_id A GO id, e.g. `"GO:0006955"`.
#' @return Logical vector, one per row; all `FALSE` without GO columns.
#' @export
features_in_go_term <- function(df, go_id) {
  hit <- rep(FALSE, nrow(df))
  for (col in go_annotation_columns(df)) {
    x <- as.character(df[[col]])
    hit <- hit | (!is.na(x) & grepl(go_id, x, fixed = TRUE))
  }
  hit
}

#' batch_correction_design
#' Build a design matrix for limma::removeBatchEffect that preserves biological
#' signal. Without a design matrix, removeBatchEffect assumes a one-group
#' experiment and can remove biology that correlates with batch.
#'
#' Never throws — returns a tagged result so callers can notify the user and
#' decide whether to skip batch correction or proceed without a design matrix.
#'
#' @param dataset Sample metadata data.frame
#' @param designColumn Column name to build the design from. "None" / empty
#'   means "no design" (batch correction still runs, without preserving any
#'   biology).
#' @param batchColumn Column name of the batch factor — required for the
#'   collinearity check.
#' @return list(design = <matrix or NULL>, skip = <logical>, msg = <chr or NULL>)
#'   - skip = TRUE  → caller must abort batch correction entirely
#'   - skip = FALSE, design = NULL → run removeBatchEffect without a design
#'   - skip = FALSE, design = <matrix> → run with the design
#' @export
batch_correction_design <- function(dataset, designColumn, batchColumn) {
  null_design <- list(design = NULL, skip = FALSE, msg = NULL)
  skipped     <- function(msg) list(design = NULL, skip = TRUE,  msg = msg)
  warned      <- function(msg) list(design = NULL, skip = FALSE, msg = msg)

  if (is.null(designColumn) || !nzchar(designColumn) || designColumn == "None") {
    return(null_design)
  }
  if (!designColumn %in% colnames(dataset)) {
    return(skipped(sprintf(
      "Design column '%s' missing from sample metadata — batch correction skipped.",
      designColumn
    )))
  }

  grp <- factor(dataset[[designColumn]])
  if (nlevels(grp) < 2) {
    return(warned(sprintf(
      "Design factor '%s' has only one level — batch correction applied without a design matrix.",
      designColumn
    )))
  }

  design <- tryCatch(stats::model.matrix(~grp), error = function(e) e)
  if (inherits(design, "error")) {
    return(skipped(sprintf(
      "Could not build design matrix from '%s': %s — batch correction skipped.",
      designColumn, conditionMessage(design)
    )))
  }

  if (!is.null(batchColumn) && nzchar(batchColumn) && batchColumn != "None" &&
      batchColumn %in% colnames(dataset)) {
    batch <- factor(dataset[[batchColumn]])
    if (nlevels(batch) >= 2) {
      batch_mm <- tryCatch(stats::model.matrix(~batch), error = function(e) NULL)
      if (!is.null(batch_mm)) {
        combined <- cbind(batch_mm, design[, -1, drop = FALSE])
        if (qr(combined)$rank < ncol(combined)) {
          return(warned(sprintf(
            "Batch ('%s') and design ('%s') are confounded — batch correction applied without a design matrix.",
            batchColumn, designColumn
          )))
        }
      }
    }
  }

  list(design = design, skip = FALSE, msg = NULL)
}

#' apply_batch_correction
#' Graceful-failure wrapper around limma::removeBatchEffect. Never throws:
#' on any validation or limma error, returns the original counts unchanged
#' and surfaces a Shiny notification.
#'
#' @param counts Numeric matrix (samples in columns)
#' @param dataset Sample metadata (nrow = ncol(counts))
#' @param batchCol Column name of batch factor, or "None" to skip entirely
#' @param designCol Column name of design factor, or "None"
#' @return The (possibly batch-corrected) counts matrix
#' @export
apply_batch_correction <- function(counts, dataset, batchCol, designCol) {
  if (is.null(batchCol) || !nzchar(batchCol) || batchCol == "None") return(counts)

  notify <- function(msg, type = "warning") {
    if (requireNamespace("shiny", quietly = TRUE) &&
        !is.null(shiny::getDefaultReactiveDomain())) {
      shiny::showNotification(msg, type = type, duration = 8)
    } else {
      warning(msg, call. = FALSE)
    }
  }

  if (!batchCol %in% colnames(dataset) ||
      nlevels(factor(dataset[[batchCol]])) < 2) {
    notify(sprintf("Batch factor '%s' has <2 levels — batch correction skipped.", batchCol))
    return(counts)
  }

  info <- batch_correction_design(dataset, designCol, batchCol)
  if (!is.null(info$msg)) {
    notify(info$msg, type = if (isTRUE(info$skip)) "error" else "warning")
  }
  if (isTRUE(info$skip)) return(counts)

  args <- list(x = counts, batch = dataset[[batchCol]])
  if (!is.null(info$design)) args$design <- info$design

  tryCatch(
    do.call(limma::removeBatchEffect, args),
    error = function(e) {
      notify(sprintf("Batch correction failed: %s — showing uncorrected counts.",
                     conditionMessage(e)),
             type = "error")
      counts
    }
  )
}

#' standardise_columns
#' Fuzzy-match incoming column names to canonical names.
#' Upstream data may call columns anything; this maps them to a known schema
#' so downstream code can rely on consistent names. Patterns come from the one
#' column dictionary (`inst/dictionaries/column_patterns.json`); DE concepts
#' are renamed to the legacy names the SE-backed app code reads
#' (`log2Ratio` / `pValue` / `fdr`). A column is only renamed when the
#' deciding pattern matches exactly one column.
#' @param df A data.frame (typically from rowData contrast)
#' @return The same data.frame with standardised column names
#' @export
standardise_columns <- function(df) {
  # dictionary concept -> legacy column name it is renamed to
  rowdata_concepts <- c("protein_Id", "peptide_Id", "site", "gene_name", "description",
                        "fasta.id", "modelName", "nrPeptides")
  de_concepts <- c(effect_log2fc = "log2Ratio", score_pvalue = "pValue", score_fdr = "fdr")

  orig    <- colnames(df)
  current <- orig
  rename <- function(patterns, to) {
    src <- eds_match_col(orig, patterns, unique = TRUE)
    if (!is.null(src)) current[match(src, orig)] <<- to
  }
  for (concept in rowdata_concepts) rename(eds_rowdata_column_patterns[[concept]], concept)
  for (concept in names(de_concepts)) rename(eds_de_column_patterns[[concept]], de_concepts[[concept]])
  colnames(df) <- current
  df
}

#' build_feature_id
#' Build a unique, human-readable feature identifier for proteomics data.
#' Uses canonical column names (after standardisation).
#' @param df A standardised data.frame with canonical column names
#' @return Character vector of unique feature IDs
#' @export
build_feature_id <- function(df) {
  pid <- df$protein_Id

  # If protein_Id is already unique, use it directly
  if (!anyDuplicated(pid)) return(pid)

  # Try known sub-feature columns (canonical names) in priority order
  # Always prepend protein_Id so users see which protein a sub-feature belongs to
  for (col in intersect(c("site", "peptide_Id"), colnames(df))) {
    combined <- paste(pid, gsub("~", "_", df[[col]]), sep = "_")
    if (!anyDuplicated(combined)) return(combined)
  }

  # Fallback: cleaned raw rownames
  gsub("~[^~]+~", "_", rownames(df))
}

#' collapse_guides_to_genes
#' Collapse sgRNA-level count matrix to gene-level using row sums.
#' Filters out negative control genes (matching ^negCtrl).
#' @param counts A numeric matrix with sgRNA rows and sample columns
#' @param sgRNA_to_gene Named character vector: names are sgRNA feature IDs, values are gene names
#' @return A numeric matrix with gene-level rows
#' @export
collapse_guides_to_genes <- function(counts, sgRNA_to_gene, exclude = NULL) {
  common <- intersect(rownames(counts), names(sgRNA_to_gene))
  counts <- counts[common, , drop = FALSE]
  collapsed <- rowsum(counts, group = sgRNA_to_gene[common], reorder = FALSE)
  if (length(exclude) > 0) {
    collapsed <- collapsed[!rownames(collapsed) %in% exclude, , drop = FALSE]
  }
  collapsed
}

#' convert_genomics_se
#' @export
convert_genomics_se <- function(se, dataDir, param = NULL) {
  # Delegates to the ExploreDataset hierarchy + legacy_gStoreList_from_eds()
  # shim (byte-equivalent to the pre-2026-05 converter it replaced; that
  # converter and its equivalence tests are archived in
  # plans/archive/legacy_converters/).
  eds <- new_explore_dataset(se)
  return(legacy_gStoreList_from_eds(eds, se = se, dataDir = dataDir, param = param))
}

#' convert_proteomics_se
#' @export
convert_proteomics_se <- function(se) {
  # Phase 1 (2026-05-15): body delegates to the ExploreDataset hierarchy +
  # legacy_gStoreList_from_eds() shim. See convert_genomics_se() for context.
  eds <- new_explore_dataset(se)
  # eds$backend() is the SE the dataset was built from (after any layout
  # conversion in new_explore_dataset()), so the shim reads the same frames.
  return(legacy_gStoreList_from_eds(eds, se = eds$backend()))
}

#' generateParams
#' @export
generateParams <- function(
  sampleGroup = "diff", refGroup = "undiff", groupingName = "Condition",
  comparison = "diff--over--undiff", featureLevel = "gene",
  runGO = TRUE, fdrThreshORA = 0.05, fdrThreshGSEA = 0.05, pValThreshGO = 0.01,
  log2RatioThreshGO = 0, backgroundExpression = 10, refBuild = "GRCh38"
) {
  return(list(
    sampleGroup = sampleGroup,
    refGroup = refGroup,
    groupingName = groupingName,
    comparison = comparison,
    featureLevel = featureLevel,
    runGO = runGO,
    fdrThreshORA = fdrThreshORA,
    fdrThreshGSEA = fdrThreshGSEA,
    pValThreshGO = pValThreshGO,
    log2RatioThreshGO = log2RatioThreshGO,
    backgroundExpression = backgroundExpression,
    refBuild = refBuild
  ))
}

#' is_url
#' Test whether a dataDir string is a URL (http/https) rather than a local path.
#' @param dataDir Character string
#' @return Logical scalar
#' @export
is_url <- function(dataDir) {
  return(grepl("^https?://", dataDir))
}


# =============================================================================
# Boxplot / barplot helpers (from server-boxplot.R, server-multi-boxplot.R)
# =============================================================================

#' dataSummary
#' Summarise counts by mean / SD / SE across grouping variables.
#' @param data data.frame containing at least `varname` and the `groupnames` columns
#' @param varname Column to summarise
#' @param groupnames Character vector of grouping columns
#' @return data.frame with one row per group combination
#' @export
dataSummary <- function(data, varname, groupnames) {
  summary_func <- function(x, col) {
    c(
      mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE),
      se = sd(x[[col]], na.rm = TRUE) / sqrt(sum(!is.na(x[[col]])))
    )
  }
  data_sum <- plyr::ddply(
    data, groupnames,
    .fun = summary_func, varname
  )
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# =============================================================================
# Heatmap helpers (from server-heatmap2.R)
# =============================================================================

#' select_genes_by_sd
#' Pick the top-N most variable features based on row-wise standard deviation.
#' Useful when ranking CountQC features without DE statistics.
#' @param counts Numeric matrix (features x samples)
#' @param top_n Number of features to keep (default 2000)
#' @return Character vector of rownames for the top features
#' @export
select_genes_by_sd <- function(counts, top_n = 2000) {
  gene_sds <- apply(counts, 1, sd, na.rm = TRUE)
  top_genes <- names(head(sort(gene_sds, decreasing = TRUE), n = top_n))
  return(top_genes)
}


# =============================================================================
# Feature bucket API (from server-featureBucket.R)
# NOTE: bucket_add/remove/clear and resolve_features reference the Shiny
# reactives `featureBucket`, `seqAnnoReactive`, and `inputDataReactive()`.
# They will error if called outside a Shiny session.
# =============================================================================

#' bucket_add
#' Add features to the central feature bucket reactive store.
#' @param new_features Character vector of feature IDs (or gene-level names in CrisprCountQC collapsed mode)
#' @param source Label for telemetry / debugging
#' @export
bucket_add <- function(new_features, source = "unknown") {
  valid_ids <- seqAnnoReactive$sa$feature_id
  # Also accept gene-level names for CrisprCountQC collapsed mode
  if (!is.null(inputDataReactive()$sgRNA_to_gene)) {
    valid_ids <- c(valid_ids, unique(inputDataReactive()$sgRNA_to_gene))
  }
  valid <- new_features[new_features %in% valid_ids & new_features != ""]
  current <- isolate(featureBucket$features)
  to_add <- setdiff(valid, current)
  if (length(to_add) > 0) {
    featureBucket$features <- c(current, to_add)
    featureBucket$version <- isolate(featureBucket$version) + 1L
  }
}

#' bucket_remove
#' Remove features from the central feature bucket.
#' @param to_remove Character vector of feature IDs
#' @param source Label for telemetry
#' @export
bucket_remove <- function(to_remove, source = "unknown") {
  current <- isolate(featureBucket$features)
  remaining <- current[!current %in% to_remove]
  if (length(remaining) != length(current)) {
    featureBucket$features <- remaining
    featureBucket$version <- isolate(featureBucket$version) + 1L
  }
}

#' bucket_clear
#' Empty the central feature bucket.
#' @param source Label for telemetry
#' @export
bucket_clear <- function(source = "unknown") {
  featureBucket$features <- character(0)
  featureBucket$version <- isolate(featureBucket$version) + 1L
}

#' parse_gene_text
#' Split a textarea string of gene names by whitespace / comma / slash / newline / tab.
#' @param text Character scalar
#' @return Character vector of non-empty tokens
#' @export
parse_gene_text <- function(text) {
  if (is.null(text) || text == "") return(character(0))
  genes <- strsplit(text, "\\ +|\\/+|\\,+|\n+|\r+|\n\r+|\t+")[[1]]
  genes[genes != ""]
}

#' resolve_features
#' Resolve gene symbols / protein IDs to feature_ids when not in feature_id mode.
#' @param input_genes Character vector of user input tokens
#' @param search_mode One of "feature_id", "gene_name", "protein_Id"
#' @return Character vector of feature_ids
#' @export
resolve_features <- function(input_genes, search_mode = "feature_id") {
  if (search_mode == "gene_name" && !is.null(inputDataReactive()$gene_lookup)) {
    expanded <- unlist(inputDataReactive()$gene_lookup[input_genes], use.names = FALSE)
    direct <- input_genes[input_genes %in% seqAnnoReactive$sa$feature_id]
    return(unique(c(expanded, direct)))
  }
  if (search_mode == "protein_Id" && !is.null(inputDataReactive()$protein_lookup)) {
    expanded <- unlist(inputDataReactive()$protein_lookup[input_genes], use.names = FALSE)
    direct <- input_genes[input_genes %in% seqAnnoReactive$sa$feature_id]
    return(unique(c(expanded, direct)))
  }
  return(input_genes)
}


# =============================================================================
# Replicability helpers (from server-replicability.R)
# NOTE: reference Shiny `input` from the enclosing scope.
# =============================================================================

#' make_rep_heatmap
#' Render the replicability correlation heatmap.
#' @param data List containing `corr` (correlation matrix) and `dataset` (sample metadata)
#' @param title_suffix String appended to heatmap title
#' @export
make_rep_heatmap <- function(data, title_suffix = "") {
    req(data)

    # Prepare annotations
    ha_top <- NULL
    ha_left <- NULL
    row_names_gp <- gpar(fontsize = input$replicabilityTextSize)
    column_names_gp <- gpar(fontsize = input$replicabilityTextSize)

    if (!is.null(input$replicabilityFactors)) {
        dataset <- data$dataset
        sample_names <- rownames(data$corr)
        dataset <- dataset[match(sample_names, dataset$Name), ]

        # Get colors
        colourList <- lapply(input$replicabilityFactors, function(f) {
            levs <- levels(as.factor(dataset[[f]]))
            colors <- vapply(levs, function(l) {
                col_input <- input[[paste0("GroupColour", f, l)]]
                if (is.null(col_input)) {
                    return("#808080FF")
                }
                tryCatch(paste0(col2hex(col_input), "FF"), error = function(e) "#808080FF")
            }, character(1))
            names(colors) <- levs
            return(colors)
        }) %>% setNames(input$replicabilityFactors)

        # Color labels based on the first selected factor
        factor_col <- input$replicabilityFactors[1]
        factor_levels <- levels(as.factor(dataset[[factor_col]]))
        colors <- vapply(factor_levels, function(l) {
            col_input <- input[[paste0("GroupColour", factor_col, l)]]
            if (is.null(col_input)) {
                return("#808080FF")
            }
            tryCatch(paste0(col2hex(col_input), "FF"), error = function(e) "#808080FF")
        }, character(1))

        sample_colors <- colors[as.character(dataset[[factor_col]])]
        names(sample_colors) <- dataset$Name

        row_names_gp <- gpar(fontsize = input$replicabilityTextSize, col = sample_colors)
        column_names_gp <- gpar(fontsize = input$replicabilityTextSize, col = sample_colors)

        anno_df <- dataset[, input$replicabilityFactors, drop = FALSE]

        # Legend shows the Colours-tab display labels (plain text)
        annoLegendParam <- lapply(input$replicabilityFactors, function(f) {
            levs <- names(colourList[[f]])
            list(at = levs, labels = unname(level_labels(input, f, levs, mode = "plain")))
        }) %>% setNames(input$replicabilityFactors)

        ha_top <- HeatmapAnnotation(
            df = anno_df,
            col = colourList,
            show_legend = TRUE,
            annotation_legend_param = annoLegendParam,
            annotation_name_gp = gpar(fontsize = input$replicabilityTextSize, fontface = "bold"),
            gp = gpar(col = input$replicabilityBorderColor1, lwd = input$replicabilityBorderSize1 / 10)
        )

        ha_left <- rowAnnotation(
            df = anno_df,
            col = colourList,
            show_legend = FALSE,
            annotation_name_gp = gpar(fontsize = input$replicabilityTextSize, fontface = "bold"),
            gp = gpar(col = input$replicabilityBorderColor1, lwd = input$replicabilityBorderSize1 / 10)
        )
    }

    # Get color breaks
    at <- c(input$replicabilityAtLow, input$replicabilityAtMid, input$replicabilityAtHigh)

    # Determine color scheme - use palette if selected, otherwise use manual colors
    if (!is.null(input$replicabilityPalette) && input$replicabilityPalette != "None") {
        palette_colours <- brewer.pal(9, input$replicabilityPalette)
        selected_colours <- palette_colours[c(1, 5, 9)]
        if (input$replicabilityPaletteRev) {
            heatmap_colours <- colorRamp2(at, selected_colours)
        } else {
            heatmap_colours <- colorRamp2(at, rev(selected_colours))
        }
    } else {
        heatmap_colours <- colorRamp2(at, c(input$replicabilityColorLow, input$replicabilityColorMid, input$replicabilityColorHigh))
    }

    hm <- ComplexHeatmap::Heatmap(
        data$corr,
        name = "Correlation",
        col = heatmap_colours,
        top_annotation = ha_top,
        left_annotation = ha_left,
        show_row_names = input$replicabilityShowNames,
        show_column_names = input$replicabilityShowNames,
        row_names_gp = row_names_gp,
        column_names_gp = column_names_gp,
        heatmap_legend_param = list(title_gp = gpar(fontsize = input$replicabilityTextSize, fontface = "bold")),
        rect_gp = gpar(col = input$replicabilityBorderColor2, lwd = input$replicabilityBorderSize2 / 10),
        column_title = paste0("Correlation Heatmap: ", title_suffix),
        column_title_gp = gpar(fontsize = input$replicabilityTextSize + 2, fontface = "bold")
    )

    draw(hm, merge_legend = TRUE, padding = unit(c(10, 10, 10, 10), "mm"))
}

#' make_rep_dendrogram
#' Render the replicability dendrogram with coloured labels / coloured bars.
#' @param data List containing `hclust` (hclust object) and `dataset` (sample metadata)
#' @param title_suffix String appended to title
#' @export
make_rep_dendrogram <- function(data, title_suffix = "") {
    req(data)
    dend <- as.dendrogram(data$hclust)

    if (!is.null(input$replicabilityFactors)) {
        dataset <- data$dataset
        factor_col <- input$replicabilityFactors[1]
        factor_levels <- levels(as.factor(dataset[[factor_col]]))
        colors <- vapply(factor_levels, function(l) {
            col_input <- input[[paste0("GroupColour", factor_col, l)]]
            if (is.null(col_input)) {
                return("#808080FF")
            }
            tryCatch(paste0(col2hex(col_input), "FF"), error = function(e) "#808080FF")
        }, character(1))

        sample_colors <- colors[as.character(dataset[[factor_col]])]
        names(sample_colors) <- dataset$Name
        labels_colors(dend) <- sample_colors[labels(dend)]

        colors_mat <- as.data.frame(dataset[, input$replicabilityFactors, drop = FALSE])
        for (f in input$replicabilityFactors) {
            colors_mat[[f]] <- vapply(colors_mat[[f]], function(l) {
                col_input <- input[[paste0("GroupColour", f, l)]]
                if (is.null(col_input)) {
                    return("#808080FF")
                }
                tryCatch(paste0(col2hex(col_input), "FF"), error = function(e) "#808080FF")
            }, character(1))
        }
        colors_mat <- colors_mat[match(labels(dend), dataset$Name), , drop = FALSE]
    }

    par(mar = c(20, 8, 6, 4))
    plot(dend,
        main = paste0("Cluster Dendrogram (", input$replicabilityDistMethod, "/", input$replicabilityClustMethod, "): ", title_suffix),
        xlab = "", sub = "", cex = input$replicabilityTextSize / 10
    )

    if (!is.null(input$replicabilityFactors)) {
        colored_bars(colors = colors_mat, dend = dend, sort_by_labels_order = FALSE)
    }
}


# =============================================================================
# Summary / methods-text helpers (from server-summary.R)
# =============================================================================

#' count_result_summary
#' SUSHI RNA-seq settings table for the Summary tab. A port of
#' `ezRun::makeCountResultSummary()`, so exploreDE does not need ezRun for it.
#' @param param SUSHI parameter list (`metadata(se)$param`).
#' @param se The SUSHI deResult `SummarizedExperiment`.
#' @return A one-column data.frame of settings.
#' @export
count_result_summary <- function(param, se) {
  specified <- function(x) !is.null(x) && length(x) > 0 && !is.na(x[1]) && nzchar(x[1])
  md <- S4Vectors::metadata(se)
  settings <- character()
  settings["Analysis:"] <- md$analysis
  settings["Reference:"] <- param$refBuild
  settings["Feature level:"] <- md$featureLevel
  settings["Data Column Used:"] <- md$countName
  settings["Method:"] <- md$method
  if (specified(param$sampleGroupBaseline) && specified(param$refGroupBaseline)) {
    settings["Baseline correction:"] <- paste(param$sampleGroupBaseline, "and",
                                              param$refGroupBaseline)
  }
  settings["Comparison:"] <- param$comparison
  if (!is.null(param$normMethod)) settings["Normalization:"] <- param$normMethod
  if (!is.null(param$deTest)) settings["Differential expression test:"] <- param$deTest
  if (isTRUE(param$useSigThresh)) {
    settings["Log2 signal threshold:"] <- signif(log2(param$sigThresh), digits = 4)
    settings["Linear signal threshold:"] <- signif(param$sigThresh, digits = 4)
  }
  as.data.frame(settings)
}

#' detectAligner
#' Inspect input_dataset.tsv to guess whether STAR / Kallisto / FeatureCounts was used.
#' @param dataDir Analysis directory containing ../input_dataset.tsv
#' @return Character: one of "STAR", "Kallisto", "FeatureCounts", "Unknown"
#' @export
detectAligner <- function(dataDir) {
  tryCatch({
    dataset_path <- file.path(dataDir, "../input_dataset.tsv")
    if (file.exists(dataset_path)) {
      dataset <- readr::read_tsv(dataset_path, show_col_types = FALSE)

      # Check Count [File] column for patterns
      count_col_name <- grep("Count.*File", names(dataset), value = TRUE, ignore.case = TRUE)
      if (length(count_col_name) > 0) {
        count_col <- dataset[[count_col_name[1]]][1]

        if (grepl("Kallisto", count_col, ignore.case = TRUE)) {
          return("Kallisto")
        } else if (grepl("FeatureCounts", count_col, ignore.case = TRUE)) {
          return("FeatureCounts")
        }
      }

      # Check if BAM column exists (indicates STAR was used)
      if (any(grepl("BAM.*File", names(dataset), ignore.case = TRUE))) {
        return("STAR")
      }
    }
    return("Unknown")
  }, error = function(e) {
    return("Unknown")
  })
}

#' detectDEMethod
#' Guess DE method (DESeq2 / edgeR) from dataDir string and param fields.
#' @param dataDir Analysis directory string
#' @param param Parameter list
#' @return Character: "DESeq2" or "edgeR"
#' @export
detectDEMethod <- function(dataDir, param) {
  # First check directory path - most reliable
  if (grepl("DESeq2", dataDir, ignore.case = TRUE)) {
    return("DESeq2")
  } else if (grepl("EdgeR|edgeR", dataDir, ignore.case = FALSE)) {
    return("edgeR")
  }

  # Check normalization method
  if (!is.null(param$normMethod)) {
    norm_method <- tolower(param$normMethod)
    if (grepl("tmm|rle", norm_method)) {
      return("edgeR")
    }
  }

  # Check param$testMethod if it exists
  if (!is.null(param$testMethod)) {
    # Extract first word before "using" if present
    test_method <- tolower(param$testMethod)
    test_method <- strsplit(test_method, " using ")[[1]][1]
    test_method <- trimws(test_method)

    # edgeR test methods
    if (test_method %in% c("glm", "lrt", "exacttest", "ql", "quasi-likelihood")) {
      return("edgeR")
    }
    # DESeq2 test methods
    if (test_method %in% c("wald")) {
      return("DESeq2")
    }
  }

  return("DESeq2")  # Default assumption
}

#' formatRefBuild
#' Pretty-print a reference build path like "Homo_sapiens/GENCODE/GRCh38.p13/.../Release_42-..."
#' @param ref_build Character reference build path
#' @return Human-readable string (e.g. "GENCODE human genome build GRCh38.p13 (Release 42)")
#' @export
formatRefBuild <- function(ref_build) {
  # If ref_build is already formatted nicely, return as-is
  if (!grepl("/", ref_build)) {
    return(ref_build)
  }

  # Parse path like "Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30"
  parts <- strsplit(ref_build, "/")[[1]]

  # Extract components
  species <- parts[1]
  database <- parts[2]
  build <- parts[3]

  # Format species
  species_formatted <- if (grepl("Homo_sapiens", species, ignore.case = TRUE)) {
    "human"
  } else if (grepl("Mus_musculus", species, ignore.case = TRUE)) {
    "mouse"
  } else {
    gsub("_", " ", species)
  }

  # Extract release from the last part if it exists
  release <- ""
  if (length(parts) >= 5) {
    release_part <- parts[5]
    # Extract release number from "Release_42-2023-01-30" or "Release_M31-2023-01-30"
    if (grepl("Release_", release_part)) {
      release_num <- sub("Release_([^-]+).*", "\\1", release_part)
      release <- paste0(" (Release ", release_num, ")")
    }
  }

  # Format final string
  formatted <- paste0(database, " ", species_formatted, " genome build ", build, release)
  return(formatted)
}

#' getFastpVersion
#' Scan preprocessing logs in `dirname(dataDir)` for the fastp version string.
#' @param dataDir Analysis directory
#' @return Version string or NULL
#' @export
getFastpVersion <- function(dataDir) {
  tryCatch({
    parent_dir <- dirname(dataDir)
    log_files <- list.files(parent_dir, pattern = "_preprocessing\\.log$", full.names = TRUE)

    if (length(log_files) > 0) {
      log_content <- readLines(log_files[1], n = 100, warn = FALSE)
      version_line <- grep("fastp v", log_content, value = TRUE)
      if (length(version_line) > 0) {
        version <- sub(".*fastp v([0-9.]+).*", "\\1", version_line[1])
        return(version)
      }
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

#' getPackageVersions
#' Look up installed versions of the given R packages (used for methods-text).
#' @param packages Character vector of package names
#' @return Named character vector of versions (or "version unavailable")
#' @export
getPackageVersions <- function(packages) {
  versions <- sapply(packages, function(pkg) {
    tryCatch({
      as.character(packageVersion(pkg))
    }, error = function(e) {
      "version unavailable"
    })
  })
  return(versions)
}

#' buildCitationsList
#' Build a named list of citation strings for methods section of the report.
#' @param aligner Aligner name from detectAligner()
#' @param de_method DE method name from detectDEMethod()
#' @param run_go Logical; include clusterProfiler citation
#' @param ref_build Reference build path (optional; drives GENCODE citation)
#' @return Named list of citation strings
#' @export
buildCitationsList <- function(aligner, de_method, run_go = FALSE, ref_build = NULL) {
  citations <- list()

  # Always cite SUSHI
  citations$SUSHI <- "Hatakeyama M, Opitz L, Russo G, Qi W, Schlapbach R, Rehrauer H (2016) SUSHI: an exquisite recipe for fully documented, reproducible and reusable NGS data analysis. BMC Bioinformatics 17, 228. https://doi.org/10.1186/s12859-016-1104-8"

  # FastQC
  citations$FastQC <- "Andrews S (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data. Available online at: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"

  # MultiQC
  citations$MultiQC <- "Ewels P, Magnusson M, Lundin S, Käller M (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354"

  # fastp
  citations$fastp <- "Chen S, Zhou Y, Chen Y, Gu J (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884–i890. https://doi.org/10.1093/bioinformatics/bty560"

  # Aligner-specific citations
  if (aligner == "STAR") {
    citations$STAR <- "Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013) STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21. https://doi.org/10.1093/bioinformatics/bts635"
  } else if (aligner == "Kallisto") {
    citations$Kallisto <- "Bray NL, Pimentel H, Melsted P, Pachter L (2016) Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology 34, 525–527. https://doi.org/10.1038/nbt.3519"
  } else if (aligner == "FeatureCounts") {
    citations$Rsubread <- "Liao Y, Smyth GK, Shi W (2013) The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research 41, e108. https://doi.org/10.1093/nar/gkt214"
    citations$STAR <- "Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013) STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21. https://doi.org/10.1093/bioinformatics/bts635"
  }

  # GENCODE citation if detected in reference build
  if (!is.null(ref_build) && grepl("GENCODE", ref_build, ignore.case = TRUE)) {
    citations$GENCODE <- "Frankish A, Diekhans M, Jungreis I, Lagarde J, Loveland JE, Mudge JM, Sisu C, Wright JC, Armstrong J, Barnes I, et al. (2021) GENCODE: reference annotation for the human and mouse genomes. Nucleic Acids Research 49, D916–D923. https://doi.org/10.1093/nar/gkaa1087"
  }

  # DE method citations
  if (de_method == "DESeq2") {
    citations$DESeq2 <- "Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550. https://doi.org/10.1186/s13059-014-0550-8"
  } else if (de_method == "edgeR") {
    citations$edgeR <- "Robinson MD, McCarthy DJ, Smyth GK (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139–140. https://doi.org/10.1093/bioinformatics/btp616"
  }

  # GO analysis citation if applicable
  if (run_go) {
    citations$clusterProfiler <- "Yu G, Wang LG, Han Y, He QY (2012) clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16, 284–287. https://doi.org/10.1089/omi.2011.0118"
  }

  # R and exploreDE citations
  citations$R <- "R Core Team (2020) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org"
  citations$exploreDE <- "Leary P, Rehrauer H (2024) exploreDE Interactive Shiny App. Zenodo. https://doi.org/10.5281/zenodo.13927692"

  return(citations)
}

#' generateMethodsText
#' Build the named list of methods-section paragraphs (trimming, mapping, DE analysis, GO, processing).
#' @param param Parameter list from the SE
#' @param aligner Aligner name from detectAligner()
#' @param de_method DE method name from detectDEMethod()
#' @param fastp_version Version string from getFastpVersion() (or NULL)
#' @param pkg_versions Named character vector from getPackageVersions()
#' @param se SummarizedExperiment (for library-size info)
#' @param raw_reads_text Optional prefix describing raw-read counts
#' @return Named list of methods paragraphs
#' @export
generateMethodsText <- function(param, aligner, de_method, fastp_version, pkg_versions, se, raw_reads_text = "") {
  sections <- list()

  # Final Mapped/Counted Reads (from SE)
  final_reads_info <- ""
  if (!is.null(se)) {
    tryCatch({
      lib_sizes <- colSums(assay(se, "counts"))
      min_size <- round(min(lib_sizes) / 1e6, 1)
      max_size <- round(max(lib_sizes) / 1e6, 1)
      final_reads_info <- paste0("Final library sizes ranged from ", min_size, " million to ", max_size, " million reads. ")
    }, error = function(e) {
      final_reads_info <- ""
    })
  }

  # Read trimming section
  fastp_ver_text <- if (!is.null(fastp_version)) paste0(" (version ", fastp_version, ")") else ""
  sections$trimming <- paste0(
    raw_reads_text,
    "Read quality was inspected using FastQC (Andrews, 2010) and MultiQC (Ewels et al., 2016), and sequencing adaptors removed using fastp",
    fastp_ver_text,
    " (Chen et al., 2018). ",
    "Quality trimming was performed to remove low-quality bases, with reads shorter than 18 bases discarded. ",
    "Poly-X tails were trimmed with a minimum length threshold of 10 bases."
  )

  # Mapping and quantification section
  ref_build_raw <- if (!is.null(param$refBuild)) param$refBuild else "the reference genome"
  ref_build <- formatRefBuild(ref_build_raw)
  feature_level <- if (!is.null(param$featureLevel)) param$featureLevel else "gene"
  transcript_types <- if (!is.null(param$transcriptTypes)) param$transcriptTypes else "protein-coding"

  if (aligner == "Kallisto") {
    sections$mapping <- paste0(
      "Transcript-level abundance was estimated using Kallisto ",
      "(Bray et al., 2016), which performs pseudo-alignment of reads to the ", ref_build, " reference transcriptome. ",
      "Kallisto uses an expectation-maximisation algorithm to quantify transcript abundances. ",
      "Only ", transcript_types, " transcripts were retained for downstream analysis."
    )
  } else if (aligner == "STAR") {
    sections$mapping <- paste0(
      "Alignment of the RNA-seq reads was performed using the STAR aligner ",
      "(Dobin et al., 2013) with the ", ref_build, " reference genome as the reference. ",
      "STAR employs a seed-and-extend approach optimised for RNA-seq data, ",
      "enabling accurate detection of splice junctions. ",
      "Only ", transcript_types, " genes were retained for downstream analysis."
    )
  } else if (aligner == "FeatureCounts") {
    sections$mapping <- paste0(
      "Alignment of the RNA-seq reads was performed using the STAR aligner ",
      "(Dobin et al., 2013) with the ", ref_build, " reference genome as the reference. ",
      "The counting of ", feature_level, "-level expression values was performed using the 'featureCounts' function ",
      "of the R package Rsubread (Liao et al., 2013). ",
      "Only ", transcript_types, " genes were retained for downstream analysis."
    )
  } else {
    sections$mapping <- paste0(
      "Reads were aligned to the ", ref_build, " reference genome and ",
      "quantified at the ", feature_level, " level using standard RNA-seq ",
      "alignment and quantification methods, ",
      "with only ", transcript_types, " genes retained for downstream analysis."
    )
  }

  # Append final library size info
  sections$mapping <- paste0(sections$mapping, " ", final_reads_info)

  # Filter, Normalisation and Differential Expression section

  deseq2_ver <- if (pkg_versions["DESeq2"] != "version unavailable") paste0(" (version ", pkg_versions["DESeq2"], ")") else ""
  edger_ver <- if (pkg_versions["edgeR"] != "version unavailable") paste0(" (version ", pkg_versions["edgeR"], ")") else ""

  if (de_method == "DESeq2") {
    sections$de_analysis <- paste0(
      "Count data were normalised using DESeq2's median-of-ratios method, which accounts for differences in library size ",
      "and RNA composition bias between samples. ",
      "Differential expression was performed using the DESeq2 Bioconductor R package", deseq2_ver, " ",
      "(Love et al., 2014), which implements a generalised linear model with a negative binomial distribution ",
      "to account for technical and biological variability. ",
      "Dispersion and log2 fold change estimates were moderated using shrinkage methods to improve stability. ",
      "Statistical significance was assessed using the Wald test, with p-values adjusted for multiple testing ",
      "using the Benjamini-Hochberg method to control the false discovery rate (FDR)."
    )
  } else if (de_method == "edgeR") {
    sections$de_analysis <- paste0(
      "Count data were normalised using the trimmed mean of M-values (TMM) method implemented in edgeR, ",
      "which corrects for differences in library size and RNA composition between samples. ",
      "Differential expression was performed using the edgeR Bioconductor R package", edger_ver, " ",
      "(Robinson et al., 2010), which implements a generalised linear model with a negative binomial distribution ",
      "to account for technical and biological variability. ",
      "Dispersion estimates were moderated using empirical Bayes shrinkage to improve reliability for genes with low counts. ",
      "Statistical significance was assessed using the likelihood ratio test (LRT), ",
      "with p-values adjusted for multiple testing using the Benjamini-Hochberg method to control the false discovery rate (FDR)."
    )
  }

  # Gene ontology section (if applicable)
  if (!is.null(param$runGO) && param$runGO) {
    cluster_ver <- if (pkg_versions["clusterProfiler"] != "version unavailable") paste0(" (version ", pkg_versions["clusterProfiler"], ")") else ""
    fdr_ora <- if (!is.null(param$fdrThreshORA)) param$fdrThreshORA else 0.05
    fdr_gsea <- if (!is.null(param$fdrThreshGSEA)) param$fdrThreshGSEA else 0.05

    sections$go_analysis <- paste0(
      "Gene Ontology (GO) term pathway analysis was performed using both the hypergeometric over-representation test ",
      "via the 'enricher' function, and gene-set enrichment analysis via the 'GSEA' function, ",
      "of the clusterProfiler Bioconductor R package", cluster_ver, " ",
      "(Yu et al., 2012), ",
      "with FDR thresholds of ", fdr_ora, " and ", fdr_gsea, ", respectively."
    )
  }

  # Data processing section
  # Get R version info if available
  r_version <- paste0("R version ", R.version$major, ".", R.version$minor)

  # Get Bioconductor version if available
  bioc_version <- ""
  tryCatch({
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      bv <- as.character(BiocManager::version())
      bioc_version <- paste0(" and Bioconductor version ", bv)
    }
  }, error = function(e) {
    bioc_version <- ""
  })

  sections$data_processing <- paste0(
    "All analyses were conducted using the SUSHI framework (Hatakeyama et al., 2016), ",
    "which provides a reproducible and documented workflow for NGS data analysis. ",
    "All R functions were executed on ", r_version, bioc_version, " (R Core Team, 2020). ",
    "Figures were generated using the exploreDE Shiny application using the data generated from the above steps (Leary and Rehrauer, 2024)."
  )

  return(sections)
}


# =============================================================================
# Ordination helpers (from server-ordination.R)
# NOTE: plot_ordination_samples and plot_ordination_biplot take a Shiny
# `input` object; step_batch_correction calls showNotification(). These
# will error if called outside a Shiny session.
# =============================================================================

#' calculate_pca
#' Run base R prcomp on a sample-by-feature matrix and return scores / loadings / variance.
#' @param counts Matrix with samples as rows, features as columns
#' @param center Logical
#' @param scale Logical
#' @return List with scores, loadings, variance, method, component_names, removed_samples
#' @export
calculate_pca <- function(counts, center = TRUE, scale = FALSE) {
  # Remove samples with missing values
  completeCases <- complete.cases(counts)
  counts_clean <- counts[completeCases, ]

  if (nrow(counts_clean) < 3 || ncol(counts_clean) < 3) {
    return(NULL)
  }

  # Run PCA
  pca_result <- tryCatch(
    {
      prcomp(counts_clean, center = center, scale. = scale)
    },
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(pca_result)) {
    return(NULL)
  }

  # Calculate variance explained
  variance <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

  list(
    scores = tibble::as_tibble(pca_result$x, rownames = "Name"),
    loadings = tibble::as_tibble(pca_result$rotation, rownames = "gene"),
    variance = variance,
    method = "PCA",
    component_names = paste0("PC", 1:length(variance)),
    removed_samples = sum(!completeCases)
  )
}

#' calculate_pcoa
#' Run PCoA (via cmdscale) or NMDS (via vegan::metaMDS) on a distance matrix.
#' @param counts Matrix with samples as rows, features as columns
#' @param method Distance method passed to vegan::vegdist
#' @param nmds Logical, run NMDS instead of PCoA
#' @return List with scores, loadings (NULL), variance / stress, method, component_names
#' @export
calculate_pcoa <- function(counts, method = "euclidean", nmds = FALSE) {
  # Remove samples with missing values
  completeCases <- complete.cases(counts)
  counts_clean <- counts[completeCases, ]

  if (nrow(counts_clean) < 3 || ncol(counts_clean) < 3) {
    return(NULL)
  }

  # Calculate distance matrix
  dist_mat <- tryCatch(
    {
      vegan::vegdist(counts_clean, method = method)
    },
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(dist_mat)) {
    return(NULL)
  }

  if (nmds) {
    # NMDS
    nmds_result <- tryCatch(
      {
        vegan::metaMDS(dist_mat, k = 2, trymax = 100, trace = 0)
      },
      error = function(e) {
        return(NULL)
      }
    )

    if (is.null(nmds_result)) {
      return(NULL)
    }

    scores_mat <- vegan::scores(nmds_result)

    # Name the columns before converting to tibble
    component_names <- paste0("NMDS", 1:ncol(scores_mat))
    colnames(scores_mat) <- component_names

    list(
      scores = tibble::as_tibble(scores_mat, rownames = "Name"),
      loadings = NULL,
      variance = NULL,
      stress = nmds_result$stress,
      method = "NMDS",
      component_names = component_names,
      removed_samples = sum(!completeCases)
    )
  } else {
    # Classical PCoA
    pcoa_result <- tryCatch(
      {
        cmdscale(dist_mat, k = min(10, nrow(counts_clean) - 1), eig = TRUE)
      },
      error = function(e) {
        return(NULL)
      }
    )

    if (is.null(pcoa_result)) {
      return(NULL)
    }

    # Name the columns before converting to tibble
    n_dims <- ncol(pcoa_result$points)
    component_names <- paste0("PCo", 1:n_dims)
    colnames(pcoa_result$points) <- component_names

    # Calculate variance explained - only for the dimensions we kept
    # cmdscale returns all eigenvalues, but we only use the first n_dims
    variance_explained <- (pcoa_result$eig[1:n_dims] / sum(abs(pcoa_result$eig))) * 100

    list(
      scores = tibble::as_tibble(pcoa_result$points, rownames = "Name"),
      loadings = NULL,
      variance = variance_explained,
      stress = NULL,
      method = "PCoA",
      component_names = component_names,
      removed_samples = sum(!completeCases)
    )
  }
}

#' calculate_glmpca
#' Run GLM-PCA on raw counts using the glmpca package.
#' @param counts Matrix with features as rows, samples as columns (RAW counts)
#' @param L Number of latent dimensions
#' @param family One of "poi", "nb", "binom"
#' @param X Optional covariate matrix
#' @return List with scores, loadings (NULL), variance (NULL), method, component_names
#' @export
calculate_glmpca <- function(counts, L = 10, family = "nb", X = NULL) {
  if (nrow(counts) < 3 || ncol(counts) < 3) {
    return(NULL)
  }

  # Validate family parameter
  if (is.null(family) || !is.character(family) || length(family) != 1) {
    warning("GLM-PCA: Invalid family parameter (must be a single character string)")
    return(NULL)
  }

  valid_families <- c("poi", "nb", "binom")
  if (!family %in% valid_families) {
    warning("GLM-PCA: family must be one of: poi, nb, binom. Got: ", family)
    return(NULL)
  }

  # Remove genes with all zeros (required by glmpca)
  row_sums <- rowSums(counts)
  non_zero_rows <- row_sums > 0
  if (sum(non_zero_rows) < 3) {
    warning("GLM-PCA: Too few non-zero genes remaining")
    return(NULL)
  }
  counts <- counts[non_zero_rows, ]

  # GLM-PCA needs raw/normalized counts (not log-transformed)
  glmpca_result <- tryCatch(
    {
      glmpca::glmpca(counts, L = L, fam = family, X = X)
    },
    error = function(e) {
      warning("GLM-PCA error: ", e$message)
      return(NULL)
    }
  )

  if (is.null(glmpca_result)) {
    return(NULL)
  }

  # Name the columns before converting to tibble
  component_names <- paste0("dim", 1:L)
  colnames(glmpca_result$factors) <- component_names

  list(
    scores = tibble::as_tibble(glmpca_result$factors, rownames = "Name"),
    loadings = NULL, # GLM-PCA loadings are on a different scale, biplot not supported
    variance = NULL, # GLM-PCA doesn't have simple variance explained like PCA
    method = "GLM-PCA",
    component_names = component_names,
    removed_samples = 0 # GLM-PCA handles missing values internally
  )
}

#' prepare_biplot_data
#' Pick which features to overlay as arrows on a biplot based on selection method.
#' @param ordination_result Result from calculate_* function
#' @param top_n Number of top features
#' @param custom_genes Manual selection
#' @param pc_x PC index or name for X
#' @param pc_y PC index or name for Y
#' @param selection_method One of "overall", "per_pc", "manual"
#' @return List with sample_data, gene_data, pc_x_name, pc_y_name
#' @export
prepare_biplot_data <- function(ordination_result, top_n = 10,
                                custom_genes = NULL, pc_x = 1, pc_y = 2,
                                selection_method = "overall") {
  if (is.null(ordination_result$loadings)) {
    # Return NULL gracefully instead of stopping - biplot not available for distance-based methods
    return(NULL)
  }

  loadings <- ordination_result$loadings

  # Get component names (pc_x and pc_y can be either indices or component name strings)
  pc_x_name <- if (is.numeric(pc_x)) ordination_result$component_names[pc_x] else pc_x
  pc_y_name <- if (is.numeric(pc_y)) ordination_result$component_names[pc_y] else pc_y

  # Calculate loading magnitude
  loadings_subset <- loadings %>%
    dplyr::mutate(magnitude = sqrt(.data[[pc_x_name]]^2 + .data[[pc_y_name]]^2))

  # Select genes based on selection method
  if (selection_method == "manual") {
    # Manual selection: use custom_genes
    if (!is.null(custom_genes) && length(custom_genes) > 0) {
      genes_to_plot <- loadings_subset %>%
        dplyr::filter(gene %in% custom_genes)
    } else {
      genes_to_plot <- loadings_subset[0, ] # Empty tibble if no genes specified
    }
  } else if (selection_method == "per_pc") {
    # Per-PC selection: top N from each PC separately
    top_x <- loadings_subset %>%
      dplyr::arrange(dplyr::desc(abs(.data[[pc_x_name]]))) %>%
      head(top_n) %>%
      dplyr::mutate(pc_source = pc_x_name)

    top_y <- loadings_subset %>%
      dplyr::arrange(dplyr::desc(abs(.data[[pc_y_name]]))) %>%
      head(top_n) %>%
      dplyr::mutate(pc_source = pc_y_name)

    # Combine and deduplicate - keep first occurrence (prioritize PC X if gene is in both)
    genes_to_plot <- dplyr::bind_rows(top_x, top_y) %>%
      dplyr::distinct(gene, .keep_all = TRUE)
  } else {
    # Overall selection: top N by combined magnitude (default behavior)
    genes_to_plot <- loadings_subset %>%
      dplyr::arrange(dplyr::desc(magnitude)) %>%
      head(top_n)
  }

  list(
    sample_data = ordination_result$scores,
    gene_data = genes_to_plot,
    pc_x_name = pc_x_name,
    pc_y_name = pc_y_name
  )
}

#' step_group_subset
#' Cascade step: subset counts and dataset to the specified groups.
#' @export
step_group_subset <- function(counts, dataset, factor_col, groups) {
  dataset <- dataset[dataset[[factor_col]] %in% groups, ]
  counts <- counts[, rownames(dataset), drop = FALSE]
  list(counts = counts, dataset = dataset)
}

#' step_peptide_filter
#' Cascade step: keep features with >= min_peptides entries (proteomics only).
#' @export
step_peptide_filter <- function(counts, dataset, seqAnnoList, min_peptides) {
  if (!is.null(seqAnnoList) && !is.null(min_peptides)) {
    keep <- seqAnnoList[[1]]$feature_id[seqAnnoList[[1]]$nrPeptides >= min_peptides]
    counts <- counts[rownames(counts) %in% keep, , drop = FALSE]
  }
  list(counts = counts, dataset = dataset)
}

#' step_log2
#' Cascade step: log2(counts + 1) for count types that are not already log-scale.
#' @export
step_log2 <- function(counts, dataset, counts_type) {
  if (counts_type %in% c("TPM", "FPKM", "Raw", "Normalised")) {
    counts <- log2(counts + 1)
  }
  list(counts = counts, dataset = dataset)
}

#' step_batch_correction
#' Cascade step: apply limma::removeBatchEffect via apply_batch_correction().
#' Calls showNotification() on failure / NA batch; never throws.
#' @export
step_batch_correction <- function(counts, dataset, batch_column, design_column = "None") {
  if (!is.null(batch_column) && nzchar(batch_column) && batch_column != "None") {
    if (any(is.na(dataset[[batch_column]]))) {
      showNotification("Batch correction skipped: batch column contains NA values",
        type = "warning", duration = 8)
    } else {
      counts <- apply_batch_correction(
        counts    = counts,
        dataset   = dataset,
        batchCol  = batch_column,
        designCol = design_column
      )
    }
  }
  list(counts = counts, dataset = dataset)
}

#' step_top_n_sd
#' Cascade step: keep top-N features by row-wise SD.
#' @export
step_top_n_sd <- function(counts, dataset, top_n) {
  top_n <- min(top_n, nrow(counts))
  if (nrow(counts) > top_n) {
    sds <- apply(counts, 1, sd, na.rm = FALSE)
    top_features <- names(head(sort(sds, decreasing = TRUE), n = top_n))
    counts <- counts[top_features, , drop = FALSE]
  }
  list(counts = counts, dataset = dataset)
}

#' plot_ordination_samples
#' Render the sample-only ordination scatter plot (no gene arrows).
#' @export
plot_ordination_samples <- function(ordination_result, pc_scores, input, coloursPCA, levelLabels = NULL) {
  fillLabels  <- if (is.null(levelLabels$fill)) ggplot2::waiver() else levelLabels$fill
  shapeLabels <- if (is.null(levelLabels$shape)) ggplot2::waiver() else levelLabels$shape
  pc_x <- input$pcaX
  pc_y <- input$pcaY

  # Base plot
  p <- ggplot(pc_scores, aes(x = .data[[pc_x]], y = .data[[pc_y]], fill = .data[[input$pcaFactor1]]))

  if (input$pcaAddEllipses) {
    p <- p + stat_ellipse(aes(colour = .data[[input$pcaFactor1]], fill = .data[[input$pcaFactor1]]),
      geom = "polygon", alpha = input$pcaEllipsesAlpha, show.legend = FALSE
    )
    p <- p + scale_colour_manual(breaks = names(coloursPCA), values = as.character(coloursPCA))
  }

  if (input$pcaFactor2 == "None") {
    p <- p + geom_point(size = input$dotSizePCA, alpha = input$alphaPCA, shape = 21, stroke = input$dotBorderPCA)
  } else {
    p <- p + geom_point(aes(shape = .data[[input$pcaFactor2]]),
      size = input$dotSizePCA,
      alpha = input$alphaPCA, stroke = input$dotBorderPCA
    )
    p <- p + scale_shape_manual(values = c(21, 22, 23, 24, 25, 8, 3, 4), labels = shapeLabels)
  }

  if (input$pcaShowNames) {
    p <- p + geom_label_repel(aes(label = Name, fill = .data[[input$pcaFactor1]]),
      force = 2, nudge_y = as.numeric(input$geneLabelNudgePCAY) / 10,
      nudge_x = as.numeric(input$geneLabelNudgePCAX) / 10,
      size = (input$textSizePCA / 3), max.overlaps = input$pcaLabelMaxOverlap,
      fontface = "bold", show.legend = FALSE, na.rm = TRUE
    )
  }

  # Axis labels
  if (!is.null(ordination_result$variance)) {
    pc_x_num <- as.numeric(gsub("\\D", "", pc_x))
    pc_y_num <- as.numeric(gsub("\\D", "", pc_y))
    x_lab <- paste0(pc_x, ": ", round(ordination_result$variance[pc_x_num]), "% variance")
    y_lab <- paste0(pc_y, ": ", round(ordination_result$variance[pc_y_num]), "% variance")
  } else if (!is.null(ordination_result$stress)) {
    x_lab <- paste0(pc_x, " (stress: ", round(ordination_result$stress, 3), ")")
    y_lab <- pc_y
  } else {
    x_lab <- pc_x
    y_lab <- pc_y
  }

  p <- p + labs(x = x_lab, y = y_lab, fill = input$pcaFactor1)

  if (input$pcaAxesProp && !is.null(ordination_result$variance)) {
    pc_x_num <- as.numeric(gsub("\\D", "", pc_x))
    pc_y_num <- as.numeric(gsub("\\D", "", pc_y))
    p <- p + coord_fixed(ratio = ordination_result$variance[pc_y_num] / ordination_result$variance[pc_x_num])
  }

  p <- p + scale_fill_manual(breaks = names(coloursPCA), values = as.character(coloursPCA), labels = fillLabels)
  p <- p + theme_bw(base_size = input$textSizePCA)
  face <- switch(input$boldPCA,
    "TRUE" = "bold",
    "FALSE" = "plain"
  )
  p <- p + theme(axis.text = element_text(face = face), axis.title = element_text(face = face))
  if (isTRUE(levelLabels$rich)) p <- p + rich_label_theme(legend = TRUE)

  # Conditionally hide grid lines
  if (!input$showLinesPCA) {
    p <- p + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  }

  # Conditionally hide axes
  if (!input$showAxesPCA) {
    p <- p + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank()
    )
  }

  return(p)
}

#' plot_ordination_biplot
#' Render the ordination biplot (samples + feature arrows).
#' @export
plot_ordination_biplot <- function(ordination_result, pc_scores, input, coloursPCA, biplot_data, arrow_scale = 1, arrow_colour = "red", selection_method = "overall", levelLabels = NULL) {
  fillLabels <- if (is.null(levelLabels$fill)) ggplot2::waiver() else levelLabels$fill
  pc_x <- biplot_data$pc_x_name
  pc_y <- biplot_data$pc_y_name

  # Scale arrows
  sample_range <- max(abs(c(biplot_data$sample_data[[pc_x]], biplot_data$sample_data[[pc_y]])), na.rm = TRUE)
  gene_range <- max(abs(c(biplot_data$gene_data[[pc_x]], biplot_data$gene_data[[pc_y]])), na.rm = TRUE)
  scale_factor <- (sample_range / gene_range) * arrow_scale

  biplot_data$gene_data <- biplot_data$gene_data %>%
    dplyr::mutate(pc_x_scaled = .data[[pc_x]] * scale_factor, pc_y_scaled = .data[[pc_y]] * scale_factor)

  # Add arrow color column for per_pc mode
  if (selection_method == "per_pc" && "pc_source" %in% names(biplot_data$gene_data)) {
    biplot_data$gene_data <- biplot_data$gene_data %>%
      dplyr::mutate(arrow_color = ifelse(pc_source == pc_x, arrow_colour$x, arrow_colour$y))
  } else {
    # Single color for all arrows
    biplot_data$gene_data$arrow_color <- if (is.list(arrow_colour)) arrow_colour$x else arrow_colour
  }

  p <- ggplot()
  p <- p + geom_point(
    data = pc_scores, aes(x = .data[[pc_x]], y = .data[[pc_y]], fill = .data[[input$pcaFactor1]]),
    size = input$dotSizePCA, alpha = input$alphaPCA, shape = 21, stroke = input$dotBorderPCA
  )
  p <- p + geom_segment(
    data = biplot_data$gene_data, aes(x = 0, y = 0, xend = pc_x_scaled, yend = pc_y_scaled, colour = arrow_color),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 0.5, alpha = 0.7
  )
  p <- p + geom_text_repel(
    data = biplot_data$gene_data, aes(x = pc_x_scaled, y = pc_y_scaled, label = gene, colour = arrow_color),
    size = input$textSizePCA / 4, fontface = "bold", max.overlaps = 20
  )
  p <- p + scale_colour_identity() # Use the actual color values

  # Axes
  if (!is.null(ordination_result$variance)) {
    pc_x_num <- as.numeric(gsub("\\D", "", pc_x))
    pc_y_num <- as.numeric(gsub("\\D", "", pc_y))
    x_lab <- paste0(pc_x, ": ", round(ordination_result$variance[pc_x_num]), "% variance")
    y_lab <- paste0(pc_y, ": ", round(ordination_result$variance[pc_y_num]), "% variance")
  } else {
    x_lab <- pc_x
    y_lab <- pc_y
  }

  p <- p + labs(x = x_lab, y = y_lab, fill = input$pcaFactor1)
  p <- p + scale_fill_manual(breaks = names(coloursPCA), values = as.character(coloursPCA), labels = fillLabels)
  p <- p + theme_bw(base_size = input$textSizePCA)
  face <- switch(input$boldPCA,
    "TRUE" = "bold",
    "FALSE" = "plain"
  )
  p <- p + theme(axis.text = element_text(face = face), axis.title = element_text(face = face))
  if (isTRUE(levelLabels$rich)) p <- p + rich_label_theme(legend = TRUE)

  # Conditionally hide grid lines
  if (!input$showLinesPCA) {
    p <- p + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  }

  # Conditionally hide axes
  if (!input$showAxesPCA) {
    p <- p + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank()
    )
  }

  return(p)
}

#' make_qc_scatter
#' Build a QC scatter plot in the style of ezRun's ezXYScatter.2:
#' log10 axes, three typed colour layers, an optional 1:1 solid guide line,
#' and optional N-fold dashed guide lines. Used by the DE QC and Count QC
#' tabs for between-group and intra-group scatter plots.
#' @param x Numeric vector (linear scale, NOT log2)
#' @param y Numeric vector (linear scale, same length as x)
#' @param types_vec Character vector of category labels, same length as x/y
#' @param types_colours Named character vector mapping level -> colour
#' @param xlab,ylab Axis labels
#' @param title Plot title (NULL for none)
#' @param point_size,alpha Point aesthetics
#' @param base_size Font size passed to ggprism::theme_prism
#' @param draw_order Optional character vector of level names defining z-order
#'   (last level = drawn on top). Defaults to names(types_colours).
#' @param guide_colour Colour for the 1:1 and fold-change guide lines.
#'   Default "blue".
#' @param guide_fold Fold-change threshold for the dashed guide lines.
#'   Default 2 (gives ±log10(2) dashed lines).
#' @param show_identity Whether to draw the solid 1:1 identity line.
#'   Default TRUE.
#' @param show_fold_bands Whether to draw the dashed fold-change guide lines.
#'   Default TRUE.
#' @return ggplot object
#' @export
make_qc_scatter <- function(x, y, types_vec, types_colours,
                            xlab, ylab, title = NULL,
                            point_size = 1, alpha = 0.6,
                            base_size = 12, draw_order = NULL,
                            guide_colour = "blue",
                            guide_fold = 2,
                            show_identity = TRUE,
                            show_fold_bands = TRUE,
                            log_scale = TRUE) {
  # log_scale = TRUE: x/y are LINEAR values (e.g. RNA-seq group-mean counts)
  #   plotted on log10 axes; a `guide_fold`-fold band is log10(guide_fold) apart.
  # log_scale = FALSE: x/y are ALREADY in log space (e.g. proteomics log2
  #   intensities, centred near 0, can be negative) — plot on linear axes and
  #   space the fold band by log2(guide_fold). Applying log10 to already-log
  #   data drops negatives and mislocates the identity/fold guides.
  if (is.null(draw_order)) draw_order <- names(types_colours)
  df <- tibble::tibble(x = x, y = y, types = types_vec)
  df$types <- factor(df$types, levels = draw_order)
  df <- df[order(df$types), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = types)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::scale_color_manual(values = types_colours)
  if (isTRUE(log_scale)) {
    p <- p + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
  }

  if (isTRUE(show_identity)) {
    p <- p + ggplot2::geom_abline(intercept = 0, colour = guide_colour)
  }
  if (isTRUE(show_fold_bands) && is.numeric(guide_fold) && guide_fold > 1) {
    off <- if (isTRUE(log_scale)) log10(guide_fold) else log2(guide_fold)
    p <- p +
      ggplot2::geom_abline(intercept =  off, colour = guide_colour, linetype = "dashed") +
      ggplot2::geom_abline(intercept = -off, colour = guide_colour, linetype = "dashed")
  }

  p +
    ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )) +
    ggplot2::labs(x = xlab, y = ylab, title = title)
}

#' build_qc_scatter_types
#' Build the `types_vec`, `types_colours` and `draw_order` used by
#' make_qc_scatter() for the three shared colour-by modes (Standard /
#' GC / Width). Used by both the DE QC between-group scatter and the
#' Count QC intra-group scatter.
#'
#' @param rd Row-data data frame aligned with the x/y vectors (must have
#'   columns `feature_id`, optionally `usedInTest`, `log2Ratio`, `pValue`,
#'   `gc`, `featWidth`).
#' @param mode One of "Standard", "GC content", "Gene width".
#' @param present_mask Logical vector, same length as nrow(rd): TRUE means
#'   the feature is "Present" for Standard colouring. Caller resolves this
#'   from dataType (RNASeq: usedInTest; CountQC: presentFlag mean).
#' @param p_thresh,lfc_thresh Significance thresholds for Standard mode.
#'   Only applied when `rd` has `pValue` and `log2Ratio`.
#' @param gc_lo,gc_hi,w_lo,w_hi Bin thresholds for GC / Width modes.
#' @param std_cols Named list with `Absent`, `Present`, `Significants` colours.
#' @param gc_cols Named list with `low`, `mid`, `high` colours for GC mode.
#' @param width_cols Named list with `low`, `mid`, `high` colours for Width mode.
#' @return A list with `types` (character vector), `colours` (named char),
#'   and `draw_order` (level order for z-ordering).
#' @export
build_qc_scatter_types <- function(rd, mode, present_mask,
                                   p_thresh = 0.05, lfc_thresh = 1,
                                   gc_lo = 0.4, gc_hi = 0.6,
                                   w_lo = 500, w_hi = 5000,
                                   std_cols = list(Absent = "grey70",
                                                   Present = "black",
                                                   Significants = "#A6CEE3"),
                                   gc_cols = list(low = "#1F78B4",
                                                  mid = "black",
                                                  high = "#FF7F00"),
                                   width_cols = list(low = "#B15928",
                                                     mid = "black",
                                                     high = "#33A02C")) {
  n <- nrow(rd)
  if (is.null(present_mask)) present_mask <- rep(TRUE, n)

  if (mode == "Standard") {
    types <- ifelse(present_mask, "Present", "Absent")
    if ("log2Ratio" %in% colnames(rd) && "pValue" %in% colnames(rd)) {
      is_sig <- !is.na(rd$pValue) & rd$pValue < p_thresh &
                !is.na(rd$log2Ratio) & abs(rd$log2Ratio) >= lfc_thresh
      types[is_sig] <- "Significants"
    }
    colours <- c(Absent = std_cols$Absent,
                 Present = std_cols$Present,
                 Significants = std_cols$Significants)
    return(list(types = types, colours = colours,
                draw_order = c("Absent", "Present", "Significants")))
  }

  if (mode == "GC content") {
    mid_lab <- sprintf("%.2f-%.2f", gc_lo, gc_hi)
    lo_lab  <- sprintf("GC < %.2f", gc_lo)
    hi_lab  <- sprintf("GC > %.2f", gc_hi)
    types <- rep(mid_lab, n)
    types[!is.na(rd$gc) & rd$gc < gc_lo] <- lo_lab
    types[!is.na(rd$gc) & rd$gc > gc_hi] <- hi_lab
    colours <- stats::setNames(
      c(gc_cols$low, gc_cols$mid, gc_cols$high),
      c(lo_lab, mid_lab, hi_lab)
    )
    return(list(types = types, colours = colours,
                draw_order = c(mid_lab, lo_lab, hi_lab)))
  }

  if (mode == "Gene width") {
    mid_lab <- sprintf("%d-%d nt", as.integer(w_lo), as.integer(w_hi))
    lo_lab  <- sprintf("< %d nt", as.integer(w_lo))
    hi_lab  <- sprintf("> %d nt", as.integer(w_hi))
    types <- rep(mid_lab, n)
    types[!is.na(rd$featWidth) & rd$featWidth < w_lo] <- lo_lab
    types[!is.na(rd$featWidth) & rd$featWidth > w_hi] <- hi_lab
    colours <- stats::setNames(
      c(width_cols$low, width_cols$mid, width_cols$high),
      c(lo_lab, mid_lab, hi_lab)
    )
    return(list(types = types, colours = colours,
                draw_order = c(mid_lab, lo_lab, hi_lab)))
  }

  list(types = rep("All", n), colours = c(All = "black"), draw_order = "All")
}

#' transform_feature_id
#' Mirror the character substitutions that convert_genomics_se applies to
#' countList rownames for isoform-level data (replacing ":", "/", "(", ")",
#' "-", and spaces with "."). Gene-level feature_ids are unaffected. Use
#' this before matching seqAnno$feature_id against signal rownames.
#' @param x Character vector of feature_ids
#' @return Character vector with init-compatible substitutions applied
#' @export
transform_feature_id <- function(x) {
  gsub(" ", ".", gsub("[:/()-]", ".", x), fixed = FALSE)
}

#' align_signal_to_seqanno
#' Drop-through alignment of a count/signal matrix to a seqAnno data frame
#' by feature_id. Falls back to the isoform-style transformation when the
#' direct match leaves NAs. Returns the matrix reordered to the seqAnno rows
#' that successfully matched, along with the subset of seqAnno kept.
#' @param signal Count/signal matrix with rownames to align
#' @param rd A data frame with a `feature_id` column
#' @return list(signal = aligned matrix, rd = subset of rd that matched)
#' @export
align_signal_to_seqanno <- function(signal, rd) {
  align <- match(rd$feature_id, rownames(signal))
  if (any(is.na(align))) {
    alt <- transform_feature_id(rd$feature_id)
    align_alt <- match(alt, rownames(signal))
    if (sum(!is.na(align_alt)) > sum(!is.na(align))) align <- align_alt
  }
  keep <- !is.na(align)
  list(
    signal = signal[align[keep], , drop = FALSE],
    rd     = rd[keep, , drop = FALSE]
  )
}
