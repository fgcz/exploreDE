# ============================================
# HELPER FUNCTIONS FOR ORDINATION ANALYSIS
# ============================================

#' Calculate PCA
#' @param counts Matrix with samples as rows, genes as columns (transposed from original)
#' @param center Logical, center the data?
#' @param scale Logical, scale the data?
#' @return List with scores, loadings, variance, method name, component names
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
    scores = as_tibble(pca_result$x, rownames = "Name"),
    loadings = as_tibble(pca_result$rotation, rownames = "gene"),
    variance = variance,
    method = "PCA",
    component_names = paste0("PC", 1:length(variance)),
    removed_samples = sum(!completeCases)
  )
}

#' Calculate PCoA/MDS
#' @param counts Matrix with samples as rows, genes as columns
#' @param method Distance method
#' @param nmds Use NMDS instead of PCoA?
#' @return List with scores, loadings (NULL), variance/stress, method name
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
      scores = as_tibble(scores_mat, rownames = "Name"),
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
      scores = as_tibble(pcoa_result$points, rownames = "Name"),
      loadings = NULL,
      variance = variance_explained,
      stress = NULL,
      method = "PCoA",
      component_names = component_names,
      removed_samples = sum(!completeCases)
    )
  }
}

#' Calculate GLM-PCA
#' @param counts Matrix with genes as rows, samples as columns (RAW COUNTS)
#' @param L Number of dimensions
#' @param family Family (nb, poi, binom)
#' @return List with scores, loadings, variance (deviance explained)
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
    scores = as_tibble(glmpca_result$factors, rownames = "Name"),
    loadings = NULL, # GLM-PCA loadings are on a different scale, biplot not supported
    variance = NULL, # GLM-PCA doesn't have simple variance explained like PCA
    method = "GLM-PCA",
    component_names = component_names,
    removed_samples = 0 # GLM-PCA handles missing values internally
  )
}

#' Prepare biplot data
#' @param ordination_result Result from calculate_* function
#' @param top_n Number of top genes to include (for overall and per_pc methods)
#' @param custom_genes Vector of specific genes to include (for manual method)
#' @param pc_x X axis component number or name
#' @param pc_y Y axis component number or name
#' @param selection_method Gene selection method: "overall", "per_pc", or "manual"
#' @return List with sample_data and gene_data for plotting
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
        filter(gene %in% custom_genes)
    } else {
      genes_to_plot <- loadings_subset[0, ] # Empty tibble if no genes specified
    }
  } else if (selection_method == "per_pc") {
    # Per-PC selection: top N from each PC separately
    top_x <- loadings_subset %>%
      dplyr::arrange(desc(abs(.data[[pc_x_name]]))) %>%
      head(top_n) %>%
      dplyr::mutate(pc_source = pc_x_name)

    top_y <- loadings_subset %>%
      dplyr::arrange(desc(abs(.data[[pc_y_name]]))) %>%
      head(top_n) %>%
      dplyr::mutate(pc_source = pc_y_name)

    # Combine and deduplicate - keep first occurrence (prioritize PC X if gene is in both)
    genes_to_plot <- bind_rows(top_x, top_y) %>%
      distinct(gene, .keep_all = TRUE)
  } else {
    # Overall selection: top N by combined magnitude (default behavior)
    genes_to_plot <- loadings_subset %>%
      dplyr::arrange(desc(magnitude)) %>%
      head(top_n)
  }

  list(
    sample_data = ordination_result$scores,
    gene_data = genes_to_plot,
    pc_x_name = pc_x_name,
    pc_y_name = pc_y_name
  )
}

#' Preprocessing cascade step functions
#' Each takes (counts, dataset, ...) and returns list(counts, dataset)
#' so they can be composed in any order.

step_group_subset <- function(counts, dataset, factor_col, groups) {
  dataset <- dataset[dataset[[factor_col]] %in% groups, ]
  counts <- counts[, rownames(dataset), drop = FALSE]
  list(counts = counts, dataset = dataset)
}

step_peptide_filter <- function(counts, dataset, seqAnnoList, min_peptides) {
  if (!is.null(seqAnnoList) && !is.null(min_peptides)) {
    keep <- seqAnnoList[[1]]$feature_id[seqAnnoList[[1]]$nrPeptides >= min_peptides]
    counts <- counts[rownames(counts) %in% keep, , drop = FALSE]
  }
  list(counts = counts, dataset = dataset)
}

step_log2 <- function(counts, dataset, counts_type) {
  if (counts_type %in% c("TPM", "FPKM", "Raw", "Normalised")) {
    counts <- log2(counts + 1)
  }
  list(counts = counts, dataset = dataset)
}

step_batch_correction <- function(counts, dataset, batch_column) {
  if (!is.null(batch_column)) {
    counts <- tryCatch({
      if (any(is.na(dataset[[batch_column]]))) {
        showNotification("Batch correction skipped: batch column contains NA values",
          type = "warning", duration = 3)
        counts
      } else {
        limma::removeBatchEffect(counts, dataset[[batch_column]])
      }
    }, error = function(e) {
      showNotification(paste("Batch correction failed:", e$message, "\nProceeding without."),
        type = "warning", duration = 3)
      counts
    })
  }
  list(counts = counts, dataset = dataset)
}

step_top_n_sd <- function(counts, dataset, top_n) {
  top_n <- min(top_n, nrow(counts))
  if (nrow(counts) > top_n) {
    sds <- apply(counts, 1, sd, na.rm = FALSE)
    top_features <- names(head(sort(sds, decreasing = TRUE), n = top_n))
    counts <- counts[top_features, , drop = FALSE]
  }
  list(counts = counts, dataset = dataset)
}

#' Prepare count data for ordination (legacy wrapper, kept for reference)
#' @param counts Raw count matrix (genes x samples)
#' @param dataset Dataset with sample metadata
#' @param method Ordination method
#' @param log2_transform Apply log2 transform?
#' @param counts_type Type of counts being used
#' @param batch_column Batch column name for correction (NULL for no correction)
#' @param top_n Number of top features to keep (by SD)
#' @return Processed count matrix (samples x genes for PCA/PCoA, genes x samples for GLM-PCA)
prepare_counts <- function(counts, dataset, method, log2_transform = FALSE, counts_type = NULL, batch_column = NULL, top_n = 2000) {
  # For GLM-PCA, return raw counts as-is (genes x samples)
  if (method == "glmpca") {
    return(counts)
  }

  # For PCA/PCoA: apply transformations
  if (log2_transform) {
    # Only apply log2 to non-logged data
    if (!is.null(counts_type) && counts_type %in% c("TPM", "FPKM", "Raw", "Normalised")) {
      counts <- log2(counts + 1)
    }
  }

  # Apply batch correction if requested
  if (!is.null(batch_column) && !is.null(dataset)) {
    counts <- tryCatch(
      {
        if (any(is.na(dataset[[batch_column]]))) {
          showNotification("Batch correction skipped: batch column contains NA values",
            type = "warning", duration = 3
          )
          counts
        } else {
          limma::removeBatchEffect(counts, dataset[[batch_column]])
        }
      },
      error = function(e) {
        showNotification(paste("Batch correction failed:", e$message, "\nProceeding without."),
          type = "warning", duration = 3
        )
        counts
      }
    )
  }

  # Transpose to samples x genes for PCA/PCoA
  counts <- as.data.frame(t(counts))

  # Keep top N features by standard deviation
  # Cap top_n at the number of available features to prevent issues
  top_n <- min(top_n, ncol(counts))

  if (ncol(counts) > top_n) {
    top_features <- names(head(sort(sapply(counts, sd, na.rm = FALSE), decreasing = TRUE), n = top_n))
    counts <- counts[, top_features]
  }

  return(counts)
}

#' Plot ordination (samples only)
plot_ordination_samples <- function(ordination_result, pc_scores, input, coloursPCA) {
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
    p <- p + scale_shape_manual(values = c(21, 22, 23, 24, 25, 8, 3, 4))
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

  p <- p + scale_fill_manual(breaks = names(coloursPCA), values = as.character(coloursPCA))
  p <- p + theme_bw(base_size = input$textSizePCA)
  face <- switch(input$boldPCA,
    "TRUE" = "bold",
    "FALSE" = "plain"
  )
  p <- p + theme(axis.text = element_text(face = face), axis.title = element_text(face = face))

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

#' Plot biplot (samples + genes)
plot_ordination_biplot <- function(ordination_result, pc_scores, input, coloursPCA, biplot_data, arrow_scale = 1, arrow_colour = "red", selection_method = "overall") {
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
  p <- p + scale_fill_manual(breaks = names(coloursPCA), values = as.character(coloursPCA))
  p <- p + theme_bw(base_size = input$textSizePCA)
  face <- switch(input$boldPCA,
    "TRUE" = "bold",
    "FALSE" = "plain"
  )
  p <- p + theme(axis.text = element_text(face = face), axis.title = element_text(face = face))

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

# Create reactive for design to use across all outputs and download handlers
designReactive <- reactive({
  req(inputDataReactive())
  if (inputDataReactive()$dataType == "RNASeq") {
    return(inputDataReactive()$design)
  } else if (inputDataReactive()$dataType == "proteomics") {
    if (isTruthy(input$contrastSelected)) {
      return(input$contrastSelected)
    } else {
      return("No contrast selected")
    }
  } else if (inputDataReactive()$dataType == "CountQC") {
    return(inputDataReactive()$design)
  }
  return("Unknown")
})

# Main reactive expression for ordination calculations
ordinationResults <- reactive({
  req(inputDataReactive()$dataset)

  # Comprehensive validation before proceeding
  req(
    !is.null(input$pcaFactor1),
    input$pcaFactor1 != "None",
    !is.null(input$pcaGroups),
    length(input$pcaGroups) > 0,
    input$pcaCounts %in% names(inputDataReactive()$countList),
    !is.null(input$ordinationMethod)
  )

  # Validate selected groups exist in the dataset
  available_groups <- levels(as.factor(inputDataReactive()$dataset[[input$pcaFactor1]]))
  valid_groups <- input$pcaGroups[input$pcaGroups %in% available_groups]
  req(length(valid_groups) > 0)

  # Prepare data: start with full dataset and counts
  # GLM-PCA requires raw counts; other methods use the user-selected count type
  datasetOrd <- inputDataReactive()$dataset
  if (input$ordinationMethod == "glmpca") {
    countsOrd <- inputDataReactive()$countList[["Raw"]]
    if (is.null(countsOrd)) {
      showNotification("GLM-PCA requires raw count data which is not available", type = "error", duration = 5)
      return(NULL)
    }
  } else {
    countsOrd <- inputDataReactive()$countList[[input$pcaCounts]]
  }

  # Build cascade step map
  batch_col <- if (input$pcaBatch != "None") input$pcaBatch else NULL
  step_map <- list(
    "Group subsetting" = function(cd) step_group_subset(cd$counts, cd$dataset, input$pcaFactor1, valid_groups),
    "Peptide filter" = function(cd) {
      if (inputDataReactive()$dataType == "proteomics") {
        step_peptide_filter(cd$counts, cd$dataset, inputDataReactive()$seqAnnoList, input$nrPeptidePCA)
      } else {
        cd
      }
    },
    "Log2 transform" = function(cd) {
      if (input$ordinationMethod == "glmpca") {
        cd  # GLM-PCA requires untransformed counts
      } else if (input$pcaLog2) {
        step_log2(cd$counts, cd$dataset, input$pcaCounts)
      } else {
        cd
      }
    },
    "Batch correction" = function(cd) {
      if (input$ordinationMethod == "glmpca") {
        cd  # GLM-PCA handles batch via X covariate, not removeBatchEffect
      } else {
        step_batch_correction(cd$counts, cd$dataset, batch_col)
      }
    },
    "Top N features (SD)" = function(cd) step_top_n_sd(cd$counts, cd$dataset, input$pcaTopN)
  )

  # Execute cascade in user-defined order
  cascade_order <- input$pcaCascadeOrder
  if (is.null(cascade_order) || length(cascade_order) == 0) {
    cascade_order <- c("Group subsetting", "Peptide filter", "Log2 transform", "Batch correction", "Top N features (SD)")
  }

  cd <- list(counts = countsOrd, dataset = datasetOrd)
  for (step_name in cascade_order) {
    if (step_name %in% names(step_map)) {
      cd <- step_map[[step_name]](cd)
    }
  }
  countsOrd <- cd$counts
  datasetOrd <- cd$dataset

  # Transpose for PCA/PCoA (GLM-PCA keeps genes x samples)
  if (input$ordinationMethod != "glmpca") {
    countsPrepared <- as.data.frame(t(countsOrd))
  } else {
    countsPrepared <- countsOrd
  }

  if (is.null(countsPrepared)) {
    showNotification(
      "Unable to prepare count data for ordination. Please check your data and settings.",
      type = "error",
      duration = 5
    )
    return(NULL)
  }

  # Validate we have enough features for ordination
  if (ncol(countsPrepared) < 3) {
    showNotification(
      paste0(
        "Not enough features for ordination (", ncol(countsPrepared), " features available). ",
        "Please try reducing the 'Use top n features' setting or adjust other filters."
      ),
      type = "warning",
      duration = 10
    )
    return(NULL)
  }

  # Calculate ordination based on selected method
  ord_result <- NULL

  if (input$ordinationMethod == "pca") {
    # Require PCA-specific inputs
    req(!is.null(input$pcaCentre), !is.null(input$pcaScale))

    ord_result <- calculate_pca(
      counts = countsPrepared,
      center = input$pcaCentre,
      scale = input$pcaScale
    )
  } else if (input$ordinationMethod == "pcoa") {
    # Require PCoA-specific inputs
    req(!is.null(input$pcoaDistance))

    ord_result <- calculate_pcoa(
      counts = countsPrepared,
      method = input$pcoaDistance,
      nmds = FALSE
    )
  } else if (input$ordinationMethod == "nmds") {
    # Require NMDS-specific inputs
    req(!is.null(input$pcoaDistance))

    ord_result <- calculate_pcoa(
      counts = countsPrepared,
      method = input$pcoaDistance,
      nmds = TRUE
    )
  } else if (input$ordinationMethod == "glmpca") {
    # Require GLM-PCA-specific inputs with strict validation
    req(
      !is.null(input$glmpcaFamily),
      is.character(input$glmpcaFamily),
      nchar(input$glmpcaFamily) > 0,
      input$glmpcaFamily %in% c("poi", "nb", "binom")
    )

    # GLM-PCA uses the cascade-filtered countsPrepared (genes x samples, raw counts)
    # The cascade already applied group subsetting, peptide filter, top N etc.
    # but skipped log2 and batch correction (handled here natively)
    if (is.null(countsPrepared) || nrow(countsPrepared) < 3 || ncol(countsPrepared) < 3) {
      showNotification("GLM-PCA: not enough data after preprocessing", type = "error", duration = 5)
      return(NULL)
    }

    # Build batch covariate matrix for GLM-PCA (instead of removeBatchEffect)
    X_covariates <- NULL
    if (!is.null(batch_col) && batch_col %in% colnames(datasetOrd)) {
      batch_factor <- as.factor(datasetOrd[[batch_col]])
      if (nlevels(batch_factor) > 1) {
        X_covariates <- model.matrix(~ batch_factor)[, -1, drop = FALSE]
      }
    }

    # Show waiter during GLM-PCA computation
    waiter <- waiter::Waiter$new(
      id = "pcaStatic",
      html = tagList(
        waiter::spin_fading_circles(),
        h4("Computing GLM-PCA..."),
        p("This may take 1-2 minutes for large datasets")
      ),
      color = "#527085"
    )
    waiter$show()
    on.exit(waiter$hide(), add = TRUE)

    ord_result <- calculate_glmpca(
      counts = countsPrepared,
      family = input$glmpcaFamily,
      L = 2,
      X = X_covariates
    )
  }

  if (is.null(ord_result)) {
    return(NULL)
  }

  # Add dataset information to scores
  ord_result$scores$Name <- as.character(ord_result$scores$Name)
  ord_result$scores <- left_join(ord_result$scores, datasetOrd, by = "Name")
  ord_result$dataset <- datasetOrd
  ord_result$counts_used <- countsPrepared

  return(ord_result)
})

observeEvent(inputDataReactive(), {
  if (inputDataReactive()$dataType == "RNASeq") {
    output$pcaDesign <- renderText({
      designReactive()
    })
    updateSelectInput(session = session, inputId = "pcaCounts", choices = names(inputDataReactive()$countList), selected = "Normalised + Log2")
  } else if (inputDataReactive()$dataType == "proteomics") {
    output$pcaDesign <- renderText({
      designReactive()
    })
    # Select "transformedData" if it exists, otherwise use the first available count type
    count_choices <- names(inputDataReactive()$countList)
    default_count <- if ("transformedData" %in% count_choices) "transformedData" else count_choices[1]
    updateSelectInput(session = session, inputId = "pcaCounts", choices = count_choices, selected = default_count)
  } else if (inputDataReactive()$dataType == "CountQC") {
    output$pcaDesign <- renderText({
      designReactive()
    })
    updateSelectInput(session = session, inputId = "pcaCounts", choices = names(inputDataReactive()$countList), selected = "VST")
  }
  # Automatically set batch correction if grouping2Name exists
  batchDefault <- "None"
  if (!is.null(inputDataReactive()$param$grouping2Name)) {
    if (inputDataReactive()$param$grouping2Name %in% inputDataReactive()$factors) {
      batchDefault <- inputDataReactive()$param$grouping2Name
    }
  }
  updateSelectInput(
    session = session,
    inputId = "pcaBatch",
    choices = c("None", inputDataReactive()$factors),
    selected = batchDefault
  )

  # Render cascade steps based on data type, log2 checkbox, and count method
  # Only available to staff (LDAP employee role)
  output$pcaCascadeUI <- renderUI({
    if (!exists("ldap_role") || ldap_role != "employee") {
      return(helpText("Advanced preprocessing settings are available to FGCZ staff only."))
    }
    log2_checked <- !is.null(input$pcaLog2) && input$pcaLog2
    is_linear <- !is.null(input$pcaCounts) && input$pcaCounts %in% c("TPM", "FPKM", "Raw", "Normalised")
    show_log2 <- log2_checked && is_linear

    steps <- c("Group subsetting", "Batch correction", "Top N features (SD)")
    if (inputDataReactive()$dataType == "proteomics") {
      steps <- c("Group subsetting", "Peptide filter", "Batch correction", "Top N features (SD)")
    }
    if (show_log2) {
      # Insert Log2 before Batch correction
      batch_pos <- which(steps == "Batch correction")
      steps <- append(steps, "Log2 transform", after = batch_pos - 1)
    }

    sortable::rank_list(
      text = "Preprocessing order",
      labels = as.list(steps),
      input_id = "pcaCascadeOrder",
      options = sortable::sortable_options(multiDrag = FALSE)
    )
  })
  outputOptions(output, "pcaCascadeUI", suspendWhenHidden = FALSE)

  lapply(c(1:2), function(i) {
    selected <- switch(i,
      "1" = inputDataReactive()$groupingName,
      "2" = "None"
    )
    updateSelectInput(
      session = session,
      inputId = paste0("pcaFactor", i),
      choices = c("None", inputDataReactive()$factors),
      selected = selected
    )
  })
  observeEvent(
    {
      input$pcaFactor1
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE,
    {
      updateCheckboxGroupInput(
        session = session,
        inputId = "pcaGroups",
        choices = levels(as.factor(inputDataReactive()$dataset[[input$pcaFactor1]])),
        selected = levels(as.factor(inputDataReactive()$dataset[[input$pcaFactor1]]))
      )
    }
  )
  if (inputDataReactive()$dataType == "proteomics") {
    # Handle nrPeptides slider (only if seqAnnoList exists)
    if (!is.null(inputDataReactive()$seqAnnoList)) {
      if ("nrPeptides" %in% colnames(inputDataReactive()$seqAnnoList[[1]])) {
        output$nrPeptidesPCAUI <- renderUI({
          sliderInput(inputId = "nrPeptidePCA", label = "Minimum number of peptides", min = 0, max = 20, value = 2, step = 1, width = "85%")
        })
        outputOptions(output, "nrPeptidesPCAUI", suspendWhenHidden = FALSE)
      } else {
        output$nrPeptidesPCAUI <- renderUI({
          NULL
        })
      }
    }

    # Update pcaTopN for proteomics data
    if (!is.null(inputDataReactive()$countList)) {
      n_features <- nrow(inputDataReactive()$countList[[1]])
      # Set to 500 or number of features, whichever is smaller
      default_value <- min(500, n_features)
      updateNumericInput(
        session = session,
        inputId = "pcaTopN",
        value = default_value,
        max = n_features
      )
    }
  }
})

# Show notification when GLM-PCA is selected
observeEvent(input$ordinationMethod, {
  if (!is.null(input$ordinationMethod) && input$ordinationMethod == "glmpca") {
    showNotification(
      "GLM-PCA computation may take 1-2 minutes for large datasets. Please be patient.",
      type = "message",
      duration = 8,
      closeButton = TRUE
    )
  }
})

# Separate observer to update biplot gene choices only when ordination results change
# This prevents resetting the selection when user picks genes
observeEvent(ordinationResults(),
  {
    ord_res <- ordinationResults()
    if (!is.null(ord_res) && !is.null(ord_res$loadings)) {
      updateSelectizeInput(
        session = session,
        inputId = "biplotCustomGenes",
        choices = ord_res$loadings$gene,
        server = TRUE
      )
    }
  },
  ignoreNULL = FALSE
)

observeEvent(
  {
    input$pcaFactor1
    input$pcaX
    input$pcaY
    input$pcaTopN
    input$pcaCascadeOrder
    input$pcaGroups
    input$pcaFactor2
    input$pcaCounts
    input$pcaAxesProp
    input$pcaBatch
    input$pcaCentre
    input$pcaScale
    input$textSizePCA
    input$showLinesPCA
    input$showAxesPCA
    input$boldPCA
    input$pcaShowNames
    input$pcaLog2
    input$dotSizePCA
    input$alphaPCA
    input$dotBorderPCA
    input$pcaLabelMaxOverlap
    input$geneLabelNudgePCAX
    input$geneLabelNudgePCAY
    input$geneLabelSizePCA
    input$nrPeptidePCA
    input$pcaAddEllipses
    input$pcaEllipsesAlpha
    # New ordination-specific inputs
    input$ordinationMethod
    input$ordinationPlotType
    input$pcoaDistance
    input$glmpcaFamily
    input$biplotGeneSelectionMethod
    input$biplotNGenes
    input$biplotNGenesPerPC
    input$biplotCustomGenes
    input$biplotCustomGenesText
    input$biplotArrowScale
    input$biplotArrowColour
    input$biplotArrowColourX
    input$biplotArrowColourY
    # Add all color picker inputs as dependencies so plot updates when colors change
    lapply(seq_along(inputDataReactive()$factorLevels), function(i) {
      lapply(seq_along(inputDataReactive()$factorLevels[[i]]), function(j) {
        input[[paste0("GroupColour", names(inputDataReactive()$factorLevels)[[i]], inputDataReactive()$factorLevels[[i]][[j]])]]
      })
    })
    inputDataReactive()
  },
  ignoreInit = TRUE,
  ignoreNULL = FALSE,
  {
    # Get ordination results using the new reactive
    ord_res <- ordinationResults()

    # Silently stop if ordination results aren't ready yet (prevents brief error flash)
    req(ord_res)

    # Get valid groups and colors
    available_groups <- levels(as.factor(inputDataReactive()$dataset[[input$pcaFactor1]]))
    valid_groups <- input$pcaGroups[input$pcaGroups %in% available_groups]

    coloursPCA <- setNames(lapply(valid_groups, function(k) {
      paste0(col2hex(input[[paste0("GroupColour", input$pcaFactor1, k)]]), "FF")
    }), valid_groups)

    # Extract results from ordination
    pc_scores <- ord_res$scores
    pc_loadings <- ord_res$loadings

    # Create variance/stress table
    if (!is.null(ord_res$variance)) {
      # For PCA, PCoA (methods with variance explained)
      pc_eigenvalues <- tibble(
        PC = ord_res$component_names,
        variance = ord_res$variance
      ) %>%
        mutate(pct = (variance / sum(variance)) * 100) %>%
        mutate(pct_cum = cumsum(pct))
      pc_eigenvalues$PC <- factor(pc_eigenvalues$PC, levels = pc_eigenvalues$PC)
    } else if (!is.null(ord_res$stress)) {
      # For NMDS (has stress instead of variance)
      pc_eigenvalues <- tibble(
        PC = ord_res$component_names,
        stress = ord_res$stress
      )
    } else {
      # For GLM-PCA (no variance or stress metrics available)
      pc_eigenvalues <- tibble(
        PC = ord_res$component_names
      )
      pc_eigenvalues$PC <- factor(pc_eigenvalues$PC, levels = pc_eigenvalues$PC)
    }

    # Update axis selectors with smart defaults
    # If the old selection doesn't exist in new choices, reset to first two components
    new_choices <- as.character(pc_eigenvalues$PC)
    x_selected <- if (!is.null(input$pcaX) && input$pcaX %in% new_choices) input$pcaX else new_choices[1]
    y_selected <- if (!is.null(input$pcaY) && input$pcaY %in% new_choices) input$pcaY else new_choices[min(2, length(new_choices))]

    updateSelectInput(session = session, inputId = "pcaX", choices = new_choices, selected = x_selected)
    updateSelectInput(session = session, inputId = "pcaY", choices = new_choices, selected = y_selected)

    # Store results in reactive for safe access by download handlers and brush table
    figuresDataReactive$pc_loadings <- pc_loadings
    figuresDataReactive$pc_eigenvalues <- pc_eigenvalues
    figuresDataReactive$pc_scores <- pc_scores
    figuresDataReactive$ord_method <- ord_res$method
    figuresDataReactive$x_selected <- x_selected
    figuresDataReactive$y_selected <- y_selected
    figuresDataReactive$factors <- inputDataReactive()$factors

    # Update plot type choices based on whether biplot is available
    if (!is.null(pc_loadings)) {
      # Biplot available - offer both options
      updateRadioButtons(
        session = session,
        inputId = "ordinationPlotType",
        choices = c("Samples only" = "samples", "Biplot (samples + genes)" = "biplot"),
        selected = if (!is.null(input$ordinationPlotType)) input$ordinationPlotType else "samples"
      )
    } else {
      # No loadings - biplot not available, force samples only
      updateRadioButtons(
        session = session,
        inputId = "ordinationPlotType",
        choices = c("Samples only" = "samples", "Biplot not available for this method" = "biplot_disabled"),
        selected = "samples"
      )
    }

    # Collect manual genes from dropdown and text area (for biplot)
    # Note: Gene choices are updated in a separate observer to prevent selection reset
    biplot_custom_genes <- NULL
    if (!is.null(input$biplotGeneSelectionMethod) && input$biplotGeneSelectionMethod == "manual") {
      # Get genes from dropdown
      genes_from_dropdown <- input$biplotCustomGenes

      # Get genes from text area (parse by various delimiters)
      genes_from_text <- NULL
      if (!is.null(input$biplotCustomGenesText) && nchar(input$biplotCustomGenesText) > 0) {
        genes_from_text <- strsplit(input$biplotCustomGenesText, "\\s+|/+|,+|\n+|\r+|\t+")[[1]]
        genes_from_text <- genes_from_text[genes_from_text != ""]
      }

      # Combine and deduplicate
      biplot_custom_genes <- unique(c(genes_from_dropdown, genes_from_text))

      # Validate against available genes
      if (!is.null(pc_loadings) && length(biplot_custom_genes) > 0) {
        biplot_custom_genes <- biplot_custom_genes[biplot_custom_genes %in% pc_loadings$gene]
      }
    }

    # Output the scree plot
    output$pcaScree <- renderPlot(
      {
        pc_eigenvalues2 <- pc_eigenvalues
        pc_eigenvalues2$PC <- gsub("PC|Axis|MDS|dim", "", pc_eigenvalues2$PC)
        pc_eigenvalues2$PC <- factor(pc_eigenvalues2$PC, levels = pc_eigenvalues2$PC)
        if (nrow(pc_eigenvalues2) > 15) {
          pc_eigenvalues2 <- pc_eigenvalues2[1:15, ]
        }

        if (!is.null(ord_res$variance)) {
          # Scree plot for variance-based methods (PCA, PCoA)
          pcaScree <- pc_eigenvalues2 %>%
            ggplot(aes(x = PC)) +
            geom_col(aes(y = pct)) +
            geom_line(aes(y = pct_cum, group = 1)) +
            geom_point(aes(y = pct_cum)) +
            labs(x = "Component", y = "Fraction variance explained") +
            theme_classic(base_size = as.numeric(input$textSizePCA))
        } else if (!is.null(ord_res$stress)) {
          # For NMDS, show stress as a single bar
          pcaScree <- ggplot(data = pc_eigenvalues2[1, ], aes(x = "Stress", y = stress)) +
            geom_col() +
            labs(x = "", y = "Stress") +
            theme_classic(base_size = as.numeric(input$textSizePCA))
        } else {
          # For GLM-PCA, show component list (no variance metrics available)
          pcaScree <- ggplot(data = pc_eigenvalues2, aes(x = PC, y = 1)) +
            geom_col(fill = "steelblue") +
            labs(x = "Component", y = "", title = "GLM-PCA Components") +
            theme_classic(base_size = as.numeric(input$textSizePCA)) +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        }
        pcaScree
      },
      width = 400,
      height = 400
    )

    # Output the variance table
    output$pcaVars <- DT::renderDataTable({
      if (!is.null(ord_res$variance)) {
        datatable(pc_eigenvalues[, c("PC", "pct", "pct_cum")],
          rownames = FALSE,
          colnames = c("Component", "Variance explained", "Cumulative variance")
        ) %>%
          formatRound(c("pct", "pct_cum"), digits = 3)
      } else if (!is.null(ord_res$stress)) {
        # For NMDS
        datatable(pc_eigenvalues[, c("PC", "stress")],
          rownames = FALSE,
          colnames = c("Component", "Stress")
        )
      } else {
        # For GLM-PCA (no variance metrics)
        datatable(pc_eigenvalues[, "PC", drop = FALSE],
          rownames = FALSE,
          colnames = c("Component")
        )
      }
    })

    # Output the loading table (only for methods with loadings)
    if (!is.null(pc_loadings) && nrow(pc_loadings) > 0) {
      # Use x_selected and y_selected (already validated) instead of input$pcaX/pcaY
      # This prevents crash when switching methods (input values lag behind updateSelectInput)
      top_genes <- pc_loadings %>%
        dplyr::select(all_of(c("gene", x_selected, y_selected))) %>%
        pivot_longer(matches("PC|Axis|MDS|dim|PCo"), names_to = "Component", values_to = "loading") %>%
        group_by(Component) %>%
        arrange(desc(abs(loading)))

      output$pcaLoadings <- DT::renderDataTable({
        datatable(top_genes, rownames = FALSE) %>% formatRound("loading", digits = 3)
      })
    } else {
      output$pcaLoadings <- DT::renderDataTable({
        datatable(data.frame(Message = "Loadings not available for this method"), rownames = FALSE)
      })
    }

    output$dlTable_ordinationLoadings <- downloadHandler(
      filename = function() {
        generate_filename(
          base_name = paste0(ord_res$method, "_loadings"),
          design_name = designReactive(),
          extension = "xlsx"
        )
      },
      content = function(file) {
        req(figuresDataReactive$pc_loadings)
        export_excel_standard(
          df = as.data.frame(figuresDataReactive$pc_loadings),
          filename = file,
          sheet_name = "Loadings",
          table_name = "Loadings"
        )
      }
    )

    # Generate the main ordination plot
    if (!is.null(input$ordinationPlotType) && input$ordinationPlotType == "biplot") {
      # Require biplot-specific inputs
      req(
        !is.null(input$biplotGeneSelectionMethod),
        !is.null(input$biplotArrowScale),
        !is.null(input$biplotArrowColour)
      )

      # Determine top_n based on selection method
      biplot_top_n <- if (input$biplotGeneSelectionMethod == "per_pc") {
        if (!is.null(input$biplotNGenesPerPC)) input$biplotNGenesPerPC else 5
      } else {
        if (!is.null(input$biplotNGenes)) input$biplotNGenes else 10
      }

      # Prepare biplot data (use validated axis selections)
      biplot_data <- prepare_biplot_data(
        ordination_result = ord_res,
        pc_x = x_selected,
        pc_y = y_selected,
        top_n = biplot_top_n,
        custom_genes = biplot_custom_genes,
        selection_method = input$biplotGeneSelectionMethod
      )

      # If biplot_data is NULL (no loadings available), fall back to samples-only plot
      if (is.null(biplot_data)) {
        # Create temporary input with validated axis selections
        input_validated <- reactiveValuesToList(input)
        input_validated$pcaX <- x_selected
        input_validated$pcaY <- y_selected
        pcaPlot <- plot_ordination_samples(
          ordination_result = ord_res,
          pc_scores = pc_scores,
          input = input_validated,
          coloursPCA = coloursPCA
        )
      } else {
        # Determine arrow colors based on selection method
        arrow_colour <- if (input$biplotGeneSelectionMethod == "per_pc") {
          # For per_pc mode, pass both colors
          list(x = input$biplotArrowColourX, y = input$biplotArrowColourY)
        } else {
          # For other modes, use single color
          input$biplotArrowColour
        }

        # Create temporary input with validated axis selections
        input_validated <- reactiveValuesToList(input)
        input_validated$pcaX <- x_selected
        input_validated$pcaY <- y_selected
        pcaPlot <- plot_ordination_biplot(
          ordination_result = ord_res,
          pc_scores = pc_scores,
          input = input_validated,
          coloursPCA = coloursPCA,
          biplot_data = biplot_data,
          arrow_scale = input$biplotArrowScale,
          arrow_colour = arrow_colour,
          selection_method = input$biplotGeneSelectionMethod
        )
      }
    } else {
      # Samples only plot (use validated axis selections)
      input_validated <- reactiveValuesToList(input)
      input_validated$pcaX <- x_selected
      input_validated$pcaY <- y_selected
      pcaPlot <- plot_ordination_samples(
        ordination_result = ord_res,
        pc_scores = pc_scores,
        input = input_validated,
        coloursPCA = coloursPCA
      )
    }

    figuresDataReactive$pcaStatic <- pcaPlot

    # Paired Plot
    dim_cols <- grep("^(PC|PCo|NMDS|dim)\\d+", colnames(pc_scores))
    dim_cols <- head(dim_cols, 5)  # max 5 dimensions
    req(length(dim_cols) >= 2)
    pairedPlot <- ggpairs(
      data = pc_scores,
      columns = dim_cols,
      aes(fill = .data[[input$pcaFactor1]], alpha = input$alphaPCA),
      upper = list(continuous = wrap("points", alpha = input$alphaPCA, size = input$dotSizePCA, alpha = input$alphaPCA, pch = 21, stroke = input$dotBorderPCA)),
      lower = list(continuous = wrap("points", alpha = input$alphaPCA, size = input$dotSizePCA, alpha = input$alphaPCA, pch = 21, stroke = input$dotBorderPCA)),
      diag = list(continuous = wrap("densityDiag", size = 0.3), discrete = "barDiag", na = "naDiag")
    ) +
      theme_bw(base_size = input$textSizePCA) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_colour_manual(values = as.character(coloursPCA), breaks = names(coloursPCA)) +
      scale_fill_manual(values = as.character(coloursPCA), breaks = names(coloursPCA))

    figuresDataReactive$pcaPaired <- pairedPlot

    output$dlTable_ordinationCoords <- downloadHandler(
      filename = function() {
        generate_filename(
          base_name = paste0(ord_res$method, "_coords"),
          design_name = designReactive(),
          extension = "xlsx"
        )
      },
      content = function(file) {
        export_excel_standard(
          df = figuresDataReactive$pc_scores,
          filename = file,
          sheet_name = "Coordinates",
          table_name = "Coordinates"
        )
      }
    )

    output$pcaNAWarning <- NULL
  }
)

output$pcaStatic <- renderPlot(
  {
    req(figuresDataReactive$pcaStatic)
    figuresDataReactive$pcaStatic
  },
  width = function() {
    as.numeric(input$figWidthPCA)
  },
  height = function() {
    as.numeric(input$figHeightPCA)
  }
)

output$dlPlot_ordination <- downloadHandler(
  filename = function() {
    method_name <- if (!is.null(figuresDataReactive$ord_method)) figuresDataReactive$ord_method else "ordination"
    generate_filename(
      base_name = paste0(input$filename_ordination, "_", method_name),
      design_name = designReactive(),
      extension = tolower(input$downloadFormatPCA)
    )
  },
  content = function(file) {
    req(figuresDataReactive$pcaStatic)
    export_plot_standard(
      plot_object = figuresDataReactive$pcaStatic,
      filename = file,
      format = tolower(input$downloadFormatPCA),
      width = (input$figWidthPCA / 80),
      height = (input$figHeightPCA / 80),
      dpi = as.numeric(input$dpiPCA)
    )
  }
)

output$pcaPaired <- renderPlot(
  {
    req(figuresDataReactive$pcaPaired)
    figuresDataReactive$pcaPaired
  },
  width = function() {
    as.numeric(input$figWidthPCA)
  },
  height = function() {
    as.numeric(input$figHeightPCA)
  }
)
output$dlPlot_ordinationPairs <- downloadHandler(
  filename = function() {
    method_name <- if (!is.null(figuresDataReactive$ord_method)) figuresDataReactive$ord_method else "ordination"
    generate_filename(
      base_name = paste0(input$filename_ordination, "_", method_name, "_paired"),
      design_name = designReactive(),
      extension = tolower(input$downloadFormatPCA)
    )
  },
  content = function(file) {
    req(figuresDataReactive$pcaPaired)
    export_plot_standard(
      plot_object = figuresDataReactive$pcaPaired,
      filename = file,
      format = tolower(input$downloadFormatPCA),
      width = (input$figWidthPCA / 80),
      height = (input$figHeightPCA / 80),
      dpi = as.numeric(input$dpiPCA)
    )
  }
)

output$pcaBrushTable <- DT::renderDataTable({
  pc_scores <- figuresDataReactive$pc_scores
  xvar <- figuresDataReactive$x_selected
  yvar <- figuresDataReactive$y_selected
  brushed_df <- brushedPoints(pc_scores, input$pcaBrush, xvar = xvar, yvar = yvar)

  # Select columns to display
  cols_to_show <- c("Name", xvar, yvar)
  if (!is.null(input$pcaFactor1) && input$pcaFactor1 != "None") {
    cols_to_show <- c(cols_to_show, input$pcaFactor1)
  }
  cols_to_show <- cols_to_show[cols_to_show %in% colnames(brushed_df)]

  DT::datatable(
    data = brushed_df[, cols_to_show, drop = FALSE],
    filter = "top",
    caption = "Click and drag over dots on the PCA plot to see those samples in this table.",
    rownames = FALSE
  ) %>%
    DT::formatRound(columns = c(xvar, yvar), digits = 3)
})
