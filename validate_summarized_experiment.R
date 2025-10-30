#' Validate SummarizedExperiment Against exploreDE Specification
#'
#' This function validates a SummarizedExperiment object against the exploreDE
#' data structure specification, checking for required components, naming
#' conventions, and data consistency.
#'
#' @param se A SummarizedExperiment object to validate
#' @param verbose Logical. If TRUE, prints detailed validation results
#' @return A list containing validation results with pass/fail status and details
#'
#' @examples
#' # Create a valid SummarizedExperiment
#' library(SummarizedExperiment)
#' se <- create_example_se()
#'
#' # Validate it
#' result <- validate_summarized_experiment(se, verbose = TRUE)
#'
#' @export
validate_summarized_experiment <- function(se, verbose = FALSE) {
  # Initialize results
  results <- list(
    overall_status = "PASS",
    errors = character(0),
    warnings = character(0),
    checks = list()
  )

  # Helper function to add results
  add_result <- function(check_name, status, message, details = NULL) {
    results$checks[[check_name]] <- list(
      status = status,
      message = message,
      details = details
    )

    if (status == "ERROR") {
      results$errors <- c(results$errors, paste(check_name, ":", message))
      results$overall_status <- "FAIL"
    } else if (status == "WARNING") {
      results$warnings <- c(results$warnings, paste(check_name, ":", message))
    }
  }

  # Check 1: Basic SummarizedExperiment structure
  if (!inherits(se, "SummarizedExperiment")) {
    add_result("basic_structure", "ERROR", "Object is not a SummarizedExperiment")
    return(results)
  }
  add_result("basic_structure", "PASS", "Object is a valid SummarizedExperiment")

  # Check 2: Required components exist
  rd <- rowData(se)
  rd_check <- if (inherits(rd, "DataFrame") || is.data.frame(rd)) {
    ncol(rd) > 0
  } else if (is.list(rd)) {
    length(rd) > 0
  } else {
    FALSE
  }

  components_check <- list(
    assays = length(assays(se)) > 0,
    rowData = rd_check,
    colData = nrow(colData(se)) > 0,
    metadata = length(metadata(se)) >= 0 # metadata can be empty
  )

  missing_components <- names(components_check)[!sapply(components_check, isTRUE)]

  if (length(missing_components) > 0 && !all(missing_components == "metadata")) {
    add_result(
      "required_components", "ERROR",
      paste("Missing required components:", paste(missing_components, collapse = ", "))
    )
  } else {
    add_result("required_components", "PASS", "All required components present")
  }

  # Check 3: Assays validation
  validate_assays(se, add_result)

  # Check 4: colData validation
  validate_coldata(se, add_result)

  # Check 5: rowData validation
  validate_rowdata(se, add_result)

  # Check 6: Metadata validation
  validate_metadata(se, add_result)

  # Check 7: Consistency between components
  validate_consistency(se, add_result)

  # Print results if verbose
  if (verbose) {
    print_validation_results(results)
  }

  return(results)
}

#' Validate assays component
validate_assays <- function(se, add_result) {
  assays_list <- assays(se)

  if (length(assays_list) == 0) {
    add_result("assays_presence", "ERROR", "No assays found")
    return()
  }

  # Check that at least one numeric assay exists
  add_result("assays_presence", "PASS", paste("Found", length(assays_list), "assay(s):", paste(names(assays_list), collapse = ", ")))

  # Note: Distribution information is not encoded in assay names per specification.
  # Tools should determine distribution statistically using moment-based classification.

  # Check assay dimensions consistency
  if (length(assays_list) > 1) {
    dims <- sapply(assays_list, dim)
    if (!all(apply(dims, 1, function(x) length(unique(x)) == 1))) {
      add_result("assay_dimensions", "ERROR", "Assays have inconsistent dimensions")
    } else {
      add_result("assay_dimensions", "PASS", "All assays have consistent dimensions")
    }
  }
}

#' Validate colData component
validate_coldata <- function(se, add_result) {
  col_data <- colData(se)

  if (nrow(col_data) == 0) {
    add_result("coldata_presence", "ERROR", "colData is empty")
    return()
  }

  # Check for rownames (required)
  if (is.null(rownames(col_data)) || all(rownames(col_data) == "")) {
    add_result("coldata_rownames", "ERROR", "colData missing rownames (sample identifiers)")
  } else if (any(duplicated(rownames(col_data)))) {
    add_result("coldata_rownames", "ERROR", "colData rownames are not unique")
  } else {
    add_result("coldata_rownames", "PASS", "colData rownames are present and unique")
  }

  # Check for sample_id column (optional but recommended)
  if ("sample_id" %in% colnames(col_data)) {
    if (any(duplicated(col_data$sample_id))) {
      add_result("sample_id_unique", "WARNING", "sample_id values are not unique")
    } else {
      add_result("sample_id_unique", "PASS", "sample_id column is present and unique")
    }

    # Check if sample_id matches rownames
    if (!all(rownames(col_data) == col_data$sample_id)) {
      add_result("sample_id_rowname_match", "WARNING", "rownames don't match sample_id (acceptable if both are valid)")
    }
  } else {
    add_result("sample_id_column", "INFO", "No sample_id column (using rownames as identifiers)")
  }

  # Check factor columns
  factor_cols <- colnames(col_data)[grepl("^factor_", colnames(col_data))]
  if (length(factor_cols) > 0) {
    add_result("factor_columns", "PASS", paste("Found factor columns:", paste(factor_cols, collapse = ", ")))
  } else {
    add_result("factor_columns", "WARNING", "No factor columns found")
  }
}

#' Validate rowData component
validate_rowdata <- function(se, add_result) {
  # Access rowData - could be a DataFrame or a list depending on implementation
  row_data <- rowData(se)

  # Convert to list if it's a DataFrame (standard SummarizedExperiment)
  if (inherits(row_data, "DataFrame") || is.data.frame(row_data)) {
    # Standard SummarizedExperiment has single DataFrame
    if (ncol(row_data) == 0) {
      add_result("rowdata_presence", "ERROR", "rowData is empty")
      return()
    }
    # Treat as single table
    row_data_list <- list(rowData = row_data)
  } else if (is.list(row_data)) {
    # Custom implementation with list of DataFrames
    row_data_list <- row_data
    if (length(row_data_list) == 0) {
      add_result("rowdata_presence", "ERROR", "rowData is empty")
      return()
    }
  } else {
    add_result("rowdata_presence", "ERROR", "rowData has unexpected type")
    return()
  }

  # Check for feature annotation table (gene_id, protein_id, peptide_id, etc.)
  feature_tables <- names(row_data_list)[grepl("_id$|feature_annotation", names(row_data_list))]
  if (length(feature_tables) == 0) {
    add_result("feature_annotation", "WARNING", "No feature annotation table found (expected: gene_id, protein_id, peptide_id, etc.)")
  } else {
    add_result("feature_annotation", "PASS", paste("Found feature annotation table:", paste(feature_tables, collapse = ", ")))
  }

  # Check DE tables
  de_tables <- names(row_data_list)[grepl("^DE_", names(row_data_list))]
  if (length(de_tables) == 0) {
    add_result("de_tables", "WARNING", "No DE tables found")
  } else {
    add_result("de_tables", "PASS", paste("Found DE tables:", paste(de_tables, collapse = ", ")))
  }

  # Validate each rowData table
  for (table_name in names(row_data_list)) {
    validate_rowdata_table(row_data_list[[table_name]], table_name, add_result)
  }
}

#' Validate individual rowData table
validate_rowdata_table <- function(table, table_name, add_result) {
  # Check required columns for feature annotation (tables ending with _id or containing feature_annotation)
  if (grepl("_id$|feature_annotation", table_name)) {
    if (!"description" %in% colnames(table)) {
      add_result(
        paste0("table_", table_name, "_description"), "WARNING",
        "Feature annotation table missing 'description' column"
      )
    } else {
      add_result(paste0("table_", table_name, "_description"), "PASS", "Description column present")
    }
  }

  # Check DE table structure
  if (grepl("^DE_", table_name)) {
    # Check for effect columns
    effect_cols <- colnames(table)[grepl("^effect_", colnames(table))]
    if (length(effect_cols) == 0) {
      add_result(paste0("table_", table_name, "_effect"), "ERROR", "DE table missing effect columns")
    } else {
      add_result(
        paste0("table_", table_name, "_effect"), "PASS",
        paste("Found effect columns:", paste(effect_cols, collapse = ", "))
      )
    }

    # Check for score columns
    score_cols <- colnames(table)[grepl("^score_", colnames(table))]
    if (length(score_cols) == 0) {
      add_result(paste0("table_", table_name, "_score"), "ERROR", "DE table missing score columns")
    } else {
      add_result(
        paste0("table_", table_name, "_score"), "PASS",
        paste("Found score columns:", paste(score_cols, collapse = ", "))
      )
    }
  }

  # Check column naming conventions
  valid_prefixes <- c("label_", "hierarchy_", "effect_", "score_", "factor_")
  # Check if columns follow naming conventions or are description
  valid_cols <- sapply(colnames(table), function(col) {
    grepl(paste0("^(", paste(valid_prefixes, collapse = "|"), ")"), col) ||
      col == "description"
  })

  invalid_cols <- colnames(table)[!valid_cols]

  if (length(invalid_cols) > 0) {
    add_result(
      paste0("table_", table_name, "_naming"), "WARNING",
      paste("Columns with non-standard names:", paste(invalid_cols, collapse = ", "))
    )
  }
}

#' Validate metadata component
validate_metadata <- function(se, add_result) {
  metadata_list <- metadata(se)

  if (length(metadata_list) == 0) {
    add_result("metadata_presence", "WARNING", "No metadata found")
    return()
  }

  # Check for DE test metadata
  if (!"de_tests" %in% names(metadata_list)) {
    add_result("de_tests_metadata", "WARNING", "No de_tests metadata found")
  } else {
    de_tests <- metadata_list$de_tests

    if (length(de_tests) == 0) {
      add_result("de_tests_metadata", "WARNING", "de_tests metadata is empty")
    } else {
      add_result("de_tests_metadata", "PASS", paste("Found", length(de_tests), "DE test(s)"))

      # Validate each DE test
      for (test_name in names(de_tests)) {
        validate_de_test_metadata(de_tests[[test_name]], test_name, add_result)
      }
    }
  }
}

#' Validate DE test metadata
validate_de_test_metadata <- function(de_test, test_name, add_result) {
  # Required fields per specification
  required_fields <- c("assay_used", "factor_used")
  missing_fields <- required_fields[!required_fields %in% names(de_test)]

  if (length(missing_fields) > 0) {
    add_result(
      paste0("de_test_", test_name, "_required"), "ERROR",
      paste("Missing required fields:", paste(missing_fields, collapse = ", "))
    )
  } else {
    add_result(paste0("de_test_", test_name, "_required"), "PASS", "All required fields present")
  }

  # Optional fields (check if present)
  optional_fields <- c("contrast_formula", "model")
  present_optional <- optional_fields[optional_fields %in% names(de_test)]
  if (length(present_optional) > 0) {
    add_result(
      paste0("de_test_", test_name, "_optional"), "INFO",
      paste("Optional fields present:", paste(present_optional, collapse = ", "))
    )
  }
}

#' Validate consistency between components
validate_consistency <- function(se, add_result) {
  # Check row consistency
  if (length(assays(se)) > 0) {
    assay_rows <- nrow(assays(se)[[1]])
    rd <- rowData(se)

    # Handle both DataFrame and list formats
    if (inherits(rd, "DataFrame") || is.data.frame(rd)) {
      # Standard SummarizedExperiment - single DataFrame
      if (nrow(rd) != assay_rows) {
        add_result("row_consistency", "ERROR", "rowData and assays have different row counts")
      } else {
        add_result("row_consistency", "PASS", "Row dimensions are consistent")
      }
    } else if (is.list(rd)) {
      # Custom implementation with list of DataFrames
      if (length(rd) == 0) {
        add_result("row_consistency", "ERROR", "rowData is empty but assays exist")
      } else {
        # Check if all rowData tables have same number of rows
        rowdata_lengths <- sapply(rd, nrow)
        if (length(unique(rowdata_lengths)) > 1) {
          add_result("row_consistency", "ERROR", "rowData tables have inconsistent row counts")
        } else if (rowdata_lengths[1] != assay_rows) {
          add_result("row_consistency", "ERROR", "rowData and assays have different row counts")
        } else {
          add_result("row_consistency", "PASS", "Row dimensions are consistent")
        }
      }
    }
  }

  # Check column consistency
  if (length(assays(se)) > 0) {
    assay_cols <- ncol(assays(se)[[1]])
    coldata_rows <- nrow(colData(se))

    if (coldata_rows != assay_cols) {
      add_result("col_consistency", "ERROR", "colData and assays have different column counts")
    } else {
      add_result("col_consistency", "PASS", "Column dimensions are consistent")
    }
  }

  # Check DE table names match metadata
  rd <- rowData(se)
  rd_names <- if (inherits(rd, "DataFrame") || is.data.frame(rd)) {
    colnames(rd)
  } else if (is.list(rd)) {
    names(rd)
  } else {
    character(0)
  }

  de_tables <- rd_names[grepl("^DE_", rd_names)]
  de_metadata_names <- names(metadata(se)$de_tests)

  if (length(de_tables) > 0 && length(de_metadata_names) > 0) {
    expected_de_names <- paste0("DE_", de_metadata_names)
    missing_de_tables <- expected_de_names[!expected_de_names %in% de_tables]
    extra_de_tables <- de_tables[!de_tables %in% expected_de_names]

    if (length(missing_de_tables) > 0) {
      add_result(
        "de_consistency", "ERROR",
        paste("Missing DE tables for metadata:", paste(missing_de_tables, collapse = ", "))
      )
    } else if (length(extra_de_tables) > 0) {
      add_result(
        "de_consistency", "WARNING",
        paste("Extra DE tables without metadata:", paste(extra_de_tables, collapse = ", "))
      )
    } else {
      add_result("de_consistency", "PASS", "DE tables and metadata are consistent")
    }
  }
}

#' Print validation results
print_validation_results <- function(results) {
  cat("=== SummarizedExperiment Validation Results ===\n")
  cat("Overall Status:", results$overall_status, "\n\n")

  if (length(results$errors) > 0) {
    cat("ERRORS:\n")
    for (error in results$errors) {
      cat("  [ERROR]", error, "\n")
    }
    cat("\n")
  }

  if (length(results$warnings) > 0) {
    cat("WARNINGS:\n")
    for (warning in results$warnings) {
      cat("  [WARN]", warning, "\n")
    }
    cat("\n")
  }

  cat("DETAILED RESULTS:\n")
  for (check_name in names(results$checks)) {
    check <- results$checks[[check_name]]
    status_icon <- if (check$status == "PASS") "[PASS]" else if (check$status == "WARNING") "[WARN]" else if (check$status == "INFO") "[INFO]" else "[ERROR]"
    cat("  ", status_icon, check_name, ":", check$message, "\n")
  }
}

#' Create example SummarizedExperiment for testing
create_example_se <- function() {
  library(SummarizedExperiment)

  # Define identifiers
  gene_ids <- c("GENE001", "GENE002", "GENE003")
  sample_ids <- c("sample1", "sample2", "sample3", "sample4")

  # Create assays
  counts <- matrix(
    c(
      100, 150, 200, 180,
      50, 60, 20, 25,
      1000, 1100, 1050, 1020
    ),
    nrow = 3, ncol = 4,
    dimnames = list(gene_ids, sample_ids)
  )

  # Transformed abundance (log2)
  vst <- log2(counts + 1)

  # Normalized abundance (TPM - example)
  tpm <- counts / colSums(counts) * 1e6

  # Create colData
  col_data <- DataFrame(
    row.names = sample_ids,
    sample_id = sample_ids,
    sample_label = c("Ctrl_1", "Ctrl_2", "Treat_1", "Treat_2"),
    factor_condition = factor(c("control", "control", "treated", "treated")),
    factor_batch = factor(c("A", "B", "A", "B")),
    organism = "Homo sapiens",
    feature_level = "gene"
  )

  # Create rowData
  feature_annotation <- DataFrame(
    row.names = gene_ids,
    description = c("Tumor protein p53", "BRCA1 DNA repair", "EGFR receptor"),
    label_gene_symbol = c("TP53", "BRCA1", "EGFR"),
    label_entrez_id = c("7157", "672", "1956"),
    label_global_base_mean = c(150.5, 42.3, 1050.2)
  )

  de_treated_vs_control <- DataFrame(
    row.names = gene_ids,
    effect_log2fc = c(2.5, -1.8, 0.3),
    score_p_value = c(0.001, 0.005, 0.45),
    score_fdr = c(0.01, 0.03, 0.67),
    label_de_base_mean = c(155.2, 38.1, 1048.9),
    label_t_statistic = c(8.2, -5.1, 0.8)
  )

  # Create metadata
  metadata_list <- list(
    de_tests = list(
      treated_vs_control = list(
        assay_used = "vst",
        factor_used = "factor_condition",
        contrast_formula = "factor_condition_treated - factor_condition_control",
        model = "DESeq2"
      )
    ),
    hierarchy_definition = c("gene_id"),
    ref_build = "GRCh38"
  )

  # Construct SummarizedExperiment
  se <- SummarizedExperiment(
    assays = list(
      counts = counts,
      vst = vst,
      tpm = tpm
    ),
    rowData = list(
      gene_id = feature_annotation,
      DE_treated_vs_control = de_treated_vs_control
    ),
    colData = col_data,
    metadata = metadata_list
  )

  return(se)
}
