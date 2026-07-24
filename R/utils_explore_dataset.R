# =============================================================================
# Internal helpers for the ExploreDataset R6 class hierarchy.
#
# All identifiers in this file are package-internal (not exported) and are
# prefixed with `eds_` to avoid colliding with similar helpers in
# `helperFunctions.R` (e.g. `build_feature_id`, `standardise_columns`)
# whose contracts differ.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Column-name regexes for DE tables. Map any upstream tool's column to a
# canonical AnnData-spec name (`effect_log2fc`, `score_pvalue`, `score_fdr`).
eds_de_column_patterns <- list(
  effect_log2fc = "^diff$|^log2[._]?(fold)?[._]?(change|ratio|fc)?$|^logfc$|^estimate$",
  score_pvalue  = "^p[._]?value$",
  score_fdr     = "^fdr$|^adj[._]?p[._]?value$|^padj$|^q[._]?value$|^p[._]?adjusted$"
)

# Column-name regexes for proteomics rowData / contrast DataFrames. Renames
# (not copies) the matched column to the canonical name.
eds_rowdata_column_patterns <- list(
  protein_Id  = "^protein[._]?id$",
  peptide_Id  = "^peptide[._]?id$|^peptide$",
  site        = "^site$",
  gene_name   = "^gene[._]?name$|^gene[._]?symbol$",
  description = "^description$|^desc$",
  fasta.id    = "^fasta[._]?id$",
  modelName   = "^model[._]?name$",
  nrPeptides  = "^nr[._]?peptides$|^n[._]?peptides$|^num[._]?peptides$"
)

# Add canonical DE columns to `df` (copy, not rename), so the upstream tool's
# native column names are preserved alongside. Only adds if the canonical
# column doesn't already exist.
eds_standardise_de_columns <- function(df) {
  current <- colnames(df)
  for (canon in names(eds_de_column_patterns)) {
    pattern <- eds_de_column_patterns[[canon]]
    idx <- grep(pattern, current, ignore.case = TRUE)
    if (length(idx) >= 1L && !canon %in% colnames(df)) {
      df[[canon]] <- df[[idx[1]]]
    }
  }
  df
}

# Rename proteomics rowData columns to canonical names. Each canonical name
# matches at most one source column.
eds_standardise_rowdata_columns <- function(df) {
  current <- colnames(df)
  for (canon in names(eds_rowdata_column_patterns)) {
    pattern <- eds_rowdata_column_patterns[[canon]]
    idx <- grep(pattern, current, ignore.case = TRUE)
    if (length(idx) == 1L) current[idx] <- canon
  }
  colnames(df) <- current
  df
}

# Build a unique, human-readable feature_id for proteomics features. Tries
# `protein_Id` alone, then `protein_Id + site` or `protein_Id + peptide_Id`,
# then cleaned rownames. Returns a character vector.
eds_build_proteomics_feature_id <- function(df) {
  if (!"protein_Id" %in% colnames(df)) {
    return(gsub("~[^~]+~", "_", rownames(df)))
  }
  pid <- df$protein_Id
  if (!anyDuplicated(pid)) return(pid)
  for (col in intersect(c("site", "peptide_Id"), colnames(df))) {
    combined <- paste(pid, gsub("~", "_", df[[col]]), sep = "_")
    if (!anyDuplicated(combined)) return(combined)
  }
  gsub("~[^~]+~", "_", rownames(df))
}

# SUSHI param$species is rarely set; fall back to the first slash-segment of
# refBuild ("Homo_sapiens/GENCODE/..." → "Homo sapiens").
eds_infer_organism <- function(param) {
  sp <- param$species
  if (!is.null(sp) && !is.na(sp) && nzchar(sp)) return(sp)
  rb <- param$refBuild
  if (is.null(rb) || is.na(rb) || !nzchar(rb)) return(NA_character_)
  parts <- strsplit(rb, "/", fixed = TRUE)[[1]]
  if (!length(parts)) return(NA_character_)
  gsub("_", " ", parts[1], fixed = TRUE)
}

# Recursive flatten of nested clusterProfiler results into a single named list.
# clusterProfiler outputs are typically nested as
# `list(BP = enrichResult, MF = enrichResult, ...)` or deeper for ORA with
# multiple directions (`list(BP = list(upGenes = ..., downGenes = ...), ...)`).
eds_flatten_cp_results <- function(x, prefix = "") {
  if (inherits(x, c("enrichResult", "gseaResult"))) {
    nm <- if (nzchar(prefix)) prefix else "default"
    return(setNames(list(x), nm))
  }
  if (is.list(x) && length(x)) {
    out <- list()
    nms <- names(x) %||% as.character(seq_along(x))
    for (i in seq_along(x)) {
      sub_prefix <- if (nzchar(prefix)) paste0(prefix, "/", nms[i]) else nms[i]
      out <- c(out, eds_flatten_cp_results(x[[i]], sub_prefix))
    }
    return(out)
  }
  list()
}

# Regex-safe identifier (used for contrast names built from group labels).
eds_safe <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

# Compute VST on a genomics SE — mirrors the legacy convert_genomics_se
# logic exactly: integer-cast counts; vst() for >=1000 rows, fall back to
# varianceStabilizingTransformation otherwise (and on any vst() error).
eds_compute_vst <- function(se) {
  forVST <- as.matrix(SummarizedExperiment::assay(se, "counts"))
  mode(forVST) <- "integer"
  tryCatch({
    if (nrow(forVST) < 1000) {
      DESeq2::varianceStabilizingTransformation(forVST, blind = FALSE)
    } else {
      DESeq2::vst(forVST, blind = FALSE)
    }
  }, error = function(e) {
    DESeq2::varianceStabilizingTransformation(forVST, blind = FALSE)
  })
}

# Apply fid_lookup remap to a matrix's rownames. Drops rows whose source
# rowname isn't in the lookup.
.eds_remap_rownames <- function(m, lk) {
  if (is.null(lk)) return(m)
  new_rn <- unname(lk[rownames(m)])
  keep   <- !is.na(new_rn)
  m      <- m[keep, , drop = FALSE]
  rownames(m) <- new_rn[keep]
  m
}

# Resolve a path / URL / directory to an in-memory SE. Mirrors the load
# behaviour of the existing app's SE-fetching code:
#   - http(s) URL                  → readRDS(url(...))
#   - .RData file (Mageck EzResult) → strip suffix, fall through to directory
#   - .rds file                    → readRDS
#   - .qs2 file                    → qs2::qs_read
#   - directory                    → search for the standard SUSHI artefacts
eds_load_se <- function(input) {
  if (grepl("^https?://", input)) {
    return(readRDS(url(input)))
  }
  if (grepl("\\.RData$", input, ignore.case = TRUE)) {
    dir_input <- gsub("/(result|counts)-.*-EzResult\\.RData$", "", input)
    return(eds_load_se(dir_input))
  }
  if (grepl("\\.rds$", input, ignore.case = TRUE) && file.exists(input)) {
    return(readRDS(input))
  }
  if (grepl("\\.qs2$", input, ignore.case = TRUE) && file.exists(input)) {
    return(qs2::qs_read(input, nthreads = 4))
  }
  if (dir.exists(input)) {
    for (fname in c("deResult.rds",      "deResult.qs2",
                    "rawData_final.rds", "rawData_final.qs2",
                    "rawData.rds",       "rawData.qs2")) {
      fp <- file.path(input, fname)
      if (file.exists(fp)) {
        if (grepl("\\.qs2$", fname)) return(qs2::qs_read(fp, nthreads = 4))
        return(readRDS(fp))
      }
    }
    stop("No known SE artefact in directory: ", input)
  }
  stop("Cannot interpret input: ", input)
}
