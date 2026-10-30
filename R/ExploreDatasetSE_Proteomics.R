# --- SAINTexpress transitional support --------------------------------------
# SAINTexpress AP-MS SEs reach the app via the legacy ?data=<.rds> URL path,
# but break the standard prolfqua assumptions: metadata$param is NULL and the
# per-contrast DE columns are log2_EFCs / SaintScore / BFDR with NO raw p-value
# (SaintScore is a probability where high = significant; BFDR is the FDR). This
# is a deliberately localised hack so those existing .rds files work from URL;
# it should be deleted once producers emit canonical h5ad across the board (the
# h5ad path already handles this shape declaratively via uns$exploreDE$mapping).
# Keyed on SaintScore alone: standardise_columns() (via the column dictionary's
# "log2" / "fdr" substring patterns) may already have renamed log2_EFCs/BFDR to
# log2Ratio/fdr by the time this runs, but it never touches SaintScore.
.se_is_saint <- function(contrast_df) "SaintScore" %in% colnames(contrast_df)

# Map a SAINT contrast frame's stat columns onto the canonical legacy names the
# rest of the app expects (log2Ratio / fdr) and synthesise an all-NA pValue so
# p-value-driven code (order(sa$pValue), select(pValue), ...) neither errors on a
# missing column nor mis-thresholds. Idempotent + only touches SAINT frames, so
# it is safe to call after standardise_columns() on any proteomics frame. Shared
# by the ExploreDataset backend and the legacy gStoreList glue so they can't
# drift. Delete alongside the rest of the SAINT hack once h5ad is universal.
.saint_standardise_de <- function(d) {
  if (!.se_is_saint(d)) {
    return(d)
  }
  if ("log2_EFCs" %in% colnames(d) && !"log2Ratio" %in% colnames(d)) {
    d$log2Ratio <- suppressWarnings(as.numeric(d$log2_EFCs))
  }
  if ("BFDR" %in% colnames(d) && !"fdr" %in% colnames(d)) {
    d$fdr <- suppressWarnings(as.numeric(d$BFDR))
  }
  if (!"pValue" %in% colnames(d)) {
    d$pValue <- NA_real_
  }
  d
}

# A p-value-safe ordering index for a DE/seqAnno frame: by raw p-value when it
# carries usable values, else by fdr, else the identity order. Guards against
# order(NULL) (missing column) and all-NA p-value producers (SAINTexpress).
.de_order_idx <- function(d) {
  col <- if ("pValue" %in% colnames(d) && any(!is.na(d$pValue))) {
    "pValue"
  } else if ("fdr" %in% colnames(d)) {
    "fdr"
  } else {
    NULL
  }
  if (is.null(col)) seq_len(nrow(d)) else order(d[[col]], na.last = TRUE)
}

#' Concrete subclass for FGCZ proteomics `SummarizedExperiment` outputs
#'
#' Reads from a prolfqua-derived SE. Top-level `rowData(se)` is sparse;
#' feature metadata and per-contrast DE statistics live nested inside
#' `rowData(se)$constrast_<name>` DataFrames (note the SUSHI typo is
#' preserved). The first contrast is the source of truth for feature
#' metadata. `feature_id` is **built** from `protein_Id` (+ `site` /
#' `peptide_Id` to disambiguate), and the assay matrices' rownames are
#' remapped via `.fid_lookup` since SE rownames do not match the built
#' `feature_id`s.
#'
#' Inherits the accessor contract from [ExploreDataset] via
#' [ExploreDatasetSE]; overrides `.build_config()` and `.prep()`.
#'
#' @export
ExploreDatasetSE_Proteomics <- R6::R6Class(
  classname = "ExploreDatasetSE_Proteomics",
  inherit = ExploreDatasetSE,

  private = list(
    .build_config = function(se) {
      cd <- data.frame(SummarizedExperiment::colData(se), check.names = FALSE, stringsAsFactors = FALSE)

      # Factors: legacy convert_proteomics_se exclusion regex. Names are
      # taken verbatim from colData (NO `G_` stripping) to match the legacy
      # gStoreList$factors / $groupingName values that modules consume.
      factor_cols <- grep(
        "Sample|Name|names|channel|raw\\.file",
        colnames(cd),
        invert = TRUE,
        ignore.case = TRUE,
        value = TRUE
      )
      factor_names <- factor_cols
      factor_levels <- setNames(
        lapply(factor_cols, function(fc) {
          v <- cd[[fc]]
          if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
        }),
        factor_names
      )

      # Contrasts from rowData `^constrast_` columns (typo preserved).
      rd_nms <- names(SummarizedExperiment::rowData(se))
      contrasts <- gsub("^constrast_", "", grep("^constrast_", rd_nms, value = TRUE))

      # Inspect first contrast to determine flavour (protein/phospho/peptide).
      # Uses the package-level standardise_columns() which renames both
      # rowdata-annotation columns and DE-stat columns in a single pass.
      feature_id_kind <- "protein"
      saint <- FALSE
      if (length(contrasts)) {
        df1 <- standardise_columns(
          data.frame(
            SummarizedExperiment::rowData(se)[[paste0("constrast_", contrasts[1])]],
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
        )
        if ("site" %in% colnames(df1)) {
          feature_id_kind <- "phospho"
        } else if ("peptide_Id" %in% colnames(df1)) {
          feature_id_kind <- "peptide"
        }
        saint <- .se_is_saint(df1) # standardise_columns leaves SaintScore as-is
      }

      # SAINTexpress reference = the grouping level flagged as control
      # (colData CONTROL == "C", e.g. eGFP); the contrast names carry no `_vs_`.
      reference_default <- NA_character_
      if (saint && "CONTROL" %in% colnames(cd) && (factor_names[1] %||% "") %in% colnames(cd)) {
        ctrl_lvls <- unique(as.character(cd[[factor_names[1]]])[
          toupper(trimws(as.character(cd$CONTROL))) == "C"
        ])
        if (length(ctrl_lvls) == 1L) reference_default <- ctrl_lvls
      }

      cm <- setNames(
        lapply(contrasts, function(cn) {
          grp <- strsplit(cn, "_vs_|--over--")[[1]]
          list(
            test_group = grp[1] %||% NA_character_,
            reference_group = if (saint) reference_default else (grp[2] %||% NA_character_),
            factor_name = factor_names[1] %||% NA_character_,
            method = if (saint) "SAINTexpress" else "prolfqua",
            thresholds = list(fdr = 0.05, log2fc = 1),
            # SAINT carries a second, display-only significance axis (SaintScore,
            # high = significant); mirrors the h5ad backend's score label.
            score_label = if (saint) "SaintScore" else NULL
          )
        }),
        contrasts
      )

      ass_nms <- SummarizedExperiment::assayNames(se)
      primary_key <- intersect(c("transformedData", "transformed", "intensities", "log2"), ass_nms)
      primary_key <- if (length(primary_key)) primary_key[1] else ass_nms[1]
      raw_key <- intersect(c("intensities_raw", "raw", "intensities"), ass_nms)
      raw_key <- if (length(raw_key)) raw_key[1] else primary_key

      # Derive data_type from feature_id_kind so server modules can dispatch
      # at the flavour level if they want; omics_class keeps the umbrella.
      data_type_flavour <- switch(
        feature_id_kind,
        protein = "proteomics",
        peptide = "peptidomics",
        phospho = "phosphoproteomics",
        "proteomics"
      )

      has_raw <- !identical(raw_key, primary_key) &&
        !is.null(raw_key) &&
        !is.na(raw_key)
      assay_map_l <- c(
        list(intensities = primary_key),
        if (has_raw) list(intensities_raw = raw_key)
      )
      # Display labels = original SE assayNames (matches legacy
      # convert_proteomics_se which uses them verbatim).
      assay_labels_l <- c(
        intensities = primary_key,
        if (has_raw) c(intensities_raw = raw_key)
      )

      list(
        schema_version = "exploreDE/config/1.0",
        omics_class = "proteomics",
        data_type = data_type_flavour,
        organism = NA_character_,
        ref_build = NA_character_,
        feature_id_kind = feature_id_kind,
        factors = factor_names,
        factor_levels = factor_levels,
        primary_grouping = factor_names[1] %||% NA_character_,
        secondary_grouping = NA_character_,
        assay_map = assay_map_l,
        assay_labels = assay_labels_l,
        background_expression = 1,
        raw_assay_name = "intensities",
        primary_assay_name = "intensities",
        display_name_field = NA_character_, # handled in .prep()
        description_field = "description",
        primary_contrast = if (length(contrasts)) contrasts[1] else NA_character_,
        contrast_meta = cm,
        de_display_extras = if (saint) "SaintScore" else character(),
        enrichment_kinds = character(),
        enrichment_locator = list(),
        sgrna_to_gene = NULL,
        proteomics = NULL, # filled by .prep()
        source = list(
          backend_kind = "se",
          schema_version = "exploreDE/config/1.0",
          built_at = Sys.time()
        )
      )
    },

    .prep = function() {
      se <- private$.se
      contrasts <- self$contrasts()

      if (!length(contrasts)) {
        private$.features_df <- data.frame(feature_id = character(), display_name = character())
        private$.de_dfs <- list()
        private$.fid_lookup <- setNames(rownames(se), rownames(se))
        return(invisible(NULL))
      }

      # Extract and standardise each contrast DataFrame. Uses the
      # package-level standardise_columns() so DE columns are renamed to
      # legacy `log2Ratio` / `pValue` / `fdr` (matches gStoreList$seqAnnoList).
      # Canonical AnnData-spec names are added later via
      # eds_standardise_de_columns() as aliases.
      raw_dfs <- setNames(
        lapply(contrasts, function(cn) {
          data.frame(
            SummarizedExperiment::rowData(se)[[paste0("constrast_", cn)]],
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
        }),
        contrasts
      )
      std_dfs <- lapply(raw_dfs, standardise_columns)

      # SAINTexpress: map log2_EFCs -> log2Ratio and BFDR -> fdr wherever
      # standardise_columns() hasn't already, add an all-NA pValue so the
      # modules that `select(pValue)` don't error (and has_pvalue() reports
      # FALSE, driving the FDR-only degradation), and leave SaintScore in place
      # as a display-only extra. Transitional; drops out for the h5ad path.
      std_dfs <- lapply(std_dfs, .saint_standardise_de)

      # Reference DataFrame (first contrast, NOT yet sorted by p-value) so
      # feature_id ↔ rownames(se) alignment is preserved for the fid_lookup.
      ref_df <- std_dfs[[1]]
      rownames(ref_df) <- rownames(se)[seq_len(nrow(ref_df))]
      ref_df$feature_id <- build_feature_id(ref_df)

      # fid_lookup: SE rowname → built feature_id
      private$.fid_lookup <- setNames(ref_df$feature_id, rownames(ref_df))

      # At protein-level, the built feature_id equals protein_Id — legacy
      # convert_proteomics_se drops the redundant protein_Id column in that
      # case. Replicate so eds$proteomics_choices("protein") matches
      # gs$protein_choices behaviour (NULL for protein-level data).
      if ("protein_Id" %in% colnames(ref_df) && identical(ref_df$feature_id, ref_df$protein_Id)) {
        ref_df$protein_Id <- NULL
      }

      # features_df: ref_df with display_name (gene_name → protein_Id →
      # feature_id fallback chain).
      display_priority <- c("gene_name", "protein_Id")
      display_priority <- intersect(display_priority, colnames(ref_df))
      display <- rep(NA_character_, nrow(ref_df))
      for (col in display_priority) {
        blank <- is.na(display) | !nzchar(display)
        if (!any(blank)) {
          break
        }
        vals <- as.character(ref_df[[col]])
        display[blank] <- vals[blank]
      }
      blank <- is.na(display) | !nzchar(display)
      display[blank] <- ref_df$feature_id[blank]
      ref_df$display_name <- display

      # Tidy description if present.
      if ("description" %in% colnames(ref_df)) {
        ref_df$description <- gsub(" \\[.*", "", as.character(ref_df$description))
      }
      rownames(ref_df) <- ref_df$feature_id
      private$.features_df <- ref_df

      # de_dfs: standardised columns + feature_id (built, not raw rownames),
      # sorted by pValue where possible. Canonical AnnData names are added
      # as aliases via eds_standardise_de_columns().
      de_dfs <- setNames(
        lapply(seq_along(contrasts), function(i) {
          d <- std_dfs[[i]]
          d$feature_id <- ref_df$feature_id[seq_len(nrow(d))]
          # Drop redundant protein_Id at protein level (mirrors above).
          if ("protein_Id" %in% colnames(d) && identical(d$feature_id, as.character(d$protein_Id))) {
            d$protein_Id <- NULL
          }
          d <- eds_standardise_de_columns(d)
          # usedInTest: which features entered the test. Proxy: a usable raw
          # p-value; fall back to fdr for p-value-less producers (SAINTexpress)
          # so DE views are not silently emptied. Only added when absent, so
          # frames that already carry it (standard prolfqua) are untouched.
          if (!"usedInTest" %in% colnames(d)) {
            has_p <- "pValue" %in% colnames(d) && any(!is.na(d$pValue))
            d$usedInTest <- if (has_p) {
              !is.na(d$pValue)
            } else if ("fdr" %in% colnames(d)) {
              !is.na(d$fdr)
            } else {
              TRUE
            }
          }
          if ("pValue" %in% colnames(d) && any(!is.na(d$pValue))) {
            d <- d[order(d$pValue, na.last = TRUE), , drop = FALSE]
          } else if ("fdr" %in% colnames(d)) {
            d <- d[order(d$fdr, na.last = TRUE), , drop = FALSE]
          }
          rownames(d) <- d$feature_id
          d
        }),
        contrasts
      )
      private$.de_dfs <- de_dfs

      # proteomics_choices: gene/protein lookup maps used by the
      # proteomics-only feature picker UI. protein_Id will have been dropped
      # already if it was redundant.
      sa1 <- de_dfs[[1]]
      gene_lookup <- NULL
      gene_choices <- NULL
      protein_lookup <- NULL
      protein_choices <- NULL
      if ("gene_name" %in% colnames(sa1) && !identical(sa1$feature_id, sa1$gene_name)) {
        gene_lookup <- split(sa1$feature_id, sa1$gene_name)
        gene_choices <- unique(sa1$gene_name)
      }
      if ("protein_Id" %in% colnames(sa1)) {
        protein_lookup <- split(sa1$feature_id, sa1$protein_Id)
        protein_choices <- unique(sa1$protein_Id)
      }
      private$.config$proteomics <- list(
        gene = list(choices = gene_choices, lookup = gene_lookup),
        protein = list(choices = protein_choices, lookup = protein_lookup)
      )
    }
  )
)
