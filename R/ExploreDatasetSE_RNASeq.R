#' Concrete subclass for SUSHI RNA-seq `SummarizedExperiment` outputs
#'
#' Reads from a SUSHI ezApp DESeq2 / edgeR / limma result. Top-level
#' `rowData(se)` carries feature annotations plus inline DE columns for
#' the single contrast described by `param$comparison`. Optional
#' pre-computed `enrichResult` / `gseaResult` objects live in
#' `metadata(se)`.
#'
#' Inherits the accessor contract from [ExploreDataset] via
#' [ExploreDatasetSE]; overrides only `.build_config()` and `.prep()`.
#'
#' @export
ExploreDatasetSE_RNASeq <- R6::R6Class(
  classname = "ExploreDatasetSE_RNASeq",
  inherit   = ExploreDatasetSE,

  private = list(

    .build_config = function(se) {
      param <- S4Vectors::metadata(se)$param %||% list()
      cd    <- data.frame(SummarizedExperiment::colData(se),
                          check.names = FALSE, stringsAsFactors = FALSE)

      factor_cols  <- grep("\\[Factor]", colnames(cd), value = TRUE)
      factor_names <- gsub(" \\[.*", "",  factor_cols)
      factor_names <- gsub(" ",      "_", factor_names)
      factor_levels <- setNames(
        lapply(factor_cols, function(fc) {
          v <- cd[[fc]]
          if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
        }),
        factor_names
      )

      comp <- param$comparison %||% "test--over--reference"
      parts <- strsplit(comp, "--over--", fixed = TRUE)[[1]]
      test_grp <- parts[1] %||% "test"
      ref_grp  <- parts[2] %||% "reference"
      contrast_name <- paste0(eds_safe(test_grp), "_vs_", eds_safe(ref_grp))

      cm <- list()
      cm[[contrast_name]] <- list(
        test_group      = test_grp,
        reference_group = ref_grp,
        factor_name     = param$groupingName,
        method          = "DESeq2",
        thresholds      = list(
          fdr    = param$fdrThreshORA      %||% 0.05,
          log2fc = param$log2RatioThreshGO %||% 0,
          go_p   = param$pValThreshGO      %||% 0.01
        )
      )

      enrichment_kinds <- character()
      md <- S4Vectors::metadata(se)
      if (isTRUE(param$runGO)) {
        if (!is.null(md$enrichResult))     enrichment_kinds <- c(enrichment_kinds, "ora")
        if (!is.null(md$enrichResultGSEA)) enrichment_kinds <- c(enrichment_kinds, "gsea")
      }

      norm_key <- intersect(c("signal", "xNorm", "normalizedCounts"),
                            SummarizedExperiment::assayNames(se))
      norm_key <- if (length(norm_key)) norm_key[1] else NA_character_
      has_refbuild <- !is.null(param$refBuild) && nzchar(as.character(param$refBuild))

      # Assays exposed by this subclass. NULL value = computed in
      # ExploreDatasetSE$assay() — see that method's switch().
      assay_map_l <- c(
        list(counts = "counts"),
        if (!is.na(norm_key)) list(normalised = norm_key, normalised_log2 = NULL),
        list(vst = NULL),
        if (has_refbuild) list(fpkm = NULL, tpm = NULL)
      )
      assay_labels_l <- c(
        counts          = "Raw",
        normalised      = "Normalised",
        normalised_log2 = "Normalised + Log2",
        vst             = "VST",
        fpkm            = "FPKM",
        tpm             = "TPM"
      )

      list(
        schema_version        = "exploreDE/config/1.0",
        omics_class           = "genomics",
        data_type             = "RNASeq",
        organism              = eds_infer_organism(param),
        ref_build             = param$refBuild %||% NA_character_,
        feature_id_kind       = param$featureLevel %||% "gene",
        factors               = factor_names,
        factor_levels         = factor_levels,
        primary_grouping      = param$groupingName %||%
          (if (length(factor_names)) factor_names[[1]] else NA_character_),
        secondary_grouping    = if (nzchar(param$grouping2Name %||% "")) {
          param$grouping2Name
        } else NA_character_,
        assay_map             = assay_map_l,
        assay_labels          = assay_labels_l,
        background_expression = param$backgroundExpression %||% 10,
        raw_assay_name        = "counts",
        primary_assay_name    = if (!is.na(norm_key)) "normalised_log2" else "counts",
        display_name_field    = "gene_name",
        description_field     = "description",
        primary_contrast      = contrast_name,
        contrast_meta         = cm,
        enrichment_kinds      = enrichment_kinds,
        enrichment_locator    = list(
          ora  = list(se = "enrichResult"),
          gsea = list(se = "enrichResultGSEA")
        ),
        sgrna_to_gene         = NULL,
        proteomics            = NULL,
        source = list(
          backend_kind   = "se",
          schema_version = "exploreDE/config/1.0",
          built_at       = Sys.time()
        )
      )
    },

    .prep = function() {
      se  <- private$.se
      cfg <- private$.config

      rd <- data.frame(SummarizedExperiment::rowData(se),
                       check.names = FALSE, stringsAsFactors = FALSE)
      rownames(rd) <- rownames(se)

      # Build feature_id via scuttle::uniquifyFeatureNames so the new API
      # matches the legacy convert_genomics_se() convention AND the
      # canonical AnnData spec (var_names = uniquified gene symbol).
      is_isoform <- identical(cfg$feature_id_kind, "isoform")
      rd$feature_id <- if (is_isoform && "transcript_id" %in% colnames(rd)) {
        scuttle::uniquifyFeatureNames(ID = rd$transcript_id, names = rd$gene_name)
      } else if (all(c("gene_id", "gene_name") %in% colnames(rd))) {
        scuttle::uniquifyFeatureNames(ID = rd$gene_id, names = rd$gene_name)
      } else {
        rownames(rd)
      }
      # Legacy underscore cleanup.
      leading_us <- grep("^_", rd$feature_id)
      if (length(leading_us)) {
        rd$feature_id[leading_us] <- gsub("_", "", rd$feature_id[leading_us])
      }

      rd$display_name <- if (cfg$display_name_field %in% colnames(rd)) {
        rd[[cfg$display_name_field]]
      } else {
        rd$feature_id
      }
      blank <- is.na(rd$display_name) | !nzchar(rd$display_name)
      rd$display_name[blank] <- rd$feature_id[blank]
      if ("description" %in% colnames(rd)) {
        rd$description <- gsub(" \\[.*", "", as.character(rd$description))
      }
      rownames(rd) <- rd$feature_id

      # fid_lookup: SE rowname → new feature_id. Built BEFORE sorting so the
      # vector alignment (rd[i]$feature_id ↔ rownames(se)[i]) is preserved.
      # Legacy convert_genomics_se applies a `[:/()-]` → `.` sanitisation
      # ONLY to assay rownames at isoform level, leaving seqAnno$feature_id
      # unsanitised — an internal inconsistency. The new API doesn't
      # replicate that quirk; assay() rownames match features()$feature_id,
      # so modules joining the two no longer need to sanitise manually.
      private$.fid_lookup <- setNames(rd$feature_id, rownames(se))

      # Sort by p-value where available — matches legacy convert_genomics_se,
      # which sorts gs$seqAnno (and therefore gs$countList rownames and
      # gs$genes) by pValue ascending. Users see most-significant features
      # first in autocompletes and top-N filters.
      pval_col <- grep("^p[._]?value$", colnames(rd),
                       value = TRUE, ignore.case = TRUE)
      if (length(pval_col) >= 1) {
        rd <- rd[order(rd[[pval_col[1]]], na.last = TRUE), , drop = FALSE]
      }
      private$.features_df <- rd

      # RNA-seq has DE columns inline in rowData (one contrast).
      contrast <- cfg$primary_contrast
      private$.de_dfs <- setNames(
        list(eds_standardise_de_columns(rd)),
        contrast
      )
    }
  )
)
