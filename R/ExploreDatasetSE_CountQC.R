#' Concrete subclass for SUSHI CountQC `SummarizedExperiment` outputs
#'
#' Top-level `rowData(se)` carries feature annotations; there are no
#' DE results. May include a `presentFlag` / `isPresent` assay used
#' as the QC denominator.
#'
#' Inherits the accessor contract from [ExploreDataset] via
#' [ExploreDatasetSE]; overrides only `.build_config()` and `.prep()`.
#'
#' @export
ExploreDatasetSE_CountQC <- R6::R6Class(
  classname = "ExploreDatasetSE_CountQC",
  inherit = ExploreDatasetSE,

  private = list(
    .build_config = function(se) {
      param <- S4Vectors::metadata(se)$param %||% list()
      cd <- data.frame(SummarizedExperiment::colData(se), check.names = FALSE, stringsAsFactors = FALSE)

      factor_cols <- grep("\\[Factor]", colnames(cd), value = TRUE)
      factor_names <- gsub(" \\[.*", "", factor_cols)
      factor_names <- gsub(" ", "_", factor_names)
      factor_levels <- setNames(
        lapply(factor_cols, function(fc) {
          v <- cd[[fc]]
          if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
        }),
        factor_names
      )

      ass_nms <- SummarizedExperiment::assayNames(se)
      norm_key <- intersect(c("signal", "xNorm"), ass_nms)
      norm_key <- if (length(norm_key)) norm_key[1] else NULL
      pf_key <- intersect(c("presentFlag", "isPresent"), ass_nms)
      pf_key <- if (length(pf_key)) pf_key[1] else NULL
      has_refbuild <- !is.null(param$refBuild) && nzchar(as.character(param$refBuild))

      assay_map_l <- c(
        list(counts = "counts"),
        if (!is.null(norm_key)) list(normalised = norm_key, normalised_log2 = NULL),
        if (!is.null(pf_key)) list(present_flag = pf_key),
        list(vst = NULL),
        if (has_refbuild) list(fpkm = NULL, tpm = NULL)
      )
      assay_labels_l <- c(
        counts = "Raw",
        normalised = "Normalised",
        normalised_log2 = "Normalised + Log2",
        present_flag = "Present Flag",
        vst = "VST",
        fpkm = "FPKM",
        tpm = "TPM"
      )

      list(
        schema_version = "exploreDE/config/1.0",
        omics_class = "genomics",
        data_type = "CountQC",
        organism = eds_infer_organism(param),
        ref_build = param$refBuild %||% NA_character_,
        feature_id_kind = param$featureLevel %||% "gene",
        factors = factor_names,
        factor_levels = factor_levels,
        primary_grouping = param$groupingName %||%
          (if (length(factor_names)) factor_names[[1]] else NA_character_),
        secondary_grouping = if (nzchar(param$grouping2Name %||% "")) {
          param$grouping2Name
        } else {
          NA_character_
        },
        assay_map = assay_map_l,
        assay_labels = assay_labels_l,
        background_expression = param$backgroundExpression %||% 10,
        raw_assay_name = "counts",
        primary_assay_name = "counts",
        display_name_field = "gene_name",
        description_field = "description",
        primary_contrast = NA_character_,
        contrast_meta = list(),
        enrichment_kinds = character(),
        enrichment_locator = list(),
        sgrna_to_gene = NULL,
        proteomics = NULL,
        source = list(
          backend_kind = "se",
          schema_version = "exploreDE/config/1.0",
          built_at = Sys.time()
        )
      )
    },

    .prep = function() {
      se <- private$.se
      cfg <- private$.config

      rd <- data.frame(SummarizedExperiment::rowData(se), check.names = FALSE, stringsAsFactors = FALSE)
      rownames(rd) <- rownames(se)

      # Build feature_id via scuttle::uniquifyFeatureNames where gene_id +
      # gene_name are present (genomics CountQC). CRISPR CountQC inherits
      # this .prep() but its rowData has no gene_id/gene_name, so feature_id
      # falls through to rownames(se) (guide IDs).
      rd$feature_id <- if (all(c("gene_id", "gene_name") %in% colnames(rd))) {
        scuttle::uniquifyFeatureNames(ID = rd$gene_id, names = rd$gene_name)
      } else {
        rownames(rd)
      }
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
      rownames(rd) <- rd$feature_id
      private$.features_df <- rd

      private$.de_dfs <- list() # no DE
      private$.fid_lookup <- setNames(rd$feature_id, rownames(se))
    }
  )
)
