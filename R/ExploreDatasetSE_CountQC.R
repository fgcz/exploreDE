#' Concrete subclass for SUSHI CountQC `SummarizedExperiment` outputs
#'
#' Top-level `rowData(se)` carries feature annotations; there are no
#' DE results. May include a `presentFlag` / `isPresent` assay used
#' as the QC denominator.
#'
#' Inherits the accessor contract from [ExploreDataset] via
#' [ExploreDatasetSE]; overrides only `.build_config()` and `.prep()`.
#'
#' @return An R6 class generator; `$new()` returns an
#'   `ExploreDatasetSE_CountQC` object.
#' @examples
#' library(SummarizedExperiment)
#' de_result <- readRDS(system.file("examples/genomics/deResult.rds", package = "exploreDE"))
#' # A count table without DE results, as SUSHI's CountQC writes it
#' se <- SummarizedExperiment(
#'   assays = list(counts = assay(de_result, "counts")),
#'   colData = colData(de_result),
#'   rowData = rowData(de_result)[, c("gene_id", "gene_name", "description")]
#' )
#' eds <- ExploreDatasetSE_CountQC$new(se)
#' eds
#' eds$has_de()
#' @export
ExploreDatasetSE_CountQC <- R6::R6Class(
  classname = "ExploreDatasetSE_CountQC",
  inherit = ExploreDatasetSE,

  private = list(
    .build_config = function(se) {
      param <- S4Vectors::metadata(se)$param %||% list()
      factors <- eds_genomics_factors(se)
      factor_names <- factors$names
      factor_levels <- factors$levels

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
        assay_log2 = .genomics_assay_log2(names(assay_map_l)),
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

      # CRISPR Count QC inherits this .prep(); its rowData has no gene_id /
      # gene_name, so its feature ids are the guide ids (rownames).
      rd <- eds_genomics_features(se, cfg)
      is_isoform <- identical(cfg$feature_id_kind, "isoform")
      private$.features_df <- rd

      if (is_isoform && "gene_name" %in% colnames(rd)) {
        private$.config$proteomics <- eds_isoform_gene_choices(rd$feature_id, rd$gene_name)
      }

      private$.de_dfs <- list() # no DE
      private$.fid_lookup <- setNames(rd$feature_id, rownames(se))
    }
  )
)
