#' Concrete subclass for CRISPR Mageck CountQC outputs
#'
#' Identical to [ExploreDatasetSE_CountQC] plus a sgRNA → gene mapping
#' parsed from the SE feature IDs, which are assumed to follow the
#' `GENE_<n>_<SEQUENCE>` convention.
#'
#' Inherits `.prep()` from [ExploreDatasetSE_CountQC] unchanged.
#'
#' @return An R6 class generator; `$new()` returns an
#'   `ExploreDatasetSE_CrisprCountQC` object.
#' @examples
#' library(SummarizedExperiment)
#' # Guide counts: 4 guides for each of 10 genes, named GENE_<n>_<SEQUENCE>
#' set.seed(1)
#' guides <- sprintf("GENE%d_%d_ACGTACGTACGT", rep(1:10, each = 4), rep(1:4, 10))
#' counts <- matrix(rpois(40 * 4, 50), nrow = 40, dimnames = list(guides, paste0("S", 1:4)))
#' samples <- S4Vectors::DataFrame(
#'   `Condition [Factor]` = rep(c("ctrl", "treated"), each = 2),
#'   row.names = colnames(counts),
#'   check.names = FALSE
#' )
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = samples)
#' eds <- ExploreDatasetSE_CrisprCountQC$new(se)
#' eds
#' head(eds$sgrna_to_gene())
#' @export
ExploreDatasetSE_CrisprCountQC <- R6::R6Class(
  classname = "ExploreDatasetSE_CrisprCountQC",
  inherit = ExploreDatasetSE_CountQC,

  private = list(
    .build_config = function(se) {
      cfg <- super$.build_config(se)
      cfg$data_type <- "CrisprCountQC"
      cfg$feature_id_kind <- "guide"
      rn <- rownames(se) %||% character()
      # Strip the guide suffix to recover the gene. Two conventions supported,
      # both ending in the protospacer sequence:
      #   GENE_<n>_<SEQ>              (e.g. Actb_1_ACGT… / Actb_sgrna1_ACGT…)
      #   GENE_<strand>_<locus>_<SEQ> (e.g. Ptplad2_-_88439147.23-P1P2_GATCC…)
      cfg$sgrna_to_gene <- setNames(
        gsub("_((sgrna)?[0-9]+|[+-]_[^_]+)_[ACGTN]+$", "", rn, ignore.case = TRUE),
        rn
      )
      cfg
    }
    # .prep() inherited from _CountQC unchanged
  )
)
