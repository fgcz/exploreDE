#' Concrete subclass for CRISPR Mageck CountQC outputs
#'
#' Identical to [ExploreDatasetSE_CountQC] plus a sgRNA → gene mapping
#' parsed from the SE feature IDs, which are assumed to follow the
#' `GENE_<n>_<SEQUENCE>` convention.
#'
#' Inherits `.prep()` from [ExploreDatasetSE_CountQC] unchanged.
#'
#' @export
ExploreDatasetSE_CrisprCountQC <- R6::R6Class(
  classname = "ExploreDatasetSE_CrisprCountQC",
  inherit   = ExploreDatasetSE_CountQC,

  private = list(

    .build_config = function(se) {
      cfg <- super$.build_config(se)
      cfg$data_type       <- "CrisprCountQC"
      cfg$feature_id_kind <- "guide"
      rn <- rownames(se) %||% character()
      cfg$sgrna_to_gene <- setNames(
        gsub("_(sgrna)?[0-9]+_[ACGTN]+$", "", rn, ignore.case = TRUE),
        rn
      )
      cfg
    }
    # .prep() inherited from _CountQC unchanged
  )
)
