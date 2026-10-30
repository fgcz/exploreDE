#' Abstract intermediate class for SE-backed datasets
#'
#' `ExploreDatasetSE` holds the `SummarizedExperiment` backend
#' reference and implements the accessor contract uniformly for every
#' SE shape. Concrete subclasses (`ExploreDatasetSE_RNASeq`,
#' `_CountQC`, `_CrisprCountQC`, `_Proteomics`) override two private
#' methods only:
#'
#' - `.build_config(se)` — returns the per-instance config list at
#'   `$initialize()` time
#' - `.prep()` — populates `private$.features_df`, `private$.de_dfs`,
#'   `private$.fid_lookup` from the SE
#'
#' This class is not instantiable on its own — instantiate a concrete
#' subclass instead.
#'
#' @export
ExploreDatasetSE <- R6::R6Class(
  classname = "ExploreDatasetSE",
  inherit = ExploreDataset,

  public = list(
    initialize = function(se, config = NULL) {
      stopifnot(inherits(se, "SummarizedExperiment"))
      private$.se <- se
      private$.config <- config %||% private$.build_config(se)
      private$.prep()
      invisible(self)
    },

    backend = function() private$.se,
    backend_kind = function() "se",

    # samples(): uniform colData read. SUSHI's `[Factor]` suffix is stripped
    # from colData column names (genomics convention). Proteomics G_-prefixed
    # column names are LEFT AS-IS to match the legacy convert_proteomics_se()
    # output that server modules currently consume.
    #
    # A `Name` column is added matching the legacy converters: factor-typed
    # for genomics, plain character for proteomics. Rownames are colnames(se)
    # where available, else cd$Name / cd$name fallback.
    samples = function() {
      cd <- data.frame(SummarizedExperiment::colData(private$.se), check.names = FALSE, stringsAsFactors = FALSE)
      colnames(cd) <- gsub(" \\[.*", "", colnames(cd))
      colnames(cd) <- gsub(" ", "_", colnames(cd))

      se_cn <- colnames(private$.se)
      sample_ids <- if (!is.null(se_cn) && length(se_cn) == nrow(cd)) {
        se_cn
      } else {
        id_col <- intersect(c("Sample", "Name", "sampleName", "name"), colnames(cd))
        if (length(id_col)) as.character(cd[[id_col[1]]]) else NULL
      }

      if (!is.null(sample_ids)) {
        rownames(cd) <- sample_ids
        oc <- private$.config$omics_class
        if (identical(oc, "genomics")) {
          # Legacy genomics: factor()-ize Name with explicit level order.
          # If colData lacks a Name column, create one from sample_ids.
          name_vals <- if ("Name" %in% colnames(cd)) {
            as.character(cd$Name)
          } else {
            sample_ids
          }
          cd$Name <- factor(name_vals, levels = name_vals)
        } else if (identical(oc, "proteomics")) {
          # Legacy proteomics: plain-char Name from Sample/Name/sampleName,
          # plus factor()-ize every factor column with mixedsort-derived
          # levels.
          if (!"Name" %in% colnames(cd)) {
            cd$Name <- sample_ids
          }
          factor_cols <- private$.config$factors %||% character()
          for (fc in intersect(factor_cols, colnames(cd))) {
            if (!is.factor(cd[[fc]])) {
              cd[[fc]] <- factor(
                cd[[fc]],
                levels = unique(gtools::mixedsort(cd[[fc]]))
              )
            }
          }
        }
      }
      cd
    },

    features = function() private$.features_df,

    de = function(contrast) {
      stopifnot(contrast %in% self$contrasts())
      private$.de_dfs[[contrast]] %||% data.frame()
    },

    assay_names = function() {
      backend_keys <- SummarizedExperiment::assayNames(private$.se)
      keys <- names(private$.config$assay_map %||% list())
      keys[vapply(
        keys,
        function(k) {
          v <- private$.config$assay_map[[k]]
          is.null(v) || v %in% backend_keys
        },
        logical(1)
      )]
    },

    assay = function(name, features = NULL, samples = NULL) {
      v <- private$.config$assay_map[[name]]
      m <- if (!is.null(v)) {
        # Backend assay — remap via fid_lookup.
        .eds_remap_rownames(
          SummarizedExperiment::assay(private$.se, v),
          private$.fid_lookup
        )
      } else {
        # Computed assay.
        switch(
          name,
          # normalised_log2 recurses into self$assay("normalised") which is
          # ALREADY remapped — don't re-apply the lookup.
          normalised_log2 = log2(
            self$assay("normalised") +
              (private$.config$background_expression %||% 10)
          ),
          # vst / fpkm / tpm operate on the raw SE — output has SE rownames,
          # so apply fid_lookup remap before returning.
          vst = .eds_remap_rownames(eds_compute_vst(private$.se), private$.fid_lookup),
          fpkm = .eds_remap_rownames(eds_compute_rpkm(private$.se), private$.fid_lookup),
          tpm = .eds_remap_rownames(eds_compute_tpm(private$.se), private$.fid_lookup),
          stop("Unknown computed assay: ", name)
        )
      }
      if (!is.null(features)) {
        m <- m[features, , drop = FALSE]
      }
      if (!is.null(samples)) {
        m <- m[, samples, drop = FALSE]
      }
      m
    },

    enrichment = function(contrast, kind) {
      stopifnot(kind %in% c("ora", "gsea"))
      loc <- private$.config$enrichment_locator[[kind]]$se
      if (is.null(loc)) {
        return(list())
      }
      obj <- S4Vectors::metadata(private$.se)[[loc]]
      if (is.null(obj)) {
        return(list())
      }
      eds_flatten_cp_results(obj)
    },

    enrichment_input = function(contrast) {
      ei <- S4Vectors::metadata(private$.se)$enrichInput
      if (is.null(ei)) {
        return(list(seq_anno = NULL, ranked_l2fc = NULL, selections = list()))
      }
      list(
        seq_anno = ei$seqAnno,
        ranked_l2fc = ei$log2Ratio,
        selections = ei$selections %||% list()
      )
    }
  ),

  private = list(
    .se = NULL,
    .features_df = NULL,
    .de_dfs = list(),
    .fid_lookup = NULL, # NULL → identity inside assay()

    # Subclasses must override these:
    .build_config = function(se) {
      stop("subclass of ExploreDatasetSE must implement .build_config(se)")
    },
    .prep = function() {
      stop("subclass of ExploreDatasetSE must implement .prep()")
    }
  )
)
