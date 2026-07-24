#' Abstract base class for exploreDE's storage-agnostic dataset wrapper
#'
#' `ExploreDataset` defines the accessor contract that every concrete
#' subclass must implement. It is not instantiable on its own — calling
#' `ExploreDataset$new()` raises an error. Use [new_explore_dataset()] to
#' build a concrete instance, or instantiate a subclass directly.
#'
#' Methods implemented on the base class read only from
#' `private$.config` and are inherited unchanged by every subclass.
#' Methods that touch the storage backend are virtual
#' (`stop("not implemented")`) and must be overridden.
#'
#' See `plans/anndata_migration/internal_api_spec_2026-04-30_thin.qmd`
#' for the full design and accessor contract.
#'
#' @export
ExploreDataset <- R6::R6Class(
  classname = "ExploreDataset",

  public = list(

    initialize = function(...) {
      stop("ExploreDataset is abstract. ",
           "Use new_explore_dataset() or one of the concrete subclasses.")
    },

    # ---- identity (config-only) -------------------------------------------
    omics_class      = function() private$.config$omics_class,
    data_type        = function() private$.config$data_type,
    has_de           = function() length(self$contrasts()) > 0L,
    has_enrichment   = function() {
      length(private$.config$enrichment_kinds %||% character()) > 0L
    },
    organism         = function() private$.config$organism %||% NA_character_,
    ref_build        = function() private$.config$ref_build %||% NA_character_,
    feature_id_kind  = function() private$.config$feature_id_kind,
    background_expression = function() {
      private$.config$background_expression %||% 10
    },
    backend_kind     = function() stop("not implemented"),
    backend          = function() stop("not implemented"),

    # ---- samples (samples() virtual; rest derived) ------------------------
    samples          = function() stop("not implemented"),
    sample_ids       = function() rownames(self$samples()),
    n_samples        = function() nrow(self$samples()),
    factor_names     = function() private$.config$factors %||% character(),
    factor_levels    = function(factor) {
      private$.config$factor_levels[[factor]]
    },
    primary_grouping = function() private$.config$primary_grouping,
    # NA when the dataset has no secondary grouping (e.g. proteomics, or
    # SUSHI params with grouping2Name blank).
    secondary_grouping = function() private$.config$secondary_grouping %||% NA_character_,
    sample_factor    = function(factor) self$samples()[[factor]],

    # ---- features (features() virtual; rest derived) ----------------------
    features         = function() stop("not implemented"),
    feature_ids      = function() rownames(self$features()),
    n_features       = function() nrow(self$features()),
    feature_display  = function(ids = NULL) {
      f <- self$features()
      if (is.null(ids)) ids <- rownames(f)
      setNames(f[ids, "display_name"], ids)
    },

    # ---- assays (assay_names + assay virtual; rest derived) ---------------
    assay_names        = function() stop("not implemented"),
    has_assay          = function(name) name %in% self$assay_names(),
    assay              = function(name, features = NULL, samples = NULL) {
      stop("not implemented")
    },
    # Human-readable label for an assay (used by UI dropdowns).
    # Subclass populates private$.config$assay_labels at build_config() time.
    # Falls back to the canonical name if no label is registered.
    assay_label        = function(name) {
      labels <- private$.config$assay_labels %||% list()
      if (name %in% names(labels)) return(labels[[name]])
      name
    },
    raw_assay_name     = function() private$.config$raw_assay_name,
    primary_assay_name = function() private$.config$primary_assay_name,
    raw_assay          = function(...) self$assay(self$raw_assay_name(), ...),
    primary_assay      = function(...) self$assay(self$primary_assay_name(), ...),

    # ---- DE (de() virtual; rest config-only) ------------------------------
    contrasts        = function() {
      names(private$.config$contrast_meta) %||% character()
    },
    primary_contrast = function() private$.config$primary_contrast,
    de               = function(contrast) stop("not implemented"),
    de_all           = function() {
      setNames(lapply(self$contrasts(), self$de), self$contrasts())
    },
    contrast_meta    = function(contrast) {
      stopifnot(contrast %in% self$contrasts())
      private$.config$contrast_meta[[contrast]]
    },
    # Display label for a contrast. Today the canonical contrast key is
    # already in "test_vs_reference" identifier-safe form (built that way in
    # _RNASeq$.build_config; supplied directly by prolfqua for proteomics),
    # so this helper just validates and returns it. Modules wanting prose
    # ("diff vs undiff") for plot titles can `gsub("_vs_", " vs ", x)`.
    contrast_label   = function(contrast) {
      stopifnot(contrast %in% self$contrasts())
      contrast
    },

    # ---- enrichment (all virtual except enrichment_kinds) -----------------
    enrichment_kinds = function() {
      private$.config$enrichment_kinds %||% character()
    },
    enrichment       = function(contrast, kind) stop("not implemented"),
    enrichment_input = function(contrast)        stop("not implemented"),

    # ---- domain extras (config-only) --------------------------------------
    sgrna_to_gene    = function() private$.config$sgrna_to_gene,
    proteomics_choices = function(kind = c("gene", "protein")) {
      kind <- match.arg(kind)
      private$.config$proteomics[[kind]]
    },

    # ---- GO pathways (derived from features()) ----------------------------
    # Returns character vector of "GO:XXXXXXX term name" strings for use in
    # selectize dropdowns (server-goTile.R, server-hypergeo2.R). Derived from
    # the GO BP / GO MF / GO CC columns in features(), so works for any
    # backend. Returns character(0) when those columns are absent (proteomics,
    # CountQC without GO annotations).
    all_pathways = function() {
      f       <- self$features()
      go_cols <- intersect(c("GO BP", "GO MF", "GO CC"), colnames(f))
      if (!length(go_cols)) return(character(0))
      ids <- unique(unlist(lapply(go_cols, function(col) {
        x <- f[[col]]
        x <- x[!is.na(x) & nzchar(x)]
        strsplit(x, "; ", fixed = TRUE)
      }), use.names = FALSE))
      ids <- ids[nzchar(ids)]
      if (!length(ids)) return(character(0))
      terms <- clusterProfiler::go2term(ids)
      if (!nrow(terms)) return(character(0))
      paste(terms[[1]], terms[[2]])
    },

    # ---- provenance --------------------------------------------------------
    source           = function() private$.config$source,

    # ---- lifecycle ---------------------------------------------------------
    print = function(...) {
      cat("<", class(self)[1], "  (", self$backend_kind(), ")>\n", sep = "")
      cat("  data_type:        ", self$data_type(), "\n", sep = "")
      cat("  feature_id_kind:  ", self$feature_id_kind(), "\n", sep = "")
      cat("  organism:         ", self$organism(), "\n", sep = "")
      cat("  ref_build:        ", self$ref_build(), "\n", sep = "")
      cat("  samples x features: ",
          self$n_samples(), " x ", self$n_features(), "\n", sep = "")
      cat("  factors (", length(self$factor_names()), "):\n", sep = "")
      for (f in self$factor_names()) {
        lv <- self$factor_levels(f)
        cat("    ", f,
            if (identical(f, self$primary_grouping())) "  (primary)" else "",
            ": ",
            paste(utils::head(lv, 6), collapse = ", "),
            if (length(lv) > 6) paste0(", ... (+", length(lv) - 6, ")") else "",
            "\n", sep = "")
      }
      cat("  contrasts (", length(self$contrasts()), "): ",
          paste(utils::head(self$contrasts(), 3), collapse = ", "),
          if (length(self$contrasts()) > 3) ", ..." else "",
          "\n", sep = "")
      cat("  primary_contrast: ",
          self$primary_contrast() %||% "(none)", "\n", sep = "")
      cat("  assays:           ",
          paste(self$assay_names(), collapse = ", "), "\n", sep = "")
      cat("  enrichment:       ",
          if (self$has_enrichment())
            paste(self$enrichment_kinds(), collapse = ", ")
          else "(none)",
          "\n", sep = "")
      invisible(self)
    }
  ),

  private = list(
    .config = NULL
  )
)
