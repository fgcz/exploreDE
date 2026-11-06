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
#' @return An R6 class generator. `$new()` always errors; the concrete
#'   subclasses return an object inheriting from `ExploreDataset`.
#' @examples
#' se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
#' eds <- new_explore_dataset(se)
#' inherits(eds, "ExploreDataset")
#' eds$contrasts()
#' head(eds$de(eds$primary_contrast()))
#'
#' # The base class itself cannot be instantiated
#' try(ExploreDataset$new())
#' @export
ExploreDataset <- R6::R6Class(
  classname = "ExploreDataset",

  public = list(
    #' @description Always errors: the base class is abstract. Use [new_explore_dataset()] or a concrete subclass.
    #' @param ... Ignored.
    initialize = function(...) {
      stop("ExploreDataset is abstract. ", "Use new_explore_dataset() or one of the concrete subclasses.")
    },

    # ---- identity (config-only) -------------------------------------------
    #' @description Omics class from the config: `"genomics"` or `"proteomics"`.
    omics_class = function() private$.config$omics_class,
    #' @description Dataset type from the config, e.g. `"RNASeq"`, `"CountQC"`, `"CrisprCountQC"` or `"proteomics"`.
    data_type = function() private$.config$data_type,
    #' @description `TRUE` when the dataset has at least one contrast.
    has_de = function() length(self$contrasts()) > 0L,
    #' @description `TRUE` when the config lists at least one enrichment kind.
    has_enrichment = function() {
      length(private$.config$enrichment_kinds %||% character()) > 0L
    },
    #' @description Organism from the config, or `NA` when unknown.
    organism = function() private$.config$organism %||% NA_character_,
    #' @description Reference genome build from the config, or `NA` when unknown.
    ref_build = function() private$.config$ref_build %||% NA_character_,
    #' @description What a feature ID denotes, from the config, e.g. `"gene"`, `"isoform"`, `"protein"`,
    #'   `"peptide"`, `"phospho"` or `"guide"`.
    feature_id_kind = function() private$.config$feature_id_kind,
    #' @description Low-signal offset added to expression values before taking log2; from the config, default 10.
    background_expression = function() {
      private$.config$background_expression %||% 10
    },
    #' @description Short tag naming the storage backend (`"se"`, `"h5ad"` or `"parquet"`). Implemented by subclasses.
    backend_kind = function() stop("not implemented"),
    #' @description The underlying storage object. Implemented by subclasses.
    backend = function() stop("not implemented"),

    # ---- samples (samples() virtual; rest derived) ------------------------
    #' @description Sample metadata as a data.frame, one row per sample, rownames are sample IDs.
    #'   Implemented by subclasses.
    samples = function() stop("not implemented"),
    #' @description Sample IDs: the rownames of `samples()`.
    sample_ids = function() rownames(self$samples()),
    #' @description Number of samples: the number of rows of `samples()`.
    n_samples = function() nrow(self$samples()),
    #' @description Names of the experimental factors from the config; empty when there are none.
    factor_names = function() private$.config$factors %||% character(),
    #' @description Levels of one experimental factor, from the config.
    #' @param factor Factor name, one of `factor_names()`.
    factor_levels = function(factor) {
      private$.config$factor_levels[[factor]]
    },
    #' @description Name of the factor used as the primary sample grouping, from the config.
    primary_grouping = function() private$.config$primary_grouping,
    # NA when the dataset has no secondary grouping (e.g. proteomics, or
    # SUSHI params with grouping2Name blank).
    #' @description Name of the secondary grouping factor, or `NA` when there is none.
    secondary_grouping = function() private$.config$secondary_grouping %||% NA_character_,
    #' @description One column of `samples()`: the values of a factor, one per sample.
    #' @param factor Column name in `samples()`, usually one of `factor_names()`.
    sample_factor = function(factor) self$samples()[[factor]],

    # ---- features (features() virtual; rest derived) ----------------------
    #' @description Feature annotation as a data.frame, one row per feature, rownames are feature IDs, with a
    #'   `display_name` column. Implemented by subclasses.
    features = function() stop("not implemented"),
    #' @description Feature IDs: the rownames of `features()`.
    feature_ids = function() rownames(self$features()),
    #' @description Number of features: the number of rows of `features()`.
    n_features = function() nrow(self$features()),
    #' @description Display names of features, named by feature ID.
    #' @param ids Feature IDs to look up; `NULL` (default) returns all features.
    feature_display = function(ids = NULL) {
      f <- self$features()
      if (is.null(ids)) {
        ids <- rownames(f)
      }
      setNames(f[ids, "display_name"], ids)
    },

    # ---- assays (assay_names + assay virtual; rest derived) ---------------
    #' @description Canonical names of the assays this dataset provides. Implemented by subclasses.
    assay_names = function() stop("not implemented"),
    #' @description `TRUE` when `name` is one of `assay_names()`.
    #' @param name Canonical assay name.
    has_assay = function(name) name %in% self$assay_names(),
    #' @description One assay as a features x samples matrix, optionally subset. Implemented by subclasses.
    #' @param name Canonical assay name, one of `assay_names()`.
    #' @param features Feature IDs (rows) to keep; `NULL` keeps all.
    #' @param samples Sample IDs (columns) to keep; `NULL` keeps all.
    assay = function(name, features = NULL, samples = NULL) {
      stop("not implemented")
    },
    # Human-readable label for an assay (used by UI dropdowns).
    # Subclass populates private$.config$assay_labels at build_config() time.
    # Falls back to the canonical name if no label is registered.
    #' @description Human-readable label for an assay, from the config; `name` itself when no label is registered.
    #' @param name Canonical assay name.
    assay_label = function(name) {
      labels <- private$.config$assay_labels %||% list()
      if (name %in% names(labels)) {
        return(labels[[name]])
      }
      name
    },
    #' @description Canonical name of the raw assay, from the config.
    raw_assay_name = function() private$.config$raw_assay_name,
    #' @description Canonical name of the primary assay, from the config.
    primary_assay_name = function() private$.config$primary_assay_name,
    #' @description The raw assay: shortcut for `assay(raw_assay_name(), ...)`.
    #' @param ... Passed to `assay()`: `features`, `samples`.
    raw_assay = function(...) self$assay(self$raw_assay_name(), ...),
    #' @description The primary assay: shortcut for `assay(primary_assay_name(), ...)`.
    #' @param ... Passed to `assay()`: `features`, `samples`.
    primary_assay = function(...) self$assay(self$primary_assay_name(), ...),

    # ---- DE (de() virtual; rest config-only) ------------------------------
    #' @description Names of the DE contrasts, from the config; empty when the dataset has no DE results.
    contrasts = function() {
      names(private$.config$contrast_meta) %||% character()
    },
    #' @description Name of the primary contrast, from the config; `NA` when there are no contrasts.
    primary_contrast = function() private$.config$primary_contrast,
    #' @description DE result table for one contrast, as a data.frame. Implemented by subclasses.
    #' @param contrast Contrast name, one of `contrasts()`.
    de = function(contrast) stop("not implemented"),
    #' @description DE result tables of all contrasts, as a list named by contrast.
    de_all = function() {
      setNames(lapply(self$contrasts(), self$de), self$contrasts())
    },
    #' @description Metadata list of one contrast from the config: test and reference group, factor, method,
    #'   thresholds. Errors on an unknown contrast.
    #' @param contrast Contrast name, one of `contrasts()`.
    contrast_meta = function(contrast) {
      stopifnot(contrast %in% self$contrasts())
      private$.config$contrast_meta[[contrast]]
    },
    # Display label for a contrast. Today the canonical contrast key is
    # already in "test_vs_reference" identifier-safe form (built that way in
    # _RNASeq$.build_config; supplied directly by prolfqua for proteomics),
    # so this helper just validates and returns it. Modules wanting prose
    # ("diff vs undiff") for plot titles can `gsub("_vs_", " vs ", x)`.
    #' @description Display label of a contrast: currently the contrast name itself. Errors on an unknown contrast.
    #' @param contrast Contrast name, one of `contrasts()`.
    contrast_label = function(contrast) {
      stopifnot(contrast %in% self$contrasts())
      contrast
    },

    # Extra DE columns to surface for display only (beyond log2Ratio/pValue/fdr)
    # — e.g. SAINTexpress's SaintScore. Default none; backends declaring a
    # mapping (H5AD) override via config. Names are also the DE-table labels.
    #' @description Names of extra DE columns shown for display only (e.g. `SaintScore`), from the config;
    #'   empty by default.
    de_display_extras = function() private$.config$de_display_extras %||% character(),

    # Does any contrast carry a usable RAW p-value? FALSE for FDR-only producers
    # (SAINTexpress), which the p-value-type selectors use to degrade to
    # FDR-only rather than plotting all-NA -log10(pValue).
    #' @description `TRUE` when the DE table of at least one contrast has a non-missing raw `pValue`;
    #'   `FALSE` for FDR-only producers such as SAINTexpress.
    has_pvalue = function() {
      if (!self$has_de()) {
        return(FALSE)
      }
      for (cc in self$contrasts()) {
        d <- tryCatch(self$de(cc), error = function(e) NULL)
        if (!is.null(d) && "pValue" %in% colnames(d) && any(!is.na(d$pValue))) {
          return(TRUE)
        }
      }
      FALSE
    },

    # ---- enrichment (all virtual except enrichment_kinds) -----------------
    #' @description Enrichment kinds available, from the config: any of `"ora"` and `"gsea"`.
    enrichment_kinds = function() {
      private$.config$enrichment_kinds %||% character()
    },
    #' @description Enrichment results of one contrast and kind, as a flat named list of clusterProfiler results.
    #'   Implemented by subclasses.
    #' @param contrast Contrast name, one of `contrasts()`.
    #' @param kind `"ora"` or `"gsea"`.
    enrichment = function(contrast, kind) stop("not implemented"),
    #' @description Inputs for enrichment of one contrast: `list(seq_anno, ranked_l2fc, selections)`.
    #'   Implemented by subclasses.
    #' @param contrast Contrast name, one of `contrasts()`.
    enrichment_input = function(contrast) stop("not implemented"),

    # Tidy enrichment as a plain data.frame — the form that does NOT require
    # clusterProfiler/DOSE to read, and the one exploreDElight uses (DOSE has
    # no wasm build). Every column the app displays lives here; only the three
    # set-level plots (cnetplot/gseaplot2/ridgeplot) still need enrichment().
    #
    # Default implementation derives it from enrichment() so it works on EVERY
    # file ever written, including those with no uns$enrichment block. Backends
    # holding a stored tidy mirror should override and prefer it.
    # See plans/exploredelight/enrichment_detach_plan.md.
    #' @description Enrichment results of one contrast and kind as one tidy data.frame, derived from
    #'   `enrichment()`; an empty data.frame when there are none.
    #' @param contrast Contrast name, one of `contrasts()`.
    #' @param kind `"ora"` or `"gsea"`.
    enrichment_table = function(contrast, kind) {
      objs <- self$enrichment(contrast, kind)
      if (!length(objs)) {
        return(data.frame())
      }
      # enrichment() names its leaves with the kind already prefixed
      # ("ora/BP/upGenes"); the writer's stored mirror does not ("BP/upGenes").
      # Strip it so both paths agree -- otherwise the same dataset yields
      # different `ontology` values depending on whether it happens to carry a
      # tidy block. `method` already says which kind it is.
      tabs <- lapply(names(objs), function(nm) {
        .enr_tidy_from_native(objs[[nm]], method = kind, prefix = sub(paste0("^", kind, "/"), "", nm))
      })
      tabs <- Filter(Negate(is.null), tabs)
      if (!length(tabs)) {
        return(data.frame())
      }
      plyr::rbind.fill(tabs)
    },

    # The p.adjust cut applied to the stored tidy mirror, or NA when the tables
    # were derived live and nothing was filtered. Callers that let a user raise
    # a GO threshold need this: without it an absent row is ambiguous between
    # "cut by the writer", "not significant" and "never run".
    #' @description The `p.adjust` cut-off applied to the stored enrichment tables, from the config;
    #'   `NA` when none was applied.
    enrichment_filter = function() {
      private$.config$enrichment_filter %||% NA_real_
    },

    # ---- domain extras (config-only) --------------------------------------
    #' @description Named character vector mapping sgRNA IDs to gene names (CRISPR datasets); `NULL` otherwise.
    sgrna_to_gene = function() private$.config$sgrna_to_gene,
    #' @description Choices and lookup for the proteomics feature picker, `list(choices, lookup)`; `NULL` when the
    #'   dataset has none.
    #' @param kind `"gene"` (default) or `"protein"`.
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
    #' @description GO terms named in the `GO BP`, `GO MF` and `GO CC` columns of `features()`, as
    #'   `"GO:ID term name"` strings; empty when those columns are absent.
    all_pathways = function() {
      f <- self$features()
      go_cols <- intersect(c("GO BP", "GO MF", "GO CC"), colnames(f))
      if (!length(go_cols)) {
        return(character(0))
      }
      ids <- unique(unlist(
        lapply(go_cols, function(col) {
          x <- f[[col]]
          x <- x[!is.na(x) & nzchar(x)]
          strsplit(x, "; ", fixed = TRUE)
        }),
        use.names = FALSE
      ))
      ids <- ids[nzchar(ids)]
      if (!length(ids)) {
        return(character(0))
      }
      terms <- clusterProfiler::go2term(ids)
      if (!nrow(terms)) {
        return(character(0))
      }
      paste(terms[[1]], terms[[2]])
    },

    # ---- provenance --------------------------------------------------------
    #' @description Provenance list from the config: backend kind, schema version and build time.
    source = function() private$.config$source,

    # ---- lifecycle ---------------------------------------------------------
    #' @description Print a summary of the dataset (type, dimensions, factors, contrasts, assays, enrichment)
    #'   and return the object invisibly.
    #' @param ... Ignored.
    print = function(...) {
      cat("<", class(self)[1], "  (", self$backend_kind(), ")>\n", sep = "")
      cat("  data_type:        ", self$data_type(), "\n", sep = "")
      cat("  feature_id_kind:  ", self$feature_id_kind(), "\n", sep = "")
      cat("  organism:         ", self$organism(), "\n", sep = "")
      cat("  ref_build:        ", self$ref_build(), "\n", sep = "")
      cat("  samples x features: ", self$n_samples(), " x ", self$n_features(), "\n", sep = "")
      cat("  factors (", length(self$factor_names()), "):\n", sep = "")
      for (f in self$factor_names()) {
        lv <- self$factor_levels(f)
        cat(
          "    ",
          f,
          if (identical(f, self$primary_grouping())) "  (primary)" else "",
          ": ",
          paste(utils::head(lv, 6), collapse = ", "),
          if (length(lv) > 6) paste0(", ... (+", length(lv) - 6, ")") else "",
          "\n",
          sep = ""
        )
      }
      cat(
        "  contrasts (",
        length(self$contrasts()),
        "): ",
        paste(utils::head(self$contrasts(), 3), collapse = ", "),
        if (length(self$contrasts()) > 3) ", ..." else "",
        "\n",
        sep = ""
      )
      cat("  primary_contrast: ", self$primary_contrast() %||% "(none)", "\n", sep = "")
      cat("  assays:           ", paste(self$assay_names(), collapse = ", "), "\n", sep = "")
      cat(
        "  enrichment:       ",
        if (self$has_enrichment()) {
          paste(self$enrichment_kinds(), collapse = ", ")
        } else {
          "(none)"
        },
        "\n",
        sep = ""
      )
      invisible(self)
    }
  ),

  private = list(
    .config = NULL
  )
)
