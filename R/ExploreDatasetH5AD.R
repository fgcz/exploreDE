# =============================================================================
#  ExploreDatasetH5AD  —  ExploreDataset backend over a canonical AnnData
#
#  Peer of the ExploreDatasetSE_* subclasses (see the AnnData migration
#  README): the single backend for user-supplied / externalised data. Reads a
#  canonical exploreDE/1.0 .h5ad (as produced by build_explore_h5ad() or the
#  SUSHI -> AnnData converter) and satisfies the same accessor contract every
#  server module reads through.
#
#  Canonical layout consumed here:
#    X / layers$*        samples x features   (axis-flipped vs an SE!)
#    obs                 sample metadata; experimental factors prefixed "factor_"
#    var                 feature annotation; rownames = feature_id
#    varm$de_<slot>      effect_log2fc / score_pvalue / score_fdr / label_*
#    uns$omics_type / hierarchy / design / de / enrichment (+ enrichment_support) / ...
#
#  Cross-API equivalence with ExploreDatasetSE_RNASeq is asserted in
#  tests/testthat/test-explore_dataset-h5ad.R (assays agree to max|d| = 0).
# =============================================================================

# ---- internal helpers -------------------------------------------------------

# Strip the canonical "factor_" prefix so the backend exposes clean factor
# names ("condition") that stay consistent with samples() columns.
.h5ad_strip_factor <- function(x) sub("^factor_", "", x)

# Regex-safe identifier — mirrors eds_safe() so H5AD contrast keys match the
# SE subclass's "test_vs_reference" convention.
.h5ad_safe <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

# Spec section 7.1: schema_version, omics_type and feature_id_kind are MUST.
# Refuse a file without them instead of guessing: the old defaults loaded a
# proteomics file with a missing or misspelt omics_type as RNA-seq (intensities
# treated as counts) and any AnnData as exploreDE/1.0.
.h5ad_require_identity <- function(uns) {
  fail <- function(what) {
    stop(
      "not an exploreDE/1.0 .h5ad: ",
      what,
      " (see inst/spec/exploreDE-h5ad-1.0.md, section 7.1)",
      call. = FALSE
    )
  }
  one <- function(key) {
    v <- uns[[key]]
    if (is.null(v) || !length(v) || is.na(v[[1]]) || !nzchar(v[[1]])) {
      fail(sprintf("uns[\"%s\"] is missing", key))
    }
    as.character(v[[1]])
  }
  sv <- one("schema_version")
  if (sv != "exploreDE/1.0") {
    fail(sprintf("uns[\"schema_version\"] is \"%s\"", sv))
  }
  ot <- one("omics_type")
  if (!ot %in% c("proteomics", "rnaseq_bulk", "crispr_countqc")) {
    fail(sprintf("uns[\"omics_type\"] is \"%s\", not proteomics, rnaseq_bulk or crispr_countqc", ot))
  }
  one("feature_id_kind")
  invisible(TRUE)
}

# Feature (var) columns per role, the spec's label_* name first; files from
# other writers may carry the plain names. .h5ad_role_col() returns the first
# one present, or NULL.
.h5ad_role_cols <- list(
  symbol = c("label_gene_symbol", "gene_name", "gene_symbol"),
  gene_id = c("label_gene_id", "gene_id"),
  protein_id = c("label_protein_id", "protein_id", "protein_Id")
)
.h5ad_role_col <- function(cols, role) {
  hit <- intersect(.h5ad_role_cols[[role]], cols)
  if (length(hit)) hit[[1]] else NULL
}

# Find a per-feature gene-length column for FPKM/TPM, tolerant of source (SUSHI
# featWidth, featureCounts Length, a GTF-derived gene_length, …). Preference-
# ordered; matches with or without a canonical "label_" prefix. NULL if none.
.h5ad_find_length_col <- function(cols) {
  norm <- gsub("^label_", "", tolower(cols))
  for (p in c(
    "^feat[._]?width$",
    "^gene[._]?length$",
    "^exonic[._]?length$",
    "^eff(ective)?[._]?length$",
    "^length$",
    "^width$"
  )) {
    hit <- which(grepl(p, norm))
    if (length(hit)) return(cols[hit[1]])
  }
  NULL
}

# Read a .h5ad into memory (anndataR >= 1.2.3, see DESCRIPTION).
.h5ad_read <- function(path) {
  anndataR::read_h5ad(path, as = "InMemoryAnnData")
}

# Pull the concrete pieces out of an anndataR object as plain R structures, so
# the rest of the class never touches the AnnData API. Layers are materialised
# as samples x features matrices with obs/var dimnames.
.h5ad_extract <- function(adata) {
  obs <- as.data.frame(adata$obs, check.names = FALSE, stringsAsFactors = FALSE)
  var <- as.data.frame(adata$var, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(obs) <- adata$obs_names
  rownames(var) <- adata$var_names

  layers <- adata$layers %||% list()
  layers <- lapply(layers, as.matrix)
  if (!"counts" %in% names(layers) && !is.null(adata$X)) {
    layers[["counts"]] <- as.matrix(adata$X)
  }
  for (k in names(layers)) {
    dimnames(layers[[k]]) <- list(adata$obs_names, adata$var_names)
  }

  varm <- adata$varm %||% list()
  varm <- lapply(varm, function(m) {
    d <- as.data.frame(m, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(d) <- adata$var_names
    d
  })

  list(obs = obs, var = var, layers = layers, varm = varm, uns = adata$uns %||% list())
}

# ---- class ------------------------------------------------------------------

#' Concrete `ExploreDataset` backend over a canonical AnnData (`.h5ad`)
#'
#' Reads a canonical exploreDE/1.0 AnnData and implements the accessor
#' contract inherited from [ExploreDataset]. Accepts a path to a `.h5ad`
#' file or a live `anndataR` object. This is the backend for user-supplied
#' data ("bring your own data") and the eventual single backend once
#' producers emit canonical `.h5ad` directly.
#'
#' @return An R6 class generator; `$new()` returns an `ExploreDatasetH5AD`
#'   object.
#' @examples
#' library(SummarizedExperiment)
#' se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
#' de <- as.data.frame(rowData(se)$constrast_B_vs_A)
#' path <- tempfile(fileext = ".h5ad")
#' build_explore_h5ad(
#'   x = list(raw = assay(se, "rawData"), transformed = assay(se, "transformedData")),
#'   out = path,
#'   omics = "proteomics",
#'   coldata = colData(se),
#'   rowdata = de,
#'   factors = "G_",
#'   de = list(B_vs_A = de)
#' )
#' eds <- ExploreDatasetH5AD$new(path)
#' eds
#' head(eds$de("B_vs_A"))
#' @export
ExploreDatasetH5AD <- R6::R6Class(
  classname = "ExploreDatasetH5AD",
  inherit = ExploreDataset,

  public = list(
    #' @description Read a canonical `.h5ad` or take a live AnnData, build the config from it, and prepare the
    #'   feature and DE tables.
    #' @param input An `anndataR` AnnData object, or a path to an existing `.h5ad` file.
    initialize = function(input) {
      adata <- if (inherits(input, "AbstractAnnData")) {
        input
      } else if (is.character(input) && length(input) == 1L && file.exists(input)) {
        .h5ad_read(input)
      } else {
        stop("ExploreDatasetH5AD: input must be an AnnData or a path to a .h5ad")
      }
      p <- .h5ad_extract(adata)
      private$.obs <- p$obs
      private$.var <- p$var
      private$.layers <- p$layers
      private$.varm <- p$varm
      private$.uns <- p$uns
      private$.config <- private$.build_config()
      private$.prep()
      invisible(self)
    },

    #' @description The AnnData parts read at construction, as `list(obs, var, layers, varm, uns)` of plain R
    #'   objects.
    backend = function() {
      list(obs = private$.obs, var = private$.var, layers = private$.layers, varm = private$.varm, uns = private$.uns)
    },
    #' @description Returns `"h5ad"`.
    backend_kind = function() "h5ad",

    #' @description Sample metadata from `obs`, with the `factor_` prefix stripped from column names and a `Name`
    #'   column added when missing.
    samples = function() {
      cd <- private$.obs
      colnames(cd) <- .h5ad_strip_factor(colnames(cd))
      # Mirror the SE backend: a `Name` column must always exist — ordination
      # and other modules left_join on it. Genomics gets a factor-ordered Name;
      # every other omics (proteomics, crispr) a plain-character Name from the
      # sample ids (rownames).
      if (!"Name" %in% colnames(cd)) {
        cd$Name <- if (identical(private$.config$omics_class, "genomics")) {
          factor(rownames(cd), levels = rownames(cd))
        } else {
          rownames(cd)
        }
      }
      cd
    },

    #' @description Feature annotation from `var`, with `feature_id` and `display_name` columns, ordered by the
    #'   p-value of the primary contrast.
    features = function() private$.features_df,

    #' @description DE table of one contrast from `varm`, with the legacy aliases `log2Ratio`, `pValue` and `fdr`
    #'   added; an empty data.frame when the slot is missing. Errors on an unknown contrast.
    #' @param contrast Contrast name, one of `contrasts()`.
    de = function(contrast) {
      stopifnot(contrast %in% self$contrasts())
      private$.de_dfs[[contrast]] %||% data.frame()
    },

    #' @description Canonical assay names from the config's assay map whose layer exists; computed assays are
    #'   always included.
    assay_names = function() {
      keys <- names(private$.config$assay_map %||% list())
      keys[vapply(
        keys,
        function(k) {
          v <- private$.config$assay_map[[k]]
          is.null(v) || v %in% names(private$.layers) # NULL = computed
        },
        logical(1)
      )]
    },

    #' @description One assay as a features x samples matrix (the layer transposed); `normalised`,
    #'   `normalised_log2`, `vst`, `fpkm` and `tpm` are computed when the config maps them to no layer.
    #' @param name Canonical assay name, one of `assay_names()`.
    #' @param features Feature IDs (rows) to keep; `NULL` keeps all.
    #' @param samples Sample IDs (columns) to keep; `NULL` keeps all.
    assay = function(name, features = NULL, samples = NULL) {
      layer_key <- private$.config$assay_map[[name]]
      m <- if (!is.null(layer_key)) {
        t(private$.layers[[layer_key]]) # samples x features -> features x samples
      } else {
        switch(
          name,
          normalised = private$.compute_normalised(),
          normalised_log2 = log2(
            self$assay("normalised") +
              (private$.config$background_expression %||% 10)
          ),
          vst = private$.compute_vst(),
          fpkm = private$.compute_fpkm(),
          tpm = private$.compute_tpm(),
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

    #' @description Enrichment results of one contrast and kind, rebuilt from `uns$enrichment` and
    #'   `uns$enrichment_support`.
    #' @param contrast Contrast name, one of `contrasts()`.
    #' @param kind `"ora"` or `"gsea"`.
    enrichment = function(contrast, kind) {
      stopifnot(kind %in% c("ora", "gsea"))
      stopifnot(contrast %in% self$contrasts())
      # The enrichment blocks are keyed by the BARE contrast slot; contrast_slot
      # carries the "de_" prefix (that's the varm key), so strip it.
      # Files written before 2026-09-25 also carry serialised R objects in
      # uns$enrichment_native. They are never read: unserialize() on file content
      # can run code.
      slot <- sub("^de_", "", private$.contrast_slot[[contrast]] %||% contrast)
      private$.rebuild_enrichment(slot, kind)
    },

    # Prefer the stored tidy mirror; fall back to the base class, which derives
    # it from the rebuilt enrichment objects.
    #' @description Tidy enrichment table of one contrast and kind: the table stored in `uns$enrichment` when
    #'   present, otherwise derived from `enrichment()` as in the base class.
    #' @param contrast Contrast name, one of `contrasts()`.
    #' @param kind `"ora"` or `"gsea"`.
    enrichment_table = function(contrast, kind) {
      stopifnot(kind %in% c("ora", "gsea"))
      slot <- sub("^de_", "", private$.contrast_slot[[contrast]] %||% contrast)
      stored <- ((private$.uns[["enrichment"]] %||% list())[[slot]] %||% list())[[kind]]
      if (is.data.frame(stored) && nrow(stored)) {
        return(stored)
      }
      super$enrichment_table(contrast, kind)
    },

    # Canonical spec drops enrichInput; reconstruct on demand from varm + var
    # so the (accessor-migrated) pathway modules have what they need.
    #' @description Enrichment inputs rebuilt from the DE table: `features()`, log2 fold changes sorted decreasing,
    #'   and up/down feature selections at the contrast's FDR and log2FC thresholds.
    #' @param contrast Contrast name, one of `contrasts()`.
    enrichment_input = function(contrast) {
      de <- tryCatch(self$de(contrast), error = function(e) NULL)
      if (is.null(de) || !nrow(de)) {
        return(list(seq_anno = NULL, ranked_l2fc = NULL, selections = list()))
      }
      l2fc <- setNames(de$effect_log2fc, de$feature_id)
      cm <- private$.config$contrast_meta[[contrast]]
      fdr_t <- cm$thresholds$fdr %||% 0.05
      l2fc_t <- cm$thresholds$log2fc %||% 0
      sig <- !is.na(de$score_fdr) & de$score_fdr < fdr_t
      list(
        seq_anno = self$features(),
        ranked_l2fc = sort(l2fc[!is.na(l2fc)], decreasing = TRUE),
        selections = list(
          upGenes = de$feature_id[sig & de$effect_log2fc > l2fc_t],
          downGenes = de$feature_id[sig & de$effect_log2fc < -l2fc_t]
        )
      )
    }
  ),

  private = list(
    .obs = NULL,
    .var = NULL,
    .layers = NULL,
    .varm = NULL,
    .uns = NULL,
    .features_df = NULL,
    .de_dfs = list(),
    .contrast_slot = list(),

    # Rebuild clusterProfiler objects from the tidy table + support block
    # (uns$enrichment / uns$enrichment_support), named like the native path
    # ("ora/BP/upGenes", "gsea/BP", or just "ora" for a bare object).
    .rebuild_enrichment = function(slot, kind) {
      sup_all <- private$.uns$enrichment_support[[slot]] %||% list()
      tidy_all <- private$.uns[["enrichment"]][[slot]] %||% list()
      pick <- if (kind == "gsea") "gsea" else setdiff(names(sup_all), "gsea")
      out <- list()
      for (nm in intersect(pick, names(sup_all))) {
        objs <- .enr_rebuild(tidy_all[[nm]], sup_all[[nm]])
        if (!length(objs)) {
          next
        }
        names(objs) <- ifelse(names(objs) == "default", nm, paste0(nm, "/", names(objs)))
        out <- c(out, objs)
      }
      out
    },

    # Matrix-based VST — mirrors eds_compute_vst() (which takes an SE).
    # Size-factor-normalised counts, for files built from raw counts alone
    # (SUSHI results ship their own normalised layer). DESeq2 median-of-ratios;
    # when every feature has a zero somewhere (no geometric mean), fall back to
    # library-size scaling.
    .compute_normalised = function() {
      cts <- self$assay("counts") # features x samples
      sf <- tryCatch(DESeq2::estimateSizeFactorsForMatrix(cts), error = function(e) NULL)
      if (is.null(sf) || any(!is.finite(sf)) || any(sf <= 0)) {
        sf <- colSums(cts) / mean(colSums(cts))
      }
      sweep(cts, 2, sf, "/")
    },

    .compute_vst = function() {
      forVST <- t(as.matrix(private$.layers[["counts"]])) # features x samples
      mode(forVST) <- "integer"
      tryCatch(
        {
          if (nrow(forVST) < 1000) {
            DESeq2::varianceStabilizingTransformation(forVST, blind = FALSE)
          } else {
            DESeq2::vst(forVST, blind = FALSE)
          }
        },
        error = function(e) {
          DESeq2::varianceStabilizingTransformation(forVST, blind = FALSE)
        }
      )
    },

    # Per-feature gene lengths (featWidth), aligned to feature_id.
    .feat_width = function() {
      col <- private$.config$featwidth_col
      setNames(as.numeric(private$.features_df[[col]]), rownames(private$.features_df))
    },

    # FPKM / TPM from raw counts + featWidth — mirrors ezRun::getRpkm/getTpm.
    .compute_fpkm = function() {
      cts <- self$assay("counts") # features x samples
      fw <- private$.feat_width()[rownames(cts)]
      sweep(sweep(cts * 1e9, 1, fw, "/"), 2, colSums(cts), "/")
    },
    .compute_tpm = function() {
      cts <- self$assay("counts")
      fw <- private$.feat_width()[rownames(cts)]
      rate <- sweep(cts, 1, fw, "/")
      sweep(rate * 1e6, 2, colSums(rate), "/")
    },

    .build_config = function() {
      uns <- private$.uns
      .h5ad_require_identity(uns)
      omics_type <- uns[["omics_type"]]
      omics_class <- if (omics_type == "proteomics") "proteomics" else "genomics"
      feature_id_kind <- uns[["feature_id_kind"]]
      data_type <- switch(
        omics_type,
        rnaseq_bulk = "RNASeq",
        crispr_countqc = "CrisprCountQC",
        proteomics = switch(feature_id_kind, peptide = "peptidomics", phospho = "phosphoproteomics", "proteomics")
      )

      factor_cols <- grep("^factor_", colnames(private$.obs), value = TRUE)
      factor_names <- .h5ad_strip_factor(factor_cols)
      factor_levels <- setNames(
        lapply(factor_cols, function(fc) {
          v <- private$.obs[[fc]]
          if (is.factor(v)) levels(v) else sort(unique(as.character(v)))
        }),
        factor_names
      )

      grouping <- .h5ad_strip_factor(
        uns$design$grouping %||% (if (length(factor_cols)) factor_cols[[1]] else NA_character_)
      )

      # Declared significance-semantics mapping (uns$exploreDE$mapping): the
      # display label + direction of an optional extra score, and the fdr
      # direction. Resolved per-contrast as default (+) per-contrast override.
      de_map <- ((uns$exploreDE %||% list())$mapping %||% list())$de %||% list()
      de_default <- de_map$default %||% list()

      # Contrasts: one per uns$de entry, keyed test_vs_reference to match the
      # SE subclass; remember which varm slot each maps to.
      cm <- list()
      slot_map <- list()
      # [[ ]], not $: Without a `de` key, uns$de partial-matches uns$design.
      de_meta <- uns[["de"]] %||% list()
      for (slot in names(de_meta)) {
        d <- de_meta[[slot]]
        test <- d$test_group %||% NA_character_
        ref <- d$reference_group %||% NA_character_
        # Prefer the caller's original contrast name (contrast_formula) so the
        # key matches what they typed; else rebuild an identifier from groups.
        key <- d$contrast_formula %||%
          (if (!is.na(test) && !is.na(ref)) {
            paste0(.h5ad_safe(test), "_vs_", .h5ad_safe(ref))
          } else {
            slot
          })
        sem <- utils::modifyList(de_default, de_map$contrasts[[slot]] %||% list())
        cm[[key]] <- list(
          test_group = test,
          reference_group = ref,
          factor_name = grouping,
          method = d$model %||% "DESeq2",
          thresholds = list(fdr = 0.05, log2fc = 0, go_p = 0.01),
          # Display-only extra score (e.g. SAINTexpress SaintScore): its column
          # label. NULL when the file declares none. (Older files may also carry
          # score_kind / fdr_kind tokens; nothing used them, so they are ignored.)
          score_label = sem$score_label
        )
        slot_map[[key]] <- paste0("de_", slot)
      }
      private$.contrast_slot <- slot_map
      # Union of extra-score display labels across contrasts (for the DE table).
      de_display_extras <- unique(Filter(nzchar, vapply(cm, function(x) x$score_label %||% "", character(1))))

      # Assays are driven by whichever canonical layers are present (raw is not
      # privileged). A NULL map value = computed on the fly in assay(); built
      # via c(list(...)) because `lst$k <- NULL` would DROP the key.
      have <- names(private$.layers)
      # A per-feature gene-length column (featWidth / gene_length / Length / …)
      # enables computed FPKM/TPM — no ezRun, no GTF lookup, just base R over
      # lengths already in the data.
      fw_col <- .h5ad_find_length_col(colnames(private$.var))
      can_len <- !is.null(fw_col) && "counts" %in% have
      if (omics_class == "proteomics") {
        # Proteomics carries intensities, not RNA-seq counts. Expose the
        # intensity layers directly and never synthesise RNA-seq-only derived
        # assays (vst/tpm/fpkm). Mirrors ExploreDatasetSE_Proteomics, which maps
        # a transformed "intensities" (primary) + optional "intensities_raw".
        prim <- intersect(c("intensities", "normalised", "log2", "normalised_log2"), have)
        prim <- if (length(prim)) prim[[1]] else have[[1]]
        raw <- setdiff(intersect(c("intensities_raw", "counts", "raw"), have), prim)
        assay_map <- c(
          setNames(list(prim), "intensities"),
          if (length(raw)) setNames(list(raw[[1]]), "intensities_raw")
        )
        assay_labels <- c(intensities = "Intensities", intensities_raw = "Raw intensities")
        raw_layer_name <- "intensities"
      } else {
        assay_map <- c(
          if ("counts" %in% have) list(counts = "counts"),
          if ("normalised" %in% have) {
            list(normalised = "normalised")
          } else if ("counts" %in% have) {
            list(normalised = NULL)
          }, # computed
          if ("log2" %in% have) {
            list(normalised_log2 = "log2")
          } else if (any(c("normalised", "counts") %in% have)) {
            list(normalised_log2 = NULL)
          }, # computed
          if ("vst" %in% have) {
            list(vst = "vst")
          } else if ("counts" %in% have) {
            list(vst = NULL)
          }, # computed
          if ("tpm" %in% have) {
            list(tpm = "tpm")
          } else if (can_len) {
            list(tpm = NULL)
          },
          if ("fpkm" %in% have) {
            list(fpkm = "fpkm")
          } else if (can_len) {
            list(fpkm = NULL)
          }
        )
        if (!length(assay_map)) {
          assay_map <- setNames(list(have[1]), have[1])
        }
        assay_labels <- c(
          counts = "Raw",
          normalised = "Normalised",
          normalised_log2 = "Normalised + Log2",
          vst = "VST",
          tpm = "TPM",
          fpkm = "FPKM"
        )
        raw_layer_name <- if ("counts" %in% names(assay_map)) "counts" else names(assay_map)[1]
      }

      ref_build <- uns$reference_genome %||% (private$.obs$organism %||% NA_character_)[1]
      if (is.null(ref_build) || is.na(ref_build) || !nzchar(as.character(ref_build))) {
        ref_build <- NA_character_
      }

      # Enrichment kinds come from the support block only (spec §7.4): a method
      # without one has no rebuildable objects, so it does not count. Any
      # contrast counts: a file whose first contrast has no ORA may still carry
      # ORA for a later one.
      enrich_kinds <- character()
      en <- uns$enrichment_support %||% list()
      methods <- unique(unlist(lapply(en, names), use.names = FALSE))
      if (any(grepl("ora", methods))) {
        enrich_kinds <- c(enrich_kinds, "ora")
      }
      if ("gsea" %in% methods) {
        enrich_kinds <- c(enrich_kinds, "gsea")
      }

      list(
        schema_version = uns[["schema_version"]],
        omics_class = omics_class,
        data_type = data_type,
        organism = (private$.obs$organism %||% NA_character_)[1],
        ref_build = ref_build,
        feature_id_kind = feature_id_kind,
        factors = factor_names,
        factor_levels = factor_levels,
        primary_grouping = grouping,
        secondary_grouping = if (length(factor_names) > 1) factor_names[[2]] else NA_character_,
        assay_map = assay_map,
        assay_labels = assay_labels,
        # Proteomics intensities sit on a log2 scale; the low-signal floor the
        # SE proteomics backend uses is 1, not the RNA-seq count floor of 10.
        background_expression = if (omics_class == "proteomics") 1 else 10,
        featwidth_col = fw_col,
        raw_assay_name = raw_layer_name,
        primary_assay_name = if ("normalised_log2" %in% names(assay_map)) {
          "normalised_log2"
        } else {
          raw_layer_name
        },
        display_name_field = .h5ad_role_col(colnames(private$.var), "symbol") %||% NA_character_,
        description_field = "description",
        primary_contrast = if (length(cm)) names(cm)[[1]] else NA_character_,
        contrast_meta = cm,
        de_display_extras = de_display_extras,
        enrichment_kinds = enrich_kinds,
        enrichment_filter = (uns$enrichment_filter %||% list())$padj_max %||% NA_real_,
        enrichment_locator = list(),
        sgrna_to_gene = NULL,
        proteomics = NULL,
        source = list(
          backend_kind = "h5ad",
          schema_version = uns[["schema_version"]],
          built_at = Sys.time()
        )
      )
    },

    .prep = function() {
      var <- private$.var
      cfg <- private$.config

      var$feature_id <- rownames(var)

      dfield <- cfg$display_name_field
      var$display_name <- if (dfield %in% colnames(var)) as.character(var[[dfield]]) else var$feature_id
      blank <- is.na(var$display_name) | !nzchar(var$display_name)
      var$display_name[blank] <- var$feature_id[blank]

      # Canonical GO label columns -> the "GO BP/MF/CC" the app's
      # all_pathways() and menu has_go check expect.
      go_ren <- c(label_go_bp = "GO BP", label_go_mf = "GO MF", label_go_cc = "GO CC")
      for (old in names(go_ren)) {
        if (old %in% colnames(var)) colnames(var)[colnames(var) == old] <- go_ren[[old]]
      }

      private$.features_df <- var

      # proteomics_choices: the gene / protein lookup maps that drive the
      # proteomics-only feature picker (server-featureBucket / server-correlation).
      # Mirrors ExploreDatasetSE_Proteomics$.prep(): gene symbols come from the
      # symbol column; a protein-id set is exposed only when it differs from the
      # feature_id (so protein-level data — where feature_id IS the protein_Id —
      # yields NULL protein choices, matching the SE backend).
      if (identical(cfg$omics_class, "proteomics")) {
        fid <- as.character(var$feature_id)
        gcol <- .h5ad_role_col(colnames(var), "symbol")
        gene_lookup <- gene_choices <- NULL
        if (!is.null(gcol)) {
          gene <- as.character(var[[gcol]])
          valid <- !is.na(gene) & nzchar(gene) & gene != fid
          if (any(valid)) {
            gene_lookup <- split(fid[valid], gene[valid])
            gene_choices <- unique(gene[valid])
          }
        }
        pcol <- .h5ad_role_col(colnames(var), "protein_id")
        protein_lookup <- protein_choices <- NULL
        if (!is.null(pcol)) {
          pid <- as.character(var[[pcol]])
          valid <- !is.na(pid) & nzchar(pid) & pid != fid
          if (any(valid)) {
            protein_lookup <- split(fid[valid], pid[valid])
            protein_choices <- unique(pid[valid])
          }
        }
        private$.config$proteomics <- list(
          gene = list(choices = gene_choices, lookup = gene_lookup),
          protein = list(choices = protein_choices, lookup = protein_lookup)
        )
      } else if (identical(cfg$feature_id_kind, "isoform")) {
        # Isoform RNA-seq: gene symbol -> transcripts, as for the SE backends.
        gcol <- .h5ad_role_col(colnames(var), "symbol")
        if (!is.null(gcol)) {
          private$.config$proteomics <- eds_isoform_gene_choices(var$feature_id, var[[gcol]])
        }
      }

      # de_dfs: canonical varm + the feature columns the app modules expect;
      # eds_finalise_de() adds the legacy aliases and usedInTest, and sorts.
      de_dfs <- list()
      for (key in names(cfg$contrast_meta)) {
        slot <- private$.contrast_slot[[key]]
        d <- private$.varm[[slot]]
        if (is.null(d)) {
          next
        }
        d$feature_id <- rownames(d)
        # Optional extra score -> its declared display label (e.g. SaintScore),
        # so the DE table can surface it (display-only, like nrPeptides).
        if ("score" %in% colnames(d)) {
          sl <- cfg$contrast_meta[[key]]$score_label %||% "score"
          colnames(d)[colnames(d) == "score"] <- sl
        }
        fv <- private$.features_df
        d$gene_name <- fv[d$feature_id, "display_name"]
        # Columns several app modules (DE table, volcano, DE QC) expect from a
        # SUSHI seqAnno; supply equivalents so those tabs work for uploaded data.
        gid_col <- .h5ad_role_col(colnames(fv), "gene_id")
        if (!is.null(gid_col) && !"gene_id" %in% colnames(d)) {
          d$gene_id <- fv[d$feature_id, gid_col]
        }
        if ("description" %in% colnames(fv) && !"description" %in% colnames(d)) {
          d$description <- fv[d$feature_id, "description"]
        } else if (!"description" %in% colnames(d)) {
          # As for a SUSHI rowData without one: the volcano and DE table select it by name.
          d$description <- d$gene_name
        }
        # Feature length and GC content for the DE QC bias plots, which read
        # featWidth / gc from de(). build_explore_h5ad() stores them as
        # label_featwidth / label_gc.
        if (!is.null(cfg$featwidth_col) && !"featWidth" %in% colnames(d)) {
          d$featWidth <- as.numeric(fv[d$feature_id, cfg$featwidth_col])
        }
        gc_col <- intersect(c("label_gc", "gc"), colnames(fv))
        if (length(gc_col) && !"gc" %in% colnames(d)) {
          d$gc <- as.numeric(fv[d$feature_id, gc_col[[1]]])
        }
        de_dfs[[key]] <- eds_finalise_de(d)
      }
      private$.de_dfs <- de_dfs

      # Order features by the primary contrast's p-value so autocompletes and
      # top-N filters surface the most-significant features first — mirrors
      # ExploreDatasetSE_RNASeq, which sorts features() by pValue. The h5ad var
      # carries no inline p-value, so borrow the ordering from the primary
      # contrast's DE table (itself already p-value sorted); features absent
      # from / NA in that contrast keep their original relative order at the end.
      pc <- cfg$primary_contrast
      if (!is.null(pc) && !is.na(pc) && !is.null(de_dfs[[pc]]) && "feature_id" %in% colnames(de_dfs[[pc]])) {
        ord <- as.character(de_dfs[[pc]]$feature_id)
        ord <- ord[ord %in% rownames(private$.features_df)]
        rest <- setdiff(rownames(private$.features_df), ord)
        private$.features_df <- private$.features_df[c(ord, rest), , drop = FALSE]
      }
    }
  )
)
