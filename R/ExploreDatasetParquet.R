# =============================================================================
#  ExploreDatasetParquet  —  ExploreDataset backend over a parquet directory
#
#  The exploreDElight (Tauri) backend. A native Rust process reads the user's
#  canonical .h5ad and writes the slices exploreDE needs as parquet into a temp
#  directory; this class reads that directory. See
#  plans/exploredelight/exploredelight_tauri_plan.md §2 for why the bridge exists at all
#  (hdf5r's read path is broken under wasm, so R cannot open the .h5ad itself).
#
#  It deliberately INHERITS FROM ExploreDatasetH5AD rather than from
#  ExploreDataset. That class's .build_config() and .prep() never touch the
#  AnnData API — they read private$.obs / .var / .layers / .varm / .uns, which
#  are plain data.frames, matrices and lists. Only initialize() is h5ad-shaped.
#  So the whole of the config derivation, the DE-table assembly, the legacy
#  column aliases and the proteomics lookups are inherited unchanged, and any
#  future fix to them lands here for free. Overriding more than the loader
#  would fork logic that has no business being forked.
#
#  Directory layout consumed here — this is the CONTRACT the Rust writer must
#  satisfy, and `plans/exploredelight/spike_h5ad_bridge/h5ad_to_parquet_dir.R` is its
#  executable specification:
#
#    manifest.json               schema, source identity, orientation, factor levels
#    obs.parquet                 `_index` = sample id, + sample metadata columns
#    var.parquet                 `_index` = feature id, + feature annotation columns
#    uns.json                    uns, MINUS the enrichment payloads
#    layers/<name>.parquet       `_index` = feature id, one column per sample
#    varm/<slot>.parquet         `_index` = feature id, + DE stat columns
#    enrichment/<slot>__<kind>.b64   hoisted base64 payloads, read on demand
#
# =============================================================================

# ---- internal helpers -------------------------------------------------------

# Read one parquet file as a plain data.frame.
#
# nanoparquet is the reader that matters: it is the one with a wasm build, so
# it is what actually runs inside webR. arrow is a cluster-side convenience for
# testing this class against a real dataset without a browser — it must never
# be the only thing that works, or the tests would be grading a reader we do
# not ship. Proven interoperable with arrow-rs-written parquet 2026-08-28
# (plans/exploredelight/spike_h5ad_bridge/README.md).
.pq_read <- function(path) {
  if (!file.exists(path)) {
    stop("ExploreDatasetParquet: missing file: ", path)
  }
  if (requireNamespace("nanoparquet", quietly = TRUE)) {
    return(as.data.frame(nanoparquet::read_parquet(path), check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (requireNamespace("arrow", quietly = TRUE)) {
    return(as.data.frame(arrow::read_parquet(path), check.names = FALSE, stringsAsFactors = FALSE))
  }
  stop("ExploreDatasetParquet needs nanoparquet (shipped) or arrow (testing)")
}

# Split the `_index` column off and promote it to rownames, which is how the
# inherited code identifies samples and features. Parquet has no concept of
# rownames, so the writer carries them as an ordinary column.
.pq_read_indexed <- function(path) {
  d <- .pq_read(path)
  if (!"_index" %in% colnames(d)) {
    stop("ExploreDatasetParquet: ", basename(path), " has no `_index` column")
  }
  idx <- as.character(d[["_index"]])
  d <- d[, setdiff(colnames(d), "_index"), drop = FALSE]
  rownames(d) <- idx
  d
}

# Read a layer as a features x samples numeric matrix, then TRANSPOSE it.
#
# The parquet is written features-as-rows because that is the long axis (25k
# features vs ~20 samples), which keeps the column count sane. But
# private$.layers must hold samples x features to match the AnnData
# orientation the inherited assay() expects — it transposes back on the way
# out. Two transposes of a 25k x 20 double is a few milliseconds, and paying
# it buys an exactly-inherited accessor rather than a forked one.
.pq_read_layer <- function(path) {
  d <- .pq_read_indexed(path)
  m <- as.matrix(d)
  storage.mode(m) <- "double"
  t(m)
}

# Restore factor columns from the manifest.
#
# Parquet has no factor type, so a factor column arrives as character and the
# inherited .build_config() would fall back to sort(unique(x)) for its levels
# — silently ALPHABETISING the groups. That reorders every boxplot, every
# legend and every contrast direction relative to the same data read from the
# .h5ad. The writer records the levels; this puts them back.
.pq_restore_factors <- function(df, levels_list) {
  for (col in intersect(names(levels_list), colnames(df))) {
    lv <- as.character(levels_list[[col]])
    if (length(lv)) df[[col]] <- factor(as.character(df[[col]]), levels = lv)
  }
  df
}

# ---- class ------------------------------------------------------------------

#' Concrete `ExploreDataset` backend over a parquet directory
#'
#' The exploreDElight backend. Reads the parquet directory a native Rust
#' process writes from a canonical exploreDE/1.0 `.h5ad`, and satisfies the
#' same accessor contract as every other backend. Inherits all config
#' derivation and DE assembly from [ExploreDatasetH5AD]; only the loader
#' differs.
#'
#' @param input Path to the parquet directory.
#' @param config Optional pre-built config list; derived from `uns.json` when
#'   omitted.
#'
#' @keywords internal
ExploreDatasetParquet <- R6::R6Class(
  classname = "ExploreDatasetParquet",
  inherit = ExploreDatasetH5AD,

  public = list(
    initialize = function(input, config = NULL) {
      if (!is.character(input) || length(input) != 1L || !dir.exists(input)) {
        stop("ExploreDatasetParquet: input must be a path to a parquet directory")
      }

      manifest <- jsonlite::fromJSON(
        file.path(input, "manifest.json"),
        simplifyVector = TRUE,
        simplifyDataFrame = FALSE,
        simplifyMatrix = FALSE
      )

      obs <- .pq_read_indexed(file.path(input, "obs.parquet"))
      var <- .pq_read_indexed(file.path(input, "var.parquet"))
      obs <- .pq_restore_factors(obs, manifest$factor_levels %||% list())

      layer_dir <- file.path(input, "layers")
      layers <- list()
      for (nm in (manifest$layers %||% character())) {
        layers[[nm]] <- .pq_read_layer(file.path(layer_dir, paste0(nm, ".parquet")))
      }

      varm_dir <- file.path(input, "varm")
      varm <- list()
      for (nm in (manifest$varm %||% character())) {
        varm[[nm]] <- .pq_read_indexed(file.path(varm_dir, paste0(nm, ".parquet")))
      }

      # simplifyDataFrame/Matrix off: uns is a nest of scalars, character
      # vectors and named lists, and letting jsonlite coerce a two-field list
      # into a data.frame would change its shape under the inherited config
      # reader. simplifyVector stays ON so scalars arrive as scalars rather
      # than length-1 lists.
      uns_path <- file.path(input, "uns.json")
      uns <- if (file.exists(uns_path)) {
        jsonlite::fromJSON(uns_path, simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
      } else {
        list()
      }

      private$.obs <- obs
      private$.var <- var
      private$.layers <- layers
      private$.varm <- varm
      private$.uns <- uns
      private$.dir <- input
      private$.manifest <- manifest
      private$.config <- config %||% private$.build_config()
      private$.prep()
      private$.warn_if_stale()
      invisible(self)
    },

    backend_kind = function() "parquet",

    # Hydrate-then-delegate. The payloads live in their own files so that
    # loading a dataset does not pay for enrichment nobody asked for (measured:
    # 13 MB of a 29 MB export). The first call for a contrast pulls the blobs
    # in, after which the inherited implementation runs against exactly the
    # structure it expects — no forked decode logic.
    enrichment = function(contrast, kind) {
      private$.hydrate_enrichment(contrast)
      super$enrichment(contrast, kind)
    },

    # The manifest as written, for diagnostics.
    manifest = function() private$.manifest,

    # Where the Rust process put the files. Useful for diagnostics and for the
    # app to clean up its temp dir on exit.
    source_dir = function() private$.dir,

    # Full content check against the source .h5ad, as opposed to the cheap
    # size/mtime comparison done at load. Returns TRUE when the recorded md5
    # still matches, FALSE when it does not, and NA when there is nothing to
    # check against (a live AnnData was exported, or the source has moved).
    verify_source = function(path = NULL) {
      src <- private$.manifest$source %||% list()
      p <- path %||% src$path
      if (is.null(p) || is.null(src$md5) || !nzchar(src$md5) || !file.exists(p)) {
        return(NA)
      }
      identical(unname(tools::md5sum(p)), src$md5)
    },

    # clusterProfiler has no wasm build, so the inherited implementation's
    # go2term() call would error rather than degrade. Return the GO terms when
    # the package is there (this class also runs on the cluster, under test)
    # and an empty vector when it is not, which is what the callers already
    # treat as "this dataset has no GO annotation".
    all_pathways = function() {
      if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        return(character(0))
      }
      super$all_pathways()
    }
  ),

  private = list(
    .dir = NULL,
    .manifest = NULL,
    .hydrated = character(),

    # Pull the hoisted base64 payloads for one contrast into .uns, once.
    .hydrate_enrichment = function(contrast) {
      slot <- sub("^de_", "", private$.contrast_slot[[contrast]] %||% contrast)
      if (slot %in% private$.hydrated) {
        return(invisible(NULL))
      }
      files <- (private$.manifest$enrichment_files %||% list())[[slot]]
      for (kind in names(files %||% list())) {
        f <- file.path(private$.dir, "enrichment", files[[kind]])
        if (!file.exists(f)) {
          next
        }
        private$.uns$enrichment_native[[slot]][[kind]] <-
          paste(readLines(f, warn = FALSE), collapse = "")
      }
      private$.hydrated <- c(private$.hydrated, slot)
      invisible(NULL)
    },

    # A parquet directory is a derived cache. If the file it came from is still
    # where the manifest says and has since changed, every number in this object
    # is out of date — warn rather than let it pass as current. size and mtime
    # only, because hashing the source on every load would cost more than the
    # load itself; verify_source() is there when certainty is wanted.
    .warn_if_stale = function() {
      src <- private$.manifest$source %||% list()
      if (is.null(src$path) || !nzchar(src$path) || !file.exists(src$path)) {
        return(invisible(NULL))
      }
      fi <- file.info(src$path)
      changed <- (!is.null(src$size) && !isTRUE(all.equal(as.numeric(fi$size), as.numeric(src$size)))) ||
        (!is.null(src$mtime) &&
          !identical(format(fi$mtime, "%Y-%m-%dT%H:%M:%S"), src$mtime))
      if (changed) {
        warning(
          "This parquet directory was exported from '",
          src$basename,
          "' on ",
          private$.manifest$written_at,
          ", but that file has changed since. The data you are looking at ",
          "is out of date -- re-run export_h5ad_to_parquet().",
          call. = FALSE
        )
      }
      invisible(NULL)
    },

    # DESeq2 is Bioconductor and has no wasm build, so a computed VST cannot
    # run in Light. This is not a limitation worth hiding behind a confusing
    # "package not found" error: VST is a COMPUTATION, and Light is explicitly
    # pure-visualisation. When the producer already stored a vst layer the
    # inherited assay_map points straight at it and this is never reached.
    .compute_vst = function() {
      if (requireNamespace("DESeq2", quietly = TRUE)) {
        return(super$.compute_vst())
      }
      stop(
        "VST is not available in this build: it must be computed with DESeq2, ",
        "which exploreDElight does not ship. Ask the producer to write a ",
        "`vst` layer into the .h5ad, or use a different assay."
      )
    }
  )
)
