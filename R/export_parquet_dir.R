# =============================================================================
#  export_h5ad_to_parquet()  —  reference implementation of the parquet
#  directory contract that ExploreDatasetParquet reads.
#
#  In the exploreDElight Tauri app this job belongs to a native Rust process
#  (plans/exploredelight_tauri_plan.md §2): R inside webR cannot open a .h5ad,
#  because hdf5r's wasm read path is broken. Keeping an R implementation here
#  earns its place three times over:
#
#    * it PINS THE CONTRACT — Rust has a reference to diff against rather than
#      prose, and any drift shows up as a failing equivalence test;
#    * it makes ExploreDatasetParquet testable without cargo or a browser;
#    * it is a genuine fallback — a Light bundle can be produced from the
#      cluster with no Rust involved at all.
# =============================================================================

# Write one data.frame to parquet, preferring the reader we actually ship.
.pq_write <- function(df, path) {
  if (requireNamespace("nanoparquet", quietly = TRUE))
    return(invisible(nanoparquet::write_parquet(df, path)))
  if (requireNamespace("arrow", quietly = TRUE))
    return(invisible(arrow::write_parquet(df, path, compression = "snappy")))
  stop("export_h5ad_to_parquet() needs nanoparquet or arrow")
}

# Rownames travel as an ordinary `_index` column; parquet has no other way to
# carry them. Non-atomic columns are DROPPED rather than coerced — a
# list-column silently flattened into a string is worse than an absent one,
# and nothing in the accessor contract reads one.
.pq_prep_df <- function(df, idx, what = "table", quiet = TRUE) {
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- vapply(df, is.atomic, logical(1))
  if (any(!keep) && !quiet)
    message(sprintf("  %s: dropping %d non-atomic column(s): %s", what,
                    sum(!keep), paste(colnames(df)[!keep], collapse = ", ")))
  df <- df[, keep, drop = FALSE]
  cbind(`_index` = as.character(idx), df, stringsAsFactors = FALSE)
}

#' Export a canonical `.h5ad` to the parquet directory exploreDElight reads
#'
#' Writes the directory layout consumed by [ExploreDatasetParquet]:
#' `manifest.json`, `obs.parquet`, `var.parquet`, `uns.json`,
#' `layers/<name>.parquet` and `varm/<slot>.parquet`.
#'
#' @param input Path to a canonical exploreDE/1.0 `.h5ad`, or a live
#'   `anndataR` object.
#' @param out Output directory; created if absent.
#' @param quiet Suppress per-layer progress messages.
#'
#' @return `out`, invisibly.
#' @export
export_h5ad_to_parquet <- function(input, out, quiet = TRUE) {
  adata <- if (inherits(input, "AbstractAnnData")) {
    input
  } else if (is.character(input) && length(input) == 1L && file.exists(input)) {
    .h5ad_read(input)
  } else {
    stop("export_h5ad_to_parquet(): input must be an AnnData or a path to a .h5ad")
  }

  dir.create(file.path(out, "layers"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "varm"),   recursive = TRUE, showWarnings = FALSE)

  obs_names <- adata$obs_names
  var_names <- adata$var_names
  obs <- as.data.frame(adata$obs, check.names = FALSE, stringsAsFactors = FALSE)
  var <- as.data.frame(adata$var, check.names = FALSE, stringsAsFactors = FALSE)

  # Record factor LEVELS before parquet flattens them to character. Without
  # this the reader falls back to sort(unique(x)) for its levels and silently
  # ALPHABETISES every experimental group — reordering boxplots, legends and
  # contrast directions against the same data read from the .h5ad.
  factor_levels <- list()
  for (cn in colnames(obs))
    if (is.factor(obs[[cn]])) factor_levels[[cn]] <- levels(obs[[cn]])

  .pq_write(.pq_prep_df(obs, obs_names, "obs", quiet), file.path(out, "obs.parquet"))
  .pq_write(.pq_prep_df(var, var_names, "var", quiet), file.path(out, "var.parquet"))

  # AnnData stores samples x features. Parquet gets features-as-ROWS because
  # features is the long axis (25k vs ~25), keeping the column count sane;
  # ExploreDatasetParquet transposes once on load.
  layers <- adata$layers %||% list()
  if (!"counts" %in% names(layers) && !is.null(adata$X)) layers[["counts"]] <- adata$X
  layer_names <- character()
  for (nm in names(layers)) {
    m <- as.matrix(layers[[nm]])
    dimnames(m) <- list(obs_names, var_names)
    fs <- t(m)
    d  <- cbind(`_index` = rownames(fs),
                as.data.frame(fs, check.names = FALSE, stringsAsFactors = FALSE))
    .pq_write(d, file.path(out, "layers", paste0(nm, ".parquet")))
    layer_names <- c(layer_names, nm)
    if (!quiet)
      message(sprintf("  layer %-16s %d features x %d samples", nm, nrow(fs), ncol(fs)))
  }

  # varm rows align to var POSITIONALLY — varm/<slot>/_index in the .h5ad holds
  # row NUMBERS, not feature names. The join is by position and the real
  # feature ids are written into `_index` here, so the reader never has to know.
  varm <- adata$varm %||% list()
  varm_names <- character()
  for (nm in names(varm)) {
    d <- as.data.frame(varm[[nm]], check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(d) != length(var_names))
      stop(sprintf("varm/%s has %d rows but var has %d - positional join unsafe",
                   nm, nrow(d), length(var_names)))
    d <- d[, setdiff(colnames(d), "_index"), drop = FALSE]
    d <- d[, vapply(d, is.atomic, logical(1)), drop = FALSE]
    .pq_write(cbind(`_index` = as.character(var_names), d),
              file.path(out, "varm", paste0(nm, ".parquet")))
    varm_names <- c(varm_names, nm)
    if (!quiet)
      message(sprintf("  varm  %-16s %d rows x %d cols", nm, nrow(d), ncol(d)))
  }

  # uns goes to JSON — except the enrichment payloads, which are hoisted out to
  # their own files and read lazily.
  #
  # Measured on my_study.h5ad: uns.json was 13 MB of a 29 MB export, i.e. 45% of
  # the whole transport, and jsonlite parsed all of it before the first plot
  # could render — for data that is only needed if someone opens a pathway tab.
  # The blobs are base64(zlib(R serialize)), so base64 adds 33% on top and JSON
  # cannot compress it back.
  #
  # The KEY STRUCTURE stays in uns.json with empty-string values: the inherited
  # .build_config() derives enrichment_kinds from names(uns$enrichment_native),
  # so dropping the keys would make has_enrichment() report FALSE. Empty strings
  # also mean an un-hydrated read degrades to an empty list rather than an
  # error, because .h5ad_decode_native() already returns NULL on !nzchar().
  uns <- adata$uns %||% list()
  native <- uns$enrichment_native %||% list()
  enrichment_files <- list()
  if (length(native)) {
    dir.create(file.path(out, "enrichment"), recursive = TRUE, showWarnings = FALSE)
    for (slot in names(native)) {
      for (kind in names(native[[slot]])) {
        b64 <- native[[slot]][[kind]]
        if (!is.character(b64) || !length(b64) || !nzchar(b64[[1]])) next
        fn <- paste0(slot, "__", kind, ".b64")
        writeLines(b64[[1]], file.path(out, "enrichment", fn))
        enrichment_files[[slot]][[kind]] <- fn
        uns$enrichment_native[[slot]][[kind]] <- ""      # keep the key, drop the payload
      }
    }
  }
  write(jsonlite::toJSON(uns, auto_unbox = TRUE, null = "null", digits = NA),
        file.path(out, "uns.json"))

  # Source identity. A parquet directory is a DERIVED CACHE, so without this a
  # stale temp dir left over from a previous file reads as current and the user
  # sees the wrong dataset with no indication anything is wrong. size and mtime
  # are free to compare at load; the md5 is there for verify_source() when you
  # want certainty rather than a cheap check.
  src <- if (is.character(input) && length(input) == 1L && file.exists(input)) {
    fi <- file.info(input)
    list(path     = normalizePath(input, mustWork = FALSE),
         basename = basename(input),
         size     = as.numeric(fi$size),
         mtime    = format(fi$mtime, "%Y-%m-%dT%H:%M:%S"),
         md5      = unname(tools::md5sum(input)))
  } else {
    # A live AnnData has no file behind it, so there is nothing to be stale
    # against. Say so rather than inventing an identity.
    list(path = NULL, basename = NULL, size = NULL, mtime = NULL, md5 = NULL)
  }

  manifest <- list(
    schema_version    = uns$schema_version %||% "exploreDE/1.0",
    writer            = "export_h5ad_to_parquet() [R reference implementation]",
    written_at        = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    source            = src,
    n_obs             = length(obs_names),
    n_var             = length(var_names),
    layer_orientation = "features_x_samples",
    layers            = layer_names,
    varm              = varm_names,
    enrichment_files  = enrichment_files,
    factor_levels     = factor_levels
  )
  write(jsonlite::toJSON(manifest, auto_unbox = TRUE, null = "null"),
        file.path(out, "manifest.json"))

  invisible(out)
}
