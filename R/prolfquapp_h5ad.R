# =============================================================================
# prolfquapp DEA AnnData (.h5ad) -> proteomics SummarizedExperiment
#
# The prolfquapp DEA (B-Fabric A414) writes AnnData.h5ad next to
# SummarizedExperiment.rds. Its layout (namespace uns$prolfquapp) is not the
# canonical exploreDE/1.0 AnnData, so instead of extending ExploreDatasetH5AD
# this rebuilds the SummarizedExperiment prolfquapp wrote alongside it;
# .prolfquapp_legacy_layout() then turns that (and any prolfquapp >= 2.10 rds)
# into the SE shape ExploreDatasetSE_Proteomics already reads: one rowData
# `constrast_<name>` frame per contrast, holding the feature annotation
# followed by the DE statistics.
#
# Mirrors prolfquapp's anndata_to_summarized_experiment(): layers become assays,
# obs becomes colData, each varm matrix plus uns$prolfquapp$varm_columns /
# varm_annotations / varm_column_order becomes a nested rowData frame, and the
# rest of uns$prolfquapp becomes metadata.
# =============================================================================

# uns$prolfquapp keys that describe the AnnData reshape, not the SE metadata.
.prolfquapp_reshape_keys <- c(
  "layer_names",
  "uns_table_columns",
  "varm_columns",
  "varm_annotations",
  "varm_column_order",
  "varm_key_order",
  "uns_list_order"
)

#' Is this an AnnData written by the prolfquapp DEA?
#' @param adata an anndataR AnnData
#' @keywords internal
is_prolfquapp_anndata <- function(adata) {
  !is.null(adata$uns[["prolfquapp"]])
}

.prolfquapp_flat <- function(x) as.vector(unlist(x, use.names = FALSE))

# Slot keys in the order prolfquapp recorded them (HDF5 groups are unordered).
.prolfquapp_ordered_keys <- function(recorded, present, label) {
  recorded <- .prolfquapp_flat(recorded)
  if (!setequal(recorded, present)) {
    stop(
      "prolfquapp AnnData ",
      label,
      " does not list the keys the file carries: ",
      paste(union(setdiff(recorded, present), setdiff(present, recorded)), collapse = ", ")
    )
  }
  recorded
}

# HDF5 hands vectors back as one-dimensional arrays; drop that shape.
.prolfquapp_uns_restore <- function(value) {
  if (is.data.frame(value)) {
    value[] <- lapply(value, .prolfquapp_uns_restore)
    return(value)
  }
  if (is.list(value)) {
    return(lapply(value, .prolfquapp_uns_restore))
  }
  if (length(dim(value)) == 1L) {
    return(as.vector(value))
  }
  value
}

# HDF5 hands list names back in its own order; put them back in the order
# prolfquapp recorded (uns_list_order; absent in files written before it was).
.prolfquapp_restore_list_order <- function(value, order) {
  if (is.null(order) || !is.list(value) || is.data.frame(value)) {
    return(value)
  }
  recorded <- .prolfquapp_flat(order$names)
  value <- value[c(intersect(recorded, names(value)), setdiff(names(value), recorded))]
  for (name in intersect(names(order$children), names(value))) {
    value[[name]] <- .prolfquapp_restore_list_order(value[[name]], order$children[[name]])
  }
  value
}

# uns$prolfquapp -> SE metadata: drop the reshape keys, rebuild the tables
# prolfquapp stores column by column (listed in uns_table_columns) and restore
# the recorded list order.
.prolfquapp_uns_metadata <- function(uns) {
  values <- .prolfquapp_uns_restore(uns[setdiff(names(uns), .prolfquapp_reshape_keys)])
  tables <- uns$uns_table_columns
  for (name in intersect(names(tables), names(values))) {
    columns <- .prolfquapp_ordered_keys(tables[[name]], names(values[[name]]), paste0("uns_table_columns$", name))
    values[[name]] <- as.data.frame(values[[name]][columns], stringsAsFactors = FALSE)
  }
  .prolfquapp_restore_list_order(values, uns$uns_list_order)
}

# One nested rowData frame from a varm matrix: numeric columns from the matrix,
# the non-numeric ones from varm_annotations, in the recorded column order.
.prolfquapp_varm_frame <- function(values, columns, annotations, column_order, var_names) {
  values <- as.matrix(values)
  columns <- .prolfquapp_flat(columns)
  if (ncol(values) != length(columns)) {
    stop("prolfquapp AnnData varm matrix has ", ncol(values), " columns but varm_columns names ", length(columns), ".")
  }
  frame <- as.data.frame(values)
  names(frame) <- columns
  for (name in names(annotations)) {
    frame[[name]] <- .prolfquapp_flat(annotations[[name]])
  }
  frame <- frame[, .prolfquapp_flat(column_order), drop = FALSE]
  rownames(frame) <- var_names
  frame
}

# prolfquapp >= 2.10 SE -> the SE shape ExploreDatasetSE_Proteomics reads.
#
# prolfquapp 2.10 stores the feature annotation once, in rowData(se)$annotation,
# and the contrast frames hold the DE statistics (2.10.0-2.10.4 without the
# feature keys). Every rowData frame is indexed by the feature ids, so each
# contrast frame gets the annotation joined on those ids, in front of its
# statistics; key columns a contrast frame carries are checked against the
# annotation. The declared factors are moved to the front of colData, modelled
# factor first, because the proteomics backend groups by the first
# non-sample column. An SE without rowData$annotation (prolfquapp < 2.10)
# is returned unchanged.
.prolfquapp_legacy_layout <- function(se) {
  rd <- SummarizedExperiment::rowData(se)
  if (!"annotation" %in% names(rd)) {
    return(se)
  }
  md <- S4Vectors::metadata(se)

  cd <- SummarizedExperiment::colData(se)
  factors <- intersect(names(md$analysis_configuration$factors), colnames(cd))
  SummarizedExperiment::colData(se) <- cd[, c(factors, setdiff(colnames(cd), factors)), drop = FALSE]

  annotation <- data.frame(rd$annotation, check.names = FALSE, stringsAsFactors = FALSE)
  keys <- intersect(.prolfquapp_flat(md$feature_keys), colnames(annotation))
  for (name in grep("^constrast_", names(rd), value = TRUE)) {
    frame <- data.frame(rd[[name]], check.names = FALSE, stringsAsFactors = FALSE)
    idx <- match(rownames(frame), rownames(annotation))
    if (anyNA(idx)) {
      stop("prolfquapp SE: '", name, "' has feature ids missing from rowData$annotation.")
    }
    ann <- annotation[idx, , drop = FALSE]
    for (k in intersect(keys, colnames(frame))) {
      known <- !is.na(frame[[k]])
      if (!identical(as.character(frame[[k]][known]), as.character(ann[[k]][known]))) {
        stop("prolfquapp SE: the feature keys of '", name, "' do not match rowData$annotation.")
      }
    }
    res <- cbind(ann, frame[, setdiff(colnames(frame), colnames(ann)), drop = FALSE])
    rownames(res) <- rownames(frame)
    SummarizedExperiment::rowData(se)[[name]] <- res
  }
  se
}

#' Rebuild the prolfquapp DEA SummarizedExperiment from its AnnData
#'
#' @param adata an anndataR AnnData written by the prolfquapp DEA
#'   (`uns$prolfquapp` present)
#' @return the `SummarizedExperiment` prolfquapp wrote alongside the AnnData:
#'   `rowData(se)$annotation` from var, one nested rowData frame per varm entry
#' @export
se_from_prolfquapp_anndata <- function(adata) {
  uns <- adata$uns[["prolfquapp"]]
  if (is.null(uns)) {
    stop("AnnData has no uns$prolfquapp; it was not written by prolfquapp.")
  }

  obs <- as.data.frame(adata$obs, check.names = FALSE, stringsAsFactors = FALSE)
  var <- as.data.frame(adata$var, check.names = FALSE, stringsAsFactors = FALSE)
  obs_names <- adata$obs_names
  var_names <- adata$var_names
  rownames(obs) <- obs_names
  rownames(var) <- var_names

  layer_keys <- .prolfquapp_ordered_keys(uns$layer_names, names(adata$layers), "layer_names")
  assays <- lapply(layer_keys, function(key) {
    values <- t(as.matrix(adata$layers[[key]]))
    dimnames(values) <- list(var_names, obs_names)
    values
  })
  names(assays) <- layer_keys

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(obs, check.names = FALSE),
    metadata = .prolfquapp_uns_metadata(uns)
  )
  SummarizedExperiment::rowData(se)[["annotation"]] <- var

  varm_keys <- .prolfquapp_ordered_keys(uns$varm_key_order, names(adata$varm), "varm_key_order")
  for (key in varm_keys) {
    SummarizedExperiment::rowData(se)[[utils::URLdecode(key)]] <-
      .prolfquapp_varm_frame(
        adata$varm[[key]],
        uns$varm_columns[[key]],
        uns$varm_annotations[[key]],
        uns$varm_column_order[[key]],
        var_names
      )
  }
  se
}
