# =============================================================================
#  Raw APB (apb2) results — read through apbCatalogR, converted to exploreDE/1.0
#
#  An apb2 result keeps the quantification software's own column names and
#  records what they mean in uns$apb. apbCatalogR (Suggests) reads the file and
#  answers what each column means. This file turns one level of the result into
#  a canonical exploreDE/1.0 AnnData in memory, mapped as APB Export's
#  `explorede` target maps it (rules v1_0), so a raw result and its export load
#  alike. The rest of the app sees an ordinary ExploreDatasetH5AD.
#
#  Every apbCatalogR call sits in the adapter section below: An API change in
#  apbCatalogR (it is pre-release) touches only that section.
# =============================================================================

# ---- adapter: the only place apbCatalogR is called -------------------------

.apb_available <- function() requireNamespace("apbCatalogR", quietly = TRUE)

.apb_require <- function() {
  if (!.apb_available()) {
    stop(
      "This file is a raw APB result. Reading it needs the apbCatalogR package, which is not ",
      "installed here. Export the result with APB Export's exploreDE output instead.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# The result, one anndataR AnnData per level (`$levels`), from its file.
.apb_read <- function(path) {
  .apb_require()
  apbCatalogR::read_parsed_levels(path)
}

# The `apb` catalogue set (APB's own fields) bound to the result; NULL when
# apbCatalogR can't bind it, so the caller falls back to the file's records.
.apb_catalog <- function(result) {
  tryCatch(apbCatalogR::DataCatalog$new(result), error = function(e) NULL)
}

# One `apb` var field of a level, by key: `value` is the vector, or NULL when
# the software doesn't report it; `error` says why the catalogue couldn't answer
# (an unreviewed rule, a stale field, an ambiguous answer), else NULL.
.apb_var <- function(catalog, level, key) {
  if (is.null(catalog)) {
    return(list(value = NULL, error = "the catalogue could not be bound to this result"))
  }
  tryCatch(
    list(value = catalog$levels[[level]]$var[[key]], error = NULL),
    error = function(e) list(value = NULL, error = conditionMessage(e))
  )
}

.apb_versions <- function() {
  cv <- tryCatch(apbCatalogR::catalog_version(), error = function(e) list())
  list(
    apbCatalogR = as.character(utils::packageVersion("apbCatalogR")),
    catalog_bundle = cv[["version"]] %||% NA_character_
  )
}

# ---- detection ---------------------------------------------------------------

# The names under /uns of an .h5ad, read without loading the file; NULL when it
# can't be read (the normal reader then reports the problem).
.h5ad_uns_names <- function(path) {
  tryCatch(
    {
      ls <- rhdf5::h5ls(path, recursive = 2, datasetinfo = FALSE)
      ls$name[ls$group == "/uns"]
    },
    error = function(e) NULL
  )
}

# Is this a raw APB result? apb2 writes its root record to uns["apb"]. A
# canonical file, APB Export's `explorede` output included, also carries
# uns["schema_version"] and is read as canonical; a prolfquapp DEA AnnData
# carries uns["prolfquapp"].
.h5ad_is_apb_raw <- function(uns_names) {
  "apb" %in% uns_names && !any(c("schema_version", "prolfquapp") %in% uns_names)
}

# ---- conversion ----------------------------------------------------------------

# The levels exploreDE reads, in the order it picks one, with what each becomes.
.APB_LEVELS <- list(
  protein = list(kind = "protein", unit_key = NULL, hierarchy = list(levels = "protein", unit = "protein")),
  peptide = list(
    kind = "peptide",
    unit_key = "peptide",
    hierarchy = list(levels = c("protein", "peptide"), unit = "peptide")
  ),
  peptidoform_ion = list(
    kind = "precursor",
    unit_key = "peptidoform_ion",
    hierarchy = list(levels = c("protein", "peptide", "precursor"), unit = "precursor")
  )
)

# FASTA columns, per level, in the order the label falls back through them (as
# the `explorede` rules read them), plus the description column.
.APB_FASTA <- list(
  protein = list(label = "fasta_gene_names", description = "fasta_descriptions"),
  peptide = list(
    label = c("reported_leading_gene_name", "reported_leading_id"),
    description = "reported_leading_description"
  ),
  peptidoform_ion = list(
    label = c("reported_leading_gene_name", "reported_leading_id"),
    description = "reported_leading_description"
  )
)

.apb_nonempty <- function(x) !is.na(x) & nzchar(x)

# One FASTA column of a level from varm, as apb-export looks it up: The
# `fasta_validation` table, then `fasta`. NULL when the result has no FASTA.
.apb_fasta_col <- function(adata, col) {
  varm <- adata$varm %||% list()
  for (tbl in c("fasta_validation", "fasta")) {
    d <- varm[[tbl]]
    if (!is.null(d) && col %in% colnames(d)) {
      return(as.character(as.data.frame(d)[[col]]))
    }
  }
  NULL
}

# A short display name, as apb-export's `short_name`: Per row, the first input
# with a value, cut to its first `;`-separated member and, for a `db|ACC|NAME`
# member, to NAME; then `_` and `append`, when given.
.apb_short_name <- function(inputs, append = NULL, separator = ";") {
  inputs <- Filter(Negate(is.null), inputs)
  n <- length(inputs[[1]])
  first <- rep(NA_character_, n)
  for (x in inputs) {
    x <- as.character(x)
    take <- is.na(first) & .apb_nonempty(x)
    first[take] <- x[take]
  }
  member <- vapply(strsplit(first, separator, fixed = TRUE), function(p) if (length(p)) p[[1]] else NA_character_, "")
  fields <- strsplit(member, "|", fixed = TRUE)
  short <- ifelse(lengths(fields) == 3L, vapply(fields, function(f) f[3], ""), member)
  if (!is.null(append)) {
    short <- ifelse(is.na(short), NA_character_, paste(short, append, sep = "_"))
  }
  short
}

# obs annotation columns renamed as the `explorede` rules do: SDRF
# `factor_value_<x>` and prolfquapp `group*` / `bait*` columns become
# `factor_<name>` (lower-case, [a-z0-9_]), text with sorted categories; other
# columns pass through under their own names.
.apb_obs_annotations <- function(obs, keys) {
  out <- list()
  for (nm in setdiff(colnames(obs), keys)) {
    v <- obs[[nm]]
    new <- if (grepl("^factor_value_(.+)$", nm)) {
      sub("^factor_value_(.+)$", "factor_\\1", nm)
    } else if (grepl("^((group|bait).*)$", nm, ignore.case = TRUE)) {
      paste0("factor_", nm)
    } else {
      NULL
    }
    if (!is.null(new)) {
      new <- gsub("[^a-z0-9_]", "_", tolower(new))
      text <- as.character(v)
      v <- factor(text, levels = sort(unique(text[!is.na(text)]), method = "radix"))
      nm <- new
    }
    if (nm %in% names(out)) {
      stop("The APB result's sample annotation gives two columns named '", nm, "'.", call. = FALSE)
    }
    out[[nm]] <- v
  }
  out
}

#' Convert one level of a raw APB result to an `exploreDE/1.0` AnnData
#'
#' Maps as APB Export's `explorede` target (rules v1_0): `X`, the rule's primary
#' layer, becomes log2 `intensities` and linear `intensities_raw`, with zeros
#' missing; SDRF factor values and prolfquapp `group` / `bait` columns become
#' `factor_*`; `feature_id` is the protein group, or `<protein>_<peptide>` /
#' `<protein>_<ProForma ion>`. Unlike that target, decoys are dropped by
#' default, and contaminants are flagged in `var$is_contaminant`.
#'
#' Fields come from apbCatalogR's `apb` catalogue. Where it can't answer (an
#' unreviewed rule, for instance), the file's own records are used instead, and
#' `uns$provenance$catalogue_fallbacks` says which.
#'
#' @param result A `ParsedLevels`, as `apbCatalogR::read_parsed_levels()` returns.
#' @param level The level to convert; by default the first of `protein`,
#'   `peptide`, `peptidoform_ion` the result holds.
#' @param drop_decoys Drop decoy features (default `TRUE`).
#' @return An in-memory anndataR `AnnData` in the `exploreDE/1.0` layout.
#' @noRd
anndata_from_apb_raw <- function(result, level = NULL, drop_decoys = TRUE) {
  have <- intersect(names(.APB_LEVELS), names(result$levels))
  if (!length(have)) {
    stop(
      "This APB result holds none of the levels exploreDE reads (protein, peptide, peptidoform_ion). ",
      "It holds: ",
      paste(names(result$levels), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  level <- level %||% have[[1]]
  if (!level %in% have) {
    stop("This APB result has no level '", level, "'. It has: ", paste(have, collapse = ", "), ".", call. = FALSE)
  }
  spec <- .APB_LEVELS[[level]]
  adata <- result$levels[[level]]
  record <- adata$uns[["apb"]] %||% list()
  storage <- record[["storage"]] %||% list()
  roles <- (record[["roles"]] %||% list())[["columns"]] %||% list()

  obs <- as.data.frame(adata$obs, check.names = FALSE, stringsAsFactors = FALSE)
  var <- as.data.frame(adata$var, check.names = FALSE, stringsAsFactors = FALSE)
  x <- as.matrix(adata$X)

  # Fields by meaning; the file's own records where the catalogue can't answer.
  catalog <- .apb_catalog(result)
  fallbacks <- character()
  field <- function(key, fallback_col) {
    got <- .apb_var(catalog, level, key)
    if (!is.null(got$value)) {
      return(as.vector(got$value))
    }
    if (length(fallback_col) && !is.na(fallback_col) && fallback_col %in% colnames(var)) {
      if (!is.null(got$error)) {
        fallbacks <<- c(fallbacks, paste0(key, " from column ", fallback_col, ": ", got$error))
      }
      return(as.vector(var[[fallback_col]]))
    }
    NULL
  }
  var_key <- (storage[["var_key_columns"]] %||% NA_character_)[[1]]
  protein <- field("protein_assignment", c(roles[["protein_assignment"]], var_key)[[1]])
  if (is.null(protein)) {
    stop("This APB result does not say which protein each feature belongs to.", call. = FALSE)
  }
  protein <- as.character(protein)
  unit <- if (!is.null(spec$unit_key)) as.character(field(spec$unit_key, var_key)) else NULL
  if (!is.null(spec$unit_key) && length(unit) != length(protein)) {
    stop("This APB result does not identify its ", spec$kind, "s.", call. = FALSE)
  }
  decoy <- as.logical(field("decoy", "apb_Decoy") %||% rep(FALSE, nrow(var)))
  contaminant <- as.logical(field("contaminant", "apb_Contaminant") %||% rep(FALSE, nrow(var)))
  decoy[is.na(decoy)] <- FALSE
  contaminant[is.na(contaminant)] <- FALSE

  # Feature ids and labels.
  # Non-empty parts joined with `_`, as apb-export's `join_nonempty`.
  feature_id <- if (is.null(unit)) {
    protein
  } else {
    has_p <- .apb_nonempty(protein)
    has_u <- .apb_nonempty(unit)
    ifelse(has_p & has_u, paste(protein, unit, sep = "_"), ifelse(has_u, unit, protein))
  }
  fasta <- .APB_FASTA[[level]]
  label <- .apb_short_name(
    c(lapply(fasta$label, function(col) .apb_fasta_col(adata, col)), list(protein)),
    append = unit
  )
  description <- .apb_fasta_col(adata, fasta$description)

  # Which features to keep: drop decoys and, below protein level, features
  # without a protein (as the `explorede` rules do).
  not_decoy <- if (isTRUE(drop_decoys)) !decoy else rep(TRUE, length(decoy))
  no_protein <- !.apb_nonempty(protein)
  keep <- if (is.null(unit)) not_decoy else not_decoy & !no_protein
  if (!any(keep)) {
    stop("No features are left in this APB result after removing decoys.", call. = FALSE)
  }
  if (anyDuplicated(feature_id[keep])) {
    dup <- unique(feature_id[keep][duplicated(feature_id[keep])])
    stop(
      "Feature ids repeat in this APB result (e.g. ",
      paste(utils::head(dup, 3), collapse = ", "),
      ").",
      call. = FALSE
    )
  }

  new_var <- data.frame(
    label_protein_id = protein,
    label_gene_symbol = label,
    stringsAsFactors = FALSE
  )
  if (!is.null(description)) {
    new_var$description <- description
  }
  new_var$is_contaminant <- contaminant
  if (!isTRUE(drop_decoys)) {
    new_var$is_decoy <- decoy
  }
  new_var <- new_var[keep, , drop = FALSE]
  rownames(new_var) <- feature_id[keep]

  # Samples: The obs key column gives the sample ids; the other columns are
  # annotation. Multiplexed designs key samples by several columns.
  obs_keys <- storage[["obs_key_columns"]] %||% character()
  if (length(obs_keys) > 1L) {
    stop(
      "This APB result keys its samples by several columns (",
      paste(obs_keys, collapse = ", "),
      "), as multiplexed designs do. exploreDE reads one sample per run.",
      call. = FALSE
    )
  }
  sample_id <- if (length(obs_keys) == 1L && obs_keys %in% colnames(obs)) {
    as.character(obs[[obs_keys]])
  } else {
    as.character(adata$obs_names)
  }
  if (anyDuplicated(sample_id)) {
    stop("Sample ids repeat in this APB result.", call. = FALSE)
  }
  new_obs <- data.frame(sample_id = sample_id, sample_label = sample_id, stringsAsFactors = FALSE)
  annotations <- .apb_obs_annotations(obs, obs_keys)
  for (nm in names(annotations)) {
    new_obs[[nm]] <- annotations[[nm]]
  }
  rownames(new_obs) <- sample_id

  # Matrices, samples x features: Zeros are missing; log2 of what's left.
  raw <- x[, keep, drop = FALSE]
  storage.mode(raw) <- "double"
  raw[!is.na(raw) & raw == 0] <- NA_real_
  intensities <- raw
  intensities[!is.na(intensities) & intensities <= 0] <- NA_real_
  intensities <- log2(intensities)
  dimnames(raw) <- dimnames(intensities) <- list(sample_id, rownames(new_var))

  rule <- tryCatch(
    jsonlite::fromJSON(((record[["parse"]] %||% list())[["provenance"]] %||% list())[["rule_json"]] %||% "null"),
    error = function(e) NULL
  )
  uns <- list(
    schema_version = "exploreDE/1.0",
    omics_type = "proteomics",
    feature_id_kind = spec$kind,
    hierarchy = spec$hierarchy,
    de = list(),
    enrichment = list(),
    enrichment_support = list(),
    provenance = c(
      list(
        source_format = "apb2",
        reader = "exploreDE raw APB reader",
        level = level,
        primary_layer = storage[["primary_layer"]] %||% NA_character_,
        software = rule[["software_name"]] %||% NA_character_,
        decoys_dropped = sum(!not_decoy),
        features_without_protein_dropped = if (is.null(unit)) 0L else sum(no_protein & not_decoy),
        catalogue_fallbacks = fallbacks
      ),
      .apb_versions()
    ),
    apb = record
  )

  anndataR::AnnData(
    X = intensities,
    obs = new_obs,
    var = new_var,
    layers = list(intensities = intensities, intensities_raw = raw),
    uns = uns
  )
}
