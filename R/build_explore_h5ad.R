# =============================================================================
#  build_explore_h5ad()  —  source-agnostic "bring your own data" constructor
#
#  Assembles a canonical exploreDE/1.0 .h5ad from a user's possibly-disparate
#  pieces, per the finalised contract (plans/jambon/byo_data_import_plan.md §0c;
#  memory project_exploreDE_byo_constructor_spec):
#
#    Required : >=1 count matrix (any kind; raw NOT privileged) + colData + omics
#    Optional : rowData, DE (>=1 contrast), pathway statistics
#    Pathways : (A) per-feature annotation -> rowData (tile/heatmap);
#               (B) statistics -> contrast-bound objects (ORA/GSEA tabs)
#    Identity : ONE feature_id, identical + same order across counts & rowData;
#               DE a subset within them; colData rows == counts columns.
#    Routes in: a rich object (DESeqDataSet/edgeR/limma/prolfqua) is the
#               PREFERRED starting point but never a dependency; it and loose
#               pieces COMPOSE.
#    Validation: lenient/permissive first (warn-and-degrade), tighten later.
#
#  .assemble_parts() is the pure generic core the Phase B upload wizard reuses.
# =============================================================================

# ---- small helpers ----------------------------------------------------------

# User-facing omics label -> canonical uns$omics_type.
.omics_type_of <- function(omics) {
  switch(
    tolower(omics),
    rnaseq = "rnaseq_bulk",
    rna = "rnaseq_bulk",
    rnaseq_bulk = "rnaseq_bulk",
    proteomics = "proteomics",
    crispr = "crispr_countqc",
    crispr_countqc = "crispr_countqc",
    stop("build_explore_h5ad: unknown omics '", omics, "' (expected rnaseq / proteomics / crispr)")
  )
}

# Canonicalise a user's count-matrix name to a known layer key. Proteomics
# intensity layers keep their own canonical names (intensities/intensities_raw)
# so the reader can expose them as intensities rather than RNA-seq "counts".
# Unknown names pass through unchanged (map is a NAMED VECTOR, so `[[` on an
# absent key errors — look up by membership, never `map[[key]]` blind).
.canon_count_name <- function(nm) {
  key <- gsub("[^a-z0-9]+", "", tolower(nm))
  map <- c(
    counts = "counts",
    count = "counts",
    raw = "counts",
    rawcounts = "counts",
    normalised = "normalised",
    normalized = "normalised",
    norm = "normalised",
    signal = "normalised",
    xnorm = "normalised",
    log2 = "log2",
    logcounts = "log2",
    log2counts = "log2",
    lognorm = "log2",
    normalisedlog2 = "log2",
    normalizedlog2 = "log2",
    vst = "vst",
    tpm = "tpm",
    fpkm = "fpkm",
    rpkm = "fpkm",
    intensities = "intensities",
    intensity = "intensities",
    transformeddata = "intensities",
    transformed = "intensities",
    intensitiesraw = "intensities_raw",
    rawintensities = "intensities_raw",
    rawdata = "intensities_raw"
  )
  if (key %in% names(map)) map[[key]] else nm
}

# The rowData column recognised as dictionary `concept` (gene_name, protein_Id,
# go_bp, ...; see inst/dictionaries/column_patterns.json), or NULL.
.rd_col <- function(rd, concept) {
  if (is.null(rd)) {
    return(NULL)
  }
  eds_match_col(colnames(rd), eds_rowdata_column_patterns[[concept]])
}

# as.data.frame() IGNORES check.names for an S4 DataFrame (colData/rowData), so
# "Condition [Factor]" gets mangled to "Condition..Factor.". Convert, then
# restore the true colnames from the source object. Safe for plain data.frames.
.as_df_keep_names <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  df <- as.data.frame(x, optional = TRUE, stringsAsFactors = FALSE)
  cn <- colnames(x)
  if (!is.null(cn) && ncol(df) == length(cn)) {
    colnames(df) <- cn
  }
  df
}

# Does `x` contain a real clusterProfiler S4 object anywhere (itself, or
# nested inside lists at any depth)? Decides whether a pathway entry also gets
# an uns$enrichment_support block (what the reader needs to rebuild the
# objects), or is data.frames all the way down.
.contains_native_cp <- function(x) {
  if (inherits(x, c("enrichResult", "gseaResult"))) {
    return(TRUE)
  }
  if (is.data.frame(x)) {
    return(FALSE)
  }
  if (is.list(x) && length(x)) {
    return(any(vapply(x, .contains_native_cp, logical(1))))
  }
  FALSE
}

# The overridable DE columns: log2 fold change, p-value, FDR. Anything a caller
# doesn't name here is resolved by the column dictionary. The optional extra
# `score` column (e.g. SaintScore) is dictionary-only, never overridden.
.DE_COLUMN_TARGETS <- c("effect_log2fc", "score_pvalue", "score_fdr")
.DE_MAPPING_KEYS <- .DE_COLUMN_TARGETS
.ROWDATA_MAPPING_KEYS <- c("id_col", "symbol_col")

# Warn (once per unknown key) about override keys we don't act on, so a typo
# or a retired key (e.g. `score`, `description_col`, `score_kind`, `fdr_kind`)
# is never silently ignored.
.check_mapping_keys <- function(mapping, allowed, arg) {
  bad <- setdiff(names(mapping), allowed)
  if (length(bad)) {
    warning(
      sprintf(
        "build_explore_h5ad: ignoring unknown %s key(s): %s (allowed: %s)",
        arg,
        paste(bad, collapse = ", "),
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Copy explicit-mapping columns to canonical names, then fill the rest from the
# column dictionary. attr(, "de_sources") records which source column
# fed each canonical name, explicit or detected.
.canon_de_columns <- function(df, mapping = NULL) {
  df <- as.data.frame(df, optional = TRUE, stringsAsFactors = FALSE)
  sources <- character()
  if (!is.null(mapping)) {
    for (canon in intersect(names(mapping), .DE_COLUMN_TARGETS)) {
      src <- mapping[[canon]]
      if (!is.null(src) && src %in% colnames(df)) {
        df[[canon]] <- df[[src]]
        sources[[canon]] <- src
      }
    }
  }
  attr(df, "de_sources") <- sources
  eds_standardise_de_columns(df) # only adds a canon col if absent -> mapping wins
}

# Turn a DE column map into the declared uns$exploreDE$mapping$de entry: WHICH
# source column fed each overridden target (effect_log2fc_col/score_pvalue_col/
# score_fdr_col). This is the dictionary itself, kept in the file so a reader (human or code) never has to
# re-guess "log2fc vs log2 ratio vs log2FoldChange" for THIS dataset — see
# plans/jambon/anndata_manifest_design.md. NULL entries drop.
.de_semantics <- function(mapping) {
  if (is.null(mapping)) {
    return(list())
  }
  keep <- list()
  for (canon in .DE_COLUMN_TARGETS) {
    src <- mapping[[canon]]
    if (!is.null(src)) keep[[paste0(canon, "_col")]] <- src
  }
  keep
}

# Explicit rowdata_columns override wins (if the named column actually exists
# in rd); else the column dictionary's `concept`. Same resolution order as DE
# columns: explicit mapping > dictionary.
.rowdata_col <- function(rd, explicit, concept) {
  if (!is.null(explicit) && !is.null(rd) && explicit %in% colnames(rd)) {
    return(explicit)
  }
  .rd_col(rd, concept)
}

# ---- feature_level: detection-with-override --------------------------------
# One knob spanning omics: rnaseq {gene,isoform}, proteomics
# {protein,peptide,phospho}, crispr {sgrna}. Detection mirrors
# ExploreDatasetSE_Proteomics (site -> phospho, peptide_Id -> peptide) and the
# RNA-seq isoform rule (one transcript_id per row). Explicit `override` always wins.

.FEATURE_LEVELS <- list(
  rnaseq_bulk = c("gene", "isoform"),
  proteomics = c("protein", "peptide", "phospho"),
  crispr_countqc = "sgrna"
)

.resolve_feature_level <- function(omics_type, rowdata, override = NULL) {
  choices <- .FEATURE_LEVELS[[omics_type]] %||% "gene"
  if (!is.null(override)) {
    if (!override %in% choices) {
      stop(sprintf(
        "feature_level '%s' invalid for omics '%s' (expected: %s)",
        override,
        omics_type,
        paste(choices, collapse = ", ")
      ))
    }
    return(override)
  }
  rd <- as.data.frame(rowdata %||% data.frame())
  if (omics_type == "proteomics") {
    if (!is.null(.rd_col(rd, "site"))) {
      return("phospho")
    }
    if (!is.null(.rd_col(rd, "peptide_Id"))) {
      return("peptide")
    }
    return("protein")
  }
  # Isoform only when each row IS one transcript. Gene-level tables often carry
  # a transcript_id column listing the gene's transcripts ("ENST...; ENST..."),
  # which must not flip the level.
  tcol <- if (omics_type == "rnaseq_bulk") .rd_col(rd, "transcript_id")
  if (!is.null(tcol) && !any(grepl("[;,]", rd[[tcol]]))) {
    return("isoform")
  }
  choices[[1]]
}

# Build the canonical feature_id for a resolved level. base_ids = the join key
# (count-matrix rownames); gid/sym = id/symbol; rd = rowData (may be NULL);
# id_col = the resolved (explicit-override-or-detected) protein/feature id
# column for the protein/peptide/phospho levels — an explicit rowdata_columns
# override wins here too, same resolution order as everywhere else.
.build_fid_for_level <- function(level, base_ids, gid, sym, rd, id_col = NULL) {
  scu <- requireNamespace("scuttle", quietly = TRUE)
  # A source id column can carry NA / empty entries (e.g. SAINT's protein_Id has
  # gaps where the clean rowname does not) — coalesce those to the join key so we
  # never build an NA feature_id (which would blow up data.frame(row.names = .)).
  coalesce_ids <- function(ids) {
    ids <- as.character(ids)
    bad <- is.na(ids) | !nzchar(ids)
    if (any(bad)) {
      ids[bad] <- base_ids[bad]
    }
    ids
  }
  prot_plus <- function(sub_concept) {
    # protein_Id [+ site/peptide]
    pcol <- id_col %||% .rd_col(rd, "protein_Id")
    if (is.null(pcol)) {
      return(make.unique(base_ids))
    }
    pid <- coalesce_ids(rd[[pcol]])
    scol <- .rd_col(rd, sub_concept)
    if (is.null(scol)) {
      return(make.unique(pid))
    }
    make.unique(paste(pid, gsub("~", "_", as.character(rd[[scol]])), sep = "_"))
  }
  fid <- switch(
    level,
    gene = if (scu) scuttle::uniquifyFeatureNames(gid, sym) else make.unique(as.character(sym)),
    isoform = {
      tcol <- .rd_col(rd, "transcript_id")
      tid <- if (!is.null(tcol)) as.character(rd[[tcol]]) else base_ids
      if (scu) scuttle::uniquifyFeatureNames(tid, sym) else make.unique(as.character(tid))
    },
    protein = {
      pcol <- id_col %||% .rd_col(rd, "protein_Id")
      make.unique(coalesce_ids(if (!is.null(pcol)) rd[[pcol]] else base_ids))
    },
    peptide = prot_plus("peptide_Id"),
    phospho = prot_plus("site"),
    sgrna = make.unique(base_ids),
    make.unique(as.character(gid))
  )
  # Final safety net: any residual NA / empty id -> its join key, re-uniquified.
  fid <- as.character(fid)
  bad <- is.na(fid) | !nzchar(fid)
  if (any(bad)) {
    fid[bad] <- base_ids[bad]
  }
  fid <- make.unique(fid)
  fid[grep("^_", fid)] <- gsub("_", "", fid[grep("^_", fid)])
  as.character(fid)
}

# Canonical uns$hierarchy for a level.
.hierarchy_for_level <- function(level) {
  switch(
    level,
    isoform = list(levels = c("gene", "transcript"), unit = "transcript"),
    phospho = list(levels = c("protein", "site"), unit = "site"),
    peptide = list(levels = c("protein", "peptide"), unit = "peptide"),
    protein = list(levels = "protein", unit = "protein"),
    sgrna = list(levels = "sgrna", unit = "sgrna"),
    list(levels = "gene", unit = "gene")
  )
}

# ---- generic assembler (pure) ----------------------------------------------

# counts   : named list of feature x sample matrices (>=1). A single matrix is
#            wrapped as list(counts = .).
# coldata  : data.frame of sample metadata (rows aligned to count columns).
# omics    : "rnaseq" | "proteomics" | "crispr".
# rowdata  : optional per-feature annotation (rows aligned to count features);
#            pathway-annotation columns (GO/KEGG) land in var.
# de       : optional named list keyed by contrast of DE tables.
# de_columns : optional explicit column mapping (global, or per-contrast list).
# pathways : optional named list keyed by contrast; each a named list of
#            clusterProfiler objects and/or data.frames.
# factors  : optional colData columns that are experimental factors (else all
#            character/factor-typed columns).
# enrichment_padj_max is one number, or a list keyed by pathway method name
# (ora / gsea); a method missing from the list falls back to the largest value.
.padj_max_for <- function(padj_max, method) {
  if (!is.list(padj_max)) {
    return(padj_max)
  }
  padj_max[[method]] %||% max(unlist(padj_max))
}

.assemble_parts <- function(
  counts,
  coldata,
  omics = "rnaseq",
  feature_level = NULL,
  rowdata = NULL,
  de = list(),
  de_columns = NULL,
  pathways = list(),
  enrichment_padj_max = 0.1,
  factors = NULL,
  gene_symbol = NULL,
  feature_id = NULL,
  rowdata_columns = NULL,
  organism = NA_character_,
  ref_build = NA_character_
) {
  omics_type <- .omics_type_of(omics)

  # -- counts: normalise to a named list keyed by canonical layer names -------
  if (is.matrix(counts) || inherits(counts, "Matrix")) {
    counts <- list(counts = counts)
  }
  if (!length(counts)) {
    stop("build_explore_h5ad: at least one count matrix is required")
  }
  names(counts) <- vapply(names(counts) %||% "counts", .canon_count_name, character(1))
  counts <- lapply(counts, as.matrix)

  primary <- if ("counts" %in% names(counts)) counts[["counts"]] else counts[[1]]
  base_ids <- rownames(primary) %||% as.character(seq_len(nrow(primary)))
  sample_ids <- colnames(primary) %||% rownames(as.data.frame(coldata))
  if (anyDuplicated(base_ids)) {
    warning("duplicate feature ids in counts \u2014 made unique")
    base_ids <- make.unique(base_ids)
  }

  # -- feature_id (the join key): symbol/id from args or rowData --------------
  rd <- .as_df_keep_names(rowdata)
  if (!is.null(rd)) {
    if (nrow(rd) != length(base_ids)) {
      warning("rowData rows != feature count \u2014 aligning by rowname where possible")
    }
    rn <- rownames(rd)
    if (!is.null(rn) && all(base_ids %in% rn)) rd <- rd[base_ids, , drop = FALSE]
  }
  .check_mapping_keys(rowdata_columns, .ROWDATA_MAPPING_KEYS, "rowdata_columns")
  sym_col <- .rowdata_col(rd, rowdata_columns$symbol_col, "gene_name")
  id_col <- .rowdata_col(rd, rowdata_columns$id_col, "gene_id")
  # An EXPLICIT id_col override (not merely auto-detected via the gene-id-like
  # patterns above) also drives protein/peptide/phospho feature-id construction
  # below — e.g. a foreign proteomics producer's "UniProt_Accession" column
  # standing in for the built-in "protein_Id" detection.
  explicit_id_col <- if (!is.null(rd) && !is.null(rowdata_columns$id_col) && rowdata_columns$id_col %in% colnames(rd)) {
    rowdata_columns$id_col
  }

  gid <- feature_id %||% (if (!is.null(id_col)) as.character(rd[[id_col]]) else base_ids)
  sym <- gene_symbol %||% (if (!is.null(sym_col)) as.character(rd[[sym_col]]) else gid)

  # feature_level: detection-with-override, then level-aware feature_id.
  level <- .resolve_feature_level(omics_type, rd, feature_level)
  fid <- .build_fid_for_level(level, base_ids, gid, sym, rd, id_col = explicit_id_col)
  id2fid <- setNames(fid, base_ids) # user identifiers -> canonical fid

  # -- obs (samples) : factors get "factor_" prefix; others carried plain -----
  obs <- .as_df_keep_names(coldata)
  if (!is.null(rownames(obs)) && all(sample_ids %in% rownames(obs))) {
    obs <- obs[sample_ids, , drop = FALSE]
  }
  if (is.null(factors)) {
    is_cat <- vapply(obs, function(x) is.character(x) || is.factor(x), logical(1))
    # Keep only INFORMATIVE categoricals as experimental factors. Exclude:
    #  * identifier-like columns (one distinct value per sample: sample names,
    #    raw-file paths, channels) — labels, not factors; and
    #  * constant columns (a single distinct value, e.g. proteomics
    #    isotope.label = "light") — no contrast, and would otherwise surface as
    #    a spurious grouping/batch factor.
    # Mirrors the proteomics Name/raw.file exclusion without hardcoding names.
    # Fall back to all categoricals only if this leaves none (degenerate case).
    n <- nrow(obs)
    nuniq <- vapply(obs, function(x) length(unique(as.character(x))), integer(1))
    informative <- nuniq > 1L & nuniq < n
    factors <- colnames(obs)[is_cat & informative]
    if (!length(factors)) factors <- colnames(obs)[is_cat]
  }
  new_names <- colnames(obs)
  for (i in seq_along(new_names)) {
    if (colnames(obs)[i] %in% factors) {
      # Strip a trailing SUSHI " [Factor]" tag before sanitising, so
      # "Condition [Factor]" -> factor_condition (not factor_condition_factor_).
      clean <- sub("[ _]*\\[factor\\]$", "", colnames(obs)[i], ignore.case = TRUE)
      new_names[i] <- paste0("factor_", gsub("[^A-Za-z0-9]+", "_", tolower(clean)))
      # Level ORDER is not cosmetic: it drives boxplot x-axis order, legend
      # order and group colours. A bare factor(as.character(x)) sorts levels
      # alphabetically, which turns a "0h, 6h, 24h" time course into
      # "0h, 24h, 6h" without erroring — the plot just quietly misrepresents
      # the design.
      v <- as.character(obs[[i]])
      obs[[i]] <- if (is.factor(obs[[i]])) {
        # An explicit level order is the user's stated design; keep it, but
        # drop levels no sample uses so empty groups don't appear.
        factor(v, levels = intersect(levels(obs[[i]]), unique(v)))
      } else {
        # Derived from plain character: mixedsort, not sort, so "0h/6h/24h"
        # and "S1/S2/S10" come out the way a human means them. Matches
        # ExploreDatasetSE.R:80 and helperFunctions.R:484, so the h5ad and SE
        # backends agree on level order rather than diverging.
        factor(v, levels = unique(gtools::mixedsort(v)))
      }
    }
  }
  colnames(obs) <- new_names
  rownames(obs) <- sample_ids
  obs$sample_label <- sample_ids
  if (!is.na(organism)) {
    obs$organism <- organism
  }
  factor_prefixed <- colnames(obs)[grepl("^factor_", colnames(obs))]

  # -- var (features) : ids/symbols + pathway ANNOTATION columns from rowData --
  var <- data.frame(row.names = fid, label_gene_id = gid, label_gene_symbol = sym, stringsAsFactors = FALSE)
  desc_col <- NULL
  if (!is.null(rd)) {
    desc_col <- .rd_col(rd, "description")
    if (!is.null(desc_col)) {
      var$description <- sub(" \\[.*", "", as.character(rd[[desc_col]]))
    }
    # Pathway annotation -> canonical label_go_* (reader renames to "GO BP" etc.)
    for (which in c("bp", "mf", "cc")) {
      col <- .rd_col(rd, paste0("go_", which))
      if (!is.null(col)) var[[paste0("label_go_", which)]] <- rd[[col]]
    }
    # Carry any remaining rowData columns through as label_* (annotation).
    used <- unique(c(sym_col, id_col, desc_col))
    for (col in setdiff(colnames(rd), used)) {
      dest <- if (grepl("^GO ", col) || grepl("^label_", col)) {
        col
      } else {
        paste0("label_", gsub("[^A-Za-z0-9]+", "_", tolower(col)))
      }
      if (!dest %in% colnames(var)) var[[dest]] <- rd[[col]]
    }
  }

  # -- layers : samples x features, canonical names, logical -> integer -------
  flip <- function(m) {
    m <- as.matrix(m)
    ids <- rownames(m) %||% base_ids
    keep <- match(base_ids, ids) # align every matrix to primary
    m <- m[keep, , drop = FALSE]
    x <- t(m)
    dimnames(x) <- list(sample_ids, fid)
    if (is.logical(x)) {
      storage.mode(x) <- "integer"
    }
    x
  }
  layers <- lapply(counts, flip)
  # derive log2 if we have normalised but not log2
  if ("normalised" %in% names(layers) && !"log2" %in% names(layers)) {
    layers$log2 <- log2(layers$normalised + 10)
  }
  X <- layers[["counts"]] %||% layers[[1]]

  # -- DE : per contrast, canonicalise columns, index to full fid space -------
  # Beyond the three canonical stat columns, a contrast may carry an OPTIONAL
  # extra `score` (a second thresholdable/displayable significance axis, e.g.
  # SAINTexpress's SaintScore) — see the declared-mapping block below. It is
  # display-only this milestone (like nrPeptides), so it rides in varm but is
  # never wired into the FDR threshold path.
  per_contrast_map <- !is.null(de_columns) &&
    length(intersect(names(de_columns), names(de))) > 0
  if (per_contrast_map) {
    for (cn in intersect(names(de_columns), names(de))) {
      .check_mapping_keys(de_columns[[cn]], .DE_MAPPING_KEYS, sprintf("de_columns[['%s']]", cn))
    }
  } else {
    .check_mapping_keys(de_columns, .DE_MAPPING_KEYS, "de_columns")
  }
  varm <- list()
  uns_de <- list()
  de_slot <- list()
  score_src <- list()
  for (cn in names(de)) {
    mapping <- if (per_contrast_map) de_columns[[cn]] else de_columns
    orig_cols <- colnames(de[[cn]])
    d <- .canon_de_columns(de[[cn]], mapping)
    # The optional extra score (e.g. SaintScore) is found by the dictionary
    # only; remember which column it came from, as that is its display label.
    srcs <- attr(d, "de_sources")
    if ("score" %in% names(srcs)) {
      score_src[[cn]] <- srcs[["score"]]
    }
    # effect_log2fc (magnitude) and score_fdr (the primary significance axis,
    # §6 of the design doc) are the two columns the app actually depends on.
    # score_pvalue is legitimately optional (SAINTexpress has none) and
    # degrades gracefully via has_pvalue() -- never warn on that alone. If
    # NEITHER the mapping NOR the pattern dictionary resolved one of the two
    # required columns, this contrast would otherwise go all-NA with no
    # signal anywhere that anything was wrong.
    missing_stats <- setdiff(c("effect_log2fc", "score_fdr"), colnames(d))
    if (length(missing_stats)) {
      warning(sprintf(
        "DE '%s': could not resolve %s from columns (%s) \u2014 values will be NA. Pass de_columns= for this file, or add a pattern to inst/dictionaries/column_patterns.json if this naming will recur.",
        cn,
        paste(missing_stats, collapse = " and "),
        paste(orig_cols, collapse = ", ")
      ))
    }
    slot <- gsub("[^a-z0-9]+", "_", tolower(cn))
    de_slot[[cn]] <- slot
    can <- data.frame(
      effect_log2fc = rep(NA_real_, length(fid)),
      score_pvalue = rep(NA_real_, length(fid)),
      score_fdr = rep(NA_real_, length(fid)),
      row.names = fid
    )
    tgt <- id2fid[rownames(d)] # DE feature ids -> canonical fid
    ok <- !is.na(tgt)
    if (any(!ok)) {
      warning(sprintf("DE '%s': %d feature(s) not in counts/rowData \u2014 dropped", cn, sum(!ok)))
    }
    pos <- match(tgt[ok], fid)
    stat_cols <- c("effect_log2fc", "score_pvalue", "score_fdr", if ("score" %in% colnames(d)) "score")
    for (col in stat_cols) {
      if (!col %in% colnames(can)) {
        can[[col]] <- rep(NA_real_, length(fid))
      }
      if (col %in% colnames(d)) can[[col]][pos] <- suppressWarnings(as.numeric(d[[col]]))[ok]
    }
    # Whether each feature entered the test (e.g. SUSHI's usedInTest), when the
    # DE table says so; found by the dictionary only. Features the table
    # doesn't list were not tested. Without it, the reader counts a feature as
    # tested when its p-value (or FDR) isn't NA.
    if ("used_in_test" %in% colnames(d)) {
      can$used_in_test <- rep(FALSE, length(fid))
      can$used_in_test[pos] <- eds_as_flag(d$used_in_test)[ok]
    }
    varm[[paste0("de_", slot)]] <- can
    grp <- strsplit(cn, "_vs_|--over--")[[1]]
    uns_de[[slot]] <- list(
      assay_used = "counts",
      factor_used = factor_prefixed[1] %||% NA_character_,
      test_group = grp[1] %||% NA_character_,
      reference_group = grp[2] %||% NA_character_,
      contrast_formula = cn,
      model = "DESeq2"
    )
  }

  # -- declared DE mapping (uns$exploreDE$mapping$de) --------------------------
  # The dictionary itself: which source column fed each overridden target
  # (effect_log2fc/score_pvalue/score_fdr), plus what the canonical stat arrays
  # alone can't express (the extra score's source column + display label).
  # A shared de_columns map => one
  # `de$default`; per-contrast maps => `de$contrasts`. A detected extra score
  # is always recorded per contrast: without its label the reader would not
  # surface it at all (no DE-table column, no volcano score axis).
  map_de <- list()
  if (!is.null(de_columns)) {
    if (per_contrast_map) {
      ctr <- list()
      for (cn in names(de)) {
        s <- .de_semantics(de_columns[[cn]])
        if (length(s)) ctr[[de_slot[[cn]]]] <- s
      }
      if (length(ctr)) map_de$contrasts <- ctr
    } else {
      d0 <- .de_semantics(de_columns)
      if (length(d0)) map_de$default <- d0
    }
  }
  for (cn in names(de)) {
    src <- score_src[[cn]]
    if (is.null(src)) {
      next
    }
    slot <- de_slot[[cn]]
    map_de$contrasts[[slot]] <- c(map_de$contrasts[[slot]], list(score_col = src, score_label = src))
  }

  # -- declared feature-identity mapping (uns$exploreDE$mapping$feature) ------
  # Only recorded when rowdata_columns explicitly named a column AND it was
  # actually used (i.e. it existed in rd) — auto-detected columns stay
  # implicit, same "written only when declared" lean-file principle as `de`.
  map_feature <- list()
  if (!is.null(rowdata_columns)) {
    if (!is.null(rowdata_columns$id_col) && identical(id_col, rowdata_columns$id_col)) {
      map_feature$id_col <- id_col
    }
    if (!is.null(rowdata_columns$symbol_col) && identical(sym_col, rowdata_columns$symbol_col)) {
      map_feature$symbol_col <- sym_col
    }
  }

  # -- Pathway STATISTICS (B) : contrast-bound; keys must be a subset of DE.
  # Each named entry under a contrast is stored AS-IS, no leaf-walking:
  #   * a data.frame, OR a nested list containing only data.frames (e.g. a
  #     hand-built BP/MF/CC x up/down/both tree with no clusterProfiler object
  #     anywhere) -> flattened straight to the tidy uns$enrichment table.
  #   * anything else (a single clusterProfiler object, OR a nested list that
  #     contains at least one — e.g. the SUSHI BP/MF/CC x up/down/both shape)
  #     -> the tidy table (filtered at enrichment_padj_max) PLUS a small
  #     uns$enrichment_support block: object type, and for GSEA the gene sets
  #     of the stored terms only, the ranked gene list and `exponent`. That is
  #     everything the reader needs to rebuild the objects for the network /
  #     ridge / running-score plots. The objects themselves are not stored:
  #     >90% of them is the full gene-set collection, which nothing displays
  #     (plans/exploredelight/enrichment_detach_plan.md).
  enr_support <- list()
  enr_tidy <- list()
  bad_keys <- setdiff(names(pathways), names(de))
  if (length(bad_keys)) {
    warning(sprintf(
      "pathway stats for non-existent contrast(s): %s \u2014 kept but unbound",
      paste(bad_keys, collapse = ", ")
    ))
  }
  for (cn in names(pathways)) {
    slot <- de_slot[[cn]] %||% gsub("[^a-z0-9]+", "_", tolower(cn))
    results <- pathways[[cn]]
    # A single object / data.frame passed directly (not wrapped in a named list).
    if (!is.list(results) || is.data.frame(results) || inherits(results, c("enrichResult", "gseaResult"))) {
      results <- list(result = results)
    }
    for (rn in names(results)) {
      obj <- results[[rn]]
      if (is.data.frame(obj)) {
        enr_tidy[[slot]] <- c(enr_tidy[[slot]], setNames(list(obj), rn))
      } else if (!.contains_native_cp(obj)) {
        tidy <- .enr_tidy_from_native(obj, method = rn)
        if (!is.null(tidy)) {
          enr_tidy[[slot]] <- c(enr_tidy[[slot]], setNames(list(tidy), rn))
        }
      } else {
        tidy <- .enr_tidy_from_native(obj, method = rn, padj_max = .padj_max_for(enrichment_padj_max, rn))
        if (!is.null(tidy)) {
          enr_tidy[[slot]] <- c(enr_tidy[[slot]], setNames(list(tidy), rn))
          sup <- .enr_support_from_native(obj, tidy)
          if (!is.null(sup)) {
            enr_support[[slot]] <- c(enr_support[[slot]], setNames(list(sup), rn))
          }
        }
      }
    }
  }

  uns_list <- list(
    schema_version = "exploreDE/1.0",
    omics_type = omics_type,
    feature_id_kind = level,
    hierarchy = .hierarchy_for_level(level),
    design = list(factors = factor_prefixed, grouping = factor_prefixed[1] %||% NA_character_),
    de = uns_de,
    enrichment = enr_tidy,
    enrichment_support = enr_support,
    # The tidy mirror is FILTERED. Record the threshold, or a
    # consumer cannot tell "cut by this filter" from "not
    # significant" from "enrichment was never run" — three very
    # different things that all look like an absent row.
    enrichment_filter = list(
      # one number, or a list keyed by method (ora / gsea)
      padj_max = if (length(enr_tidy)) enrichment_padj_max else NULL,
      column = "p.adjust"
    ),
    reference_genome = as.character(ref_build),
    provenance = list(source_format = "build_explore_h5ad")
  )
  # `mapping` stays a DISTINCT key from `design` (design = what the experiment
  # is; mapping = where to find things / how to read the stats) -- this is the
  # declared dictionary itself (plans/jambon/anndata_manifest_design.md): which raw
  # column/matrix supplies each concept exploreDE cares about, kept in the file
  # rather than discarded once used at build time. Written only when a
  # declared mapping actually exists, so canonical files stay lean.
  mapping <- list()
  if (length(map_de)) {
    mapping$de <- map_de
  }
  if (length(map_feature)) {
    mapping$feature <- map_feature
  }
  if (length(mapping)) {
    uns_list$exploreDE <- list(mapping = mapping)
  }

  list(X = X, layers = layers, obs = obs, var = var, varm = varm, uns = uns_list)
}

# ---- rich-object extraction (PREFERRED starting point, never required) ------

.pieces_from_deseq2 <- function(dds) {
  counts <- list(counts = SummarizedExperiment::assay(dds, "counts"))
  nd <- tryCatch(DESeq2::counts(dds, normalized = TRUE), error = function(e) NULL)
  if (!is.null(nd)) {
    counts$normalised <- nd
  }
  rd <- .as_df_keep_names(SummarizedExperiment::rowData(dds))
  if (!ncol(rd)) {
    rd <- NULL
  }
  fcols <- tryCatch(all.vars(DESeq2::design(dds)), error = function(e) NULL)
  list(
    counts = counts,
    coldata = .as_df_keep_names(SummarizedExperiment::colData(dds)),
    rowdata = rd,
    factors = intersect(fcols, colnames(SummarizedExperiment::colData(dds)))
  )
}

#' Build a canonical exploreDE `.h5ad` from user-supplied data
#'
#' The source-agnostic "bring your own data" constructor. `x` is either a rich
#' upstream object (a `DESeqDataSet` — the preferred starting point) or the
#' counts themselves (a matrix, or a named list of matrices such as
#' `list(counts=, normalised=, vst=, tpm=)`). Any parts the object doesn't
#' carry — DE tables, pathway statistics, extra count matrices — are supplied
#' as arguments and merged. Only counts + `coldata` + `omics` are required.
#'
#' @param x A `DESeqDataSet`, a counts matrix, or a named list of count
#'   matrices.
#' @param out Output `.h5ad` path.
#' @param omics Required: `"rnaseq"`, `"proteomics"`, or `"crispr"`.
#' @param feature_level Optional override of the feature granularity within the
#'   omics — rnaseq `"gene"`/`"isoform"`, proteomics
#'   `"protein"`/`"peptide"`/`"phospho"`. Omit to auto-detect from `rowdata`
#'   (`site` -> phospho, `peptide_Id` -> peptide, `transcript_id` -> isoform).
#' @param coldata Sample metadata `data.frame`. Required unless `x` is a rich
#'   object that carries it.
#' @param rowdata Optional per-feature annotation (feature names/symbols/ids and
#'   per-feature pathway annotation columns, e.g. `GO BP`).
#' @param de Optional named list of DE tables keyed by contrast
#'   (`DESeqResults` / `data.frame`). Columns are recognised by the column
#'   dictionary (`inst/dictionaries/column_patterns.json`), overridable with
#'   `de_columns`. An optional extra score column (e.g. SAINTexpress's
#'   `SaintScore`) is picked up by the dictionary alone and shown for display
#'   under its own column name.
#' @param de_columns Optional DE column overrides:
#'   `list(effect_log2fc = "<col>", score_pvalue = "<col>", score_fdr = "<col>")`,
#'   applied to all contrasts, or a named list of such lists keyed by contrast.
#'   Named columns are used as given; any not named are still found by the
#'   dictionary. A producer with no raw p-value (e.g. SAINTexpress) simply has
#'   none; the app degrades to FDR-only. Unknown keys are ignored with a
#'   warning. Overrides are persisted at `uns$exploreDE$mapping$de`, so the
#'   file stays self-describing.
#' @param rowdata_columns Optional feature column overrides:
#'   `list(id_col = "<col>", symbol_col = "<col>")` — which `rowdata` column
#'   supplies the feature id and the display symbol. Anything not named is
#'   still found by the dictionary. Unknown keys are ignored with a warning.
#'   Overrides are persisted at `uns$exploreDE$mapping$feature`.
#' @param pathways Optional pathway statistics: a named list keyed by contrast
#'   (keys a subset of `de`); each entry a named list of clusterProfiler objects
#'   and/or `data.frame`s.
#' @param enrichment_padj_max p.adjust cut-off applied to clusterProfiler
#'   results before they are stored (`p.adjust <= enrichment_padj_max`). Either
#'   one number for every method, or a named list keyed by the method names used
#'   in `pathways` (e.g. `list(ora = 0.05, gsea = 0.05)`) so that ORA and GSEA can
#'   follow their own FDR thresholds; a method not named in the list uses the
#'   largest value given. The value(s) are recorded in `uns$enrichment_filter`
#'   and the app uses a per-method value as its significance threshold.
#' @param factors Optional colData columns that are experimental factors
#'   (default: all character/factor-typed columns).
#' @param gene_symbol,organism,ref_build Optional labels.
#' @return The output path, invisibly.
#' @export
build_explore_h5ad <- function(
  x,
  out,
  omics,
  feature_level = NULL,
  coldata = NULL,
  rowdata = NULL,
  de = list(),
  de_columns = NULL,
  rowdata_columns = NULL,
  pathways = list(),
  enrichment_padj_max = 0.1,
  factors = NULL,
  gene_symbol = NULL,
  organism = NA_character_,
  ref_build = NA_character_
) {
  if (missing(omics)) {
    stop("build_explore_h5ad: `omics` is required (rnaseq / proteomics / crispr)")
  }

  # Rich object -> extract the parts it carries; explicit args fill/override.
  if (inherits(x, "DESeqDataSet")) {
    pc <- .pieces_from_deseq2(x)
    counts <- pc$counts
    coldata <- coldata %||% pc$coldata
    rowdata <- rowdata %||% pc$rowdata
    factors <- factors %||% pc$factors
  } else if (inherits(x, "SummarizedExperiment")) {
    counts <- list(counts = SummarizedExperiment::assay(x, 1))
    coldata <- coldata %||% .as_df_keep_names(SummarizedExperiment::colData(x))
    rd <- .as_df_keep_names(SummarizedExperiment::rowData(x))
    rowdata <- rowdata %||% (if (ncol(rd)) rd else NULL)
  } else if (is.matrix(x) || inherits(x, "Matrix") || is.list(x)) {
    counts <- x
  } else {
    stop("build_explore_h5ad: `x` must be a DESeqDataSet, a counts matrix, or a named list of matrices")
  }
  if (is.null(coldata)) {
    stop("build_explore_h5ad: `coldata` is required (no rich object supplied it)")
  }

  de <- lapply(de, .as_df_keep_names)

  parts <- .assemble_parts(
    counts,
    coldata,
    omics = omics,
    feature_level = feature_level,
    rowdata = rowdata,
    de = de,
    de_columns = de_columns,
    rowdata_columns = rowdata_columns,
    pathways = pathways,
    enrichment_padj_max = enrichment_padj_max,
    factors = factors,
    gene_symbol = gene_symbol,
    organism = organism,
    ref_build = ref_build
  )
  adata <- anndataR::AnnData(
    X = parts$X,
    obs = parts$obs,
    var = parts$var,
    layers = parts$layers,
    varm = parts$varm,
    uns = parts$uns
  )
  anndataR::write_h5ad(adata, out, mode = "w") # "w" allows overwrite (default is "w-")
  invisible(out)
}
