# =============================================================================
# Internal helpers for the ExploreDataset R6 class hierarchy.
#
# All identifiers in this file are package-internal (not exported) and are
# prefixed with `eds_` to avoid colliding with similar helpers in
# `helperFunctions.R` (e.g. `build_feature_id`, `standardise_columns`)
# whose contracts differ.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# The AUTOMATIC half of exploreDE's dictionary (plans/jambon/anndata_manifest_design.md):
# when no explicit de_columns/rowdata_columns mapping is given, THIS is what
# figures out which column is which for an arbitrary incoming data.frame. Lives
# as data (inst/dictionaries/column_patterns.json, validated by
# inst/schema/exploreDE_column_patterns.schema.json), not hardcoded in R, so
# covering a new pipeline's spelling of "log2 ratio" is a JSON edit, not a
# code change. Loaded once per session; each concept's pattern array is
# collapsed to one alternation regex for eds_standardise_*_columns() below.
eds_load_column_patterns <- function() {
  path <- system.file("dictionaries", "column_patterns.json", package = "exploreDE")
  if (!nzchar(path)) {
    # exploreDE is not always installed as a package -- it's routinely
    # source()'d directly from an arbitrary working directory: app.R's own
    # "dev mode -- no install required" startup, or an ad hoc test script
    # that does `source(file.path(PKG, "R", "utils_explore_dataset.R"))`
    # with PKG an absolute path and getwd() somewhere else entirely. Anchor
    # to THIS FILE'S OWN LOCATION rather than guessing from getwd(): source()
    # records the path it was given as `ofile` in its own evaluation frame,
    # which is still on the call stack right now (this call happens at the
    # top level of this very file, i.e. during that source() call).
    ofile <- NULL
    for (fr in sys.frames()) {
      cand <- fr$ofile
      if (!is.null(cand) && grepl("utils_explore_dataset\\.R$", cand)) { ofile <- cand; break }
    }
    if (!is.null(ofile)) {
      pkg_root <- dirname(dirname(normalizePath(ofile, mustWork = FALSE)))
      cand_path <- file.path(pkg_root, "inst", "dictionaries", "column_patterns.json")
      if (file.exists(cand_path)) path <- cand_path
    }
  }
  if (!nzchar(path)) {
    # Last resort: a few plausible getwd()-relative locations (devtools::
    # load_all()/testthat with cwd == package root or tests/testthat/, or
    # app.R's own inst/app runtime cwd) for the rare case source() left no
    # `ofile` on the stack (e.g. code pasted directly into the console).
    candidates <- c(
      file.path("inst", "dictionaries", "column_patterns.json"),
      file.path("..", "dictionaries", "column_patterns.json"),
      file.path("..", "..", "inst", "dictionaries", "column_patterns.json")
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) path <- hit[[1]]
  }
  if (!nzchar(path) || !file.exists(path))
    stop("exploreDE: inst/dictionaries/column_patterns.json not found ",
        "(checked the installed package location, this file's own directory, ",
        "and dev-mode relative paths). Only the column-standardising helpers ",
        "need it; if you are running a cut-down build that does not import ",
        "arbitrary column names, nothing should be reaching this.")
  raw <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  # Keep each concept's patterns as an ORDERED vector, not one alternation:
  # eds_match_col() tries them in priority order, so the anchored patterns
  # listed first always beat the broad substring catch-alls listed last.
  ordered <- function(group) lapply(group, function(pats) as.character(unlist(pats)))
  list(de = ordered(raw$de %||% list()), rowdata = ordered(raw$rowdata %||% list()))
}

# THE column matcher. Every piece of exploreDE that has to work out "which of
# these columns is the FDR / the gene symbol / the protein id" goes through
# here, against the one JSON dictionary. Patterns are tried in the order the
# dictionary lists them; the first pattern that hits anything decides. So
# `^padj$` wins over the bare substring `fdr`, and a decoy such as
# `log2_mean_intensity` can only be picked when no anchored log2fc pattern
# matched at all.
#   unique = FALSE : take the first matching column (constructor / aliases).
#   unique = TRUE  : if the deciding pattern hits >1 column, it is ambiguous;
#                    return NULL rather than guess (standardise_columns()).
eds_match_col <- function(cols, patterns, unique = FALSE) {
  for (p in patterns) {
    hit <- cols[grepl(p, cols, ignore.case = TRUE)]
    if (length(hit) > 1L && unique) return(NULL)
    if (length(hit)) return(hit[[1]])
  }
  NULL
}

# LAZY on purpose (delayedAssign, not plain assignment). These used to be
# evaluated at source() time, which made merely LOADING this file depend on
# finding the JSON on disk -- and eds_flatten_cp_results() and friends live in
# the same file. exploreDElight sources these four R files into webR, where
# there is no installed package and no inst/ tree, so the whole file failed to
# load and took the enrichment helpers down with it, for a dictionary Light
# never consults.
#
# Deliberately still an ERROR rather than a silent empty default: the loader
# keeps its stop(), it just fires when something actually asks for a pattern
# instead of when the file is read. Loud where it matters, absent where it does
# not.
#
# Note for the installed package: R CMD INSTALL forces these while building the
# lazy-load DB, so an installed exploreDE behaves exactly as before -- the
# laziness only shows up in source() mode, which is the mode that needed it.
delayedAssign(".eds_column_patterns", eds_load_column_patterns())

# Column-name regexes for DE tables, as ordered vectors per concept. Map any
# upstream tool's column to a canonical AnnData-spec name (`effect_log2fc`,
# `score_pvalue`, `score_fdr`, `score`).
delayedAssign("eds_de_column_patterns", .eds_column_patterns$de)

# Column-name regexes for per-feature annotation (ids, symbols, description,
# GO, proteomics protein/peptide/site ...). Read by build_explore_h5ad() and by
# standardise_columns() (the legacy SE proteomics path).
delayedAssign("eds_rowdata_column_patterns", .eds_column_patterns$rowdata)

# Add canonical DE columns to `df` (copy, not rename), so the upstream tool's
# native column names are preserved alongside. Only adds if the canonical
# column doesn't already exist.
#
# Matching runs against the ORIGINAL column names only, so a canonical column
# added earlier in the loop (score_pvalue, score_fdr) can never be picked up
# by a later concept's substring pattern (e.g. "score").
#
# The source column chosen for each canonical name is recorded in
# attr(df, "de_sources"), so the constructor can persist which column fed
# what (e.g. SaintScore -> score, used as its display label).
eds_standardise_de_columns <- function(df) {
  current <- colnames(df)
  sources <- attr(df, "de_sources") %||% character()
  for (canon in names(eds_de_column_patterns)) {
    if (canon %in% colnames(df)) next
    src <- eds_match_col(current, eds_de_column_patterns[[canon]])
    if (!is.null(src)) {
      df[[canon]] <- df[[src]]
      sources[[canon]] <- src
    }
  }
  attr(df, "de_sources") <- sources
  df
}

# SUSHI param$species is rarely set; fall back to the first slash-segment of
# refBuild ("Homo_sapiens/GENCODE/..." → "Homo sapiens").
eds_infer_organism <- function(param) {
  sp <- param$species
  if (!is.null(sp) && !is.na(sp) && nzchar(sp)) return(sp)
  rb <- param$refBuild
  if (is.null(rb) || is.na(rb) || !nzchar(rb)) return(NA_character_)
  parts <- strsplit(rb, "/", fixed = TRUE)[[1]]
  if (!length(parts)) return(NA_character_)
  gsub("_", " ", parts[1], fixed = TRUE)
}

# Recursive flatten of nested clusterProfiler results into a single named list.
# clusterProfiler outputs are typically nested as
# `list(BP = enrichResult, MF = enrichResult, ...)` or deeper for ORA with
# multiple directions (`list(BP = list(upGenes = ..., downGenes = ...), ...)`).
#
# A plain data.frame is ALSO a valid leaf (a BYO producer's own tidy ORA/GSEA
# table nested the same BP/MF/CC x up/down shape, with no clusterProfiler
# object anywhere) -- checked before the generic is.list() recursion because a
# data.frame is itself a list-of-columns in R, so without this branch a
# data.frame leaf gets silently walked into its columns and vanishes instead
# of flattening to itself. See plans/exploredelight/enrichment_detach_plan.md.
eds_flatten_cp_results <- function(x, prefix = "") {
  if (inherits(x, c("enrichResult", "gseaResult")) || is.data.frame(x)) {
    nm <- if (nzchar(prefix)) prefix else "default"
    return(setNames(list(x), nm))
  }
  if (is.list(x) && length(x)) {
    out <- list()
    nms <- names(x) %||% as.character(seq_along(x))
    for (i in seq_along(x)) {
      sub_prefix <- if (nzchar(prefix)) paste0(prefix, "/", nms[i]) else nms[i]
      out <- c(out, eds_flatten_cp_results(x[[i]], sub_prefix))
    }
    return(out)
  }
  list()
}

# Regex-safe identifier (used for contrast names built from group labels).
eds_safe <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

# Compute VST on a genomics SE — mirrors the legacy convert_genomics_se
# logic exactly: integer-cast counts; vst() for >=1000 rows, fall back to
# varianceStabilizingTransformation otherwise (and on any vst() error).
eds_compute_vst <- function(se) {
  forVST <- as.matrix(SummarizedExperiment::assay(se, "counts"))
  mode(forVST) <- "integer"
  tryCatch({
    if (nrow(forVST) < 1000) {
      DESeq2::varianceStabilizingTransformation(forVST, blind = FALSE)
    } else {
      DESeq2::vst(forVST, blind = FALSE)
    }
  }, error = function(e) {
    DESeq2::varianceStabilizingTransformation(forVST, blind = FALSE)
  })
}

# Compute RPKM/TPM on a genomics SE from raw counts + rowData()$featWidth.
# Self-contained (bulk-data formula only) so exploreDE doesn't need to depend
# on ezRun just for these two functions -- mirrors ezRun::getRpkm()/getTpm()'s
# bulk-data path exactly (ezRun's own single-cell/10x TPM branch is irrelevant
# here: exploreDE never sets metadata(se)$param$scProtocol). Also mirrors
# ExploreDatasetH5AD.R's private$.compute_fpkm()/.compute_tpm(), which do the
# same sweep on a plain matrix + featWidth vector instead of an SE.
eds_compute_rpkm <- function(se) {
  counts <- SummarizedExperiment::assay(se, "counts")
  featWidth <- SummarizedExperiment::rowData(se)$featWidth
  libSize <- Matrix::colSums(counts)
  rpkm <- sweep(counts * 1e9, MARGIN = 1, STATS = featWidth, FUN = "/")
  sweep(rpkm, MARGIN = 2, STATS = libSize, FUN = "/")
}

eds_compute_tpm <- function(se) {
  counts <- SummarizedExperiment::assay(se, "counts")
  featWidth <- SummarizedExperiment::rowData(se)$featWidth
  tpm <- sweep(counts, MARGIN = 1, STATS = featWidth, FUN = "/")
  sweep(tpm * 1e6, MARGIN = 2, STATS = Matrix::colSums(tpm), FUN = "/")
}

# Apply fid_lookup remap to a matrix's rownames. Drops rows whose source
# rowname isn't in the lookup.
.eds_remap_rownames <- function(m, lk) {
  if (is.null(lk)) return(m)
  new_rn <- unname(lk[rownames(m)])
  keep   <- !is.na(new_rn)
  m      <- m[keep, , drop = FALSE]
  rownames(m) <- new_rn[keep]
  m
}

# Resolve a path / URL / directory to an in-memory SE. Mirrors the load
# behaviour of the existing app's SE-fetching code:
#   - http(s) URL                  → readRDS(url(...))
#   - .RData file (Mageck EzResult) → strip suffix, fall through to directory
#   - .rds file                    → readRDS
#   - .qs2 file                    → qs2::qs_read
#   - directory                    → search for the standard SUSHI artefacts
eds_load_se <- function(input) {
  if (grepl("^https?://", input)) {
    return(readRDS(url(input)))
  }
  if (grepl("\\.RData$", input, ignore.case = TRUE)) {
    dir_input <- gsub("/(result|counts)-.*-EzResult\\.RData$", "", input)
    return(eds_load_se(dir_input))
  }
  if (grepl("\\.rds$", input, ignore.case = TRUE) && file.exists(input)) {
    return(readRDS(input))
  }
  if (grepl("\\.qs2$", input, ignore.case = TRUE) && file.exists(input)) {
    return(qs2::qs_read(input, nthreads = 4))
  }
  if (dir.exists(input)) {
    for (fname in c("deResult.rds",      "deResult.qs2",
                    "rawData_final.rds", "rawData_final.qs2",
                    "rawData.rds",       "rawData.qs2")) {
      fp <- file.path(input, fname)
      if (file.exists(fp)) {
        if (grepl("\\.qs2$", fname)) return(qs2::qs_read(fp, nthreads = 4))
        return(readRDS(fp))
      }
    }
    stop("No known SE artefact in directory: ", input)
  }
  stop("Cannot interpret input: ", input)
}

# Flatten a (possibly nested) clusterProfiler result into ONE tidy data.frame.
#
# This table (uns$enrichment) is the stored form of enrichment results: the
# clusterProfiler objects are not written, because they are unreadable without
# DOSE (no wasm build) and ~16x larger than the tables they carry: measured on
# a real dataset, >90% of the payload is geneSets/geneList that the app never
# displays. .enr_support_from_native() keeps the small part the set-level
# plots need. See plans/exploredelight/enrichment_detach_plan.md.
#
# Column names are kept EXACTLY as clusterProfiler emits them (`p.adjust`,
# `GeneRatio`, `enrichmentScore`) rather than the snake_case the 2026-04-17
# AnnData spec proposed. Deliberate deviation: the goal is dropping the
# dependency, not renaming the vocabulary, and every consuming module already
# reads these names.
#
# `ontology` carries the leaf name eds_flatten_cp_results() synthesises
# ("BP", "MF", "BP/upGenes", ...), so the nesting survives as a column.
#
# Reads the result table via attributes(), never `@` or as.data.frame(): both
# of those are fine on the cluster but as.data.frame() is S4 dispatch and dies
# where DOSE is absent, and this helper must work in both places.
.enr_tidy_from_native <- function(obj, method = NA_character_, padj_max = NULL,
                                  prefix = "") {
  # `prefix` matters when the caller already flattened: eds_flatten_cp_results()
  # on a bare leaf names it "default", which would throw away the ontology the
  # caller knows (e.g. "BP/upGenes"). Pass it back in.
  leaves <- eds_flatten_cp_results(obj, prefix)
  if (!length(leaves)) return(NULL)
  tabs <- lapply(names(leaves), function(nm) {
    leaf <- leaves[[nm]]
    # A leaf is either a real clusterProfiler S4 object (result lives in its
    # `result` slot, read via attributes() so this works without DOSE) or
    # already a plain tidy data.frame (a BYO producer's own table) -- use it
    # directly rather than looking for a non-existent "result" attribute on it.
    r <- if (is.data.frame(leaf)) leaf else attributes(leaf)[["result"]]
    if (is.null(r) || !is.data.frame(r) || !nrow(r)) return(NULL)
    r$method   <- method
    r$ontology <- nm
    r
  })
  tabs <- Filter(Negate(is.null), tabs)
  if (!length(tabs)) return(NULL)
  # rbind.fill, not rbind: ORA and GSEA result frames carry different columns
  # (GeneRatio/Count vs NES/enrichmentScore) and both must survive.
  out <- plyr::rbind.fill(tabs)
  if (!is.null(padj_max) && "p.adjust" %in% colnames(out))
    out <- out[!is.na(out$p.adjust) & out$p.adjust <= padj_max, , drop = FALSE]
  if (!nrow(out)) return(NULL)
  rownames(out) <- NULL
  out
}

# What a reader needs, beyond the tidy table, to rebuild clusterProfiler
# objects for the network / ridge / running-score plots. Written alongside
# the tidy table at uns$enrichment_support[[slot]][[kind]]:
#   class     : "enrichResult" or "gseaResult"
#   GSEA only:
#   gene_sets : data.frame(ontology, ID, gene), the members of the STORED terms
#               only (0.3 MB vs 10 MB for the full collection on real data)
#   gene_list : data.frame(ontology, gene, value), the ranked list in rank
#               order; ontology "*" when every leaf shares one list (usual)
#   exponent  : fgsea weighting exponent (the one param the plots read)
# ORA needs nothing extra: its network plot reads geneID from the table.
# Class is read with oldClass(), and slots via attributes(), so building this
# never needs DOSE loaded.
.enr_support_from_native <- function(obj, tidy) {
  leaves <- Filter(function(x) any(oldClass(x) %in% c("enrichResult", "gseaResult")),
                   eds_flatten_cp_results(obj))
  if (!length(leaves)) return(NULL)
  cls <- intersect(oldClass(leaves[[1]]), c("enrichResult", "gseaResult"))[1]
  out <- list(class = cls)
  if (cls != "gseaResult") return(out)

  sets <- lapply(names(leaves), function(ont) {
    ids <- tidy$ID[tidy$ontology == ont]
    gs  <- attributes(leaves[[ont]])[["geneSets"]][ids]
    gs  <- gs[!vapply(gs, is.null, logical(1))]
    if (!length(gs)) return(NULL)
    data.frame(ontology = ont, ID = rep(names(gs), lengths(gs)),
               gene = unlist(gs, use.names = FALSE), stringsAsFactors = FALSE)
  })
  out$gene_sets <- do.call(rbind, Filter(Negate(is.null), sets))

  lists <- lapply(leaves, function(l) attributes(l)[["geneList"]])
  lists <- Filter(function(gl) length(gl) > 0L, lists)
  if (length(lists)) {
    shared <- all(vapply(lists, identical, logical(1), lists[[1]]))
    use <- if (shared) setNames(lists[1], "*") else lists
    out$gene_list <- do.call(rbind, lapply(names(use), function(ont)
      data.frame(ontology = ont, gene = names(use[[ont]]), value = unname(use[[ont]]),
                 stringsAsFactors = FALSE)))
  }
  out$exponent <- (attributes(leaves[[1]])[["params"]] %||% list())$exponent %||% 1
  out
}

# Rebuild clusterProfiler objects from a stored tidy table + its support block
# (see .enr_support_from_native()). Returns a list of objects named by leaf
# ("BP/upGenes", "BP", or "default" for a bare object), or list() when DOSE is
# not installed (exploreDElight / wasm), in which case the three set-level
# plots are simply unavailable and everything else runs off the tidy table.
.enr_rebuild <- function(tidy, support) {
  if (is.null(tidy) || !nrow(tidy) || is.null(support$class)) return(list())
  if (!requireNamespace("DOSE", quietly = TRUE)) return(list())
  gsea <- identical(support$class, "gseaResult")
  onts <- unique(tidy$ontology)
  setNames(lapply(onts, function(ont) {
    r <- tidy[tidy$ontology == ont, setdiff(colnames(tidy), c("method", "ontology")), drop = FALSE]
    r <- r[, vapply(r, function(col) !all(is.na(col)), logical(1)), drop = FALSE]  # rbind.fill padding
    rownames(r) <- r$ID
    if (!gsea) {
      return(methods::new("enrichResult", result = r, readable = FALSE,
                          keytype = "UNKNOWN", ontology = "UNKNOWN", organism = "UNKNOWN",
                          pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1))
    }
    gs <- support$gene_sets
    gs <- gs[gs$ontology == ont, , drop = FALSE]
    gl <- support$gene_list
    gl <- gl[gl$ontology %in% c(ont, "*"), , drop = FALSE]
    methods::new("gseaResult", result = r, readable = FALSE,
                 geneSets = split(gs$gene, factor(gs$ID, levels = unique(gs$ID))),
                 geneList = stats::setNames(gl$value, gl$gene),
                 params = list(exponent = support$exponent %||% 1, pAdjustMethod = "BH"),
                 keytype = "UNKNOWN", setType = "UNKNOWN", organism = "UNKNOWN")
  }), onts)
}
