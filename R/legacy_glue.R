# =============================================================================
# Phase 1 compatibility shim: build a legacy `gStoreList` from an
# `ExploreDataset` instance.
#
# Goal: produce byte-equivalent output to the existing
# `convert_genomics_se()` / `convert_proteomics_se()` so the running app
# can be switched to the new infrastructure under the hood while every
# server module continues to read `inputDataReactive()$X` unchanged.
#
# Phase 2 migrates each server module to call `eds$method()` directly;
# once every module is migrated, this shim becomes dead code.
# =============================================================================

#' Build a legacy `gStoreList` from an `ExploreDataset` instance
#'
#' @param eds An `ExploreDataset` instance.
#' @param se The raw `SummarizedExperiment` (defaults to `eds$backend()`).
#' @param dataDir Optional path / URL — recorded as an attribute on the
#'   returned list and used for the path-based dataType heuristic in the
#'   genomics branch (`ExploreMageckCounts` / `CountQC` / RNA-seq).
#' @param param Optional SUSHI param list; defaults to
#'   `metadata(se)$param`.
#' @return A named list with the legacy gStoreList shape.
#' @export
legacy_gStoreList_from_eds <- function(eds, se = NULL, dataDir = NULL, param = NULL) {
  stopifnot(inherits(eds, "ExploreDataset"))
  if (is.null(se))    se    <- eds$backend()
  if (is.null(param)) param <- S4Vectors::metadata(se)$param %||% list()

  if (eds$omics_class() == "genomics") {
    .legacy_gStoreList_genomics(eds, se, dataDir, param)
  } else if (eds$omics_class() == "proteomics") {
    .legacy_gStoreList_proteomics(eds, se)
  } else {
    stop("legacy_gStoreList_from_eds: unsupported omics_class '",
         eds$omics_class(), "'")
  }
}

# ---------------------------------------------------------------------------
# Genomics path (RNA-seq / CountQC / CrisprCountQC)
# ---------------------------------------------------------------------------

.legacy_gStoreList_genomics <- function(eds, se, dataDir, param) {

  # ---- seqAnno + feature_id (legacy scuttle::uniquifyFeatureNames) --------
  rd <- suppressWarnings(data.frame(
    SummarizedExperiment::rowData(se),
    row.names   = rownames(se),
    check.names = FALSE
  ))
  if ("description" %in% colnames(rd)) {
    rd$description <- gsub(" \\[..*", "", rd$description)
  } else {
    rd$description <- rd$gene_name
  }
  rd$feature_id <- .legacy_feature_id_genomics(rd, param)

  # ---- dataset (with `input_dataset.tsv` fallback) ------------------------
  dataset <- .legacy_dataset_genomics(eds, se, dataDir, param)

  # Factors with [Factor] suffix already stripped via eds$samples() ordering;
  # rebuild explicitly here to preserve the legacy ordering (factors come from
  # the `[Factor]`-tagged colData columns, not from eds$factor_names() which
  # has already cleaned the names).
  colData_raw <- data.frame(SummarizedExperiment::colData(se), check.names = FALSE)
  factors_raw <- colnames(colData_raw)[grep("\\[Factor]", colnames(colData_raw))]
  factors     <- gsub(" \\[.*", "", factors_raw)
  factors     <- gsub(" ", "_", factors)
  factorLevels <- setNames(
    lapply(factors, function(f) levels(as.factor(dataset[[f]]))),
    factors
  )

  if (is.null(param$groupingName) || !nzchar(param$groupingName)) {
    param$groupingName <- factors[[1]]
  }

  # ---- countList ----------------------------------------------------------
  countList <- list()
  countList[["Raw"]] <- SummarizedExperiment::assay(se, "counts")
  normCounts <- grep("^signal$|^xNorm$", names(SummarizedExperiment::assays(se)), value = TRUE)
  if (length(normCounts)) {
    countList[["Normalised"]]          <- SummarizedExperiment::assay(se, normCounts)
    countList[["Normalised + Log2"]]   <- log2(SummarizedExperiment::assay(se, normCounts) + (param$backgroundExpression %||% 10))
  }
  if (!is.null(param$refBuild)) {
    countList[["FPKM"]] <- ezRun::getRpkm(se)
    countList[["TPM"]]  <- ezRun::getTpm(se)
  }
  forVST <- as.matrix(SummarizedExperiment::assay(se, "counts"))
  mode(forVST) <- "integer"
  countList[["VST"]] <- tryCatch({
    if (nrow(forVST) < 1000) {
      DESeq2::varianceStabilizingTransformation(forVST, blind = FALSE)
    } else {
      DESeq2::vst(forVST, blind = FALSE)
    }
  }, error = function(e) {
    DESeq2::varianceStabilizingTransformation(forVST, blind = FALSE)
  })

  # Align countList rownames to the legacy feature_id
  for (i in seq_along(countList)) {
    if (isTRUE(param$featureLevel == "isoform")) {
      idOrder <- intersect(rd$transcript_id, rownames(countList[[i]]))
      countList[[i]] <- countList[[i]][idOrder, ]
      rownames(countList[[i]]) <- rd$feature_id %>%
        gsub("[:/()-]", ".", .) %>%
        gsub(" ", ".", .)
    } else {
      idOrder <- intersect(rd$gene_id, rownames(countList[[i]]))
      countList[[i]] <- countList[[i]][idOrder, ]
      rownames(countList[[i]]) <- rd$feature_id
    }
  }

  # ---- GO pathways --------------------------------------------------------
  allPathwaysDEG <- c(
    unique(unlist(strsplit(x = rd$`GO BP`, "; "))),
    unique(unlist(strsplit(x = rd$`GO MF`, "; "))),
    unique(unlist(strsplit(x = rd$`GO CC`, "; ")))
  )
  goTerms     <- clusterProfiler::go2term(allPathwaysDEG)
  allPathways <- paste(goTerms[, 1], goTerms[, 2])

  # ---- dataType from path -------------------------------------------------
  # eds$data_type() already returns the legacy "RNASeq"/"CountQC"/"CrisprCountQC"
  # strings, so the path-based heuristic only matters when dataDir is given
  # (and could in principle disagree with what the SE looks like — defer to
  # the path).
  dataTypeGStore <- if (!is.null(dataDir) && grepl("ExploreMageckCounts", dataDir)) {
    "CrisprCountQC"
  } else if (!is.null(dataDir) && grepl("CountQC", dataDir)) {
    "CountQC"
  } else {
    eds$data_type()
  }

  genes <- rd$feature_id

  gStoreList <- list(
    dataType     = dataTypeGStore,
    se           = se,
    param        = param,
    seqAnno      = rd,
    dataset      = dataset,
    factors      = factors,
    factorLevels = factorLevels,
    groupingName = param$groupingName,
    countList    = countList,
    genes        = genes,
    allPathways  = allPathways
  )

  # ---- DE / contrast extras (RNA-seq only) -------------------------------
  if (dataTypeGStore %in% c("CountQC", "CrisprCountQC")) {
    gStoreList$design       <- dataTypeGStore
    gStoreList$designLevels <- NULL
    if (dataTypeGStore == "CrisprCountQC") {
      gStoreList$sgRNA_to_gene <- setNames(
        gsub("_((sgrna)?[0-9]+|[+-]_[^_]+)_[ACGTN]+$", "", rd$feature_id, ignore.case = TRUE),
        rd$feature_id
      )
    }
  } else if (dataTypeGStore == "RNASeq") {
    pValueColumn <- grep("^p[-_]?value$", colnames(rd), value = TRUE, ignore.case = TRUE)[[1]]
    rd           <- rd[order(rd[[pValueColumn]]), ]
    genes        <- rd$feature_id[order(rd$pValue)]

    dataset[[param$groupingName]] <- factor(
      dataset[[param$groupingName]],
      levels = c(
        param$sampleGroup,
        param$refGroup,
        unique(dataset[[param$groupingName]])[
          !unique(dataset[[param$groupingName]]) %in% c(param$sampleGroup, param$refGroup)
        ]
      )
    )

    # Per-group means / SDs appended to seqAnno
    for (count in names(countList)) {
      cts <- countList[[count]]
      rd[[paste0(count, "_mean")]] <- rowMeans(cts)[rd$feature_id]
      rd[[paste0(count, "_SD")]]   <- matrixStats::rowSds(cts)[rd$feature_id]
      isSample <- dataset[[param$groupingName]] == param$sampleGroup
      rd[[paste0(count, "_", param$sampleGroup, "_Mean")]] <-
        rowMeans(cts[, isSample, drop = FALSE])[rd$feature_id]
      isRef <- dataset[[param$groupingName]] == param$refGroup
      rd[[paste0(count, "_", param$refGroup, "_Mean")]] <-
        rowMeans(cts[, isRef, drop = FALSE])[rd$feature_id]
    }

    design       <- param$comparison
    designLevels <- strsplit(design, split = "--over--")[[1]]

    goObjects <- list(runGO = param$runGO)
    if (isTRUE(param$runGO)) {
      goObjects$ora  <- S4Vectors::metadata(se)$enrichResult
      goObjects$gsea <- S4Vectors::metadata(se)$enrichResultGSEA
    }

    gStoreList$goObjects    <- goObjects
    gStoreList$groupingName <- param$groupingName
    gStoreList$design       <- design
    gStoreList$designLevels <- designLevels
    gStoreList$seqAnno      <- rd
    gStoreList$dataset      <- dataset
    gStoreList$genes        <- genes
  }

  attr(gStoreList, "dataDir") <- dataDir
  gStoreList
}

.legacy_feature_id_genomics <- function(rd, param) {
  fid <- if (isTRUE(param$featureLevel == "isoform")) {
    scuttle::uniquifyFeatureNames(ID = rd$transcript_id, names = rd$gene_name)
  } else {
    scuttle::uniquifyFeatureNames(ID = rd$gene_id, names = rd$gene_name)
  }
  fid[grep("^_", fid)] <- gsub("_", "", fid[grep("^_", fid)])
  fid
}

.legacy_dataset_genomics <- function(eds, se, dataDir, param) {
  ds <- if (!is.null(dataDir) && file.exists(file.path(dataDir, "../input_dataset.tsv"))) {
    data.frame(
      readr::read_tsv(file.path(dataDir, "../input_dataset.tsv"), show_col_types = FALSE),
      check.names = FALSE
    )
  } else {
    data.frame(SummarizedExperiment::colData(se), check.names = FALSE)
  }
  if (length(param$samples) > 0) {
    ds <- ds[which(ds$Name %in% param$samples), ]
  }
  ds[["Name"]] <- factor(ds[["Name"]], levels = ds[["Name"]])
  rownames(ds) <- ds$Name
  colnames(ds) <- gsub(" \\[.*", "", colnames(ds))
  colnames(ds) <- gsub(" ", "_", colnames(ds))

  # Deduplicate identically-valued duplicate columns (legacy quirk)
  dups <- which(duplicated(names(ds)) | duplicated(names(ds), fromLast = TRUE))
  dup_names <- unique(names(ds)[dups])
  for (dup in dup_names) {
    cols <- which(names(ds) == dup)
    keep <- cols[1]
    drop_these <- integer(0)
    for (i in 2:length(cols)) {
      if (all(ds[[cols[i]]] == ds[[keep]], na.rm = TRUE)) {
        drop_these <- c(drop_these, cols[i])
      } else {
        names(ds)[cols[i]] <- paste0(dup, ".", i - 1)
      }
    }
    if (length(drop_these) > 0) ds <- ds[, -drop_these, drop = FALSE]
  }
  ds
}

# ---------------------------------------------------------------------------
# Multi-DE per-contrast objects from an ExploreDataset (h5ad uploads)
# ---------------------------------------------------------------------------
#
# The Multi-DE subsystem was built to compare several *separate* single-contrast
# SUSHI deResult SEs, one gStoreList per contrast (see
# `convert_genomics_se()` in server-multi-data-selection.R). An uploaded
# multi-contrast `.h5ad` is a *single* ExploreDataset holding N contrasts, so
# there is no per-contrast SE to convert. `legacy_multi_contrast_from_eds()`
# reconstructs the same per-contrast gStoreList shape those modules consume,
# but sourced entirely from `eds` accessors — one call per contrast. This is
# the genomics analogue of the proteomics auto-load in
# server-multi-data-selection.R (which already fans a single proteomics SE out
# into one pseudo-object per contrast).

# Reshape the flat enrichment accessor list — keys "ont/sig" (SE backend) or
# "ora/ont/sig" / "gsea/ont" (h5ad backend) — back into the nested
# ont -> sig (ORA) / ont (GSEA) shape the multi-DE modules read from legacy
# SUSHI metadata (`metadata(se)$enrichResult` / `$enrichResultGSEA`).
.reshape_enrichment_nested <- function(flat, kind = c("ora", "gsea")) {
  kind <- match.arg(kind)
  if (!length(flat)) return(list())
  out <- list()
  for (nm in names(flat)) {
    key   <- sub("^(ora|gsea)/", "", nm)
    parts <- strsplit(key, "/", fixed = TRUE)[[1]]
    if (kind == "gsea") {
      out[[parts[[1]]]] <- flat[[nm]]
    } else {
      ont <- parts[[1]]
      sig <- if (length(parts) >= 2L) parts[[2]] else "bothGenes"
      if (is.null(out[[ont]])) out[[ont]] <- list()
      out[[ont]][[sig]] <- flat[[nm]]
    }
  }
  out
}

# Append the `Label` column (ID + truncated Description) to each clusterProfiler
# result's @result table — mirrors the inline loop the SUSHI genomics loader
# runs, so the Upset / Tile / Pathway modules find the labels they plot.
# Tolerant of the real S4 enrichResult/gseaResult (slot access) and the S3
# stand-ins used in tests (list access); a no-op when there is no Description.
.label_enrichment_results <- function(nested) {
  add_label <- function(er) {
    if (is.null(er)) return(er)
    res   <- tryCatch(er@result, error = function(e) NULL)
    is_s4 <- !is.null(res)
    if (is.null(res)) res <- tryCatch(er[["result"]], error = function(e) NULL)
    if (is.null(res) || is.null(res$Description) || nrow(res) < 2L) return(er)
    res$Label <- vapply(seq_along(res$Description), function(j) {
      desc <- res$Description[j]
      if (nchar(desc) > 35) desc <- substr(desc, 1, 45)
      paste0(res$ID[j], "\n", desc)
    }, character(1))
    if (is_s4) er@result <- res else er[["result"]] <- res
    er
  }
  lapply(nested, function(x) if (methods::is(x, "list")) lapply(x, add_label) else add_label(x))
}

#' Build a legacy per-contrast Multi-DE object from an `ExploreDataset`
#'
#' Reconstructs the single-contrast gStoreList shape the Multi-DE modules
#' consume (`seqAnno`, `dataset`, `countList`, `param`, `design`,
#' `ora`/`gsea`, ...) for one contrast of a (possibly multi-contrast)
#' genomics dataset — used to fan an uploaded `.h5ad` out into the Multi-DE
#' data-selection list, the way the proteomics auto-load already does.
#'
#' @param eds An `ExploreDataset` instance (genomics).
#' @param contrast A single contrast name (must be in `eds$contrasts()`).
#' @return A named list with the legacy per-contrast Multi-DE shape.
#' @export
legacy_multi_contrast_from_eds <- function(eds, contrast) {
  stopifnot(inherits(eds, "ExploreDataset"), contrast %in% eds$contrasts())
  cm       <- eds$contrast_meta(contrast)
  grouping <- eds$primary_grouping()
  test     <- cm$test_group
  ref      <- cm$reference_group

  # ---- dataset (samples) --------------------------------------------------
  dataset <- eds$samples()
  if (!"Name" %in% colnames(dataset)) dataset$Name <- rownames(dataset)
  # Order the grouping factor so the contrast's test/reference groups come
  # first (legacy convention used by the boxplot bracket logic).
  if (!is.na(grouping) && grouping %in% colnames(dataset)) {
    lv      <- unique(as.character(dataset[[grouping]]))
    ordered <- c(intersect(c(test, ref), lv), setdiff(lv, c(test, ref)))
    dataset[[grouping]] <- factor(as.character(dataset[[grouping]]), levels = ordered)
  }

  factors      <- eds$factor_names()
  factorLevels <- setNames(lapply(factors, function(f) eds$factor_levels(f)), factors)

  # ---- seqAnno = feature annotation + this contrast's DE stats ------------
  de    <- eds$de(contrast)          # feature_id, log2Ratio, pValue, fdr, gene_name, ...
  feats <- eds$features()            # + GO BP/MF/CC, label_*, display_name
  extra <- setdiff(colnames(feats), colnames(de))
  seqAnno <- de
  if (length(extra))
    seqAnno <- cbind(seqAnno, feats[de$feature_id, extra, drop = FALSE])
  if ("pValue" %in% colnames(seqAnno))
    seqAnno <- seqAnno[order(seqAnno$pValue, na.last = TRUE), , drop = FALSE]
  genes <- seqAnno$feature_id

  # ---- countList (legacy display names) -----------------------------------
  assay_alias <- c(counts = "Raw", normalised = "Normalised",
                   normalised_log2 = "Normalised + Log2",
                   vst = "VST", fpkm = "FPKM", tpm = "TPM")
  countList <- list()
  for (a in eds$assay_names()) {
    countList[[assay_alias[[a]] %||% a]] <- eds$assay(a)
  }

  # ---- param --------------------------------------------------------------
  grouping2 <- eds$secondary_grouping()
  if (is.na(grouping2)) grouping2 <- ""
  param <- list(
    comparison           = paste0(test, "--over--", ref),
    sampleGroup          = test,
    refGroup             = ref,
    groupingName         = grouping,
    grouping2Name        = grouping2,
    featureLevel         = if (identical(eds$feature_id_kind(), "isoform")) "isoform" else "gene",
    runGO                = eds$has_enrichment(),
    refBuild             = eds$ref_build(),
    ezRef                = NULL,
    backgroundExpression = eds$background_expression()
  )

  out <- list(
    dataType     = eds$data_type(),
    se           = NULL,
    param        = param,
    seqAnno      = seqAnno,
    dataset      = dataset,
    factors      = factors,
    factorLevels = factorLevels,
    groupingName = grouping,
    countList    = countList,
    genes        = genes,
    allPathways  = eds$all_pathways(),
    design       = paste0(test, "--over--", ref),
    designLevels = c(test, ref)
  )

  # ---- enrichment (nested ont -> sig / ont, with Label) -------------------
  if (eds$has_enrichment()) {
    ora  <- .reshape_enrichment_nested(eds$enrichment(contrast, "ora"),  "ora")
    gsea <- .reshape_enrichment_nested(eds$enrichment(contrast, "gsea"), "gsea")
    if (length(ora))  out$ora  <- .label_enrichment_results(ora)
    if (length(gsea)) out$gsea <- .label_enrichment_results(gsea)
  }

  out
}

# ---------------------------------------------------------------------------
# Proteomics path
# ---------------------------------------------------------------------------

.legacy_gStoreList_proteomics <- function(eds, se) {
  contrasts <- names(SummarizedExperiment::rowData(se))[
    grep("^constrast_", names(SummarizedExperiment::rowData(se)))
  ]
  contrasts <- gsub("constrast_", "", contrasts)

  contrast_dfs <- setNames(
    lapply(contrasts, function(con) {
      as.data.frame(SummarizedExperiment::rowData(se)[[paste0("constrast_", con)]])
    }),
    contrasts
  )

  seqAnnoList <- setNames(
    lapply(contrasts, function(con) {
      # .saint_standardise_de maps SAINT's log2_EFCs/BFDR -> log2Ratio/fdr and
      # adds an all-NA pValue (no-op for standard prolfqua frames); .de_order_idx
      # sorts by pValue when usable, else fdr — guarding order(NULL) for SAINT.
      sa <- .saint_standardise_de(standardise_columns(contrast_dfs[[con]]))
      sa <- sa[.de_order_idx(sa), ]
      sa$feature_id <- build_feature_id(sa)
      if ("protein_Id" %in% colnames(sa) && identical(sa$feature_id, sa$protein_Id)) {
        sa$protein_Id <- NULL
      }
      stats_to_keep <- c("log2Ratio", "pValue", "fdr", "nrPeptides")
      non_numeric   <- colnames(sa)[!vapply(sa, is.numeric, logical(1))]
      keep <- unique(c("feature_id", non_numeric, intersect(stats_to_keep, colnames(sa))))
      sa <- sa[, keep, drop = FALSE]
      if (!"description" %in% colnames(sa)) {
        sa$description <- if ("gene_name" %in% colnames(sa)) sa$gene_name else sa$feature_id
      }
      sa
    }),
    contrasts
  )

  # ---- dataset ------------------------------------------------------------
  dataset <- as.data.frame(SummarizedExperiment::colData(se))
  name_col <- colnames(dataset)[
    grep("^Sample$|^Name$|^sampleName$", colnames(dataset), ignore.case = TRUE)
  ]
  dataset$Name <- dataset[[name_col]]
  rownames(dataset) <- dataset[[name_col]]

  factors <- colnames(dataset)[
    grep("Sample|Name|names|channel|raw\\.file", colnames(dataset),
         invert = TRUE, ignore.case = TRUE)
  ]
  factorLevels <- setNames(
    lapply(factors, function(f) levels(as.factor(dataset[[f]]))),
    factors
  )
  for (x in factors) {
    dataset[[x]] <- factor(dataset[[x]],
                           levels = unique(gtools::mixedsort(dataset[[x]])))
  }

  # ---- countList (assay rownames remapped via fid_lookup) -----------------
  unsorted_sa <- standardise_columns(contrast_dfs[[contrasts[1]]])
  unsorted_sa$feature_id <- build_feature_id(unsorted_sa)
  raw_rn     <- rownames(contrast_dfs[[contrasts[1]]])
  fid_lookup <- setNames(unsorted_sa$feature_id, raw_rn)

  countList <- setNames(
    lapply(names(SummarizedExperiment::assays(se)), function(ass) {
      x <- SummarizedExperiment::assay(se, ass)
      mapped <- fid_lookup[rownames(x)]
      valid  <- !is.na(mapped)
      x      <- x[valid, , drop = FALSE]
      rownames(x) <- unname(mapped[valid])
      x
    }),
    names(SummarizedExperiment::assays(se))
  )

  genes <- unique(unlist(lapply(seqAnnoList, `[[`, "feature_id")))

  # ---- gene / protein lookups --------------------------------------------
  sa1 <- seqAnnoList[[1]]
  gene_lookup    <- NULL
  gene_choices   <- NULL
  protein_lookup <- NULL
  protein_choices <- NULL
  if ("gene_name" %in% colnames(sa1) && !identical(sa1$feature_id, sa1$gene_name)) {
    gene_lookup  <- split(sa1$feature_id, sa1$gene_name)
    gene_choices <- unique(sa1$gene_name)
  }
  if ("protein_Id" %in% colnames(sa1)) {
    protein_lookup  <- split(sa1$feature_id, sa1$protein_Id)
    protein_choices <- unique(sa1$protein_Id)
  }

  list(
    dataType        = "proteomics",
    se              = se,
    seqAnnoList     = seqAnnoList,
    contrasts       = contrasts,
    dataset         = dataset,
    factors         = factors,
    countList       = countList,
    genes           = genes,
    factorLevels    = factorLevels,
    groupingName    = factors[[1]],
    gene_lookup     = gene_lookup,
    protein_lookup  = protein_lookup,
    gene_choices    = gene_choices,
    protein_choices = protein_choices
  )
}
