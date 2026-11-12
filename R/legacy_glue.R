# =============================================================================
# Multi DE records: the per-contrast list shape the Multi DE modules read
# (`seqAnno`, `dataset`, `countList`, `param`, `ora` / `gsea`, ...), built
# from an `ExploreDataset` by legacy_multi_contrast_from_eds().
# =============================================================================

# ---------------------------------------------------------------------------
# Multi-DE per-contrast objects from an ExploreDataset (h5ad uploads)
# ---------------------------------------------------------------------------
#
# The Multi-DE subsystem was built to compare several *separate* single-contrast
# SUSHI deResult SEs, one gStoreList per contrast. A dataset holding N
# contrasts has no per-contrast SE to convert. `legacy_multi_contrast_from_eds()`
# reconstructs the same per-contrast gStoreList shape those modules consume,
# but sourced entirely from `eds` accessors — one call per contrast.

# Reshape the flat enrichment accessor list — keys "ont/sig" (SE backend) or
# "ora/ont/sig" / "gsea/ont" (h5ad backend) — back into the nested
# ont -> sig (ORA) / ont (GSEA) shape the multi-DE modules read from legacy
# SUSHI metadata (`metadata(se)$enrichResult` / `$enrichResultGSEA`).
.reshape_enrichment_nested <- function(flat, kind = c("ora", "gsea")) {
  kind <- match.arg(kind)
  if (!length(flat)) {
    return(list())
  }
  out <- list()
  for (nm in names(flat)) {
    key <- sub("^(ora|gsea)/", "", nm)
    parts <- strsplit(key, "/", fixed = TRUE)[[1]]
    if (kind == "gsea") {
      out[[parts[[1]]]] <- flat[[nm]]
    } else {
      ont <- parts[[1]]
      sig <- if (length(parts) >= 2L) parts[[2]] else "bothGenes"
      if (is.null(out[[ont]])) {
        out[[ont]] <- list()
      }
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
    # Only clusterProfiler objects carry a result table; data.frame leaves stay
    if (!isS4(er)) {
      return(er)
    }
    res <- er@result
    if (is.null(res$Description) || nrow(res) < 2L) {
      return(er)
    }
    res$Label <- vapply(
      seq_along(res$Description),
      function(j) {
        desc <- res$Description[j]
        if (nchar(desc) > 35) {
          desc <- substr(desc, 1, 45)
        }
        paste0(res$ID[j], "\n", desc)
      },
      character(1)
    )
    er@result <- res
    er
  }
  lapply(nested, function(x) if (methods::is(x, "list")) lapply(x, add_label) else add_label(x))
}

# The count methods of a dataset under their legacy display names (proteomics
# intensities keep their key), and which of them are on the log2 scale. Built
# once per dataset: every contrast's Multi DE record holds the same list.
multi_count_list <- function(eds) {
  assay_alias <- c(
    counts = "Raw",
    normalised = "Normalised",
    normalised_log2 = "Normalised + Log2",
    vst = "VST",
    fpkm = "FPKM",
    tpm = "TPM"
  )
  countList <- list()
  countLog2 <- logical()
  for (a in eds$assay_names()) {
    label <- if (a %in% names(assay_alias)) assay_alias[[a]] else a
    countList[[label]] <- eds$assay(a)
    countLog2[[label]] <- eds$assay_is_log2(a)
  }
  list(countList = countList, countLog2 = countLog2)
}

#' Build a legacy per-contrast Multi-DE object from an `ExploreDataset`
#'
#' Reconstructs the single-contrast gStoreList shape the Multi-DE modules
#' consume (`seqAnno`, `dataset`, `countList`, `param`, `design`,
#' `ora`/`gsea`, ...) for one contrast of a (possibly multi-contrast)
#' dataset, genomics or proteomics, from any backend. The Multi-DE auto-load
#' fans a loaded dataset out into one record per contrast with it.
#'
#' @param eds An `ExploreDataset` instance.
#' @param contrast A single contrast name (must be in `eds$contrasts()`).
#' @param counts The dataset's count methods, `multi_count_list(eds)`; the
#'   records of one dataset share them.
#' @return A named list with the legacy per-contrast Multi-DE shape.
#' @noRd
legacy_multi_contrast_from_eds <- function(eds, contrast, counts = multi_count_list(eds)) {
  stopifnot(inherits(eds, "ExploreDataset"), contrast %in% eds$contrasts())
  cm <- eds$contrast_meta(contrast)
  grouping <- eds$primary_grouping()
  test <- cm$test_group
  ref <- cm$reference_group

  # ---- dataset (samples) --------------------------------------------------
  dataset <- eds$samples()
  if (!"Name" %in% colnames(dataset)) {
    dataset$Name <- rownames(dataset)
  }
  # Order the grouping factor so the contrast's test/reference groups come
  # first (legacy convention used by the boxplot bracket logic).
  if (!is.na(grouping) && grouping %in% colnames(dataset)) {
    lv <- unique(as.character(dataset[[grouping]]))
    ordered <- c(intersect(c(test, ref), lv), setdiff(lv, c(test, ref)))
    dataset[[grouping]] <- factor(as.character(dataset[[grouping]]), levels = ordered)
  }

  factors <- eds$factor_names()
  factorLevels <- stats::setNames(lapply(factors, function(f) eds$factor_levels(f)), factors)

  # ---- seqAnno = feature annotation + this contrast's DE stats ------------
  de <- eds$de(contrast) # feature_id, log2Ratio, pValue, fdr, gene_name, ...
  feats <- eds$features() # + GO BP/MF/CC, label_*, display_name
  extra <- setdiff(colnames(feats), colnames(de))
  seqAnno <- de
  if (length(extra)) {
    seqAnno <- cbind(seqAnno, feats[de$feature_id, extra, drop = FALSE])
  }
  genes <- seqAnno$feature_id # in eds$de() order (eds_finalise_de())

  # ---- param --------------------------------------------------------------
  grouping2 <- eds$secondary_grouping()
  if (is.na(grouping2)) {
    grouping2 <- ""
  }
  param <- list(
    comparison = paste0(test, "--over--", ref),
    sampleGroup = test,
    refGroup = ref,
    groupingName = grouping,
    grouping2Name = grouping2,
    featureLevel = if (identical(eds$feature_id_kind(), "isoform")) "isoform" else "gene",
    runGO = eds$has_enrichment(),
    refBuild = eds$ref_build(),
    ezRef = NULL,
    backgroundExpression = eds$background_expression()
  )

  out <- list(
    # The Multi-DE modules test for "proteomics", whatever the level
    # (data_type() can be "phosphoproteomics" or "peptidomics").
    dataType = if (eds$omics_class() == "proteomics") "proteomics" else eds$data_type(),
    se = NULL,
    param = param,
    seqAnno = seqAnno,
    dataset = dataset,
    factors = factors,
    factorLevels = factorLevels,
    groupingName = grouping,
    countList = counts$countList,
    countLog2 = counts$countLog2,
    genes = genes,
    allPathways = eds$all_pathways(),
    design = paste0(test, "--over--", ref),
    designLevels = c(test, ref)
  )

  # ---- enrichment (nested ont -> sig / ont, with Label) -------------------
  if (eds$has_enrichment()) {
    ora <- .reshape_enrichment_nested(eds$enrichment(contrast, "ora"), "ora")
    gsea <- .reshape_enrichment_nested(eds$enrichment(contrast, "gsea"), "gsea")
    if (length(ora)) {
      out$ora <- .label_enrichment_results(ora)
    }
    if (length(gsea)) out$gsea <- .label_enrichment_results(gsea)
  }

  out
}
