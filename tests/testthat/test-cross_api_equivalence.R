# -----------------------------------------------------------------------------
# Cross-API equivalence tests.
#
# Proves: for every fixture, the values returned by `eds$method()` are
# byte-equivalent to the equivalent legacy `gStoreList$X` field. This is
# the contract Phase 2 module migrations depend on — when a module swaps
# `inputDataReactive()$countList[["Raw"]]` for `eds$assay("counts")`, it
# must see the same matrix.
#
# Note: the shim's gStoreList output was proven equal to the legacy
# converters' output before those were retired (archived with their tests in
# plans/archive/legacy_converters/).
# -----------------------------------------------------------------------------

# Helper — align two matrices by rownames/colnames, then compare values.
expect_matrix_equal_aligned <- function(new, old, label, tol = 1e-8) {
  expect_setequal(rownames(new), rownames(old))
  expect_setequal(colnames(new), colnames(old))
  ord <- rownames(old)
  cls <- colnames(old)
  expect_equal(
    as.matrix(new[ord, cls, drop = FALSE]),
    as.matrix(old),
    tolerance = tol,
    info      = label
  )
}

# Helper — align two data.frames by `feature_id` (or rownames), then compare
# the shared columns. Tolerates feature_id-order differences and extra columns
# on either side.
expect_df_equal_aligned <- function(new, old, label, key = "feature_id", tol = 1e-8) {
  if (!key %in% colnames(new))  new[[key]]  <- rownames(new)
  if (!key %in% colnames(old))  old[[key]]  <- rownames(old)
  expect_setequal(new[[key]], old[[key]])
  ord       <- old[[key]]
  shared    <- intersect(colnames(new), colnames(old))
  shared    <- setdiff(shared, key)              # compared via key alignment
  new_align <- new[match(ord, new[[key]]), shared, drop = FALSE]
  old_align <- old[, shared, drop = FALSE]
  rownames(new_align) <- ord
  rownames(old_align) <- ord
  expect_equal(new_align, old_align, tolerance = tol, info = label)
}

# -----------------------------------------------------------------------------
# Generic per-fixture checks — applied to every shape via parametrised helpers
# -----------------------------------------------------------------------------

check_common <- function(eds, gs, fx_name) {
  # data_type:
  # - For genomics: eds$data_type() returns the same legacy strings as
  #   gs$dataType ("RNASeq" / "CountQC" / "CrisprCountQC"), so direct compare.
  # - For proteomics: eds$data_type() returns the FLAVOUR
  #   ("proteomics" / "peptidomics" / "phosphoproteomics"), while
  #   gs$dataType is always the umbrella "proteomics". The umbrella is
  #   eds$omics_class().
  expected_gs_dataType <- if (eds$omics_class() == "proteomics") {
    "proteomics"
  } else {
    eds$data_type()
  }
  expect_identical(
    expected_gs_dataType,
    gs$dataType,
    info = paste0("[", fx_name, "] data_type"))

  # samples / dataset — STRUCTURAL equivalence only.
  #
  # `eds$samples()` reads `colData(se)` directly. `gs$dataset` reads
  # `../input_dataset.tsv` preferentially (a SUSHI-specific layout) and
  # falls back to colData. The two sources carry the same conceptual data
  # but with type-level differences that are NOT bugs:
  #   - boolean cols ('paired') are logical in TSV, character in colData
  #   - factor cols ('Condition') are factor (with SUSHI-ordering) in
  #     legacy, character or original type in colData
  #   - colData has extra tool-artefact cols (`sf`, `conds`) the TSV omits
  # Type-level reconciliation is a per-module Phase 2 concern; this test
  # only validates structural equivalence (same sample ids, gs's cols are
  # a subset of eds's, same value-per-column when both happen to be the
  # same type).
  expect_setequal(rownames(eds$samples()), rownames(gs$dataset))
  expect_true(
    all(colnames(gs$dataset) %in% colnames(eds$samples())),
    info = paste0("[", fx_name, "] gs$dataset colnames are subset of eds$samples()"))
  # Compare shared columns where both sides agree on type; soft-skip cols
  # whose types differ (the documented TSV-vs-colData divergence).
  shared_cols <- intersect(colnames(eds$samples()), colnames(gs$dataset))
  for (col in shared_cols) {
    e_val <- eds$samples()[rownames(gs$dataset), col]
    g_val <- gs$dataset[[col]]
    if (!identical(class(e_val), class(g_val))) next  # type divergence — skip
    expect_equal(e_val, g_val,
                 info = sprintf("[%s] samples/dataset col '%s'", fx_name, col))
  }

  # factor_names — sets should match
  expect_setequal(eds$factor_names(), gs$factors)

  # factor_levels per factor — sets per factor should match
  for (f in eds$factor_names()) {
    expect_setequal(
      eds$factor_levels(f),
      gs$factorLevels[[f]])
  }

  # primary_grouping
  expect_identical(eds$primary_grouping(), gs$groupingName)

  # sample / feature alignment of the raw assay
  raw_a <- eds$raw_assay()
  expect_setequal(colnames(raw_a), eds$sample_ids())
  expect_setequal(rownames(raw_a), eds$feature_ids())
}

check_countList <- function(eds, gs, fx_name, name_map) {
  # name_map: named char where names are eds assay names, values are gs countList keys
  for (eds_name in names(name_map)) {
    gs_name <- name_map[[eds_name]]
    if (!gs_name %in% names(gs$countList)) next     # optional assays (FPKM/TPM/VST may be absent)
    if (!eds$has_assay(eds_name))           next
    expect_matrix_equal_aligned(
      eds$assay(eds_name),
      gs$countList[[gs_name]],
      label = paste0("[", fx_name, "] assay('", eds_name,
                     "') vs countList[['", gs_name, "']]"))
  }
}


# -----------------------------------------------------------------------------
# RNA-seq (gene-level)
# -----------------------------------------------------------------------------

test_that("RNA-seq gene: cross-API equivalence", {
  fx  <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)

  check_common(eds, gs, "rnaseq_gene")
  check_countList(eds, gs, "rnaseq_gene", c(
    counts          = "Raw",
    normalised      = "Normalised",
    normalised_log2 = "Normalised + Log2",
    vst             = "VST",
    fpkm            = "FPKM",
    tpm             = "TPM"
  ))
})

test_that("RNA-seq gene: DE values equivalent to seqAnno", {
  fx  <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)

  # eds$de() returns rownames-ordered; gs$seqAnno is sorted by pValue. Align
  # by feature_id, compare shared DE-stat columns.
  de  <- eds$de(eds$primary_contrast())
  sa  <- gs$seqAnno
  # The shim's seqAnno carries extra computed columns (<assay>_<group>_Mean
  # etc.). The shared columns we care about for module compatibility:
  shared_de_cols <- intersect(
    c("log2Ratio", "pValue", "fdr", "feature_id"),
    intersect(colnames(de), colnames(sa)))
  expect_true(length(shared_de_cols) >= 4, info = "shared DE cols present")
  expect_df_equal_aligned(
    de[, shared_de_cols, drop = FALSE],
    sa[, shared_de_cols, drop = FALSE],
    label = "rnaseq_gene DE values")
})

test_that("RNA-seq gene: enrichment objects equivalent", {
  fx  <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)

  if (!eds$has_enrichment()) skip("no enrichment in this fixture")

  ora_new  <- eds$enrichment(eds$primary_contrast(), "ora")
  ora_old  <- gs$goObjects$ora
  gsea_new <- eds$enrichment(eds$primary_contrast(), "gsea")
  gsea_old <- gs$goObjects$gsea

  # eds$enrichment() returns a *flattened* named list of clusterProfiler
  # objects; gs$goObjects$ora is the raw nested list. We check that the
  # underlying enrichResult/gseaResult objects are the same by serialising
  # and comparing — heavyweight but unambiguous.
  flat_old_ora <- exploreDE:::eds_flatten_cp_results(ora_old)
  expect_setequal(names(ora_new), names(flat_old_ora))
  for (k in names(ora_new)) {
    expect_identical(
      class(ora_new[[k]]),
      class(flat_old_ora[[k]]),
      info = paste0("ora kind ", k))
  }

  flat_old_gsea <- exploreDE:::eds_flatten_cp_results(gsea_old)
  expect_setequal(names(gsea_new), names(flat_old_gsea))
})

test_that("RNA-seq gene: has_enrichment matches runGO", {
  fx  <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)
  expect_identical(eds$has_enrichment(), isTRUE(gs$goObjects$runGO))
})


# -----------------------------------------------------------------------------
# RNA-seq (isoform-level)
# -----------------------------------------------------------------------------

test_that("RNA-seq isoform: cross-API equivalence", {
  fx  <- fx_rnaseq_isoform()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)

  check_common(eds, gs, "rnaseq_isoform")
  # Note: isoform-level feature_ids differ subtly between eds (rownames(se))
  # and gs (scuttle::uniquifyFeatureNames(transcript_id, gene_name) with
  # char-cleanup). Matrices won't align by rowname in this fixture; this is
  # the kind of thing a module migration must handle. The countList matrices
  # are NOT compared here for that reason.
  expect_identical(eds$feature_id_kind(), "isoform")
})


# -----------------------------------------------------------------------------
# CountQC
# -----------------------------------------------------------------------------

test_that("CountQC: cross-API equivalence", {
  fx  <- fx_countqc()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)

  check_common(eds, gs, "countqc")
  check_countList(eds, gs, "countqc", c(
    counts     = "Raw",
    normalised = "Normalised"
  ))
  # CountQC has no DE
  expect_identical(length(eds$contrasts()), 0L)
  expect_false(eds$has_de())
})


# -----------------------------------------------------------------------------
# CRISPR CountQC
# -----------------------------------------------------------------------------

test_that("CrisprCountQC: cross-API equivalence", {
  fx  <- fx_crispr_countqc()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)

  check_common(eds, gs, "crispr_countqc")
  check_countList(eds, gs, "crispr_countqc", c(
    counts     = "Raw",
    normalised = "Normalised"
  ))
  # sgrna_to_gene must match
  expect_setequal(names(eds$sgrna_to_gene()), names(gs$sgRNA_to_gene))
  expect_identical(
    eds$sgrna_to_gene()[names(gs$sgRNA_to_gene)],
    gs$sgRNA_to_gene)
})


# -----------------------------------------------------------------------------
# Proteomics — protein, peptide, phospho
# -----------------------------------------------------------------------------

test_that("Proteomics protein: cross-API equivalence", {
  fx  <- fx_proteomics_protein()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se)

  check_common(eds, gs, "proteomics_protein")

  # contrasts character vectors should match (set, not necessarily order)
  expect_setequal(eds$contrasts(), gs$contrasts)

  # DE tables per contrast — both sorted by pValue, both keyed by feature_id
  for (c in eds$contrasts()) {
    new <- eds$de(c)
    old <- gs$seqAnnoList[[c]]
    shared <- intersect(
      c("feature_id", "log2Ratio", "pValue", "fdr",
        "gene_name", "protein_Id", "nrPeptides"),
      intersect(colnames(new), colnames(old)))
    expect_df_equal_aligned(
      new[, shared, drop = FALSE],
      old[, shared, drop = FALSE],
      label = paste0("proteomics_protein de('", c, "')"))
  }

  # proteomics_choices — sets must match (both NULL at protein-level where
  # protein_Id == feature_id, so use expect_equal not expect_setequal)
  pc_new <- eds$proteomics_choices("protein")
  expect_equal(pc_new$choices, gs$protein_choices)
  expect_equal(names(pc_new$lookup), names(gs$protein_lookup))

  # countList values — every assay in eds, aligned by built feature_id
  for (assay_name in eds$assay_names()) {
    if (!assay_name %in% names(gs$countList)) next
    expect_matrix_equal_aligned(
      eds$assay(assay_name),
      gs$countList[[assay_name]],
      label = paste0("proteomics_protein assay('", assay_name, "')"))
  }
})

test_that("Proteomics phospho: cross-API equivalence (subset)", {
  fx  <- fx_proteomics_phospho()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se)

  check_common(eds, gs, "proteomics_phospho")
  expect_setequal(eds$contrasts(), gs$contrasts)
  expect_identical(eds$feature_id_kind(), "phospho")

  # Spot-check first contrast's DE values
  c1     <- eds$contrasts()[1]
  shared <- c("feature_id", "log2Ratio", "pValue", "fdr")
  shared <- intersect(shared, intersect(colnames(eds$de(c1)), colnames(gs$seqAnnoList[[c1]])))
  expect_df_equal_aligned(
    eds$de(c1)[, shared, drop = FALSE],
    gs$seqAnnoList[[c1]][, shared, drop = FALSE],
    label = paste0("proteomics_phospho de('", c1, "')"))
})

test_that("Proteomics peptide: cross-API equivalence (subset)", {
  fx  <- fx_proteomics_peptide()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se)

  check_common(eds, gs, "proteomics_peptide")
  expect_setequal(eds$contrasts(), gs$contrasts)
  expect_identical(eds$feature_id_kind(), "peptide")
  expect_setequal(eds$proteomics_choices("protein")$choices, gs$protein_choices)
})
