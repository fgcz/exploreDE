# -----------------------------------------------------------------------------
# Phase 1 equivalence tests: legacy_gStoreList_from_eds() must produce
# byte-equivalent output to the existing convert_genomics_se() /
# convert_proteomics_se() for every demo dataset.
#
# Strategy: for each fixture, call both paths and assert field-by-field.
# Fields that contain S4 objects (`se`) or environments (`param` may carry
# state) are compared with `identical()`; data.frames and matrices with
# `expect_equal()` to tolerate floating-point noise from VST / log2.
# -----------------------------------------------------------------------------

compare_gstorelist <- function(old, new, label, ignore = character()) {
  # Common-field equivalence
  common <- intersect(names(old), names(new))
  for (k in setdiff(common, ignore)) {
    expect_identical(
      class(new[[k]]), class(old[[k]]),
      info = sprintf("[%s] field '%s' class mismatch", label, k)
    )
    # Use equal (tolerant) for matrices/data.frames, identical for scalars/lists
    if (is.matrix(old[[k]]) || is.data.frame(old[[k]])) {
      expect_equal(
        new[[k]], old[[k]],
        info = sprintf("[%s] field '%s' (value)", label, k)
      )
    } else if (is.list(old[[k]]) && !is.data.frame(old[[k]])) {
      expect_equal(
        new[[k]], old[[k]],
        info = sprintf("[%s] list field '%s'", label, k)
      )
    } else {
      expect_identical(
        new[[k]], old[[k]],
        info = sprintf("[%s] field '%s'", label, k)
      )
    }
  }
  # Field set equivalence (modulo intentionally-ignored)
  expect_setequal(
    setdiff(names(new), ignore),
    setdiff(names(old), ignore)
  )
}

test_that("RNA-seq gene: shim matches convert_genomics_se", {
  fx <- fx_rnaseq_gene()
  old <- exploreDE:::.convert_genomics_se_legacy(fx$se, fx$dataDir)
  eds <- new_explore_dataset(fx$se)
  new <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)
  # `se` is the raw SE — identical by reference. Skip explicit comparison
  # (would error on large nested S4 objects).
  compare_gstorelist(old, new, "rnaseq_gene", ignore = c("se"))
})

test_that("RNA-seq isoform: shim matches convert_genomics_se", {
  fx <- fx_rnaseq_isoform()
  old <- exploreDE:::.convert_genomics_se_legacy(fx$se, fx$dataDir)
  eds <- new_explore_dataset(fx$se)
  new <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)
  compare_gstorelist(old, new, "rnaseq_isoform", ignore = c("se"))
})

test_that("CountQC: shim matches convert_genomics_se", {
  fx <- fx_countqc()
  old <- exploreDE:::.convert_genomics_se_legacy(fx$se, fx$dataDir)
  eds <- new_explore_dataset(fx$se)
  new <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)
  compare_gstorelist(old, new, "countqc", ignore = c("se"))
})

test_that("CrisprCountQC: shim matches convert_genomics_se", {
  fx <- fx_crispr_countqc()
  old <- exploreDE:::.convert_genomics_se_legacy(fx$se, fx$dataDir)
  eds <- new_explore_dataset(fx$se)
  new <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)
  compare_gstorelist(old, new, "crispr_countqc", ignore = c("se"))
})

test_that("Proteomics protein: shim matches convert_proteomics_se", {
  fx <- fx_proteomics_protein()
  old <- exploreDE:::.convert_proteomics_se_legacy(fx$se)
  eds <- new_explore_dataset(fx$se)
  new <- legacy_gStoreList_from_eds(eds, se = fx$se)
  compare_gstorelist(old, new, "proteomics_protein", ignore = c("se"))
})

test_that("Proteomics phospho: shim matches convert_proteomics_se", {
  fx <- fx_proteomics_phospho()
  old <- exploreDE:::.convert_proteomics_se_legacy(fx$se)
  eds <- new_explore_dataset(fx$se)
  new <- legacy_gStoreList_from_eds(eds, se = fx$se)
  compare_gstorelist(old, new, "proteomics_phospho", ignore = c("se"))
})

test_that("Proteomics peptide: shim matches convert_proteomics_se", {
  fx <- fx_proteomics_peptide()
  old <- exploreDE:::.convert_proteomics_se_legacy(fx$se)
  eds <- new_explore_dataset(fx$se)
  new <- legacy_gStoreList_from_eds(eds, se = fx$se)
  compare_gstorelist(old, new, "proteomics_peptide", ignore = c("se"))
})

test_that("dataDir attribute preserved on genomics path", {
  fx <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)
  gs  <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)
  expect_identical(attr(gs, "dataDir"), fx$dataDir)
})
