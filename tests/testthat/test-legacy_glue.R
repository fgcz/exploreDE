# -----------------------------------------------------------------------------
# legacy_gStoreList_from_eds(): the gStoreList shape the multi-DE modules still
# consume. (Its byte-equivalence tests against the pre-2026-05 converters are
# archived in plans/archive/legacy_converters/, with those converters.)
# -----------------------------------------------------------------------------

test_that("dataDir attribute preserved on genomics path", {
  fx <- fx_rnaseq_gene()
  eds <- new_explore_dataset(fx$se)
  gs <- legacy_gStoreList_from_eds(eds, se = fx$se, dataDir = fx$dataDir)
  expect_identical(attr(gs, "dataDir"), fx$dataDir)
})

test_that("SAINT proteomics: legacy shim builds without error (order(NULL) guard)", {
  # The legacy gStoreList glue re-parses the SE and standardise_columns() does
  # not recognise SAINT's log2_EFCs/SaintScore/BFDR, so sa$pValue was NULL and
  # order(NULL) threw — aborting the load observer and hanging the app spinner.
  # No old-converter equivalence: the old converter breaks on SAINT too.
  fx <- fx_proteomics_saint()
  eds <- new_explore_dataset(fx$se)
  gs <- expect_error(legacy_gStoreList_from_eds(eds, se = fx$se), NA) # NA = no error
  sa <- gs$seqAnnoList[[1]]
  expect_true(all(c("log2Ratio", "pValue", "fdr") %in% colnames(sa)))
  expect_true(all(is.na(sa$pValue))) # SAINT has no raw p-value
  expect_true(any(!is.na(sa$fdr))) # but a real (BFDR) fdr
})
