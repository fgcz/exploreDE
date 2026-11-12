# -----------------------------------------------------------------------------
# ExploreDatasetSE_Proteomics tests (protein / peptide / phospho flavours)
# -----------------------------------------------------------------------------

test_that("Proteomics protein-level: subclass + populated accessors", {
  fx <- fx_proteomics_protein()
  eds <- new_explore_dataset(fx$se)

  expect_s3_class(eds, "ExploreDatasetSE_Proteomics")
  expect_identical(eds$omics_class(), "proteomics")
  expect_identical(eds$data_type(), "proteomics")
  expect_identical(eds$feature_id_kind(), "protein")
  expect_gt(length(eds$contrasts()), 0L)
  expect_true(eds$has_de())

  f <- eds$features()
  expect_true("feature_id" %in% colnames(f))
  expect_true("display_name" %in% colnames(f))
  # At protein-level, feature_id == protein_Id, so the column is dropped
  # (matches legacy convert_proteomics_se behaviour).
  expect_false("protein_Id" %in% colnames(f))
})

test_that("Proteomics protein-level: de() returns populated rows", {
  fx <- fx_proteomics_protein()
  eds <- new_explore_dataset(fx$se)

  c1 <- eds$primary_contrast()
  de <- eds$de(c1)
  expect_s3_class(de, "data.frame")
  expect_gt(nrow(de), 0L)
  expect_true("effect_log2fc" %in% colnames(de))
  expect_true("score_pvalue" %in% colnames(de))
  expect_true("score_fdr" %in% colnames(de))
  expect_true(any(!is.na(de$effect_log2fc)))
})

test_that("Proteomics protein-level: assay rownames remap via fid_lookup", {
  fx <- fx_proteomics_protein()
  eds <- new_explore_dataset(fx$se)

  m <- eds$assay("intensities")
  expect_true(all(rownames(m) %in% eds$feature_ids()))
})

test_that("Proteomics protein-level: proteomics_choices behaves correctly", {
  fx <- fx_proteomics_protein()
  eds <- new_explore_dataset(fx$se)

  # At protein-level: protein_Id == feature_id, so protein_choices is NULL
  # (legacy convention — no separate lookup needed).
  prot_choices <- eds$proteomics_choices("protein")
  expect_null(prot_choices$choices)
  expect_null(prot_choices$lookup)

  # Gene lookups should still be populated if gene_name is present.
  gene_choices <- eds$proteomics_choices("gene")
  expect_gt(length(gene_choices$choices), 0L)
  expect_true(is.list(gene_choices$lookup))
})

test_that("Proteomics protein-level: validator passes", {
  fx <- fx_proteomics_protein()
  eds <- new_explore_dataset(fx$se)
  expect_silent(validate_explore_dataset(eds))
})

test_that("Proteomics phospho: feature_id_kind + data_type", {
  fx <- fx_proteomics_phospho()
  eds <- new_explore_dataset(fx$se)
  expect_identical(eds$omics_class(), "proteomics")
  expect_identical(eds$data_type(), "phosphoproteomics")
  expect_identical(eds$feature_id_kind(), "phospho")

  c1 <- eds$primary_contrast()
  de <- eds$de(c1)
  expect_gt(nrow(de), 0L)
  expect_silent(validate_explore_dataset(eds))
})

test_that("Proteomics peptide: feature_id_kind + data_type", {
  fx <- fx_proteomics_peptide()
  eds <- new_explore_dataset(fx$se)
  expect_identical(eds$omics_class(), "proteomics")
  expect_identical(eds$data_type(), "peptidomics")
  expect_identical(eds$feature_id_kind(), "peptide")
  expect_silent(validate_explore_dataset(eds))
})

# --- SAINTexpress transitional SE-backend support ----------------------------
# The p42111 / wu_349205 AP-MS file reaches the app via the legacy ?data=<.rds>
# URL path. metadata$param is NULL; DE cols are log2_EFCs / SaintScore / BFDR
# with no raw p-value. These lock in the localised SE-backend hack that makes
# such files work until producers emit canonical h5ad.

test_that("SAINT SE: constructs as proteomics despite NULL metadata$param", {
  fx <- fx_proteomics_saint()
  expect_true(is.null(S4Vectors::metadata(fx$se)$param))
  eds <- new_explore_dataset(fx$se)
  expect_s3_class(eds, "ExploreDatasetSE_Proteomics")
  expect_identical(eds$omics_class(), "proteomics")
  expect_identical(eds$data_type(), "proteomics")
  expect_gt(length(eds$contrasts()), 0L)
})

test_that("SAINT SE: log2_EFCs/BFDR map to log2Ratio/fdr, pValue all-NA", {
  eds <- new_explore_dataset(fx_proteomics_saint()$se)
  d <- eds$de(eds$primary_contrast())
  expect_true(all(c("log2Ratio", "fdr", "pValue") %in% colnames(d)))
  expect_true(any(!is.na(d$log2Ratio)))
  expect_true(any(!is.na(d$fdr)))
  expect_true(all(is.na(d$pValue)))
  expect_false(eds$has_pvalue())
})

test_that("SAINT SE: usedInTest falls back to fdr, not emptied to all-FALSE", {
  eds <- new_explore_dataset(fx_proteomics_saint()$se)
  d <- eds$de(eds$primary_contrast())
  expect_true("usedInTest" %in% colnames(d))
  expect_true(any(d$usedInTest))
  expect_equal(d$usedInTest, !is.na(d$fdr))
})

test_that("SAINT SE: SaintScore surfaces as a display-only extra", {
  eds <- new_explore_dataset(fx_proteomics_saint()$se)
  expect_identical(eds$de_display_extras(), "SaintScore")
  d <- eds$de(eds$primary_contrast())
  expect_true("SaintScore" %in% colnames(d))
  expect_false(isTRUE(all.equal(d$SaintScore, d$fdr)))
})

test_that("SAINT SE: reference group resolves to the CONTROL level", {
  eds <- new_explore_dataset(fx_proteomics_saint()$se)
  cm <- eds$contrast_meta(eds$primary_contrast())
  expect_identical(cm$method, "SAINTexpress")
  expect_identical(cm$reference_group, "eGFP")
  expect_identical(cm$score_label, "SaintScore")
})

test_that("SAINT SE: standard prolfqua proteomics is unaffected (no SAINT tagging)", {
  eds <- new_explore_dataset(fx_proteomics_protein()$se)
  expect_identical(eds$de_display_extras(), character())
  expect_true(eds$has_pvalue())
  cm <- eds$contrast_meta(eds$primary_contrast())
  expect_identical(cm$method, "prolfqua")
  expect_null(cm$score_label)
})

test_that("a prolfqua SummarizedExperiment without contrasts opens as proteomics with all its features", {
  se <- readRDS(system.file("examples", "proteomics", "SummarizedExperiment.rds", package = "exploreDE"))
  SummarizedExperiment::rowData(se)$constrast_B_vs_A <- NULL
  expect_identical(detect_se_shape(se), "proteomics")
  eds <- new_explore_dataset(se)
  expect_identical(eds$omics_class(), "proteomics")
  expect_length(eds$contrasts(), 0)
  expect_identical(eds$features()$feature_id, rownames(se))
  expect_true(validate_explore_dataset(eds, strict = FALSE)$ok)
})

test_that("the raw assay of proteomics data is the raw intensities when there are any", {
  se <- readRDS(system.file("examples", "proteomics", "SummarizedExperiment.rds", package = "exploreDE"))
  # rawData is not offered: the transformed intensities are the only assay
  expect_identical(new_explore_dataset(se)$raw_assay_name(), "intensities")
  SummarizedExperiment::assayNames(se)[SummarizedExperiment::assayNames(se) == "rawData"] <- "raw"
  eds <- new_explore_dataset(se)
  expect_identical(eds$raw_assay_name(), "intensities_raw")
  expect_identical(eds$primary_assay_name(), "intensities")

  skip_if_not_installed("anndataR")
  pr <- make_toy_prot("protein")
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(
    list(raw = 2^pr$intensities, transformed = pr$intensities),
    out = out,
    omics = "proteomics",
    coldata = pr$coldata,
    rowdata = pr$rowdata
  )
  eds <- new_explore_dataset(out)
  expect_identical(eds$raw_assay_name(), "intensities_raw")
  expect_identical(eds$primary_assay_name(), "intensities")
  expect_identical(eds$raw_assay(), eds$assay("intensities_raw"))
})

test_that("eds_proteomics_choices maps gene symbols and protein ids to feature ids", {
  ch <- eds_proteomics_choices(
    feature_id = c("P1~a", "P1~b", "P2~c", "P3"),
    gene = c("G1", "G1", NA, "P3"),
    protein = c("P1", "P1", "P2", "P3")
  )
  expect_identical(ch$gene, list(choices = "G1", lookup = list(G1 = c("P1~a", "P1~b"))))
  expect_identical(ch$protein, list(choices = c("P1", "P2"), lookup = list(P1 = c("P1~a", "P1~b"), P2 = "P2~c")))
  # Protein-level data: the protein ids are the feature ids
  none <- list(choices = NULL, lookup = NULL)
  expect_identical(eds_proteomics_choices(c("P1", "P2"), protein = c("P1", "P2"))$protein, none)
  expect_identical(eds_proteomics_choices(c("P1", "P2"))$gene, none)
})
