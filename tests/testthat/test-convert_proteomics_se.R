# Shared assertions helper — structural checks that apply to every proteomics
# conversion regardless of feature level.
expect_proteomics_shape <- function(out) {
  expect_type(out, "list")
  expect_equal(out$dataType, "proteomics")
  expect_true(all(c("se", "seqAnnoList", "contrasts", "dataset", "factors",
                    "countList", "genes", "factorLevels", "groupingName",
                    "gene_lookup", "protein_lookup",
                    "gene_choices", "protein_choices") %in% names(out)))

  # --- Contrasts --------------------------------------------------------------
  expect_type(out$contrasts, "character")
  expect_gt(length(out$contrasts), 0)
  expect_setequal(names(out$seqAnnoList), out$contrasts)

  # --- seqAnnoList: each entry is a sorted data.frame with feature_id --------
  for (con in out$contrasts) {
    sa <- out$seqAnnoList[[con]]
    expect_s3_class(sa, "data.frame")
    expect_true("feature_id" %in% colnames(sa))
    expect_equal(anyDuplicated(sa$feature_id), 0L)
    # Sorted by raw pValue (ascending), NAs allowed
    expect_false(is.unsorted(sa$pValue, na.rm = TRUE))
  }

  # --- countList: each matrix has the right shape ----------------------------
  expect_type(out$countList, "list")
  expect_gt(length(out$countList), 0)
  for (m in names(out$countList)) {
    mat <- out$countList[[m]]
    # Rownames must be a subset of feature_ids in at least one contrast
    all_fids <- unique(unlist(lapply(out$seqAnnoList, `[[`, "feature_id")))
    expect_true(all(rownames(mat) %in% all_fids),
      info = paste("count matrix", m, "has rownames not in any feature_id"))
    # Columns match sample dataset
    expect_equal(ncol(mat), nrow(out$dataset),
      info = paste("count matrix", m, "column count != dataset rows"))
  }

  # --- genes vector ----------------------------------------------------------
  expect_type(out$genes, "character")
  expect_equal(anyDuplicated(out$genes), 0L)

  # --- Dataset / factors -----------------------------------------------------
  expect_s3_class(out$dataset, "data.frame")
  expect_true("Name" %in% colnames(out$dataset))
  expect_true(out$groupingName %in% out$factors)
}


# ============================================================================
# Protein-level: feature_id == protein_Id, so protein_Id dropped as redundant
# ============================================================================
test_that("convert_proteomics_se on protein-level fixture has clean protein-level structure", {
  fx <- fx_proteomics_protein()
  out <- convert_proteomics_se(fx$se)

  expect_proteomics_shape(out)

  sa1 <- out$seqAnnoList[[1]]

  # protein_Id is dropped when it equals feature_id (protein-level data)
  expect_false("protein_Id" %in% colnames(sa1))
  # Correspondingly: protein_lookup and protein_choices are NULL
  expect_null(out$protein_lookup)
  expect_null(out$protein_choices)

  # gene_lookup is populated whenever gene_name differs from feature_id
  if ("gene_name" %in% colnames(sa1) && !identical(sa1$feature_id, sa1$gene_name)) {
    expect_type(out$gene_lookup, "list")
    expect_type(out$gene_choices, "character")
    expect_gt(length(out$gene_choices), 0)
    # Every looked-up name maps to ≥1 feature_id
    for (gene in head(names(out$gene_lookup), 5)) {
      expect_true(all(out$gene_lookup[[gene]] %in% sa1$feature_id))
    }
  }
})


# ============================================================================
# Peptide-level: feature_id is protein_Id + peptide suffix; protein_Id retained
# ============================================================================
test_that("convert_proteomics_se on peptide-level fixture constructs composite feature_ids", {
  fx <- fx_proteomics_peptide()
  out <- convert_proteomics_se(fx$se)

  expect_proteomics_shape(out)

  sa1 <- out$seqAnnoList[[1]]

  # feature_id != protein_Id (peptide-level), so protein_Id column is preserved
  expect_true("protein_Id" %in% colnames(sa1))
  # ... and protein_lookup is populated
  expect_type(out$protein_lookup, "list")
  expect_type(out$protein_choices, "character")
  expect_gt(length(out$protein_choices), 0)

  # Each feature_id starts with a known protein_Id
  sample_fids <- head(sa1$feature_id, 20)
  sample_pids <- head(sa1$protein_Id, 20)
  expect_true(all(startsWith(sample_fids, sample_pids)))

  # protein_lookup[pid] returns the feature_ids for that protein
  sample_protein <- out$protein_choices[[1]]
  expect_true(all(out$protein_lookup[[sample_protein]] %in% sa1$feature_id))
})


# ============================================================================
# Phospho-level: feature_id is protein_Id + site
# ============================================================================
test_that("convert_proteomics_se on phospho-level fixture builds protein_Id_site feature_ids", {
  fx <- fx_proteomics_phospho()
  out <- convert_proteomics_se(fx$se)

  expect_proteomics_shape(out)

  sa1 <- out$seqAnnoList[[1]]

  # feature_id composite keeps protein_Id prefix
  expect_true("protein_Id" %in% colnames(sa1))
  # And site column is present (phospho-specific)
  expect_true("site" %in% colnames(sa1))

  # protein_lookup populated
  expect_type(out$protein_lookup, "list")
  expect_gt(length(out$protein_choices), 0)

  # Every feature_id starts with a protein_Id
  sample_fids <- head(sa1$feature_id, 20)
  sample_pids <- head(sa1$protein_Id, 20)
  expect_true(all(startsWith(sample_fids, sample_pids)))
})


# ============================================================================
# Cross-cutting: countList rownames align with feature_ids (same fix we shipped
# to guard against the sort-order alignment bug in the lookup).
# ============================================================================
test_that("convert_proteomics_se countList rownames match unique feature_ids from seqAnnoList", {
  for (key in c("proteomics_protein", "proteomics_peptide", "proteomics_phospho")) {
    fx <- load_fixture(key)
    out <- convert_proteomics_se(fx$se)

    all_fids <- unique(unlist(lapply(out$seqAnnoList, `[[`, "feature_id")))
    for (m in names(out$countList)) {
      rn <- rownames(out$countList[[m]])
      expect_equal(anyDuplicated(rn), 0L,
        info = paste(key, "/", m, ": duplicate feature_ids in count matrix"))
      expect_true(all(rn %in% all_fids),
        info = paste(key, "/", m, ": rownames not in seqAnnoList feature_ids"))
    }
  }
})
