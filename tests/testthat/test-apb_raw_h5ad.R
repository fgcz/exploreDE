# Raw APB (apb2) results: detection, and conversion to exploreDE/1.0 through
# apbCatalogR. The conversion tests use the synthetic results apbCatalogR ships
# in its test bundle (Parquet, so they need arrow); apbCatalogR reads them into
# the same ParsedLevels as an .h5ad.

apb_results_dir <- function() {
  skip_if_not_installed("apbCatalogR")
  skip_if_not_installed("arrow")
  archive <- system.file("extdata", "test-bundle", "results.tar.gz", package = "apbCatalogR")
  skip_if(!nzchar(archive), "apbCatalogR ships no test results")
  dir <- file.path(tempdir(), "apb_test_results")
  if (!dir.exists(dir)) {
    utils::untar(archive, exdir = dir)
  }
  dir
}

read_apb <- function(name) {
  apbCatalogR::read_parsed_levels(file.path(apb_results_dir(), paste0(name, ".parquet")))
}

# A minimal .h5ad whose uns looks like a raw APB result's (root record, no
# schema_version). Only for detection: Its contents aren't an apb2 layout.
write_fake_apb_h5ad <- function(path, uns = list(apb = list(parse = list(note = "root record")))) {
  ad <- anndataR::AnnData(
    X = matrix(1:4, 2, 2, dimnames = list(c("S0", "S1"), c("F0", "F1"))),
    obs = data.frame(row.names = c("S0", "S1")),
    var = data.frame(row.names = c("F0", "F1")),
    uns = uns
  )
  anndataR::write_h5ad(ad, path, mode = "w")
  path
}

test_that(".apb_short_name takes the first input with a value, its first member, and NAME from db|ACC|NAME", {
  expect_equal(
    .apb_short_name(list(c("GENE1;GENE2", "", NA), c("x", "sp|P12345|ENO_HUMAN;sp|Q1|X_HUMAN", "P99999"))),
    c("GENE1", "ENO_HUMAN", "P99999")
  )
  expect_equal(.apb_short_name(list(c("A", NA)), append = c("PEPTIDE", "IONS/2")), c("A_PEPTIDE", NA))
  expect_equal(.apb_short_name(list(NULL, c("P1", "P2"))), c("P1", "P2"))
})

test_that(".apb_obs_annotations renames SDRF and prolfquapp factors, sorted, and passes other columns through", {
  obs <- data.frame(
    Run = c("r1", "r2", "r3"),
    factor_value_Treatment = c("trt", "ctrl", "trt"),
    `Group Name` = c("b", "a", "b"),
    characteristics_biological_replicate = c(1, 2, 3),
    check.names = FALSE
  )
  out <- .apb_obs_annotations(obs, keys = "Run")
  expect_named(out, c("factor_treatment", "factor_group_name", "characteristics_biological_replicate"))
  expect_equal(levels(out$factor_treatment), c("ctrl", "trt"))
  expect_equal(as.character(out$factor_group_name), c("b", "a", "b"))
  expect_equal(out$characteristics_biological_replicate, c(1, 2, 3))
})

test_that(".apb_contaminants: the vendor marking, or the FASTA check's contaminant match", {
  marked <- c(TRUE, FALSE, FALSE, NA)
  no_fasta <- .apb_contaminants(list(varm = list()), marked)
  expect_equal(no_fasta$value, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(no_fasta$sources, "apb contaminant")

  with_fasta <- list(varm = list(fasta_validation = data.frame(fasta_matches_contaminant = c(FALSE, TRUE, NA, FALSE))))
  out <- .apb_contaminants(with_fasta, marked)
  expect_equal(out$value, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(out$sources, c("apb contaminant", "FASTA fasta_matches_contaminant"))

  # apb-export reads the match from fasta_validation only.
  other_table <- list(varm = list(fasta = data.frame(fasta_matches_contaminant = c(TRUE, TRUE, TRUE, TRUE))))
  expect_equal(.apb_contaminants(other_table, marked)$value, c(TRUE, FALSE, FALSE, FALSE))
})

test_that("raw APB detection: uns['apb'] without schema_version or prolfquapp", {
  expect_true(.h5ad_is_apb_raw(c("apb", "protein")))
  expect_false(.h5ad_is_apb_raw(c("apb", "schema_version", "omics_type"))) # APB Export's explorede output
  expect_false(.h5ad_is_apb_raw(c("apb", "prolfquapp")))
  expect_false(.h5ad_is_apb_raw(c("schema_version", "de")))
  expect_false(.h5ad_is_apb_raw(NULL))

  path <- write_fake_apb_h5ad(tempfile(fileext = ".h5ad"))
  on.exit(unlink(path))
  expect_true("apb" %in% .h5ad_uns_names(path))
  expect_null(.h5ad_uns_names(tempfile(fileext = ".h5ad"))) # no such file
})

test_that("protein level: log2 intensities with zeros missing, decoys dropped, contaminants flagged", {
  r <- read_apb("diann_v2-protein")
  lvl <- r$levels$protein
  decoy <- as.logical(lvl$var$apb_Decoy)
  ad <- anndata_from_apb_raw(r)

  expect_equal(ad$uns$schema_version, "exploreDE/1.0")
  expect_equal(ad$uns$feature_id_kind, "protein")
  expect_equal(ad$uns$provenance$primary_layer, "PG_MaxLFQ")
  expect_equal(ad$uns$provenance$decoys_dropped, sum(decoy))
  expect_equal(ad$var_names, lvl$var$Protein_Group[!decoy])
  expect_equal(ad$var$label_protein_id, lvl$var$Protein_Group[!decoy])
  expect_equal(ad$var$is_contaminant, as.logical(lvl$var$apb_Contaminant[!decoy]))
  expect_equal(ad$uns$provenance$contaminant_sources, "apb contaminant")

  x <- as.matrix(lvl$X)[, !decoy, drop = FALSE]
  raw <- ad$layers$intensities_raw
  expect_equal(unname(raw[x != 0]), unname(x[x != 0]))
  expect_true(all(is.na(raw[x == 0])))
  expect_equal(unname(ad$layers$intensities[x > 0]), log2(unname(x[x > 0])))

  eds <- ExploreDatasetH5AD$new(ad)
  expect_equal(eds$data_type(), "proteomics")
  expect_setequal(eds$assay_names(), c("intensities", "intensities_raw"))
  expect_true(isTRUE(validate_explore_dataset(eds, strict = FALSE)$ok))
})

test_that("drop_decoys = FALSE keeps decoys and marks them", {
  r <- read_apb("diann_v2-protein")
  ad <- anndata_from_apb_raw(r, drop_decoys = FALSE)
  expect_equal(ad$n_vars(), r$levels$protein$n_vars())
  expect_equal(ad$var$is_decoy, as.logical(r$levels$protein$var$apb_Decoy))
  expect_equal(ad$uns$provenance$decoys_dropped, 0L)
})

test_that("precursor level: feature_id is <protein>_<ProForma ion>, the label carries the ion", {
  r <- read_apb("diann_v2-peptidoform_ion")
  lvl <- r$levels$peptidoform_ion
  keep <- !as.logical(lvl$var$apb_Decoy)
  ad <- anndata_from_apb_raw(r)

  expect_equal(ad$uns$feature_id_kind, "precursor")
  expect_equal(ad$uns$hierarchy$unit, "precursor")
  expect_equal(ad$var_names, paste(lvl$var$Protein_Group[keep], lvl$var$ProForma_ion[keep], sep = "_"))
  expect_equal(ad$var$label_gene_symbol, paste(lvl$var$Protein_Group[keep], lvl$var$ProForma_ion[keep], sep = "_"))
  expect_equal(ExploreDatasetH5AD$new(ad)$feature_id_kind(), "precursor")
})

test_that("peptide level becomes feature_id_kind peptide (peptidomics)", {
  ad <- anndata_from_apb_raw(read_apb("maxquant_wide-peptide"))
  expect_equal(ad$uns$feature_id_kind, "peptide")
  expect_equal(ExploreDatasetH5AD$new(ad)$data_type(), "peptidomics")
})

test_that("SDRF factor values become sorted factor_* columns; the run key gives the sample ids", {
  r <- read_apb("edge-sdrf")
  ad <- anndata_from_apb_raw(r)
  obs <- ad$obs
  expect_equal(ad$obs_names, as.character(r$levels[[1]]$obs$Run))
  expect_equal(obs$sample_id, ad$obs_names)
  expect_true(all(c("factor_treatment", "factor_dose", "characteristics_biological_replicate") %in% colnames(obs)))
  expect_false("Run" %in% colnames(obs))
  expect_equal(levels(obs$factor_treatment), sort(unique(as.character(r$levels[[1]]$obs$factor_value_treatment))))
  expect_equal(ExploreDatasetH5AD$new(ad)$factor_names()[[1]], "treatment")
})

test_that("an unreviewed rule falls back to the file's own records, and says so", {
  ad <- anndata_from_apb_raw(read_apb("edge-unreviewed"))
  fb <- ad$uns$provenance$catalogue_fallbacks
  expect_true(any(grepl("^protein_assignment from column Protein_IDs", fb)))
  expect_true(any(grepl("^decoy from column apb_Decoy", fb)))
  expect_true(isTRUE(validate_explore_dataset(ExploreDatasetH5AD$new(ad), strict = FALSE)$ok))
})

test_that("level choice: protein first; another level on request; unknown levels refused", {
  r <- read_apb("edge-two-levels")
  expect_equal(anndata_from_apb_raw(r)$uns$feature_id_kind, "protein")
  expect_equal(anndata_from_apb_raw(r, level = "peptidoform_ion")$uns$feature_id_kind, "precursor")
  expect_error(anndata_from_apb_raw(r, level = "peptide"), "has no level 'peptide'")
  expect_error(anndata_from_apb_raw(read_apb("maxquant_wide-peptidoform")), "none of the levels exploreDE reads")
})

test_that("new_explore_dataset() routes a raw APB .h5ad through apbCatalogR", {
  r <- read_apb("diann_v2-protein")
  path <- write_fake_apb_h5ad(tempfile(fileext = ".h5ad"))
  on.exit(unlink(path))
  local_mocked_bindings(.apb_read = function(path) r)
  eds <- new_explore_dataset(path)
  expect_s3_class(eds, "ExploreDatasetH5AD")
  expect_equal(eds$feature_id_kind(), "protein")
})

test_that("a raw APB .h5ad without apbCatalogR gives a clear message", {
  path <- write_fake_apb_h5ad(tempfile(fileext = ".h5ad"))
  on.exit(unlink(path))
  local_mocked_bindings(.apb_available = function() FALSE)
  expect_error(new_explore_dataset(path), "needs the apbCatalogR package")
})

test_that("a live AnnData that is a raw APB result is refused with a hint", {
  ad <- anndataR::AnnData(
    X = matrix(1:4, 2, 2, dimnames = list(c("S0", "S1"), c("F0", "F1"))),
    obs = data.frame(row.names = c("S0", "S1")),
    var = data.frame(row.names = c("F0", "F1")),
    uns = list(apb = list(parse = list(note = "root record")))
  )
  expect_error(new_explore_dataset(ad), "Pass the file's path")
})
