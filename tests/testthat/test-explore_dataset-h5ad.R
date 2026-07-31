# -----------------------------------------------------------------------------
# ExploreDatasetH5AD + build_explore_h5ad (source-agnostic constructor) tests
#
# Covers the finalised contract (plans/byo_data_import_plan.md §0c):
#   - rich object (DESeqDataSet) route
#   - loose-pieces route (matrix + colData + omics, no rich object)
#   - counts as a named collection of matrices
#   - DE with non-canonical columns via explicit de_columns mapping
#   - rowData pathway annotation (A) -> tile/heatmap capability
#   - pathway statistics (B) keyed by contrast
#   - cross-API equivalence vs ExploreDatasetSE_RNASeq (real data)
# -----------------------------------------------------------------------------

skip_if_no_anndataR <- function() testthat::skip_if_not_installed("anndataR")

# ---- 1. rich object (DESeqDataSet) route ------------------------------------

test_that("rich object: DESeqDataSet round-trips into a valid eds", {
  skip_if_no_anndataR()
  dds <- make_toy_dds()
  res <- DESeq2::results(dds, contrast = c("condition", "treated", "control"))
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(dds, out = out, omics = "rnaseq",
                     de = list(treated_vs_control = res), organism = "Homo sapiens")

  eds <- new_explore_dataset(out)
  expect_s3_class(eds, "ExploreDatasetH5AD")
  expect_identical(eds$data_type(), "RNASeq")
  expect_silent(validate_explore_dataset(eds))
  expect_identical(eds$factor_names(), "condition")      # not sizeFactor
  expect_true("treated_vs_control" %in% eds$contrasts())
  expect_true(eds$has_assay("normalised"))

  de <- eds$de("treated_vs_control")
  ref <- setNames(res$log2FoldChange, rownames(dds))
  common <- intersect(rownames(de), names(ref))
  expect_lt(max(abs(de[common, "effect_log2fc"] - ref[common]), na.rm = TRUE), 1e-6)
})

test_that("omics is required", {
  skip_if_no_anndataR()
  dds <- make_toy_dds()
  out <- withr::local_tempfile(fileext = ".h5ad")
  expect_error(build_explore_h5ad(dds, out = out), "omics")
})

# ---- 2. loose-pieces route (no rich object) ---------------------------------

test_that("loose pieces: counts + colData + omics, no DE, builds a valid eds", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata)

  eds <- new_explore_dataset(out)
  expect_silent(validate_explore_dataset(eds))
  expect_identical(eds$data_type(), "RNASeq")
  expect_false(eds$has_de())                             # counts-only degrades gracefully
  expect_identical(eds$factor_names(), "condition")
  expect_identical(dim(eds$assay("counts")), dim(p$counts))
})

test_that("counts as a named collection of matrices exposes each assay", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  tpm <- p$normalised * 1.5                              # stand-in tpm matrix
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(list(counts = p$counts, normalised = p$normalised, tpm = tpm),
                     out = out, omics = "rnaseq", coldata = p$coldata)

  eds <- new_explore_dataset(out)
  an <- eds$assay_names()
  expect_true(all(c("counts", "normalised", "normalised_log2", "vst", "tpm") %in% an))
  a <- eds$assay("tpm")                                 # features x samples
  expect_lt(max(abs(a - tpm[rownames(a), colnames(a)]), na.rm = TRUE), 1e-6)
})

# ---- 3. DE with non-canonical columns via explicit mapping ------------------

test_that("de_columns maps non-standard DE column names", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()                                 # de_df has LFC / P / Q
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata,
                     de = list(treated_vs_control = p$de_df),
                     de_columns = list(effect_log2fc = "LFC", score_pvalue = "P", score_fdr = "Q"))

  eds <- new_explore_dataset(out)
  expect_true(eds$has_de())
  de <- eds$de("treated_vs_control")
  ref <- setNames(p$de_df$LFC, rownames(p$de_df))
  common <- intersect(rownames(de), names(ref))
  expect_gt(length(common), 0L)
  expect_lt(max(abs(de[common, "effect_log2fc"] - ref[common]), na.rm = TRUE), 1e-6)
})

# ---- 4. pathway ANNOTATION (A) in rowData -----------------------------------

test_that("rowData pathway annotation surfaces as GO columns (tile/heatmap capability)", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata,
                     rowdata = p$rowdata)

  eds <- new_explore_dataset(out)
  expect_true("GO BP" %in% colnames(eds$features()))
  expect_true(any(!is.na(eds$features()[["GO BP"]])))
})

# ---- 5. pathway STATISTICS (B), contrast-keyed ------------------------------

test_that("pathway stats (data.frame) are stored under their contrast", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  ora_df <- data.frame(ID = c("GO:1", "GO:2"), Description = c("a", "b"),
                       pvalue = c(0.01, 0.2), stringsAsFactors = FALSE)
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata,
                     de = list(treated_vs_control = p$de_df),
                     de_columns = list(effect_log2fc = "LFC", score_pvalue = "P", score_fdr = "Q"),
                     pathways = list(treated_vs_control = list(my_ora = ora_df)))

  eds <- new_explore_dataset(out)
  enr <- eds$backend()$uns$enrichment
  expect_true("treated_vs_control" %in% names(enr) || "treated_vs_control" %in% names(eds$backend()$uns$enrichment))
})

test_that("nested SUSHI-shaped pathway stats round-trip and flatten", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  # Lightweight cP stand-ins: S3 objects that inherit the class the flattener
  # keys on. Serialise/deserialise like the real S4 objects for this purpose.
  er <- function(n) structure(list(result = data.frame(ID = paste0("GO:", seq_len(n)))),
                              class = "enrichResult")
  gr <- function(n) structure(list(result = data.frame(ID = paste0("GO:", seq_len(n)))),
                              class = "gseaResult")
  # SUSHI shapes: ORA nested ontology x direction; GSEA nested by ontology.
  ora  <- list(BP = list(upGenes = er(3), downGenes = er(2), bothGenes = er(5)),
               MF = list(upGenes = er(1), downGenes = er(1), bothGenes = er(2)),
               CC = list(upGenes = er(4), downGenes = er(1), bothGenes = er(6)))
  gsea <- list(BP = gr(4), MF = gr(2), CC = gr(6))
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata,
                     de = list(treated_vs_control = p$de_df),
                     de_columns = list(effect_log2fc = "LFC", score_pvalue = "P", score_fdr = "Q"),
                     pathways = list(treated_vs_control = list(ora = ora, gsea = gsea)))

  eds <- new_explore_dataset(out)
  expect_setequal(eds$enrichment_kinds(), c("ora", "gsea"))

  ora_flat <- eds$enrichment("treated_vs_control", "ora")
  expect_length(ora_flat, 9L)                              # 3 ontologies x 3 directions
  expect_true("ora/BP/upGenes" %in% names(ora_flat))
  expect_s3_class(ora_flat[["ora/BP/upGenes"]], "enrichResult")

  gsea_flat <- eds$enrichment("treated_vs_control", "gsea")
  expect_length(gsea_flat, 3L)                             # BP / MF / CC
  expect_true("gsea/CC" %in% names(gsea_flat))
  expect_s3_class(gsea_flat[["gsea/BP"]], "gseaResult")
})

test_that("pathway stats for a non-existent contrast warn", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  out <- withr::local_tempfile(fileext = ".h5ad")
  expect_warning(
    build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata,
                       de = list(treated_vs_control = p$de_df),
                       de_columns = list(effect_log2fc = "LFC", score_pvalue = "P", score_fdr = "Q"),
                       pathways = list(nonexistent_contrast = list(x = data.frame(ID = "GO:1")))),
    "non-existent contrast")
})

# ---- 5b. Multi-DE per-contrast fan-out (legacy_multi_contrast_from_eds) ------

test_that("legacy_multi_contrast_from_eds fans a multi-contrast h5ad out per contrast", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  # Two contrasts against the same reference so both test/reference groups are
  # real levels of the grouping factor.
  de2 <- p$de_df; de2$LFC <- -de2$LFC
  er <- function(n) structure(
    list(result = data.frame(ID = paste0("GO:", seq_len(n)),
                             Description = paste("term", seq_len(n)))),
    class = "enrichResult")
  gr <- function(n) structure(
    list(result = data.frame(ID = paste0("GO:", seq_len(n)),
                             Description = paste("term", seq_len(n)))),
    class = "gseaResult")
  ora  <- list(BP = list(upGenes = er(3), downGenes = er(2), bothGenes = er(5)),
               MF = list(upGenes = er(1), downGenes = er(1), bothGenes = er(2)))
  gsea <- list(BP = gr(4), MF = gr(2))
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(
    p$counts, out = out, omics = "rnaseq", coldata = p$coldata,
    de = list(treated_vs_control = p$de_df, control_vs_treated = de2),
    de_columns = list(effect_log2fc = "LFC", score_pvalue = "P", score_fdr = "Q"),
    pathways = list(treated_vs_control = list(ora = ora, gsea = gsea)))

  eds <- new_explore_dataset(out)
  expect_setequal(eds$contrasts(), c("treated_vs_control", "control_vs_treated"))

  o <- legacy_multi_contrast_from_eds(eds, "treated_vs_control")

  # Legacy per-contrast gStoreList shape the Multi-DE modules consume.
  expect_identical(o$dataType, "RNASeq")
  expect_identical(o$design, "treated--over--control")
  expect_identical(o$designLevels, c("treated", "control"))
  expect_identical(o$param$sampleGroup, "treated")
  expect_identical(o$param$refGroup, "control")
  expect_true(o$groupingName %in% colnames(o$dataset))
  expect_true("Name" %in% colnames(o$dataset))

  # seqAnno carries feature_id + DE stats, sorted by p-value, aligned to genes.
  expect_true(all(c("feature_id", "log2Ratio", "pValue", "fdr") %in% colnames(o$seqAnno)))
  expect_false(is.unsorted(o$seqAnno$pValue, na.rm = TRUE))
  expect_identical(o$genes, o$seqAnno$feature_id)

  # countList uses the legacy display names and feature_id rownames.
  expect_true("Raw" %in% names(o$countList))
  expect_setequal(rownames(o$countList[["Raw"]]), o$seqAnno$feature_id)

  # Enrichment reshaped back to nested ont -> sig (ORA) / ont (GSEA). (The
  # `Label` column added by .label_enrichment_results is exercised against real
  # S4 clusterProfiler objects in the manual round-trip check — the S3
  # stand-ins here collide with clusterProfiler's own S4 `[[`/`@` methods.)
  expect_setequal(names(o$ora), c("BP", "MF"))
  expect_setequal(names(o$ora$BP), c("upGenes", "downGenes", "bothGenes"))
  expect_s3_class(o$ora$BP$upGenes, "enrichResult")
  expect_setequal(names(o$gsea), c("BP", "MF"))
  expect_s3_class(o$gsea$BP, "gseaResult")

  # The second contrast has no enrichment attached (pathways only for the first).
  o2 <- legacy_multi_contrast_from_eds(eds, "control_vs_treated")
  expect_null(o2$ora)
  expect_null(o2$gsea)
  expect_identical(o2$param$sampleGroup, "control")
})

# ---- 6. feature_level: detection-with-override ------------------------------

test_that("rnaseq isoform is detected from a transcript_id column", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  rd <- p$rowdata; rd$transcript_id <- paste0("ENST", seq_len(nrow(rd)))
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata, rowdata = rd)
  eds <- new_explore_dataset(out)
  expect_identical(eds$feature_id_kind(), "isoform")
  expect_identical(eds$data_type(), "RNASeq")
})

test_that("proteomics flavour: protein & phospho detected, peptide via override", {
  skip_if_no_anndataR()
  build_prot <- function(pr, ...) {
    out <- withr::local_tempfile(fileext = ".h5ad", .local_envir = parent.frame())
    build_explore_h5ad(pr$intensities, out = out, omics = "proteomics",
                       coldata = pr$coldata, rowdata = pr$rowdata, ...)
    new_explore_dataset(out)
  }

  e_prot <- build_prot(make_toy_prot("protein"))
  expect_identical(e_prot$omics_class(), "proteomics")
  expect_identical(e_prot$feature_id_kind(), "protein")
  expect_identical(e_prot$data_type(), "proteomics")
  expect_silent(validate_explore_dataset(e_prot))

  e_ph <- build_prot(make_toy_prot("phospho"))
  expect_identical(e_ph$feature_id_kind(), "phospho")
  expect_identical(e_ph$data_type(), "phosphoproteomics")
  expect_true(all(grepl("_S[0-9]+$", e_ph$feature_ids())))   # protein_Id_site, unique

  e_pep <- build_prot(make_toy_prot("peptide"), feature_level = "peptide")
  expect_identical(e_pep$feature_id_kind(), "peptide")
  expect_identical(e_pep$data_type(), "peptidomics")

  # Intensity layers surface as intensities (not RNA-seq "counts") and no
  # RNA-seq-only derived assays (vst/tpm/fpkm) are synthesised.
  expect_true("intensities" %in% e_prot$assay_names())
  expect_false(any(c("counts", "vst", "tpm", "fpkm") %in% e_prot$assay_names()))
  expect_identical(e_prot$primary_assay_name(), "intensities")

  # proteomics_choices: gene lookup populated (from gene_name); protein choices
  # populate only when the protein id differs from the feature id — i.e. at
  # peptide/phospho level, NOT at protein level (matches the SE backend).
  gp <- e_prot$proteomics_choices("gene")
  expect_true(length(gp$choices) > 0)
  expect_true(all(unlist(gp$lookup) %in% e_prot$feature_ids()))
  # At peptide level feature_id = protein_Id_peptide, so protein choices populate
  # and each resolves to real features.
  pp <- e_pep$proteomics_choices("protein")
  expect_true(length(pp$choices) > 0)
  expect_true(all(unlist(pp$lookup) %in% e_pep$feature_ids()))

  # Grouping falls back to a real sample-annotation column (no SUSHI "G_").
  expect_true(e_prot$primary_grouping() %in% colnames(e_prot$samples()))
})

test_that("an invalid feature_level for the omics errors", {
  skip_if_no_anndataR()
  pr <- make_toy_prot("protein")
  out <- withr::local_tempfile(fileext = ".h5ad")
  expect_error(
    build_explore_h5ad(pr$intensities, out = out, omics = "proteomics",
                       coldata = pr$coldata, rowdata = pr$rowdata, feature_level = "isoform"),
    "invalid")
})

# ---- 6b. S4 DataFrame colData with a "[Factor]" column doesn't mangle -------

test_that("S4 DataFrame colData preserves '[Factor]' column names and factors cleanly", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  # colData as an S4 DataFrame with a SUSHI-style "[Factor]" tag + a junk col.
  cd <- S4Vectors::DataFrame(
    "Condition [Factor]" = rep(c("ctrl", "trt"), each = 3),
    "Count [File]"       = paste0("f", 1:6, ".txt"),
    row.names = colnames(p$counts), check.names = FALSE)
  out <- withr::local_tempfile(fileext = ".h5ad")

  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = cd,
                     factors = "Condition [Factor]")   # name must survive to match

  eds <- new_explore_dataset(out)
  expect_silent(validate_explore_dataset(eds))
  expect_identical(eds$factor_names(), "condition")    # [Factor] stripped, not mangled
  expect_false(is.na(eds$primary_grouping()))
  expect_setequal(as.character(eds$samples()$condition), c("ctrl", "trt"))
})

# ---- 7. cross-API equivalence vs the SE backend (real data) -----------------

test_that("ExploreDatasetH5AD matches ExploreDatasetSE_RNASeq accessor-for-accessor", {
  skip_if_no_anndataR()
  fx  <- fx_rnaseq_gene()                                # skips if gstore unreachable
  eds_se <- new_explore_dataset(fx$se)

  out <- withr::local_tempfile(fileext = ".h5ad")
  se_to_canonical_h5ad(fx$se, out)
  eds_h5 <- new_explore_dataset(out)

  expect_silent(validate_explore_dataset(eds_h5))
  expect_identical(eds_se$data_type(), eds_h5$data_type())
  expect_identical(eds_se$n_features(), eds_h5$n_features())
  expect_setequal(eds_se$feature_ids(), eds_h5$feature_ids())
  expect_setequal(tolower(eds_se$factor_names()), tolower(eds_h5$factor_names()))
  expect_setequal(eds_se$contrasts(), eds_h5$contrasts())

  c1 <- eds_se$primary_contrast()
  dse <- eds_se$de(c1); dh5 <- eds_h5$de(c1)
  common <- intersect(rownames(dse), rownames(dh5))
  expect_gt(length(common), 0L)
  for (col in c("effect_log2fc", "score_pvalue", "score_fdr"))
    expect_lt(max(abs(dse[common, col] - dh5[common, col]), na.rm = TRUE), 1e-6)

  # incl. fpkm/tpm — h5ad computes from featWidth, SE from ezRun::getRpkm/getTpm.
  for (nm in c("counts", "normalised", "normalised_log2", "vst", "fpkm", "tpm")) {
    expect_true(eds_se$has_assay(nm) && eds_h5$has_assay(nm))
    a <- eds_se$assay(nm); b <- eds_h5$assay(nm)
    fs <- intersect(rownames(a), rownames(b)); ss <- intersect(colnames(a), colnames(b))
    expect_lt(max(abs(a[fs, ss] - b[fs, ss]), na.rm = TRUE), 1e-6)
  }
})

test_that("computed TPM/FPKM from a generic (non-SUSHI) gene-length column", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  rd <- p$rowdata
  rd$gene_length <- sample(500:5000, nrow(rd), replace = TRUE)   # NOT 'featWidth'
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata, rowdata = rd)

  eds <- new_explore_dataset(out)
  expect_true(all(c("tpm", "fpkm") %in% eds$assay_names()))   # detected despite generic name
  expect_equal(unname(colSums(eds$assay("tpm"))), rep(1e6, eds$n_samples()), tolerance = 1e-6)
  # FPKM matches the closed form counts*1e9 / length / libSize
  cts <- eds$assay("counts"); fw <- setNames(rd$gene_length, rownames(rd))[rownames(cts)]
  ref <- sweep(sweep(cts * 1e9, 1, fw, "/"), 2, colSums(cts), "/")
  expect_lt(max(abs(eds$assay("fpkm") - ref)), 1e-6)
})

test_that("no gene-length column -> FPKM/TPM gracefully absent (no ezRun needed)", {
  skip_if_no_anndataR()
  p <- make_toy_pieces()
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(p$counts, out = out, omics = "rnaseq", coldata = p$coldata)  # no lengths
  eds <- new_explore_dataset(out)
  expect_false(any(c("tpm", "fpkm") %in% eds$assay_names()))
})
