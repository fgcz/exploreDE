test_that("convert_genomics_se on RNASeq gene-level fixture returns expected structure", {
  fx <- fx_rnaseq_gene()
  out <- convert_genomics_se(fx$se, fx$dataDir)

  # --- Top-level shape -------------------------------------------------------
  expect_type(out, "list")
  expect_equal(out$dataType, "RNASeq")
  expect_true(all(
    c(
      "se",
      "param",
      "seqAnno",
      "dataset",
      "factors",
      "factorLevels",
      "groupingName",
      "countList",
      "genes",
      "design",
      "designLevels",
      "goObjects",
      "allPathways"
    ) %in%
      names(out)
  ))

  # --- seqAnno ---------------------------------------------------------------
  expect_s3_class(out$seqAnno, "data.frame")
  expect_true("feature_id" %in% colnames(out$seqAnno))
  expect_true("pValue" %in% colnames(out$seqAnno))
  expect_true("log2Ratio" %in% colnames(out$seqAnno))
  # Sorted by raw pValue (ascending)
  expect_false(is.unsorted(out$seqAnno$pValue, na.rm = TRUE))
  # feature_ids are unique
  expect_equal(anyDuplicated(out$seqAnno$feature_id), 0L)

  # --- countList -------------------------------------------------------------
  expect_type(out$countList, "list")
  expect_true(all(c("Raw", "Normalised", "Normalised + Log2", "VST") %in% names(out$countList)))
  # TPM / FPKM present when refBuild is set
  if (!is.null(out$param$refBuild)) {
    expect_true(all(c("TPM", "FPKM") %in% names(out$countList)))
  }
  # Rownames align with seqAnno$feature_id
  for (m in names(out$countList)) {
    expect_setequal(rownames(out$countList[[m]]), out$seqAnno$feature_id)
  }

  # --- Sample metadata -------------------------------------------------------
  expect_s3_class(out$dataset, "data.frame")
  expect_true("Name" %in% colnames(out$dataset))
  # countList columns line up with dataset rows
  for (m in names(out$countList)) {
    expect_equal(ncol(out$countList[[m]]), nrow(out$dataset))
  }

  # --- Design ----------------------------------------------------------------
  expect_type(out$design, "character")
  expect_type(out$designLevels, "character")
  expect_length(out$designLevels, 2)

  # --- Genes vector ----------------------------------------------------------
  expect_type(out$genes, "character")
  expect_setequal(out$genes, out$seqAnno$feature_id)
})

test_that("convert_genomics_se on CountQC fixture returns expected structure", {
  fx <- fx_countqc()
  out <- convert_genomics_se(fx$se, fx$dataDir)

  expect_equal(out$dataType, "CountQC")
  # CountQC should have an empty/placeholder design (no DE results)
  expect_equal(out$design, "CountQC")
  expect_null(out$designLevels)

  # Still has the core bits
  expect_s3_class(out$seqAnno, "data.frame")
  expect_true("feature_id" %in% colnames(out$seqAnno))
  expect_true(all(c("Raw", "Normalised", "Normalised + Log2", "VST") %in% names(out$countList)))
})

test_that("convert_genomics_se on CrisprCountQC fixture adds sgRNA_to_gene lookup", {
  fx <- fx_crispr_countqc()
  out <- convert_genomics_se(fx$se, fx$dataDir)

  expect_equal(out$dataType, "CrisprCountQC")
  expect_true("sgRNA_to_gene" %in% names(out))
  expect_type(out$sgRNA_to_gene, "character")
  # sgRNA_to_gene maps every feature_id → a gene name
  expect_setequal(names(out$sgRNA_to_gene), out$seqAnno$feature_id)
})

test_that("convert_genomics_se on RNASeq isoform fixture uses transcript_id for feature ordering", {
  fx <- fx_rnaseq_isoform()
  out <- convert_genomics_se(fx$se, fx$dataDir)

  expect_equal(out$dataType, "RNASeq")
  expect_equal(out$param$featureLevel, "isoform")
  expect_true("transcript_id" %in% colnames(out$seqAnno))
  # Isoform feature_ids are cleaned to replace special characters with "."
  expect_false(any(grepl("[:/()-]| ", rownames(out$countList[[1]]))))
})
