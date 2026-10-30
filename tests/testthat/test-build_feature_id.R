test_that("build_feature_id returns protein_Id directly when it is unique (protein-level data)", {
  df <- data.frame(
    protein_Id = c("P1", "P2", "P3"),
    stringsAsFactors = FALSE
  )
  expect_equal(build_feature_id(df), c("P1", "P2", "P3"))
})

test_that("build_feature_id falls back to protein_Id_site for phospho data", {
  df <- data.frame(
    protein_Id = c("P1", "P1", "P2"),
    site = c("S1", "S2", "S1"),
    stringsAsFactors = FALSE
  )
  expect_equal(build_feature_id(df), c("P1_S1", "P1_S2", "P2_S1"))
})

test_that("build_feature_id falls back to protein_Id_peptide_Id for peptide data", {
  df <- data.frame(
    protein_Id = c("P1", "P1", "P2"),
    peptide_Id = c("AAA", "BBB", "AAA"),
    stringsAsFactors = FALSE
  )
  expect_equal(build_feature_id(df), c("P1_AAA", "P1_BBB", "P2_AAA"))
})

test_that("build_feature_id prefers site over peptide_Id when both are present and site disambiguates", {
  df <- data.frame(
    protein_Id = c("P1", "P1", "P2"),
    site = c("S1", "S2", "S1"),
    peptide_Id = c("AAA", "BBB", "AAA"),
    stringsAsFactors = FALSE
  )
  # site is listed first in the intersect() loop, so if it disambiguates, it wins
  expect_equal(build_feature_id(df), c("P1_S1", "P1_S2", "P2_S1"))
})

test_that("build_feature_id strips tildes from sub-feature columns when combining", {
  df <- data.frame(
    protein_Id = c("P1", "P1"),
    site = c("S1~lfq~extra", "S2~lfq~extra"),
    stringsAsFactors = FALSE
  )
  expect_equal(build_feature_id(df), c("P1_S1_lfq_extra", "P1_S2_lfq_extra"))
})

test_that("build_feature_id falls through to cleaned rownames when no column disambiguates", {
  df <- data.frame(
    protein_Id = c("P1", "P1"),
    stringsAsFactors = FALSE
  )
  rownames(df) <- c("P1~lfq~foo", "P1~lfq~bar")
  out <- build_feature_id(df)
  # gsub("~[^~]+~", "_", rownames(df)) collapses the middle "~lfq~" to "_"
  expect_equal(out, c("P1_foo", "P1_bar"))
})

test_that("build_feature_id handles single-row data frame", {
  df <- data.frame(protein_Id = "P1", stringsAsFactors = FALSE)
  expect_equal(build_feature_id(df), "P1")
})
