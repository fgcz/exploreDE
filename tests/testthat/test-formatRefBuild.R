test_that("formatRefBuild formats a GENCODE human reference build path", {
  out <- formatRefBuild("Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30")
  expect_match(out, "GENCODE")
  expect_match(out, "human")
  expect_match(out, "GRCh38\\.p13")
  expect_match(out, "Release 42")
})

test_that("formatRefBuild formats a GENCODE mouse reference build path", {
  out <- formatRefBuild("Mus_musculus/GENCODE/GRCm39/Annotation/Release_M31-2023-01-30")
  expect_match(out, "GENCODE")
  expect_match(out, "mouse")
  expect_match(out, "GRCm39")
  expect_match(out, "Release M31")
})

test_that("formatRefBuild returns a non-path string unchanged", {
  expect_equal(formatRefBuild("GRCh38.p13"), "GRCh38.p13")
  expect_equal(formatRefBuild("the reference genome"), "the reference genome")
})

test_that("formatRefBuild converts species underscores to spaces for unknown species", {
  out <- formatRefBuild("Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Release_104")
  expect_match(out, "Rattus norvegicus")
})

test_that("formatRefBuild omits the release suffix when not present", {
  out <- formatRefBuild("Homo_sapiens/GENCODE/GRCh38.p13/Annotation")
  expect_false(grepl("Release", out))
})
