# --- is_url -----------------------------------------------------------------
test_that("is_url recognises http / https URLs and rejects paths", {
  expect_true(is_url("http://example.com/file.rds"))
  expect_true(is_url("https://fgcz-ms.uzh.ch/public/pStore/foo.rds"))
  expect_false(is_url("/srv/gstore/projects/p3009/file.rds"))
  expect_false(is_url("relative/path.rds"))
  expect_false(is_url(""))
  # "http" in middle of a path must not trigger
  expect_false(is_url("/tmp/http_dump/file"))
})


# --- parse_gene_text --------------------------------------------------------
test_that("parse_gene_text splits by whitespace, comma, slash, tab, newline", {
  expect_equal(parse_gene_text("GENE1 GENE2"),          c("GENE1", "GENE2"))
  expect_equal(parse_gene_text("GENE1,GENE2,GENE3"),    c("GENE1", "GENE2", "GENE3"))
  expect_equal(parse_gene_text("GENE1/GENE2"),          c("GENE1", "GENE2"))
  expect_equal(parse_gene_text("GENE1\tGENE2"),         c("GENE1", "GENE2"))
  expect_equal(parse_gene_text("GENE1\nGENE2"),         c("GENE1", "GENE2"))
  expect_equal(parse_gene_text("A, B/C\tD\nE"),         c("A", "B", "C", "D", "E"))
})

test_that("parse_gene_text returns character(0) for empty / NULL input", {
  expect_equal(parse_gene_text(""),   character(0))
  expect_equal(parse_gene_text(NULL), character(0))
})

test_that("parse_gene_text collapses runs of separators", {
  expect_equal(parse_gene_text("GENE1   GENE2"),    c("GENE1", "GENE2"))
  expect_equal(parse_gene_text("GENE1,,,GENE2"),    c("GENE1", "GENE2"))
})


# --- getPackageVersions -----------------------------------------------------
test_that("getPackageVersions returns a valid version string for installed packages", {
  vers <- getPackageVersions(c("base", "stats"))
  expect_length(vers, 2)
  expect_named(vers, c("base", "stats"))
  expect_true(all(grepl("^[0-9]", vers)))
})

test_that("getPackageVersions returns 'version unavailable' for unknown packages", {
  vers <- getPackageVersions("this_package_definitely_does_not_exist_42")
  expect_equal(unname(vers), "version unavailable")
})
