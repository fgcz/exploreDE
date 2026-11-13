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
  expect_equal(parse_gene_text("GENE1 GENE2"), c("GENE1", "GENE2"))
  expect_equal(parse_gene_text("GENE1,GENE2,GENE3"), c("GENE1", "GENE2", "GENE3"))
  expect_equal(parse_gene_text("GENE1/GENE2"), c("GENE1", "GENE2"))
  expect_equal(parse_gene_text("GENE1\tGENE2"), c("GENE1", "GENE2"))
  expect_equal(parse_gene_text("GENE1\nGENE2"), c("GENE1", "GENE2"))
  expect_equal(parse_gene_text("A, B/C\tD\nE"), c("A", "B", "C", "D", "E"))
})

test_that("parse_gene_text returns character(0) for empty / NULL input", {
  expect_equal(parse_gene_text(""), character(0))
  expect_equal(parse_gene_text(NULL), character(0))
})

test_that("parse_gene_text collapses runs of separators", {
  expect_equal(parse_gene_text("GENE1   GENE2"), c("GENE1", "GENE2"))
  expect_equal(parse_gene_text("GENE1,,,GENE2"), c("GENE1", "GENE2"))
})


# --- match_feature_pattern --------------------------------------------------
pattern_candidates <- c(
  "Hspa1a",
  "Hspa1b",
  "Hsp90aa1",
  "HSPB1",
  "Gm12345",
  "1700001A01Rik",
  "Actb",
  "sp|P02768|ALBU_HUMAN",
  "sp|P00001|X_HUMAN"
)

test_that("match_feature_pattern takes * and ? as wildcards and everything else literally, ignoring case", {
  expect_equal(
    as.character(match_feature_pattern("hsp*", pattern_candidates)),
    c("Hspa1a", "Hspa1b", "Hsp90aa1", "HSPB1")
  )
  expect_equal(as.character(match_feature_pattern("*rik", pattern_candidates)), "1700001A01Rik")
  expect_equal(as.character(match_feature_pattern("Hspa1?", pattern_candidates)), c("Hspa1a", "Hspa1b"))
  # A bare name is that name only, not a substring
  expect_equal(as.character(match_feature_pattern("actb", pattern_candidates)), "Actb")
  expect_equal(as.character(match_feature_pattern("Hsp", pattern_candidates)), character(0))
  # Regular expression characters in a name are literal
  expect_equal(as.character(match_feature_pattern("sp|P02768|ALBU_HUMAN", pattern_candidates)), "sp|P02768|ALBU_HUMAN")
  expect_equal(
    as.character(match_feature_pattern("sp|*_HUMAN", pattern_candidates)),
    c("sp|P02768|ALBU_HUMAN", "sp|P00001|X_HUMAN")
  )
  expect_equal(as.character(match_feature_pattern("Hsp(", pattern_candidates)), character(0))
})

test_that("match_feature_pattern takes a regular expression when asked", {
  rx <- function(p) as.character(match_feature_pattern(p, pattern_candidates, regex = TRUE))
  expect_equal(rx("^Hsp[a9]"), c("Hspa1a", "Hspa1b", "Hsp90aa1"))
  expect_equal(rx("rik$"), "1700001A01Rik")
  expect_equal(rx("Actb|HSPB1"), c("HSPB1", "Actb"))
  # unanchored, like grep
  expect_equal(rx("hsp"), c("Hspa1a", "Hspa1b", "Hsp90aa1", "HSPB1"))
})

test_that("match_feature_pattern returns every match and their count", {
  m <- match_feature_pattern("hsp*", pattern_candidates)
  expect_equal(as.character(m), c("Hspa1a", "Hspa1b", "Hsp90aa1", "HSPB1"))
  expect_identical(attr(m, "n_matched"), 4L)
  expect_identical(attr(match_feature_pattern("zzz*", pattern_candidates), "n_matched"), 0L)
})

test_that("match_feature_pattern matches nothing for an empty, NULL or invalid pattern", {
  for (p in list("", "   ", NULL)) {
    m <- match_feature_pattern(p, pattern_candidates)
    expect_equal(as.character(m), character(0))
    expect_identical(attr(m, "n_matched"), 0L)
  }
  # TRE warns before it errors on the pattern; neither reaches the user
  expect_no_warning(bad <- match_feature_pattern("Hsp(", pattern_candidates, regex = TRUE))
  expect_equal(as.character(bad), character(0))
  expect_identical(attr(bad, "n_matched"), NA_integer_)
  expect_identical(attr(match_feature_pattern("hsp*", character(0)), "n_matched"), 0L)
})

test_that("cap_pattern_features caps the features the matched names resolve to", {
  # Protein mode: each name stands for 5 peptides
  resolve <- function(x) as.vector(outer(paste0(x, "_pep"), 1:5, paste0))
  names <- sprintf("P%04d", 1:2000)
  ids <- cap_pattern_features(names, resolve)
  expect_length(ids, 1000)
  expect_identical(attr(ids, "n_features"), 10000L)
  few <- cap_pattern_features(names[1:3], resolve)
  expect_length(few, 15)
  expect_identical(attr(few, "n_features"), 15L)
  expect_identical(attr(cap_pattern_features(character(0), stop), "n_features"), 0L)
})

test_that("fb_pattern_message describes the matches", {
  msg <- function(p, ...) fb_pattern_message(p, match_feature_pattern(p, pattern_candidates, ...))
  expect_identical(msg(""), "")
  expect_identical(msg("actb"), "1 name matches.")
  expect_identical(msg("hsp*"), "4 names match.")
  capped <- structure(sprintf("f%d", 1:1000), n_features = 1250L)
  expect_identical(
    fb_pattern_message("hsp*", match_feature_pattern("hsp*", pattern_candidates), capped),
    "4 names match (1250 features); the first 1000 features were added."
  )
  expect_identical(msg("Hsp(", regex = TRUE), "Not a valid pattern.")
  # a plain vector without the attribute counts its elements
  expect_identical(fb_pattern_message("hsp*", c("Hspa1a", "Hspa1b")), "2 names match.")
})

test_that("fb_pattern_ui builds the pattern box, the regex tick box and the info line of a tab", {
  html <- as.character(fb_pattern_ui("boxplot"))
  expect_match(html, 'id="fb_pattern_boxplot"', fixed = TRUE)
  expect_match(html, 'id="fb_pattern_regex_boxplot"', fixed = TRUE)
  expect_match(html, 'id="fb_pattern_info_boxplot"', fixed = TRUE)
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

test_that("match_feature_pattern cannot be made to backtrack for long", {
  set.seed(1)
  n <- 500
  acc <- replicate(n, paste0(sample(c(LETTERS, 0:9), 6, TRUE), collapse = ""))
  ids <- paste0("sp|", acc, "|", replicate(n, paste0(sample(LETTERS, 5, TRUE), collapse = "")), "_HUMAN")
  expect_lt(system.time(match_feature_pattern("(.+)+X$", ids, regex = TRUE))[["elapsed"]], 1)
  expect_lt(system.time(match_feature_pattern("*?*?*?*?*?*?*?*?*?*?*Q", ids))[["elapsed"]], 1)
  # Regular expressions still match as before
  expect_identical(
    as.character(match_feature_pattern("^sp\\|A|_HUMAN$", c("sp|A1|X_HUMAN", "tr|B|Y_MOUSE"), regex = TRUE)),
    "sp|A1|X_HUMAN"
  )
})
