# -----------------------------------------------------------------------------
# Multi-DE dataset links: multi_resolve_dataset(), ldap_allows_project(),
# multi_load_de_result() (R/utils_multi_de.R)
# -----------------------------------------------------------------------------

# Two fake data roots with SUSHI-shaped result directories.
fake_roots <- function(env = parent.frame()) {
  base <- withr::local_tempdir(.local_envir = env)
  gstore <- file.path(base, "gstore", "projects")
  course <- file.path(base, "course", "projects")
  mk <- function(dir, file = "deResult.rds") {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(file)) file.create(file.path(dir, file))
  }
  mk(file.path(gstore, "p1234", "DESeq2_KO--over--WT_2024", "KO--over--WT"))
  mk(file.path(gstore, "p9999", "DESeq2_A--over--B_2024", "A--over--B"))
  mk(file.path(gstore, "p1234", "DESeq2_qs2_2024", "X--over--Y"), "deResult.qs2")
  mk(file.path(gstore, "p1234", "DESeq2_empty_2024", "E--over--F"), NULL)
  mk(file.path(gstore, "log", "DESeq2_x", "x"))
  mk(file.path(course, "p1000", "DESeq2_SDY--over--Control", "SDY--over--Control"))
  mk(file.path(base, "gstore", "projects_evil", "p1234", "DESeq2_z", "z"))
  list(base = base, gstore = gstore, course = course, roots = c(gstore, course))
}

test_that("a gStore-relative link resolves to its result directory and project", {
  r <- fake_roots()
  res <- multi_resolve_dataset("p1234/DESeq2_KO--over--WT_2024/KO--over--WT/", r$roots)
  expect_equal(res$path, file.path(r$gstore, "p1234/DESeq2_KO--over--WT_2024/KO--over--WT"))
  expect_equal(res$root, r$gstore)
  expect_equal(res$project, "p1234")
})

test_that("an exploreDE URL and stray spaces are reduced to the ?data= path", {
  r <- fake_roots()
  url <- "https://fgcz-shiny.uzh.ch/exploreDEG/?data=p1234/DESeq2_KO--over--WT_2024/KO--over--WT/"
  expect_equal(multi_resolve_dataset(url, r$roots)$project, "p1234")
  expect_equal(multi_resolve_dataset(" p1234/DESeq2_KO--over--WT_2024/ KO--over--WT ", r$roots)$project, "p1234")
})

test_that("'..' is refused, even when the target exists", {
  r <- fake_roots()
  expect_null(multi_resolve_dataset("p1234/../p9999/DESeq2_A--over--B_2024/A--over--B", r$roots))
  expect_null(multi_resolve_dataset("../projects/p9999/DESeq2_A--over--B_2024/A--over--B", r$roots))
  expect_null(multi_resolve_dataset(
    "https://x/?data=p1234/../p9999/DESeq2_A--over--B_2024/A--over--B",
    r$roots
  ))
})

test_that("absolute paths must lie under a root, not merely share its prefix", {
  r <- fake_roots()
  ok <- file.path(r$gstore, "p1234/DESeq2_KO--over--WT_2024/KO--over--WT")
  expect_equal(multi_resolve_dataset(ok, r$roots)$path, ok)
  expect_null(multi_resolve_dataset(file.path(r$base, "gstore/projects_evil/p1234/DESeq2_z/z"), r$roots))
  expect_null(multi_resolve_dataset("/etc", r$roots))
})

test_that("the first component must be a project, and a DE result must be present", {
  r <- fake_roots()
  expect_null(multi_resolve_dataset("log/DESeq2_x/x", r$roots))
  expect_null(multi_resolve_dataset("p1234/DESeq2_empty_2024/E--over--F", r$roots))
  expect_null(multi_resolve_dataset("p1234/does/not/exist", r$roots))
  expect_equal(multi_resolve_dataset("p1234/DESeq2_qs2_2024/X--over--Y", r$roots)$project, "p1234")
})

test_that("the second root is tried when the first holds nothing", {
  r <- fake_roots()
  res <- multi_resolve_dataset("p1000/DESeq2_SDY--over--Control/SDY--over--Control", r$roots)
  expect_equal(res$root, r$course)
})

test_that("non-links are rejected", {
  r <- fake_roots()
  for (bad in list(NULL, NA_character_, "", "   ", 42, c("a", "b"))) {
    expect_null(multi_resolve_dataset(bad, r$roots))
  }
})

test_that("ldap_allows_project mirrors the ?data= rule", {
  expect_true(ldap_allows_project(list(roles = c("R_2", "R_3"), projects = character()), "p1"))
  expect_true(ldap_allows_project(list(roles = "R_3", projects = c("p1", "p2")), "p2"))
  expect_false(ldap_allows_project(list(roles = "R_3", projects = c("p1", "p2")), "p3"))
  expect_false(ldap_allows_project(list(roles = character(), projects = "p1"), "p1"))
  expect_false(ldap_allows_project(list(), "p1"))
})

test_that("multi_load_de_result returns NULL without a DE result", {
  r <- fake_roots()
  expect_null(multi_load_de_result(file.path(r$gstore, "p1234/DESeq2_empty_2024/E--over--F")))
})

test_that("multi_load_de_result loads a real SUSHI DESeq2 result", {
  path <- "/srv/GT/analysis/course_sushi/public/gstore/projects/p1000/DESeq2_SDY--over--Control_2025-01-29--14-26-27/SDY--over--Control"
  skip_if_not(file.exists(file.path(path, "deResult.rds")), "course data not mounted")
  res <- multi_resolve_dataset(
    sub("^.*/projects/", "", path),
    c("/srv/gstore/projects", "/srv/GT/analysis/course_sushi/public/gstore/projects")
  )
  expect_equal(res$path, path)
  out <- multi_load_de_result(res$path)
  expect_true("log2Ratio" %in% colnames(out$seqAnno))
})

test_that("multi_load_de_result gives the demo result's record, samples included", {
  path <- system.file("examples/genomics", package = "exploreDE")
  rec <- multi_load_de_result(path)
  expect_identical(rec$dataType, "RNASeq")
  expect_true(all(c("feature_id", "log2Ratio", "pValue", "fdr", "usedInTest") %in% colnames(rec$seqAnno)))
  expect_equal(nrow(rec$dataset), 9L)
  expect_true(all(rec$factors %in% colnames(rec$dataset)))
  expect_match(rec$param$comparison, "--over--", fixed = TRUE)
  expect_true(all(c("Raw", "Normalised + Log2") %in% names(rec$countList)))
  expect_true("Label" %in% colnames(.res_df(rec$ora[[1]][[1]])))
})

# -----------------------------------------------------------------------------
# Uploaded files in Multi DE: multi_upload_label(), multi_records_from_dataset()
# -----------------------------------------------------------------------------

test_that("multi_upload_label is the file name made id-safe and unique", {
  expect_identical(multi_upload_label("/tmp/upload/GSE1234 condition.1.h5ad"), "GSE1234_condition_1")
  expect_identical(multi_upload_label("study-2.h5ad"), "study-2")
  expect_identical(multi_upload_label("study.h5ad", existing = c("study", "other")), "study_1")
  expect_identical(multi_upload_label("study.h5ad", existing = c("study", "study_1")), "study_2")
  expect_identical(multi_upload_label("(((.h5ad"), "upload")
})

test_that("multi_records_from_dataset gives one record per contrast, named and sourced by origin", {
  eds <- new_explore_dataset(readRDS(system.file("examples/genomics/deResult.rds", package = "exploreDE")))
  con <- eds$contrasts()
  own <- multi_records_from_dataset(eds)
  expect_identical(names(own), con)
  expect_identical(own[[1]]$source, list(kind = "contrast", contrast = con[1]))
  expect_match(own[[1]]$param$comparison, "--over--", fixed = TRUE)

  up <- multi_records_from_dataset(eds, label = "study2")
  expect_identical(names(up), paste0("study2_", con))
  expect_identical(up[[1]]$source, list(kind = "upload", file = "study2", contrast = con[1]))

  labelled <- multi_records_from_dataset(eds, label = "main", uploaded = FALSE)
  expect_identical(names(labelled), paste0("main_", con))
  expect_identical(labelled[[1]]$source, list(kind = "contrast", contrast = con[1]))
  expect_identical(attr(labelled, "failed"), character())
})

test_that("multi_records_from_dataset leaves out a contrast it cannot build and names it", {
  se <- readRDS(system.file("examples/genomics/deResult.rds", package = "exploreDE"))
  Broken <- R6::R6Class(
    "Broken",
    inherit = ExploreDatasetSE_RNASeq,
    public = list(
      contrasts = function() c(super$contrasts(), "nope"),
      contrast_meta = function(contrast) {
        if (contrast == "nope") {
          stop("no such contrast")
        }
        super$contrast_meta(contrast)
      }
    )
  )
  broken <- Broken$new(se)
  con <- setdiff(broken$contrasts(), "nope")
  recs <- multi_records_from_dataset(broken, label = "x")
  expect_identical(names(recs), paste0("x_", con))
  expect_identical(names(attr(recs, "failed")), "nope")
  expect_match(attr(recs, "failed")[["nope"]], "no such contrast")
})
