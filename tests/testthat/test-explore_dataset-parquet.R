# -----------------------------------------------------------------------------
# ExploreDatasetParquet + export_h5ad_to_parquet() tests
#
# The exploreDElight backend (plans/exploredelight_tauri_plan.md). The whole
# point of this backend is that it must be INDISTINGUISHABLE from the h5ad
# backend through the accessor contract — Light runs the same server modules —
# so the central test is accessor-for-accessor equivalence on one source file,
# in the same spirit as test-cross_api_equivalence.R.
#
# Everything here goes through export_h5ad_to_parquet(), the R reference
# implementation of the directory contract. The Rust writer in the Tauri app
# has to match it; when it drifts, these are the tests that should fail.
# -----------------------------------------------------------------------------

skip_if_no_parquet <- function() {
  testthat::skip_if_not_installed("anndataR")
  if (!requireNamespace("nanoparquet", quietly = TRUE) &&
      !requireNamespace("arrow", quietly = TRUE)) {
    testthat::skip("neither nanoparquet nor arrow is installed")
  }
}

# A toy canonical .h5ad + its parquet export. Built once per test that needs it
# (cheap: 300 features x 6 samples).
toy_pair <- function(env = parent.frame()) {
  dds <- make_toy_dds()
  res <- DESeq2::results(dds, contrast = c("condition", "treated", "control"))
  h5  <- withr::local_tempfile(fileext = ".h5ad", .local_envir = env)
  build_explore_h5ad(dds, out = h5, omics = "rnaseq",
                     de = list(treated_vs_control = res),
                     organism = "Homo sapiens")
  pq <- withr::local_tempdir(.local_envir = env)
  export_h5ad_to_parquet(h5, pq)
  list(h5 = h5, pq = pq)
}

# ---- 1. the central claim: accessor equivalence -----------------------------

test_that("parquet backend is accessor-identical to the h5ad backend", {
  skip_if_no_parquet()
  p <- toy_pair()
  a <- ExploreDatasetH5AD$new(p$h5)
  b <- ExploreDatasetParquet$new(p$pq)

  # identity / config
  expect_identical(a$omics_class(),        b$omics_class())
  expect_identical(a$data_type(),          b$data_type())
  expect_identical(a$feature_id_kind(),    b$feature_id_kind())
  expect_identical(a$organism(),           b$organism())
  expect_identical(a$factor_names(),       b$factor_names())
  expect_identical(a$primary_grouping(),   b$primary_grouping())
  expect_identical(a$secondary_grouping(), b$secondary_grouping())
  expect_identical(a$contrasts(),          b$contrasts())
  expect_identical(a$primary_contrast(),   b$primary_contrast())
  expect_identical(a$assay_names(),        b$assay_names())
  expect_identical(a$enrichment_kinds(),   b$enrichment_kinds())
  expect_identical(a$has_pvalue(),         b$has_pvalue())

  # dimensions and ids
  expect_identical(a$n_samples(),   b$n_samples())
  expect_identical(a$n_features(),  b$n_features())
  expect_identical(a$sample_ids(),  b$sample_ids())
  expect_identical(a$feature_ids(), b$feature_ids())

  # tables
  expect_equal(a$samples(),  b$samples())
  expect_equal(a$features(), b$features())
  for (cc in a$contrasts()) expect_equal(a$de(cc), b$de(cc))
})

test_that("every assay is numerically identical, not merely close", {
  skip_if_no_parquet()
  p <- toy_pair()
  a <- ExploreDatasetH5AD$new(p$h5)
  b <- ExploreDatasetParquet$new(p$pq)

  expect_gt(length(b$assay_names()), 0L)
  for (nm in b$assay_names()) {
    ma <- a$assay(nm); mb <- b$assay(nm)
    expect_identical(dim(ma), dim(mb), info = nm)
    expect_identical(dimnames(ma), dimnames(mb), info = nm)
    # max|d| == 0: a float round-trip through parquet must not perturb a
    # single value. "Close enough" here would hide a real precision loss.
    expect_identical(max(abs(ma - mb), na.rm = TRUE), 0, info = nm)
  }
})

# ---- 2. the trap the manifest exists to prevent -----------------------------

test_that("factor levels are whatever the h5ad backend says, not re-sorted", {
  skip_if_no_parquet()
  p <- toy_pair()
  a <- ExploreDatasetH5AD$new(p$h5)
  b <- ExploreDatasetParquet$new(p$pq)
  # Fidelity is the contract: the parquet backend must agree with the h5ad
  # backend about level order, whatever that order happens to be.
  for (f in a$factor_names())
    expect_identical(a$factor_levels(f), b$factor_levels(f), info = f)
})

test_that("a non-alphabetical factor order survives h5ad -> parquet", {
  skip_if_no_parquet()
  # Parquet has no factor type, so without the manifest's factor_levels block
  # the reader falls back to sort(unique(x)) and silently alphabetises every
  # experimental group -- turning a "0h, 6h, 24h" time course into
  # "0h, 24h, 6h". Levels chosen so alphabetical, mixedsort and the intended
  # order are all different, so only honouring the declared levels passes.
  lv  <- c("zeta", "beta", "alpha")
  cd  <- data.frame(timepoint = factor(rep(lv, each = 2), levels = lv))
  m   <- matrix(stats::rpois(50 * 6, 100), nrow = 50,
                dimnames = list(sprintf("G%02d", 1:50), sprintf("S%d", 1:6)))
  rownames(cd) <- colnames(m)
  h5 <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(m, out = h5, omics = "rnaseq", coldata = cd,
                     organism = "Homo sapiens")
  pq <- withr::local_tempdir()
  export_h5ad_to_parquet(h5, pq)

  a <- ExploreDatasetH5AD$new(h5)
  b <- ExploreDatasetParquet$new(pq)
  expect_identical(a$factor_levels("timepoint"), lv)
  expect_identical(b$factor_levels("timepoint"), lv)
  expect_false(identical(b$factor_levels("timepoint"), sort(lv)))
})

# ---- 3. the directory contract ----------------------------------------------

test_that("export writes the documented layout and a readable manifest", {
  skip_if_no_parquet()
  p <- toy_pair()
  expect_true(file.exists(file.path(p$pq, "manifest.json")))
  expect_true(file.exists(file.path(p$pq, "obs.parquet")))
  expect_true(file.exists(file.path(p$pq, "var.parquet")))
  expect_true(file.exists(file.path(p$pq, "uns.json")))
  expect_true(dir.exists(file.path(p$pq, "layers")))

  m <- jsonlite::fromJSON(file.path(p$pq, "manifest.json"), simplifyVector = TRUE,
                          simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
  expect_identical(m$layer_orientation, "features_x_samples")
  expect_true("counts" %in% m$layers)
  expect_gt(m$n_var, 0L)
  expect_gt(m$n_obs, 0L)
  # Every declared layer/varm slot must actually be on disk — a manifest that
  # promises a file the reader then cannot find is the failure mode most likely
  # to appear when the Rust writer drifts from this reference.
  for (nm in m$layers)
    expect_true(file.exists(file.path(p$pq, "layers", paste0(nm, ".parquet"))))
  for (nm in (m$varm %||% character()))
    expect_true(file.exists(file.path(p$pq, "varm", paste0(nm, ".parquet"))))
})

test_that("layers are written features-as-rows with an _index column", {
  skip_if_no_parquet()
  p <- toy_pair()
  b <- ExploreDatasetParquet$new(p$pq)
  d <- exploreDE:::.pq_read(file.path(p$pq, "layers", "counts.parquet"))
  expect_true("_index" %in% colnames(d))
  # rows = features, data columns = samples (the long axis goes down the file)
  expect_identical(nrow(d), b$n_features())
  expect_identical(ncol(d) - 1L, b$n_samples())
})

# ---- 4. dispatch and identity -----------------------------------------------

test_that("new_explore_dataset() routes a parquet directory to the backend", {
  skip_if_no_parquet()
  p <- toy_pair()
  eds <- new_explore_dataset(p$pq)
  expect_s3_class(eds, "ExploreDatasetParquet")
  expect_identical(eds$backend_kind(), "parquet")
  expect_identical(eds$source_dir(), p$pq)
  expect_silent(validate_explore_dataset(eds))
})

test_that("the parquet backend inherits the h5ad backend's behaviour", {
  skip_if_no_parquet()
  p <- toy_pair()
  b <- ExploreDatasetParquet$new(p$pq)
  # Inheritance is the design: config derivation and DE assembly are shared,
  # so a change to ExploreDatasetH5AD must reach this backend too.
  expect_s3_class(b, "ExploreDatasetH5AD")
  expect_s3_class(b, "ExploreDataset")
})

# ---- 5. failure modes -------------------------------------------------------

test_that("bad input is rejected with a message that says what was wrong", {
  skip_if_no_parquet()
  expect_error(ExploreDatasetParquet$new("/no/such/directory"),
               "parquet directory")
  expect_error(ExploreDatasetParquet$new(42), "parquet directory")

  # A directory that exists but is missing the files.
  empty <- withr::local_tempdir()
  expect_error(ExploreDatasetParquet$new(empty))
})

test_that("a parquet table without _index is rejected rather than mis-read", {
  skip_if_no_parquet()
  d <- withr::local_tempdir()
  f <- file.path(d, "bad.parquet")
  if (requireNamespace("nanoparquet", quietly = TRUE)) {
    nanoparquet::write_parquet(data.frame(a = 1:3), f)
  } else {
    arrow::write_parquet(data.frame(a = 1:3), f)
  }
  # Silently inventing rownames here would misalign every feature downstream,
  # so this must be an error and not a warning.
  expect_error(exploreDE:::.pq_read_indexed(f), "_index")
})

# ---- 6. lazy enrichment -----------------------------------------------------

test_that("enrichment payloads are hoisted out of uns.json, not inlined", {
  skip_if_no_parquet()
  p <- toy_pair()
  # The toy fixture carries no enrichment, so assert the shape of the rule
  # rather than the presence of files: uns.json must not be carrying megabytes
  # of base64. Measured on my_study.h5ad it was 13 MB of a 29 MB export,
  # parsed before the first plot could render.
  expect_lt(file.info(file.path(p$pq, "uns.json"))$size, 1e6)
  m <- jsonlite::fromJSON(file.path(p$pq, "manifest.json"), simplifyVector = TRUE,
                          simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
  expect_true("enrichment_files" %in% names(m))
})

test_that("hoisting keeps the enrichment KEYS so has_enrichment() still works", {
  skip_if_no_parquet()
  # .build_config() derives enrichment_kinds from names(uns$enrichment_native).
  # The writer empties the payloads but must keep the keys, or a dataset with
  # enrichment would silently report having none.
  d   <- withr::local_tempdir()
  uns <- list(schema_version = "exploreDE/1.0",
              enrichment_native = list(mycontrast = list(ora = "AAAA", gsea = "BBBB")))
  ad  <- anndataR::AnnData(
    # AnnData is samples x features, not the SE orientation.
    X = matrix(as.double(1:6), nrow = 2,
               dimnames = list(c("s1","s2"), c("a","b","c"))),
    obs = data.frame(row.names = c("s1","s2"), factor_g = factor(c("x","y"))),
    var = data.frame(row.names = c("a","b","c"), description = c("A","B","C")),
    uns = uns)
  export_h5ad_to_parquet(ad, d)

  written <- jsonlite::fromJSON(file.path(d, "uns.json"), simplifyVector = TRUE,
                                simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
  expect_identical(names(written$enrichment_native$mycontrast), c("ora", "gsea"))
  expect_identical(written$enrichment_native$mycontrast$ora, "")   # key kept, payload gone
  expect_true(file.exists(file.path(d, "enrichment", "mycontrast__ora.b64")))
  expect_identical(readLines(file.path(d, "enrichment", "mycontrast__ora.b64"),
                             warn = FALSE), "AAAA")
})

# ---- 7. source identity / staleness -----------------------------------------

test_that("the manifest records what the export came from", {
  skip_if_no_parquet()
  p   <- toy_pair()
  m   <- jsonlite::fromJSON(file.path(p$pq, "manifest.json"), simplifyVector = TRUE,
                            simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
  expect_identical(m$source$basename, basename(p$h5))
  expect_gt(m$source$size, 0)
  expect_true(nzchar(m$source$md5))
  expect_true(ExploreDatasetParquet$new(p$pq)$verify_source())
})

test_that("a stale export warns instead of passing off old data as current", {
  skip_if_no_parquet()
  p <- toy_pair()
  expect_silent(ExploreDatasetParquet$new(p$pq))       # fresh: no warning

  # Rewrite the source so size and mtime both move. A derived cache that reads
  # as current after its source changed is the failure mode worth shouting
  # about: every number would be silently out of date.
  writeBin(as.raw(rep(0L, 4096)), p$h5)
  expect_warning(ExploreDatasetParquet$new(p$pq), "out of date")
  expect_false(isTRUE(suppressWarnings(ExploreDatasetParquet$new(p$pq))$verify_source()))
})

test_that("exporting a live AnnData records no source rather than a fake one", {
  skip_if_no_parquet()
  d  <- withr::local_tempdir()
  ad <- anndataR::AnnData(
    # AnnData is samples x features, not the SE orientation.
    X = matrix(as.double(1:6), nrow = 2,
               dimnames = list(c("s1","s2"), c("a","b","c"))),
    obs = data.frame(row.names = c("s1","s2"), factor_g = factor(c("x","y"))),
    var = data.frame(row.names = c("a","b","c"), description = c("A","B","C")),
    uns = list(schema_version = "exploreDE/1.0"))
  export_h5ad_to_parquet(ad, d)
  # Nothing to be stale against, so verify_source() must say "unknown" (NA),
  # not invent a pass or a fail.
  expect_true(is.na(ExploreDatasetParquet$new(d)$verify_source()))
})
