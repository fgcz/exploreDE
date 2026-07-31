# ============================================================================
# Save/restore session-state — pure helpers (R/utils_savestate.R)
# Covers the input skip-filter, the exact-dataset fingerprint + matcher, and the
# restore_default seeding logic. Shiny runtime (retry burst, renderUI seeding,
# widget restoration) is out of scope here — it needs a running app.
# ============================================================================

# ---- save_state_filter_inputs ----------------------------------------------

test_that("save_state_filter_inputs keeps restorable ids and drops the rest", {
  keep_ids <- c("boxplotFactor1", "pcaGroups", "heatmapCounts", "boxKeepBucket",
                "GroupColourGenotypeWT", "volcanoColourNotSignificant", "tabs",
                "contrastSelected", "filename_boxplot", "replicabilityAtLow", "oraType")
  drop_ids <- c("load_state_file", "importColourFile", "multiBoxplotGenes",
                "multiContrastLink1", "degTable_rows_selected",
                "selectedTable_ORA_rows_selected", "boxplotBrush", "dlPlot_boxplot",
                "download_state", "heatmapGoGo", "fb_reset_boxplot", "fb_add_brushed",
                "deQcScatterBrush", "sidebarMenu", "resetDEGTableSelectRows")
  x <- setNames(as.list(seq_along(c(keep_ids, drop_ids))), c(keep_ids, drop_ids))
  kept <- names(save_state_filter_inputs(x))
  expect_true(all(keep_ids %in% kept))
  expect_false(any(drop_ids %in% kept))
})

test_that("contrastSelected is kept but *_rows_selected echoes are dropped", {
  x <- list(contrastSelected = "A", selectedTable_ORA_rows_selected = 3L,
            degTable_rows_selected = 1L)
  kept <- names(save_state_filter_inputs(x))
  expect_true("contrastSelected" %in% kept)
  expect_false("selectedTable_ORA_rows_selected" %in% kept)
  expect_false("degTable_rows_selected" %in% kept)
})

test_that("skip regex has no empty-alternative over-match (José || bug guard)", {
  # Ordinary ids must NEVER be dropped by a stray empty alternative.
  ordinary <- c("boxplotFactor1", "heatmapZScore", "pcaCounts", "textSizeORA")
  x <- setNames(as.list(ordinary), ordinary)
  expect_equal(sort(names(save_state_filter_inputs(x))), sort(ordinary))
})

# ---- save_state_fingerprint -------------------------------------------------

fake_eds <- function() list(
  data_type    = function() "RNASeq",
  omics_class  = function() "genomics",
  n_samples    = function() 6L,
  feature_ids  = function() paste0("g", 1:100),
  factor_names = function() c("Sex", "Genotype"),   # deliberately unsorted
  contrasts    = function() c("A_vs_B")
)

test_that("save_state_fingerprint builds a sorted, complete identity", {
  fp <- save_state_fingerprint(fake_eds(), "/gstore/pX/oY/")
  expect_type(fp, "list")
  expect_equal(sort(names(fp)),
               sort(c("dataDir", "dataType", "omics", "nSamples", "nFeatures",
                      "factors", "contrasts")))
  expect_equal(fp$factors, c("Genotype", "Sex"))   # sorted -> order-independent
  expect_equal(fp$nFeatures, 100L)
  expect_equal(fp$nSamples, 6L)
  expect_equal(fp$dataDir, "/gstore/pX/oY/")
})

test_that("save_state_fingerprint is NULL-safe and never throws on a bad accessor", {
  expect_null(save_state_fingerprint(NULL))
  bad <- list(data_type = function() stop("boom"), omics_class = function() "g",
              n_samples = function() 3L, feature_ids = function() letters,
              factor_names = function() "F", contrasts = function() character())
  expect_silent(fp <- save_state_fingerprint(bad, NULL))
  expect_true(is.na(fp$dataType))    # errored accessor -> NA, not a crash
  expect_true(is.na(fp$dataDir))     # NULL dataDir -> NA
})

# ---- save_state_fp_match (exact-dataset gate) -------------------------------

base_fp <- list(dataDir = "/g/pX/", dataType = "RNASeq", omics = "genomics",
                nSamples = 6L, nFeatures = 100L, factors = c("Genotype", "Sex"),
                contrasts = "A_vs_B")

test_that("identical fingerprints match; any field difference refuses", {
  expect_true(save_state_fp_match(base_fp, base_fp))
  mut <- function(...) modifyList(base_fp, list(...))
  expect_false(save_state_fp_match(base_fp, mut(dataDir  = "/other/")))
  expect_false(save_state_fp_match(base_fp, mut(nSamples = 7L)))
  expect_false(save_state_fp_match(base_fp, mut(nFeatures = 99L)))
  expect_false(save_state_fp_match(base_fp, mut(factors  = "Genotype")))
  expect_false(save_state_fp_match(base_fp, mut(contrasts = character())))
  expect_false(save_state_fp_match(base_fp, mut(dataType = "CountQC")))
})

test_that("fp_match is NULL-safe; NA dataDir falls back to content match", {
  expect_false(save_state_fp_match(NULL, base_fp))
  expect_false(save_state_fp_match(base_fp, NULL))
  na_dd <- modifyList(base_fp, list(dataDir = NA_character_))
  expect_true(save_state_fp_match(na_dd, base_fp))   # content still matches
})

test_that("CountQC-style fingerprint (empty contrasts) self-matches", {
  cqc <- list(dataDir = "/g/c/", dataType = "CountQC", omics = "genomics",
              nSamples = 4L, nFeatures = 50L, factors = "Batch",
              contrasts = character())
  expect_true(save_state_fp_match(cqc, cqc))
})

# ---- restore_default_value --------------------------------------------------

test_that("restore_default_value returns default unless restoring with a seed", {
  restoring <- shiny::reactiveVal(FALSE)
  seed <- shiny::reactiveValues()
  shiny::isolate({
    expect_equal(restore_default_value("a", "DEF", restoring = restoring, restore_seed = seed), "DEF")
    restoring(TRUE)
    expect_equal(restore_default_value("a", "DEF", restoring = restoring, restore_seed = seed), "DEF")
    seed[["a"]] <- "SAVED"
    expect_equal(restore_default_value("a", "DEF", restoring = restoring, restore_seed = seed), "SAVED")
  })
})

test_that("restore_default_value valid= intersects, falling back when all stale", {
  restoring <- shiny::reactiveVal(TRUE)
  seed <- shiny::reactiveValues()
  shiny::isolate({
    seed[["b"]] <- c("x", "y", "z")
    expect_equal(
      restore_default_value("b", c("x", "y"), valid = c("x", "y"),
                            restoring = restoring, restore_seed = seed),
      c("x", "y"))
    seed[["c"]] <- c("old1", "old2")   # none valid for the current factor
    expect_equal(
      restore_default_value("c", c("n1", "n2"), valid = c("n1", "n2"),
                            restoring = restoring, restore_seed = seed),
      c("n1", "n2"))
  })
})
