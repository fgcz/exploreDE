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

test_that("contrastSelected is kept but its per-tab mirrors are dropped", {
  x <- list(contrastSelected = "A", contrastSelected_volcano = "A",
            contrastSelected_heatmap = "A", contrastSelected_ora = "A")
  kept <- names(save_state_filter_inputs(x))
  expect_equal(kept, "contrastSelected")
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

# ---- JSON session files (exploreDE-session/2) -------------------------------

# A bundle shaped like the app's build_state_bundle(), with the value types a
# real session holds (548 inputs from a demo session were character, integer,
# numeric, logical, NULL, empty, and two named lists for bucket_list groups).
state_bundle <- function(contrasts = "diff_vs_undiff") {
  list(
    meta = list(app = "exploreDE", version = 1L,
                saved_at = as.POSIXct("2026-10-29 11:51:23"), dataUrl = "",
                fingerprint = list(dataDir = "/data/genomics", dataType = "RNASeq", omics = "genomics",
                                   nSamples = 9L, nFeatures = 21057L, factors = "Condition",
                                   contrasts = contrasts)),
    inputs = list(heatmapColourRed = "#8B1A1A", featureCalcLFC = 0.5, boxplotNCol = 3L,
                  pcaEllipses = TRUE, emptySelect = character(0), sidebarItemExpanded = NULL,
                  multiPick = c("a", "b"),
                  bucket_list_group = list(boxKeepBucket = c("diff", "undiff"), boxExcludeBucket = character(0)),
                  someButton = structure(2L, class = c("shinyActionButtonValue", "integer"))),
    rv = list(featureBucket = c("GAPDH", "ACTB"), goHeatmap = list(genes = NULL, name = NULL, term = NULL),
              oraSelected = c(1L, 3L), gseaSelected = NULL, gseaType = "BP"),
    multi = list(
      datasets = list(list(name = "KO over WT", source = "link", link = "p1000/DESeq2_x/KO--over--WT"),
                      list(name = "B_vs_A", source = "contrast", contrast = "B_vs_A")),
      inputs = list(multiBoxplotNCol = 2L, multiGroupColourKO_over_WT = "#1B9E77")))
}

write_json_session <- function(bundle, env = parent.frame()) {
  f <- withr::local_tempfile(fileext = ".json", .local_envir = env)
  writeLines(save_state_to_json(bundle), f)
  f
}

test_that("a session round-trips through JSON with the same values", {
  b <- state_bundle()
  r <- save_state_from_json(write_json_session(b))
  expect_identical(r$meta$format, "exploreDE-session/2")
  for (k in setdiff(names(b$inputs), c("someButton", "sidebarItemExpanded")))
    expect_equal(r$inputs[[k]], b$inputs[[k]], info = k)
  expect_null(r$inputs$sidebarItemExpanded)
  expect_identical(r$rv$featureBucket, b$rv$featureBucket)
  expect_identical(r$rv$oraSelected, b$rv$oraSelected)
  expect_identical(r$rv$gseaType, "BP")
  expect_equal(as.numeric(r$meta$saved_at), as.numeric(b$meta$saved_at))
})

test_that("the fingerprint still matches after JSON, also with no contrasts", {
  for (con in list("diff_vs_undiff", character())) {
    b <- state_bundle(contrasts = con)
    r <- save_state_from_json(write_json_session(b))
    expect_true(save_state_fp_match(r$meta$fingerprint, b$meta$fingerprint))
  }
})

test_that("action-button values are not written", {
  json <- save_state_to_json(state_bundle())
  expect_false(grepl("someButton", json, fixed = TRUE))
  expect_true(grepl('"format": ["exploreDE-session/2"]', json, fixed = TRUE))
})

test_that("the Multi DE block round-trips: datasets in order, and their widgets", {
  r <- save_state_from_json(write_json_session(state_bundle()))
  expect_equal(vapply(r$multi$datasets, `[[`, "", "name"), c("KO over WT", "B_vs_A"))
  expect_identical(r$multi$datasets[[1]], list(name = "KO over WT", source = "link", link = "p1000/DESeq2_x/KO--over--WT"))
  expect_identical(r$multi$datasets[[2]], list(name = "B_vs_A", source = "contrast", contrast = "B_vs_A"))
  expect_identical(r$multi$inputs$multiGroupColourKO_over_WT, "#1B9E77")
})

test_that("Multi DE entries with an unknown source or no reference are dropped", {
  b <- state_bundle()
  b$multi$datasets <- c(b$multi$datasets,
                        list(list(name = "x", source = "file", link = "/etc/passwd"),
                             list(name = "y", source = "link")))
  r <- save_state_from_json(write_json_session(b))
  expect_equal(length(r$multi$datasets), 2L)
})

test_that("the reader refuses what isn't a session", {
  f <- withr::local_tempfile(fileext = ".json")
  bad <- list(
    not_json   = "{not json",
    array      = "[1, 2, 3]",
    other_app  = '{"meta": {"app": ["other"], "format": ["exploreDE-session/2"]}}',
    old_format = '{"meta": {"app": ["exploreDE"], "format": ["exploreDE-session/1"]}}',
    no_meta    = '{"inputs": {}}',
    too_deep   = paste0('{"meta": {"app": ["exploreDE"], "format": ["exploreDE-session/2"]}, "inputs": ',
                        strrep('{"a": ', 8), "1", strrep("}", 8), "}"))
  for (nm in names(bad)) {
    writeLines(bad[[nm]], f)
    expect_error(save_state_from_json(f), info = nm)
  }
})

test_that("the reader refuses files over the size limit", {
  f <- write_json_session(state_bundle())
  expect_error(save_state_from_json(f, max_bytes = 100), "larger than")
})

test_that("code-like content stays text, and R objects cannot be written", {
  b <- state_bundle()
  b$inputs$heatmapColourRed <- 'system("id")'
  r <- save_state_from_json(write_json_session(b))
  expect_identical(r$inputs$heatmapColourRed, 'system("id")')
  e <- new.env()
  expect_error(save_state_to_json(list(meta = list(app = "exploreDE"), inputs = list(x = e))))
})

test_that("a retired .rds session renamed to .json is refused, never deserialised", {
  f <- tempfile(fileext = ".json")
  saveRDS(state_bundle(), f)
  expect_error(save_state_from_json(f), "not valid JSON")
})

# ---- Multi DE helpers ------------------------------------------------------

test_that("save_state_multi_inputs keeps multi widgets and drops the rest", {
  x <- list(multiBoxplotNCol = 2L, multiGroupColourA = "#000000", multiContrastLink1 = "p1/x",
            multiImportFile = list(), multiDatasetsCheckbox = "/srv/x", multiContrastOrder = c("a", "b"),
            multiDlColourTemplate = 1L, multiSelectedDatasetsTable_rows_selected = 1L,
            multiApplyOrder = structure(1L, class = c("shinyActionButtonValue", "integer")),
            boxplotFactor1 = "Condition")
  expect_setequal(names(save_state_multi_inputs(x)),
                  c("multiBoxplotNCol", "multiGroupColourA", "multiContrastLink1"))
})

test_that("save_state_multi_datasets lists sources in display order", {
  res <- list(`My KO` = list(source = list(kind = "link", link = "p1/D/KO")),
              B_vs_A = list(source = list(kind = "contrast", contrast = "B_vs_A")),
              nosource = list(seqAnno = data.frame()))
  out <- save_state_multi_datasets(res)
  expect_equal(vapply(out, `[[`, "", "name"), c("My KO", "B_vs_A"))
  expect_equal(out[[1]]$link, "p1/D/KO")
})
