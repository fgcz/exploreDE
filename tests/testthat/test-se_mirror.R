# -----------------------------------------------------------------------------
# build_explore_h5ad(mirror = se) stores the whole SummarizedExperiment and
# se_from_h5ad() gives it back identical(): Every assay, colData and rowData
# column, and metadata element at its path, with its class and attributes.
# -----------------------------------------------------------------------------

skip_if_not_installed("anndataR")

mirror_write <- function(se, factors, env = parent.frame()) {
  out <- withr::local_tempfile(fileext = ".h5ad", .local_envir = env)
  build_explore_h5ad(
    x = list(raw = SummarizedExperiment::assay(se, "counts"), normalised = SummarizedExperiment::assay(se, "xNorm")),
    out = out,
    omics = "rnaseq",
    coldata = SummarizedExperiment::colData(se),
    rowdata = SummarizedExperiment::rowData(se),
    factors = factors,
    de = list(trt_vs_ctrl = SummarizedExperiment::rowData(se)),
    mirror = se
  )
  out
}

# Every awkward case met in SUSHI twoGroups results (DESeq2, edgeR, limma):
# NA / NaN in layers, logical matrices, S3 lists with their own `[[`
# (package_version), S4 NULL slots, gene-set lists, captured environments with
# cycles, closures with missing-argument formals, formulas.
awkward_se <- function() {
  G <- 30L
  S <- 4L
  genes <- sprintf("ENSG%011d", seq_len(G))
  samples <- paste0("S", seq_len(S))
  counts <- matrix(rpois(G * S, 20), G, S, dimnames = list(genes, samples))
  storage.mode(counts) <- "double"
  xnorm <- counts / 2
  xnorm[1, 1] <- NA_real_
  xnorm[2, 2] <- NaN
  present <- counts > 15
  present[3, 3] <- NA
  cd <- S4Vectors::DataFrame(
    `Condition [Factor]` = c("ctrl", "ctrl", "trt", "trt"),
    `Read Count` = c(10L, 12L, 11L, 9L),
    sf = c(1, 1.1, 0.9, 1.05),
    row.names = samples,
    check.names = FALSE
  )
  rd <- S4Vectors::DataFrame(
    gene_id = genes,
    gene_name = paste0("G", seq_len(G)),
    isPresentProbe = rep(c(TRUE, FALSE), length.out = G),
    log2Ratio = seq(-2, 2, length.out = G),
    pValue = seq(0.001, 1, length.out = G),
    fdr = seq(0.01, 1, length.out = G),
    usedInTest = rep(c(TRUE, TRUE, FALSE), length.out = G),
    row.names = genes,
    check.names = FALSE
  )
  # A function frame like the ones ezRun and DESeq2 leave behind; its parent
  # is fixed, so the comparison does not wander into testthat's environments.
  frame <- new.env(parent = baseenv())
  local(
    {
      sampleGroup <- "trt"
      dispersion <- function(q, scale) coefs[1] + coefs[2] / q
      coefs <- c(asymptDisp = 0.1, extraPois = 4)
    },
    envir = frame
  )
  frame$design <- stats::as.formula("~grouping", env = frame) # a cycle: frame -> formula -> frame
  gene_sets <- stats::setNames(lapply(seq_len(25), function(i) genes[seq_len(i %% 7)]), sprintf("GO:%07d", 1:25))
  md <- list(
    isLog = FALSE,
    summary = c(Name = "trt--over--ctrl", Normalization = "none"),
    param = list(
      cores = 4L,
      grouping = stats::setNames(c("ctrl", "ctrl", "trt", "trt"), samples),
      seqNames = character(0),
      enrichrDatabase = NULL,
      ratio = 1 / 3
    ),
    nativeResult = list(
      baseMean = runif(G),
      sf = stats::setNames(c(1, 1.1, 0.9, 1.05), samples),
      coefficients = matrix(rnorm(G * 2), G, 2, dimnames = list(NULL, c("Intercept", "trt"))),
      design = matrix(c(1, 1, 1, 1, 0, 0, 1, 1), S, 2),
      fit = frame$dispersion,
      formula = frame$design,
      version = packageVersion("stats"),
      when = as.Date("2026-11-27"),
      level = factor(c("b", "a", NA), levels = c("b", "a"))
    ),
    enrichInput = list(
      table = data.frame(term = c("x", "y"), n = c(1L, NA), row.names = c("r1", "r2")),
      empty = list()
    ),
    geneSets = gene_sets,
    s4 = S4Vectors::DataFrame(a = 1:3)
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, presentFlag = present, xNorm = xnorm),
    colData = cd,
    rowData = rd,
    metadata = md
  )
}

test_that("the demo SUSHI result comes back identical()", {
  se <- readRDS(system.file("examples", "genomics", "deResult.rds", package = "exploreDE"))
  param <- S4Vectors::metadata(se)$param
  out <- withr::local_tempfile(fileext = ".h5ad")
  md <- S4Vectors::metadata(se)
  build_explore_h5ad(
    x = list(raw = SummarizedExperiment::assays(se)$counts, normalised = SummarizedExperiment::assays(se)$xNorm),
    out = out,
    omics = "rnaseq",
    feature_level = param$featureLevel,
    ref_build = param$refBuild,
    coldata = SummarizedExperiment::colData(se),
    rowdata = SummarizedExperiment::rowData(se),
    factors = paste(param$groupingName, "[Factor]"),
    de = stats::setNames(list(SummarizedExperiment::rowData(se)), param$comparison),
    pathways = stats::setNames(list(list(ora = md$enrichResult, gsea = md$enrichResultGSEA)), param$comparison),
    enrichment_fdr = list(ora = param$fdrThreshORA, gsea = param$fdrThreshGSEA),
    mirror = se
  )
  expect_identical(se_from_h5ad(out), se)
  # exploreDE still reads it as before
  eds <- new_explore_dataset(out)
  expect_identical(eds$contrasts(), new_explore_dataset(se)$contrasts())
})

test_that("awkward values come back identical(), environments by contents", {
  set.seed(1)
  se <- awkward_se()
  back <- se_from_h5ad(mirror_write(se, "Condition [Factor]"))
  # Environments cannot keep their address through any copy; everything else
  # is identical().
  expect_true(isTRUE(.se_mirror_compare(back, se)))
  md <- S4Vectors::metadata(se)
  mb <- S4Vectors::metadata(back)
  for (k in setdiff(names(md), "nativeResult")) {
    expect_identical(mb[[k]], md[[k]], info = k)
  }
  for (k in setdiff(names(md$nativeResult), c("fit", "formula"))) {
    expect_identical(mb$nativeResult[[k]], md$nativeResult[[k]], info = k)
  }
  expect_identical(SummarizedExperiment::assays(back), SummarizedExperiment::assays(se))
  expect_identical(SummarizedExperiment::colData(back), SummarizedExperiment::colData(se))
  expect_identical(SummarizedExperiment::rowData(back), SummarizedExperiment::rowData(se))
  # the rebuilt function works, from its rebuilt environment
  expect_identical(mb$nativeResult$fit(2), md$nativeResult$fit(2))
  expect_identical(environment(mb$nativeResult$formula)$design, mb$nativeResult$formula)
})

test_that("values are placed by shape, and the DE columns live in the DE table only", {
  set.seed(1)
  se <- awkward_se()
  ad <- anndataR::read_h5ad(mirror_write(se, "Condition [Factor]"), as = "InMemoryAnnData")
  expect_true("assays.presentFlag" %in% names(ad$layers))
  expect_true("nativeResult.coefficients" %in% names(ad$varm))
  expect_true("nativeResult.design" %in% names(ad$obsm))
  expect_false(any(c("label_log2ratio", "label_pvalue", "label_fdr", "label_usedintest") %in% colnames(ad$var)))
  expect_true("label_ispresentprobe" %in% colnames(ad$var))
  m <- ad$uns$se$manifest
  expect_identical(m$how[m$path == "se@assays@data@listData$counts"], "atomic")
  expect_match(m$location[m$path == "se@assays@data@listData$counts"], "^layers\\|counts$")
  expect_identical(m$how[m$path == "se@metadata$geneSets"], "long")
})

test_that("a file without a mirror says so", {
  set.seed(1)
  se <- awkward_se()
  out <- withr::local_tempfile(fileext = ".h5ad")
  build_explore_h5ad(
    x = list(raw = SummarizedExperiment::assay(se, "counts")),
    out = out,
    omics = "rnaseq",
    coldata = SummarizedExperiment::colData(se),
    factors = "Condition [Factor]"
  )
  expect_error(se_from_h5ad(out), "no SummarizedExperiment")
})
