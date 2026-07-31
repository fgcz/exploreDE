# -----------------------------------------------------------------------------
# Canonical test fixtures
#
# Loads real SE objects from gstore (local paths) or pStore (URLs) and caches
# them for the life of the R session. Each test references a fixture via the
# `fx_*()` accessors below; if the underlying file / URL is unreachable, the
# test skips gracefully (does not error).
#
# Bioconductor / ezRun generics are used unqualified inside R/helperFunctions.R
# (metadata, assay, rowData, colData, getRpkm, getTpm, ...). The installed app
# attaches these via library() at startup; tests need to mirror that, since
# NAMESPACE does not currently importFrom them.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(shiny)               # isolate() / reactiveValues() used unqualified
  library(SummarizedExperiment)
  library(DESeq2)
  library(ezRun)
  library(ggrepel)             # ggrepel::geom_*_repel used unqualified
})

.fixture_cache <- new.env(parent = emptyenv())

fixture_specs <- list(
  rnaseq_gene = list(
    kind    = "genomics",
    dataDir = "/srv/gstore/projects/p3009/o5638_DESeq2_diff--over--undiff_2024-11-20--12-29-40/diff--over--undiff",
    file    = "deResult.rds"
  ),
  rnaseq_isoform = list(
    kind    = "genomics",
    dataDir = "/srv/gstore/projects/p34429/o35733_DESeq2_Palbo_day6--over--Palbo_day2_2024-08-22--11-21-44/Palbo_day6--over--Palbo_day2",
    file    = "deResult.rds"
  ),
  countqc = list(
    kind    = "genomics",
    dataDir = "/srv/gstore/projects/p3009/o5638_CountQC_2025-11-14--15-14-04/Count_QC",
    file    = "rawData_final.rds"
  ),
  crispr_countqc = list(
    kind    = "genomics",
    dataDir = "/srv/gstore/projects/p36911/o39664_o37783_o39124_o38710_o41353_ExploreMageckCounts_2026-04-09--15-07-24/ExploreMageckCounts",
    file    = "rawData_final.rds"
  ),
  proteomics_protein = list(
    kind = "proteomics",
    url  = "https://fgcz-ms.uzh.ch/public/pStore/p23078/bfabric/Proteomics/exploreDE/2026/2026-03/2026-03-10/workunit_342063/3106962.rds"
  ),
  proteomics_peptide = list(
    kind = "proteomics",
    url  = "https://fgcz-ms.uzh.ch/public/pStore/p23078/bfabric/Proteomics/exploreDE/2026/2026-03/2026-03-10/workunit_342065/3106968.rds"
  ),
  proteomics_phospho = list(
    kind = "proteomics",
    url  = "https://fgcz-ms.uzh.ch/public/pStore/p39819/bfabric/Proteomics/exploreDE/2025/2025-10/2025-10-22/workunit_335485/2997996.rds"
  ),
  # SAINTexpress AP-MS: metadata$param NULL, DE cols log2_EFCs/SaintScore/BFDR,
  # no raw p-value. Exercises the transitional SE-backend SAINT support.
  proteomics_saint = list(
    kind = "proteomics",
    url  = "https://fgcz-ms.uzh.ch/public/pStore/p42111/bfabric/Proteomics/exploreDE/2026/2026-07/2026-07-28/workunit_349205/3206509.rds"
  )
)

load_fixture <- function(key) {
  spec <- fixture_specs[[key]]
  if (is.null(spec)) stop("Unknown fixture key: ", key)

  cached <- .fixture_cache[[key]]
  if (!is.null(cached)) return(cached)

  result <- if (spec$kind == "genomics") {
    path <- file.path(spec$dataDir, spec$file)
    if (!file.exists(path)) {
      testthat::skip(sprintf("Fixture '%s' not found at %s", key, path))
    }
    list(se = readRDS(path), dataDir = spec$dataDir, kind = "genomics", key = key)
  } else if (spec$kind == "proteomics") {
    se <- tryCatch(
      readRDS(url(spec$url)),
      error = function(e) {
        testthat::skip(sprintf("Fixture '%s' unreachable: %s", key, conditionMessage(e)))
      }
    )
    list(se = se, url = spec$url, kind = "proteomics", key = key)
  } else {
    stop("Unknown fixture kind: ", spec$kind)
  }

  .fixture_cache[[key]] <- result
  result
}

fx_rnaseq_gene        <- function() load_fixture("rnaseq_gene")
fx_rnaseq_isoform     <- function() load_fixture("rnaseq_isoform")
fx_countqc            <- function() load_fixture("countqc")
fx_crispr_countqc     <- function() load_fixture("crispr_countqc")
fx_proteomics_protein <- function() load_fixture("proteomics_protein")
fx_proteomics_peptide <- function() load_fixture("proteomics_peptide")
fx_proteomics_phospho <- function() load_fixture("proteomics_phospho")
fx_proteomics_saint   <- function() load_fixture("proteomics_saint")
