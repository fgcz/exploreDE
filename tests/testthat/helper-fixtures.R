# -----------------------------------------------------------------------------
# Canonical test fixtures
#
# Loads SE objects and caches them for the life of the R session. Each test
# references a fixture via the `fx_*()` accessors below. The genomics fixtures
# are the shipped demo and SUSHI outputs built from it, so they run
# everywhere; the proteomics ones are read from pStore (URLs) and the test
# skips when the URL is unreachable.
#
# Bioconductor generics are used unqualified inside R/ (e.g. feature_ids.R)
# (metadata, assay, rowData, colData, ...). The installed app
# attaches these via library() at startup; tests need to mirror that, since
# NAMESPACE does not currently importFrom them.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(shiny) # isolate() / reactiveValues() used unqualified
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggrepel) # ggrepel::geom_*_repel used unqualified
})

.fixture_cache <- new.env(parent = emptyenv())

fixture_specs <- list(
  # The package's demo dataset: this SUSHI comparison
  # (p3009/o5638_DESeq2_diff--over--undiff), trimmed by
  # inst/scripts/trim_demo_deResult.R. Shipped, so its tests run everywhere.
  rnaseq_gene = list(
    kind = "genomics",
    package_file = "examples/genomics/deResult.rds"
  ),
  # SUSHI outputs built from the demo by the build_*() functions below. Their
  # layout follows the real gstore files (read on 2026-11-13): the p3009 o5638
  # DESeq2 directory and Count_QC, p34429 o35733 DESeq2 at isoform level, and
  # p36911 ExploreMageckCounts.
  rnaseq_gene_dir = list(kind = "genomics", build = function() build_rnaseq_gene_dir()),
  rnaseq_isoform = list(kind = "genomics", build = function() build_rnaseq_isoform()),
  countqc = list(kind = "genomics", build = function() build_countqc()),
  crispr_countqc = list(kind = "genomics", build = function() build_crispr_countqc()),
  proteomics_protein = list(
    kind = "proteomics",
    url = "https://fgcz-ms.uzh.ch/public/pStore/p23078/bfabric/Proteomics/exploreDE/2026/2026-03/2026-03-10/workunit_342063/3106962.rds"
  ),
  proteomics_peptide = list(
    kind = "proteomics",
    url = "https://fgcz-ms.uzh.ch/public/pStore/p23078/bfabric/Proteomics/exploreDE/2026/2026-03/2026-03-10/workunit_342065/3106968.rds"
  ),
  proteomics_phospho = list(
    kind = "proteomics",
    url = "https://fgcz-ms.uzh.ch/public/pStore/p39819/bfabric/Proteomics/exploreDE/2025/2025-10/2025-10-22/workunit_335485/2997996.rds"
  ),
  # SAINTexpress AP-MS: metadata$param NULL, DE cols log2_EFCs/SaintScore/BFDR,
  # no raw p-value. Exercises the transitional SE-backend SAINT support.
  proteomics_saint = list(
    kind = "proteomics",
    url = "https://fgcz-ms.uzh.ch/public/pStore/p42111/bfabric/Proteomics/exploreDE/2026/2026-07/2026-07-28/workunit_349205/3206509.rds"
  )
)

load_fixture <- function(key) {
  spec <- fixture_specs[[key]]
  if (is.null(spec)) {
    stop("Unknown fixture key: ", key)
  }

  cached <- .fixture_cache[[key]]
  if (!is.null(cached)) {
    return(cached)
  }

  result <- if (!is.null(spec$build)) {
    c(spec$build(), kind = spec$kind, key = key)
  } else if (spec$kind == "genomics") {
    path <- system.file(spec$package_file, package = "exploreDE", mustWork = TRUE)
    list(se = readRDS(path), dataDir = dirname(path), kind = "genomics", key = key)
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

fx_rnaseq_gene <- function() load_fixture("rnaseq_gene")
fx_rnaseq_gene_dir <- function() load_fixture("rnaseq_gene_dir")
fx_rnaseq_isoform <- function() load_fixture("rnaseq_isoform")
fx_countqc <- function() load_fixture("countqc")
fx_crispr_countqc <- function() load_fixture("crispr_countqc")
fx_proteomics_protein <- function() load_fixture("proteomics_protein")
fx_proteomics_peptide <- function() load_fixture("proteomics_peptide")
fx_proteomics_phospho <- function() load_fixture("proteomics_phospho")
fx_proteomics_saint <- function() load_fixture("proteomics_saint")

# ---- SUSHI outputs built from the demo --------------------------------------

# The DESeq2 result directory with SUSHI's ../input_dataset.tsv next to it
# (Name plus the sample table; the demo's Count [File] paths name Kallisto).
build_rnaseq_gene_dir <- function() {
  se <- load_fixture("rnaseq_gene")$se
  root <- tempfile("sushi_deseq2_")
  data_dir <- file.path(root, "diff--over--undiff")
  dir.create(data_dir, recursive = TRUE)
  cd <- data.frame(SummarizedExperiment::colData(se), check.names = FALSE)
  cd$sf <- NULL
  utils::write.table(
    data.frame(Name = colnames(se), cd, check.names = FALSE),
    file.path(root, "input_dataset.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  list(se = se, dataDir = data_dir)
}

# DESeq2 result at isoform level: 1 to 3 transcripts per gene, each with its
# gene's statistics, keyed on Ensembl transcript ids; featureLevel "isoform"
# in metadata and param, as SUSHI writes it.
build_rnaseq_isoform <- function() {
  se <- load_fixture("rnaseq_gene")$se[1:300, ]
  se <- se[rep(seq_len(nrow(se)), rep_len(1:3, nrow(se))), ]
  tx <- sprintf("ENST%011d", seq_len(nrow(se)))
  SummarizedExperiment::rowData(se)$transcript_id <- tx
  rownames(se) <- tx
  md <- S4Vectors::metadata(se)
  md$featureLevel <- "isoform"
  md$param$featureLevel <- "isoform"
  S4Vectors::metadata(se) <- md
  list(se = se)
}

# Count QC of the same samples: no DE columns and no comparison; assays
# counts, presentFlag and signal (SUSHI also writes log2signal and
# log2signalCentered, which exploreDE does not read); rowData flag isValid.
build_countqc <- function() {
  se <- load_fixture("rnaseq_gene")$se
  SummarizedExperiment::assays(se) <- list(
    counts = SummarizedExperiment::assay(se, "counts"),
    presentFlag = SummarizedExperiment::assay(se, "presentFlag"),
    signal = SummarizedExperiment::assay(se, "xNorm")
  )
  rd <- SummarizedExperiment::rowData(se)
  rd <- rd[, setdiff(colnames(rd), c("isPresentProbe", "log2Ratio", "pValue", "fdr", "usedInTest"))]
  rd$isValid <- rep(TRUE, nrow(rd))
  SummarizedExperiment::rowData(se) <- rd
  se$sf <- NULL
  se$conds <- se$`Condition [Factor]`
  md <- S4Vectors::metadata(se)[c("isLog", "featureLevel", "type", "countName", "param")]
  md$param[c("comparison", "sampleGroup", "refGroup", "grouping")] <- ""
  md$param$normMethod <- "logMean"
  S4Vectors::metadata(se) <- md
  list(se = se)
}

# CRISPR guide counts as ExploreMageckCounts writes them: 3 guides for each of
# 44 genes (132 rows, as in p36911), ids GENE_<n>_<protospacer>, and gene_id,
# transcript_id and gene_name all holding the guide id.
build_crispr_countqc <- function() {
  se <- load_fixture("countqc")$se[1:44, ]
  genes <- SummarizedExperiment::rowData(se)$gene_name
  se <- se[rep(seq_len(nrow(se)), each = 3), ]
  # 20-nt protospacer: the guide's index in base 4 over ACGT
  protospacer <- vapply(
    seq_len(nrow(se)),
    function(i) paste(c("A", "C", "G", "T")[(i %/% 4^(19:0)) %% 4 + 1], collapse = ""),
    character(1)
  )
  guides <- paste(rep(genes, each = 3), rep(1:3, length(genes)), protospacer, sep = "_")
  rownames(se) <- guides
  for (col in c("gene_id", "transcript_id", "gene_name")) {
    SummarizedExperiment::rowData(se)[[col]] <- guides
  }
  list(se = se)
}
