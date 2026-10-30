# -----------------------------------------------------------------------------
# Helpers for the ExploreDatasetH5AD / build_explore_h5ad tests.
#
#   make_toy_dds()          — a small, self-contained DESeqDataSet (no gstore),
#                             so the writer+reader+validator loop always runs.
#   se_to_canonical_h5ad()  — convert a SUSHI RNA-seq SE to a canonical .h5ad,
#                             so the h5ad backend can be diffed accessor-for-
#                             accessor against the SE backend (equivalence).
#
# These mirror plans/archive/anndata_migration/poc_h5ad_roundtrip.R, and fold in the
# three anndataR gotchas (logical layers -> integer; grep [Factor] on the S4
# colData; write_h5ad mode = "w").
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# A tiny two-condition RNA-seq DESeqDataSet with unique gene symbols.
make_toy_dds <- function(n_genes = 300, n_per_group = 3, seed = 42) {
  set.seed(seed)
  n <- n_per_group * 2L
  mu <- matrix(stats::rnbinom(n_genes * n, mu = 100, size = 5), nrow = n_genes)
  # inject a real signal into the first 30 genes for the treated group
  treated <- rep(c(FALSE, TRUE), each = n_per_group)
  mu[1:30, treated] <- mu[1:30, treated] * 4L
  rownames(mu) <- sprintf("Gene%03d", seq_len(n_genes))
  colnames(mu) <- sprintf("S%02d", seq_len(n))
  cd <- S4Vectors::DataFrame(
    condition = factor(rep(c("control", "treated"), each = n_per_group), levels = c("control", "treated")),
    row.names = colnames(mu)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(mu, cd, design = ~condition)
  DESeq2::DESeq(dds, quiet = TRUE)
}

# Loose, disparate pieces (as a non-DESeq2 user would have them): a counts
# matrix, a normalised matrix, colData, a DE table with NON-canonical column
# names, and rowData carrying a per-feature pathway annotation.
make_toy_pieces <- function() {
  dds <- make_toy_dds()
  cnt <- DESeq2::counts(dds)
  nrm <- DESeq2::counts(dds, normalized = TRUE)
  res <- as.data.frame(DESeq2::results(dds, contrast = c("condition", "treated", "control")))
  de_df <- data.frame(
    row.names = rownames(res), # deliberately odd names
    LFC = res$log2FoldChange,
    P = res$pvalue,
    Q = res$padj
  )
  n <- nrow(cnt)
  rowdata <- data.frame(
    row.names = rownames(cnt),
    gene_name = rownames(cnt),
    `GO BP` = ifelse(seq_len(n) <= 20L, "GO:0006955", NA_character_),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  list(
    counts = cnt,
    normalised = nrm,
    coldata = as.data.frame(SummarizedExperiment::colData(dds), check.names = FALSE),
    de_df = de_df,
    rowdata = rowdata,
    contrast = "treated_vs_control"
  )
}

# Toy proteomics pieces (intensities + colData + rowData) for a given flavour.
# Protein ids are duplicated (2 features per protein) so that protein-level
# uniquifies but phospho/peptide disambiguate via site/peptide_Id.
make_toy_prot <- function(flavour = c("protein", "phospho", "peptide")) {
  flavour <- match.arg(flavour)
  n <- 40L
  ns <- 6L
  m <- matrix(seq_len(n * ns) + 0.5, nrow = n, ncol = ns)
  rownames(m) <- sprintf("F%03d", seq_len(n))
  colnames(m) <- sprintf("S%02d", seq_len(ns))
  rd <- data.frame(
    row.names = rownames(m),
    protein_Id = sprintf("P%03d", rep(seq_len(n / 2L), each = 2L)),
    gene_name = sprintf("Gene%03d", rep(seq_len(n / 2L), each = 2L)),
    stringsAsFactors = FALSE
  )
  if (flavour == "phospho") {
    rd$site <- sprintf("S%d", seq_len(n))
  }
  if (flavour == "peptide") {
    rd$peptide_Id <- sprintf("PEP%d", seq_len(n))
  }
  coldata <- data.frame(row.names = colnames(m), group = rep(c("A", "B"), each = ns / 2L))
  list(intensities = m, coldata = coldata, rowdata = rd)
}

# Faithful SUSHI RNA-seq SE -> canonical exploreDE/1.0 .h5ad.
se_to_canonical_h5ad <- function(se, out) {
  md <- S4Vectors::metadata(se)
  param <- md$param %||% list()

  # grep on the S4 colData colnames — as.data.frame() mangles "[Factor]".
  cn <- colnames(SummarizedExperiment::colData(se))
  fo <- grep("\\[Factor\\]", cn, value = TRUE)
  fn <- paste0("factor_", gsub("[^a-z0-9]+", "_", tolower(sub(" \\[Factor\\]", "", fo))))
  obs <- data.frame(row.names = colnames(se))
  for (i in seq_along(fo)) {
    obs[[fn[i]]] <- factor(as.character(SummarizedExperiment::colData(se)[[fo[i]]]))
  }
  obs$sample_label <- rownames(obs)
  if (!is.null(param$refBuild)) {
    obs$organism <- sub("/.*", "", param$refBuild)
  }

  rd <- as.data.frame(SummarizedExperiment::rowData(se), check.names = FALSE)
  fid <- scuttle::uniquifyFeatureNames(rd$gene_id, rd$gene_name)
  fid[grep("^_", fid)] <- gsub("_", "", fid[grep("^_", fid)])
  var <- data.frame(
    row.names = fid,
    description = sub(" \\[.*", "", as.character(rd$description %||% rd$gene_name)),
    label_gene_id = rd$gene_id,
    label_gene_symbol = rd$gene_name,
    stringsAsFactors = FALSE
  )
  go_map <- c(
    "GO BP" = "label_go_bp",
    "GO.BP" = "label_go_bp",
    "GO MF" = "label_go_mf",
    "GO.MF" = "label_go_mf",
    "GO CC" = "label_go_cc",
    "GO.CC" = "label_go_cc"
  )
  for (g in names(go_map)) {
    if (g %in% colnames(rd)) var[[go_map[[g]]]] <- rd[[g]]
  }
  if ("featWidth" %in% colnames(rd)) {
    var$label_featwidth <- rd$featWidth
  } # -> computed FPKM/TPM

  amap <- c(
    counts = "counts",
    xNorm = "normalised",
    signal = "normalised",
    presentFlag = "present_flag",
    isPresent = "present_flag"
  )
  layers <- list()
  for (a in SummarizedExperiment::assayNames(se)) {
    key <- amap[[a]] %||% a
    m <- as.matrix(SummarizedExperiment::assay(se, a))
    rownames(m) <- fid
    if (is.logical(m)) {
      storage.mode(m) <- "integer"
    } # anndataR mangles logical
    layers[[key]] <- t(m)
  }
  if ("normalised" %in% names(layers)) {
    layers$log2 <- log2(layers$normalised + (param$backgroundExpression %||% 10))
  }
  X <- layers$counts

  slot <- gsub("[^a-z0-9]+", "_", tolower(param$comparison %||% "de"))
  varm <- list()
  if (all(c("log2Ratio", "pValue", "fdr") %in% colnames(rd))) {
    varm[[paste0("de_", slot)]] <- data.frame(
      row.names = fid,
      effect_log2fc = as.numeric(rd$log2Ratio),
      score_pvalue = as.numeric(rd$pValue),
      score_fdr = as.numeric(rd$fdr)
    )
  }
  uns_de <- setNames(
    list(list(
      assay_used = "counts",
      factor_used = fn[1],
      test_group = as.character(param$sampleGroup %||% NA),
      reference_group = as.character(param$refGroup %||% NA),
      model = md$method %||% "DESeq2"
    )),
    slot
  )

  adata <- anndataR::AnnData(
    X = X,
    obs = obs,
    var = var,
    layers = layers,
    varm = varm,
    uns = list(
      schema_version = "exploreDE/1.0",
      omics_type = "rnaseq_bulk",
      hierarchy = list(levels = "gene", unit = "gene"),
      design = list(factors = fn, grouping = fn[1]),
      de = uns_de,
      enrichment = list(),
      enrichment_native = list(),
      reference_genome = as.character(param$refBuild %||% NA)
    )
  )
  anndataR::write_h5ad(adata, out, mode = "w")
  out
}

# Real (tiny) clusterProfiler objects with everything the support block reads:
# significant terms, geneID/core_enrichment, and for GSEA gene sets (including
# one term that is NOT in the result, which must not be stored) + a ranked list.
.mini_cp <- function(cls, n) {
  ids <- sprintf("GO:%04d", seq_len(n))
  res <- data.frame(
    ID = ids,
    Description = paste("term", seq_len(n)),
    pvalue = 0.001,
    p.adjust = 0.01,
    stringsAsFactors = FALSE,
    row.names = ids
  )
  if (cls == "enrichResult") {
    res$geneID <- "G1/G2"
    res$Count <- 2L
    return(methods::new(cls, result = res))
  }
  res$NES <- 1.5
  res$enrichmentScore <- 0.5
  res$setSize <- 3L
  res$core_enrichment <- "G1/G2"
  sets <- c(setNames(rep(list(c("G1", "G2", "G3")), n), ids), list("GO:9999" = c("G4", "G5")))
  methods::new(
    cls,
    result = res,
    geneSets = sets,
    geneList = c(G1 = 3, G2 = 2, G3 = 1, G4 = -1, G5 = -2),
    params = list(exponent = 1)
  )
}
