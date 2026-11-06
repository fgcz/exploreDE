# Trim the demo dataset inst/examples/genomics/deResult.rds to the parts
# exploreDE reads, so the file stays under Bioconductor's 5 MB limit.
#
# The input is a SUSHI DESeq2 result (RNA-seq, diff--over--undiff, human
# GENCODE GRCh38 release 42). Kept unchanged: all assays, colData, rowData and
# the result tables of the ORA and GSEA objects. The output is a plain
# SummarizedExperiment. Removed:
#   - rowRanges (unused; its Seqinfo needs GenomeInfoDb)
#   - metadata$nativeResult (the DESeqDataSet; the app never reads it)
#   - metadata$param$ezRef (an ezRun object; param$refBuild holds the same
#     reference build, and exploreDE does not depend on ezRun)
#   - metadata$enrichInput fields other than seqAnno, log2Ratio and selections
#     (the ones ExploreDatasetSE$enrichment_input() returns), and seqAnno
#     columns other than gene_id, gene_name, pValue and fdr (the ones the ORA
#     and GSEA tabs read; seqAnno repeats rowData, which keeps every column)
#   - gene sets of GO terms that are not in an ORA/GSEA result table; the
#     network, ridge and running-score plots read only those of the stored terms
#
# Usage: Rscript inst/scripts/trim_demo_deResult.R <in.rds> <out.rds>

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DOSE)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
se <- readRDS(args[1])

keep_result_gene_sets <- function(x) {
  if (is(x, "enrichResult") || is(x, "gseaResult")) {
    x@geneSets <- x@geneSets[intersect(x@result$ID, names(x@geneSets))]
    return(x)
  }
  if (is.list(x)) {
    return(lapply(x, keep_result_gene_sets))
  }
  x
}

md <- metadata(se)
md$nativeResult <- NULL
md$param$ezRef <- NULL
md$enrichInput <- md$enrichInput[c("seqAnno", "log2Ratio", "selections")]
md$enrichInput$seqAnno <- md$enrichInput$seqAnno[, c("gene_id", "gene_name", "pValue", "fdr")]
md$enrichResult <- keep_result_gene_sets(md$enrichResult)
md$enrichResultGSEA <- keep_result_gene_sets(md$enrichResultGSEA)

# A plain SummarizedExperiment: the RangedSummarizedExperiment's rowRanges
# (unused by exploreDE) needs GenomeInfoDb, which exploreDE does not depend on.
se <- SummarizedExperiment(
  assays = assays(se),
  colData = colData(se),
  rowData = rowData(se),
  metadata = md
)

saveRDS(se, args[2], compress = "xz")
cat(sprintf("%s: %.1f MB\n", args[2], file.size(args[2]) / 1e6))
