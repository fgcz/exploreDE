---
name: exploreDE-byod
description: Write an R script that turns a user's own omics data (counts, optionally a differential-expression table and pathway/enrichment results) into an exploreDE .h5ad file with build_explore_h5ad(). Use when someone wants to load their own RNA-seq, proteomics or CRISPR screen data into the exploreDE Shiny app.
---

# Building an exploreDE file from a user's data

The user may know nothing about exploreDE's internals. Your job is a short,
correct R script that calls `build_explore_h5ad()` on *their* data, then
checks the result. Read this whole file before writing code.

## Rules

- exploreDE reads exactly one format: An `.h5ad` in the `exploreDE/1.0`
  layout. In R, always write it with `build_explore_h5ad()`; never hand-write
  `.h5ad` structure or invent another format.
- If the user's pipeline isn't in R (Python, an exporter app), don't write R
  for them: Point them to the format spec,
  `system.file("spec/exploreDE-h5ad-1.0.md", package = "exploreDE")`, which
  has every required key and a minimal Python writer.
- Some files need no script at all: APB Export's exploreDE output and the
  prolfquapp DEA's `AnnData.h5ad` (prolfquapp 2.10.x) open as they are. Check
  with `new_explore_dataset()` (see "Verify") before writing anything.
- Look at the user's actual data first (`dim()`, `colnames()`,
  `head(rownames())`, `class()`). Never guess column names, and never write
  a mapping to a column that isn't there.
- Don't rename or pre-process DE or annotation columns. The column
  dictionary recognises common names; renaming defeats it and can hide
  problems.
- Only use arguments the user's data needs. Don't invent placeholder data.
- Finish with the check in "Verify". "No error" is not success: an
  unrecognised DE column comes through as `NA`, not as an error.

## The function

```r
build_explore_h5ad(
  x, out, omics,
  feature_level = NULL, coldata = NULL, rowdata = NULL,
  de = list(), de_columns = NULL, rowdata_columns = NULL,
  pathways = list(), enrichment_padj_max = 0.1, factors = NULL,
  gene_symbol = NULL, organism = NA, ref_build = NA
)
```

Required:

- `x`: counts, **features as rows, samples as columns**. A matrix, a named
  list of matrices (`list(counts = , normalised = )`), or a `DESeqDataSet`.
  If unsure of orientation, the sample axis is the one whose names match
  `coldata`'s rownames.
- `coldata`: one row per sample, rownames equal to `colnames(counts)`.
  Not needed if `x` is a `DESeqDataSet`.
- `omics`: exactly `"rnaseq"`, `"proteomics"`, or `"crispr"`.
- `out`: output path ending in `.h5ad`.

Optional:

- `rowdata`: per-feature annotation (IDs, symbols, description, GO columns).
  Every column is kept.
- `de`: named list, one DE table (data.frame / `DESeqResults`) per contrast,
  e.g. `list(Treat_vs_Ctrl = res)`. Only the recognised statistics are kept
  per contrast (log2 fold change, p-value, FDR, optional extra score, optional
  `usedInTest` flag); other
  DE columns are dropped, so put per-feature information in `rowdata`.
- `pathways`: named by contrast (a subset of `de`'s names), each a list with
  `ora` and/or `gsea`. Values: clusterProfiler `enrichResult`/`gseaResult`,
  plain data.frames, or nested lists of either. Data.frames must use
  clusterProfiler column names (`ID`, `Description`, `p.adjust`,
  `GeneRatio`, `Count`, `geneID`; GSEA adds `enrichmentScore`, `NES`,
  `setSize`, `core_enrichment`). clusterProfiler objects are filtered to
  `p.adjust <= enrichment_padj_max`; data.frames are stored as given.
- `factors`: which `coldata` columns are experimental factors (default: all
  character/factor columns).
- `organism`, `ref_build`: labels, e.g. `"Homo sapiens"`.
- `feature_level`: normally detected (see below); override only if wrong.

A counts-only file is complete and valid. DE and pathways are optional.

## Examples (real, runnable on the data shipped with exploreDE)

Counts only:

```r
library(exploreDE); library(SummarizedExperiment)
se <- readRDS(system.file("examples/genomics/deResult.rds", package = "exploreDE"))
build_explore_h5ad(assay(se, "counts"), out = "counts_only.h5ad",
                   omics = "rnaseq", coldata = colData(se))
```

Proteomics (prolfqua; the DE table is nested in `rowData`):

```r
se_prot <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds",
                               package = "exploreDE"))
de_prot <- as.data.frame(rowData(se_prot)$constrast_B_vs_A)   # diff, p.value, FDR
build_explore_h5ad(
  x       = list(raw = assay(se_prot, "rawData"),
                 transformed = assay(se_prot, "transformedData")),
  out     = "proteomics.h5ad", omics = "proteomics",
  coldata = colData(se_prot), rowdata = de_prot, factors = "G_",
  de      = list(B_vs_A = de_prot)
)
```

RNA-seq with DE and pathways (SUSHI/DESeq2 result; DE stats live in
`rowData`, clusterProfiler results in `metadata`):

```r
build_explore_h5ad(
  x         = list(counts = assay(se, "counts"), normalised = assay(se, "xNorm")),
  out       = "rnaseq_full.h5ad", omics = "rnaseq",
  organism  = "Homo sapiens", ref_build = metadata(se)$param$refBuild,
  coldata   = colData(se), rowdata = rowData(se),
  factors   = "Condition [Factor]",
  de        = list(diff_vs_undiff = rowData(se)),
  pathways  = list(diff_vs_undiff = list(ora  = metadata(se)$enrichResult,
                                          gsea = metadata(se)$enrichResultGSEA))
)
```

Adapt these to the user's objects; don't copy the example names into their
script.

## How columns are recognised

A dictionary (`system.file("dictionaries/column_patterns.json", package =
"exploreDE")`) holds ordered, case-insensitive regex patterns per concept;
the first pattern that matches any column wins.

| Concept | Recognised by (in order) |
|---|---|
| log2 fold change | `diff`, `log2FoldChange`/`log2FC`/`log2Ratio`/…, `logFC`, `estimate`, then anything containing `log2` |
| p-value | `pvalue`, `p.value`, `P.Value`, `p_value` |
| FDR | `fdr`, `adj.P.Val`, `padj`, `qvalue`, `p.adjusted`, then anything containing `fdr` |
| extra score | anything containing `score` (e.g. `SaintScore`) |
| tested flag | `usedInTest`, `used_in_test` (optional; without it, "tested" = p-value not `NA`) |

So DESeq2, edgeR, limma, prolfqua, and SAINTexpress tables work as they are.
SAINTexpress has `log2_EFCs`, `BFDR`, `SaintScore`, and no p-value: That's
fine, since the app works from the FDR, and the score is shown in the DE table and
as a volcano y-axis (display only; thresholds use the FDR). A table with an
FDR but no raw p-value gives `has_pvalue() == FALSE`, which is correct.

`rowdata` columns are recognised the same way: gene ID (`gene_id`,
`ensembl…`), symbol (`gene_name`, `gene_symbol`, `symbol`), `description`,
`transcript_id`, `GO BP/MF/CC`, protein ID (`protein_Id`, `uniprot`,
`accession`), `peptide_Id`, `site`.

**Feature level** is detected from `rowdata`: `site` → phospho,
`peptide_Id` → peptide, a `transcript_id` with one transcript per row →
isoform, otherwise gene / protein.

### Organisms without an annotation package

Nothing needs an organism database. Leave `organism`/`ref_build` out and use
the user's own gene ids as the counts' rownames.

- KEGG tabs appear only for human, mouse, and rat; don't try to add KEGG.
- GO tabs appear when `rowdata` has any of `GO BP`, `GO MF`, `GO CC`, each a
  `"; "`-separated string of GO ids (`NA` for none). One ontology is fine.
  If the user has one GO list per gene (e.g. eggNOG-mapper's comma-separated
  `GOs`), split it with `AnnotationDbi::Ontology()` (GO.db) into those three
  columns; the vignette's "Organisms without an annotation package" section
  has the code.
- The app never runs enrichment. For ORA/GSEA, run
  `clusterProfiler::enricher()` / `GSEA()` with a `TERM2GENE` built from the
  user's annotation, gene ids = the counts' rownames, and pass the results as
  `pathways = list(<contrast> = list(ora = list(BP = list(upGenes = res))))`
  (ontology, then `upGenes`/`downGenes`/`bothGenes`).

### When the builder warns that it couldn't resolve a column

The warning names the contrast, the missing field, and the columns it saw.

1. **One file:** map the columns explicitly, using names that exist in the
   user's table:

   ```r
   de_columns      = list(effect_log2fc = "<col>", score_pvalue = "<col>",
                          score_fdr = "<col>")   # leave out what doesn't exist
   rowdata_columns = list(id_col = "<col>", symbol_col = "<col>")
   ```

   These five keys are the only ones. Unnamed fields are still found by the
   dictionary. `de_columns` applies to all contrasts, or give a list of such
   lists named by contrast.
2. **A naming convention the user will see again:** suggest adding a pattern
   to `inst/dictionaries/column_patterns.json` in the exploreDE source (a
   JSON edit; anchored `^…$` patterns before loose substrings). Don't do this
   for a one-off file.

## Verify

```r
eds <- new_explore_dataset("my_data.h5ad")
eds                                  # summary: omics, feature level, factors, contrasts, assays
eds$has_de()                         # TRUE if `de` was supplied
eds$has_pvalue()                     # TRUE only if a raw p-value was found
d <- eds$de(eds$primary_contrast())
all(is.na(d$effect_log2fc))          # must be FALSE when `de` was supplied
all(is.na(d$score_fdr))              # must be FALSE when `de` was supplied
```

If either `all(is.na(...))` is `TRUE`, a column wasn't recognised: go back
to the warning, fix it, rebuild. Also check that the printed feature level
and factors are what the user expects.

## Opening the file

```r
library(exploreDE)
exploreDE("my_data.h5ad")   # opens straight into the file
# or: exploreDE(), then upload the .h5ad on the landing page (up to 1 GB)
```

On a server, `?data=<percent-encoded https link to the .h5ad>` opens a file
from an allowed source (APB Export jobs and public pStore by default; the
environment variable `EXPLOREDE_URL_SOURCES` sets the list).

Install with `BiocManager::install("fgcz/exploreDE")`. The full guide is
`vignette("bring-your-own-data", package = "exploreDE")`; the argument
reference is `?build_explore_h5ad`.
