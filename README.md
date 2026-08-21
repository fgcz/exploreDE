# exploreDE

## How to install

```
BiocManager::install(c("fgsea", "DESeq2", "ComplexHeatmap", "clusterProfiler"))
devtools::install_github("uzh/ezRun")
```

## Bring your own data

exploreDE reads one canonical, self-describing file format: an `.h5ad`
(AnnData) file built by `build_explore_h5ad()`. This is the supported way to
get your own data into the app — for any omics type it knows about
(currently `"rnaseq"`, `"proteomics"`, `"crispr"`), from any pipeline,
regardless of that pipeline's column-naming conventions.

Only a count matrix, sample metadata, and an omics type are required. A DE
table and pathway/enrichment results are both optional — the app adapts to
whatever is present. The three examples below span the full range, from the
bare floor to every argument at once.

### Minimal — counts only

The floor: no DE, no rowdata, nothing else. This is a completely valid file —
the app just won't show DE-dependent tabs.

```r
build_explore_h5ad(
  counts,                # a plain matrix, features x samples
  out     = "minimal.h5ad",
  omics   = "rnaseq",
  coldata = coldata      # one row per sample, rownames = counts' colnames
)
```

### Typical — counts + annotation + one contrast

```r
counts  <- read.csv("counts.csv", row.names = 1) |> as.matrix()
coldata <- data.frame(row.names = colnames(counts),
                      group = sub("_[0-9]+$", "", colnames(counts)))
rowdata <- read.csv("gene_annotation.csv", row.names = 1)   # symbols, GO terms, etc.
de      <- read.csv("de_results.csv", row.names = "gene")

build_explore_h5ad(
  counts,
  out     = "my_data.h5ad",
  omics   = "rnaseq",
  coldata = coldata,
  rowdata = rowdata,
  de      = list(Treat_vs_Ctrl = de)
)
```

### Maximal — every argument at once

A real example, from a SUSHI/DESeq2-based RNA-seq analysis with two
contrasts: multiple count layers, organism/ref_build/feature_level/factors,
and per-contrast pathway results (ORA + GSEA).

```r
se1 <- qs2::qs_read("<path to contrast 1's deResult.qs2>")
se2 <- qs2::qs_read("<path to contrast 2's deResult.qs2>")

build_explore_h5ad(
  x             = list(raw = assay(se1, "counts"), normalised = assay(se1, "xNorm")),
  out           = "my_study.h5ad",
  omics         = "rnaseq",
  organism      = "Mus musculus",
  ref_build     = metadata(se1)$param$refBuild,
  coldata       = colData(se1),                   # study metadata
  rowdata       = rowData(se1),                   # feature IDs, symbols, pathway annotations
  feature_level = "gene",                         # gene, isoform, protein, peptide, phospho, sgrna
  factors       = "Condition [Factor]",            # can be a vector; all must exist in coldata
  de = list(
    "KO-inf--over--WT-inf"  = rowData(se1),        # DESeq2/SUSHI keeps DE stats in rowData
    "WT-untr--over--WT-inf" = rowData(se2)
  ),
  pathways = list(
    "KO-inf--over--WT-inf"  = list(ora = metadata(se1)$enrichResult, gsea = metadata(se1)$enrichResultGSEA),
    "WT-untr--over--WT-inf" = list(ora = metadata(se2)$enrichResult, gsea = metadata(se2)$enrichResultGSEA)
  )
)
```

If a column name doesn't match anything in the automatic dictionary (see
below), add `de_columns=`/`rowdata_columns=` — for a real worked case, see
the SAINTexpress producer in `plans/anndata_manifest_design.md`.

You do **not** need to tell it which column is the p-value, FDR, or log2
fold-change — a built-in dictionary (`inst/dictionaries/column_patterns.json`)
recognises common naming conventions automatically (`logFC`, `log2FoldChange`,
`padj`, `adj.P.Val`, `BFDR`, ...). If a column genuinely doesn't match
anything, `build_explore_h5ad()` warns loudly (never silently drops data) and
tells you what to do next — either pass an explicit mapping for that one file
via `de_columns=`/`rowdata_columns=` (see `?build_explore_h5ad`), or add a
pattern to `column_patterns.json` if you expect to see that naming again.

Once built, load and sanity-check the file:

```r
eds <- new_explore_dataset("my_data.h5ad")
eds$has_de()        # TRUE if any contrast was supplied
eds$has_pvalue()     # TRUE only if a raw (non-adjusted) p-value was resolved
eds$de(eds$primary_contrast())
```

See `?build_explore_h5ad` for the full argument reference, and `SKILLS.md`
for an LLM-oriented guide to writing this kind of import script for a user.

