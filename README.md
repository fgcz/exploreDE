# exploreDE

## Install and run

```r
install.packages("BiocManager")
BiocManager::install("fgcz/exploreDE")   # pulls in the Bioconductor dependencies

library(exploreDE)
exploreDE()                   # landing page: upload a file or load the demo
exploreDE("my_data.h5ad")     # or open straight into your own file
```

With no file, the app opens on its landing page, where you can upload an
exploreDE `.h5ad` (see below) or load the demo dataset that ships with the
package. At FGCZ the app is also reachable from SUSHI with
`?data=<result directory>`.

## Bring your own data

exploreDE reads one file format: an `.h5ad` (AnnData) file built by
`build_explore_h5ad()`, for `"rnaseq"`, `"proteomics"` or `"crispr"` data
from any pipeline. Only counts, sample metadata and the omics type are
required; a DE table and pathway results are optional, and the app shows
the tabs your data supports.

```r
build_explore_h5ad(
  counts,                        # features x samples matrix
  out     = "my_data.h5ad",
  omics   = "rnaseq",
  coldata = coldata,             # one row per sample, rownames = colnames(counts)
  rowdata = annotation,          # optional: gene IDs, symbols, GO terms, ...
  de      = list(Treat_vs_Ctrl = de_table)   # optional: one table per contrast
)
```

You don't need to rename columns: a built-in dictionary recognises the usual
DESeq2, edgeR, limma, prolfqua and SAINTexpress column names, and warns when
it can't place one.

- **Full guide:** `vignette("bring-your-own-data", package = "exploreDE")`
  covers counts-only, proteomics and RNA-seq + pathway examples (all runnable
  on the data shipped with the package), checking a file, and what to do when
  a column isn't recognised.
- **Using an AI assistant:** give it
  `system.file("skills/exploreDE-byod/SKILL.md", package = "exploreDE")`.
  It's written so an assistant can write the import script for your data.
- **Argument reference:** `?build_explore_h5ad`.

## Citing exploreDE

See `CITATION.cff`.
