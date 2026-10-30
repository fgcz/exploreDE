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
package. On a server, `?data=` opens a dataset from the link:

- `?data=<result directory>`: A SUSHI result on FGCZ's gStore (LDAP-checked).
- `?data=<percent-encoded https link to an .h5ad>`: A file from an allowed
  source (by default APB Export jobs and public pStore; set
  `EXPLOREDE_URL_SOURCES` to change the list). Plain `https://` only, no
  redirects, up to 1 GB.

## Bring your own data

exploreDE reads one file format: an `.h5ad` (AnnData) file in the
`exploreDE/1.0` layout, for `"rnaseq"`, `"proteomics"`, or `"crispr"` data
from any pipeline. In R, `build_explore_h5ad()` writes it. Only counts,
sample metadata, and the omics type are required; a DE table and pathway
results are optional, and the app shows the tabs your data supports.

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

You don't need to rename columns: A built-in dictionary recognises the usual
DESeq2, edgeR, limma, prolfqua, and SAINTexpress column names, and warns when
it can't place one.

- **Full guide:** `vignette("bring-your-own-data", package = "exploreDE")`
  covers counts-only, proteomics, and RNA-seq + pathway examples (all runnable
  on the data shipped with the package), checking a file, and what to do when
  a column isn't recognised.
- **Using an AI assistant:** give it
  `system.file("skills/exploreDE-byod/SKILL.md", package = "exploreDE")`.
  It's written so an assistant can write the import script for your data.
- **Argument reference:** `?build_explore_h5ad`.
- **Writing files from another tool (Python, an exporter app):**
  `system.file("spec/exploreDE-h5ad-1.0.md", package = "exploreDE")` is the
  written `.h5ad` format spec, including the proteomics profile and a minimal
  Python writer. APB Export's exploreDE output follows it.
- **Also opens as it is:** The prolfquapp DEA's `AnnData.h5ad` (prolfquapp
  2.10.x; spec §10).

## Saving your settings

The **Save / Load Session** tab downloads your settings across every tab
(widgets, the feature bucket, pathway selections, and Multi DE datasets) as
one `.json` file, and restores them on the same dataset later. It never holds
data. The format is `exploreDE-session/2`:
`system.file("spec/exploreDE-session-2.md", package = "exploreDE")`. The
earlier `.rds` session files were retired in 3.9.0.

## Citing exploreDE

See `CITATION.cff`.
