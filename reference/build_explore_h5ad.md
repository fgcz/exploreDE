# Build a canonical exploreDE `.h5ad` from user-supplied data

The source-agnostic "bring your own data" constructor. `x` is a
`DESeqDataSet` (its counts, `colData`, `rowData` and factors are used),
a `SummarizedExperiment` (every assay under its name, its `colData` and
its `rowData`), or the counts themselves (a matrix, or a named list of
matrices such as `list(counts=, normalised=, vst=, tpm=)`). Assay and
matrix names become layer names: `rawData` and `raw` are raw intensities
in proteomics, `transformedData` the analysis-ready ones. DE tables and
pathway statistics are supplied as arguments. Only counts + `coldata` +
`omics` are required; `coldata` and `rowdata`, when given, replace those
of `x`.

## Usage

``` r
build_explore_h5ad(
  x,
  out,
  omics,
  feature_level = NULL,
  coldata = NULL,
  rowdata = NULL,
  de = list(),
  de_columns = NULL,
  rowdata_columns = NULL,
  pathways = list(),
  enrichment_padj_max = 0.1,
  factors = NULL,
  organism = NA_character_,
  ref_build = NA_character_
)
```

## Arguments

- x:

  A `DESeqDataSet`, a `SummarizedExperiment` with named assays, a counts
  matrix, or a named list of count matrices.

- out:

  Output `.h5ad` path.

- omics:

  Required: `"rnaseq"`, `"proteomics"`, or `"crispr"`.

- feature_level:

  Optional override of the feature granularity within the omics — rnaseq
  `"gene"`/`"isoform"`, proteomics `"protein"`/`"peptide"`/`"phospho"`.
  Omit to auto-detect from `rowdata` (`site` -\> phospho, `peptide_Id`
  -\> peptide, `transcript_id` -\> isoform).

- coldata:

  Sample metadata `data.frame`. Required unless `x` is a rich object
  that carries it.

- rowdata:

  Optional per-feature annotation (feature names/symbols/ids and
  per-feature pathway annotation columns, e.g. `GO BP`).

- de:

  Optional named list of DE tables keyed by contrast (`DESeqResults` /
  `data.frame`). Name each contrast `<test>_vs_<reference>`, with both
  groups levels of one `coldata` factor (an error otherwise). Columns
  are recognised by the column dictionary
  (`inst/dictionaries/column_patterns.json`), overridable with
  `de_columns`. An optional extra score column (e.g. SAINTexpress's
  `SaintScore`) is picked up by the dictionary alone and shown for
  display under its own column name.

- de_columns:

  Optional DE column overrides:
  `list(effect_log2fc = "<col>", score_pvalue = "<col>", score_fdr = "<col>")`,
  applied to all contrasts, or a named list of such lists keyed by
  contrast. Named columns are used as given; any not named are still
  found by the dictionary. A producer with no raw p-value (e.g.
  SAINTexpress) simply has none; the app degrades to FDR-only. Unknown
  keys are ignored with a warning. Overrides are persisted at
  `uns$exploreDE$mapping$de`, so the file stays self-describing.

- rowdata_columns:

  Optional feature column overrides:
  `list(id_col = "<col>", symbol_col = "<col>")` — which `rowdata`
  column supplies the feature id and the display symbol. Anything not
  named is still found by the dictionary. Unknown keys are ignored with
  a warning. Overrides are persisted at `uns$exploreDE$mapping$feature`.

- pathways:

  Optional pathway statistics: a named list keyed by contrast (keys a
  subset of `de`); each entry a named list of clusterProfiler objects
  and/or `data.frame`s.

- enrichment_padj_max:

  p.adjust cut-off applied to clusterProfiler results before they are
  stored (`p.adjust <= enrichment_padj_max`). Either one number for
  every method, or a named list keyed by the method names used in
  `pathways` (e.g. `list(ora = 0.05, gsea = 0.05)`) so that ORA and GSEA
  can follow their own FDR thresholds; a method not named in the list
  uses the largest value given. The value(s) are recorded in
  `uns$enrichment_filter` and the app uses a per-method value as its
  significance threshold.

- factors:

  Optional colData columns that are experimental factors (default: all
  character/factor-typed columns).

- organism, ref_build:

  Optional labels.

## Value

The output path, invisibly.

## Examples

``` r
library(SummarizedExperiment)
# A prolfqua proteomics result: one contrast, nested in rowData
se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
de <- as.data.frame(rowData(se)$constrast_B_vs_A)
path <- build_explore_h5ad(
  x = list(raw = assay(se, "rawData"), transformed = assay(se, "transformedData")),
  out = tempfile(fileext = ".h5ad"),
  omics = "proteomics",
  coldata = colData(se),
  rowdata = de,
  factors = "G_",
  de = list(B_vs_A = de)
)
file.exists(path)
#> [1] TRUE

# The SummarizedExperiment itself: rawData and transformedData become the
# layers intensities_raw and intensities
se_path <- build_explore_h5ad(
  se,
  out = tempfile(fileext = ".h5ad"),
  omics = "proteomics",
  rowdata = de,
  factors = "G_",
  de = list(B_vs_A = de)
)

# Counts only
counts_only <- build_explore_h5ad(
  assay(se, "rawData"),
  out = tempfile(fileext = ".h5ad"),
  omics = "proteomics",
  coldata = colData(se)
)
```
