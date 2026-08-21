# exploreDE: building a compatible data file — instructions for an LLM assistant

You are helping a user (who may have no knowledge of exploreDE's internals)
turn their own data — counts, a DE results table, maybe pathway results —
into a file exploreDE can load. This document is everything you need to do
that correctly. Read it in full before writing any code.

## What you're building

exploreDE reads exactly one file format: a canonical `.h5ad` (AnnData) file,
produced by the R function `build_explore_h5ad()`. There is no other
supported way to get user data into the app. Do not hand-write `.h5ad`
structure yourself, and do not invent an alternative format — always go
through this function.

## The one function you need

```r
build_explore_h5ad(
  x, out, omics,
  feature_level = NULL, coldata = NULL, rowdata = NULL,
  de = list(), de_columns = NULL, rowdata_columns = NULL,
  pathways = list(), factors = NULL,
  gene_symbol = NULL, organism = NA, ref_build = NA
)
```

Run `?build_explore_h5ad` for the full reference if it's available in the
user's R session; the summary below is enough for the common case.

**Required:**
- `x` — either a counts matrix (features × samples, or samples × features —
  check which axis is longer / has gene-like rownames before transposing), a
  named list of matrices (e.g. `list(counts=, vst=, tpm=)`), or a rich object
  the function already understands (`DESeqDataSet`, `SummarizedExperiment`).
  If `x` is one of those rich objects it may already carry `coldata`/`rowdata`
  — don't ask the user to re-supply what the object already has.
- `out` — where to write the `.h5ad`.
- `omics` — one of `"rnaseq"`, `"proteomics"`, `"crispr"`. This is a **closed
  set** right now. If the user's data is something else (metabolomics,
  lipidomics, etc.), do not invent a new omics string — tell the user this
  isn't supported yet rather than guessing a value that will error or silently
  mis-tag the data.
- `coldata` — a data.frame of sample metadata, one row per sample, rownames
  matching the counts matrix's sample columns. Skip only if `x` already
  supplies it.

**Optional — the app adapts to what's present, nothing else needs to change:**
- `rowdata` — per-feature annotation (gene symbol, protein ID, description,
  GO terms, ...).
- `de` — a **named list** of DE result tables, one entry per contrast. The
  name of each list entry becomes the contrast's identifier in the app
  (convention: `"treated_vs_control"`).
- `pathways` — enrichment/ORA/GSEA results, named list keyed by the same
  contrast names as `de`.

If the user only has counts, that's a valid, complete file — DE and pathways
are genuinely optional, not "usually required."

## Worked examples — minimal, typical, maximal

Match the script you write to what the user actually has. Don't add
arguments they don't need, and don't invent placeholder data to fill in ones
they didn't mention.

### Minimal — counts only

The floor. No DE, no rowdata. Still a completely valid file — the app just
won't show DE-dependent tabs.

```r
build_explore_h5ad(
  counts,                # a plain matrix, features x samples
  out     = "minimal.h5ad",
  omics   = "rnaseq",
  coldata = coldata      # one row per sample, rownames = counts' colnames
)
```

### Typical — counts + annotation + one contrast

The common case: a user with a counts table, some feature annotation, and
one DE result table from an ordinary tool (DESeq2, limma, edgeR, ...).

```r
counts  <- read.csv("counts.csv", row.names = 1) |> as.matrix()
coldata <- data.frame(row.names = colnames(counts),
                      group = sub("_[0-9]+$", "", colnames(counts)))
rowdata <- read.csv("gene_annotation.csv", row.names = 1)   # symbols, GO terms, etc.
de      <- read.csv("de_results.csv", row.names = "gene")   # e.g. logFC, P.Value, adj.P.Val

build_explore_h5ad(
  counts,
  out     = "my_data.h5ad",
  omics   = "rnaseq",
  coldata = coldata,
  rowdata = rowdata,
  de      = list(Treat_vs_Ctrl = de)
)
```

### Maximal — every argument, from a real SUSHI/DESeq2 analysis

Multiple count layers, organism/ref_build/feature_level/factors, and two
contrasts each with pathway results. Note DESeq2/SUSHI's own convention of
keeping DE stats *inside* `rowData()` rather than a separate table — don't
"fix" that, it's the correct real-world shape for that producer.

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

If the user's data instead needs an explicit column-name override
(`de_columns=`/`rowdata_columns=`), see the next section — and ground any
example you write in the user's *actual* column names, never a fabricated
naming convention.

## The column-detection dictionary — the part that means you don't hardcode anything

This is the single most important thing to understand: **you do not need to
tell `build_explore_h5ad()` which column is the p-value, which is the FDR,
which is the log2 fold-change, or which is the gene/protein ID** for common
naming conventions. A dictionary at `inst/dictionaries/column_patterns.json`
already recognises a wide range of real-world spellings — `logFC`,
`log2FoldChange`, `diff`, `pvalue`, `p.value`, `padj`, `adj.P.Val`, `qval`,
`BFDR`, `protein_Id`, `UniProt`, `gene_name`, `nrPeptides`, and more — by
regex, including substring matches (so `BFDR` is recognised via `fdr`,
`log2_EFCs` via `log2`, `SaintScore` via `score`).

**So: for a normal DE table (DESeq2, limma, edgeR, prolfqua, etc.), just pass
it as-is.** Do not manually rename columns before calling
`build_explore_h5ad()`. Let it do the recognition.

```r
de <- read.csv("de_results.csv", row.names = "gene")   # e.g. logFC, P.Value, adj.P.Val
build_explore_h5ad(counts, out = "data.h5ad", omics = "rnaseq",
                   coldata = coldata, de = list(treated_vs_control = de))
```

### When it can't figure something out

If a column name carries no recognisable textual signal at all (a genuinely
opaque or idiosyncratic name), `build_explore_h5ad()` **warns loudly** —
never silently — naming the contrast, which canonical fields it couldn't
resolve, and the exact source columns it saw. If you see that warning:

1. **One-off fix** (this file only): pass an explicit mapping.
   ```r
   de_columns = list(effect_log2fc = "<source column>",
                     score_fdr     = "<source column>",
                     score_pvalue  = "<source column>")   # omit if genuinely absent
   ```
   `rowdata_columns = list(id_col = "<col>", symbol_col = "<col>",
   description_col = "<col>")` does the same for feature identity.
2. **Permanent fix** (this naming convention will recur): tell the user this
   requires adding a pattern to `inst/dictionaries/column_patterns.json` — a
   JSON edit, not an R code change. Don't do this yourself for a one-off file;
   only suggest it when the user says they'll see this naming again.

Never invent a de_columns mapping to a column that doesn't actually exist in
the user's data.frame — check `colnames()` first.

### Extra, domain-specific columns (SaintScore, nrPeptides, etc.)

These don't need any mapping at all — every column in `de`/`rowdata` that
doesn't match a canonical target rides through automatically and is available
for display. Only use `de_columns$score=`/`score_kind=` if the extra column is
a genuine *second significance axis* you want the app to know the direction
of (e.g. SAINTexpress's `SaintScore`, where high = significant):

```r
de_columns = list(effect_log2fc = "log2_EFCs", score_fdr = "BFDR",
                  score = "SaintScore", score_kind = "probability_high_sig")
```
(`score_kind` / `fdr_kind` enum: `"p_value_low_sig"`, `"fdr_low_sig"`,
`"probability_high_sig"`.)

## Always verify the output before telling the user you're done

```r
eds <- new_explore_dataset("data.h5ad")
eds$has_de()                              # sanity: TRUE if de= was supplied
eds$has_pvalue()                          # TRUE only if a RAW p-value resolved (FDR alone -> FALSE, this is correct, not a bug)
d <- eds$de(eds$primary_contrast())
all(is.na(d$effect_log2fc))               # should be FALSE if a DE table was supplied
all(is.na(d$score_fdr))                   # should be FALSE if a DE table was supplied
```

If either of those `all(is.na(...))` checks comes back `TRUE`, something
didn't resolve — go back and check the warning message, don't tell the user
the file is ready.

## Common mistakes to avoid

- Don't manually rename or pre-process DE/rowdata columns before calling
  `build_explore_h5ad()` — that defeats the automatic dictionary and can mask
  real problems.
- Don't guess at an `omics=` value outside the three supported strings.
- Don't fabricate a `de_columns=`/`rowdata_columns=` mapping to a column name
  that isn't actually in the user's data.
- Don't skip the verification step above — "no error was thrown" is not the
  same as "the DE columns resolved correctly" (an unresolved column produces
  `NA`, not an error).
- Counts orientation matters: exploreDE expects features (genes/proteins) as
  rows, samples as columns. If unsure, check which axis has sample-like names
  (matching `coldata` rownames) and which has feature-like names.

## Reference

- Full argument documentation: `?build_explore_h5ad` (or read the roxygen
  comment block directly above `build_explore_h5ad` in
  `R/build_explore_h5ad.R` if help isn't available in your context).
- Design rationale for the dictionary mechanism:
  `plans/anndata_manifest_design.md`.
- The pattern dictionary itself: `inst/dictionaries/column_patterns.json`.
