# exploreDE `.h5ad` format — `exploreDE/1.0`

| | |
|---|---|
| **Status** | Current. Describes what exploreDE 3.9.0 reads and writes. |
| **Date** | 2026-10-30 |
| **Supersedes** | The 2026-04-17 draft (`plans/archive/anndata_migration/anndata_spec_2026-04-17.qmd`), which the implementation has since departed from. |
| **Reference reader** | `R/ExploreDatasetH5AD.R` (via `new_explore_dataset()`) |
| **Reference writer** | `build_explore_h5ad()` in `R/build_explore_h5ad.R` |
| **Machine-readable parts** | `inst/schema/exploreDE_mapping.schema.json` (the `uns["exploreDE"]["mapping"]` block) and `inst/dictionaries/column_patterns.json` (column dictionary used by the writer) |

This document is for anyone producing exploreDE files outside
`build_explore_h5ad()` — another R package, a Python pipeline, or an exporter
app. If you are in R and can call `build_explore_h5ad()`, do that instead: It
applies every rule below for you.

The key words MUST, SHOULD, and MAY are used as in RFC 2119. "Reader" means
exploreDE's `ExploreDatasetH5AD`.

---

## 1. Container

- A standard AnnData file (`.h5ad`, HDF5), in the on-disk encoding written by
  Python `anndata` ≥ 0.8 or R `anndataR`. Both have been checked against the
  reader (Python `anndata` 0.10.5 and 0.12, R `anndataR` 1.2.3 and 1.3.1).
  The deployed app (ShinyProxy `fgcz_base_image`) reads with `anndataR` 1.2.3.
- Text MAY be stored in any standard encoding: `string-array`, `categorical`,
  or `nullable-string-array` (what Python `anndata` ≥ 0.12 writes for pandas
  string columns). Reading `nullable-string-array` needs `anndataR` ≥ 1.2.2.
  Older `anndataR` versions only warn, and the file loads without its sample
  names, factors, and feature ids. To open a file there, cast text to `object`
  before writing:

  ```python
  def plain_strings(df):
      df = df.astype({c: object for c in df.columns if isinstance(df[c].dtype, pd.StringDtype)})
      df.index = df.index.astype(object)
      return df
  adata.obs, adata.var = plain_strings(adata.obs), plain_strings(adata.var)
  ```
- **One modality per file.** MuData (`.h5mu`) is not read.
- **Axes are AnnData's:** `obs` = samples, `var` = features. Every matrix is
  samples × features (the transpose of a SummarizedExperiment).
- Matrices SHOULD be dense. Missing values are `NaN` (proteomics intensities
  are expected to have them).

```text
X                          samples × features   (a copy of the primary matrix, see §4)
layers[<name>]             samples × features   every matrix the app should use
obs                        sample table; experimental factors are "factor_*" columns
var                        feature table; index = feature_id
varm["de_<slot>"]          one data frame per contrast (optional)
uns["schema_version"]      "exploreDE/1.0"
uns["omics_type"]          "proteomics" | "rnaseq_bulk" | "crispr_countqc"
uns["feature_id_kind"]     what one var row is
uns["hierarchy"]           feature levels
uns["design"]              factors and default grouping
uns["de"]                  per-contrast metadata, keyed by <slot> (may be empty)
uns["enrichment"]          pathway tables per contrast (may be empty)
uns["enrichment_support"]  what the app needs to show / rebuild them
uns["enrichment_filter"]   the cut-off the stored tables were filtered at
uns["exploreDE"]["mapping"] which source column fed which statistic (optional)
uns["reference_genome"], uns["provenance"]
```

## 2. `obs` — samples

| Column | Status | Type | Meaning |
|---|---|---|---|
| index (`obs_names`) | MUST | unique strings | Sample ids. |
| `factor_<name>` | SHOULD (≥ 1) | categorical | Experimental factors. The reader exposes them without the prefix (`factor_condition` → `condition`). Category order is the display order (boxplot x-axis, legend, colours), so set it deliberately. |
| `sample_label` | SHOULD | string | Display name. The reference writer copies the sample id. |
| `organism` | MAY | string | E.g. `"Homo sapiens"`. The reader takes the first value. |
| anything else | MAY | any | Carried through as sample annotation. |

Rules:

- Only `factor_*` columns are treated as factors. Do not prefix identifier
  columns (one value per sample, e.g. raw file names) or constant columns
  (e.g. `isotopeLabel = "light"`): They give a useless grouping.
- The default grouping is `uns["design"]["grouping"]`, else the first
  `factor_*` column. The second `factor_*` column, if any, becomes the
  secondary grouping.

## 3. `var` — features

| Column | Status | Type | Meaning |
|---|---|---|---|
| index (`var_names`) | MUST | unique strings | `feature_id`: The join key for every matrix, every DE table, and the app's feature pickers. See §6 for the proteomics conventions. |
| `label_gene_symbol` | SHOULD | string | Display name in plots, tables, and search. Empty or missing values fall back to `feature_id`. |
| `description` | SHOULD | string | Free text; searchable. |
| `label_gene_id` | MAY | string | Stable gene id (RNA-seq: Ensembl). |
| `label_protein_id` | proteomics: SHOULD | string | Protein accession. See §6. |
| `label_go_bp`, `label_go_mf`, `label_go_cc` | MAY | string | Per-feature GO annotation; the reader renames them to `GO BP` / `GO MF` / `GO CC`. |
| a length column (`featWidth`, `gene_length`, `length`, …) | RNA-seq: MAY | number | Enables computed TPM and FPKM. |
| any other `label_*` | MAY | any | Carried through as annotation. |

## 4. `X` and `layers` — matrices

The reader takes every matrix from `layers`. `X` is only used when `layers`
has no `counts` entry, in which case `X` is read as `counts`. So:

- Every matrix the app should use MUST be in `layers`, under the names below.
- `X` SHOULD be a copy of the primary matrix (other AnnData tools expect
  one); the reader ignores it whenever `layers` has `counts`.

### Proteomics (`uns["omics_type"] = "proteomics"`)

| Layer | Status | Shown as | Content |
|---|---|---|---|
| `intensities` | MUST | Intensities | The analysis-ready matrix: Transformed (log2-scale), normalised intensities, `NaN` where missing. |
| `intensities_raw` | SHOULD | Raw intensities | Untransformed intensities. |

- If `intensities` is absent, the reader falls back to `normalised`, `log2`,
  or `normalised_log2`, then to whichever layer comes first. Raw falls back to
  `counts` or `raw`. New writers SHOULD use the two canonical names.
- RNA-seq derived assays (VST, TPM, and FPKM) are never computed for
  proteomics. The low-signal floor is 1 (log2 scale), not 10.
- If you have no raw intensities, omit `intensities_raw` and set `X` to the
  same values as `intensities`. (The reader then shows the `X` copy as
  "Raw intensities" as well — harmless, but known.)

### RNA-seq (`"rnaseq_bulk"`)

| Layer | Status | Content |
|---|---|---|
| `counts` | MUST | Raw counts. |
| `normalised` | MAY | Library-size normalised counts; computed (DESeq2 size factors) if absent. |
| `log2` | MAY | `log2(normalised + pseudocount)`; computed if absent. |
| `vst` | MAY | Variance-stabilised counts; computed if absent. |
| `tpm`, `fpkm` | MAY | Computed if absent and `var` has a length column. |

### CRISPR screen counts (`"crispr_countqc"`)

`counts` = raw sgRNA counts; `feature_id_kind = "sgrna"`.

## 5. `varm["de_<slot>"]` — differential expression, one per contrast

Optional as a whole. A file with no contrasts loads and the DE and pathway
tabs are hidden.

Each entry MUST be a data frame (`encoding-type: dataframe`) with one row per
feature, in `var` order, and its index equal to `var_names`. Features that were
not tested are `NaN`, or are marked `false` in `used_in_test`.

| Column | Status | Meaning |
|---|---|---|
| `effect_log2fc` | MUST | log2 fold change (volcano x-axis). |
| `score_fdr` | MUST | Adjusted p-value / q-value / FDR. **This is the significance axis**: Every threshold in the app uses it. |
| `score_pvalue` | MAY | Raw p-value (volcano y-axis). Leave it out if the method has none (e.g. SAINTexpress); the app then works from the FDR alone. |
| `score` | MAY | One extra display-only score (e.g. SaintScore). Shown in the DE table and offered as a volcano y-axis, never used for thresholds. Its display label comes from `uns["exploreDE"]["mapping"]` (§7.6); without a label it is shown as `score`. |
| `used_in_test` | MAY | Boolean: Whether the feature entered the test (e.g. SUSHI's `usedInTest`, the expression filter). The DE views show only tested features. Without this column, a feature counts as tested when its p-value (or, with no p-values, its FDR) is not `NaN`. Readers also accept 0/1. |
| other columns | MAY | Passed through with the DE table; not interpreted. |

`<slot>` SHOULD match `[a-z0-9_]+` (the reference writer lower-cases the
contrast name and replaces anything else with `_`). The same `<slot>` keys
`uns["de"]`, `uns["enrichment"]`, `uns["enrichment_support"]`, and
`uns["exploreDE"]["mapping"]["de"]["contrasts"]`.

## 6. Proteomics profile

### Feature levels

| `uns["feature_id_kind"]` | `uns["hierarchy"]` | One `var` row is | Recommended `feature_id` |
|---|---|---|---|
| `protein` | `{levels: ["protein"], unit: "protein"}` | a protein (group) | the protein accession |
| `peptide` | `{levels: ["protein", "peptide"], unit: "peptide"}` | a peptide | `<protein>_<peptide>` |
| `phospho` | `{levels: ["protein", "site"], unit: "site"}` | a modification site | `<protein>_<site>` |
| `precursor` | `{levels: ["protein", "peptide", "precursor"], unit: "precursor"}` | a precursor ion | `<protein>_<ProForma ion>`, e.g. `Q5VYK3_AAAAASASQDELNQLER/2` |

- `feature_id_kind` sets the app's data type: `protein` and `precursor` →
  proteomics, `peptide` → peptidomics, `phospho` → phosphoproteomics.
  `build_explore_h5ad()` does not offer `precursor` yet; files written
  elsewhere load (checked on DIA-NN, MaxQuant, and FragPipe precursor exports).
- Any unique string is a valid `feature_id`; the recommendations match what
  the reference writer builds (`~` in a site is replaced by `_`), so files from
  either route look alike.
- For `peptide`, `phospho`, and `precursor`, `var` SHOULD carry
  `label_protein_id` (and `label_gene_symbol`): The feature picker then offers
  search by protein and by gene, each mapping to all rows of that protein or
  gene.
- `label_gene_symbol` is the label on every plot and table. When there is no
  gene name, a short name such as the UniProt entry name of the first protein
  (`ENOA_HUMAN`) reads far better than a whole protein group.

### Without a DEA

`obs` + `var` + `layers["intensities"]` (+ `intensities_raw`) + `uns` with
`uns["de"] = {}` and `uns["enrichment"] = {}`. Leave `varm` empty.

### With a DEA

As above, plus one `varm["de_<slot>"]` and one `uns["de"][<slot>]` per
contrast. A prolfqua contrast table maps as `diff` → `effect_log2fc`,
`p.value` → `score_pvalue`, `FDR` → `score_fdr`.

## 7. `uns` — metadata

### 7.1 Identity — all MUST

| Key | Value |
|---|---|
| `schema_version` | `"exploreDE/1.0"` |
| `omics_type` | `"proteomics"`, `"rnaseq_bulk"`, or `"crispr_countqc"` |
| `feature_id_kind` | proteomics: `protein`, `peptide`, `phospho`, `precursor`; RNA-seq: `gene`, `isoform`; CRISPR: `sgrna` |
| `hierarchy` | `{levels: [...], unit: "..."}` as in §6 (RNA-seq: `["gene"]`, or `["gene", "transcript"]` with unit `transcript`) |

### 7.2 `design` — SHOULD

```yaml
design:
  factors:  ["factor_condition", "factor_batch"]   # informational
  grouping: "factor_condition"                     # default grouping
```

### 7.3 `de` — SHOULD be present, MAY be empty

One entry per contrast, keyed by `<slot>` (same as `varm["de_<slot>"]`).
Write `{}` when there are none. exploreDE ≥ 3.9.0 reads a file without the
`de` key as having no contrasts; exploreDE 3.8.0 does not load it, because
there `uns$de` partial-matches `uns$design`.

| Field | Status | Meaning |
|---|---|---|
| `test_group` | MUST | A level of the grouping factor. Drives the boxplot group order and the per-group means. |
| `reference_group` | MUST | A level of the grouping factor. |
| `contrast_formula` | SHOULD | The contrast's display name and its key inside the app. Default: `<test>_vs_<reference>`, else `<slot>`. |
| `factor_used` | SHOULD | The `factor_*` column the contrast is on. |
| `model` | SHOULD | `"prolfqua"`, `"DESeq2"`, `"limma"`, … (informational). |
| `assay_used` | MAY | Which layer the test ran on (informational). |

```yaml
de:
  trt_vs_ctrl:
    test_group: "trt"
    reference_group: "ctrl"
    contrast_formula: "trt_vs_ctrl"
    factor_used: "factor_condition"
    model: "prolfqua"
    assay_used: "intensities"
```

### 7.4 `enrichment`, `enrichment_support`, `enrichment_filter` — pathways

Optional. Write `enrichment = {}` and `enrichment_support = {}` when there are
none. The pathway tabs (ORA, GSEA, tile plots, KEGG) are shown only for
RNA-seq files with DE (§9); for proteomics, pathway results are stored and
read but not displayed.

**`enrichment[<slot>][<method>]`** — one tidy data frame per method, where
`<method>` is exactly `ora` or `gsea`. Columns are those of clusterProfiler's
result tables, under clusterProfiler's names; typically:

| Method | Columns |
|---|---|
| `ora` | `ID`, `Description`, `GeneRatio`, `BgRatio`, `pvalue`, `p.adjust`, `qvalue`, `geneID` (`/`-separated), `Count` |
| `gsea` | `ID`, `Description`, `setSize`, `enrichmentScore`, `NES`, `pvalue`, `p.adjust`, `qvalue`, `rank`, `leading_edge`, `core_enrichment` (`/`-separated) |
| both | `ontology`: `"<domain>/<direction>"` for ORA (e.g. `BP/upGenes`, `BP/downGenes`), `"<domain>"` for GSEA (e.g. `BP`); `method`: `ora` or `gsea` |

**`enrichment_support[<slot>][<method>]`** — SHOULD be present for every method
in `enrichment`; without it the dataset reports no enrichment
(`has_enrichment()` is `FALSE`) and no clusterProfiler objects are rebuilt:

| Field | Status | Meaning |
|---|---|---|
| `class` | MUST | `"enrichResult"` (ORA) or `"gseaResult"` (GSEA). |
| `gene_sets` | GSEA: SHOULD | Data frame `ontology`, `ID`, `gene`: Members of the stored terms only. |
| `gene_list` | GSEA: SHOULD | Data frame `ontology`, `gene`, `value`: The ranked list, in rank order; `ontology = "*"` when all domains share one list. |
| `exponent` | GSEA: SHOULD | fgsea weighting exponent (usually 1). |

With these, the reader rebuilds clusterProfiler objects for the network, ridge,
and running-score plots (when DOSE is installed); the tables alone drive
everything else.

**`enrichment_filter`** — SHOULD be present when tables were filtered:
`{padj_max: <number> | {ora: <number>, gsea: <number>}, column: "p.adjust"}`.
The app uses `padj_max` as its pathway threshold, and it tells a consumer that
a missing term was cut, not untested.

### 7.5 Provenance — SHOULD

- `reference_genome`: String, e.g. `"Homo_sapiens/GENCODE/GRCh38.p13"`.
- `provenance`: Dict; at least `source_format` (the producing tool).

### 7.6 `exploreDE.mapping` — MAY

Records which source column fed each canonical statistic, so the file
describes itself. Validated by `inst/schema/exploreDE_mapping.schema.json`:

```yaml
exploreDE:
  mapping:
    de:
      default:   {effect_log2fc_col: "diff", score_pvalue_col: "p.value", score_fdr_col: "FDR"}
      contrasts:
        trt_vs_ctrl: {score_col: "SaintScore", score_label: "SaintScore"}
    feature: {id_col: "UniProt_Accession", symbol_col: "Gene"}
```

The reader currently uses only `score_label` (merged as `default`, then the
contrast's entry). Everything else is provenance.

## 8. Checking a file

In R, with exploreDE installed:

```r
library(exploreDE)
eds <- new_explore_dataset("my_data.h5ad")
validate_explore_dataset(eds)        # stops on structural problems
eds$data_type(); eds$feature_id_kind(); eds$assay_names()
eds$has_de(); eds$has_pvalue(); eds$has_enrichment(); eds$contrasts()
d <- eds$de(eds$primary_contrast())
c(all(is.na(d$effect_log2fc)), all(is.na(d$score_fdr)))   # both FALSE
```

Then open it: `exploreDE("my_data.h5ad")`.

### Minimal proteomics file in Python

Checked against the reader (Python `anndata` 0.10.5): It loads as proteomics,
protein level, with one contrast. Drop the `varm` entry and set `uns["de"]` to
`{}` for the no-DEA case.

```python
import numpy as np, pandas as pd, anndata as ad

samples = [f"S{i}" for i in range(1, 7)]
proteins = [f"P{i:05d}" for i in range(1, 51)]
rng = np.random.default_rng(1)
log2_int = rng.normal(20, 2, (6, 50))
log2_int[0, :3] = np.nan                                 # missing values stay NaN

obs = pd.DataFrame({
    "sample_label": samples,
    "factor_condition": pd.Categorical(["ctrl"] * 3 + ["trt"] * 3,
                                       categories=["ctrl", "trt"]),
}, index=samples)
var = pd.DataFrame({
    "label_gene_symbol": [f"G{i}" for i in range(1, 51)],
    "label_protein_id": proteins,
    "description": [f"protein {p}" for p in proteins],
}, index=proteins)
de = pd.DataFrame({
    "effect_log2fc": rng.normal(0, 1, 50),
    "score_pvalue": rng.uniform(0, 1, 50),
    "score_fdr": rng.uniform(0, 1, 50),
}, index=proteins)                                       # index MUST equal var_names

adata = ad.AnnData(
    X=log2_int.copy(), obs=obs, var=var,
    layers={"intensities": log2_int.copy(), "intensities_raw": 2 ** log2_int},
    varm={"de_trt_vs_ctrl": de},
    uns={
        "schema_version": "exploreDE/1.0",
        "omics_type": "proteomics",
        "feature_id_kind": "protein",
        "hierarchy": {"levels": ["protein"], "unit": "protein"},
        "design": {"factors": ["factor_condition"], "grouping": "factor_condition"},
        "de": {"trt_vs_ctrl": {
            "test_group": "trt", "reference_group": "ctrl",
            "contrast_formula": "trt_vs_ctrl", "factor_used": "factor_condition",
            "model": "prolfqua", "assay_used": "intensities"}},
        "enrichment": {}, "enrichment_support": {},
        "provenance": {"source_format": "my_exporter"},
    },
)
adata.write_h5ad("my_proteomics.h5ad")
```

## 9. Reader limitations (exploreDE 3.9.0)

These are reader behaviours, not format rules; the MUSTs above already work
around them.

1. **Pathway tabs are RNA-seq only.** The menu shows them for
   `omics_type = "rnaseq_bulk"` with DE; it does not consult the file's
   enrichment. `has_enrichment()` counts only methods with an
   `enrichment_support` entry.
2. **`X` doubles as `counts`** when `layers` has no `counts` (§4).
3. **R-written files may not open in Python.** `anndataR` before 1.2.3 (and
   the 1.3.1 development build) writes each `varm` data frame's `_index` as
   `1..n` instead of the var names, and Python `anndata` refuses the file
   ("value.index does not match parent's var names"). The exploreDE reader
   doesn't mind. Fixed in `anndataR` 1.2.3 (scverse/anndataR#525).
4. **Fixed significance cut-offs.** DE thresholds default to FDR 0.05 and
   |log2FC| 0 in the app; they are not read from the file.

## 10. Accepted dialect — prolfquapp DEA `AnnData.h5ad`

The prolfquapp DEA (B-Fabric A414) writes `AnnData.h5ad` next to
`SummarizedExperiment.rds`. It is **not** `exploreDE/1.0`: Its own namespace,
`uns["prolfquapp"]`, holds the SummarizedExperiment metadata, and its `varm`
holds the SE's nested rowData frames. exploreDE detects it by the presence of
`uns["prolfquapp"]` (this wins over everything else in the file), rebuilds the
SummarizedExperiment (`se_from_prolfquapp_anndata()` in `R/prolfquapp_h5ad.R`),
and reads that with the proteomics SE backend. A new producer SHOULD write
`exploreDE/1.0` instead; this section exists so prolfquapp output keeps working
and so the two can be told apart.

### Shared by both prolfquapp layouts

| Slot | Content |
|---|---|
| `layers` | The SE assays, by name: `rawData`, `transformedData`, and, depending on version and model, `nr_children`, `imputedData` (`lm_impute` model), `ibaq`. Order recorded in `uns["prolfquapp"]["layer_names"]`. |
| `obs` | `colData`: Sample name, `raw.file`, the grouping column (e.g. `G_`), `CONTROL`, … The factors in `uns["prolfquapp"]["analysis_configuration"]["factors"]`, modelled factor first, become the groupings. |
| `var` | `rowData$annotation`: `protein_Id` (+ `peptide_Id` / `site`), `description`, `nrPeptides`, `fasta.id`, … |
| `varm[<key>]` | One per nested rowData frame: `constrast_<contrast>` (sic — prolfquapp's spelling; the contrast name follows), `stats_normalized_wide`, `stats_raw_wide`, `imputation`. Keys are URL-encoded; order in `uns["prolfquapp"]["varm_key_order"]`. Contrast frames hold prolfqua's columns (`diff`, `std.error`, `avgAbd`, `statistic`, `df`, `p.value`, `conf.low`, `conf.high`, `sigma`, `FDR`, `modelName`, …). |
| `uns["prolfquapp"]` | The SE metadata, plus bookkeeping keys (`layer_names`, `uns_table_columns`, `varm_key_order`, …) recording what HDF5 does not keep (key order, tables stored column by column). |

### prolfquapp 2.10.0–2.10.4 — read by exploreDE

Each `varm` entry is a **numeric matrix**. Its column names are in
`uns["prolfquapp"]["varm_columns"][<key>]`, its non-numeric columns in
`varm_annotations[<key>]`, and the final column order in
`varm_column_order[<key>]`. Contrast frames carry no feature keys; exploreDE
joins `var` onto each on the feature id.

### prolfquapp ≥ 2.11.0 — not yet read by exploreDE

- Each `varm` entry is a **data frame** (with its key columns), and
  `varm_columns` / `varm_annotations` / `varm_column_order` are gone.
- `uns["prolfquapp"]["uns_list_order"]` records metadata list order, and
  `uns["X_layer_name"] = "transformedData"` names the `X` copy.
- `uns["prolfquapp"]` MUST hold `artifact_type = "dea_results"`,
  `schema_version` (prolfquapp's own; `"2.1.0"` at 2.11), `source_software`,
  `analysis_configuration`, `layer_names`, `contrasts`, `formula`, and
  `provenance`.
- A DEA of an AnnData input (e.g. an apb-export file) also writes
  `MuData.h5mu` (modalities `lfqdata` and `dea`). exploreDE reads neither the
  `.h5mu` nor an `artifact_type = "lfqdata"` AnnData.

exploreDE 3.9.0 stops on this layout with "varm matrix has N columns but
varm_columns names 0".

---

## Version history

| Date | Version | Change |
|---|---|---|
| 2026-04-17 | draft 0.1 | Planning draft (archived). |
| 2026-10-29 | 1.0 | First written spec of what the reader and writer actually do; proteomics profile; prolfquapp dialect. |
| 2026-10-30 | 1.0 | Format unchanged. `nullable-string-array` text is read (the deployed app has `anndataR` 1.2.3); `uns["de"]` relaxed from MUST to SHOULD (exploreDE 3.9.0 reader fix); optional `used_in_test` DE column. |
