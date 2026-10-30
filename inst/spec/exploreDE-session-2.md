# exploreDE saved-session file — `exploreDE-session/2`

| | |
|---|---|
| **Status** | Current. The only session format since exploreDE 3.9.0. |
| **Date** | 2026-10-30 |
| **Written by** | Save / Load Session tab, "Save session (.json)" (`save_state_to_json()`) |
| **Read by** | Save / Load Session tab, any `.json` upload (`save_state_from_json()`) |
| **Code** | `R/utils_savestate.R`, `inst/app/server-saveState.R` |

A saved session records a user's settings for one dataset: Every tab's
widgets, the feature bucket, the GO/KEGG and pathway selections, and the Multi
DE datasets with their settings. It never holds data. Datasets are loaded
again from their source on restore.

The file is JSON, so it can only hold text, numbers, true/false, null, arrays,
and objects. Reading it can't run code. It opens in any text editor and can be
diffed or edited by hand.

## 1. Encoding

- UTF-8 JSON, one top-level object.
- **Every value is an array**, even a single one: `"gseaType": ["BP"]`. A
  single value and a one-item multi-select read back the same way. An empty
  selection is `[]`; an unset value is `null`.
- Named groups of values are objects (`bucket_list_group` below).
- Numbers keep full precision. Times are ISO 8601 local time,
  `"2026-10-29T12:54:03"`.

## 2. Layout

```json
{
  "meta":   { ... },   // MUST
  "inputs": { ... },   // MUST (may be {})
  "rv":     { ... },   // SHOULD
  "multi":  { ... }    // only when Multi DE datasets were loaded
}
```

Unknown top-level keys are ignored.

### 2.1 `meta`

| Key | Status | Meaning |
|---|---|---|
| `app` | MUST | `["exploreDE"]` |
| `format` | MUST | `["exploreDE-session/2"]`. Any other value is refused. |
| `version` | SHOULD | `[1]`, the bundle version |
| `saved_at` | SHOULD | When it was saved. |
| `dataUrl` | SHOULD | The `?data=` query the dataset was opened with (`""` for uploads and the demo). |
| `fingerprint` | MUST | Identifies the dataset; the session is refused on any other. |

`fingerprint` holds `dataDir`, `dataType`, `omics`, `nSamples`, `nFeatures`,
`factors` (sorted), and `contrasts` (sorted; `[]` when there are none). A
restore requires every field to match the loaded dataset, and `dataDir` too
when both sides have one (`save_state_fp_match()`).

### 2.2 `inputs`

Every restorable widget of the single-dataset tabs, keyed by its Shiny input
id. `save_state_filter_inputs()` decides which: It leaves out uploads, buttons,
table-selection and plot-brush echoes, downloads, the per-tab copies of
`contrastSelected`, and everything starting `multi`.

```json
"inputs": {
  "heatmapColourRed": ["#123456"],
  "bucket_list_group": { "boxKeepBucket": ["diff", "undiff"], "boxExcludeBucket": [] }
}
```

Restore applies values as their widgets appear, re-trying for two minutes
after the upload whenever a tab is opened. Values for widgets this dataset
doesn't have are dropped.

### 2.3 `rv`

App state that isn't a widget value:

| Key | Meaning |
|---|---|
| `featureBucket` | Feature ids in the feature bucket. |
| `goHeatmap` | `{genes, name, term}` of the GO heatmap. |
| `keggHeatmap` | `{genes, name}` of the KEGG heatmap. |
| `goTileKegg` | `{genes, keggId, name}` of the tile-plot KEGG selection. |
| `oraSelected`, `gseaSelected` | Selected row numbers in the ORA / GSEA tables. |
| `gseaType` | The GSEA ontology shown. |

### 2.4 `multi`

Written only when Multi DE datasets are loaded.

```json
"multi": {
  "datasets": [
    { "name": ["courseB"], "source": ["link"],
      "link": ["p1000/DESeq2_Treated--over--Control_2025-11-25--15-08-43/Treated--over--Control"] },
    { "name": ["courseA"], "source": ["link"],
      "link": ["p1000/DESeq2_SDY--over--Control_2025-01-29--14-26-27/SDY--over--Control"] }
  ],
  "inputs": {
    "multiBoxplotNCol": [1],
    "multiGroupColourcourseA": ["#009E73"]
  }
}
```

`datasets` lists the loaded datasets **in display order**, under their
display names (after any renaming):

| `source` | Reference | Restored by |
|---|---|---|
| `link` | `link`: A result directory relative to a gStore root (`p<project>/<analysis>/<contrast>`) | Re-reading it through `multi_resolve_dataset()` and the same project check as the Data Selection tab |
| `contrast` | `contrast`: A contrast of the loaded dataset (proteomics and uploaded `.h5ad` datasets fill Multi DE with their own contrasts) | Taking the matching contrast the loaded dataset provides |

Datasets that can't be found or aren't accessible to the user are skipped,
and the restore message names them. Entries with any other `source` are
dropped when the file is read.

`inputs` holds the Multi DE widgets (`save_state_multi_inputs()`): Every input
starting `multi`, except uploads, the project browser's checkboxes, the order
bucket (the order is in `datasets`), table-selection and download echoes, and
buttons. Group colours are stored as the resolved hex value, so a random
palette comes back the same.

## 3. Limits on reading

The reader refuses, with a short message:

- files larger than 5 MB (a full session with Multi DE is about 25 KB);
- anything that isn't valid JSON or isn't an object;
- `meta.app` other than `exploreDE`, or `meta.format` other than
  `exploreDE-session/2`;
- nesting deeper than 6 levels (a session uses 4).

Strings are only ever used as widget values, feature ids, or dataset links, and
links pass the same checks as links typed into the Data Selection tab.

## 4. The retired `.rds` format

Before this format, a session was saved as an `.rds` file: The same bundle as
an R object. Reading an `.rds` file can run code, so exploreDE 3.9.0 removed it.
An uploaded `.rds` is refused with a message asking the user to save the
session again as `.json`. Nothing converts old files: Open the dataset, set it
up, and save again.

## Version history

| Date | Format | Change |
|---|---|---|
| 2026-09 | `.rds`, bundle version 1 | Save / Load Session (beta), single dataset. |
| 2026-10-29 | `exploreDE-session/2` | JSON format; Multi DE datasets and settings. |
| 2026-10-30 | `exploreDE-session/2` | `.rds` removed (exploreDE 3.9.0); Save / Load Session out of beta. |
