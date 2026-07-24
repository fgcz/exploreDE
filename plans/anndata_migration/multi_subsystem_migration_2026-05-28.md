# Multi-* subsystem migration — exploreDE AnnData Phase 2 + Phase 3

**Date:** 2026-05-28
**Branch:** `anndata-conversion`
**Related:** `internal_api_spec_2026-04-30_thin.qmd`, `module_inventory_2026-04-30.md`, `anndata_spec_2026-04-17.qmd`

---

## Context

The multi-* tab subsystem in exploreDE is the second of two "data pipelines"
in the app. The single-dataset path was migrated to `ExploreDataset` via
the Phase 1 legacy shim and is well underway through Phase 2. The multi-*
subsystem still runs on a parallel, gStoreList-shaped pipeline:
`multiDeReactive$results[[name]]` holds objects produced by directly
calling the legacy `convert_genomics_se()`, with a separate manually-built
`pseudoObj` for proteomics. As a result, the multi-* tab is hard-coupled
to gStore-resident `.rds`/`.qs2` SE files — researchers with their own data
cannot use it at all.

This migration does two things at once:

1. **Moves multi-* onto the `ExploreDataset` (`eds`) class**, completing
   the Phase 2 migration for the multi-* subsystem.
2. **Lands the `ExploreDatasetH5AD` subclass** (spec §1.5, originally
   slated for Phase 3) and adds a real BYOD entry point to multi-*. After
   this PR, a researcher with no FGCZ access can drop in an `.h5ad` (or
   `.rds`/`.qs2`) and use the multi-* tabs end-to-end.

Decisions locked in with Peter (this session, 2026-05-28):

- **Data sources**: full BYOD — file upload widget + generalised URL + existing Browse Projects, all in this PR.
- **Per-slot shape**: `list(eds = eds, contrast = "contrast_name")`. Proteomics: N slots sharing one `eds`. Genomics: one slot per loaded dataset.
- **Mixed-omics policy**: type-pure. First load establishes `omics_class`; subsequent loads with different `omics_class` are rejected with `showNotification(type="error")`.
- **DE column names in module code**: canonical (`effect_log2fc` / `score_pvalue` / `score_fdr`).
- **Proteomics auto-load**: preserved — when `inputDataReactive()` resolves to proteomics, slots auto-populate from `eds$contrasts()`.
- **h5ad scope**: `ExploreDatasetH5AD` ships in this PR.
- **PR shape**: one big PR.

---

## Architectural additions to the package

### 1. `ExploreDatasetH5AD` subclass — `R/ExploreDatasetH5AD.R`

New peer of `ExploreDatasetSE`, reading from an in-memory
`anndataR::AnnData` that conforms to the canonical layout
(`exploreDE/plans/anndata_migration/anndata_spec_2026-04-17.qmd`).

Implements every accessor in the contract (`ExploreDataset.R` §2) by
walking the canonical AnnData slots — see spec §1.5 for the
slot-by-method table. Critical mappings:

- `samples()` ← `adata$obs` (strip `factor_` prefix on factor cols)
- `features()` ← `adata$var` + computed `display_name` from `label_*` fallback chain
- `assay(name)` ← `t(adata$X)` if `name == uns$primary_assay_name`, else `t(adata$layers[[name]])`
- `de(c)` ← `as.data.frame(adata$varm[[paste0("de_", c)]])`
- `contrast_meta(c)` ← `adata$uns$de[[c]]`
- `enrichment(c, kind)` ← `adata$uns$enrichment_native[[c]]` filtered by kind, flattened via `eds_flatten_cp_results()`
- `omics_class()` / `data_type()` ← `adata$uns$omics_class` / `adata$uns$omics_type`
- `backend()` returns the AnnData; `backend_kind()` returns `"h5ad"`

### 2. `validate_anndata_schema()` — `R/validate_anndata_schema.R`

Per spec §4.2. First-cut scope: required `uns` fields
(`schema_version`, `omics_type`), required `varm[de_*]` columns,
factor-column convention. Strict full validation can come later — the
runtime `validate_explore_dataset()` already catches most cross-reference
issues on a live `eds`.

### 3. Factory routing — `R/new_explore_dataset.R`

Extend `new_explore_dataset(input)`:

```r
if (is_h5ad_input(input)) {              # path ending .h5ad, or AnnData
  adata <- if (is.character(input)) anndataR::read_h5ad(input) else input
  validate_anndata_schema(adata)
  return(ExploreDatasetH5AD$new(adata))
}
# existing SE path unchanged
```

### 4. Missing accessor: `eds$secondary_grouping()`

`server-multi-boxplot.R` (and `server-correlation.R`, where it's
already noted as the one legacy read remaining) needs the `param$grouping2Name`
equivalent. Add to base class:

```r
secondary_grouping = function() private$.config$secondary_grouping  # may be NULL
```

Populate in each subclass's `.build_config()` from `param$grouping2Name`.

### 5. Tests

- `tests/testthat/test-subclass-h5ad.R` — H5AD subclass accessors against synthetic AnnData fixtures.
- `tests/testthat/test-cross-backend.R` — run `build_example_anndata.R` on each demo SE, load both ways, assert every accessor returns equivalent results (modulo documented divergences).
- `tests/testthat/test-multi-loader.R` — type-pure enforcement, list(eds, contrast) shape on Browse / Manual / Upload / Proteomics-auto.

Existing `tests/testthat/helper-fixtures.R` gains 1-2 canonical AnnData fixtures (built via the converter functions in `plans/anndata_migration/build_example_anndata.R`).

---

## Multi-* subsystem rewrites

### Loader — `inst/app/server-multi-data-selection.R`

This file is the heart of the change. Full rewrite of the loader logic
(not the colour-picker / rename / reorder logic — those are gStoreList-shape-agnostic).

**Replace** every `convert_genomics_se(se, fullPath)` call with:

```r
eds <- new_explore_dataset(input)   # accepts SE, path, or AnnData
.eds_assert_type_pure(eds, multiDeReactive$results)  # showNotification + return on mismatch

if (eds$omics_class() == "proteomics") {
  for (c in eds$contrasts())
    multiDeReactive$results[[paste0(name, " | ", c)]] <- list(eds = eds, contrast = c)
} else {
  multiDeReactive$results[[name]] <- list(eds = eds,
                                          contrast = eds$primary_contrast())
}
```

Branches affected:

- **Browse Projects** loop (line 280–365 today): swap loader call.
- **Manual Import** loop (line 444–545 today): swap loader call, generalise URL accept (drop `urlDataRoot` whitelist — let `new_explore_dataset()` accept any URL or local path).
- **File Upload** (new): `fileInput` accepting `.rds`, `.qs2`, `.h5ad`. Loader handles file path the same way.
- **Proteomics auto-load** (line 734–840 today): collapses to ~10 lines — when `inputDataReactive()$.eds$omics_class() == "proteomics"`, produce N slots via `for (c in eds$contrasts())`. The 20-field `pseudoObj` construction disappears entirely.

Type-pure enforcement: `.eds_assert_type_pure()` checks the first existing slot's `omics_class()` against the incoming `eds`. Mismatch → `showNotification("Cannot mix omics. Clear existing datasets first.", type='error')`, return without adding.

### UI — `inst/app/ui-tab-multi-data-selection.R`

Add a third tab inside the existing `tabBox`:

```r
tabPanel(
  title = "Upload File", icon = icon("upload"),
  fileInput("multiFileUpload", "Choose a file (.rds, .qs2, .h5ad)",
            accept = c(".rds", ".qs2", ".h5ad"), multiple = TRUE),
  helpText("Drop one or more datasets...")
)
```

Existing Browse Projects + Manual Import tabs unchanged in UI; loader logic flows through the same code path.

### Downstream module rewrites

Same mechanical pattern across all 7 modules. For each:

```r
# Old:
seqAnno <- results[[n]]$seqAnno
counts  <- results[[n]]$countList[[input$assay]]
dataset <- results[[n]]$dataset

# New:
eds      <- results[[n]]$eds
contrast <- results[[n]]$contrast
seqAnno  <- eds$de(contrast)        # canonical cols: effect_log2fc, score_pvalue, score_fdr
counts   <- eds$assay(input$assay)
dataset  <- eds$samples()
```

Per-file specifics:

- **`server-multi-de-table.R`**: replace fallback column-name probing with direct `effect_log2fc`/`score_pvalue`/`score_fdr` reads. Use `eds$feature_display()` for gene name column.
- **`server-multi-boxplot.R`**: proxy `multiBoxplotDataReactive()` now returns `results[[1]]$eds`. Reads of `d$factors`, `d$groupingName`, `d$param$grouping2Name` go through `eds$factor_names()`, `eds$primary_grouping()`, `eds$secondary_grouping()`. The contrast-parsing logic at line 339+ uses `eds$contrast_meta(c)$test_group` / `$reference_group` directly, killing the multi-way `--over--` / `_vs_` / `vs ` regex.
- **`server-multi-tile.R`**: kill the seqAnno / seqAnnoList / degResultsTable fallback chain — `eds$de(contrast)` is unconditional. GO columns still read from `eds$features()` (which keeps them).
- **`server-multi-upset.R`**: enrichment access becomes `eds$enrichment(contrast, "ora")` / `(..., "gsea")`. Returns flat keys (`"BP/upGenes"`, `"BP"`). Update `upsetGeneList` to use flat keys — the dataToPlotHeatmap / dataToPlotBubble lookup tables in the consumer modules (and their UI choices) move to flat keys (`"BP/upGenes"` instead of `"BPupGenes"`).
- **`server-multi-bubble.R`** + **`server-multi-heatmap.R`**: consume `upsetGeneList` with flat keys. Their nested ORA/GSEA accesses (`results[[con]]$ora[[ont]][[sig]]`) become `eds$enrichment(contrast, "ora")[["BP/upGenes"]]` etc.
- **`server-multi-kegg.R`**: `results[[n]]$param$refBuild` → `eds$ref_build()`. `entrez_id` column still lives in `eds$features()`. DE filtering uses canonical column names.

### Removed dead code

- The `convert_genomics_se()` / `convert_proteomics_se()` calls inside multi-* (they remain in the single-dataset path until Phase 2.7).
- The 20-field proteomics `pseudoObj` construction.
- The `urlDataRoot` hardcoded path whitelist.
- The seqAnno / seqAnnoList / degResultsTable fallback chains in multi-tile / multi-de-table.

---

## Spec updates

- Move "Phase 3" status of `ExploreDatasetH5AD` to "delivered" in `internal_api_spec_2026-04-30_thin.qmd` §7.5.
- Add `eds$secondary_grouping()` to the §2 accessor contract table.
- Update §6.3 migration map: note `param$grouping2Name` → `eds$secondary_grouping()`.
- Update spec §3 config schema to list `secondary_grouping`.
- Add to README: BYOD now works for multi-*.

---

## Implementation order within the single PR

1. `R/ExploreDatasetH5AD.R` + `R/validate_anndata_schema.R` + factory extension + `eds$secondary_grouping()` accessor.
2. Generate canonical `.h5ad` fixtures from existing demo SE files via `build_example_anndata.R`.
3. Tests: `test-subclass-h5ad.R`, `test-cross-backend.R`, `test-multi-loader.R`. All must pass before module rewrites begin.
4. Rewrite `server-multi-data-selection.R` (loader + auto-load + type-pure).
5. Add Upload tab to `ui-tab-multi-data-selection.R`.
6. Rewrite the 7 downstream modules — order: de-table, boxplot, kegg, tile, upset, bubble, heatmap (upset before bubble/heatmap because the latter consume `upsetGeneList`).
7. Update spec + bump existing testthat snapshots where DE column names changed.
8. Smoke test in app (see Verification).

---

## Verification

End-to-end manual testing before declaring done:

1. **Existing genomics demo**: Browse Projects → load 2 contrasts → exercise every multi-* tab (de-table / boxplot / tile / upset / bubble / heatmap / kegg). Confirm no visual regressions vs. main.
2. **Existing proteomics demo**: open the app — auto-load should populate multi-* slots from `eds$contrasts()`. Exercise every multi-* tab (the ones that apply to proteomics).
3. **BYOD path — h5ad**: produce a canonical h5ad via `build_example_anndata.R` from one of the demo SEs. Upload via the new file widget. Confirm loads, populates correctly, all tabs work.
4. **BYOD path — local .rds/.qs2**: same as 3 but with an SE file uploaded from disk rather than a gStore path.
5. **Type-pure rejection**: load genomics, then try to load proteomics. Confirm error notification, existing slots preserved.
6. **Mixed loading + reordering**: load 3 datasets, reorder via drag bucket, rename, remove one. Confirm slot order propagates to every downstream module (most read `names(multiDeReactive$results)`).

Test suite:

```bash
R -e 'devtools::test()'    # whole testthat suite must pass
R -e 'devtools::check()'   # R CMD check clean
```

Expect test count to rise from current 951 by ~30 (h5ad subclass + cross-backend + multi-loader).

---

## What this PR does NOT include

- Migration of any other Phase 2 server module (volcano, DEGTable, deQc, heatmap2, ordination, gsea2, hypergeo2, summary, initInputData). Those continue under the legacy shim.
- Removal of the legacy shim (`legacy_gStoreList_from_eds()`). It's retired in Phase 2.7 only after every single-dataset server module migrates.
- CSV → h5ad helper (Phase 4a per spec).
- Friendly schema-validator error UI (Phase 4a per spec).
- Producer-side `.h5ad` emission (Phase 4b per spec — separate workstream).
