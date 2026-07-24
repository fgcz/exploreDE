# exploreDE AnnData migration

Started 2026-04-17 with the on-disk schema; expanded 2026-04-30 with
the in-memory API + module migration plan.

## What this is

Moving exploreDE off the current SUSHI / FGCZ-proteomics
SummarizedExperiment outputs and onto a single canonical AnnData
(`.h5ad`) layout. Two intertwined wins:

1. **Internal harmonisation** — one accessor API for genomics,
   proteomics, future metabolomics; one on-disk format.
2. **Externalisation** — exploreDE stops being FGCZ-internal.
   Researchers anywhere can point the app at their own data (canonical
   `.h5ad` they built themselves, or counts CSV + sample sheet via a
   helper) and get the full app experience without any FGCZ
   infrastructure.

This affects three independent pieces, all designed together:

| Front | Today | After |
|---|---|---|
| **Upstream producers** | SUSHI ezApps emit `deResult.rds` / `rawData_final.rds` (genomics, varied shapes); FGCZ proteomics workflow emits `<id>.rds` on pStore (different shape again, nested `constrast_*` rowData) | All producers emit canonical `.h5ad` per `anndata_spec_2026-04-17.qmd` |
| **Canonical format** | — (no shared contract) | `anndata_spec_2026-04-17.qmd` defines the on-disk shape: `obs` / `var` / `X` / `layers` / `varm["de_<contrast>"]` / `uns["enrichment_native"]` / `uns["omics_type"]` etc. |
| **Downstream consumer** (exploreDE) | Reads SE shapes via `convert_genomics_se()` / `convert_proteomics_se()` into a `gStoreList` plain list; modules read fields directly | Reads through a small R6 class hierarchy. Today: per-shape `ExploreDatasetSE_*` subclasses wrapping the legacy SEs. Future: a single `ExploreDatasetH5AD` subclass once producers emit canonical `.h5ad`. Modules call accessors (`eds$samples()`, `eds$de(c)`, …) — same API regardless of backend. |

## Design choice (2026-04-30)

Per-shape SE subclasses now, single h5ad subclass later. Class
hierarchy:

```
ExploreDataset                (abstract base)
├── ExploreDatasetSE          (abstract intermediate)
│   ├── ExploreDatasetSE_RNASeq        (DESeq2/EdgeR/Limma)
│   ├── ExploreDatasetSE_CountQC
│   ├── ExploreDatasetSE_CrisprCountQC (inherits from _CountQC; sgrna_to_gene flag)
│   └── ExploreDatasetSE_Proteomics
└── ExploreDatasetH5AD        (peer; lands Phase 3, one class for every omics)
```

A four-SE-subclass POC validates the design end-to-end on six real
datasets; see `poc_explore_dataset_2026-04-30.R`.

## Status (2026-05-15)

- **On-disk format**: spec drafted, validated against six real
  datasets via the converter in `build_example_anndata.R` (every
  example dataset successfully converted to in-memory AnnData
  conforming to the spec).
- **In-memory API**: design complete, subclass POC working on the six
  demo datasets.
- **Module migration plan**: full per-module audit (29 files) in
  `module_inventory_2026-04-30.md`. Phased migration described in
  `internal_api_spec_2026-04-30_thin.qmd`.
- **Implementation**: not started. Phase 0 is next — lift the
  subclass POC into `R/` package code, add validators + tests.
- **Upstream producer migration**: separate workstream tracked in
  `producer_migration_2026-04-30.md`. Not yet pursued; consumer-side
  works without it via the in-memory SE subclasses.

## Artefacts in this directory

### On-disk format (front 2)

- **`anndata_spec_2026-04-17.qmd`** — the canonical AnnData on-disk
  schema. Required slots, optional slots, naming conventions
  (`factor_*`, `label_*`, `hierarchy_*`, `de_<contrast>`,
  `enrichment_native`), per-omics specifics (RNA-seq / CRISPR /
  proteomics / metabolomics+lipidomics future). Every front
  references this. **Frozen.**

### SE → AnnData converter (the bridge / documentation)

- **`build_example_anndata.R`** — narrated SE → canonical AnnData
  converter. `convert_rnaseq()` / `convert_crispr_countqc()` /
  `convert_proteomics()` each map every SE slot to its AnnData
  destination with inline `.mapline()` narration. Doubles as the
  reference for upstream producers who want to emit `.h5ad` from
  their existing SE pipelines. Reading this file IS part of the
  documentation.

### In-memory API + consumer (front 3)

- **`internal_api_spec_2026-04-30_thin.qmd`** ⭐ — the design spec.
  Describes the abstract `ExploreDataset` class, the per-shape SE
  subclass hierarchy, the accessor contract, and the migration plan.
  Read this first if you want to know what's getting built.
- **`poc_explore_dataset_2026-04-30.R`** ⭐ — working POC of the
  subclass design (1051 lines, 5 classes). Tested on six real
  datasets. Phase 0 lifts this into `R/` package code.
- **`module_inventory_2026-04-30.md`** — full per-module audit of
  exploreDE (29 files). Drives the Phase 2 PR plan: which module
  migrates when, with what risk, against what API surface.

### Upstream producers (front 1)

- **`producer_migration_2026-04-30.md`** — stub doc tracking the
  separate workstream of getting SUSHI ezApps + FGCZ proteomics
  workflow to emit canonical `.h5ad` directly. Not yet pursued;
  tracked here so it doesn't get lost.

## How to read this

If you're new to the project, read in this order:

1. This README (you're here).
2. `internal_api_spec_2026-04-30_thin.qmd` — the design spec.
3. `poc_explore_dataset_2026-04-30.R` — the working code.
4. `anndata_spec_2026-04-17.qmd` — when you need on-disk format
   details.
5. `build_example_anndata.R` — when you need to see the SE → AnnData
   mapping spelled out.
6. `module_inventory_2026-04-30.md` — when you're ready to plan the
   per-module migration PRs.

## What lands first (Phase 0)

Files added to the exploreDE R package, lifted from the subclass POC:

- `DESCRIPTION` — add `anndataR`, `R6` to `Imports:`
- `R/ExploreDataset.R` — abstract base
- `R/ExploreDatasetSE.R` — abstract intermediate
- `R/ExploreDatasetSE_RNASeq.R` — concrete RNA-seq subclass
- `R/ExploreDatasetSE_CountQC.R` — concrete CountQC subclass
- `R/ExploreDatasetSE_CrisprCountQC.R` — concrete CRISPR CountQC subclass (inherits from `_CountQC`)
- `R/ExploreDatasetSE_Proteomics.R` — concrete proteomics subclass
- `R/new_explore_dataset.R` — factory (detects SE shape, picks the right subclass)
- `R/utils_explore_dataset.R` — shared helpers (`build_feature_id`, `standardise_de_columns`, etc.)
- `R/validate_explore_dataset.R` — runtime validator on a live `eds`
- `tests/testthat/test-subclass-*.R`, `test-factory.R`,
  `test-validator-*.R` against the six demo datasets

Phase 0 is **not stubs** — it ships working code from the POC,
package-quality, with tests covering the six demo datasets. The
running Shiny app is unchanged (the new code isn't wired in yet —
that's Phase 1's job, via the legacy `gStoreList` shim).
