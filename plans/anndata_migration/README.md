# exploreDE AnnData migration

Started 2026-04-17 with the on-disk schema; expanded 2026-04-30 with the
in-memory API + module migration plan.

## What this is

Moving exploreDE off the current SUSHI / FGCZ-proteomics
SummarizedExperiment outputs and onto a single canonical AnnData
(`.h5ad`) layout. Two intertwined wins:

1. **Internal harmonisation** — one accessor API for genomics,
   proteomics, future metabolomics; one in-memory class; one on-disk
   format.
2. **Externalisation** — exploreDE stops being FGCZ-internal.
   Researchers anywhere can point the app at their own data
   (canonical `.h5ad` they built themselves, or counts CSV + sample
   sheet via a helper) and get the full app experience without any
   FGCZ infrastructure.

This affects three independent pieces, all designed together:

| Front | Today | After |
|---|---|---|
| **Upstream producers** | SUSHI ezApps emit `deResult.rds` / `rawData_final.rds` (genomics, varied shapes); FGCZ proteomics workflow emits `<id>.rds` on pStore (different shape again, nested `constrast_*` rowData) | All producers emit canonical `.h5ad` per `anndata_spec_2026-04-17.qmd` |
| **Canonical format** | — (no shared contract) | `anndata_spec_2026-04-17.qmd` defines the on-disk shape: `obs` / `var` / `X` / `layers` / `varm["de_<contrast>"]` / `uns["enrichment_native"]` / `uns["omics_type"]` etc. |
| **Downstream consumer** (exploreDE) | Reads SE shapes via `convert_genomics_se()` / `convert_proteomics_se()` into a `gStoreList` plain list; modules read fields directly | Reads canonical `.h5ad` into a single `ExploreDatasetH5AD` R6 instance; modules use the accessor API (`eds$samples()`, `eds$de(c)`, etc.) |

The **converters** (`build_example_anndata.R`, `se_to_canonical_anndata()`
inside `poc_canonical_2026-04-30.R`) are the transitional bridge: they
turn today's SEs into canonical AnnData so the new exploreDE works
*now* even before upstream producers migrate. As producers adopt the
canonical format directly, the converters become obsolete — which is
the win condition.

## Status (2026-04-30)

- **On-disk format**: spec drafted, validated against six real
  datasets via the canonical POC's converter (the converter wrote
  every example dataset successfully into in-memory AnnData
  conforming to the spec).
- **In-memory API**: design complete, two POCs built, comparison
  done; **canonical-converter approach chosen** — see `pocs_comparison_2026-04-30.md`.
- **Module migration plan**: full per-module audit (29 files) in
  `module_inventory_2026-04-30.md`. Phased migration described in
  `internal_api_spec_2026-04-30_thin.qmd` §10 (chosen plan,
  authoritative), with canonical-design details in §15.
- **Implementation**: not started. Phase 0 is next — lift the
  canonical POC into `R/` package code, add validators + tests.
- **Upstream producer migration**: separate workstream tracked in
  `producer_migration_2026-04-30.md`. Not yet pursued; consumer-side
  works without it via the SE converters.

## Artefacts in this directory

### On-disk format (front 2)

- **`anndata_spec_2026-04-17.qmd`** — the canonical AnnData on-disk
  schema. Required slots, optional slots, naming conventions
  (`factor_*`, `label_*`, `hierarchy_*`, `de_<contrast>`,
  `enrichment_native`), per-omics specifics (RNA-seq / CRISPR /
  proteomics / metabolomics+lipidomics future). Every front references
  this.

### Converter (the bridge)

- **`build_example_anndata.R`** — earlier (2026-04-17) verbose
  converter exploring the SE → AnnData mapping for the canonical
  example datasets. Predates the in-memory API work; informed it.
- **`se_to_canonical_anndata()`** — inline in
  `poc_canonical_2026-04-30.R`. Production-shaped converter covering
  all four legacy SE shapes (RNA-seq, CountQC, CrisprCountQC,
  proteomics). Also doubles as the SE → `.h5ad` export pipeline:
  `anndataR::write_h5ad(eds$backend(), path)`.

### In-memory API + consumer (front 3)

- **`internal_api_spec_2026-04-30_thin.qmd`** ⭐ — the chosen design
  spec (1367 lines). Describes the abstract `ExploreDataset` class,
  the accessor contract, the migration plan, and (in §15) the
  canonical-converter implementation strategy that won. Read this
  first if you want to know what's getting built.
- **`internal_api_spec_2026-04-30.qmd`** — Draft 1, the materialised-
  internal-object alternative. Documented but not chosen; preserved
  for context.
- **`poc_canonical_2026-04-30.R`** ⭐ — the chosen POC (883 lines).
  One `ExploreDatasetH5AD` reader + four per-shape converter
  functions producing in-memory `anndataR::AnnData`. Tested on six
  real datasets.
- **`poc_explore_dataset_2026-04-30.R`** — alternative POC (1051
  lines): four SE subclasses + virtual h5ad. Works but had bugs in
  the proteomics path that the canonical converter exposed and
  fixed. Preserved for reference and as a fallback if the AnnData
  direction is reconsidered.
- **`pocs_comparison_2026-04-30.md`** — head-to-head comparison of
  the two POCs with the architectural argument and the one minor
  memory caveat. Decision rationale lives here.
- **`module_inventory_2026-04-30.md`** + `.html` — full per-module
  audit of exploreDE (29 files). Drives the Phase 2 PR plan: which
  module migrates when, with what risk, against what API surface.

### Upstream producers (front 1)

- **`producer_migration_2026-04-30.md`** — stub doc tracking the
  separate workstream of getting SUSHI ezApps + FGCZ proteomics
  workflow to emit canonical `.h5ad` directly. Not yet pursued; tracked
  here so it doesn't get lost.

## How to read this

If you're new to the project, read in this order:

1. This README (you're here).
2. `pocs_comparison_2026-04-30.md` — quick orientation on the two
   approaches and what won.
3. `internal_api_spec_2026-04-30_thin.qmd` — start at the ⚠ callout
   at the top, then jump to §15. The earlier sections (§1–§14)
   describe the four-subclass alternative we considered.
4. `poc_canonical_2026-04-30.R` — the working code.
5. `anndata_spec_2026-04-17.qmd` — when you need the on-disk format
   details.
6. `module_inventory_2026-04-30.md` — when you're ready to plan the
   per-module migration PRs.

## What lands first (Phase 0)

Per the spec §10.2 / §13:

- `DESCRIPTION` — add `anndataR`, `R6` to Imports
- `R/ExploreDataset.R` — abstract base (lifted from POC)
- `R/ExploreDatasetH5AD.R` — concrete reader (lifted from POC)
- `R/se_to_canonical_anndata.R` — converter (lifted from POC)
- `R/new_explore_dataset.R` — factory (lifted from POC)
- `R/utils_explore_dataset.R` — shared helpers
- `R/validate_explore_dataset.R` — runtime validator
- `R/validate_anndata_schema.R` — `.h5ad` schema validator
- `tests/testthat/test-converter-*.R`, `test-reader-shape.R`,
  `test-round-trip.R`, `test-validator-*.R`

Phase 0 is **not stubs** — it ships working code from the POC,
package-quality, with tests covering the six demo datasets. The
running Shiny app is unchanged (the new code isn't wired in yet —
that's Phase 1's job, via the legacy `gStoreList` shim).
