# Producer-side migration to canonical AnnData

**Status**: not yet pursued. Stub doc, 2026-04-30.

This is the upstream-producer arm of the AnnData migration (front 1 of
the three-fronts framing in `README.md`). Out of scope for the
exploreDE consumer-side spec
(`internal_api_spec_2026-04-30_thin.qmd`), tracked here separately
because:

- The stakeholders are different teams (SUSHI / NGS, FGCZ proteomics)
  rather than just exploreDE / Peter.
- The timeline is multi-quarter and depends on those teams' priorities.
- The relevant code lives in *their* repositories, not in exploreDE.
- exploreDE works fine without this — the per-shape SE subclasses in
  `poc_explore_dataset_2026-04-30.R` read the legacy SEs directly, and
  the SE→canonical AnnData mapping is documented in
  `build_example_anndata.R`. This producer work just makes the SE
  bridge obsolete.

## Producers in scope

### 1. SUSHI ezApps (FGCZ NGS team)

Currently emit `deResult.rds` / `deResult.qs2` (RNA-seq DESeq2),
`rawData_final.rds` / `.qs2` (CountQC, CrisprCountQC) on gStore.
Should emit canonical `.h5ad` (per
`anndata_spec_2026-04-17.qmd`) alongside or instead of the SE.

Affected ezApps (rough inventory — confirm with NGS team):
- DESeq2 (RNA-seq, all variants — single contrast, multi-factor,
  isoform)
- ExploreMageckCounts (CRISPR Mageck CountQC)
- Any other DESeq2-shaped or CountQC-shaped ezApp that produces an
  exploreDE-readable SE today.

### 2. FGCZ proteomics workflow (B-Fabric / pStore)

Currently emits `<id>.rds` (prolfqua-derived SE with nested
`constrast_*` rowData) on pStore. Should emit canonical `.h5ad`.

Affected workflows:
- Protein-level prolfqua
- Phospho-site prolfqua
- Peptide-level prolfqua
- Any future workflow targeting exploreDE.

## Why coordinate this

- **Eliminates the SE → canonical AnnData converter step** at
  exploreDE load time (small but real performance win, more
  importantly removes a maintenance burden as legacy SE shapes
  evolve).
- **Producer teams own the format conversion**, which is the right
  place for it — they understand their own pipelines best and can
  incorporate the conversion at the appropriate stage.
- **The canonical `.h5ad` is portable** — Python users, downstream
  ML pipelines, and external collaborators can consume it directly
  without going through exploreDE's R idioms.
- **One archive format for all FGCZ omics public-data releases**.

## Reference implementation

`build_example_anndata.R` is the reference: its
`convert_rnaseq()` / `convert_crispr_countqc()` / `convert_proteomics()`
functions produce the canonical AnnData layout from each of the legacy
SE shapes, with every SE-slot → AnnData-slot mapping narrated inline.
Producer teams can:

1. Run their existing pipeline as-is, producing the SE.
2. Pass the SE through the appropriate `convert_*()` function.
3. Write the result via `anndataR::write_h5ad()` to the same gStore /
   pStore path alongside (or instead of) the `.rds`.

That's the simplest adoption path and requires no changes inside the
pipelines themselves. Eventually, producers can choose to build the
AnnData directly without going through an SE.

## When this is done

The exploreDE-side cleanup (per
`internal_api_spec_2026-04-30_thin.qmd`):

- Retire the per-subclass `ExploreDatasetSE_*` classes (only the
  `ExploreDatasetH5AD` subclass remains).
- Drop `convert_genomics_se()` / `convert_proteomics_se()` legacy
  entry points.
- Remove `SummarizedExperiment` from `Imports:`.
- exploreDE talks `.h5ad` only.

## Status

Not started. Open this document to discuss with the relevant teams
when there's bandwidth on both sides. Until then, exploreDE handles
the legacy SEs via the converters and works fine.
