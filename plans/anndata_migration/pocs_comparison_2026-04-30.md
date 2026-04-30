# POC comparison — 2026-04-30

Two POCs were built side-by-side to validate the storage-agnostic API
design from `internal_api_spec_2026-04-30_thin.qmd`. Both load the same
seven datasets through the same accessor contract; they differ in *how
the legacy SE shapes are reconciled* with the new API.

| | Subclass POC | Canonical POC ⭐ |
|---|---|---|
| **File** | `poc_explore_dataset_2026-04-30.R` | `poc_canonical_2026-04-30.R` |
| **Lines** | 1051 | 883 |
| **Class hierarchy** | 5 classes (1 abstract, 1 SE intermediate, 4 concrete SE subclasses, +h5ad later) | 2 classes (1 abstract, 1 concrete reader) |
| **SE shape handling** | Per-shape subclass with `.build_config` + `.prep` overrides | One per-shape converter function producing canonical AnnData |
| **AnnData support** | Future `ExploreDatasetH5AD` subclass | Already present — same reader |
| **`.h5ad` export** | Separate work item (Phase 3) | Free — `anndataR::write_h5ad(eds$backend(), …)` |
| **Spec compliance** | Subclasses produce per-subclass-shape data | Canonical layout enforced at conversion time |
| **External deps** | R6, SummarizedExperiment | + anndataR |

**Recommendation: canonical POC.**

## Why canonical wins

### 1. Fewer concepts in the reader

The reader (`ExploreDatasetH5AD`) reads from a uniform layout. There's no
"is this proteomics, do I need to remap rownames?" branching anywhere
in any read method. All shape-specific knowledge lives in one place
(`se_to_canonical_anndata.R`), runs once at load, and produces the same
in-memory object regardless of source.

The subclass POC distributes shape knowledge across four subclasses;
each one is small, but the contributor needs to read four files to
understand what "proteomics features" really means.

### 2. The proteomics flow is materially better

This is the empirical finding from running both POCs on the same data:

| Surface | Subclass POC | Canonical POC |
|---|---|---|
| `eds$de(c)` for proteomics | empty (off-by-prefix bug) | populated |
| `eds$features()$label_gene_symbol` for proteomics | not exposed | `NUDT4B`, `GATD3B`, … |
| `eds$features()$label_uniprot_id` for proteomics | not exposed | populated |
| Phospho `feature_id_kind` | `"protein"` (wrong) | `"site"` (correct) |
| Proteomics `description` (full UniProt) | not exposed | full descriptions |
| Proteomics gene/protein dropdowns (`proteomics_choices`) | needed special-case method | derived from `var$label_*` columns uniformly |

The canonical converter is forced to lay out *all* proteomics data in a
documented layout, which surfaces issues the subclass POC's ad-hoc
extraction had silently dropped.

### 3. The migration's Phase 3 (SE → .h5ad export) drops out for free

The migration plan calls for an SE → `.h5ad` exporter for public-data
release. With the canonical approach, this is one line:

```r
anndataR::write_h5ad(eds$backend(), "/path/to/output.h5ad")
```

…because `eds$backend()` is already the canonical AnnData. With the
subclass approach, Phase 3 is a separate ~350-line work item.

### 4. New backends don't fan out

If we ever need a Parquet or Zarr backend, the subclass tree would be
tempted to grow `ExploreDatasetParquet_RNASeq`, `_CountQC`, etc.
(repeating SE's legacy fan-out). With the canonical approach, the new
backend writes a converter — `parquet_to_canonical_anndata()` — and the
existing reader handles it. One new file per backend, ever.

### 5. Spec compliance is mechanically enforced

Every SE that enters the system passes through the canonical
normaliser. If the layout is wrong, conversion fails loudly. With the
subclass approach, each subclass interprets the SE differently and
nothing checks that the result conforms to a documented schema —
schema compliance is a code-review concern.

## The one downside (minor, irrelevant in practice)

**Brief peak memory during loading.** During `new_explore_dataset(path)`
for an SE-derived dataset, the SE and the converted AnnData both exist
in memory until the SE is dropped (which happens when
`se_to_canonical_anndata()` returns and the local `se` variable goes
out of scope). For a 170k-feature isoform-level dataset this is a
short-lived ~1-2 GB peak.

In practice: irrelevant. Our biggest dataset is the isoform RNA-seq
(170k features × 12 samples) — peak at ~1 GB extra, lasting maybe a
second. The Shiny worker has many GB of headroom.

If it ever became a concern, the converter could nullify the SE
reference internally as it builds the AnnData, freeing memory before
returning. ~5 lines of change.

## Other things flagged but not problems

- **anndataR dependency.** It's pure R, no compiled code, in CRAN, two
  small dependencies. Trivial to install. The migration is committing
  to `.h5ad` as the primary format anyway.
- **Spec coupling.** The canonical layout becomes the in-memory contract.
  We *want* this — it's the migration's whole point. If the spec
  changes, the converter changes; modules don't notice.
- **Round-trip fidelity.** Verified: same DE values, same feature
  counts, same enrichment objects on the teaching dataset. The
  canonical approach actually produced *more* correct data than the
  subclass POC (proteomics DE bug).

## What changes in the implementation plan

The migration plan in `internal_api_spec_2026-04-30_thin.qmd` §10
remains valid in shape. With the canonical approach, the per-phase
content shifts:

| Phase | Subclass approach | Canonical approach |
|---|---|---|
| 0 | Abstract base + factory shell + validators | Same |
| 1 | `ExploreDatasetSE_RNASeq` + SUSHI config builder | `se_to_canonical_anndata()` for RNA-seq + `ExploreDatasetH5AD` reader |
| (Phase 1 add) | `ExploreDatasetSE_CountQC` + `_CrisprCountQC` + `_Proteomics` | `se_to_canon_count_qc()` + `_crispr_count_qc()` + `_proteomics()` |
| 2 | Module migration (unchanged) | Module migration (unchanged) |
| 3 | New `ExploreDatasetH5AD` subclass | **Already done in Phase 1** — Phase 3 becomes "ship the .h5ad export utility (`anndataR::write_h5ad`) and remove the SE adapter step" |
| 4 | BYOD config builders | BYOD = user supplies a canonical `.h5ad`; converter not needed |

Net effort: similar total lines of code; same migration sequence;
Phase 3 collapses into Phase 1.

## Decision

Going forward with the canonical POC as the primary implementation
target. The subclass POC is preserved for reference (it documents the
alternative design we considered) and as a fallback in case the
AnnData direction is ever reconsidered.
