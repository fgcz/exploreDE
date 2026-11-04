# exploreDE: Reading APB results natively (2026-10-30)

Goal: exploreDE opens a raw APB result directly, as produced by apb2 or APB Export's plain `.h5ad` output (not the `explorede` target). It recognises the file by `uns["apb"]`, and Witold's apbCatalogR reads it and finds columns by their meaning, not by vendor name.

**Status (17:15):** Witold's agent and ours agree on reading through apbCatalogR ("option 2"). **Peter's sign-off pending.** Discussion: `/scratch/talahon_talk/talahon_talk.md`, entries from 10:15 to 17:02 on 30 Oct. Nothing is implemented in exploreDE yet. apbCatalogR hasn't been installed or run here: Auto mode blocks installing downloaded code, so Peter has to install it or allow the install.

## Tested (17:45, apbCatalogR `2d61d04`, installed by Peter in `/scratch/peterleary/Rlib_apbCatalogR`, Dev/R/4.6.0)

Scripts: `/scratch/peterleary/apb_tests/test_apbCatalogR_exports.R` and `test_apbCatalogR_raw_bundle.R`.

- **Loads without `arrow`:** `library(apbCatalogR)` and reading a `.h5ad` leave `arrow` unloaded.
- **`.h5ad` reader:** All eight `explorede_exports/` files and Witold's job file read in 0.6–1.7 s, with SDRF metadata. Because they're `explorede` exports (layers `X` / `intensities` / `intensities_raw`, columns `feature_id` / `factor_*`), every catalogue lookup returns `NULL`. That's correct for them, and it's why canonical detection must come first.
- **Raw results** (the synthetic Parquet results from the bundle's `results.tar.gz`, 2 × 3 each; the archive is no longer in the bundle at `2d61d04`):
  - **Abundance lookups by tool:**

    | Result | `intensity` | `intensity`, `none` | `maxlfq` |
    |---|---|---|---|
    | DIA-NN 2 protein | `NULL` | `NULL` | `PG_MaxLFQ` |
    | DIA-NN 2 precursor | `UnresolvedField` (two layers) | `Precursor_Quantity` | — |
    | DIA-NN 1.8 protein | `UnresolvedField` (two layers) | `PG_Quantity` | `PG_MaxLFQ` |
    | MaxQuant protein | `Intensity` | `Intensity` | `LFQ_Intensity` |
    | FragPipe precursor | `Intensity` | `UnresolvedField` (normalisation unknown) | — |
    | Spectronaut protein | `UnresolvedField` (kind unknown) | — | `UnresolvedField` |

    Asking for a specific kind on Spectronaut protein raises an error rather than returning `NULL`.
  - **`apb` set:** `protein_assignment`, `decoy`, and `contaminant` resolve on every reviewed result.
  - **SDRF:** `obs(level, "condition")` is **ambiguous** when the SDRF has more than one `factor_value_*` column (`edge-sdrf`: treatment and dose). It resolves for prolfquapp.
  - **Unreviewed rule** (`edge-unreviewed`, no software version in `rule_json`): Every lookup raises `UnresolvedField`. The file's own records are still there: `roles$columns` (`protein_assignment = Protein_IDs`, `fasta_accessions = Proteins`) and the `apb_Decoy` / `apb_Contaminant` columns.
- **Consequences for the converter:**
  - Default to `X`, so the default view needs no catalogue lookup.
  - Wrap every catalogue lookup in `tryCatch`, with fallback to the file's own `roles` and `apb_*` columns.
  - Take factors from all `factor_value_*` / `group*` / `bait*` obs columns directly, as the `explorede` target does, not from the single `condition` lookup.
- **Still untested:** Real raw apb2 `.h5ad` files (Witold's question 8).

## Decision: Read through apbCatalogR, not our own parser

- apbCatalogR `2d61d04` (GitHub `anndata-omics-bridge/apbCatalogR`, `main`): `read_parsed_levels(path)` takes a `.h5ad` file, a `.duckdb` file, or a `.parquet` directory and returns a `ParsedLevels`. For `.h5ad` (`R/h5ad.R`), it maps apb2's layout itself: The level record at `uns[<level>]$apb`, the root record at `uns$apb`, `storage` as JSON text, `json_values`, and layers stored under `physical_name`. It accepts only `apb2-parsed-levels` format version 8, and passes `sdrf` / `fasta` through as metadata, so SDRF factors survive.
- **Dependencies:** Imports `anndataR (>= 1.2.2)`, `jsonlite`, `Matrix`, `R6`, `stats`, `tools`, `utils`, all in `fgcz_base_image`. `arrow` (Parquet only) and `duckdb` / `DBI` are in Suggests. Witold's agent reports `make check` clean, and that neither loading the package nor reading a `.h5ad` loads `arrow`. I checked DESCRIPTION, NAMESPACE, and `read_table()`.
- **Why not our own parser:** The layout knowledge stays with apb2's authors, and exploreDE gets kinds and normalisation from the catalogue. Our parser's advantage, no new packages in the base image, also holds for apbCatalogR now (see Deployment).
- **Not yet:** `.h5mu`, a release tag, and the compatibility rule for files from a newer apb2.

## Which layer exploreDE shows

Agreed with Witold's agent: **`X`, the rule's primary layer, by default.** Every rule has one, and it is what the `explorede` export writes.

```r
r <- read_parsed_levels("job.h5ad")
lvl <- r$levels$protein
t(lvl$X)                                # proteins x samples
lvl$uns$apb$storage$primary_layer       # its name, e.g. PG_MaxLFQ
```

- `X` means different things across tools: DIA-NN's is normalised (`Precursor_Normalised`, `PG_MaxLFQ`), MaxQuant's is raw `Intensity`. **Show the layer's name in the app** (Results Summary) and record it in `uns["provenance"]`.
- **For specific requests** (e.g. unnormalised precursors), go through the catalogue with a kind, plus a normalisation where a tool has several:
  ```r
  abundance <- DataCatalog$new(r, "abundance")
  abundance$kinds("peptidoform_ion", "layers")                            # "intensity" "ms1_intensity"
  abundance$layer("peptidoform_ion", "intensity", normalization = "none")  # Precursor_Quantity
  ```
- At protein level, no single kind fits every tool: MaxQuant `intensity`, DIA-NN `maxlfq` (re-catalogued in apb-catalog `811caaa`), AlphaDIA `directlfq`, Spectronaut and MSAngel `unknown`. That's why `X` is the default.
- As in the `explorede` target: Take log2 and turn ≤ 0 into NaN, keep linear values as `intensities_raw`.

## Implementation plan in exploreDE

1. **Detection**, in `new_explore_dataset()` (`R/new_explore_dataset.R:58`), after `.h5ad_read()`:
   1. prolfquapp (`is_prolfquapp_anndata()`)
   2. canonical: `uns[["schema_version"]] == "exploreDE/1.0"`
   3. raw APB: `uns[["apb"]][["storage"]]` parses as JSON with `format == "apb2-parsed-levels"`, and there is no `schema_version`
   4. otherwise, canonical as today

   Canonical must come before raw, because `explorede` files carry the same `uns$apb` marker. The `[[` matters: `$` partial-matches (the 3.9.0 `uns$de` bug).
2. **Read:** `apbCatalogR::read_parsed_levels(path)`. It needs a path, not the `adata` already read, so the file is read twice; ask Witold for an `AnnData` entry point if that's slow on big files.
3. **Convert to canonical in memory:** Build an `exploreDE/1.0` `AnnData` from the level and pass it to `ExploreDatasetH5AD$new()`, which already accepts an `AbstractAnnData`. The rest of the app stays unchanged. Mapping, as in the `explorede` target:
   - **Level:** protein → `protein`, peptide → `peptide`, `peptidoform_ion` → `precursor` (`feature_id_kind`, `hierarchy`). A single-level `.h5ad` holds exactly one.
   - **`feature_id`:** The `protein_assignment` role (catalogue `apb`). At peptide and precursor level, `_` plus the peptide / ion id.
   - **`label_gene_symbol`:** FASTA gene name (`varm["fasta_validation"]` or `varm["fasta"]`), else the leading id, else `protein_assignment`; at peptide and precursor level, with the peptide / ion appended. `description` comes from the FASTA columns.
   - **`obs`:** The obs keys become `sample_id` / `sample_label`. `factor_value_<x>` → `factor_<x>`, and `group*` / `bait*` → `factor_<name>`, with names in `[a-z0-9_]` and sorted categories.
   - **Decoys:** Dropped (`apb_Decoy` / catalogue `apb` kind `decoy`).
   - **Contaminants:** Kept and flagged, waiting on Witold's question 4 for the column (`is_contaminant` proposed, to go into the spec).
   - **`uns`:** `schema_version`, `omics_type = "proteomics"`, `de = {}`, `enrichment = {}`, `enrichment_support = {}`, and `provenance` with `source_format = "apb2"` plus the source layer and its normalisation. Keep `uns["apb"]`.
4. **apbCatalogR optional:** Put it in exploreDE's Suggests. Without it, a raw APB file gets a clear message: "Reading raw APB files needs apbCatalogR; or export with APB Export's exploreDE output."
5. **Tests:** Raw files against the `explorede` exports from the same jobs (Witold's question 8). They should have the same `feature_id`s, factors, and values for the same layer.

## Deployment

- apbCatalogR has no compiled code, and everything it imports is in `fgcz_base_image`. So install it into a library next to exploreDE rather than into the base image, e.g. `exploreDE/rlib/` (git-ignored), with `.libPaths()` extended at the top of `inst/app/app.R`. Peter can then update it on each Witold release without a rebuild by José Luis.
- Pin a commit (`2d61d04` today) until Witold tags releases. Every apb2 format bump needs a new apbCatalogR, because the reader accepts only format version 8.

## What the `explorede` target does today (findings, 30 Oct)

- **Values:** `X` / `intensities` = log2 of the primary layer (the same choice as above). The delivered DIA-NN and Spectronaut files are vendor-normalised. The layer used isn't recorded in the file; asked for in `uns["provenance"]`.
- **Decoys are kept:** 334 of 6,896 MaxQuant protein groups are led by `REV__`, and 421 of 46,560 peptides. Contaminants are kept, with no column marking them. Neither the spec nor the reader has a decoy or contaminant field.
- **Version label:** `target_version_pattern` is `^3\.(8|9)\.` since apb-export `97ebda1`.

## Open with Witold

- Release tag, and the compatibility rule for newer apb2 files.
- Question 4: Drop decoys in `explorede`, and the contaminant column.
- Question 5: Are the FASTA column names stable?
- Question 6: `.h5mu` (later); does APB Export serve Parquet?
- Question 7: Does a format bump always mean a layout change?
- Question 8: Raw test files in `/scratch/talahon_talk/apb_raw/`.

## Source

- Read-only copies on fgcz-c-051, in `/scratch/peterleary/apb_review_src/`: apbCatalogR (`HEAD` of the morning, `ad895fa`, `2d61d04`), apb-catalog, apb-export, and apb2. Delete when done.
