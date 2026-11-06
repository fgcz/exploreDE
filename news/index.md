# Changelog

## exploreDE 3.10.0

- Security: an `.h5ad` can no longer run code when its pathway tabs
  open. The reader decoded `uns$enrichment_native` with
  [`unserialize()`](https://rdrr.io/r/base/serialize.html), which runs
  code embedded in the file; that block is now ignored, and
  [`export_h5ad_to_parquet()`](https://fgcz.github.io/exploreDE/reference/export_h5ad_to_parquet.md)
  no longer copies it. Files written before 2026-09-25 that carry their
  ORA/GSEA results only in that block show no enrichment; rebuild them
  with
  [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md).
- Security: a `?data=` link whose path contains “zip” no longer loads a
  fixed 2023 proteomics dataset over `ssh` before the access check; that
  leftover debugging branch is removed.
- Security: Multi DE “Browse Projects” lists a project’s analysis
  folders only for users with access to it. A link the user was denied
  no longer leaves its project browsable after loading the demo or an
  upload.
- Security: the Multi DE Upset/Venn tab shows the features of a selected
  intersection as text. It inserted the feature ids as HTML, so a
  dataset whose ids contained markup could run script in the viewer’s
  browser. The “Analyse at Enrichr website” link is gone; it called a
  function the app never defined and did nothing.
- [`exploreDE()`](https://fgcz.github.io/exploreDE/reference/exploreDE.md)
  no longer turns on
  [`shiny::reactiveConsole()`](https://rdrr.io/pkg/shiny/man/reactiveConsole.html),
  which stayed on in the calling R session after the app closed and hid
  errors from reading reactive values outside a reactive context.
- [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  aligns every count matrix to the samples by column name. A second
  matrix (for example `tpm`) with its samples in another order was
  stored under the wrong sample labels; a matrix with different samples
  is now an error.
- [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  keeps the values of duplicated feature ids. Every duplicate after the
  first was stored as an all-NA row, although the warning said the ids
  were made unique.
- [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  stops when two count matrices map to the same layer (for example `raw`
  and `counts`) or two contrast names to the same storage key (for
  example `B_vs_A` and `b-vs-a`). Before, one of them was silently
  dropped.
- DE Table: the significance calculator and the “selected” and
  “filtered” downloads always use the selected contrast. Each contrast
  change registered another copy of their handlers, so after switching,
  a copy for an earlier contrast could answer, and selected row numbers
  could be applied to the earlier contrast’s table.
- Breaking: an `.h5ad` must carry `uns` `schema_version`
  (`"exploreDE/1.0"`), `omics_type` and `feature_id_kind`, as the spec
  requires; a file without them is refused with a message naming the
  missing key. Before, defaults were filled in: a proteomics file with a
  missing or misspelt `omics_type` loaded as RNA-seq (VST and
  normalisation computed on log2 intensities), and any AnnData loaded as
  exploreDE/1.0. Files written by
  [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  carry all three.
- Breaking:
  [`new_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/new_explore_dataset.md)
  and the constructors of the `ExploreDataset*` classes no longer take a
  `config` argument. Nothing in the package passed one; for an `.h5ad`
  or a Parquet directory a supplied config made loading fail and skipped
  the check that the file is exploreDE/1.0.
- Breaking:
  [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  requires contrast names of the form `<test>_vs_<reference>`, with both
  groups levels of one `coldata` factor, and records that factor as
  `factor_used`. Before, it stored any name split at `_vs_` as the
  groups (`contrast1` gave the reference group “NA”).
- [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  writes proteomics layers under the spec’s names: a single matrix as
  `intensities`, raw ones (`raw`, `counts`) as `intensities_raw`.
  Before, the single matrix became `counts`. `uns$de` records
  `model = "DESeq2"` only for DESeq2 results, instead of for every
  table, and `enrichment_filter$padj_max` only for pathway tables that
  were filtered.
- [`validate_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/validate_explore_dataset.md)
  also fails on DE tables without `effect_log2fc` or `score_fdr`, and on
  a raw or primary assay that is not offered. With `strict = FALSE` it
  no longer throws when
  [`samples()`](https://rdrr.io/pkg/Biobase/man/class.MIAME.html) fails,
  and records a failing `assay_names()` instead of skipping the check.
- Multi DE Manual Import loads the datasets when all twelve link boxes
  are filled; before, nothing loaded. An import file whose first cell is
  empty no longer ends the session.
- Proteomics results whose `protein_Id` has gaps (SAINT, for example)
  load: a missing id is replaced by the row name. Before, one missing id
  stopped the dataset from loading, and two or more could replace every
  feature id by its row name. A table without a `protein_Id` column is
  refused with a message that lists its columns.
- Batch correction is skipped, with a message, when the batch or design
  column has missing values. A missing value in the design column
  stopped the tab with “number of rows of matrices must match”.
- Ordination (PCA, PCoA, NMDS) leaves out features with missing values
  instead of samples. Proteomics data has missing values in every
  sample, so when “Use top n features” covered all features (the default
  for datasets with fewer than 500), no plot appeared.
- Ordination: FGCZ staff see the “Preprocessing order” settings again;
  everybody saw the staff-only message instead.
- The README says how the `?data=` LDAP check identifies the user
  (`SHINYPROXY_USERNAME`, one ShinyProxy container per user), and that
  without it every gStore link is refused.
- Every backend builds its DE table the same way: SAINT `.h5ad` results
  are sorted by FDR, as SAINT `SummarizedExperiment`s are (they kept
  feature order), and an RNA-seq `SummarizedExperiment` without a
  `usedInTest` column marks as tested the features with a p-value, as
  `.h5ad` files do (it had no such column).
- Multi DE opens canonical proteomics `.h5ad` files and Parquet
  directories: their contrasts are loaded into the Multi DE tabs, as for
  a prolfquapp result; before, the tabs stayed empty. All proteomics
  data now goes through the same converter, so the Multi-Contrast
  Boxplot offers the intensities as `intensities` / `intensities_raw`
  instead of the result’s own assay names (`transformedData`,
  `imputedData`, …).
- The Multi DE “Integrated DE Table” shows each contrast’s p-value. It
  looked for a `pvalue` column, while the results carry `pValue`, so the
  column was never there.
- Multi DE results added through Browse Projects or Manual Import are
  read like the dataset in the main tabs, with the sample table taken
  from the result itself. Before, it came from `../input_dataset.tsv`,
  and a result without that file got an empty sample table.
- An RNA-seq `SummarizedExperiment` whose p-value column is named
  `P Value` lists its features by p-value, as one named `pValue` or
  `P.Value` does.
- An `.h5ad` written by another tool, with gene symbols in a `gene_name`
  or `gene_symbol` column instead of `label_gene_symbol`, labels its
  features with those symbols; before, the labels were the feature ids,
  while the gene picker offered the symbols.
- DE QC offers the GC content and gene width plots for `.h5ad` files
  that carry feature lengths and GC content. Files written by
  [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  store them as `label_featwidth` and `label_gc`, which the DE table
  never picked up, so only the standard view was offered.
- The DESCRIPTION declares every package the app uses: pathview,
  enrichplot, KEGGREST, ggvenn, cowplot, colorspace, scales and tidyr
  were missing, so on a fresh install the app failed at startup
  (pathview) or when reaching the Venn diagram (ggvenn). The two
  `shinyalert` pop-ups are now Shiny’s own modal and notification, so
  shinyalert is no longer needed.
- Heatmap and Ordination treat peptide and phospho datasets like protein
  datasets when choosing the batch-correction defaults: with two
  factors, the design is locked to the group factor, as Boxplot and
  Correlation already did. Before, the design defaulted to “None”, so
  batch correction did not protect the group effect. All four tabs now
  use one rule.
- Batch correction is no longer applied to raw (linear) proteomics
  intensities as if they were log2: the tabs say that log-scale data is
  needed, as they do for raw RNA-seq counts. Only the transformed
  intensities count as log2.
- The proteomics Volcano download file names carry the contrast, and the
  colour-template file name the dataset’s design label. Both used to
  take a value set by the ORA tab.
- The version follows Bioconductor’s numbering, where an even middle
  number marks a release.
- anndataR is installed from Bioconductor: the DESCRIPTION no longer has
  a `Remotes:` field pointing to GitHub, and requires anndataR \>= 1.2.3
  (Bioconductor 3.23), the version the deployed app reads with. A file
  that is not valid HDF5 is read once before the error is shown, not
  three times.
- The demo dataset (`inst/examples/genomics/deResult.rds`) is 4.4 MB
  instead of 34 MB. It shows the same samples, features, DE table and
  ORA/GSEA results; removed are parts the app never reads: the DESeq2
  object, the gene sets of GO terms that are in no result table, and
  duplicate annotation columns. `inst/scripts/trim_demo_deResult.R`
  records the trimming.
- Every exported function and class has a help-page example that runs on
  the shipped example datasets, and says what it returns. The
  DESCRIPTION gives the package’s URL, bug tracker, Bioconductor
  categories, and the maintainer’s ORCID.
- exploreDE has a hex sticker (`man/figures/logo.png`, drawn by
  `inst/scripts/make_hex_logo.R`): a volcano plot whose point density
  draws a cat face. The README opens with what exploreDE is.
- The hex sticker now has dotted cat eyes, ears and whiskers, a
  background-coloured mouth, and subdued grey “re” lettering. The README
  links the sticker to the full-size PNG.
- The sticker colours every point, the cat’s eyes, ears and whiskers
  included, by the volcano’s cutoffs: blue or orange beyond a log2 fold
  change of ±1 and above the FDR 0.05 line, grey elsewhere. The eyes,
  ears and whiskers are jittered, so they read as scattered points like
  the rest of the volcano.
- Documentation website at <https://fgcz.github.io/exploreDE/>
  (pkgdown): the README, the help pages, the bring-your-own-data
  vignette and this changelog, rebuilt on every push to `main`.
- Help pages show code and links as such: backticks and `[function()]`
  references in the documentation were printed literally before.
- `ExploreDataset` objects print the same summary as before, now from a
  [`format()`](https://rdrr.io/r/base/format.html) method: `format(eds)`
  returns its lines, and the `eds$print()` method is gone (`print(eds)`
  works as before).

## exploreDE 3.9.0

### Breaking changes

- Saved sessions are `.json` only (format `exploreDE-session/2`). The
  `.rds` session format is removed: An uploaded `.rds` is refused with a
  message asking the user to save again as `.json`. Reading an `.rds`
  can run code; reading JSON can’t. Save / Load Session is out of beta.
- Only the documented API is exported:
  [`exploreDE()`](https://fgcz.github.io/exploreDE/reference/exploreDE.md),
  [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md),
  [`new_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/new_explore_dataset.md),
  [`validate_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/validate_explore_dataset.md),
  and the `ExploreDataset*` classes. The 89 app helpers
  (`generateParams()`, `step_*()`, `banana_*()`, `save_state_*()`, …)
  are internal, so `exploreDE::<helper>()` no longer works; the app gets
  them by sourcing `R/` or, when installed, from the package namespace.

### New features

- Open an `.h5ad` from a link: `?data=<percent-encoded https URL>`, from
  an allowed source (APB Export jobs and public pStore by default; set
  `EXPLOREDE_URL_SOURCES` to change them). Plain `https://` only, no
  redirects, up to 1 GB, then the same checks as an upload.
- JSON session files also save and restore the Multi DE datasets, in
  order, with their settings.
- Per-tab contrast selectors throughout the app, kept in sync.
- `.h5ad` DE tables can say which features entered the test: Optional
  `used_in_test` column (spec §5).
  [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  picks up a `usedInTest` / `used_in_test` column from the DE table
  through the column dictionary, so files built from SUSHI results show
  the same tested features as the SUSHI original. Without it, “tested”
  still means a non-`NA` p-value.
- Reads the prolfquapp DEA’s `AnnData.h5ad` and prolfquapp ≥ 2.10
  `SummarizedExperiment` results.
- Container RAM/CPU badge in the header.
- GO behaves like KEGG: The GO parts of the app (tile plot selector,
  heatmap GO sub-tab) appear only when the features carry GO annotation
  columns, and any subset of `GO BP` / `GO MF` / `GO CC` works. The
  Pathway Tile Plots tab shows when GO annotation or a supported KEGG
  species is available, so human/mouse/rat data without GO columns now
  get KEGG tiles. New exported helpers `go_annotation_columns()` and
  `features_in_go_term()`.

### Bug fixes

- An `.h5ad` without `uns["de"]` loads as counts-only; before, `uns$de`
  partial-matched `uns$design` and the file failed to load. Pathway
  tables are read with `[[` for the same reason.
- Multi DE: Dataset paths outside the data roots (or containing `..`)
  are refused, and each dataset gets its own LDAP project check.
- An empty feature bucket no longer prints “display:none;”.
- GO tile plots and the GO heatmap work on `.h5ad` files: GO membership
  is read from the features, not the DE table (which `.h5ad` files don’t
  give GO columns). Before, they failed with “object ‘GO BP’ not found”.
- GO tile plots skip features without a log2 fold change and draw no
  stars for an `NA` p-value, instead of failing (“NA/NaN/Inf in foreign
  function call”, “missing value where TRUE/FALSE needed”).
- GO tabs follow the changed `twoGroups` output of ezRun; small ORA/GSEA
  and volcano fixes.
- `rhdf5` is a declared dependency, so `.h5ad` files load on a fresh
  install; before, anndataR failed with “HDF5AnnData requires the rhdf5
  package”.
- The app gets the package’s functions, internal helpers included, when
  it runs from the installed package, however it is started
  ([`exploreDE::exploreDE()`](https://fgcz.github.io/exploreDE/reference/exploreDE.md),
  or [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) on
  `system.file("app", package = "exploreDE")`). Before, it needed
  [`library(exploreDE)`](https://fgcz.github.io/exploreDE/) first
  (‘could not find function “legend_offset_ui”’), and internal helpers
  such as the prolfquapp `.h5ad` reader were never available.
- The ORA and GSEA tabs keep their settings apart. Before, for an
  uploaded `.h5ad`, the GSEA tab overwrote the ORA FDR threshold, and
  the ORA table failed with “length(cuts) must be equal to
  length(values) - 1”.
- The Volcano and DE Table tabs work for an `.h5ad` whose features have
  no `description` column (the description falls back to the gene name,
  as for SUSHI results). Before, the volcano failed with “undefined
  columns selected”.
- Species without KEGG support (anything but human, mouse and rat):
  Changing the contrast no longer ends the session with “EXPR must be a
  length 1 vector”. An `.h5ad` with two or more contrasts hung on load.
- Level labels with Greek letters, or with a space between two tags,
  keep their formatting when a PDF device is open on Linux. Before, the
  label check measured them on that device, failed, and fell back to
  plain text.

### Documentation

- Written format specs in `inst/spec/`: `exploreDE-h5ad-1.0.md` (the
  `.h5ad` layout, proteomics profile including precursor level,
  prolfquapp dialects, and a minimal Python writer) and
  `exploreDE-session-2.md` (session files).
- README, the bring-your-own-data vignette, and the AI skill now cover
  writing files outside R, opening files from a link, and saving
  sessions.
- The help pages of the `ExploreDataset*` classes describe every public
  method and its arguments.

### Removed

- The App Data tab (deprecated).
- `devtools` as a dependency: installing exploreDE no longer installs
  it. `dendextend` and `rlang`, which the package code calls, move from
  Suggests to Imports.

## exploreDE 3.8.0

- Project Jambon: One column dictionary with ordered patterns;
  enrichment rebuilt from tidy tables on load;
  [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  documented in the bring-your-own-data vignette and an AI skill;
  `exploreDE(file)` opens a file directly, and the demo ships with the
  package; ezRun dropped.
- Level labels, heatmap legend sliders, feature bucket fixes, and “DEG”
  → “DE” throughout the interface.
