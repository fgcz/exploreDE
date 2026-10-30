# exploreDE 3.9.0

## Breaking changes

- Saved sessions are `.json` only (format `exploreDE-session/2`). The `.rds`
  session format is removed: An uploaded `.rds` is refused with a message
  asking the user to save again as `.json`. Reading an `.rds` can run code;
  reading JSON can't. Save / Load Session is out of beta.

## New features

- Open an `.h5ad` from a link: `?data=<percent-encoded https URL>`, from an
  allowed source (APB Export jobs and public pStore by default; set
  `EXPLOREDE_URL_SOURCES` to change them). Plain `https://` only, no
  redirects, up to 1 GB, then the same checks as an upload.
- JSON session files also save and restore the Multi DE datasets, in order,
  with their settings.
- Per-tab contrast selectors throughout the app, kept in sync.
- `.h5ad` DE tables can say which features entered the test: Optional
  `used_in_test` column (spec §5). `build_explore_h5ad()` picks up a
  `usedInTest` / `used_in_test` column from the DE table through the column
  dictionary, so files built from SUSHI results show the same tested features
  as the SUSHI original. Without it, "tested" still means a non-`NA` p-value.
- Reads the prolfquapp DEA's `AnnData.h5ad` and prolfquapp ≥ 2.10
  `SummarizedExperiment` results.
- Container RAM/CPU badge in the header.
- GO behaves like KEGG: The GO parts of the app (tile plot selector, heatmap
  GO sub-tab) appear only when the features carry GO annotation columns, and
  any subset of `GO BP` / `GO MF` / `GO CC` works. The Pathway Tile Plots tab
  shows when GO annotation or a supported KEGG species is available, so
  human/mouse/rat data without GO columns now get KEGG tiles. New exported
  helpers `go_annotation_columns()` and `features_in_go_term()`.

## Bug fixes

- An `.h5ad` without `uns["de"]` loads as counts-only; before, `uns$de`
  partial-matched `uns$design` and the file failed to load. Pathway tables
  are read with `[[` for the same reason.
- Multi DE: Dataset paths outside the data roots (or containing `..`) are
  refused, and each dataset gets its own LDAP project check.
- An empty feature bucket no longer prints "display:none;".
- GO tile plots and the GO heatmap work on `.h5ad` files: GO membership is
  read from the features, not the DE table (which `.h5ad` files don't give
  GO columns). Before, they failed with "object 'GO BP' not found".
- GO tile plots skip features without a log2 fold change and draw no stars
  for an `NA` p-value, instead of failing ("NA/NaN/Inf in foreign function
  call", "missing value where TRUE/FALSE needed").
- GO tabs follow the changed `twoGroups` output of ezRun; small ORA/GSEA and
  volcano fixes.

## Documentation

- Written format specs in `inst/spec/`: `exploreDE-h5ad-1.0.md` (the `.h5ad`
  layout, proteomics profile including precursor level, prolfquapp dialects,
  and a minimal Python writer) and `exploreDE-session-2.md` (session files).
- README, the bring-your-own-data vignette, and the AI skill now cover
  writing files outside R, opening files from a link, and saving sessions.

## Removed

- The App Data tab (deprecated).

# exploreDE 3.8.0

- Project Jambon: One column dictionary with ordered patterns; enrichment
  rebuilt from tidy tables on load; `build_explore_h5ad()` documented in the
  bring-your-own-data vignette and an AI skill; `exploreDE(file)` opens a file
  directly, and the demo ships with the package; ezRun dropped.
- Level labels, heatmap legend sliders, feature bucket fixes, and "DEG" →
  "DE" throughout the interface.
