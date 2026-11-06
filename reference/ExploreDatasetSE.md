# Abstract intermediate class for SE-backed datasets

`ExploreDatasetSE` holds the `SummarizedExperiment` backend reference
and implements the accessor contract uniformly for every SE shape.
Concrete subclasses (`ExploreDatasetSE_RNASeq`, `_CountQC`,
`_CrisprCountQC`, `_Proteomics`) override two private methods only:

- `.build_config(se)` — returns the per-instance config list at
  `$initialize()` time

- `.prep()` — populates `private$.features_df`, `private$.de_dfs`,
  `private$.fid_lookup` from the SE

This class is not instantiable on its own — instantiate a concrete
subclass instead.

## Value

An R6 class generator; the concrete subclasses' `$new()` returns an
object inheriting from `ExploreDatasetSE`.

## Super class

[`ExploreDataset`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md)
-\> `ExploreDatasetSE`

## Methods

### Public methods

- [`ExploreDatasetSE$new()`](#method-ExploreDatasetSE-initialize)

- [`ExploreDatasetSE$backend()`](#method-ExploreDatasetSE-backend)

- [`ExploreDatasetSE$backend_kind()`](#method-ExploreDatasetSE-backend_kind)

- [`ExploreDatasetSE$samples()`](#method-ExploreDatasetSE-samples)

- [`ExploreDatasetSE$features()`](#method-ExploreDatasetSE-features)

- [`ExploreDatasetSE$de()`](#method-ExploreDatasetSE-de)

- [`ExploreDatasetSE$assay_names()`](#method-ExploreDatasetSE-assay_names)

- [`ExploreDatasetSE$assay()`](#method-ExploreDatasetSE-assay)

- [`ExploreDatasetSE$enrichment()`](#method-ExploreDatasetSE-enrichment)

- [`ExploreDatasetSE$enrichment_input()`](#method-ExploreDatasetSE-enrichment_input)

- [`ExploreDatasetSE$clone()`](#method-ExploreDatasetSE-clone)

Inherited methods

- [`ExploreDataset$all_pathways()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-all_pathways)
- [`ExploreDataset$assay_label()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-assay_label)
- [`ExploreDataset$background_expression()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-background_expression)
- [`ExploreDataset$contrast_label()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-contrast_label)
- [`ExploreDataset$contrast_meta()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-contrast_meta)
- [`ExploreDataset$contrasts()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-contrasts)
- [`ExploreDataset$data_type()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-data_type)
- [`ExploreDataset$de_all()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-de_all)
- [`ExploreDataset$de_display_extras()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-de_display_extras)
- [`ExploreDataset$enrichment_filter()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-enrichment_filter)
- [`ExploreDataset$enrichment_kinds()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-enrichment_kinds)
- [`ExploreDataset$enrichment_table()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-enrichment_table)
- [`ExploreDataset$factor_levels()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-factor_levels)
- [`ExploreDataset$factor_names()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-factor_names)
- [`ExploreDataset$feature_display()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-feature_display)
- [`ExploreDataset$feature_id_kind()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-feature_id_kind)
- [`ExploreDataset$feature_ids()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-feature_ids)
- [`ExploreDataset$format()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-format)
- [`ExploreDataset$has_assay()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-has_assay)
- [`ExploreDataset$has_de()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-has_de)
- [`ExploreDataset$has_enrichment()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-has_enrichment)
- [`ExploreDataset$has_pvalue()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-has_pvalue)
- [`ExploreDataset$n_features()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-n_features)
- [`ExploreDataset$n_samples()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-n_samples)
- [`ExploreDataset$omics_class()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-omics_class)
- [`ExploreDataset$organism()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-organism)
- [`ExploreDataset$primary_assay()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-primary_assay)
- [`ExploreDataset$primary_assay_name()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-primary_assay_name)
- [`ExploreDataset$primary_contrast()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-primary_contrast)
- [`ExploreDataset$primary_grouping()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-primary_grouping)
- [`ExploreDataset$proteomics_choices()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-proteomics_choices)
- [`ExploreDataset$raw_assay()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-raw_assay)
- [`ExploreDataset$raw_assay_name()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-raw_assay_name)
- [`ExploreDataset$ref_build()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-ref_build)
- [`ExploreDataset$sample_factor()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-sample_factor)
- [`ExploreDataset$sample_ids()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-sample_ids)
- [`ExploreDataset$secondary_grouping()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-secondary_grouping)
- [`ExploreDataset$sgrna_to_gene()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-sgrna_to_gene)
- [`ExploreDataset$source()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-source)

------------------------------------------------------------------------

### `ExploreDatasetSE$new()`

Wrap a `SummarizedExperiment`: store it, build the config from it unless
one is given, and prepare the feature and DE tables.

#### Usage

    ExploreDatasetSE$new(se, config = NULL)

#### Arguments

- `se`:

  A `SummarizedExperiment`.

- `config`:

  Optional pre-built config list; built from `se` by the subclass when
  `NULL`.

------------------------------------------------------------------------

### `ExploreDatasetSE$backend()`

The wrapped `SummarizedExperiment`.

#### Usage

    ExploreDatasetSE$backend()

------------------------------------------------------------------------

### `ExploreDatasetSE$backend_kind()`

Returns `"se"`.

#### Usage

    ExploreDatasetSE$backend_kind()

------------------------------------------------------------------------

### `ExploreDatasetSE$samples()`

Sample metadata from
[`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
with the SUSHI ` [Factor]` suffix stripped from column names, rownames
set to sample IDs, and a `Name` column added.

#### Usage

    ExploreDatasetSE$samples()

------------------------------------------------------------------------

### `ExploreDatasetSE$features()`

Feature annotation table, prepared from
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
by the subclass.

#### Usage

    ExploreDatasetSE$features()

------------------------------------------------------------------------

### `ExploreDatasetSE$de()`

DE table of one contrast, prepared from
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
by the subclass; an empty data.frame when none was prepared. Errors on
an unknown contrast.

#### Usage

    ExploreDatasetSE$de(contrast)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

------------------------------------------------------------------------

### `ExploreDatasetSE$assay_names()`

Canonical assay names from the config's assay map whose SE assay exists;
computed assays are always included.

#### Usage

    ExploreDatasetSE$assay_names()

------------------------------------------------------------------------

### `ExploreDatasetSE$assay()`

One assay as a features x samples matrix with rows renamed to feature
IDs; `normalised_log2`, `vst`, `fpkm` and `tpm` are computed when the
config maps them to no SE assay.

#### Usage

    ExploreDatasetSE$assay(name, features = NULL, samples = NULL)

#### Arguments

- `name`:

  Canonical assay name, one of `assay_names()`.

- `features`:

  Feature IDs (rows) to keep; `NULL` keeps all.

- `samples`:

  Sample IDs (columns) to keep; `NULL` keeps all.

------------------------------------------------------------------------

### `ExploreDatasetSE$enrichment()`

Enrichment results of one kind, read from `metadata(se)` and flattened
to a named list; an empty list when absent. `contrast` is not used.

#### Usage

    ExploreDatasetSE$enrichment(contrast, kind)

#### Arguments

- `contrast`:

  Contrast name; not used by this backend.

- `kind`:

  `"ora"` or `"gsea"`.

------------------------------------------------------------------------

### `ExploreDatasetSE$enrichment_input()`

Enrichment inputs read from `metadata(se)$enrichInput`, with `NULL`
fields and no selections when absent. `contrast` is not used.

#### Usage

    ExploreDatasetSE$enrichment_input(contrast)

#### Arguments

- `contrast`:

  Contrast name; not used by this backend.

------------------------------------------------------------------------

### `ExploreDatasetSE$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ExploreDatasetSE$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
eds <- new_explore_dataset(se)
inherits(eds, "ExploreDatasetSE")
#> [1] TRUE
eds$assay_names()
#> [1] "intensities"
dim(eds$raw_assay())
#> [1] 1384    6
```
