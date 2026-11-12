# Concrete `ExploreDataset` backend over a parquet directory

The exploreDElight backend. Reads the parquet directory a native Rust
process writes from a canonical exploreDE/1.0 `.h5ad`, and satisfies the
same accessor contract as every other backend. Inherits all config
derivation and DE assembly from
[ExploreDatasetH5AD](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.md);
only the loader differs.

## Super classes

[`ExploreDataset`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md)
-\>
[`ExploreDatasetH5AD`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.md)
-\> `ExploreDatasetParquet`

## Methods

### Public methods

- [`ExploreDatasetParquet$new()`](#method-ExploreDatasetParquet-initialize)

- [`ExploreDatasetParquet$backend_kind()`](#method-ExploreDatasetParquet-backend_kind)

- [`ExploreDatasetParquet$manifest()`](#method-ExploreDatasetParquet-manifest)

- [`ExploreDatasetParquet$source_dir()`](#method-ExploreDatasetParquet-source_dir)

- [`ExploreDatasetParquet$verify_source()`](#method-ExploreDatasetParquet-verify_source)

- [`ExploreDatasetParquet$all_pathways()`](#method-ExploreDatasetParquet-all_pathways)

- [`ExploreDatasetParquet$clone()`](#method-ExploreDatasetParquet-clone)

Inherited methods

- [`ExploreDataset$assay_is_log2()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-assay_is_log2)
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
- [`ExploreDatasetH5AD$assay()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.html#method-assay)
- [`ExploreDatasetH5AD$assay_names()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.html#method-assay_names)
- [`ExploreDatasetH5AD$backend()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.html#method-backend)
- [`ExploreDatasetH5AD$de()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.html#method-de)
- [`ExploreDatasetH5AD$enrichment()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.html#method-enrichment)
- [`ExploreDatasetH5AD$enrichment_input()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.html#method-enrichment_input)
- [`ExploreDatasetH5AD$enrichment_table()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.html#method-enrichment_table)
- [`ExploreDatasetH5AD$features()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.html#method-features)
- [`ExploreDatasetH5AD$samples()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.html#method-samples)

------------------------------------------------------------------------

### `ExploreDatasetParquet$new()`

Read the parquet directory (manifest, obs, var, layers, varm, uns),
build the config, prepare the feature and DE tables, and warn when the
source `.h5ad` has changed since export.

#### Usage

    ExploreDatasetParquet$new(input)

#### Arguments

- `input`:

  Path to the parquet directory.

------------------------------------------------------------------------

### `ExploreDatasetParquet$backend_kind()`

Returns `"parquet"`.

#### Usage

    ExploreDatasetParquet$backend_kind()

------------------------------------------------------------------------

### `ExploreDatasetParquet$manifest()`

The parsed `manifest.json` of the parquet directory.

#### Usage

    ExploreDatasetParquet$manifest()

------------------------------------------------------------------------

### `ExploreDatasetParquet$source_dir()`

Path of the parquet directory.

#### Usage

    ExploreDatasetParquet$source_dir()

------------------------------------------------------------------------

### `ExploreDatasetParquet$verify_source()`

Compare the md5 of the source `.h5ad` with the one recorded in the
manifest: `TRUE` when it matches, `FALSE` when not, `NA` when there is
nothing to check against.

#### Usage

    ExploreDatasetParquet$verify_source(path = NULL)

#### Arguments

- `path`:

  Path to the source `.h5ad`; `NULL` (default) uses the path recorded in
  the manifest.

------------------------------------------------------------------------

### `ExploreDatasetParquet$all_pathways()`

GO terms as in the base class, but an empty vector when clusterProfiler
is not installed.

#### Usage

    ExploreDatasetParquet$all_pathways()

------------------------------------------------------------------------

### `ExploreDatasetParquet$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ExploreDatasetParquet$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
