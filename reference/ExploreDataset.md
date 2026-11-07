# Abstract base class for exploreDE's storage-agnostic dataset wrapper

`ExploreDataset` defines the accessor contract that every concrete
subclass must implement. It is not instantiable on its own — calling
`ExploreDataset$new()` raises an error. Use
[`new_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/new_explore_dataset.md)
to build a concrete instance, or instantiate a subclass directly.

Methods implemented on the base class read only from `private$.config`
and are inherited unchanged by every subclass. Methods that touch the
storage backend are virtual (`stop("not implemented")`) and must be
overridden.

## Value

An R6 class generator. `$new()` always errors; the concrete subclasses
return an object inheriting from `ExploreDataset`.

## Methods

### Public methods

- [`ExploreDataset$new()`](#method-ExploreDataset-initialize)

- [`ExploreDataset$omics_class()`](#method-ExploreDataset-omics_class)

- [`ExploreDataset$data_type()`](#method-ExploreDataset-data_type)

- [`ExploreDataset$has_de()`](#method-ExploreDataset-has_de)

- [`ExploreDataset$has_enrichment()`](#method-ExploreDataset-has_enrichment)

- [`ExploreDataset$organism()`](#method-ExploreDataset-organism)

- [`ExploreDataset$ref_build()`](#method-ExploreDataset-ref_build)

- [`ExploreDataset$feature_id_kind()`](#method-ExploreDataset-feature_id_kind)

- [`ExploreDataset$background_expression()`](#method-ExploreDataset-background_expression)

- [`ExploreDataset$backend_kind()`](#method-ExploreDataset-backend_kind)

- [`ExploreDataset$backend()`](#method-ExploreDataset-backend)

- [`ExploreDataset$samples()`](#method-ExploreDataset-samples)

- [`ExploreDataset$sample_ids()`](#method-ExploreDataset-sample_ids)

- [`ExploreDataset$n_samples()`](#method-ExploreDataset-n_samples)

- [`ExploreDataset$factor_names()`](#method-ExploreDataset-factor_names)

- [`ExploreDataset$factor_levels()`](#method-ExploreDataset-factor_levels)

- [`ExploreDataset$primary_grouping()`](#method-ExploreDataset-primary_grouping)

- [`ExploreDataset$secondary_grouping()`](#method-ExploreDataset-secondary_grouping)

- [`ExploreDataset$sample_factor()`](#method-ExploreDataset-sample_factor)

- [`ExploreDataset$features()`](#method-ExploreDataset-features)

- [`ExploreDataset$feature_ids()`](#method-ExploreDataset-feature_ids)

- [`ExploreDataset$n_features()`](#method-ExploreDataset-n_features)

- [`ExploreDataset$feature_display()`](#method-ExploreDataset-feature_display)

- [`ExploreDataset$assay_names()`](#method-ExploreDataset-assay_names)

- [`ExploreDataset$has_assay()`](#method-ExploreDataset-has_assay)

- [`ExploreDataset$assay()`](#method-ExploreDataset-assay)

- [`ExploreDataset$assay_label()`](#method-ExploreDataset-assay_label)

- [`ExploreDataset$raw_assay_name()`](#method-ExploreDataset-raw_assay_name)

- [`ExploreDataset$primary_assay_name()`](#method-ExploreDataset-primary_assay_name)

- [`ExploreDataset$raw_assay()`](#method-ExploreDataset-raw_assay)

- [`ExploreDataset$primary_assay()`](#method-ExploreDataset-primary_assay)

- [`ExploreDataset$contrasts()`](#method-ExploreDataset-contrasts)

- [`ExploreDataset$primary_contrast()`](#method-ExploreDataset-primary_contrast)

- [`ExploreDataset$de()`](#method-ExploreDataset-de)

- [`ExploreDataset$de_all()`](#method-ExploreDataset-de_all)

- [`ExploreDataset$contrast_meta()`](#method-ExploreDataset-contrast_meta)

- [`ExploreDataset$contrast_label()`](#method-ExploreDataset-contrast_label)

- [`ExploreDataset$de_display_extras()`](#method-ExploreDataset-de_display_extras)

- [`ExploreDataset$has_pvalue()`](#method-ExploreDataset-has_pvalue)

- [`ExploreDataset$enrichment_kinds()`](#method-ExploreDataset-enrichment_kinds)

- [`ExploreDataset$enrichment()`](#method-ExploreDataset-enrichment)

- [`ExploreDataset$enrichment_input()`](#method-ExploreDataset-enrichment_input)

- [`ExploreDataset$enrichment_table()`](#method-ExploreDataset-enrichment_table)

- [`ExploreDataset$enrichment_filter()`](#method-ExploreDataset-enrichment_filter)

- [`ExploreDataset$sgrna_to_gene()`](#method-ExploreDataset-sgrna_to_gene)

- [`ExploreDataset$proteomics_choices()`](#method-ExploreDataset-proteomics_choices)

- [`ExploreDataset$all_pathways()`](#method-ExploreDataset-all_pathways)

- [`ExploreDataset$source()`](#method-ExploreDataset-source)

- [`ExploreDataset$format()`](#method-ExploreDataset-format)

- [`ExploreDataset$clone()`](#method-ExploreDataset-clone)

------------------------------------------------------------------------

### `ExploreDataset$new()`

Always errors: the base class is abstract. Use
[`new_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/new_explore_dataset.md)
or a concrete subclass.

#### Usage

    ExploreDataset$new(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `ExploreDataset$omics_class()`

Omics class from the config: `"genomics"` or `"proteomics"`.

#### Usage

    ExploreDataset$omics_class()

------------------------------------------------------------------------

### `ExploreDataset$data_type()`

Dataset type from the config, e.g. `"RNASeq"`, `"CountQC"`,
`"CrisprCountQC"` or `"proteomics"`.

#### Usage

    ExploreDataset$data_type()

------------------------------------------------------------------------

### `ExploreDataset$has_de()`

`TRUE` when the dataset has at least one contrast.

#### Usage

    ExploreDataset$has_de()

------------------------------------------------------------------------

### `ExploreDataset$has_enrichment()`

`TRUE` when the config lists at least one enrichment kind.

#### Usage

    ExploreDataset$has_enrichment()

------------------------------------------------------------------------

### `ExploreDataset$organism()`

Organism from the config, or `NA` when unknown.

#### Usage

    ExploreDataset$organism()

------------------------------------------------------------------------

### `ExploreDataset$ref_build()`

Reference genome build from the config, or `NA` when unknown.

#### Usage

    ExploreDataset$ref_build()

------------------------------------------------------------------------

### `ExploreDataset$feature_id_kind()`

What a feature ID denotes, from the config, e.g. `"gene"`, `"isoform"`,
`"protein"`, `"peptide"`, `"phospho"` or `"guide"`.

#### Usage

    ExploreDataset$feature_id_kind()

------------------------------------------------------------------------

### `ExploreDataset$background_expression()`

Low-signal offset added to expression values before taking log2; from
the config, default 10.

#### Usage

    ExploreDataset$background_expression()

------------------------------------------------------------------------

### `ExploreDataset$backend_kind()`

Short tag naming the storage backend (`"se"`, `"h5ad"` or `"parquet"`).
Implemented by subclasses.

#### Usage

    ExploreDataset$backend_kind()

------------------------------------------------------------------------

### `ExploreDataset$backend()`

The underlying storage object. Implemented by subclasses.

#### Usage

    ExploreDataset$backend()

------------------------------------------------------------------------

### `ExploreDataset$samples()`

Sample metadata as a data.frame, one row per sample, rownames are sample
IDs. Implemented by subclasses.

#### Usage

    ExploreDataset$samples()

------------------------------------------------------------------------

### `ExploreDataset$sample_ids()`

Sample IDs: the rownames of `samples()`.

#### Usage

    ExploreDataset$sample_ids()

------------------------------------------------------------------------

### `ExploreDataset$n_samples()`

Number of samples: the number of rows of `samples()`.

#### Usage

    ExploreDataset$n_samples()

------------------------------------------------------------------------

### `ExploreDataset$factor_names()`

Names of the experimental factors from the config; empty when there are
none.

#### Usage

    ExploreDataset$factor_names()

------------------------------------------------------------------------

### `ExploreDataset$factor_levels()`

Levels of one experimental factor, from the config.

#### Usage

    ExploreDataset$factor_levels(factor)

#### Arguments

- `factor`:

  Factor name, one of `factor_names()`.

------------------------------------------------------------------------

### `ExploreDataset$primary_grouping()`

Name of the factor used as the primary sample grouping, from the config.

#### Usage

    ExploreDataset$primary_grouping()

------------------------------------------------------------------------

### `ExploreDataset$secondary_grouping()`

Name of the secondary grouping factor, or `NA` when there is none.

#### Usage

    ExploreDataset$secondary_grouping()

------------------------------------------------------------------------

### `ExploreDataset$sample_factor()`

One column of `samples()`: the values of a factor, one per sample.

#### Usage

    ExploreDataset$sample_factor(factor)

#### Arguments

- `factor`:

  Column name in `samples()`, usually one of `factor_names()`.

------------------------------------------------------------------------

### `ExploreDataset$features()`

Feature annotation as a data.frame, one row per feature, rownames are
feature IDs, with a `display_name` column. Implemented by subclasses.

#### Usage

    ExploreDataset$features()

------------------------------------------------------------------------

### `ExploreDataset$feature_ids()`

Feature IDs: the rownames of `features()`.

#### Usage

    ExploreDataset$feature_ids()

------------------------------------------------------------------------

### `ExploreDataset$n_features()`

Number of features: the number of rows of `features()`.

#### Usage

    ExploreDataset$n_features()

------------------------------------------------------------------------

### `ExploreDataset$feature_display()`

Display names of features, named by feature ID.

#### Usage

    ExploreDataset$feature_display(ids = NULL)

#### Arguments

- `ids`:

  Feature IDs to look up; `NULL` (default) returns all features.

------------------------------------------------------------------------

### `ExploreDataset$assay_names()`

Canonical names of the assays this dataset provides. Implemented by
subclasses.

#### Usage

    ExploreDataset$assay_names()

------------------------------------------------------------------------

### `ExploreDataset$has_assay()`

`TRUE` when `name` is one of `assay_names()`.

#### Usage

    ExploreDataset$has_assay(name)

#### Arguments

- `name`:

  Canonical assay name.

------------------------------------------------------------------------

### `ExploreDataset$assay()`

One assay as a features x samples matrix, optionally subset. Implemented
by subclasses.

#### Usage

    ExploreDataset$assay(name, features = NULL, samples = NULL)

#### Arguments

- `name`:

  Canonical assay name, one of `assay_names()`.

- `features`:

  Feature IDs (rows) to keep; `NULL` keeps all.

- `samples`:

  Sample IDs (columns) to keep; `NULL` keeps all.

------------------------------------------------------------------------

### `ExploreDataset$assay_label()`

Human-readable label for an assay, from the config; `name` itself when
no label is registered.

#### Usage

    ExploreDataset$assay_label(name)

#### Arguments

- `name`:

  Canonical assay name.

------------------------------------------------------------------------

### `ExploreDataset$raw_assay_name()`

Canonical name of the raw assay, from the config.

#### Usage

    ExploreDataset$raw_assay_name()

------------------------------------------------------------------------

### `ExploreDataset$primary_assay_name()`

Canonical name of the primary assay, from the config.

#### Usage

    ExploreDataset$primary_assay_name()

------------------------------------------------------------------------

### `ExploreDataset$raw_assay()`

The raw assay: shortcut for `assay(raw_assay_name(), ...)`.

#### Usage

    ExploreDataset$raw_assay(...)

#### Arguments

- `...`:

  Passed to `assay()`: `features`, `samples`.

------------------------------------------------------------------------

### `ExploreDataset$primary_assay()`

The primary assay: shortcut for `assay(primary_assay_name(), ...)`.

#### Usage

    ExploreDataset$primary_assay(...)

#### Arguments

- `...`:

  Passed to `assay()`: `features`, `samples`.

------------------------------------------------------------------------

### `ExploreDataset$contrasts()`

Names of the DE contrasts, from the config; empty when the dataset has
no DE results.

#### Usage

    ExploreDataset$contrasts()

------------------------------------------------------------------------

### `ExploreDataset$primary_contrast()`

Name of the primary contrast, from the config; `NA` when there are no
contrasts.

#### Usage

    ExploreDataset$primary_contrast()

------------------------------------------------------------------------

### `ExploreDataset$de()`

DE result table for one contrast, as a data.frame. Implemented by
subclasses.

#### Usage

    ExploreDataset$de(contrast)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

------------------------------------------------------------------------

### `ExploreDataset$de_all()`

DE result tables of all contrasts, as a list named by contrast.

#### Usage

    ExploreDataset$de_all()

------------------------------------------------------------------------

### `ExploreDataset$contrast_meta()`

Metadata list of one contrast from the config: test and reference group,
factor, method, thresholds. Errors on an unknown contrast.

#### Usage

    ExploreDataset$contrast_meta(contrast)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

------------------------------------------------------------------------

### `ExploreDataset$contrast_label()`

Display label of a contrast: currently the contrast name itself. Errors
on an unknown contrast.

#### Usage

    ExploreDataset$contrast_label(contrast)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

------------------------------------------------------------------------

### `ExploreDataset$de_display_extras()`

Names of extra DE columns shown for display only (e.g. `SaintScore`),
from the config; empty by default.

#### Usage

    ExploreDataset$de_display_extras()

------------------------------------------------------------------------

### `ExploreDataset$has_pvalue()`

`TRUE` when the DE table of at least one contrast has a non-missing raw
`pValue`; `FALSE` for FDR-only producers such as SAINTexpress.

#### Usage

    ExploreDataset$has_pvalue()

------------------------------------------------------------------------

### `ExploreDataset$enrichment_kinds()`

Enrichment kinds available, from the config: any of `"ora"` and
`"gsea"`.

#### Usage

    ExploreDataset$enrichment_kinds()

------------------------------------------------------------------------

### `ExploreDataset$enrichment()`

Enrichment results of one contrast and kind, as a flat named list of
clusterProfiler results. Implemented by subclasses.

#### Usage

    ExploreDataset$enrichment(contrast, kind)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

- `kind`:

  `"ora"` or `"gsea"`.

------------------------------------------------------------------------

### `ExploreDataset$enrichment_input()`

Inputs for enrichment of one contrast:
`list(seq_anno, ranked_l2fc, selections)`. Implemented by subclasses.

#### Usage

    ExploreDataset$enrichment_input(contrast)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

------------------------------------------------------------------------

### `ExploreDataset$enrichment_table()`

Enrichment results of one contrast and kind as one tidy data.frame,
derived from `enrichment()`; an empty data.frame when there are none.

#### Usage

    ExploreDataset$enrichment_table(contrast, kind)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

- `kind`:

  `"ora"` or `"gsea"`.

------------------------------------------------------------------------

### `ExploreDataset$enrichment_filter()`

The `p.adjust` cut-off applied to the stored enrichment tables, from the
config; `NA` when none was applied.

#### Usage

    ExploreDataset$enrichment_filter()

------------------------------------------------------------------------

### `ExploreDataset$sgrna_to_gene()`

Named character vector mapping sgRNA IDs to gene names (CRISPR
datasets); `NULL` otherwise.

#### Usage

    ExploreDataset$sgrna_to_gene()

------------------------------------------------------------------------

### `ExploreDataset$proteomics_choices()`

Choices and lookup for the proteomics feature picker,
`list(choices, lookup)`; `NULL` when the dataset has none.

#### Usage

    ExploreDataset$proteomics_choices(kind = c("gene", "protein"))

#### Arguments

- `kind`:

  `"gene"` (default) or `"protein"`.

------------------------------------------------------------------------

### `ExploreDataset$all_pathways()`

GO terms named in the `GO BP`, `GO MF` and `GO CC` columns of
`features()`, as `"GO:ID term name"` strings; empty when those columns
are absent.

#### Usage

    ExploreDataset$all_pathways()

------------------------------------------------------------------------

### `ExploreDataset$source()`

Provenance list from the config: backend kind, schema version and build
time.

#### Usage

    ExploreDataset$source()

------------------------------------------------------------------------

### `ExploreDataset$format()`

Summary of the dataset (type, dimensions, factors, contrasts, assays,
enrichment), one line per element.
[`print()`](https://rdrr.io/r/base/print.html) on the object shows these
lines.

#### Usage

    ExploreDataset$format(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `ExploreDataset$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ExploreDataset$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
#> Loading required namespace: SummarizedExperiment
eds <- new_explore_dataset(se)
inherits(eds, "ExploreDataset")
#> [1] TRUE
eds$contrasts()
#> [1] "B_vs_A"
head(eds$de(eds$primary_contrast()))
#>        nrPeptides IDcolumn              fasta.id           description
#> P0AFG3         25   P0AFG3  sp|P0AFG3|ODO1_ECOLI  sp_P0AFG3_ODO1_ECOLI
#> P0AC41         15   P0AC41  sp|P0AC41|SDHA_ECOLI  sp_P0AC41_SDHA_ECOLI
#> P0A7W1          5   P0A7W1   sp|P0A7W1|RS5_ECOLI   sp_P0A7W1_RS5_ECOLI
#> P0A7V0          8   P0A7V0   sp|P0A7V0|RS2_ECOLI   sp_P0A7V0_RS2_ECOLI
#> A1L429          1   A1L429 sp|A1L429|GG12C_HUMAN sp_A1L429_GG12C_HUMAN
#> P0AFG6         10   P0AFG6  sp|P0AFG6|ODO2_ECOLI  sp_P0AFG6_ODO2_ECOLI
#>        protein_length nr_tryptic_peptides   CON   REV              modelName
#> P0AFG3            933                  44 FALSE FALSE Linear_Model_moderated
#> P0AC41            588                  31 FALSE FALSE Linear_Model_moderated
#> P0A7W1            167                  11 FALSE FALSE Linear_Model_moderated
#> P0A7V0            241                  14 FALSE FALSE Linear_Model_moderated
#> A1L429            117                   2 FALSE FALSE Imputed_Mean_moderated
#> P0AFG6            405                  20 FALSE FALSE Linear_Model_moderated
#>        contrast log2Ratio  std.error   avgAbd statistic       df       pValue
#> P0AFG3   B_vs_A  2.076169 0.02592986 19.27265  26.95160 5.832013 2.407749e-07
#> P0AC41   B_vs_A  2.207672 0.05267450 19.16701  25.70541 5.832013 3.166968e-07
#> P0A7W1   B_vs_A  2.232960 0.05541211 21.40682  25.64938 5.832013 3.207217e-07
#> P0A7V0   B_vs_A  2.290956 0.06851577 20.84410  24.57190 5.832013 4.111082e-07
#> A1L429   B_vs_A -7.134870 0.06665272 21.11332 -68.04320 3.798723 5.272144e-07
#> P0AFG6   B_vs_A  2.690083 0.11613502 20.00267  22.16970 5.832013 7.449241e-07
#>         conf.low conf.high      sigma          fdr feature_id effect_log2fc
#> P0AFG3  1.843691  2.308647 0.09434608 0.0001329935     P0AFG3      2.076169
#> P0AC41  1.948485  2.466860 0.10518546 0.0001329935     P0AC41      2.207672
#> P0A7W1  1.970231  2.495689 0.10662272 0.0001329935     P0A7W1      2.232960
#> P0A7V0  2.009583  2.572328 0.11418883 0.0001329935     P0A7V0      2.290956
#> A1L429 -7.499009 -6.770731 0.12842423 0.0001459329     A1L429     -7.134870
#> P0AFG6  2.323891  3.056276 0.14861122 0.0001562791     P0AFG6      2.690083
#>        score_pvalue    score_fdr usedInTest
#> P0AFG3 2.407749e-07 0.0001329935       TRUE
#> P0AC41 3.166968e-07 0.0001329935       TRUE
#> P0A7W1 3.207217e-07 0.0001329935       TRUE
#> P0A7V0 4.111082e-07 0.0001329935       TRUE
#> A1L429 5.272144e-07 0.0001459329       TRUE
#> P0AFG6 7.449241e-07 0.0001562791       TRUE

# The base class itself cannot be instantiated
try(ExploreDataset$new())
#> Error in initialize(...) : 
#>   ExploreDataset is abstract. Use new_explore_dataset() or one of the concrete subclasses.
```
