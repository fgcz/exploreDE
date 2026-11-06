# Concrete `ExploreDataset` backend over a canonical AnnData (`.h5ad`)

Reads a canonical exploreDE/1.0 AnnData and implements the accessor
contract inherited from
[ExploreDataset](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md).
Accepts a path to a `.h5ad` file or a live `anndataR` object. This is
the backend for user-supplied data ("bring your own data") and the
eventual single backend once producers emit canonical `.h5ad` directly.

## Value

An R6 class generator; `$new()` returns an `ExploreDatasetH5AD` object.

## Super class

[`ExploreDataset`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md)
-\> `ExploreDatasetH5AD`

## Methods

### Public methods

- [`ExploreDatasetH5AD$new()`](#method-ExploreDatasetH5AD-initialize)

- [`ExploreDatasetH5AD$backend()`](#method-ExploreDatasetH5AD-backend)

- [`ExploreDatasetH5AD$backend_kind()`](#method-ExploreDatasetH5AD-backend_kind)

- [`ExploreDatasetH5AD$samples()`](#method-ExploreDatasetH5AD-samples)

- [`ExploreDatasetH5AD$features()`](#method-ExploreDatasetH5AD-features)

- [`ExploreDatasetH5AD$de()`](#method-ExploreDatasetH5AD-de)

- [`ExploreDatasetH5AD$assay_names()`](#method-ExploreDatasetH5AD-assay_names)

- [`ExploreDatasetH5AD$assay()`](#method-ExploreDatasetH5AD-assay)

- [`ExploreDatasetH5AD$enrichment()`](#method-ExploreDatasetH5AD-enrichment)

- [`ExploreDatasetH5AD$enrichment_table()`](#method-ExploreDatasetH5AD-enrichment_table)

- [`ExploreDatasetH5AD$enrichment_input()`](#method-ExploreDatasetH5AD-enrichment_input)

- [`ExploreDatasetH5AD$clone()`](#method-ExploreDatasetH5AD-clone)

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

### `ExploreDatasetH5AD$new()`

Read a canonical `.h5ad` or take a live AnnData, build the config from
it unless one is given, and prepare the feature and DE tables.

#### Usage

    ExploreDatasetH5AD$new(input, config = NULL)

#### Arguments

- `input`:

  An `anndataR` AnnData object, or a path to an existing `.h5ad` file.

- `config`:

  Optional pre-built config list; built from the AnnData when `NULL`.

------------------------------------------------------------------------

### `ExploreDatasetH5AD$backend()`

The AnnData parts read at construction, as
`list(obs, var, layers, varm, uns)` of plain R objects.

#### Usage

    ExploreDatasetH5AD$backend()

------------------------------------------------------------------------

### `ExploreDatasetH5AD$backend_kind()`

Returns `"h5ad"`.

#### Usage

    ExploreDatasetH5AD$backend_kind()

------------------------------------------------------------------------

### `ExploreDatasetH5AD$samples()`

Sample metadata from `obs`, with the `factor_` prefix stripped from
column names and a `Name` column added when missing.

#### Usage

    ExploreDatasetH5AD$samples()

------------------------------------------------------------------------

### `ExploreDatasetH5AD$features()`

Feature annotation from `var`, with `feature_id` and `display_name`
columns, ordered by the p-value of the primary contrast.

#### Usage

    ExploreDatasetH5AD$features()

------------------------------------------------------------------------

### `ExploreDatasetH5AD$de()`

DE table of one contrast from `varm`, with the legacy aliases
`log2Ratio`, `pValue` and `fdr` added; an empty data.frame when the slot
is missing. Errors on an unknown contrast.

#### Usage

    ExploreDatasetH5AD$de(contrast)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

------------------------------------------------------------------------

### `ExploreDatasetH5AD$assay_names()`

Canonical assay names from the config's assay map whose layer exists;
computed assays are always included.

#### Usage

    ExploreDatasetH5AD$assay_names()

------------------------------------------------------------------------

### `ExploreDatasetH5AD$assay()`

One assay as a features x samples matrix (the layer transposed);
`normalised`, `normalised_log2`, `vst`, `fpkm` and `tpm` are computed
when the config maps them to no layer.

#### Usage

    ExploreDatasetH5AD$assay(name, features = NULL, samples = NULL)

#### Arguments

- `name`:

  Canonical assay name, one of `assay_names()`.

- `features`:

  Feature IDs (rows) to keep; `NULL` keeps all.

- `samples`:

  Sample IDs (columns) to keep; `NULL` keeps all.

------------------------------------------------------------------------

### `ExploreDatasetH5AD$enrichment()`

Enrichment results of one contrast and kind: decoded from the serialised
objects in older files, otherwise rebuilt from `uns$enrichment` and
`uns$enrichment_support`.

#### Usage

    ExploreDatasetH5AD$enrichment(contrast, kind)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

- `kind`:

  `"ora"` or `"gsea"`.

------------------------------------------------------------------------

### `ExploreDatasetH5AD$enrichment_table()`

Tidy enrichment table of one contrast and kind: the table stored in
`uns$enrichment` when present, otherwise derived from `enrichment()` as
in the base class.

#### Usage

    ExploreDatasetH5AD$enrichment_table(contrast, kind)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

- `kind`:

  `"ora"` or `"gsea"`.

------------------------------------------------------------------------

### `ExploreDatasetH5AD$enrichment_input()`

Enrichment inputs rebuilt from the DE table: `features()`, log2 fold
changes sorted decreasing, and up/down feature selections at the
contrast's FDR and log2FC thresholds.

#### Usage

    ExploreDatasetH5AD$enrichment_input(contrast)

#### Arguments

- `contrast`:

  Contrast name, one of
  [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

------------------------------------------------------------------------

### `ExploreDatasetH5AD$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ExploreDatasetH5AD$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following object is masked from ‘package:utils’:
#> 
#>     data
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, scale, sequence, table,
#>     tapply, transform, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
de <- as.data.frame(rowData(se)$constrast_B_vs_A)
path <- tempfile(fileext = ".h5ad")
build_explore_h5ad(
  x = list(raw = assay(se, "rawData"), transformed = assay(se, "transformedData")),
  out = path,
  omics = "proteomics",
  coldata = colData(se),
  rowdata = de,
  factors = "G_",
  de = list(B_vs_A = de)
)
eds <- ExploreDatasetH5AD$new(path)
eds
#> <ExploreDatasetH5AD  (h5ad)>
#>   data_type:        proteomics
#>   feature_id_kind:  protein
#>   organism:         NA
#>   ref_build:        NA
#>   samples x features: 6 x 1384
#>   factors (1):
#>     g_  (primary): A, B
#>   contrasts (1): B_vs_A
#>   primary_contrast: B_vs_A
#>   assays:           intensities, intensities_raw
#>   enrichment:       (none)
head(eds$de("B_vs_A"))
#>        effect_log2fc score_pvalue    score_fdr feature_id log2Ratio
#> P0AFG3      2.076169 2.407749e-07 0.0001329935     P0AFG3  2.076169
#> P0AC41      2.207672 3.166968e-07 0.0001329935     P0AC41  2.207672
#> P0A7W1      2.232960 3.207217e-07 0.0001329935     P0A7W1  2.232960
#> P0A7V0      2.290956 4.111082e-07 0.0001329935     P0A7V0  2.290956
#> A1L429     -7.134870 5.272144e-07 0.0001459329     A1L429 -7.134870
#> P0AFG6      2.690083 7.449241e-07 0.0001562791     P0AFG6  2.690083
#>              pValue          fdr gene_name usedInTest gene_id
#> P0AFG3 2.407749e-07 0.0001329935    P0AFG3       TRUE  P0AFG3
#> P0AC41 3.166968e-07 0.0001329935    P0AC41       TRUE  P0AC41
#> P0A7W1 3.207217e-07 0.0001329935    P0A7W1       TRUE  P0A7W1
#> P0A7V0 4.111082e-07 0.0001329935    P0A7V0       TRUE  P0A7V0
#> A1L429 5.272144e-07 0.0001459329    A1L429       TRUE  A1L429
#> P0AFG6 7.449241e-07 0.0001562791    P0AFG6       TRUE  P0AFG6
#>                  description
#> P0AFG3  sp_P0AFG3_ODO1_ECOLI
#> P0AC41  sp_P0AC41_SDHA_ECOLI
#> P0A7W1   sp_P0A7W1_RS5_ECOLI
#> P0A7V0   sp_P0A7V0_RS2_ECOLI
#> A1L429 sp_A1L429_GG12C_HUMAN
#> P0AFG6  sp_P0AFG6_ODO2_ECOLI
```
