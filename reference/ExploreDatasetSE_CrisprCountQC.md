# Concrete subclass for CRISPR Mageck CountQC outputs

Identical to
[ExploreDatasetSE_CountQC](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE_CountQC.md)
plus a sgRNA → gene mapping parsed from the SE feature IDs, which are
assumed to follow the `GENE_<n>_<SEQUENCE>` convention.

## Value

An R6 class generator; `$new()` returns an
`ExploreDatasetSE_CrisprCountQC` object.

## Details

Inherits `.prep()` from
[ExploreDatasetSE_CountQC](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE_CountQC.md)
unchanged.

## Super classes

[`ExploreDataset`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md)
-\>
[`ExploreDatasetSE`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.md)
-\>
[`ExploreDatasetSE_CountQC`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE_CountQC.md)
-\> `ExploreDatasetSE_CrisprCountQC`

## Methods

### Public methods

- [`ExploreDatasetSE_CrisprCountQC$clone()`](#method-ExploreDatasetSE_CrisprCountQC-clone)

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
- [`ExploreDatasetSE$assay()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.html#method-assay)
- [`ExploreDatasetSE$assay_names()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.html#method-assay_names)
- [`ExploreDatasetSE$backend()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.html#method-backend)
- [`ExploreDatasetSE$backend_kind()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.html#method-backend_kind)
- [`ExploreDatasetSE$de()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.html#method-de)
- [`ExploreDatasetSE$enrichment()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.html#method-enrichment)
- [`ExploreDatasetSE$enrichment_input()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.html#method-enrichment_input)
- [`ExploreDatasetSE$features()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.html#method-features)
- [`ExploreDatasetSE$initialize()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.html#method-initialize)
- [`ExploreDatasetSE$samples()`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.html#method-samples)

------------------------------------------------------------------------

### `ExploreDatasetSE_CrisprCountQC$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ExploreDatasetSE_CrisprCountQC$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
library(SummarizedExperiment)
# Guide counts: 4 guides for each of 10 genes, named GENE_<n>_<SEQUENCE>
set.seed(1)
guides <- sprintf("GENE%d_%d_ACGTACGTACGT", rep(1:10, each = 4), rep(1:4, 10))
counts <- matrix(rpois(40 * 4, 50), nrow = 40, dimnames = list(guides, paste0("S", 1:4)))
samples <- S4Vectors::DataFrame(
  `Condition [Factor]` = rep(c("ctrl", "treated"), each = 2),
  row.names = colnames(counts),
  check.names = FALSE
)
se <- SummarizedExperiment(assays = list(counts = counts), colData = samples)
eds <- ExploreDatasetSE_CrisprCountQC$new(se)
eds
#> <ExploreDatasetSE_CrisprCountQC  (se)>
#>   data_type:        CrisprCountQC
#>   feature_id_kind:  guide
#>   organism:         NA
#>   ref_build:        NA
#>   samples x features: 4 x 40
#>   factors (1):
#>     Condition  (primary): ctrl, treated
#>   contrasts (0): 
#>   primary_contrast: (none)
#>   assays:           counts, vst
#>   enrichment:       (none)
head(eds$sgrna_to_gene())
#> GENE1_1_ACGTACGTACGT GENE1_2_ACGTACGTACGT GENE1_3_ACGTACGTACGT 
#>              "GENE1"              "GENE1"              "GENE1" 
#> GENE1_4_ACGTACGTACGT GENE2_1_ACGTACGTACGT GENE2_2_ACGTACGTACGT 
#>              "GENE1"              "GENE2"              "GENE2" 
```
