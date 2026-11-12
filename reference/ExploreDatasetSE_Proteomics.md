# Concrete subclass for FGCZ proteomics `SummarizedExperiment` outputs

Reads from a prolfqua-derived SE. Top-level `rowData(se)` is sparse;
feature metadata and per-contrast DE statistics live nested inside
`rowData(se)$constrast_<name>` DataFrames (note the SUSHI typo is
preserved). The first contrast is the source of truth for feature
metadata. `feature_id` is **built** from `protein_Id` (+ `site` /
`peptide_Id` to disambiguate), and the assay matrices' rownames are
remapped via `.fid_lookup` since SE rownames do not match the built
`feature_id`s.

## Value

An R6 class generator; `$new()` returns an `ExploreDatasetSE_Proteomics`
object.

## Details

Inherits the accessor contract from
[ExploreDataset](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md)
via
[ExploreDatasetSE](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.md);
overrides `.build_config()` and `.prep()`.

## Super classes

[`ExploreDataset`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md)
-\>
[`ExploreDatasetSE`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.md)
-\> `ExploreDatasetSE_Proteomics`

## Methods

### Public methods

- [`ExploreDatasetSE_Proteomics$clone()`](#method-ExploreDatasetSE_Proteomics-clone)

Inherited methods

- [`ExploreDataset$all_pathways()`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.html#method-all_pathways)
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

### `ExploreDatasetSE_Proteomics$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ExploreDatasetSE_Proteomics$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
eds <- ExploreDatasetSE_Proteomics$new(se)
eds
#> <ExploreDatasetSE_Proteomics  (se)>
#>   data_type:        proteomics
#>   feature_id_kind:  protein
#>   organism:         NA
#>   ref_build:        NA
#>   samples x features: 6 x 1384
#>   factors (3):
#>     G_  (primary): A, B
#>     CONTROL: C, T
#>     isotopeLabel: light
#>   contrasts (1): B_vs_A
#>   primary_contrast: B_vs_A
#>   assays:           intensities
#>   enrichment:       (none)
head(eds$feature_ids())
#> [1] "A0A0U1RRI6" "A1L429"     "A2RU48"     "A6NCE7"     "A6NEC2"    
#> [6] "A6NEK1"    
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
```
