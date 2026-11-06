# Concrete subclass for SUSHI RNA-seq `SummarizedExperiment` outputs

Reads from a SUSHI ezApp DESeq2 / edgeR / limma result. Top-level
`rowData(se)` carries feature annotations plus inline DE columns for the
single contrast described by `param$comparison`. Optional pre-computed
`enrichResult` / `gseaResult` objects live in `metadata(se)`.

## Value

An R6 class generator; `$new()` returns an `ExploreDatasetSE_RNASeq`
object.

## Details

Inherits the accessor contract from
[ExploreDataset](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md)
via
[ExploreDatasetSE](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.md);
overrides only `.build_config()` and `.prep()`.

## Super classes

[`ExploreDataset`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md)
-\>
[`ExploreDatasetSE`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.md)
-\> `ExploreDatasetSE_RNASeq`

## Methods

### Public methods

- [`ExploreDatasetSE_RNASeq$clone()`](#method-ExploreDatasetSE_RNASeq-clone)

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

### `ExploreDatasetSE_RNASeq$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ExploreDatasetSE_RNASeq$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
se <- readRDS(system.file("examples/genomics/deResult.rds", package = "exploreDE"))
eds <- ExploreDatasetSE_RNASeq$new(se)
eds
#> <ExploreDatasetSE_RNASeq  (se)>
#>   data_type:        RNASeq
#>   feature_id_kind:  gene
#>   organism:         Homo sapiens
#>   ref_build:        Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30
#>   samples x features: 9 x 21057
#>   factors (1):
#>     Condition  (primary): diff, undiff
#>   contrasts (1): diff_vs_undiff
#>   primary_contrast: diff_vs_undiff
#>   assays:           counts, normalised, normalised_log2, vst, fpkm, tpm
#>   enrichment:       ora, gsea
eds$contrasts()
#> [1] "diff_vs_undiff"
head(eds$de(eds$primary_contrast()))
#>                 gene_id
#> DHRS9   ENSG00000073737
#> B4GALT5 ENSG00000158470
#> GCNT3   ENSG00000140297
#> SAA2    ENSG00000134339
#> CXCL6   ENSG00000124875
#> FST     ENSG00000134363
#>                                                                                                                                                                    transcript_id
#> DHRS9                    ENST00000602501; ENST00000357546; ENST00000432060; ENST00000428522; ENST00000450153; ENST00000674881; ENST00000421653; ENST00000436483; ENST00000412271
#> B4GALT5                                                                                                                                                          ENST00000371711
#> GCNT3   ENST00000558721; ENST00000560874; ENST00000560580; ENST00000560111; ENST00000559189; ENST00000396065; ENST00000560585; ENST00000559626; ENST00000560210; ENST00000559200
#> SAA2                                                                        ENST00000414546; ENST00000530400; ENST00000528349; ENST00000529528; ENST00000256733; ENST00000526900
#> CXCL6                                                                                                                          ENST00000515050; ENST00000226317; ENST00000503446
#> FST                                                                                          ENST00000396947; ENST00000256759; ENST00000491717; ENST00000504226; ENST00000497789
#>         gene_name           type       biotypes
#> DHRS9       DHRS9 protein_coding protein_coding
#> B4GALT5   B4GALT5 protein_coding protein_coding
#> GCNT3       GCNT3 protein_coding protein_coding
#> SAA2         SAA2 protein_coding protein_coding
#> CXCL6       CXCL6 protein_coding protein_coding
#> FST           FST protein_coding protein_coding
#>                                               description     gc featWidth
#> DHRS9                           dehydrogenase/reductase 9 0.4434    1615.0
#> B4GALT5                  beta-1,4-galactosyltransferase 5 0.4803    4722.0
#> GCNT3   glucosaminyl (N-acetyl) transferase 3, mucin type 0.4790    1204.0
#> SAA2                                     serum amyloid A2 0.5030     694.3
#> CXCL6                      C-X-C motif chemokine ligand 6 0.4438     878.3
#> FST                                           follistatin 0.4793    1379.0
#>                                                                                                                                                                                                                          GO BP
#> DHRS9                                                                                                                           GO:0006629; GO:0008202; GO:0042572; GO:0030855; GO:0008209; GO:0042448; GO:0042904; GO:0002138
#> B4GALT5                                                 GO:0006629; GO:0005975; GO:0006665; GO:0006486; GO:0031647; GO:0009101; GO:0010706; GO:0016266; GO:0021955; GO:0022010; GO:0030311; GO:0040019; GO:0042551; GO:0070085
#> GCNT3                                                                                                                                       GO:0006486; GO:0048729; GO:0050892; GO:0060993; GO:0005975; GO:0006493; GO:0016266
#> SAA2                                                                                                                                                                                                                GO:0006953
#> CXCL6   GO:0006955; GO:0007165; GO:0006935; GO:0009605; GO:0006952; GO:0060326; GO:0042742; GO:0030593; GO:0071222; GO:0061844; GO:0006954; GO:0007267; GO:0070098; GO:0042119; GO:0001776; GO:0032496; GO:0032642; GO:0070951
#> FST     GO:0000122; GO:0030154; GO:0002244; GO:0030510; GO:0032926; GO:0051798; GO:0001501; GO:0007276; GO:0007389; GO:0008585; GO:0014070; GO:0030509; GO:0030857; GO:0031069; GO:0036305; GO:0042475; GO:0043616; GO:0071363
#>                                                                                                              GO MF
#> DHRS9   GO:0016491; GO:0004745; GO:0016229; GO:0047023; GO:0047035; GO:0047044; GO:0004022; GO:0016854; GO:0005515
#> B4GALT5                                     GO:0046872; GO:0016740; GO:0016757; GO:0003945; GO:0008378; GO:0008489
#> GCNT3               GO:0016757; GO:0016740; GO:0008109; GO:0047225; GO:0008375; GO:0003829; GO:0106325; GO:0106326
#> SAA2                                                                                                    GO:0005515
#> CXCL6                                                   GO:0008009; GO:0005125; GO:0005515; GO:0008201; GO:0045236
#> FST                                                                 GO:0005515; GO:0048185; GO:0038102; GO:0043395
#>                                                              GO CC entrez_id
#> DHRS9   GO:0016020; GO:0005783; GO:0005789; GO:0043231; GO:0030176     10170
#> B4GALT5 GO:0016020; GO:0016021; GO:0005794; GO:0000139; GO:0032580      9334
#> GCNT3   GO:0016020; GO:0016021; GO:0005794; GO:0000139; GO:0070062      9245
#> SAA2                            GO:0005576; GO:0034364; GO:0070062      6289
#> CXCL6                                       GO:0005576; GO:0005615      6372
#> FST                 GO:0005576; GO:0005615; GO:0005634; GO:0005737     10468
#>         isPresentProbe log2Ratio        pValue           fdr usedInTest
#> DHRS9             TRUE  6.325017 2.111083e-293 2.822940e-289       TRUE
#> B4GALT5           TRUE  3.275219 4.239871e-214 2.834778e-210       TRUE
#> GCNT3             TRUE  7.095957 1.716907e-188 7.652826e-185       TRUE
#> SAA2              TRUE  5.864014 7.179677e-188 2.400166e-184       TRUE
#> CXCL6             TRUE  7.121160 3.999797e-167 1.069706e-163       TRUE
#> FST               TRUE -5.094279 1.818724e-165 4.053330e-162       TRUE
#>         feature_id display_name effect_log2fc  score_pvalue     score_fdr
#> DHRS9        DHRS9        DHRS9      6.325017 2.111083e-293 2.822940e-289
#> B4GALT5    B4GALT5      B4GALT5      3.275219 4.239871e-214 2.834778e-210
#> GCNT3        GCNT3        GCNT3      7.095957 1.716907e-188 7.652826e-185
#> SAA2          SAA2         SAA2      5.864014 7.179677e-188 2.400166e-184
#> CXCL6        CXCL6        CXCL6      7.121160 3.999797e-167 1.069706e-163
#> FST            FST          FST     -5.094279 1.818724e-165 4.053330e-162
#>         used_in_test
#> DHRS9           TRUE
#> B4GALT5         TRUE
#> GCNT3           TRUE
#> SAA2            TRUE
#> CXCL6           TRUE
#> FST             TRUE
eds$enrichment_kinds()
#> [1] "ora"  "gsea"
```
