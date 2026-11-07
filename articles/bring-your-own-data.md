# Bring your own data: building an exploreDE file

exploreDE reads one file format: an `.h5ad` (AnnData) file in the
`exploreDE/1.0` layout. In R,
[`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
writes it for you, whatever pipeline produced your data. This vignette
covers:

1.  installing exploreDE and running it on your own machine;
2.  building a file, from counts only up to counts + DE + pathways;
3.  checking the file before you open it;
4.  what to do when your column names aren’t recognised;
5.  organisms without an annotation package (your own GO terms and
    enrichment);
6.  writing the file from outside R, and saving your settings in the
    app.

Every example below runs as written against the two datasets that ship
with the package (`inst/examples/`).

## Install and run

``` r
install.packages("BiocManager")
BiocManager::install("fgcz/exploreDE")   # pulls in the Bioconductor dependencies
```

Start the app with:

``` r
library(exploreDE)
exploreDE("my_data.h5ad")   # open straight into your file
exploreDE()                 # or start on the landing page
```

On the landing page, use **Choose an exploreDE .h5ad file** to upload
your file (up to 1 GB), or load the demo dataset. The file is read into
memory for your session only. The menu then shows only the tabs your
data supports: DE tabs need a DE table, and pathway tabs need pathway
results.

## What goes into a file

Only three things are required:

| Argument | What it is |
|----|----|
| `x` | the counts: a matrix (features × samples), a named list of matrices (e.g. `list(counts = , normalised = )`), or a `DESeqDataSet` |
| `coldata` | sample metadata, one row per sample, rownames = the counts’ column names |
| `omics` | `"rnaseq"`, `"proteomics"`, or `"crispr"` |

Everything else is optional: `rowdata` (feature annotation), `de` (one
DE table per contrast), `pathways` (ORA/GSEA results per contrast), plus
labels such as `organism`, `ref_build`, and `factors`.

``` r
library(exploreDE)
library(SummarizedExperiment)
out_dir <- tempdir()
```

### Counts only

The shipped RNA-seq example is a SUSHI/DESeq2 result. Using only its
counts and sample table gives the smallest valid file:

``` r
se <- readRDS(system.file("examples/genomics/deResult.rds", package = "exploreDE"))

build_explore_h5ad(
  assay(se, "counts"),
  out     = file.path(out_dir, "counts_only.h5ad"),
  omics   = "rnaseq",
  coldata = colData(se)
)
```

### Proteomics with one contrast

The shipped proteomics example is a prolfqua result. Its DE table sits
in `rowData` as a nested table, with prolfqua’s own column names
(`diff`, `p.value`, `FDR`). You pass it as it is; the column dictionary
(below) recognises those names.

``` r
se_prot <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds",
                               package = "exploreDE"))
de_prot <- as.data.frame(rowData(se_prot)$constrast_B_vs_A)

build_explore_h5ad(
  x       = list(raw = assay(se_prot, "rawData"),
                 transformed = assay(se_prot, "transformedData")),
  out     = file.path(out_dir, "proteomics.h5ad"),
  omics   = "proteomics",
  coldata = colData(se_prot),
  rowdata = de_prot,
  factors = "G_",
  de      = list(B_vs_A = de_prot)
)
```

### RNA-seq with DE and pathways

The same SUSHI result, this time with everything: two count layers,
feature annotation, the DE statistics (SUSHI keeps them in `rowData`)
and the clusterProfiler ORA and GSEA results from its metadata.

``` r
build_explore_h5ad(
  x         = list(counts = assay(se, "counts"), normalised = assay(se, "xNorm")),
  out       = file.path(out_dir, "rnaseq_full.h5ad"),
  omics     = "rnaseq",
  organism  = "Homo sapiens",
  ref_build = metadata(se)$param$refBuild,
  coldata   = colData(se),
  rowdata   = rowData(se),
  factors   = "Condition [Factor]",
  de        = list(diff_vs_undiff = rowData(se)),
  pathways  = list(diff_vs_undiff = list(ora  = metadata(se)$enrichResult,
                                          gsea = metadata(se)$enrichResultGSEA))
)
```

`de` and `pathways` are named by contrast; add one entry per contrast.
Name a contrast `<test>_vs_<reference>`, where both groups are levels of
one `coldata` factor (`diff_vs_undiff` above); other names are an error.

Pathway results can be clusterProfiler objects (`enrichResult`,
`gseaResult`), plain data.frames, or nested lists of either (as above:
BP/MF/CC × up/down/both). A plain data.frame must use clusterProfiler’s
column names: `ID`, `Description`, `p.adjust`, `GeneRatio`, `Count`,
`geneID` for ORA, plus `enrichmentScore`, `NES`, `setSize`,
`core_enrichment` for GSEA. clusterProfiler objects are filtered to
`p.adjust <= enrichment_padj_max` (default 0.1); data.frames are stored
as given. With clusterProfiler objects, the network, ridge and
running-score plots are available too; with data.frames, the table-based
pathway plots are.

## Check the file before you open it

``` r
eds <- new_explore_dataset(file.path(out_dir, "rnaseq_full.h5ad"))
eds
#> <ExploreDatasetH5AD  (h5ad)>
#>   data_type:        RNASeq
#>   feature_id_kind:  gene
#>   organism:         Homo sapiens
#>   ref_build:        Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30
#>   samples x features: 9 x 21057
#>   factors (1):
#>     condition  (primary): diff, undiff
#>   contrasts (1): diff_vs_undiff
#>   primary_contrast: diff_vs_undiff
#>   assays:           counts, normalised, normalised_log2, vst, tpm, fpkm
#>   enrichment:       ora, gsea

eds$has_de()                                   # TRUE: a DE table was supplied
#> [1] TRUE
eds$has_pvalue()                               # TRUE only if a raw p-value was found
#> [1] TRUE
d <- eds$de(eds$primary_contrast())
c(log2fc_missing = all(is.na(d$effect_log2fc)),
  fdr_missing    = all(is.na(d$score_fdr)))   # both should be FALSE
#> log2fc_missing    fdr_missing 
#>          FALSE          FALSE
```

If a DE column wasn’t recognised it comes through as `NA` rather than as
an error, so do this check before trusting the file.

## How columns are recognised

You don’t tell
[`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
which column is the log2 fold change, the p-value or the FDR. A
dictionary shipped with the package
(`inst/dictionaries/column_patterns.json`) holds an ordered list of
case-insensitive patterns per concept. The first pattern that matches
any column wins, so exact names beat loose substrings:

``` json
"effect_log2fc": ["^diff$", "^log2[._ ]?(fold)?[._ ]?(change|ratio|fc)?$",
                  "^logfc$", "^estimate$", "log2"],
"score_pvalue":  ["^p[._ ]?value$"],
"score_fdr":     ["^fdr$", "^adj[._ ]?p[._ ]?val(ue)?$", "^padj$",
                  "^q[._ ]?val(ue)?$", "^p[._ ]?adjusted$", "fdr"],
"score":         ["score"],
"used_in_test":  ["^used[._ ]?in[._ ]?test$"]
```

This covers DESeq2, edgeR, limma, and prolfqua out of the box. It also
covers SAINTexpress, whose tables have `log2_EFCs`, `BFDR`, and
`SaintScore` and no p-value at all: `log2_EFCs` matches `log2`, `BFDR`
matches `fdr`, and `SaintScore` becomes an extra score. That score is
shown in the DE table and offered as a volcano y-axis; significance
thresholds still use the FDR. With no p-value, the app works from the
FDR alone.

The `rowdata` side works the same way for gene IDs, symbols,
descriptions, protein and peptide IDs, phospho sites, and GO columns.

Only the recognised statistics are stored for each contrast, plus a
`usedInTest` flag when the table has one (SUSHI’s expression filter; the
DE views then show only the tested features). Without it, a feature
counts as tested when its p-value isn’t `NA`. Other DE-table columns
(e.g. a per-contrast `avgAbd`) are not kept; put per-feature annotation
in `rowdata`, where every column is kept.

### When a column isn’t recognised

[`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
warns, naming the contrast, the field it couldn’t fill, and the columns
it saw. Two fixes:

- **For one file:** name the columns yourself. Five keys exist; anything
  you don’t name is still found by the dictionary.

  ``` r
  build_explore_h5ad(...,
    de_columns      = list(effect_log2fc = "<your column>",
                           score_pvalue  = "<your column>",
                           score_fdr     = "<your column>"),
    rowdata_columns = list(id_col = "<your column>", symbol_col = "<your column>"))
  ```

  `de_columns` applies to every contrast, or pass a list of such lists
  named by contrast. The mapping is saved in the file.

- **For a naming convention you’ll see again:** add a pattern to
  `inst/dictionaries/column_patterns.json`. It’s a JSON edit, not a code
  change; put anchored patterns (`^...$`) before loose ones.

### Feature level

`feature_level` is detected from `rowdata`: a `site` column means
phospho, a `peptide_Id` column means peptide, and a `transcript_id`
column with one transcript per row means isoform. Otherwise it’s gene
(RNA-seq) or protein (proteomics). Pass `feature_level =` to override.

## Organisms without an annotation package

Nothing in the file needs an organism database. Leave `organism` and
`ref_build` out, and use your own gene ids as the rownames of the
counts. What changes is what the app can offer:

- **KEGG** tabs appear only for human, mouse, and rat.
- **GO** tile plots and the GO heatmap appear when `rowdata` has GO
  annotation columns: Any of `GO BP`, `GO MF`, and `GO CC`, each holding
  `"; "`-separated GO ids (`NA` for none). One ontology is enough.
- **ORA / GSEA** show the results you pass in `pathways`; the app never
  runs enrichment itself.

If your annotation is one list of GO ids per gene (e.g. eggNOG-mapper’s
comma-separated `GOs` column), split it by ontology with GO.db:

``` r
# annot: one row per gene, with gene_id and GOs ("GO:0008150,GO:0003674", or "-")
ids  <- strsplit(annot$GOs, ",", fixed = TRUE)
long <- data.frame(gene_id = rep(annot$gene_id, lengths(ids)), go = unlist(ids))
long <- long[grepl("^GO:", long$go), ]
long$ontology <- AnnotationDbi::Ontology(long$go)   # "BP", "MF" or "CC"
for (ont in c("BP", "MF", "CC")) {
  hits     <- long[long$ontology %in% ont, ]
  per_gene <- split(hits$go, factor(hits$gene_id, levels = annot$gene_id))
  annot[[paste("GO", ont)]] <- vapply(per_gene, function(x)
    if (length(x)) paste(x, collapse = "; ") else NA_character_, "")
}
```

For pathway results, run clusterProfiler’s `enricher()` (or `GSEA()`)
with a `TERM2GENE` built from the same annotation. The example below
does this on the shipped dataset as if it had no annotation package: It
keeps only the `GO BP` column, builds `TERM2GENE` from it, and runs ORA
on the up-regulated genes against the tested ones. Gene ids in
`TERM2GENE` must be the rownames of your counts.

``` r
rd  <- rowData(se)
bp  <- strsplit(as.character(rd[["GO BP"]]), "; ", fixed = TRUE)
t2g <- data.frame(term = unlist(bp), gene = rep(rownames(rd), lengths(bp)))
t2g <- t2g[!is.na(t2g$term), ]

tested <- rownames(rd)[rd$usedInTest]
up     <- rownames(rd)[rd$usedInTest & !is.na(rd$fdr) & rd$fdr < 0.05 & rd$log2Ratio > 0]
ora_up <- clusterProfiler::enricher(up, TERM2GENE = t2g, universe = tested)

build_explore_h5ad(
  x        = list(counts = assay(se, "counts"), normalised = assay(se, "xNorm")),
  out      = file.path(out_dir, "own_annotation.h5ad"),
  omics    = "rnaseq",
  coldata  = colData(se),
  rowdata  = rd[, c("gene_id", "gene_name", "description", "GO BP")],
  factors  = "Condition [Factor]",
  de       = list(diff_vs_undiff = rd),
  pathways = list(diff_vs_undiff = list(ora = list(BP = list(upGenes = ora_up))))
)
```

Name the nested results as the app shows them: The ontology (`BP`, `MF`,
`CC`, or your own label), then the direction (`upGenes`, `downGenes`,
`bothGenes`).

## Writing the file outside R

If your pipeline isn’t in R (a Python script, an exporter app), write
the `.h5ad` yourself following the format spec that ships with the
package:

``` r
system.file("spec/exploreDE-h5ad-1.0.md", package = "exploreDE")
```

It lists every required key, the proteomics profile (protein, peptide,
phospho, and precursor levels), and a minimal Python writer that has
been checked against the reader. Check the result with
[`new_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/new_explore_dataset.md)
as above. APB Export’s exploreDE output and the prolfquapp DEA’s
`AnnData.h5ad` (prolfquapp 2.10.x) open as they are.

## Opening a file from a link

On a server, `?data=` takes a percent-encoded `https://` link to an
`.h5ad` from an allowed source, for example
`https://shiny-public.fgcz.uzh.ch/app/exploreDE/?data=https%3A%2F%2F...%2Fmy_data.h5ad`.
The default sources are APB Export jobs and public pStore; set the
environment variable `EXPLOREDE_URL_SOURCES` (comma-separated `https://`
prefixes) to change them. The file is downloaded without redirects, up
to 1 GB, and then checked like an upload.

## Saving your settings

The **Save / Load Session** tab downloads your settings across every tab
(widgets, the feature bucket, pathway selections, and Multi DE datasets)
as one `.json` file. Upload it later on the same dataset to pick up
where you left off; a session from a different dataset is refused. The
file holds no data, and its format is described in
`system.file("spec/exploreDE-session-2.md", package = "exploreDE")`.

## Session info

``` r
sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] SummarizedExperiment_1.43.0 Biobase_2.73.2             
#>  [3] GenomicRanges_1.65.4        Seqinfo_1.3.2              
#>  [5] IRanges_2.47.5              S4Vectors_0.51.10          
#>  [7] BiocGenerics_0.59.12        generics_0.1.4             
#>  [9] MatrixGenerics_1.25.0       matrixStats_1.5.0          
#> [11] exploreDE_3.10.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.6.2               later_1.4.8                
#>   [3] ggplotify_0.1.4             tibble_3.3.1               
#>   [5] polyclip_1.10-7             enrichit_0.2.7             
#>   [7] httr2_1.3.0                 lifecycle_1.0.5            
#>   [9] rstatix_1.1.0               rprojroot_2.1.1            
#>  [11] processx_3.9.0              lattice_0.23-1             
#>  [13] MASS_7.3-66                 backports_1.5.1            
#>  [15] magrittr_2.0.5              plotly_4.12.1              
#>  [17] sass_0.4.10                 rmarkdown_2.32             
#>  [19] jquerylib_0.1.4             yaml_2.3.12                
#>  [21] httpuv_1.6.17               otel_0.2.0                 
#>  [23] ggtangle_0.1.5              zip_3.0.2                  
#>  [25] askpass_1.2.1               reticulate_1.47.0          
#>  [27] ggvenn_0.1.19               cowplot_1.2.0              
#>  [29] DBI_1.3.0                   RColorBrewer_1.1-3         
#>  [31] abind_1.4-8                 purrr_1.2.2                
#>  [33] yulab.utils_0.2.5           tweenr_2.0.3               
#>  [35] rappdirs_0.3.4              aisdk_1.4.12               
#>  [37] gdtools_0.5.3               enrichplot_1.99.8          
#>  [39] ggrepel_0.9.8               tidytree_0.4.8             
#>  [41] pkgdown_2.2.1               svglite_2.2.2              
#>  [43] codetools_0.2-20            DelayedArray_0.39.8        
#>  [45] scuttle_1.23.2              DOSE_4.7.3                 
#>  [47] DT_0.34.0                   xml2_1.6.0                 
#>  [49] ggforce_0.5.0               tidyselect_1.2.1           
#>  [51] aplot_0.3.2                 farver_2.1.2               
#>  [53] jsonlite_2.0.0              waiter_0.2.5.1             
#>  [55] Formula_1.2-6               systemfonts_1.3.2          
#>  [57] tools_4.6.2                 ggnewscale_0.5.2           
#>  [59] treeio_1.37.1               ragg_1.5.2                 
#>  [61] Rcpp_1.1.2                  glue_1.8.1                 
#>  [63] SparseArray_1.13.4          BiocBaseUtils_1.15.1       
#>  [65] xfun_0.62                   qvalue_2.45.0              
#>  [67] dplyr_1.2.1                 shinydashboard_0.7.3       
#>  [69] withr_3.0.3                 fastmap_1.2.0              
#>  [71] GGally_2.4.0                rhdf5filters_1.25.5        
#>  [73] shinyjs_2.1.1               openssl_2.5.0              
#>  [75] callr_3.8.0                 digest_0.6.39              
#>  [77] R6_2.6.1                    mime_0.13                  
#>  [79] gridGraphics_0.5-1          textshaping_1.0.5          
#>  [81] colorspace_2.1-4            GO.db_3.24.0               
#>  [83] gtools_3.9.5                RSQLite_3.53.3             
#>  [85] tidyr_1.3.2                 fontLiberation_0.1.0       
#>  [87] data.table_1.18.6.1         S4Arrays_1.13.2            
#>  [89] httr_1.4.9                  htmlwidgets_1.6.4          
#>  [91] scatterpie_0.2.6            ggstats_0.14.0             
#>  [93] pkgconfig_2.0.3             gtable_0.3.6               
#>  [95] blob_1.3.0                  S7_0.2.2                   
#>  [97] SingleCellExperiment_1.35.2 XVector_0.53.0             
#>  [99] clusterProfiler_4.21.4      htmltools_0.5.9            
#> [101] fontBitstreamVera_0.1.1     carData_3.0-6              
#> [103] scales_1.4.0                kableExtra_1.4.1           
#> [105] fresh_0.2.2                 png_0.1-9                  
#> [107] ggfun_0.2.2                 knitr_1.52                 
#> [109] rstudioapi_0.19.0           reshape2_1.4.5             
#> [111] uuid_1.2-2                  nlme_3.1-171               
#> [113] anndataR_1.3.4              cachem_1.1.0               
#> [115] rhdf5_2.57.20               stringr_1.6.0              
#> [117] parallel_4.6.2              miniUI_0.1.2               
#> [119] vipor_0.4.7                 shinycssloaders_1.1.0      
#> [121] AnnotationDbi_1.75.2        desc_1.4.3                 
#> [123] pillar_1.11.1               grid_4.6.2                 
#> [125] vctrs_0.7.3                 promises_1.5.0             
#> [127] ggpubr_1.0.0                car_3.1-5                  
#> [129] tidydr_0.0.7                beachmat_2.29.3            
#> [131] xtable_1.8-8                cluster_2.1.8.3            
#> [133] beeswarm_0.4.0              evaluate_1.0.5             
#> [135] cli_3.6.6                   compiler_4.6.2             
#> [137] rlang_1.3.0                 crayon_1.5.3               
#> [139] sortable_0.6.0              ggsignif_0.6.4             
#> [141] ps_1.9.3                    plyr_1.8.9                 
#> [143] fs_2.1.0                    ggbeeswarm_0.7.3           
#> [145] ggiraph_0.9.6               stringi_1.8.9              
#> [147] BiocParallel_1.47.0         viridisLite_0.4.3          
#> [149] assertthat_0.2.1            Biostrings_2.81.9          
#> [151] lazyeval_0.2.3              colourpicker_1.3.0         
#> [153] GOSemSim_2.39.4             fontquiver_0.2.1           
#> [155] Matrix_1.7-6                patchwork_1.3.2            
#> [157] bit64_4.8.6                 ggplot2_4.0.3              
#> [159] learnr_0.11.6               Rhdf5lib_2.1.0             
#> [161] KEGGREST_1.53.6             shiny_1.14.0               
#> [163] igraph_2.3.4                broom_1.0.13               
#> [165] memoise_2.0.1               bslib_0.12.0               
#> [167] ggtree_4.3.0                bit_4.6.0                  
#> [169] officer_0.7.7               gson_0.2.2                 
#> [171] ape_5.8-1
```
