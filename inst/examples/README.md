# Example datasets

The examples, tests and vignettes use these two files. The app's "Load a demo dataset" button opens the RNA-seq one.

## genomics/deResult.rds

- Content: a `SummarizedExperiment` with the result of a SUSHI DESeq2 two-group analysis (bulk RNA-seq, human, GENCODE GRCh38.p13 release 42, gene level); 21,057 genes × 9 samples, groups `diff` and `undiff`, one contrast `diff--over--undiff`, with its GO over-representation (ORA) and gene set enrichment (GSEA) results
- Assays: `counts`, `xNorm` (DESeq2 median-ratio normalised), `presentFlag`, `isPresent`
- Source: FGCZ SUSHI run of 2024-11-20, B-Fabric project p3009, order o5638
- How it was made: [`inst/scripts/trim_demo_deResult.R`](../scripts/trim_demo_deResult.R) trimmed the SUSHI output to the parts exploreDE reads, so that the file stays under Bioconductor's 5 MB limit
- Size: 4.2 MB

## proteomics/SummarizedExperiment.rds

- Content: a `SummarizedExperiment` with a prolfqua label-free quantification and differential expression result; 1,384 proteins × 6 runs, groups `A` and `B` (factor `G_`), one contrast `B_vs_A`, nested in `rowData` as `constrast_B_vs_A`
- Assays: `rawData` (protein intensities) and `transformedData` (log2, normalised)
- Source: the six Orbitrap all-ion-fragmentation runs `LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01..03` and `..._Condition_B_Sample_Alpha_01..03` of the public LFQ benchmark dataset PRIDE [PXD028735](https://www.ebi.ac.uk/pride/archive/projects/PXD028735) (Van Puyvelde et al. 2022, Scientific Data 9:126, [doi:10.1038/s41597-022-01216-6](https://doi.org/10.1038/s41597-022-01216-6)), quantified and analysed with prolfquapp at the FGCZ
- Size: 0.5 MB
