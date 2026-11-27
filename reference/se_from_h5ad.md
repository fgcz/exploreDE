# Rebuild the SummarizedExperiment stored in an exploreDE `.h5ad`

Files written by
[`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
with `mirror = se` hold the whole SummarizedExperiment in `uns["se"]`
(spec section 11). This rebuilds it: Every value at its original path,
with its class and attributes, including S4 objects such as
`DESeqDataSet`, `DGEGLM`, `enrichResult`, and `gseaResult` (their
packages must be installed), functions, and the environments they refer
to.

## Usage

``` r
se_from_h5ad(x)
```

## Arguments

- x:

  Path to an `.h5ad` file, or an anndataR `AnnData` read from one.

## Value

The SummarizedExperiment.
