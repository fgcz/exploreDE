# Export a canonical `.h5ad` to the parquet directory exploreDElight reads

Writes the directory layout consumed by
[ExploreDatasetParquet](https://fgcz.github.io/exploreDE/reference/ExploreDatasetParquet.md):
`manifest.json`, `obs.parquet`, `var.parquet`, `uns.json`,
`layers/<name>.parquet` and `varm/<slot>.parquet`.

## Usage

``` r
export_h5ad_to_parquet(input, out, quiet = TRUE)
```

## Arguments

- input:

  Path to a canonical exploreDE/1.0 `.h5ad`, or a live `anndataR`
  object.

- out:

  Output directory; created if absent.

- quiet:

  Suppress per-layer progress messages.

## Value

`out`, invisibly.
