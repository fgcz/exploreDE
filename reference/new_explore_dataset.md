# Construct an `ExploreDataset` instance from a path, URL, or SE

Factory function: loads the input (or accepts a live
`SummarizedExperiment`), classifies its shape via `detect_se_shape()`,
and returns the appropriate concrete subclass.

## Usage

``` r
new_explore_dataset(input, config = NULL)
```

## Arguments

- input:

  See above.

- config:

  Optional pre-built config list; advanced use only.

## Value

An object inheriting from
[ExploreDataset](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md).

## Details

Accepted inputs:

- A live `SummarizedExperiment` object.

- A character path to an `.rds` / `.qs2` file.

- A character path to a `.RData` file (Mageck EzResult convention).

- A character path to a directory containing one of the standard SUSHI
  artefacts (`deResult.rds`, `rawData_final.rds`, etc.).

- An `http(s)://` URL pointing at an `.rds` file.

## Examples

``` r
# From a live SummarizedExperiment
se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
eds <- new_explore_dataset(se)
class(eds)[1]
#> [1] "ExploreDatasetSE_Proteomics"

# From an .h5ad file written by build_explore_h5ad()
path <- build_explore_h5ad(
  SummarizedExperiment::assay(se, "rawData"),
  out = tempfile(fileext = ".h5ad"),
  omics = "proteomics",
  coldata = SummarizedExperiment::colData(se)
)
eds_h5ad <- new_explore_dataset(path)
class(eds_h5ad)[1]
#> [1] "ExploreDatasetH5AD"
eds_h5ad$has_de()
#> [1] FALSE
```
