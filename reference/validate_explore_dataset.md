# Runtime validator for a live `ExploreDataset` instance

Checks subclass-agnostic invariants by walking the accessor API:

## Usage

``` r
validate_explore_dataset(eds, strict = TRUE)
```

## Arguments

- eds:

  An `ExploreDataset` instance.

- strict:

  If `TRUE` (default), [`stop()`](https://rdrr.io/r/base/stop.html) on
  the first failure. If `FALSE`, returns `list(ok, errors)`.

## Value

Invisibly, `TRUE` on success (strict mode); or
`list(ok = logical(1), errors = character)` (non-strict mode).

## Details

1.  The runtime class is a concrete `ExploreDataset` subclass — not the
    abstract base or intermediate.

2.  Cross-references resolve: `factor_names()` are columns of
    [`samples()`](https://rdrr.io/pkg/Biobase/man/class.MIAME.html);
    every [`contrasts()`](https://rdrr.io/r/stats/contrasts.html) entry
    has a working [`de()`](https://rdrr.io/r/utils/dataentry.html);
    `primary_grouping()` is a member of `factor_names()`;
    `primary_contrast()` is in
    [`contrasts()`](https://rdrr.io/r/stats/contrasts.html) (or `NA`).

3.  Sample / feature alignment: assay rownames / colnames agree with
    `sample_ids()` / `feature_ids()`.

Cost is proportional to the number of contrasts (one `de(c)` round-trip
each) plus a small constant of other accessor calls.

## Examples

``` r
se <- readRDS(system.file("examples/proteomics/SummarizedExperiment.rds", package = "exploreDE"))
eds <- new_explore_dataset(se)
validate_explore_dataset(eds)
validate_explore_dataset(eds, strict = FALSE)
#> $ok
#> [1] TRUE
#> 
#> $errors
#> character(0)
#> 
```
