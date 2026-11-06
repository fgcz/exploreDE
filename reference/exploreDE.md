# Start the exploreDE app

Opens the app locally. With no `file`, it starts on the landing page,
where you can upload an `.h5ad` built with
[`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
or load the demo dataset. With `file`, it opens straight into that
dataset.

## Usage

``` r
exploreDE(file = NULL, ...)
```

## Arguments

- file:

  Optional path to an exploreDE `.h5ad`, a SummarizedExperiment
  `.rds`/`.qs2`, or a SUSHI result directory.

- ...:

  Passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) (e.g.
  `port`, `launch.browser`).

## Value

Called to run the app; returns the value passed to
[`shiny::stopApp()`](https://rdrr.io/pkg/shiny/man/stopApp.html) when
the app stops (usually `NULL`), as
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) does.

## Examples

``` r
if (interactive()) {
  exploreDE()
  exploreDE("my_data.h5ad")
}
```
