# exploreDE: interactive exploration of differential expression results

A Shiny app for exploring differential expression results from RNA-seq,
proteomics and CRISPR screens, and for making figures from them.

## Details

- [`exploreDE()`](https://fgcz.github.io/exploreDE/reference/exploreDE.md)
  starts the app, on a file or a `SummarizedExperiment` or AnnData
  object.

- [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  writes an exploreDE `.h5ad` from any pipeline's results, and
  [`validate_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/validate_explore_dataset.md)
  checks it.

- [`new_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/new_explore_dataset.md)
  opens a dataset in R, through the same interface the app's tabs use.

Start with
[`vignette("exploreDE")`](https://fgcz.github.io/exploreDE/articles/exploreDE.md);
[`vignette("bring-your-own-data")`](https://fgcz.github.io/exploreDE/articles/bring-your-own-data.md)
explains how to build a file from your own results.

## See also

Useful links:

- <https://fgcz.github.io/exploreDE/>

- <https://github.com/fgcz/exploreDE>

- Report bugs at <https://github.com/fgcz/exploreDE/issues>

## Author

**Maintainer**: Witold Wolski <wewolski@gmail.com>
([ORCID](https://orcid.org/0000-0002-6468-120X))

Authors:

- Peter Leary <peter.leary@uzh.ch>
  ([ORCID](https://orcid.org/0000-0003-1430-6185))

Other contributors:

- Functional Genomics Center Zurich \[funder\]
