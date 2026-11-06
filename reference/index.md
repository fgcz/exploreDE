# Package index

## Run the app

- [`exploreDE()`](https://fgcz.github.io/exploreDE/reference/exploreDE.md)
  : Start the exploreDE app

## Build and check exploreDE files

Write an `.h5ad` in the exploreDE/1.0 layout, open it as a dataset, and
check it.

- [`build_explore_h5ad()`](https://fgcz.github.io/exploreDE/reference/build_explore_h5ad.md)
  :

  Build a canonical exploreDE `.h5ad` from user-supplied data

- [`new_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/new_explore_dataset.md)
  :

  Construct an `ExploreDataset` instance from a path, URL, or SE

- [`validate_explore_dataset()`](https://fgcz.github.io/exploreDE/reference/validate_explore_dataset.md)
  :

  Runtime validator for a live `ExploreDataset` instance

## Dataset classes

The accessor API the app reads every dataset through, one class per
storage backend.

- [`ExploreDataset`](https://fgcz.github.io/exploreDE/reference/ExploreDataset.md)
  : Abstract base class for exploreDE's storage-agnostic dataset wrapper

- [`ExploreDatasetH5AD`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetH5AD.md)
  :

  Concrete `ExploreDataset` backend over a canonical AnnData (`.h5ad`)

- [`ExploreDatasetSE`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE.md)
  : Abstract intermediate class for SE-backed datasets

- [`ExploreDatasetSE_CountQC`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE_CountQC.md)
  :

  Concrete subclass for SUSHI CountQC `SummarizedExperiment` outputs

- [`ExploreDatasetSE_CrisprCountQC`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE_CrisprCountQC.md)
  : Concrete subclass for CRISPR Mageck CountQC outputs

- [`ExploreDatasetSE_Proteomics`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE_Proteomics.md)
  :

  Concrete subclass for FGCZ proteomics `SummarizedExperiment` outputs

- [`ExploreDatasetSE_RNASeq`](https://fgcz.github.io/exploreDE/reference/ExploreDatasetSE_RNASeq.md)
  :

  Concrete subclass for SUSHI RNA-seq `SummarizedExperiment` outputs
