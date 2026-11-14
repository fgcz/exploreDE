#' exploreDE: interactive exploration of differential expression results
#'
#' A Shiny app for exploring differential expression results from RNA-seq,
#' proteomics and CRISPR screens, and for making figures from them.
#'
#' - [exploreDE()] starts the app, on a file or a `SummarizedExperiment` or
#'   AnnData object.
#' - [build_explore_h5ad()] writes an exploreDE `.h5ad` from any pipeline's
#'   results, and [validate_explore_dataset()] checks it.
#' - [new_explore_dataset()] opens a dataset in R, through the same interface
#'   the app's tabs use.
#'
#' Start with `vignette("exploreDE")`; `vignette("bring-your-own-data")`
#' explains how to build a file from your own results.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom R6 R6Class
#' @importFrom rlang .data
NULL

# Packages only the Shiny app uses. inst/app attaches them with library() and calls them unqualified, which
# R CMD check does not see; one import each records the dependency. rhdf5 is never called: anndataR reads
# .h5ad files through it but lists it only in Suggests.
#' @importFrom clusterProfiler go2term
#' @importFrom colorspace lighten
#' @importFrom colourpicker colourInput
#' @importFrom cowplot plot_grid
#' @importFrom DT datatable
#' @importFrom enrichplot gseaplot2
#' @importFrom fresh create_theme
#' @importFrom GGally ggpairs
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom ggpubr stat_cor
#' @importFrom ggvenn ggvenn
#' @importFrom kableExtra kable_styling
#' @importFrom matrixStats rowMedians
#' @importFrom officer read_docx
#' @importFrom parallel detectCores
#' @importFrom plotly plot_ly
#' @importFrom rhdf5 h5ls
#' @importFrom rstatix add_significance
#' @importFrom scales hue_pal
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinydashboard dashboardPage
#' @importFrom shinyjs useShinyjs
#' @importFrom sortable bucket_list
#' @importFrom tidyr pivot_longer
#' @importFrom waiter use_waiter
NULL
