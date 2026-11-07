#' @importFrom dplyr %>%
#' @importFrom R6 R6Class
#' @importFrom rlang .data
NULL

# Packages only the Shiny app uses. inst/app attaches them with library() and calls them unqualified, which
# R CMD check does not see; one import each records the dependency. rhdf5 is never called: anndataR reads
# .h5ad files through it but lists it only in Suggests.
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
