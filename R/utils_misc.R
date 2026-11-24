# Small app utilities: URLs, boxplot summaries, colour conversion.

#' is_url
#' Test whether a dataDir string is a URL (http/https) rather than a local path.
#' @param dataDir Character string
#' @return Logical scalar
#' @noRd
is_url <- function(dataDir) {
  return(grepl("^https?://", dataDir))
}

# Telemetry app name for the `dataType` a ?data= link set (NULL without one).
telemetry_name <- function(data_type) {
  switch(
    data_type %||% "",
    proteomics = "exploreDE_proteomics",
    metabolomics = "exploreDE_metabolomics",
    "exploreDE_genomics"
  )
}

# What the UI calls feature values: metabolomics measures intensities, not counts.
.values_words <- c(
  "Counts" = "Intensities",
  "Count" = "Intensity",
  "counts" = "intensities",
  "count" = "intensity",
  "count method" = "intensity type"
)
values_word <- function(data_type, word = "Counts") {
  if (identical(data_type, "metabolomics")) .values_words[[word]] else word
}
values_class <- function(word) paste0("values-word-", gsub(" ", "-", word, fixed = TRUE))
# A word in the static UI that server-summary.R relabels once the data type is known.
values_span <- function(word) shiny::span(class = values_class(word), word, .noWS = "outside")

#' dataSummary
#' Summarise counts by mean / SD / SE across grouping variables.
#' @param data data.frame containing at least `varname` and the `groupnames` columns
#' @param varname Column to summarise
#' @param groupnames Character vector of grouping columns
#' @return data.frame with one row per group combination
#' @noRd
dataSummary <- function(data, varname, groupnames) {
  data_sum <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(groupnames))) |>
    dplyr::summarise(
      mean = mean(.data[[varname]], na.rm = TRUE),
      sd = stats::sd(.data[[varname]], na.rm = TRUE),
      se = stats::sd(.data[[varname]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[varname]]))),
      .groups = "drop"
    ) |>
    as.data.frame()
  names(data_sum)[names(data_sum) == "mean"] <- varname
  data_sum
}

#' col2hex
#' Convert R colour names or hex strings to `#RRGGBB`, dropping any alpha.
#' @param x Character vector of colours
#' @return Character vector of `#RRGGBB` strings
#' @noRd
col2hex <- function(x) {
  grDevices::rgb(t(grDevices::col2rgb(x)), maxColorValue = 255)
}
