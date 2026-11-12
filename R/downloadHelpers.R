#' Standardized Excel Export
#'
#' Exports a data frame to Excel using openxlsx with consistent formatting:
#' - TableStyleMedium9
#' - Frozen first row
#' - Auto-filter on first row
#' - Banded rows
#'
#' @param df Data frame to export
#' @param filename Output filename (including .xlsx extension)
#' @param sheet_name Name of the worksheet (default: "Data")
#' @param table_name Name of the Excel table (default: "DataTable")
#'
#' @noRd
export_excel_standard <- function(df, filename, sheet_name = "Data", table_name = "DataTable") {
  # Create workbook and add sheet
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet_name)

  # Write data as a formatted table
  openxlsx::writeDataTable(
    wb,
    sheet = sheet_name,
    x = df,
    startCol = 1,
    startRow = 1,
    tableStyle = "TableStyleMedium9",
    tableName = table_name,
    withFilter = TRUE,
    bandedRows = TRUE,
    bandedCols = FALSE
  )

  # Freeze the first row
  openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)

  # Save workbook
  openxlsx::saveWorkbook(wb, filename, overwrite = TRUE)
}

#' Standardized Plot Export
#'
#' Exports a plot object to PDF, SVG, or PNG with consistent settings.
#' Handles ggplot objects, grid objects, and base graphics (passed as function).
#'
#' @param plot_object The plot to save. Can be a ggplot object, a grob, or a function that draws a base plot.
#' @param filename Output filename (including extension)
#' @param format One of "PDF", "SVG", "PNG" (case insensitive)
#' @param dpi Resolution for PNG (default: 300)
#' @param width Width in inches
#' @param height Height in inches
#'
#' @noRd
export_plot_standard <- function(plot_object, filename, format = "pdf", dpi = 300, width = 8, height = 6) {
  format <- tolower(format)

  # Set up the graphics device. PDF stays on base pdf() (Affinity-compatible);
  # Greek in rich level labels is rewritten to the Symbol font for it.
  if (format == "pdf") {
    plot_object <- prepare_plot_for_pdf(plot_object)
    grDevices::pdf(file = filename, width = width, height = height)
  } else if (format == "svg") {
    grDevices::svg(file = filename, width = width, height = height)
  } else if (format == "png") {
    grDevices::png(filename = filename, width = width, height = height, units = "in", res = as.numeric(dpi))
  } else {
    warning("Unknown format '", format, "'. Defaulting to PDF.")
    grDevices::pdf(file = filename, width = width, height = height)
  }

  # Draw the plot. show() draws S4 plots (ComplexHeatmap, UpSet) and calls
  # print() for S3 ones (ggplot, patchwork, GGally::ggpairs).
  tryCatch(
    {
      if (is.function(plot_object)) {
        plot_object()
      } else {
        methods::show(plot_object)
      }
    },
    finally = {
      # Ensure device is closed even if printing fails
      grDevices::dev.off()
    }
  )
}

#' Generate Consistent Filename
#'
#' Creates a standardized filename in the format: `[design]_[base]_[suffix].[ext]`
#'
#' @param base_name Core name of the file (e.g., "volcano", "heatmap")
#' @param suffix Optional suffix (e.g., "counts", "stats")
#' @param extension File extension (e.g., "xlsx", "pdf")
#' @param design_name Optional design/contrast name to prepend
#'
#' @return A string with the constructed filename
#' @noRd
generate_filename <- function(base_name, suffix = NULL, extension, design_name = NULL) {
  # Clean inputs (remove NULLs and empty strings)
  parts <- c(design_name, base_name, suffix)
  parts <- parts[!is.null(parts) & parts != ""]

  # Combine with underscores
  filename <- paste(parts, collapse = "_")

  # Add extension if not already present (simple check)
  if (!grepl(paste0("\\.", extension, "$"), filename)) {
    filename <- paste0(filename, ".", extension)
  }

  return(filename)
}

#' Valid Excel worksheet names
#'
#' Excel refuses sheet names longer than 31 characters, containing
#' `: \ / ? * [ ]` or starting or ending with `'`, the name `History`, and two
#' sheets whose names differ only in case. Replaces those characters with `_`,
#' strips the apostrophes, truncates to 31 characters and makes the names
#' unique (case-insensitively, `History` counting as taken), keeping their
#' order.
#'
#' @param x Character vector of wanted sheet names
#' @return Character vector of the same length
#' @noRd
xlsx_sheet_names <- function(x) {
  x <- gsub("[\\[\\]:\\\\/?*]", "_", x, perl = TRUE)
  x <- gsub("^'+|'+$", "", x)
  x[!nzchar(x)] <- "Sheet"
  out <- substr(x, 1, 31)
  for (i in seq_along(out)) {
    k <- 1L
    while (tolower(out[i]) %in% tolower(c("History", out[seq_len(i - 1L)]))) {
      suffix <- paste0("_", k)
      out[i] <- paste0(substr(x[i], 1, 31 - nchar(suffix)), suffix)
      k <- k + 1L
    }
  }
  out
}
