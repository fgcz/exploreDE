
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
#' @export
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
#' @export
export_plot_standard <- function(plot_object, filename, format = "pdf", 
                                 dpi = 300, width = 8, height = 6) {
  format <- tolower(format)
  
  # Set up the graphics device. PDF stays on base pdf() (Affinity-compatible);
  # Greek in rich level labels is rewritten to the Symbol font for it.
  if (format == "pdf") {
    plot_object <- prepare_plot_for_pdf(plot_object)
    pdf(file = filename, width = width, height = height)
  } else if (format == "svg") {
    svg(file = filename, width = width, height = height)
  } else if (format == "png") {
    png(filename = filename, width = width, height = height, 
        units = "in", res = as.numeric(dpi))
  } else {
    warning("Unknown format '", format, "'. Defaulting to PDF.")
    pdf(file = filename, width = width, height = height)
  }
  
  # Print/Draw the plot
  tryCatch({
    if (inherits(plot_object, "ggplot") || inherits(plot_object, "grob")) {
      print(plot_object)
    } else if (is.function(plot_object)) {
      plot_object()
    } else {
      # Fallback for other objects (e.g. ComplexHeatmap) that implement draw() or print()
      print(plot_object)
    }
  }, finally = {
    # Ensure device is closed even if printing fails
    dev.off()
  })
}

#' Generate Consistent Filename
#'
#' Creates a standardized filename in the format: [design]_[base]_[suffix].[ext]
#'
#' @param base_name Core name of the file (e.g., "volcano", "heatmap")
#' @param suffix Optional suffix (e.g., "counts", "stats")
#' @param extension File extension (e.g., "xlsx", "pdf")
#' @param design_name Optional design/contrast name to prepend
#'
#' @return A string with the constructed filename
#' @export
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
