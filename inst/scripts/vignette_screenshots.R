# Screenshots of the app on the demo dataset for vignettes/exploreDE.Rmd.
#
# Run from the package root: Rscript inst/scripts/vignette_screenshots.R
# Needs shinytest2 and Chrome. Writes vignettes/figures/<tab>.png.

# AppDriver skips itself unless told it is not on CRAN
Sys.setenv(NOT_CRAN = "true")
out_dir <- file.path("vignettes", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

app <- shinytest2::AppDriver$new(
  "inst/app",
  name = "vignette",
  width = 1400,
  height = 1050,
  load_timeout = 180e3,
  timeout = 60e3
)
on.exit(app$stop(), add = TRUE)

app$click("loadDemo")
invisible(app$wait_for_value(output = "settingsTable", timeout = 180e3))
app$wait_for_idle()

shoot <- function(tab, file = tab) {
  app$set_inputs(tabs = tab)
  app$wait_for_idle(duration = 2000, timeout = 120e3)
  png_file <- file.path(out_dir, paste0(file, ".png"))
  unlink(png_file) # get_screenshot() does not overwrite
  # Scale 0.75: 1050 px wide, still above the vignette's display width
  app$get_screenshot(
    png_file,
    selector = "viewport",
    screenshot_args = list(scale = 0.75)
  )
}

shoot("ordinationTab", "ordination")
shoot("heatmapTab", "heatmap")
shoot("volcanoTab", "volcano")
