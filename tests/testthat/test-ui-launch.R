# -----------------------------------------------------------------------------
# Files handed to exploreDE(c(...)) (helper-app.R), which sets the option
# exploreDE.file to their paths: every test here opens a session on one app
# process launched with the same files. In order: a file that is not an
# .h5ad, the demo, a copy of the demo, two RNA-seq .h5ad and the proteomics
# example. The first that loads (the demo) opens; the other RNA-seq files'
# contrasts join Multi DE, the proteomics one is refused.
# -----------------------------------------------------------------------------

launch_files <- local({
  files <- NULL
  function() {
    if (is.null(files)) {
      dir <- file.path(tempdir(), "ui-launch")
      dir.create(dir, showWarnings = FALSE)
      bad <- file.path(dir, "bad.h5ad")
      writeLines("not an h5ad file", bad)
      examples <- system.file("examples", package = "exploreDE")
      demo <- file.path(examples, "genomics", "deResult.rds")
      copy <- file.path(dir, "copy.rds")
      file.copy(demo, copy, overwrite = TRUE)
      files <<- c(
        bad,
        demo,
        copy,
        toy_h5ad(dir, "design1", seed = 1),
        toy_h5ad(dir, "design2", seed = 2),
        file.path(examples, "proteomics", "SummarizedExperiment.rds")
      )
    }
    files
  }
})

launch_session <- function(name, .local_envir = parent.frame()) {
  skip_if_not_installed("anndataR")
  app <- app_session(name, options = list(exploreDE.file = launch_files()), .local_envir = .local_envir)
  app$wait_for_value(output = "settingsTable", timeout = 180e3)
  app$wait_for_idle()
  app
}

test_that("the first launch file that loads opens, and a file that fails is named", {
  app <- launch_session("launch-first-bad")
  expect_match(app$get_js("document.querySelector('.modal-body').innerText"), "bad.h5ad", fixed = TRUE)
  expect_true(all(c("volcanoTab", "oraTab", "multiDataSelectionTab") %in% app_menu_tabs(app)))
  expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})

test_that("the other RNA-seq launch files join Multi DE, the proteomics one does not", {
  app <- launch_session("launch-multi-de")
  app_open_tab(app, "multiDataSelectionTab")
  own <- c("control_vs_treated", "treated_vs_control")
  expect_setequal(
    multi_datasets(app),
    c("diff_vs_undiff", "copy_diff_vs_undiff", paste0("design1_", own), paste0("design2_", own))
  )
  expect_identical(app_output_errors(app), character(0))
})

# Bubble plots compare the pathway results of at least two datasets: the demo
# and its copy carry ORA and GSEA.
test_that("the Multi DE Bubble table shows ORA and GSEA results", {
  app <- launch_session("launch-bubble")
  expect_true("multiBubbleTab" %in% app_menu_tabs(app))
  # The pathway lists are built on the Upset/Venn tab
  app_open_tab(app, "multiUpsetPlotTab")
  app_open_tab(app, "multiBubbleTab")
  headers <- function() {
    unlist(app$get_js("Array.from(document.querySelectorAll('#BubbleTable thead th')).map(th => th.innerText)"))
  }
  rows <- function() app$get_js("document.querySelectorAll('#BubbleTable tbody tr').length")

  app$set_inputs(dataToBubble = "ORA")
  app$set_inputs(bubblePathwayNumber = 5)
  app$wait_for_idle(duration = 1500)
  expect_gt(rows(), 0)
  expect_true(all(c("pvalue", "p.adjust") %in% headers()))

  app$set_inputs(dataToBubble = "GSEA")
  app$wait_for_idle(duration = 1500)
  expect_gt(rows(), 0)
  expect_true(all(c("enrichmentScore", "NES") %in% headers()))
  expect_identical(app_output_errors(app), character(0))
})
