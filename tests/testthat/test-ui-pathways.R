# -----------------------------------------------------------------------------
# One session on an uploaded .h5ad with ORA results on one contrast and GSEA
# results on the other (helper-app.R). The ORA tab shows the results of a
# contrast other than the primary one and clears them when switching back;
# the GSEA table comes back after a contrast without results.
# -----------------------------------------------------------------------------

# A clusterProfiler ORA result with the columns the ORA table shows.
ora_result <- function() {
  ids <- c("GO:0001", "GO:0002")
  res <- data.frame(
    ID = ids,
    Description = c("term 1", "term 2"),
    GeneRatio = "2/10",
    BgRatio = "5/80",
    pvalue = 0.001,
    p.adjust = 0.01,
    qvalue = 0.01,
    geneID = "Gene001/Gene002",
    Count = 2L,
    row.names = ids
  )
  methods::new("enrichResult", result = res)
}

# The ORA table, or the message shown in its place (an htmlwidget shows it
# next to the table, in the same box).
ora_table_text <- function(app) {
  app$get_js("document.getElementById('selectedTable_ORA').parentElement.innerText")
}

gsea_pane_text <- function(app) {
  app$get_js("document.getElementById('shiny-tab-gseaTab').innerText")
}

skip_if_not_installed("anndataR")
skip_if_not_installed("DOSE")
dds <- make_toy_dds(n_genes = 80)
h5ad <- withr::local_tempfile(fileext = ".h5ad")
# The primary contrast (first as written and alphabetically) has GSEA and no
# ORA; the other has ORA and no GSEA.
build_explore_h5ad(
  dds,
  out = h5ad,
  omics = "rnaseq",
  de = list(
    control_vs_treated = DESeq2::results(dds, contrast = c("condition", "control", "treated")),
    treated_vs_control = DESeq2::results(dds, contrast = c("condition", "treated", "control"))
  ),
  pathways = list(
    control_vs_treated = list(gsea = list(BP = .mini_cp("gseaResult", 3))),
    treated_vs_control = list(ora = list(BP = list(upGenes = ora_result())))
  ),
  organism = "Homo sapiens"
)
app <- app_session("pathways", .local_envir = environment())
app$upload_file(uploadH5ad = h5ad)
app$wait_for_value(output = "settingsTable", timeout = 180e3)
app$wait_for_idle()

test_that("the ORA tab shows the results of a contrast other than the primary one", {
  expect_true("oraTab" %in% app_menu_tabs(app))

  app_open_tab(app, "oraTab")
  expect_match(
    app$get_js("document.querySelector('.modal-body').innerText"),
    "no ORA results for the selected contrast",
    fixed = TRUE
  )
  app$run_js("$('.modal').modal('hide')")

  app$set_inputs(contrastSelected_ora = "treated_vs_control")
  app$wait_for_idle()
  expect_match(ora_table_text(app), "term 1", fixed = TRUE)
  expect_match(ora_table_text(app), "term 2", fixed = TRUE)
  expect_identical(app_output_errors(app), character(0))

  app$set_inputs(contrastSelected_ora = "control_vs_treated")
  app$wait_for_idle()
  app$run_js("$('.modal').modal('hide')")
  expect_no_match(ora_table_text(app), "term 1", fixed = TRUE)
  expect_match(ora_table_text(app), "No ORA results for this contrast", fixed = TRUE)
  expect_identical(app_output_errors(app), character(0))
})

test_that("the GSEA table comes back after a contrast without results", {
  app_open_tab(app, "gseaTab")
  app$set_inputs(contrastSelected_gsea = "control_vs_treated")
  app$wait_for_idle()
  expect_match(gsea_pane_text(app), "term 1", fixed = TRUE)

  app$set_inputs(contrastSelected_gsea = "treated_vs_control")
  app$wait_for_idle()
  app$run_js("$('.modal').modal('hide')")
  expect_match(gsea_pane_text(app), "No GSEA results for this contrast", fixed = TRUE)

  app$set_inputs(contrastSelected_gsea = "control_vs_treated")
  app$wait_for_idle()
  expect_match(gsea_pane_text(app), "term 1", fixed = TRUE)
  expect_identical(app_output_errors(app), character(0))
})
