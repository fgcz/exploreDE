# Writes the two prolfquapp AnnData fixtures used by test-prolfquapp_h5ad.R
# (simulated data only). Run from the package root.
#
# - prolfquapp_2.10_AnnData.h5ad: the first 40 proteins of a DEA AnnData.h5ad
#   written by prolfquapp 2.10 (layout schema 2.0.0, read by exploreDE), from
#   a simulated dataset in the ptmbenchmark workspace (a local path).
# - prolfquapp_2.11_AnnData.h5ad: a DEA of prolfquapp's 30-protein simulated
#   example, written by prolfquapp 2.11.1 (layout schema 2.1.0, not read yet).

trim_prolfquapp_h5ad <- function(src, out, n = 40L) {
  ad <- exploreDE:::.h5ad_read(src)
  nv <- nrow(ad$var)
  keep <- seq_len(min(n, nv))
  # Every per-feature array in uns (prolfquapp keeps the non-numeric varm
  # columns there) is cut to the same features.
  cut <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
      return(lapply(x, cut))
    }
    if (is.atomic(x) && length(x) == nv) {
      return(as.vector(x)[keep])
    }
    x
  }
  trimmed <- anndataR::AnnData(
    X = if (is.null(ad$X)) NULL else ad$X[, keep, drop = FALSE],
    obs = ad$obs,
    var = ad$var[keep, , drop = FALSE],
    layers = lapply(ad$layers, function(m) m[, keep, drop = FALSE]),
    varm = lapply(ad$varm, function(m) m[keep, , drop = FALSE]),
    uns = cut(ad$uns)
  )
  anndataR::write_h5ad(trimmed, out, mode = "w")
  invisible(out)
}

sim <- "/Users/wolski/projects/prolfqua_fml/ptmbenchmark/work"
trim_prolfquapp_h5ad(
  file.path(
    sim,
    "simulations/sim1_sd02_rep02_cond2/pipeline/protein",
    "protein_20260921_PIptmbenchmark_Osim1_sd02_rep02_cond2_WUprotein_none/Results_WU_protein/AnnData.h5ad"
  ),
  "tests/testthat/fixtures/prolfquapp_2.10_AnnData.h5ad"
)
stopifnot(utils::packageVersion("prolfquapp") == "2.11.1")
dea <- prolfquapp::example_deanalyse(Nprot = 30)
se <- prolfquapp::DEAReportGenerator$new(dea, dea$prolfq_app_config, name = "")$make_SummarizedExperiment()
prolfquapp:::write_summarized_experiment_h5ad(se, "tests/testthat/fixtures/prolfquapp_2.11_AnnData.h5ad")
