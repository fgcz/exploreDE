# -----------------------------------------------------------------------------
# ?data=https://... loading: url_load_allowed(), url_load_sources(),
# url_load_fetch() (R/utils_url_load.R)
# -----------------------------------------------------------------------------

apb <- "https://apps-dev.bfabric.org/apb-export/api/jobs/1234/files/maxquant_protein.h5ad"

test_that("links to .h5ad files under an allowed prefix are accepted", {
  expect_true(url_load_allowed(apb))
  expect_true(url_load_allowed("https://fgcz-ms.uzh.ch/public/pStore/p1000/bfabric/Proteomics/x/AnnData.h5ad"))
})

test_that("everything else is refused, with a reason", {
  bad <- c(
    http = "http://apps-dev.bfabric.org/apb-export/api/jobs/1/files/a.h5ad",
    userinfo = "https://apps-dev.bfabric.org@evil.example/apb-export/api/jobs/1/files/a.h5ad",
    query = "https://apps-dev.bfabric.org/apb-export/api/jobs/1/files/a.h5ad?x=1",
    dotdot = "https://apps-dev.bfabric.org/apb-export/api/jobs/../../admin/a.h5ad",
    dotdot_enc = "https://apps-dev.bfabric.org/apb-export/api/jobs/%2e%2e/%2E%2E/a.h5ad",
    slash_enc = "https://apps-dev.bfabric.org/apb-export/api/jobs/1%2F..%2Fa.h5ad",
    other_host = "https://evil.example/apb-export/api/jobs/1/files/a.h5ad",
    prefix_sib = "https://apps-dev.bfabric.org/apb-export-other/api/jobs/1/files/a.h5ad",
    not_h5ad = "https://apps-dev.bfabric.org/apb-export/api/jobs/1/files/a.rds",
    file = "file:///etc/passwd.h5ad"
  )
  for (nm in names(bad)) {
    res <- url_load_allowed(bad[[nm]])
    expect_false(res, info = nm)
    expect_true(nzchar(attr(res, "reason")), info = nm)
  }
  for (x in list(NULL, NA_character_, c(apb, apb), 1)) {
    expect_false(url_load_allowed(x))
  }
})

test_that("EXPLOREDE_URL_SOURCES replaces the list; prefixes must be https and end in /", {
  withr::local_envvar(
    EXPLOREDE_URL_SOURCES = "https://data.example.org/share/, http://insecure.example/, https://noslash.example"
  )
  expect_true(url_load_allowed("https://data.example.org/share/a.h5ad"))
  expect_false(url_load_allowed(apb))
  expect_false(url_load_allowed("https://noslash.example/a.h5ad"))
  expect_false(url_load_allowed("https://noslash.example.evil/a.h5ad"))
})

test_that("url_load_fetch refuses a link before touching the network", {
  dest <- withr::local_tempfile(fileext = ".h5ad")
  expect_error(url_load_fetch("https://evil.example/a.h5ad", dest), "list of data sources")
  expect_false(file.exists(dest))
})

test_that("url_load_fetch downloads an allowed file, and reports HTTP errors and the size limit", {
  skip_on_cran()
  skip_if_offline("raw.githubusercontent.com")
  withr::local_envvar(EXPLOREDE_URL_SOURCES = "https://raw.githubusercontent.com/scverse/anndataR/")
  base <- "https://raw.githubusercontent.com/scverse/anndataR/main/inst/extdata/"
  dest <- withr::local_tempfile(fileext = ".h5ad")
  url_load_fetch(paste0(base, "example.h5ad"), dest)
  expect_gt(file.size(dest), 1000)
  expect_s3_class(anndataR::read_h5ad(dest), "AbstractAnnData")
  expect_error(url_load_fetch(paste0(base, "does_not_exist.h5ad"), withr::local_tempfile()), "answered 404")
  expect_error(url_load_fetch(paste0(base, "example.h5ad"), withr::local_tempfile(), max_bytes = 1000))
})
