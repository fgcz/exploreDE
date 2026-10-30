# =============================================================================
# Loading an .h5ad from a URL (?data=https://...)
# Used by inst/app/server-initInputData.R. Only .h5ad files are fetched (HDF5
# holds data, never code), only from the URL prefixes in url_load_sources(),
# without following redirects, and up to a size limit. The file then goes
# through the same checks as an upload.
# =============================================================================

#' URL prefixes exploreDE may load `.h5ad` files from
#'
#' The default list covers APB Export job files and public pStore. Set the
#' environment variable `EXPLOREDE_URL_SOURCES` (comma-separated `https://`
#' prefixes, each ending in `/`) to replace it, e.g. in the ShinyProxy app's
#' `container-env`.
#'
#' @return Character vector of URL prefixes.
#' @export
url_load_sources <- function() {
  env <- Sys.getenv("EXPLOREDE_URL_SOURCES", "")
  if (nzchar(env)) {
    return(trimws(strsplit(env, ",", fixed = TRUE)[[1]]))
  }
  c("https://apps-dev.bfabric.org/apb-export/api/jobs/", "https://fgcz-ms.uzh.ch/public/pStore/")
}

#' May exploreDE load this URL?
#'
#' Allowed: An `https://` URL with no user name or password, whose path has no
#' `..` (also not percent-encoded), ends in `.h5ad`, and starts with one of
#' `sources`. The comparison is on the URL text, so a prefix must end in `/`
#' to stop `.../apb-export-other/` matching `.../apb-export/`.
#'
#' @param url Character(1).
#' @param sources URL prefixes, by default [url_load_sources()].
#' @return `TRUE`, or `FALSE` with attribute `reason`.
#' @export
url_load_allowed <- function(url, sources = url_load_sources()) {
  no <- function(reason) structure(FALSE, reason = reason)
  if (!is.character(url) || length(url) != 1L || is.na(url)) {
    return(no("not a URL"))
  }
  m <- regmatches(url, regexec("^https://([A-Za-z0-9.-]+)(:[0-9]+)?(/[^?#]*)$", url))[[1]]
  if (!length(m)) {
    return(no("only plain https:// links to a file are accepted"))
  }
  path <- m[[4]]
  if (grepl("(^|/)\\.\\.?(/|$)|%2e|%2f|%5c|\\\\", path, ignore.case = TRUE)) {
    return(no("the link contains a relative path"))
  }
  if (!grepl("\\.h5ad$", path, ignore.case = TRUE)) {
    return(no("only .h5ad files can be loaded from a link"))
  }
  sources <- sources[startsWith(sources, "https://") & endsWith(sources, "/")]
  if (!any(startsWith(url, sources))) {
    return(no("this site is not on exploreDE's list of data sources"))
  }
  TRUE
}

#' Download an allowed `.h5ad` URL to a file
#'
#' No redirects are followed, the download stops at `max_bytes`, and anything
#' but HTTP 200 is an error.
#'
#' @param url An URL for which [url_load_allowed()] is `TRUE`.
#' @param dest Destination path.
#' @param max_bytes Size limit (1 GB, as for uploads).
#' @param timeout Seconds.
#' @return `dest`, invisibly. Stops with a short message on failure.
#' @export
url_load_fetch <- function(url, dest, max_bytes = 1000 * 1024^2, timeout = 600) {
  ok <- url_load_allowed(url)
  if (!isTRUE(ok)) {
    stop(attr(ok, "reason"))
  }
  h <- curl::new_handle(
    followlocation = FALSE,
    timeout = timeout,
    maxfilesize_large = max_bytes,
    protocols_str = "https"
  )
  res <- tryCatch(curl::curl_fetch_disk(url, dest, handle = h), error = function(e) {
    stop("download failed: ", conditionMessage(e))
  })
  if (res$status_code != 200L) {
    unlink(dest)
    stop(sprintf("the server answered %d", res$status_code))
  }
  if (file.size(dest) > max_bytes) {
    unlink(dest)
    stop("the file is larger than the limit")
  }
  invisible(dest)
}
