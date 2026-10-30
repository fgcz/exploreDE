# ---------------------------------------------------------------------------
# Banana mode: hidden easter egg that swaps p-value significance markers for a
# pile of banana images. The marker PNG is base64-EMBEDDED below (not a served
# www/ asset) so nothing is fetchable from the running app; it is decoded to
# the session tempdir on first use. Rendering is raster (grid.raster for
# ComplexHeatmap cells, ggtext <img> for ggplot brackets) so it draws the same
# on every graphics device.
#
# Activated by app.R: triple-click the footer heart (#bananaHeart) or ?banana=1.
# ---------------------------------------------------------------------------

# Encoded marker image (PNG). Opaque data.
.banana_b64 <- "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
"

.banana_env <- new.env(parent = emptyenv())

# Banana mode is a hidden easter egg — it must NEVER take down the app. Every
# function below is failure-tolerant: on any error it returns the benign "no
# banana" value (NULL / "" / invisible(NULL)), and `banana_on()` gates the
# whole feature on the optional deps actually being present. The plot-builder
# call sites run inside observeEvent()s where an uncaught error would kill the
# session, so this is load-bearing, not just tidiness.
.banana_safe <- function(expr, fallback) tryCatch(expr, error = function(e) fallback)

#' Absolute path to the decoded marker PNG (written to the session tempdir on
#' first call, cached thereafter). Returns "" if it can't be written.
#' @return length-1 character path ("" on failure)
#' @export
banana_png_path <- function() {
  .banana_safe(
    {
      p <- .banana_env$path
      if (is.null(p) || !file.exists(p)) {
        p <- file.path(tempdir(), "exploreDE-mk.png")
        writeBin(jsonlite::base64_dec(.banana_b64), p)
        .banana_env$path <- p
      }
      p
    },
    fallback = ""
  )
}

#' Significance tier (0-4) for a p-value / FDR. Shared cutpoints with the
#' asterisk ladder: <=1e-4 -> 4, <=1e-3 -> 3, <=1e-2 -> 2, <=0.05 -> 1, else 0.
#' @return integer scalar
#' @export
banana_tier <- function(p) {
  if (length(p) == 0L || is.na(p)) {
    return(0L)
  }
  if (p <= 1e-04) {
    4L
  } else if (p <= 1e-03) {
    3L
  } else if (p <= 1e-02) {
    2L
  } else if (p <= 5e-02) {
    1L
  } else {
    0L
  }
}

#' Is banana mode active AND safe to render? Requires the toggle on, the
#' embedded asset decodable, and the optional deps (ggtext, png) installed.
#' If any of that is missing, banana mode is simply off everywhere and no
#' banana code path is entered.
#' @param toggle The raw bananaMode() reactive value.
#' @return logical scalar
#' @export
banana_on <- function(toggle) {
  isTRUE(toggle) &&
    nzchar(.banana_b64) &&
    .banana_safe(requireNamespace("ggtext", quietly = TRUE), FALSE) &&
    .banana_safe(requireNamespace("png", quietly = TRUE), FALSE) &&
    nzchar(banana_png_path())
}

#' Decoded marker raster (RGBA array), cached. NULL on any failure.
#' @export
banana_raster <- function() {
  .banana_safe(
    {
      src <- banana_png_path()
      if (!nzchar(src)) {
        return(NULL)
      }
      if (!identical(.banana_env$img_src, src)) {
        .banana_env$img <- png::readPNG(src)
        .banana_env$img_src <- src
      }
      .banana_env$img
    },
    fallback = NULL
  )
}

#' <img> markdown for n markers -- for ggtext::geom_richtext labels.
#' @param n integer count; px display height in pixels.
#' @return length-1 character ("" when n < 1 or the asset is unavailable)
#' @export
banana_md <- function(n, px = 16L) {
  n <- as.integer(n)
  if (is.na(n) || n < 1L) {
    return("")
  }
  src <- banana_png_path()
  if (!nzchar(src)) {
    return("")
  }
  q <- "'"
  strrep(paste0("<img src=", q, src, q, " height=", q, max(6L, as.integer(px)), q, "/>"), n)
}

#' Draw n markers in a centred horizontal row at (x, y) -- for a ComplexHeatmap
#' cell_fun (call inside the cell viewport). Never throws.
#' @param n integer count; h glyph size (a grid unit).
#' @export
banana_grid <- function(n, x, y, h = grid::unit(4.2, "mm")) {
  .banana_safe(
    {
      n <- as.integer(n)
      img <- banana_raster()
      if (is.na(n) || n < 1L || is.null(img)) {
        return(invisible(NULL))
      }
      gap <- grid::unit(0.5, "mm")
      for (k in seq_len(n)) {
        off <- k - (n + 1) / 2
        grid::grid.raster(img, x = x + off * (h + gap), y = y, width = h, height = h, interpolate = TRUE)
      }
      invisible(NULL)
    },
    fallback = invisible(NULL)
  )
}

#' A ggtext::geom_richtext layer: a row of marker images at (x, y) per facet,
#' sitting clear above the (kept, label-blanked) p-value bracket. Returns NULL
#' when nothing qualifies -- p + NULL is safe.
#' @param data plot data (one+ rows per facet level).
#' @param x numeric x position for the row (bracket midpoint).
#' @param y numeric -- scalar or a vector aligned to nrow(data) -- bracket-line y.
#' @param fdr_col,facet_col column names in data.
#' @param px marker display height in px.
#' @export
banana_layer <- function(data, x, y, fdr_col, facet_col = "Feature", px = 16L) {
  .banana_safe(
    {
      if (is.null(data) || !all(c(fdr_col, facet_col) %in% names(data))) {
        return(NULL)
      }
      keep <- !duplicated(data[[facet_col]])
      d <- data[keep, c(facet_col, fdr_col), drop = FALSE]
      names(d) <- c(facet_col, "fdr")
      d$by <- if (length(y) == 1L) rep(y, nrow(d)) else y[keep]
      d$n <- vapply(d$fdr, banana_tier, integer(1))
      d <- d[d$n >= 1L & is.finite(d$by), , drop = FALSE]
      if (!nrow(d)) {
        return(NULL)
      }
      d$md <- vapply(d$n, banana_md, character(1), px = px)
      if (!any(nzchar(d$md))) {
        return(NULL)
      }
      d$bx <- x
      ggtext::geom_richtext(
        data = d,
        mapping = ggplot2::aes(x = bx, y = by, label = md),
        vjust = -0.1,
        inherit.aes = FALSE,
        fill = NA,
        label.colour = NA,
        label.padding = grid::unit(c(0, 0, 0, 0), "pt")
      )
    },
    fallback = NULL
  )
}
