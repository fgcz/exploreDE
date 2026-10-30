# ---------------------------------------------------------------------------
# Custom factor-level display labels (Colours tab).
#
# Users can give any factor level a display label with a small HTML subset
# (<i>, <b>, <em>, <strong>, <sup>, <sub>, <br>) plus HTML entities
# (e.g. "&Delta;"). Labels are DISPLAY-ONLY: the raw level stays the key for
# input IDs, buckets, p-value brackets, contrasts and save-state.
#
# Three renderings of every label:
#   * rich  -- ggtext markdown for ggplot (element_markdown). User text is
#              backslash-escaped so only the whitelisted tags are markup.
#   * plain -- tags stripped, <br> -> "\n". Fallback when ggtext/gridtext
#              aren't installed.
#   * ComplexHeatmap legends are rich too, but built by us at DRAW time
#              (build_rich_heatmap_legend(), called from draw_heatmap()):
#              ComplexHeatmap 2.28 sizes gt_render legend labels by their raw
#              markup length (discrete_legend_body drops the gridtext class),
#              which blows the legend out to page width, so the width is
#              re-measured with gridtext on the device actually drawing.
#   * PDF   -- base pdf() is single-byte (Latin-1), so Greek letters are
#              rewritten at export time into spans of the standard PDF Symbol
#              font (prepare_plot_for_pdf()). Keeps the file Affinity-safe,
#              which cairo_pdf is not.
# ---------------------------------------------------------------------------

.label_allowed_tags <- c("i", "b", "em", "strong", "sup", "sub", "br")
.label_tag_pattern <- "<\\s*/?\\s*([A-Za-z]+)\\s*/?\\s*>"

# Named HTML entities we decode (everything else numeric: &#916; / &#x394;).
.label_entities <- c(
  amp = "&",
  lt = "<",
  gt = ">",
  quot = "\"",
  apos = "'",
  nbsp = "\u00a0",
  micro = "\u00b5",
  deg = "\u00b0",
  plusmn = "\u00b1",
  times = "\u00d7",
  middot = "\u00b7",
  minus = "\u2212",
  ndash = "\u2013",
  mdash = "\u2014",
  hellip = "\u2026",
  rarr = "\u2192",
  larr = "\u2190",
  Alpha = "\u0391",
  Beta = "\u0392",
  Gamma = "\u0393",
  Delta = "\u0394",
  Epsilon = "\u0395",
  Zeta = "\u0396",
  Eta = "\u0397",
  Theta = "\u0398",
  Iota = "\u0399",
  Kappa = "\u039a",
  Lambda = "\u039b",
  Mu = "\u039c",
  Nu = "\u039d",
  Xi = "\u039e",
  Omicron = "\u039f",
  Pi = "\u03a0",
  Rho = "\u03a1",
  Sigma = "\u03a3",
  Tau = "\u03a4",
  Upsilon = "\u03a5",
  Phi = "\u03a6",
  Chi = "\u03a7",
  Psi = "\u03a8",
  Omega = "\u03a9",
  alpha = "\u03b1",
  beta = "\u03b2",
  gamma = "\u03b3",
  delta = "\u03b4",
  epsilon = "\u03b5",
  zeta = "\u03b6",
  eta = "\u03b7",
  theta = "\u03b8",
  iota = "\u03b9",
  kappa = "\u03ba",
  lambda = "\u03bb",
  mu = "\u03bc",
  nu = "\u03bd",
  xi = "\u03be",
  omicron = "\u03bf",
  pi = "\u03c0",
  rho = "\u03c1",
  sigmaf = "\u03c2",
  sigma = "\u03c3",
  tau = "\u03c4",
  upsilon = "\u03c5",
  phi = "\u03c6",
  chi = "\u03c7",
  psi = "\u03c8",
  omega = "\u03c9"
)

# Greek (Unicode) -> character in the Adobe Symbol font encoding.
.symbol_font_map <- c(
  "\u0391" = "A",
  "\u0392" = "B",
  "\u0393" = "G",
  "\u0394" = "D",
  "\u2206" = "D",
  "\u0395" = "E",
  "\u0396" = "Z",
  "\u0397" = "H",
  "\u0398" = "Q",
  "\u0399" = "I",
  "\u039a" = "K",
  "\u039b" = "L",
  "\u039c" = "M",
  "\u039d" = "N",
  "\u039e" = "X",
  "\u039f" = "O",
  "\u03a0" = "P",
  "\u03a1" = "R",
  "\u03a3" = "S",
  "\u03a4" = "T",
  "\u03a5" = "U",
  "\u03a6" = "F",
  "\u03a7" = "C",
  "\u03a8" = "Y",
  "\u03a9" = "W",
  "\u2126" = "W",
  "\u03b1" = "a",
  "\u03b2" = "b",
  "\u03b3" = "g",
  "\u03b4" = "d",
  "\u03b5" = "e",
  "\u03b6" = "z",
  "\u03b7" = "h",
  "\u03b8" = "q",
  "\u03b9" = "i",
  "\u03ba" = "k",
  "\u03bb" = "l",
  "\u03bc" = "m",
  "\u03bd" = "n",
  "\u03be" = "x",
  "\u03bf" = "o",
  "\u03c0" = "p",
  "\u03c1" = "r",
  "\u03c2" = "V",
  "\u03c3" = "s",
  "\u03c4" = "t",
  "\u03c5" = "u",
  "\u03c6" = "f",
  "\u03c7" = "c",
  "\u03c8" = "y",
  "\u03c9" = "w",
  "\u03d1" = "J",
  "\u03d5" = "j",
  "\u03d6" = "v"
)

# Non-Latin-1 punctuation base pdf() can't draw -> ASCII stand-ins (PDF only).
.pdf_ascii_map <- c(
  "\u2212" = "-",
  "\u2013" = "-",
  "\u2014" = "-",
  "\u2018" = "'",
  "\u2019" = "'",
  "\u201c" = "\"",
  "\u201d" = "\"",
  "\u2026" = "...",
  "\u2192" = "->",
  "\u2190" = "<-"
)

.pdf_symbol_family <- "GreekSym"

#' Backslash-escape all ASCII punctuation so text is literal in ggtext markdown
#'
#' @param x Character vector.
#' @return Character vector safe to embed in a ggtext/gridtext string.
#' @export
escape_label_markdown <- function(x) {
  gsub("([[:punct:]])", "\\\\\\1", x, perl = TRUE)
}

#' Decode HTML entities (named subset + numeric) to Unicode
#' @param x Character vector.
#' @export
decode_label_entities <- function(x) {
  vapply(
    x,
    function(s) {
      if (is.na(s) || !grepl("&", s, fixed = TRUE)) {
        return(s)
      }
      m <- gregexpr("&(#[0-9]+|#[xX][0-9A-Fa-f]+|[A-Za-z]+);", s)
      regmatches(s, m) <- list(vapply(
        regmatches(s, m)[[1]],
        function(e) {
          body <- substr(e, 2, nchar(e) - 1)
          if (startsWith(body, "#x") || startsWith(body, "#X")) {
            return(intToUtf8(strtoi(substring(body, 3), 16L)))
          }
          if (startsWith(body, "#")) {
            return(intToUtf8(as.integer(substring(body, 2))))
          }
          if (body %in% names(.label_entities)) {
            return(.label_entities[[body]])
          }
          # Any other named HTML entity (&iuml;, &eacute;, ...): let libxml2 decode
          # it; unknown names stay as typed.
          tryCatch(
            {
              d <- xml2::xml_text(xml2::read_html(paste0("<p>", e, "</p>")))
              if (nzchar(d)) d else e
            },
            error = function(err) e
          )
        },
        character(1)
      ))
      s
    },
    character(1),
    USE.NAMES = FALSE
  )
}

#' Parse one user-typed level label into rich (ggtext) and plain renderings
#'
#' Whitelisted tags are kept, everything else is literal text. A label that
#' gridtext still can't render (e.g. unbalanced tags) falls back to its plain
#' text with `ok = FALSE`, so a typo can never crash a plot.
#'
#' @param x Single string.
#' @return list(rich, plain, ok)
#' @export
parse_level_label <- function(x) {
  x <- enc2utf8(as.character(x))
  m <- gregexpr(.label_tag_pattern, x, perl = TRUE)
  parts <- regmatches(x, m, invert = NA)[[1]]
  # invert = NA interleaves text, tag, text, ...: odd = text, even = tag
  is_tag <- seq_along(parts) %% 2 == 0
  rich <- plain <- character(length(parts))
  for (i in seq_along(parts)) {
    p <- parts[[i]]
    if (is_tag[i]) {
      tag <- tolower(sub(.label_tag_pattern, "\\1", p, perl = TRUE))
      if (tag %in% .label_allowed_tags) {
        closing <- grepl("^<\\s*/", p)
        rich[i] <- if (tag == "br") "<br>" else paste0("<", if (closing) "/", tag, ">")
        plain[i] <- if (tag == "br") "\n" else ""
        next
      }
    }
    txt <- decode_label_entities(p)
    # gridtext splits text on Unicode whitespace, so a whitespace-only run
    # between two tags ("</sup> <i>") vanishes and squashes the words
    # together -- &nbsp; included. Anchoring the space between zero-width
    # spaces keeps it; symbolise_rich_labels() swaps those for a PDF-safe
    # anchor, since U+200B isn't Latin-1.
    if (i > 1 && i < length(parts) && grepl("^\\s+$", txt)) {
      rich[i] <- "\u200b \u200b"
      plain[i] <- " "
      next
    }
    rich[i] <- escape_label_markdown(txt)
    plain[i] <- txt
  }
  rich <- paste(rich, collapse = "")
  plain <- paste(plain, collapse = "")
  ok <- tryCatch(
    {
      gridtext::richtext_grob(rich)
      TRUE
    },
    error = function(e) FALSE
  )
  if (!ok) {
    rich <- gsub("\n", "<br>", escape_label_markdown(plain), fixed = TRUE)
  }
  list(rich = rich, plain = plain, ok = ok)
}

#' Can this session render rich (ggtext) labels?
#' @export
rich_labels_available <- function() {
  requireNamespace("ggtext", quietly = TRUE) && requireNamespace("gridtext", quietly = TRUE)
}

#' Display labels for a factor's levels, read from the Colours tab inputs
#'
#' @param input Shiny input (or list) holding `GroupLabel<factor><level>`.
#' @param fname Factor (column) name.
#' @param levels Raw levels, in plot order.
#' @param mode "rich" for ggplot (ggtext) or "plain" for ComplexHeatmap/base.
#' @param wrap Replace `_ . space /` in UNLABELLED levels with line breaks
#'   (the boxplot "format condition names" option).
#' @param force_rich In rich mode, return escaped markdown even when no level
#'   has a custom label -- for a second scale (e.g. shape) that shares a
#'   legend.text element_markdown() switched on by another factor.
#' @return Character vector named by the raw levels, with attribute
#'   `rich = TRUE` when it must be drawn with element_markdown().
#' @export
level_labels <- function(input, fname, levels, mode = c("rich", "plain"), wrap = FALSE, force_rich = FALSE) {
  mode <- match.arg(mode)
  levels <- as.character(levels)
  custom <- vapply(
    levels,
    function(l) {
      v <- if (is.null(fname) || is.null(input)) NULL else input[[paste0("GroupLabel", fname, l)]]
      if (is.null(v) || is.na(v[1]) || !nzchar(trimws(v[1]))) NA_character_ else trimws(v[1])
    },
    character(1),
    USE.NAMES = FALSE
  )
  has_custom <- !is.na(custom)
  # Wrap BEFORE escaping (via a placeholder), or the escapes on "_" "." "/"
  # would be left dangling in front of the <br>.
  wrap_raw <- function(x, br, escape = FALSE) {
    if (wrap) {
      x <- gsub("\\_|\\.|\\ |\\/", "\001", x)
    }
    if (escape) {
      x <- escape_label_markdown(x)
    }
    gsub("\001", br, x, fixed = TRUE)
  }

  rich <- mode == "rich" && (any(has_custom) || force_rich) && rich_labels_available()
  out <- vapply(
    seq_along(levels),
    function(i) {
      if (!has_custom[i]) {
        if (!rich) {
          return(wrap_raw(levels[i], "\n"))
        }
        return(wrap_raw(levels[i], "<br>", escape = TRUE))
      }
      p <- parse_level_label(custom[i])
      if (rich) p$rich else p$plain
    },
    character(1)
  )
  names(out) <- levels
  attr(out, "rich") <- rich
  out
}

#' Is a label vector (from level_labels()) rich markdown?
#' @export
is_rich_labels <- function(x) isTRUE(attr(x, "rich"))

#' Theme switching the given text elements to ggtext::element_markdown()
#'
#' Add AFTER the complete theme (theme_prism()/theme_bw()): the markdown
#' element merges with, and so inherits, the existing size/face/angle.
#'
#' @export
rich_label_theme <- function(axis_x = FALSE, axis_y = FALSE, legend = FALSE, strip = FALSE) {
  if (!rich_labels_available()) {
    return(ggplot2::theme())
  }
  els <- list()
  if (axis_x) {
    els$axis.text.x <- ggtext::element_markdown()
  }
  if (axis_y) {
    els$axis.text.y <- ggtext::element_markdown()
  }
  if (legend) {
    els$legend.text <- ggtext::element_markdown()
  }
  # strip.text.x/.y explicitly: theme_prism() styles those, which would
  # otherwise shadow a markdown strip.text parent.
  if (strip) {
    els$strip.text.x <- ggtext::element_markdown()
    els$strip.text.y <- ggtext::element_markdown()
  }
  do.call(ggplot2::theme, els)
}

#' Characters in labels that base pdf() cannot draw even after the Greek ->
#' Symbol-font and punctuation rewrites
#'
#' @param x Character vector of labels (plain or rich).
#' @return Unique offending characters.
#' @export
label_pdf_unsupported_chars <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(character(0))
  }
  chars <- unique(unlist(strsplit(enc2utf8(x), "")))
  chars <- setdiff(chars, c(names(.symbol_font_map), names(.pdf_ascii_map)))
  chars[is.na(iconv(chars, "UTF-8", "latin1"))]
}

#' Does any label contain Greek letters?
#' @export
label_has_greek <- function(x) {
  x <- x[!is.na(x)]
  length(x) > 0 && any(strsplit(paste(x, collapse = ""), "")[[1]] %in% names(.symbol_font_map))
}

#' Register the Symbol font as a pdf() family (idempotent)
#' @export
ensure_pdf_symbol_font <- function() {
  if (!.pdf_symbol_family %in% names(grDevices::pdfFonts())) {
    fam <- list(grDevices::Type1Font("Symbol", rep("Symbol.afm", 4), encoding = "AdobeSym.enc"))
    names(fam) <- .pdf_symbol_family
    do.call(grDevices::pdfFonts, fam)
  }
  invisible(.pdf_symbol_family)
}

#' Rewrite rich labels for base pdf(): Greek -> Symbol-font spans,
#' non-Latin-1 punctuation -> ASCII
#' @param x Character vector of rich (markdown) labels; names are kept.
#' @export
symbolise_rich_labels <- function(x) {
  nm <- names(x)
  # Zero-width space anchors (parse_level_label()) -> a near-invisible glyph
  x <- gsub("\u200b \u200b", "<span style='font-size:0.01pt'>.</span> ", x, fixed = TRUE)
  x <- gsub("\u200b", "", x, fixed = TRUE)
  for (ch in names(.pdf_ascii_map)) {
    x <- gsub(ch, escape_label_markdown(.pdf_ascii_map[[ch]]), x, fixed = TRUE)
  }
  for (ch in names(.symbol_font_map)) {
    x <- gsub(
      ch,
      sprintf("<span style='font-family:%s'>%s</span>", .pdf_symbol_family, .symbol_font_map[[ch]]),
      x,
      fixed = TRUE
    )
  }
  names(x) <- nm
  x
}

# Built at top level so the closure captures only `lab` -- a closure made
# inside prepare_plot_for_pdf() would capture the plot itself, and the
# plot -> facet -> labeller -> plot cycle overflows ggplot_build().
.symbolising_labeller <- function(lab) {
  force(lab)
  f <- function(labels) lapply(lab(labels), symbolise_rich_labels)
  class(f) <- class(lab)
  f
}

#' Make a ggplot's rich labels drawable by base pdf()
#'
#' Returns a modified COPY (scales cloned) so the on-screen plot keeps its
#' Unicode labels. Only plots using element_markdown() are touched -- plain
#' element_text can't draw the Symbol-font spans.
#'
#' @param p A plot object; non-ggplots are returned unchanged.
#' @export
prepare_plot_for_pdf <- function(p) {
  if (!inherits(p, "ggplot") || !rich_labels_available()) {
    return(p)
  }
  els <- c("axis.text.x", "axis.text.y", "legend.text", "strip.text.x", "strip.text.y")
  rich <- any(vapply(els, function(e) inherits(p$theme[[e]], "element_markdown"), logical(1)))
  if (!rich) {
    return(p)
  }
  ensure_pdf_symbol_font()
  p$scales <- p$scales$clone()
  for (s in p$scales$scales) {
    if (is.character(s$labels)) s$labels <- symbolise_rich_labels(s$labels)
  }
  lab <- p$facet$params$labeller
  strip_rich <- inherits(p$theme[["strip.text.x"]], "element_markdown") ||
    inherits(p$theme[["strip.text.y"]], "element_markdown")
  if (strip_rich && is.function(lab)) {
    params <- p$facet$params
    params$labeller <- .symbolising_labeller(lab)
    # ggproto() resolves its parent lazily: bind it first, or the new facet
    # ends up inheriting from itself once p$facet is reassigned.
    parent_facet <- p$facet
    p$facet <- ggplot2::ggproto(NULL, parent_facet, params = params)
  }
  p
}

# ---------------------------------------------------------------------------
# Rich ComplexHeatmap legends
# ---------------------------------------------------------------------------

#' Describe a rich annotation legend for draw_heatmap() to build at draw time
#'
#' Attach a list of these to a HeatmapAnnotation as
#' `attr(ha, "rich_legends")` (and give the annotation `show_legend = FALSE`);
#' draw_heatmap() then builds them on the device actually drawing, so label
#' widths are measured -- and Greek symbolised for pdf() -- correctly.
#'
#' @param title Legend title.
#' @param labels Rich labels from level_labels(), named by level.
#' @param fill Colours, named by level.
#' @param title_gp,labels_gp grid gpar for the title and labels.
#' @export
heatmap_legend_spec <- function(
  title,
  labels,
  fill,
  title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  labels_gp = grid::gpar(fontsize = 10)
) {
  at <- names(labels)
  list(
    title = title,
    at = at,
    labels = unname(labels),
    fill = unname(fill[at]),
    title_gp = title_gp,
    labels_gp = labels_gp
  )
}

#' Build a rich ComplexHeatmap legend with a correctly measured width
#'
#' @param spec From heatmap_legend_spec().
#' @param pdf Rewrite Greek to the Symbol font for base pdf().
#' @return A ComplexHeatmap `Legends` object.
#' @export
build_rich_heatmap_legend <- function(spec, pdf = FALSE) {
  labels <- spec$labels
  if (pdf) {
    ensure_pdf_symbol_font()
    labels <- symbolise_rich_labels(labels)
  }
  lgd <- ComplexHeatmap::Legend(
    at = spec$at,
    labels = ComplexHeatmap::gt_render(labels),
    legend_gp = grid::gpar(fill = spec$fill),
    title = spec$title,
    title_gp = spec$title_gp,
    labels_gp = spec$labels_gp
  )
  # Legend() measured the labels as raw markup: re-measure the rendered text
  # and correct the three places the (inflated) width is stored.
  g <- lgd@grob
  body <- which(vapply(g$children, inherits, logical(1), "legend_body"))
  if (length(body) != 1) {
    return(lgd)
  }
  lab_w <- max(do.call(
    grid::unit.c,
    lapply(labels, function(x) {
      grid::grobWidth(gridtext::richtext_grob(x, gp = spec$labels_gp))
    })
  ))
  # key (4 mm default) + ComplexHeatmap's 1 mm label padding + slack
  body_w <- grid::convertWidth(grid::unit(4, "mm") + grid::unit(2, "mm") + lab_w, "mm")
  others <- setdiff(seq_along(g$children), body)
  title_w <- if (length(others)) {
    grid::convertWidth(grid::grobWidth(g$children[[others[1]]]), "mm")
  } else {
    grid::unit(0, "mm")
  }
  w <- max(body_w, title_w)
  attr(g$children[[body]], "width") <- body_w
  attr(g, "width") <- w
  if (!is.null(g$vp)) {
    g$vp$width <- w
  }
  lgd@grob <- g
  lgd
}

# Rich (gt_render) column labels for pdf(): swap Greek for Symbol-font spans
# in a COPY of the column-name annotation (its labels live in an environment
# shared with the on-screen object, so copy the environment too).
.symbolise_heatmap_column_labels <- function(ht) {
  if (inherits(ht, "HeatmapList")) {
    for (i in seq_along(ht@ht_list)) {
      if (inherits(ht@ht_list[[i]], "Heatmap")) ht@ht_list[[i]] <- .symbolise_heatmap_column_labels(ht@ht_list[[i]])
    }
    return(ht)
  }
  if (!inherits(ht, "Heatmap")) {
    return(ht)
  }
  a <- ht@column_names_param$anno
  if (is.null(a) || !inherits(a, "AnnotationFunction")) {
    return(ht)
  }
  e <- a@var_env
  v <- e$value
  if (!inherits(v, "gridtext")) {
    return(ht)
  }
  ne <- list2env(as.list(e, all.names = TRUE), parent = parent.env(e))
  ne$value <- structure(symbolise_rich_labels(as.character(v)), class = class(v), param = attr(v, "param"))
  a@var_env <- ne
  environment(a@fun) <- ne
  ht@column_names_param$anno <- a
  ht
}

# Rich legend specs attached to a Heatmap / HeatmapList's top annotations.
.heatmap_rich_legend_specs <- function(ht) {
  hts <- if (inherits(ht, "HeatmapList")) ht@ht_list else list(ht)
  unlist(
    lapply(hts, function(h) {
      if (!inherits(h, "Heatmap") || is.null(h@top_annotation)) {
        return(NULL)
      }
      attr(h@top_annotation, "rich_legends")
    }),
    recursive = FALSE
  )
}

#' Prepare a heatmap + draw() args for rich level labels on the current device
#'
#' Called by draw_heatmap() BEFORE any off-screen grab, so `pdf` reflects the
#' real target device.
#' @return list(ht, args)
#' @export
prepare_rich_heatmap <- function(ht, args) {
  pdf <- identical(names(grDevices::dev.cur()), "pdf")
  specs <- .heatmap_rich_legend_specs(ht)
  if (length(specs)) {
    args$annotation_legend_list <- c(args$annotation_legend_list, lapply(specs, build_rich_heatmap_legend, pdf = pdf))
  }
  if (pdf && rich_labels_available()) {
    ht <- .symbolise_heatmap_column_labels(ht)
  }
  list(ht = ht, args = args)
}
