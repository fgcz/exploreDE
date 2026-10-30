# ============================================================================
# Custom factor-level display labels (R/level_labels.R)
# ============================================================================

lab_input <- list(
  GroupLabelCondcKO_treated = "<i>Il4a</i><sup>fl/fl</sup> N.b.",
  GroupLabelCondKO_d = "&Delta;<i>Il4</i> &alpha;-CD3",
  GroupLabelCondWT = "   ",
  GroupLabelCondtypo = "<u>under</u>"
)
lab_levels <- c("WT", "cKO_treated", "KO_d", "typo", "1. odd")

test_that("escape_label_markdown makes every ASCII punctuation literal", {
  expect_equal(escape_label_markdown("cKO_treated*"), "cKO\\_treated\\*")
  expect_equal(escape_label_markdown("1. a<b"), "1\\. a\\<b")
  expect_equal(escape_label_markdown("déjà"), "déjà")
})

test_that("escaped raw levels render literally through gridtext", {
  skip_if_not_installed("gridtext")
  collect <- function(g) {
    if (inherits(g, "text")) {
      return(g$label)
    }
    if (!is.null(g$children)) {
      return(unlist(lapply(g$children, collect)))
    }
    NULL
  }
  for (s in c("cKO_treated", "*WT*", "a<b & c", "1. first", "- neg", "---", "[x](y)", "<i>x</i>")) {
    rendered <- paste(collect(gridtext::richtext_grob(escape_label_markdown(s))), collapse = "")
    expect_equal(gsub(" ", "", rendered), gsub(" ", "", s), info = s)
  }
})

test_that("decode_label_entities handles named, decimal and hex entities", {
  expect_equal(decode_label_entities("&Delta;x &alpha; &#946; &#x3B3; &amp;"), "Δx α β γ &")
  expect_equal(decode_label_entities("&unknown; plain"), "&unknown; plain")
  expect_equal(decode_label_entities(c("a", NA)), c("a", NA))
})

test_that("parse_level_label keeps whitelisted tags and escapes the rest", {
  p <- parse_level_label("<I>Il4a</I><sup>fl/fl</sup> N.b.")
  expect_true(p$ok)
  expect_equal(p$rich, "<i>Il4a</i><sup>fl\\/fl</sup> N\\.b\\.")
  expect_equal(p$plain, "Il4afl/fl N.b.")

  u <- parse_level_label("<u>under</u>")
  expect_equal(u$plain, "<u>under</u>") # not whitelisted: literal text
  expect_match(u$rich, "\\\\<u\\\\>", fixed = FALSE)

  br <- parse_level_label("top<br/>bottom")
  expect_equal(br$rich, "top<br>bottom")
  expect_equal(br$plain, "top\nbottom")
})

test_that("level_labels returns raw levels untouched when nothing is labelled", {
  out <- level_labels(list(), "Cond", lab_levels)
  expect_equal(as.character(out), lab_levels)
  expect_equal(names(out), lab_levels)
  expect_false(is_rich_labels(out))
  # legacy "format condition names" wrapping still uses \n in plain mode
  wrapped <- level_labels(list(), "Cond", c("cKO_treated", "WT"), wrap = TRUE)
  expect_equal(as.character(wrapped), c("cKO\ntreated", "WT"))
})

test_that("level_labels rich mode escapes unlabelled levels and ignores blanks", {
  skip_if_not(rich_labels_available())
  out <- level_labels(lab_input, "Cond", lab_levels)
  expect_true(is_rich_labels(out))
  expect_equal(out[["WT"]], "WT") # whitespace-only label ignored
  expect_equal(out[["1. odd"]], "1\\. odd") # escaped raw level
  expect_equal(out[["cKO_treated"]], "<i>Il4a</i><sup>fl\\/fl</sup> N\\.b\\.")
  expect_equal(out[["KO_d"]], "Δ<i>Il4</i> α\\-CD3")
  wrapped <- level_labels(lab_input, "Cond", c("WT", "a_b"), wrap = TRUE, force_rich = TRUE)
  expect_equal(wrapped[["a_b"]], "a<br>b")
  expect_equal(level_labels(list(), "X", "a.b_c*", wrap = TRUE, force_rich = TRUE)[[1]], "a<br>b<br>c\\*")
})

test_that("level_labels plain mode strips markup for ComplexHeatmap", {
  out <- level_labels(lab_input, "Cond", lab_levels, mode = "plain")
  expect_false(is_rich_labels(out))
  expect_equal(out[["cKO_treated"]], "Il4afl/fl N.b.")
  expect_equal(out[["KO_d"]], "ΔIl4 α-CD3")
  expect_equal(out[["1. odd"]], "1. odd")
})

test_that("force_rich escapes a second factor sharing the markdown legend", {
  skip_if_not(rich_labels_available())
  out <- level_labels(list(), "Other", c("a*b", "c"), force_rich = TRUE)
  expect_true(is_rich_labels(out))
  expect_equal(as.character(out), c("a\\*b", "c"))
})

test_that("PDF character checks flag only what base pdf() can't draw", {
  expect_true(label_has_greek("ΔIl4"))
  expect_false(label_has_greek("café"))
  expect_equal(label_pdf_unsupported_chars(c("Δx", "café", "a−b", "中")), "中")
})

test_that("symbolise_rich_labels swaps Greek for Symbol-font spans, keeping names", {
  out <- symbolise_rich_labels(c(k = "Δ<i>Il4</i> α\\-CD3 −"))
  expect_equal(names(out), "k")
  expect_match(out[["k"]], "<span style='font-family:GreekSym'>D</span>", fixed = TRUE)
  expect_match(out[["k"]], "<span style='font-family:GreekSym'>a</span>", fixed = TRUE)
  expect_false(grepl("−", out[["k"]]))
})

test_that("prepare_plot_for_pdf rewrites a copy and leaves the on-screen plot alone", {
  skip_if_not(rich_labels_available())
  lbl <- level_labels(lab_input, "Cond", lab_levels)
  df <- data.frame(g = factor(rep(lab_levels, each = 2), levels = lab_levels), y = seq_len(10), f = rep(c("A", "B"), 5))
  p <- ggplot2::ggplot(df, ggplot2::aes(g, y, fill = g)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(cols = ggplot2::vars(f), labeller = ggplot2::as_labeller(c(A = "β one", B = "B"))) +
    ggplot2::scale_x_discrete(labels = lbl) +
    ggplot2::scale_fill_manual(values = rainbow(5), breaks = lab_levels, labels = lbl) +
    ggprism::theme_prism(axis_text_angle = 45) +
    rich_label_theme(axis_x = TRUE, legend = TRUE, strip = TRUE)
  q <- prepare_plot_for_pdf(p)
  expect_equal(p$scales$get_scales("x")$labels[["KO_d"]], lbl[["KO_d"]])
  expect_match(q$scales$get_scales("x")$labels[["KO_d"]], "GreekSym", fixed = TRUE)
  expect_match(q$scales$get_scales("fill")$labels[["KO_d"]], "GreekSym", fixed = TRUE)
  # theme_prism's axis angle survives the markdown element
  expect_equal(ggplot2::calc_element("axis.text.x", ggplot2::theme_get() + p$theme)$angle, 45)
  f <- tempfile(fileext = ".pdf")
  on.exit(unlink(f))
  expect_no_error(export_plot_standard(p, f, format = "pdf", width = 6, height = 4))
  expect_gt(file.size(f), 0)
})

test_that("prepare_plot_for_pdf is a no-op for plain plots and non-ggplots", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(factor(cyl), mpg)) + ggplot2::geom_point()
  expect_identical(prepare_plot_for_pdf(p), p)
  expect_identical(prepare_plot_for_pdf("x"), "x")
})

test_that("whitespace between two tags survives (zero-width anchors, PDF-safe swap)", {
  p <- parse_level_label("<i>Vil1</i><sup>+/+</sup> <i>Il4ra</i>")
  expect_match(p$rich, "</sup>\u200b \u200b<i>", fixed = TRUE)
  expect_false(grepl("\u200b", symbolise_rich_labels(p$rich)))
  expect_equal(p$plain, "Vil1+/+ Il4ra")
})

test_that("build_rich_heatmap_legend corrects ComplexHeatmap's inflated width", {
  skip_if_not(rich_labels_available())
  labs <- level_labels(
    list(
      GroupLabelCa = "<i>Vil1</i><sup>CreER</sup> <i>Il4ra</i><sup>fl/fl</sup> (N.b.)",
      GroupLabelCb = "&Delta;<i>Il4</i>"
    ),
    "C",
    c("a", "b")
  )
  spec <- heatmap_legend_spec("C", labs, c(a = "red", b = "blue"))
  f <- tempfile(fileext = ".pdf")
  pdf(f)
  on.exit(
    {
      dev.off()
      unlink(f)
    },
    add = TRUE
  )
  raw <- ComplexHeatmap::Legend(
    at = spec$at,
    labels = ComplexHeatmap::gt_render(spec$labels),
    legend_gp = grid::gpar(fill = spec$fill),
    title = "C"
  )
  fixed <- build_rich_heatmap_legend(spec, pdf = TRUE)
  w_raw <- grid::convertWidth(attr(raw@grob, "width"), "mm", valueOnly = TRUE)
  w_fix <- grid::convertWidth(attr(fixed@grob, "width"), "mm", valueOnly = TRUE)
  expect_lt(w_fix, w_raw / 2)
  expect_equal(grid::convertWidth(fixed@grob$vp$width, "mm", valueOnly = TRUE), w_fix)
})

test_that("draw_heatmap builds rich legends and PDF column labels without touching the object", {
  skip_if_not(rich_labels_available())
  labs <- level_labels(list(GroupLabelCa = "&alpha;<i>x</i>"), "C", c("a", "b"))
  m <- matrix(rnorm(12), 3, dimnames = list(NULL, c("a", "b", "a2", "b2")))
  ha <- ComplexHeatmap::HeatmapAnnotation(
    C = c("a", "b", "a", "b"),
    col = list(C = c(a = "red", b = "blue")),
    show_legend = FALSE
  )
  attr(ha, "rich_legends") <- list(heatmap_legend_spec("C", labs, c(a = "red", b = "blue")))
  ht <- ComplexHeatmap::Heatmap(
    m,
    name = "z",
    top_annotation = ha,
    column_labels = ComplexHeatmap::gt_render(c("α<i>x</i>", "b", "c", "d"))
  )
  f <- tempfile(fileext = ".pdf")
  on.exit(unlink(f), add = TRUE)
  expect_no_error(export_plot_standard(function() draw_heatmap(ht, merge_legend = TRUE), f, format = "pdf"))
  expect_gt(file.size(f), 0)
  expect_equal(as.character(ht@column_names_param$anno@var_env$value)[1], "α<i>x</i>")
  prep <- prepare_rich_heatmap(ht, list())
  expect_length(prep$args$annotation_legend_list, 1)
})

test_that("heatmap PDFs with rich legends have exactly one page", {
  skip_if_not(rich_labels_available())
  pages <- function(f) sum(grepl("/Type /Page[^s]|/Type /Page$", readLines(f, warn = FALSE)))
  ha <- ComplexHeatmap::HeatmapAnnotation(
    C = c("a", "b", "a", "b"),
    col = list(C = c(a = "red", b = "blue")),
    show_legend = FALSE
  )
  attr(ha, "rich_legends") <- list(heatmap_legend_spec(
    "C",
    level_labels(list(GroupLabelCa = "<i>x</i><sup>fl</sup>"), "C", c("a", "b")),
    c(a = "red", b = "blue")
  ))
  ht <- ComplexHeatmap::Heatmap(matrix(rnorm(12), 3), name = "z", top_annotation = ha)
  for (dx in c(0, 5)) {
    f <- tempfile(fileext = ".pdf")
    export_plot_standard(function() draw_heatmap(ht, merge_legend = TRUE, legend_dx = dx), f, format = "pdf")
    expect_equal(pages(f), 1, info = paste("legend_dx =", dx))
    unlink(f)
  }
})
