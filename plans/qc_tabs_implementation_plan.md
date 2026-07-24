# exploreDE: QC Tabs Implementation Plan

**Author:** Peter Leary + Claude
**Created:** 2026-03-27
**Revised:** 2026-04-17
**Status:** Implemented
**Target:** exploreDE v3.6.0

## Overview

Add two new sidebar sections to exploreDE, porting key QC figures from the static
`twoGroups.Rmd` and `CountQC.Rmd` reports into interactive Shiny equivalents.
Both sections must work for RNASeq; applicable plots must also work for proteomics
and CountQC data types.

### Non-goals (already covered elsewhere in exploreDE)

| Static-report plot | Existing exploreDE equivalent |
|---|---|
| Sample correlation heatmap | Replicability tab |
| Sample clustering dendrogram | Replicability tab |
| All-pairs NxN scatter matrix (`allPairs`) | Replicability tab covers the purpose |
| MDS 2D/3D | Ordination tab (PCA/MDS/GLM-PCA) |
| High-variance clustering heatmap + GO | Heatmaps tab |
| EnrichR topMarkers links | Low ROI, skip |

## Data Availability Matrix

| Plot                           | RNASeq | CountQC | Proteomics |
|--------------------------------|:------:|:-------:|:----------:|
| Between-group scatter          |   ✓    |    —    |     ✓      |
| P-value histogram              |   ✓    |    —    |     ✓      |
| Intra-group scatter (Std)      |   ✓    |    ✓    |     ✓      |
| Intra-group scatter (GC/Width) |   ✓    |    ✓    |     —      |
| Total reads / above threshold  |   ✓    |    ✓    |     —      |
| Dominant genes boxplot         |   ✓    |    ✓    |     —      |
| GC / length LOESS bias scatter |   ✓    |    —    |     —      |

Two forms of bias visualisation:
- **Primary (all data types with `gc`/`featWidth`):** the **GC** and **Width**
  flavours of the intra-group and between-group scatter plots. Coloured bins
  drifting off the 1:1 line indicate bias. This is what CountQC.Rmd produces.
- **Additional (RNASeq DE only):** dedicated LOESS scatters of log2FC vs GC
  content and log2FC vs log10(feature width). Requires `log2Ratio` + `pValue`
  in seqAnno, so RNASeq DE only.

Key: ✓ = data available, — = not applicable / data absent.

---

## Shared: Scatter Plot Builder

All scatter plots (between-group, intra-group) share the same visual language
from ezRun's `ezXYScatter.2`: log10 axes, three typed colour layers, 1:1 solid
blue line + 2-fold dashed blue guide lines.

### Helper function: `make_qc_scatter()`

This goes in a shared helper file (either a new `helpers-qc.R` or appended to
`helperFunctions.R`). Both tabs use it.

```r
#' Build a QC scatter plot matching ezXYScatter.2 ggplot2 style
#' @param x Numeric vector (linear scale, NOT log2)
#' @param y Numeric vector (linear scale)
#' @param types_vec Character vector same length as x/y: level names from types_colours
#' @param types_colours Named character vector: level → colour
#' @param xlab,ylab Axis labels
#' @param title Plot title (NULL for none)
#' @param point_size,alpha Aesthetics
#' @param draw_order Character vector of level names defining z-order
#'        (last = drawn on top). Default: names(types_colours).
make_qc_scatter <- function(x, y, types_vec, types_colours,
                            xlab, ylab, title = NULL,
                            point_size = 1, alpha = 0.6,
                            base_size = 12, draw_order = NULL) {
  if (is.null(draw_order)) draw_order <- names(types_colours)
  df <- tibble(x = x, y = y, types = types_vec)
  df$types <- factor(df$types, levels = draw_order)
  df <- df |> arrange(types)  # last factor level drawn on top

  ggplot(df, aes(x = x, y = y, colour = types)) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_manual(values = types_colours) +
    scale_x_log10() + scale_y_log10() +
    geom_abline(intercept = 0, colour = "blue") +
    geom_abline(intercept = log10(2), colour = "blue", linetype = "dashed") +
    geom_abline(intercept = -log10(2), colour = "blue", linetype = "dashed") +
    theme_prism(base_size = base_size) +
    theme(legend.title = element_blank()) +
    labs(x = xlab, y = ylab, title = title)
}
```

### Colour-by modes (selectable via `radioButtons`)

All three modes share a common pattern: a `types_vec` of category labels
and a `types_colours` named-vector of colours built from user-controlled
`colourpicker::colourInput`s and `numericInput`s.

```r
# 1. Standard: Absent / Present / Significants
# RNASeq: usedInTest + significance thresholds
# CountQC: usedInTest only (no DE results → no Significants)
# Proteomics: all features are "Present" + significance from seqAnnoList
std_types <- ifelse(rd$usedInTest, "Present", "Absent")
is_sig <- !is.na(rd$pValue) &
          rd$pValue < input$scatterPval &
          abs(rd$log2Ratio) >= input$scatterLFC
std_types[is_sig] <- "Significants"
std_colours <- c(
  "Absent"       = input$scatterColAbsent,       # default: "grey"
  "Present"      = input$scatterColPresent,      # default: "black"
  "Significants" = input$scatterColSig           # default: "#A6CEE3"
)

# 2. GC bins — thresholds user-configurable
gc_lo <- input$scatterGCLow   # default 0.4
gc_hi <- input$scatterGCHigh  # default 0.6
gc_types <- rep(sprintf("%.2f–%.2f", gc_lo, gc_hi), nrow(rd))
gc_types[rd$gc < gc_lo] <- sprintf("GC < %.2f", gc_lo)
gc_types[rd$gc > gc_hi] <- sprintf("GC > %.2f", gc_hi)
gc_colours <- setNames(
  c(input$scatterColGCLow,   # default: brewer "Paired" blue
    input$scatterColGCMid,   # default: "black"
    input$scatterColGCHigh), # default: brewer "Paired" orange
  c(sprintf("GC < %.2f", gc_lo),
    sprintf("%.2f–%.2f", gc_lo, gc_hi),
    sprintf("GC > %.2f", gc_hi))
)
# Draw outlier bins last so they're visible on top
draw_order_gc <- c(names(gc_colours)[2], names(gc_colours)[c(1, 3)])

# 3. Width bins — thresholds user-configurable
w_lo <- input$scatterWidthLow   # default 500
w_hi <- input$scatterWidthHigh  # default 5000
width_types <- rep(sprintf("%d–%d nt", w_lo, w_hi), nrow(rd))
width_types[rd$featWidth < w_lo] <- sprintf("< %d nt", w_lo)
width_types[rd$featWidth > w_hi] <- sprintf("> %d nt", w_hi)
width_colours <- setNames(
  c(input$scatterColWidthLow,
    input$scatterColWidthMid,
    input$scatterColWidthHigh),
  c(sprintf("< %d nt", w_lo),
    sprintf("%d–%d nt", w_lo, w_hi),
    sprintf("> %d nt", w_hi))
)
draw_order_w <- c(names(width_colours)[2], names(width_colours)[c(1, 3)])
```

### Signal computation

Scatter plots use **linear signal** (not log2), matching ezRun:

```r
# RNASeq / CountQC:
norm_assay <- grep("^signal$|^xNorm$", names(assays(se)), value = TRUE)
signal <- assay(se, norm_assay) + param$backgroundExpression

# Proteomics:
# countList assays are already in the right space
signal <- inputDataReactive()$countList[[input$scatterCountType]]
```

---

## Tab 1: Count QC

**Sidebar section:** Counts (alongside Ordination, Heatmaps, etc.)
**Tab name:** `countQcTab`
**Visibility:** All data types (conditional panels within)
**Files to create:**
- `inst/app/ui-tab-countQc.R`
- `inst/app/server-countQc.R`

### 1.1 Total Reads & Features Above Threshold bar charts

![Total Reads](plan_figures/01_total_reads.png)
![Features Present](plan_figures/02_features_present.png)

**Data type:** RNASeq, CountQC only
**Source:** `CountQC.Rmd` lines 226–253
**Visibility:** `conditionalPanel` hidden for proteomics (no raw counts / presentFlag)

**Data access:**
```r
se <- inputDataReactive()$se
raw_counts <- assay(se, "counts")
present_flag <- assay(se, "presentFlag")
dataset <- inputDataReactive()$dataset

toPlot <- tibble(
  Sample = factor(colnames(raw_counts), levels = colnames(raw_counts)),
  Total_M = colSums(raw_counts) / 1e6,
  Present = colSums(present_flag),
  Pct = paste0(signif(100 * Present / nrow(raw_counts), 3), "%"),
  Group = dataset[[input$barColourBy]][match(colnames(raw_counts), rownames(dataset))]
)
```

**Static (ggplot2):**
```r
p_total <- ggplot(toPlot, aes(x = Sample, y = Total_M, fill = Group)) +
  geom_col(colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = group_colours) +
  theme_prism(base_size = input$textSizeCountQc, axis_text_angle = 45) +
  labs(title = "Total Reads per Sample", y = "Counts (millions)", x = NULL)

p_present <- ggplot(toPlot, aes(x = Sample, y = Present, fill = Group)) +
  geom_col(colour = "black", linewidth = 0.3) +
  geom_text(aes(label = Pct), vjust = -0.3, size = 3, fontface = "bold") +
  scale_fill_manual(values = group_colours) +
  theme_prism(base_size = input$textSizeCountQc, axis_text_angle = 45) +
  labs(title = "Features Above Threshold", y = "Count", x = NULL)
```

**UI controls:**
- Colour by factor (`selectInput` from `inputDataReactive()$factors`)
- Uses app's group colour system (`input$GroupColour{Factor}{Level}`)
- Font size slider
- Download button (plot + data table as xlsx)

### 1.2 Dominant Genes Boxplot

**Data type:** RNASeq, CountQC
**Source:** `CountQC.Rmd` lines 255–283 (extended with gene-family highlighter)

**Concept:** Per-sample boxplot of `fraction of total expression` across all
genes. Highlights the top-N over-represented genes per sample (default behaviour
from CountQC.Rmd) **and** allows the user to type a pattern / pick a preset to
flag specific gene families regardless of their per-sample rank. Useful for
spotting blood contamination (HB genes), incomplete rRNA depletion (RP[LS]),
high mitochondrial content (MT-), immune contamination (IG[HKL]), etc.

**Data access:**
```r
counts <- assay(se, "counts")
gene_names <- rowData(se)$gene_name
# use unique names; fall back to rownames when gene_name is missing
labels <- ifelse(is.na(gene_names) | gene_names == "",
                 rownames(counts), gene_names)
labels <- make.unique(labels)

countFrac <- sweep(counts, 2, colSums(counts), FUN = "/")
rownames(countFrac) <- labels

df <- as.data.frame(countFrac) |>
  tibble::rownames_to_column("Gene") |>
  pivot_longer(-Gene, names_to = "Sample", values_to = "Fraction")
```

**Top-N default layer (from CountQC.Rmd):**
```r
topN <- df |>
  group_by(Sample) |>
  slice_max(Fraction, n = input$domTopN, with_ties = FALSE) |>  # default 3
  mutate(rank = factor(row_number(), levels = seq_len(input$domTopN))) |>
  ungroup()
```

**Gene-family highlight layer (new):**
```r
# input$domPattern: free-text regex, e.g. "^MT-|^HB[AB]"
# input$domPresets: multi-select of curated presets; union them with user regex
presets <- list(
  "Mitochondrial"    = "^MT-",
  "Ribosomal (cyto)" = "^RP[LS]",
  "Ribosomal (mito)" = "^MRP[LS]",
  "Haemoglobin"      = "^HB[ABEGMQZ]",
  "Immunoglobulin"   = "^IG[HKL][VDJC]",
  "TCR"              = "^TR[ABGD][VDJC]",
  "Heat shock"       = "^HSP",
  "Histone"          = "^HIST"
)
regex_parts <- c(input$domPattern, unlist(presets[input$domPresets]))
regex_parts <- regex_parts[nzchar(regex_parts)]
hit_regex <- if (length(regex_parts)) paste(regex_parts, collapse = "|") else NA_character_

family <- df |>
  filter(!is.na(hit_regex) & grepl(hit_regex, Gene, ignore.case = TRUE))
```

**Static (ggplot2):**
```r
p <- ggplot(df, aes(x = Sample, y = Fraction)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  # Top-N ranked points (CountQC.Rmd behaviour)
  geom_point(data = topN, aes(colour = rank), size = 2,
             position = position_dodge(width = 0.4)) +
  geom_text_repel(data = topN, aes(label = Gene, colour = rank),
                  size = 2.5, nudge_y = 0.02, direction = "y",
                  segment.size = 0.2, segment.alpha = 0.4,
                  max.overlaps = Inf) +
  scale_colour_manual(values = c("red", "orange", "gold"),
                      name = "Top-gene rank") +
  # Family hits (new layer — separate scale needed)
  ggnewscale::new_scale_colour() +
  geom_point(data = family, aes(colour = "Family hit"),
             size = 2, shape = 17, alpha = 0.8) +
  geom_text_repel(data = family, aes(label = Gene), colour = "#1965B0",
                  size = 2.2, max.overlaps = Inf, segment.alpha = 0.3) +
  scale_colour_manual(values = c("Family hit" = "#1965B0"), name = NULL) +
  theme_prism(base_size = input$textSizeCountQc, axis_text_angle = 90) +
  labs(y = "Fraction of total expression",
       title = "Dominant genes per sample")
```

**UI controls:**
- `numericInput("domTopN", "Top-N per sample", value = 3, min = 0, max = 20)`
  (0 disables top-N, leaves only family highlights)
- `checkboxGroupInput("domPresets", "Gene family presets", choices = names(presets))`
- `textInput("domPattern", "Custom regex (OR'd with presets)")`
- Font size, point size sliders
- Download (plot + tidy data as xlsx)

### 1.3 Intra-group Scatter Plots

Three flavours matching CountQC's `countQcScatterPlots()` → `ezScatter()`:
Standard (Absent/Present/Significants), GC bins, gene width bins.
Log10 axes, 1:1 + 2-fold blue guide lines, geometric mean on x-axis.

![Intra-group Standard](plan_figures/download (27).png)
![Intra-group GC](plan_figures/download (28).png)
![Intra-group Width](plan_figures/download (29).png)

**Data type:** All three (GC/width flavours only when `gc`/`featWidth` available)
**Source:** `CountQC.Rmd` lines 591–735 (`countQcScatterPlots` → `ezScatter`)

**Concept:** For each group level, plot each replicate's expression (y) against
the group **geometric** mean (x). Reveals outlier samples and technical bias.

**Data access:**
```r
signal <- assay(se, norm_assay) + param$backgroundExpression  # linear scale

group <- input$intraGroup
idx <- which(dataset[[grouping]] == group)
samples <- colnames(signal)[idx]

# Geometric mean (matching ezScatter)
group_geomean <- 2^rowMeans(log2(signal[, idx, drop = FALSE]))
```

**Static (ggplot2) — one panel per sample using `make_qc_scatter`:**
```r
plots <- lapply(samples, function(s) {
  make_qc_scatter(
    x = group_geomean,
    y = signal[, s],
    types_vec = current_types,        # from colour-by selector
    types_colours = current_colours,  # from colour-by selector
    xlab = paste("Avg of", group),
    ylab = s, title = s,
    point_size = input$pointSizeCountQc,
    alpha = input$alphaCountQc,
    base_size = input$textSizeCountQc,
    draw_order = current_draw_order   # outlier bins drawn last
  ) + theme(legend.position = "none")
})
cowplot::plot_grid(plotlist = plots, ncol = min(length(plots), 4))
```

**UI controls** (shared Standard/GC/Width):
- `selectInput("intraGroup", ...)` — group level to display (levels of `groupingName`)
- `radioButtons("scatterColourBy", choices = c("Standard", "GC content", "Gene width"))`
  - GC / Width options hidden when columns absent
  - For CountQC: Standard has no Significants (just Absent/Present)
  - For Proteomics: no GC/Width options
- **Standard colouring** — 3 `colourpicker::colourInput`s (Absent, Present, Significants)
  + P-value threshold, LFC threshold
- **GC colouring** — 2 `numericInput`s (low, high; defaults 0.4 / 0.6)
  + 3 `colourpicker::colourInput`s (low-bin, mid-bin, high-bin)
- **Width colouring** — 2 `numericInput`s (low, high; defaults 500 / 5000)
  + 3 `colourpicker::colourInput`s (low-bin, mid-bin, high-bin)
- Point size slider (default 0.5)
- Alpha slider (default 0.4)
- Font size slider
- Figure dimensions, download

**Edge cases:**
- Single-sample group: show notification, skip (can't compute group mean from 1 sample)
- >6 replicates: grid wraps (ncol capped at 4)
- Proteomics: `is_present` may not exist — all features treated as "Present"

---

## Tab 2: DE QC

**Sidebar section:** DE Results (alongside DE Table, Volcano)
**Tab name:** `deQcTab`
**Visibility:** `dataType != "CountQC"` (same as Volcano, DE Table)
**Files to create:**
- `inst/app/ui-tab-deQc.R`
- `inst/app/server-deQc.R`

### 2.1 Between-group Scatter Plot

Single pair (`refGroup` vs `sampleGroup`) with three flavours — Standard, GC, Width.
Log10 axes, 1:1 + 2-fold guide lines.

![Between-group Scatter — all three flavours](plan_figures/04_between_group_scatter.png)

**Data type:** RNASeq, Proteomics (GC/width flavours RNASeq only)
**Source:** `twoGroups.Rmd` lines 168–198 + `countQcScatterPlots` pattern

**Data access (RNASeq):**
```r
signal <- assay(se, norm_assay) + param$backgroundExpression
sample_idx <- which(dataset[[grouping]] == param$sampleGroup)
ref_idx <- which(dataset[[grouping]] == param$refGroup)
sample_mean <- rowMeans(signal[, sample_idx, drop = FALSE])
ref_mean <- rowMeans(signal[, ref_idx, drop = FALSE])
```

**Data access (Proteomics):**
```r
# Proteomics has no param$sampleGroup/refGroup.
# User selects two levels from the first factor.
signal <- inputDataReactive()$countList[[input$scatterCountType]]
sample_idx <- which(dataset[[grouping]] == input$scatterSampleGroup)
ref_idx <- which(dataset[[grouping]] == input$scatterRefGroup)
sample_mean <- rowMeans(signal[, sample_idx, drop = FALSE])
ref_mean <- rowMeans(signal[, ref_idx, drop = FALSE])
```

**Static (ggplot2) — using `make_qc_scatter`:**
```r
p <- make_qc_scatter(
  x = ref_mean, y = sample_mean,
  types_vec = current_types,
  types_colours = current_colours,
  xlab = param$refGroup,
  ylab = param$sampleGroup,
  title = "Comparison of average expression",
  base_size = input$textSizeDeQc,
  point_size = input$pointSizeDeQc,
  alpha = input$alphaDeQc,
  draw_order = current_draw_order
)
```

**Interactive (plotly):** Phase 3 polish — `scattergl` wrapping for large datasets,
with log10 axes, blue 1:1 + 2-fold shape lines, feature_id tooltip.

**UI controls:** Same pattern as intra-group scatter (see 1.3).
Additionally for proteomics: `selectInput`s for which two levels to compare.

### 2.2 P-value Histogram

![P-value Histogram](plan_figures/03_pvalue_histogram.png)

**Data type:** RNASeq, Proteomics
**Source:** `twoGroups.Rmd` lines 262–280

**Concept:** Histogram of raw p-values. A well-calibrated test shows a uniform
distribution with a spike near 0 (true positives). The dashed horizontal line
shows the expected count per bin under the null. Modes elsewhere indicate
model misspecification. Orange bars = features not expressed / not used in test.

**Data access (RNASeq):**
```r
seqAnno <- inputDataReactive()$seqAnno
pvals <- seqAnno$pValue
categories <- ifelse(seqAnno$usedInTest, "Used in test", "Not expressed")
n_used <- sum(seqAnno$usedInTest, na.rm = TRUE)
```

**Data access (Proteomics):**
```r
seqAnnoFilt <- inputDataReactive()$seqAnnoList[[input$contrastSelected]]
pvals <- seqAnnoFilt$pValue
# No usedInTest flag — single colour, or split by modelName if available:
if ("modelName" %in% colnames(seqAnnoFilt)) {
  categories <- ifelse(grepl("^Linear", seqAnnoFilt$modelName), "Linear model", "Imputed")
} else {
  categories <- rep("All features", length(pvals))
}
n_used <- sum(!is.na(pvals))
```

**Static (ggplot2):**
```r
df <- tibble(pValue = pvals, Category = categories) |> filter(!is.na(pValue))
nbins <- input$phistBins

p <- ggplot(df, aes(x = pValue, fill = Category)) +
  geom_histogram(bins = nbins, boundary = 0, colour = "white", linewidth = 0.15) +
  geom_hline(yintercept = n_used / nbins, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  scale_fill_manual(values = c(
    "Used in test" = input$phistColUsed,     # default "#1965B0"
    "Not expressed" = input$phistColNotUsed, # default "#FF7B00"
    "Linear model" = input$phistColUsed,
    "Imputed" = input$phistColNotUsed,
    "All features" = input$phistColUsed
  )) +
  theme_prism(base_size = input$textSizeDeQc) +
  labs(x = "P-value", y = "Count", title = "P-value Distribution",
       subtitle = "Dashed line = expected count under the null (uniform)")
```

**UI controls:**
- Number of bins (slider, min 50, max 500, default 200)
- For proteomics: contrast selector (reuses existing `input$contrastSelected`)
- Toggle raw p-value vs FDR
- 2× `colourpicker::colourInput` (Used / Not used)
- Font size
- Download

### 2.3 GC Content & Gene Length LOESS Bias Scatter (additional)

![GC Bias](plan_figures/06_gc_bias.png)
![Gene Length Bias](plan_figures/07_length_bias.png)

**Data type:** RNASeq DE only (requires `gc`, `featWidth`, `log2Ratio`,
`pValue` in seqAnno). This is *additional* to the GC/Width flavoured scatters
in §1.3 and §2.1 — those show bias by position of coloured bins; these show it
directly as a LOESS trend of log2FC vs GC (or length).
**Source:** `twoGroups.Rmd` technical bias section (verify path before
implementation — I attributed this to lines 889–966 in the original plan but
haven't confirmed).

**Concept:** Scatter of GC content (or gene length) vs log2FC, coloured by
significance, with a LOESS trend line. If the LOESS deviates from y=0, DE
results may be confounded by the technical feature.

**Data access:**
```r
seqAnno <- inputDataReactive()$seqAnno
has_gc <- "gc" %in% colnames(seqAnno) && !all(is.na(seqAnno$gc))
has_width <- "featWidth" %in% colnames(seqAnno) && !all(is.na(seqAnno$featWidth))

df <- seqAnno[seqAnno$usedInTest & !is.na(seqAnno$pValue), ]
df$types <- ifelse(
  abs(df$log2Ratio) >= input$biasLFC & df$pValue <= input$biasPval,
  "Significants", "Present"
)
df$types <- factor(df$types, levels = c("Present", "Significants"))
```

**Static (ggplot2):**
```r
bias_colours <- c(
  "Present"      = input$biasColPresent, # default "black"
  "Significants" = input$biasColSig      # default "#A6CEE3"
)

p_gc <- ggplot(df, aes(x = gc, y = log2Ratio, colour = types)) +
  geom_point(size = input$biasPointSize, alpha = input$biasAlpha) +
  scale_color_manual(values = bias_colours) +
  geom_hline(yintercept = 0, colour = "blue") +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE,
              colour = "#DC050C", fill = "#DC050C", alpha = 0.15,
              linewidth = 0.7, formula = y ~ x) +
  theme_prism(base_size = input$textSizeDeQc) +
  theme(legend.title = element_blank()) +
  labs(x = "GC Content", y = "Log2 Fold Change", title = "GC Content vs Fold Change")

# Gene length: same but x = log10(featWidth)
```

**UI controls:**
- P-value threshold for significance colouring (numericInput, default 0.05)
- LFC threshold (numericInput, default 1)
- Toggle LOESS line on/off
- 2× `colourpicker::colourInput` (Present / Significants)
- Font size, point size, alpha sliders
- Download

**Visibility:** wrap in `conditionalPanel` gated on `dataType == "RNASeq"` AND
`has_gc` / `has_width` (via server-side `uiOutput`).

---

## UI Layout Specification

### Tab 1: Count QC (`ui-tab-countQc.R`)

```
tabItem(tabName = "countQcTab",
  fluidRow(
    column(width = 3,
      box(title = "Settings", width = NULL, status = "primary", solidHeader = TRUE,
        # Bar chart section (conditional on RNASeq/CountQC)
        uiOutput("countQcBarSettingsUI"),
        hr(),
        # Dominant genes section (conditional on RNASeq/CountQC)
        h4("Dominant genes"),
        numericInput("domTopN", ..., value = 3),
        checkboxGroupInput("domPresets", ..., choices = preset_names),
        textInput("domPattern", ...),
        hr(),
        # Intra-group scatter section
        h4("Intra-group scatter"),
        selectInput("intraGroup", ...),
        radioButtons("scatterColourBy", ...),
        # Conditional panels per colour-by mode:
        conditionalPanel("input.scatterColourBy == 'Standard'",
          numericInput("scatterPval", value = 0.05), numericInput("scatterLFC", value = 1),
          colourInput("scatterColAbsent"), colourInput("scatterColPresent"), colourInput("scatterColSig")),
        conditionalPanel("input.scatterColourBy == 'GC content'",
          numericInput("scatterGCLow", value = 0.4), numericInput("scatterGCHigh", value = 0.6),
          colourInput("scatterColGCLow"), colourInput("scatterColGCMid"), colourInput("scatterColGCHigh")),
        conditionalPanel("input.scatterColourBy == 'Gene width'",
          numericInput("scatterWidthLow", value = 500), numericInput("scatterWidthHigh", value = 5000),
          colourInput("scatterColWidthLow"), colourInput("scatterColWidthMid"), colourInput("scatterColWidthHigh")),
        selectInput("intraCountType", ..., choices from countList names),
        hr(),
        # Shared aesthetics
        sliderInput("textSizeCountQc", ...),
        sliderInput("pointSizeCountQc", ...),
        sliderInput("alphaCountQc", ...),
        numericInput("figHeightCountQc", ...), numericInput("figWidthCountQc", ...)
      )
    ),
    column(width = 9,
      uiOutput("countQcBarPlotsUI"),              # bar charts (conditional)
      uiOutput("countQcDominantGenesUI"),         # dominant genes boxplot (conditional)
      box(title = "Intra-group Scatter", width = NULL, status = "primary", solidHeader = TRUE,
        plotOutput("intraGroupScatter"),
        downloadButton("dlPlot_intraScatter", ...)
      )
    )
  )
)
```

### Tab 2: DE QC (`ui-tab-deQc.R`)

```
tabItem(tabName = "deQcTab",
  fluidRow(
    column(width = 3,
      box(title = "Settings", width = NULL, status = "primary", solidHeader = TRUE,
        h4("Between-group scatter"),
        radioButtons("deQcScatterColourBy", ...),
        # Same 3 conditional panels as Count QC (Standard/GC/Width)
        uiOutput("deQcScatterColourUI"),
        uiOutput("deQcGroupSelectorsUI"),          # proteomics only
        hr(),
        h4("P-value histogram"),
        sliderInput("phistBins", ..., value = 200),
        radioButtons("phistType", choices = c("Raw" = "pValue", "FDR" = "fdr")),
        colourInput("phistColUsed"), colourInput("phistColNotUsed"),
        hr(),
        # LOESS bias (RNASeq DE only, conditional)
        uiOutput("deQcBiasSettingsUI"),
        #   Contains: numericInput("biasPval"), numericInput("biasLFC"),
        #             checkboxInput("biasLoess"),
        #             colourInput("biasColPresent"), colourInput("biasColSig"),
        #             sliderInput("biasPointSize"), sliderInput("biasAlpha")
        hr(),
        sliderInput("textSizeDeQc", ...),
        sliderInput("pointSizeDeQc", ...),
        sliderInput("alphaDeQc", ...),
        numericInput("figHeightDeQc", ...), numericInput("figWidthDeQc", ...)
      )
    ),
    column(width = 9,
      box(title = "Between-group Scatter", width = NULL, status = "primary", solidHeader = TRUE,
        tabsetPanel(
          tabPanel("Static", plotOutput("betweenGroupScatter")),
          tabPanel("Interactive", plotlyOutput("betweenGroupScatterPlotly"))
        ),
        downloadButton("dlPlot_betweenScatter", ...)
      ),
      box(title = "P-value Distribution", width = NULL, status = "primary", solidHeader = TRUE,
        plotOutput("pvalHistogram"),
        downloadButton("dlPlot_pvalHist", ...)
      ),
      uiOutput("deQcBiasUI")                      # LOESS bias plots (RNASeq DE only)
    )
  )
)
```

---

## Sidebar Menu Integration

In `app.R`, the sidebar menu is dynamically generated in `output$sidebarMenu`.

```r
# In the Counts section (always shown, alongside Ordination, Heatmaps, etc.):
menuItem(text = "Count QC", tabName = "countQcTab", icon = icon("chart-bar"))

# In the DE Results section (dataType != "CountQC"):
menuItem(text = "DE QC", tabName = "deQcTab", icon = icon("chart-line"))
```

Both UI tab files need to be sourced in `app.R`'s `body` definition alongside
existing `source("ui-tab-*.R")` calls. Both server files need to be sourced in
the `server` function alongside existing `source("server-*.R")` calls.

---

## File Structure

```
inst/app/
├── ui-tab-countQc.R      # NEW — Count QC tab UI
├── ui-tab-deQc.R         # NEW — DE QC tab UI
├── server-countQc.R      # NEW — Count QC server logic
├── server-deQc.R         # NEW — DE QC server logic
```

`make_qc_scatter()` lives in an existing helpers file (append to `R/helperFunctions.R`
in the exploreDE package, or a new `inst/app/helpers-qc.R` sourced at startup).

---

## Edge Cases

| Scenario | Handling |
|----------|----------|
| Single-sample group | Skip intra-group for that level; show notification |
| >6 replicates per group | Grid wraps (ncol capped at 4) |
| Missing `gc` column | Hide GC colour-by option and GC LOESS bias plot |
| Missing `featWidth` column | Hide Width colour-by option and length LOESS bias plot |
| Proteomics: no `presentFlag` / `isPresent` | All features = "Present" in type vector |
| Proteomics: no `sampleGroup`/`refGroup` | User selects two levels from factor dropdown |
| CountQC: no DE results | Standard colouring = Absent/Present only (no Significants) |
| Isoform-level: 170k features | `scattergl` for plotly performance; warn if >100k |
| All p-values NA | Show empty histogram with notification |
| `presentFlag` vs `isPresent` assay | Check both: `assay(se, "isPresent")` then fall back to `assay(se, "presentFlag")` |
| Gene-family regex matches 0 genes | Silent; top-N layer still shown |

---

## Implementation Phases

### Phase 1: DE QC tab — RNASeq, Standard colouring only
1. Create `ui-tab-deQc.R` with layout (Standard panel, p-value hist, LOESS bias section)
2. Implement `server-deQc.R`:
   - `make_qc_scatter()` helper
   - Between-group scatter (RNASeq, Standard colouring)
   - P-value histogram (RNASeq)
   - GC / length LOESS bias scatters (RNASeq, conditional on gc/featWidth)
3. Wire into `app.R` (sidebar + source)
4. Test with gene-level RNASeq data

### Phase 2: DE QC tab — GC/Width flavours + proteomics
1. Add GC and Width colour-by modes to between-group scatter
2. Add proteomics paths (contrast selector, group selectors, p-value histogram)
3. Test with proteomics (protein + phospho) and isoform data

### Phase 3: Count QC tab
1. Create `ui-tab-countQc.R`
2. Implement `server-countQc.R`:
   - Total reads / features above threshold bar charts (RNASeq + CountQC)
   - Dominant genes boxplot (top-N + family highlighter)
   - Intra-group scatter with 3 colour-by modes
3. Wire into `app.R`
4. Test with RNASeq, CountQC, and proteomics data

### Phase 4: Polish
- Download handlers (plot + data export as xlsx) for all plots
- Plotly interactive variants for scatter plots (`scattergl` for >100k features)
- Info boxes with interpretation guidance for each plot
- Debounce on threshold inputs
- Consistent styling audit against existing tabs

---

## Open Questions

1. **Intra-group for proteomics:** Proteomics factor levels may not map cleanly
   to "groups" in the same way. Use the first factor by default, let user pick
   via a `selectInput`.

2. **ezRun functions:** `ezXYScatter.2`, `ezScatter`, `countQcScatterPlots` are
   available. Recommendation: **reimplement in pure ggplot2** for Shiny — we need
   reactive inputs (colour-by, thresholds, font size, point aesthetics) that
   don't map to ezRun's static API. The `make_qc_scatter()` helper gives us full
   control while matching the visual style.

3. **P-value histogram for proteomics with `modelName`:** Split by Linear/Imputed
   (same colour scheme as Used/Not expressed). Single colour for datasets without
   `modelName`.

4. **Dominant genes gene-family presets:** the curated list in §1.2 covers the
   common cases (MT, rRNA, HB, Ig, TCR, HSP, histone). Add species-specific
   presets (e.g. mouse lowercase `mt-`, `Rp[ls]`, `Hb[ab]`) via a case-insensitive
   match, or keep strict capitalisation and ship human + mouse variants?

---

## Test Datasets

```r
# RNASeq (gene-level, standard)
dataDir <- "/srv/gstore/projects/p36257/o38191_o39085_o40237_DESeq2_Ang_CTRL_T2dot5--over--CTRL_T2dot5_2026-03-06--13-16-47/Ang_CTRL_T2.5--over--CTRL_T2.5"

# RNASeq (isoform-level, 170k features)
dataDir <- "/srv/gstore/projects/p34429/o35733_DESeq2_Palbo_day6--over--Palbo_day2_2024-08-22--11-21-44/Palbo_day6--over--Palbo_day2/"

# CountQC
dataDir <- "/srv/gstore/projects/p41168/o41216_CountQC_2026-03-27--12-03-12/Count_QC/"

# Proteomics (protein-level)
dataDir <- "https://fgcz-ms.uzh.ch/public/pStore/p39819/bfabric/Proteomics/exploreDE/2025/2025-10/2025-10-22/workunit_335485/2997996.rds"

# Proteomics (phospho/site-level)
dataDir <- "https://fgcz-ms.uzh.ch/public/pStore/p23078/bfabric/Proteomics/exploreDE/2026/2026-03/2026-03-10/workunit_342065/3106968.rds"
```

---

## Generating prototype figures

Prototype figures in `plan_figures/`:
- `01_total_reads.png`, `02_features_present.png` — bar charts
- `03_pvalue_histogram.png` — p-value histogram
- `04_between_group_scatter.png` — between-group, all 3 flavours
- `05a/b/c_intra_group_*.png` — early intra-group prototypes
- `06_gc_bias.png`, `07_length_bias.png` — LOESS bias scatters (RNASeq DE, additional)
- `download (27/28/29).png` — real CountQC intra-group scatter (Std/GC/Width)
   showing desired visual style with better GC<0.4/>0.6 and width<500/>5000 highlighting
- `CountQC (4).Rmd` — authoritative source for CountQC plot code

Regenerate with:
```bash
module load Dev/R/4.5.0
Rscript plan_figures/generate_figures.R
```
