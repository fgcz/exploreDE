# Draw the exploreDE hex sticker, man/figures/logo.png: a volcano plot of the
# demo dataset (inst/examples/genomics/deResult.rds) in the app's colours.
#
# Usage, from the package root: Rscript inst/scripts/make_hex_logo.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(SummarizedExperiment)
})

fill <- "#222d32" # app sidebar
border <- "#527085" # app header (adminlte light_blue)
up <- "#F39C12"
down <- "#6FA8DC"
ns <- "#6b7b84"

rd <- as.data.frame(rowData(readRDS("inst/examples/genomics/deResult.rds")))
# Square-root y axis: the p-values have a long tail (-log10 p up to ~300),
# which on a linear axis squeezes the volcano into a block.
rd <- rd[!is.na(rd$pValue) & !is.na(rd$log2Ratio) & abs(rd$log2Ratio) <= 8, ]
rd$y <- sqrt(-log10(rd$pValue))
rd <- rd[rd$y <= stats::quantile(rd$y, 0.995), ] # a few extreme points would flatten the rest
rd$x <- rd$log2Ratio
sig <- !is.na(rd$fdr) & rd$fdr < 0.05
rd$group <- ifelse(sig & rd$log2Ratio > 1, "up", ifelse(sig & rd$log2Ratio < -1, "down", "ns"))
rd <- rd[order(rd$group != "ns"), ] # significant points on top

grDevices::pdf(NULL) # ggplotGrob() needs an open device; this one writes no file
volcano <- ggplot(rd, aes(x, y, colour = group)) +
  geom_point(size = 0.1, alpha = 0.55, shape = 16) +
  scale_colour_manual(values = c(up = up, down = down, ns = ns), guide = "none") +
  scale_x_continuous(limits = c(-8, 8), expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  theme_void()

# Pointy-top hexagon of radius r, centred at the origin.
hexagon <- function(r) {
  a <- (90 + 60 * 0:5) * pi / 180
  data.frame(x = r * cos(a), y = r * sin(a))
}

sticker <- ggplot() +
  geom_polygon(aes(x, y), hexagon(1), fill = border) +
  geom_polygon(aes(x, y), hexagon(0.93), fill = fill) +
  annotation_custom(ggplotGrob(volcano), xmin = -0.64, xmax = 0.64, ymin = -0.24, ymax = 0.66) +
  annotate(
    "text",
    x = 0,
    y = -0.46,
    label = "exploreDE",
    colour = "white",
    size = 5.6,
    family = "Avenir Next",
    fontface = "bold"
  ) +
  coord_fixed(xlim = c(-sqrt(3) / 2, sqrt(3) / 2), ylim = c(-1, 1), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

dir.create("man/figures", showWarnings = FALSE)
ggsave(
  "man/figures/logo.png",
  sticker,
  width = sqrt(3),
  height = 2,
  dpi = 300,
  bg = "transparent",
  device = ragg::agg_png
)
cat("man/figures/logo.png:", file.size("man/figures/logo.png"), "bytes\n")
