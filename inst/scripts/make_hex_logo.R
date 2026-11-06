# Draw man/figures/logo.png: volcano points suggest a cat's face.
# Usage, from the package root: Rscript inst/scripts/make_hex_logo.R
# An absolute script path also works from any directory.
# This is a stylised logo, not a statistical plot: facial points are decorative.

suppressPackageStartupMessages({
  library(ggplot2)
  library(SummarizedExperiment)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
package_dir <- normalizePath(file.path(dirname(normalizePath(script_file)), "../.."))
output_dir <- file.path(package_dir, "man/figures")
data_file <- file.path(package_dir, "inst/examples/genomics/deResult.rds")

fill <- "#222d32"
border <- "#527085"
up <- "#F39C12"
down <- "#6FA8DC"
ns <- "#6b7b84"

rd <- as.data.frame(rowData(readRDS(data_file)))
rd <- rd[!is.na(rd$pValue) & !is.na(rd$log2Ratio) & abs(rd$log2Ratio) <= 8, ]
rd$y <- sqrt(-log10(rd$pValue))
rd <- rd[rd$y <= stats::quantile(rd$y, 0.995), ]
rd$x <- rd$log2Ratio / 8
rd$y <- rd$y / max(rd$y)
sig <- !is.na(rd$fdr) & rd$fdr < 0.05
rd$colour <- ifelse(sig & rd$x > 1 / 8, up, ifelse(sig & rd$x < -1 / 8, down, ns))
rd <- rd[order(rd$colour != ns), ]

# A little irregularity keeps the facial details in the same language as the
# volcano: individual points, never a solid cartoon outline.
set.seed(27)
dot_curve <- function(x, y, colour, size = 0.48, alpha = 0.9, jitter = 0.002) {
  data.frame(
    x = x + stats::rnorm(length(x), sd = jitter),
    y = y + stats::rnorm(length(y), sd = jitter),
    colour = colour,
    size = size,
    alpha = alpha
  )
}

features <- list()
for (side in c(-1, 1)) {
  colour <- if (side < 0) down else up
  eye_x <- 0.29
  eye_y <- 0.36
  eye_width <- 0.14
  eye_height <- 0.065

  # Open almond eyes, outlined by points with a dotted vertical pupil.
  eye_distance <- abs(rd$x - side * eye_x) / eye_width
  in_eye <- eye_distance < 1 & abs(rd$y - eye_y) < eye_height * (1 - eye_distance^1.6) + 0.02
  rd <- rd[!in_eye, ]
  t <- seq(-1, 1, length.out = 15)
  for (arc in c(-1, 1)) {
    features[[length(features) + 1]] <- dot_curve(
      side * eye_x + eye_width * t,
      eye_y + arc * eye_height * (1 - t^2),
      colour,
      size = 0.30
    )
  }
  features[[length(features) + 1]] <- dot_curve(
    rep(side * eye_x, 5),
    seq(eye_y - 0.043, eye_y + 0.043, length.out = 5),
    colour,
    size = 0.30
  )

  # Three whiskers fan out from each cheek; gaps preserve the scatter texture.
  t <- seq(0, 1, length.out = 22)
  for (slope in c(-0.07, 0.02, 0.13)) {
    whisker_x <- side * (0.20 + 0.61 * t)
    whisker_y <- 0.18 + slope * t + 0.012 * sin(pi * t)
    for (i in seq_along(t)) {
      near_dot <- (rd$x - whisker_x[i])^2 + ((rd$y - whisker_y[i]) * 1.5)^2 < 0.018^2
      rd <- rd[!near_dot, ]
    }
    features[[length(features) + 1]] <- dot_curve(
      whisker_x,
      whisker_y,
      colour,
      size = 0.30,
      alpha = 0.95
    )
  }

  # Dotted ear tips echo the two rising arms of the volcano.
  t <- seq(0, 1, length.out = 24)
  features[[length(features) + 1]] <- dot_curve(
    side * (0.20 + 0.33 * t),
    0.56 + 0.38 * t,
    colour,
    size = 0.26,
    alpha = 0.75,
    jitter = 0.004
  )
  features[[length(features) + 1]] <- dot_curve(
    side * (0.53 + 0.17 * t),
    0.94 - 0.37 * t,
    colour,
    size = 0.26,
    alpha = 0.75,
    jitter = 0.004
  )
}

# Leave a downward triangular nose in the background colour by removing dots.
in_nose <- rd$y >= 0.137 & rd$y <= 0.21 & abs(rd$x) <= 0.045 * (rd$y - 0.137) / (0.21 - 0.137)
rd <- rd[!in_nose, ]
# Mouth curves erase dots with the exact background colour.
t <- seq(0, 1, length.out = 10)
for (side in c(-1, 1)) {
  features[[length(features) + 1]] <- dot_curve(
    side * 0.10 * t,
    0.135 - 0.035 * sin(pi * t),
    fill,
    size = 0.27,
    alpha = 1
  )
}
face <- do.call(rbind, features)

# Keep the highest outlying data points inside the hexagon's sloping shoulders.
sticker_x <- rd$x / 2.10 * 1.44
sticker_y <- -0.25 + (rd$y + 0.03) / 1.06 * 0.92
rd <- rd[sticker_y < 0.90 - abs(sticker_x) / sqrt(3), ]

volcano <- ggplot() +
  geom_point(aes(x, y, colour = colour), rd, size = 0.13, alpha = 0.4, shape = 16) +
  geom_point(aes(x, y, colour = colour, size = size, alpha = alpha), face, shape = 16) +
  scale_colour_identity() +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-1.05, 1.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.03, 1.03), expand = c(0, 0)) +
  theme_void()

hexagon <- function(r) {
  a <- (90 + 60 * 0:5) * pi / 180
  data.frame(x = r * cos(a), y = r * sin(a))
}

grDevices::pdf(NULL)
sticker <- ggplot() +
  geom_polygon(aes(x, y), hexagon(1), fill = border) +
  geom_polygon(aes(x, y), hexagon(0.93), fill = fill) +
  annotation_custom(ggplotGrob(volcano), xmin = -0.72, xmax = 0.72, ymin = -0.25, ymax = 0.67) +
  annotate(
    ggtext::GeomRichText,
    x = 0,
    y = -0.46,
    label = paste0("explo<span style='color:", ns, "'>re</span>DE"),
    colour = "white",
    size = 5.6,
    family = "Avenir Next",
    fontface = "bold",
    fill = NA,
    label.colour = NA,
    label.padding = grid::unit(0, "pt")
  ) +
  coord_fixed(xlim = c(-sqrt(3) / 2, sqrt(3) / 2), ylim = c(-1, 1), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

dir.create(output_dir, showWarnings = FALSE)
ggsave(
  file.path(output_dir, "logo.png"),
  sticker,
  width = sqrt(3),
  height = 2,
  dpi = 600,
  bg = "transparent",
  device = ragg::agg_png
)
grDevices::dev.off()
cat("Wrote", file.path(output_dir, "logo.png"), "\n")
