# theme_marte — Marte-lab ggplot2 figure aesthetic (theme_minimal base).
# source("theme_marte.R"); marte_defaults(); ggplot(...) + theme_marte()
#
# Raincloud across groups — box + jittered caviar 'spread' on the LEFT, thin-dark-
# bordered half-violin 'cloud' on the RIGHT. Version-robust: uses gghalves for the
# cloud only (ggdist::stat_halfeye and gghalves::geom_half_point are BROKEN on
# ggplot2 >= 4.0), and precomputed jitter for the points:
#   df$xnum <- as.integer(df$g)
#   df$xj   <- df$xnum - 0.20 + runif(nrow(df), -0.055, 0.055)
#   ggplot(df, aes(xnum, y)) +
#     gghalves::geom_half_violin(aes(fill=g), side="r", colour="#1a1a1a",
#       linewidth=.4, alpha=.55, position=position_nudge(x=.05), width=.9) +
#     geom_boxplot(aes(group=g), width=.10, outlier.shape=NA, alpha=.6, fill="white",
#       position=position_nudge(x=-.03)) +
#     geom_point(aes(x=xj, fill=g), shape=21, colour="#1a1a1a", stroke=.3, size=1.9) +
#     scale_fill_manual(values=group_colors) +
#     scale_x_continuous(breaks=seq_along(levels(df$g)), labels=levels(df$g))
suppressPackageStartupMessages({ library(ggplot2) })

# ── accent palette (lab brand) ────────────────────────────────────
PURPLE  <- "#400286"; PURPLE2 <- "#7b3fa8"; ORANGE <- "#E87722"
BLUE    <- "#4393C3"; GREEN   <- "#2E7D32"; RED    <- "#C62828"
CORAL   <- "#D6604D"; LGREEN  <- "#74C476"; LORANGE<- "#FD8D3C"; GREY <- "#999999"
PURPLE_RAMP  <- c("#C6ACDF","#9970AB","#6A3D9A","#400286")
model_colors <- c(ridge="#4393C3", svr="#D6604D", rf="#74C476", xgb="#FD8D3C")
posneg_colors<- c("1"="#2E7D32", "0"="#C62828")
# ── ggsci NPG (Nature) categorical palette; groups for this project ──
NPG <- c("#E64B35","#3C5488","#00A087","#4DBBD5","#F39B7F","#8491B4","#91D1C2")
# NPG teal/navy/red — keyed for both short (L/R) and full (LH/RH) group labels
group_colors <- c(HC="#00A087", L="#3C5488", R="#E64B35",
                  LH="#3C5488", RH="#E64B35")
# (equivalently: ggsci::scale_colour_npg() / scale_fill_npg())

theme_marte <- function() {
  # Tuned for pixel-parity with mstyle.py (Python). Keep the two in lockstep:
  # keep in step with the matplotlib version of this theme
  theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      axis.ticks           = element_blank(),
      panel.grid.minor     = element_blank(),
      panel.grid.major     = element_line(color = "#EBEBEB", linewidth = 0.25),
      strip.text           = element_text(face = "bold", size = 10),
      legend.position      = "bottom",
      legend.title         = element_text(size = 8.5),
      legend.text          = element_text(size = 8.5),
      axis.title           = element_text(size = 10, colour = "#222222"),
      axis.text            = element_text(size = 9,  colour = "#666666"),
      plot.title           = element_text(face = "bold", size = 12, colour = "#222222",
                                          hjust = 0, margin = margin(b = 6)),
      plot.subtitle        = element_text(size = 9, color = "#666666", hjust = 0),
      # title/caption align above the PANEL (matches matplotlib loc="left" for parity)
      plot.title.position   = "panel",
      plot.caption.position = "panel"
    )
}

CAVIAR_EDGE <- "#1a1a1a"   # thin dark point outline (the caviar look)

# Caviar the HOUSE DEFAULT: call once after sourcing so every geom_point/jitter
# gets a thin dark outline (shape 21, fill = aesthetic). From ohmyggplot.
marte_defaults <- function() {
  # Mirrors mstyle.apply(): sets the house style GLOBALLY so every plot inherits it
  # (theme + caviar geom defaults) — no need to add + theme_marte() to each plot.
  ggplot2::theme_set(theme_marte())
  # size 2.5 / stroke 0.3 ≈ matplotlib scatter s=42 / patch.linewidth 0.8 (parity)
  ggplot2::update_geom_defaults("point",
    list(shape = 21, colour = CAVIAR_EDGE, stroke = 0.3, size = 2.5, fill = "gray35"))
  ggplot2::update_geom_defaults("jitter",
    list(shape = 21, colour = CAVIAR_EDGE, stroke = 0.3, size = 2.2))
  invisible(TRUE)
}

# Explicit bordered 'caviar' point (from ohmyggplot::geom_caviar). Use when you
# want to force the look on one layer regardless of geom defaults. Map fill for color.
geom_caviar <- function(..., size = 2.5, stroke = 0.3) {
  ggplot2::geom_point(..., shape = 21, colour = CAVIAR_EDGE, size = size, stroke = stroke)
}

# House stat strip (parity with ms.statbox): small, detailed, BORDERED corner label.
# `label` is markdown for ggtext — bold the key value, <br> for line breaks, e.g.
#   "**AUC 0.81**<br>95% CI [0.75, 0.87]<br>n = 200"
# loc in {ul, ur, ll, lr}. Requires ggtext.
statbox_marte <- function(label, loc = "lr", size = 2.7) {
  p <- switch(loc,
    ur = c( Inf,  Inf, 1.08,  1.35),  lr = c( Inf, -Inf, 1.08, -0.35),
    ul = c(-Inf,  Inf, -0.08, 1.35),  ll = c(-Inf, -Inf, -0.08, -0.35))
  ggtext::geom_richtext(
    data = data.frame(x = p[1], y = p[2], lab = label),
    aes(x = x, y = y, label = lab), hjust = p[3], vjust = p[4],
    size = size, family = "Helvetica", lineheight = 1.3,
    text.colour = "#222222", fill = "white",
    label.colour = "#222222", label.size = 0.25,          # ~0.7pt border
    label.r = unit(0.18, "lines"),
    label.padding = unit(c(.18, .3, .18, .3), "lines"),
    inherit.aes = FALSE)
}

# Opt-in Tufte range-frame ('floating axes'). Add ggthemes::geom_rangeframe() as a
# layer for the trimmed spines; this theme strips grid and unhides the axis lines.
theme_marte_floating <- function() {
  theme_marte() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "gray40", linewidth = 0.5)
    )
  # usage: p + theme_marte_floating() + ggthemes::geom_rangeframe(colour = "gray40")
}

# Horizontal colorbar on top for a continuous fill/colour (from ohmyggplot).
better_fill_legend  <- ggplot2::guides(fill  = ggplot2::guide_colorbar(
  title.position = "top", title.hjust = 0.5,
  barwidth = grid::unit(12, "lines"), barheight = grid::unit(0.4, "lines")))
better_color_legend <- ggplot2::guides(colour = ggplot2::guide_colorbar(
  title.position = "top", title.hjust = 0.5,
  barwidth = grid::unit(12, "lines"), barheight = grid::unit(0.4, "lines")))

# ── API-parity helpers: mirror mstyle.py so R and Python are called the same way ──
# Python ms.top_colorbar ↔ R better_fill_legend / better_color_legend (above).

# Python ms.panel(ax,"A","Title") ↔ panel_marte("A","Title")  (adds a labs title)
panel_marte <- function(letter, text = "")
  ggplot2::labs(title = if (nzchar(letter)) paste0(letter, ". ", text) else text)

# Python ms.subtitle ↔ subtitle_marte("...")
subtitle_marte <- function(text) ggplot2::labs(subtitle = text)

# Python ms.emph(ax,x,y,"0.97") ↔ emph_marte(x,y,"0.97") — the ONE bold answer
emph_marte <- function(x, y, label, colour = "#222222", size = 4)
  ggplot2::annotate("text", x = x, y = y, label = label, fontface = "bold",
                    colour = colour, size = size, family = "Helvetica")

# Python ms.ci_band(ax,x,lo,hi) ↔ ci_band_marte() — a styled ribbon layer
ci_band_marte <- function(mapping = NULL, data = NULL, colour = BLUE, alpha = 0.18, ...)
  ggplot2::geom_ribbon(mapping = mapping, data = data, fill = colour, alpha = alpha,
                       colour = NA, ...)

# Python ms.regline_ci(ax,x,y) ↔ regline_marte() — OLS fit + mean-prediction CI band (layer)
regline_marte <- function(colour = BLUE, ci = 0.95, ...)
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = ci,
                       colour = colour, fill = colour, alpha = 0.18, linewidth = 0.6, ...)

# Python ms.roc_ci(ax,y,s) ↔ roc_marte(y,s) — returns a complete ggplot; AUC/CI in attributes
.roc_curve <- function(y, s) {
  o <- order(-s); y <- y[o]
  P <- sum(y == 1); N <- sum(y == 0)
  tpr <- c(0, cumsum(y == 1) / P); fpr <- c(0, cumsum(y == 0) / N)
  list(fpr = fpr, tpr = tpr,
       auc = sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2))
}
roc_marte <- function(y, s, n_boot = 1000, colour = PURPLE) {
  y <- as.integer(y); s <- as.numeric(s); grid <- seq(0, 1, length.out = 101)
  rc <- .roc_curve(y, s)
  bt <- matrix(NA_real_, n_boot, length(grid)); aucs <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(length(y), replace = TRUE)
    if (length(unique(y[idx])) < 2) next
    r <- .roc_curve(y[idx], s[idx])
    bt[b, ] <- stats::approx(r$fpr, r$tpr, grid, ties = max, rule = 2)$y
    aucs[b] <- r$auc
  }
  band <- data.frame(grid = grid,
                     lo = apply(bt, 2, quantile, .025, na.rm = TRUE),
                     hi = apply(bt, 2, quantile, .975, na.rm = TRUE))
  lo <- unname(quantile(aucs, .025, na.rm = TRUE)); hi <- unname(quantile(aucs, .975, na.rm = TRUE))
  p <- ggplot2::ggplot() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         linewidth = 0.3, colour = "#999999") +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(grid, ymin = lo, ymax = hi),
                         fill = colour, alpha = 0.2) +
    ggplot2::geom_line(data = data.frame(fpr = rc$fpr, tpr = rc$tpr),
                       ggplot2::aes(fpr, tpr), colour = colour, linewidth = 0.7) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "False positive rate", y = "True positive rate") + theme_marte()
  attr(p, "auc") <- rc$auc; attr(p, "ci") <- c(lo, hi)
  message(sprintf("roc_marte: AUC %.3f [%.3f, %.3f]", rc$auc, lo, hi))
  p
}

# Python ms.raincloud(ax,...) ↔ raincloud_marte(df, x, y) — returns a complete ggplot.
# box + jittered caviar 'spread' left, thin-dark-bordered half-violin 'cloud' right.
# Version-robust (gghalves half-violin only; precomputed jitter for points).
# `interval` adds the inferential layer (house 'always plot CIs' rule) in the lane
# between box and cloud: "mean" (filled diamond + 95% bootstrap CI, DEFAULT) /
# "median" (open circle + 95% bootstrap CI, robust for skew) / "none".
raincloud_marte <- function(df, x, y, colors = NULL, dist_border = CAVIAR_EDGE,
                            interval = "mean", n_boot = 2000) {
  df <- as.data.frame(df); df$.g <- factor(df[[x]]); df$.y <- df[[y]]
  df$.xn <- as.integer(df$.g)
  df$.xj <- df$.xn - 0.20 + stats::runif(nrow(df), -0.055, 0.055)
  p <- ggplot2::ggplot(df, ggplot2::aes(.xn, .y)) +
    gghalves::geom_half_violin(ggplot2::aes(fill = .g), side = "r", colour = dist_border,
      linewidth = 0.4, alpha = 0.55, position = ggplot2::position_nudge(x = 0.05), width = 0.9) +
    ggplot2::geom_boxplot(ggplot2::aes(group = .g), width = 0.10, outlier.shape = NA,
      alpha = 0.6, fill = "white", position = ggplot2::position_nudge(x = -0.03)) +
    ggplot2::geom_point(ggplot2::aes(.xj, fill = .g), shape = 21, colour = dist_border,
      stroke = 0.3, size = 2, alpha = 0.9) +
    ggplot2::scale_x_continuous(breaks = seq_along(levels(df$.g)), labels = levels(df$.g)) +
    theme_marte()
  if (interval %in% c("mean", "median")) {
    fstat <- if (interval == "mean") mean else stats::median
    est <- do.call(rbind, lapply(levels(df$.g), function(g) {
      d <- df$.y[df$.g == g]
      bv <- vapply(seq_len(n_boot), function(i) fstat(sample(d, length(d), replace = TRUE)), numeric(1))
      data.frame(.xn = as.integer(factor(g, levels = levels(df$.g))) + 0.02,
                 m = fstat(d), lo = quantile(bv, .025), hi = quantile(bv, .975), .g = g)
    }))
    p <- p +
      ggplot2::geom_linerange(data = est, ggplot2::aes(x = .xn, ymin = lo, ymax = hi),
        colour = "#222222", linewidth = 0.6, inherit.aes = FALSE) +
      (if (interval == "mean")
        ggplot2::geom_point(data = est, ggplot2::aes(.xn, m, fill = .g), shape = 23,
          colour = "#222222", stroke = 0.4, size = 3, inherit.aes = FALSE)
       else
        ggplot2::geom_point(data = est, ggplot2::aes(.xn, m), shape = 21, fill = "white",
          colour = "#222222", stroke = 0.6, size = 2.6, inherit.aes = FALSE))
  }
  if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
  p
}

# Python ms.save(fig,path,tif=) ↔ save_marte(plot,path,...) — 300 dpi PNG (+ LZW TIFF)
save_marte <- function(plot, path, width = 6, height = 4, tif = FALSE) {
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 300, bg = "white")
  if (tif) ggplot2::ggsave(sub("\\.[^.]+$", ".tif", path), plot, width = width,
                           height = height, dpi = 300, bg = "white",
                           device = grDevices::tiff, compression = "lzw")
  invisible(path)
}
