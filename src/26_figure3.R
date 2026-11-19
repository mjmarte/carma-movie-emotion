# Figure 3: mean real-time valence rating per second for each clip in PWA and HC, with 95% confidence bands (mean +/- 1.96
# SEM across participants), clips ordered from the most positive to the most negative mean rating. Writes figure3.png.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

ratings <- readRDS(output_path("realtime_ratings.rds"))
group_colors <- c(PWA = "#CC0000", HC = "#0072B2")
group_fill <- c(PWA = "#CC000070", HC = "#0072B270")
mean_ci <- function(x) {
  m <- mean(x, na.rm = TRUE); se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  data.frame(y = m, ymin = m - 1.96 * se, ymax = m + 1.96 * se)
}
figure3 <- ggplot(ratings, aes(Timepoint, Rating, color = Group, fill = Group)) +
  stat_summary(geom = "ribbon", fun.data = mean_ci, alpha = 0.4, color = NA) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 0.5) +
  stat_summary(geom = "line", fun = mean, na.rm = TRUE, linewidth = 0.75) +
  scale_color_manual(values = group_colors) + scale_fill_manual(values = group_fill) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + coord_cartesian(ylim = c(-3, 3)) +
  labs(y = "Valence Rating", x = "Time per clip (secs)", color = "Group:", fill = "Group:") +
  facet_wrap(~ clip_abbrev, ncol = 4, scales = "free_x") +
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
        strip.background = element_rect(fill = "grey90", color = "black"), strip.text = element_text(face = "bold", size = 15),
        axis.text = element_text(color = "black", size = 14), axis.title = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 15), legend.text = element_text(size = 14), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5), panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(size = 13))
ggsave(output_path("figure3.png"), figure3, width = 16, height = 12, dpi = 300, bg = "white")
ggsave(output_path("figure3.pdf"), figure3, width = 16, height = 12, bg = "white", device = grDevices::cairo_pdf)
cat("Figure 3:", n_distinct(ratings$Participant), "participants,", nlevels(ratings$clip_abbrev), "clips\n")
