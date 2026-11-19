# Figure 5 and the Results 2.5 correlations: PC1 against the WAB-R Aphasia Quotient within PWA, against the PHQ-8 and GAD-7
# scores within each group, and AQ against the two mood scores within PWA (Pearson, complete cases). PC scores and clinical
# scores are joined on the full participant ID. Writes figure5_correlations.csv and figure5.png.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

aq <- read_csv(input_file("aphasia_quotient"), show_col_types = FALSE) %>% select(Participant, AQ)
mood <- read_csv(input_file("demographics"), show_col_types = FALSE) %>% select(Participant, PHQ_Total, GAD_Total)
stopifnot(!anyDuplicated(aq$Participant), !anyDuplicated(mood$Participant))
d <- read.csv(output_path("pc_scores.csv")) %>% select(Participant, Group, PC1) %>%
  left_join(mood, by = "Participant") %>% left_join(aq, by = "Participant")
stopifnot(nrow(d) == 100, sum(d$Group == "PWA") == 57, sum(d$Group == "HC") == 43, !anyNA(d$AQ[d$Group == "PWA"]))
pwa <- d %>% filter(Group == "PWA")
hc <- d %>% filter(Group == "HC")

ct <- function(x, y) {
  ok <- complete.cases(x, y)
  t <- cor.test(x[ok], y[ok])
  data.frame(r = unname(t$estimate), p = t$p.value, n = sum(ok))
}
correlations <- bind_rows(cbind(label = "PWA PC1~AQ", ct(pwa$PC1, pwa$AQ)),
                          cbind(label = "PWA PC1~PHQ", ct(pwa$PC1, pwa$PHQ_Total)),
                          cbind(label = "PWA PC1~GAD", ct(pwa$PC1, pwa$GAD_Total)),
                          cbind(label = "HC PC1~PHQ", ct(hc$PC1, hc$PHQ_Total)),
                          cbind(label = "HC PC1~GAD", ct(hc$PC1, hc$GAD_Total)),
                          cbind(label = "PWA AQ~PHQ", ct(pwa$AQ, pwa$PHQ_Total)),
                          cbind(label = "PWA AQ~GAD", ct(pwa$AQ, pwa$GAD_Total)))
write_csv(correlations, output_path("figure5_correlations.csv"))

# Figure 5
group_colors <- c(HC = "#0072B2", PWA = "#CC0000")
annotation <- function(label) {
  x <- correlations[correlations$label == label, ]
  sprintf("r = %.3f, p = %s", x$r, if (x$p < 0.001) "<0.001" else sprintf("%.3f", x$p))
}
fig_theme <- theme_classic(base_size = 15) +
  theme(panel.grid.major = element_line(color = "grey92", linewidth = 0.4), axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(face = "bold", size = 16), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 14, color = "gray40", lineheight = 1.2),
        panel.border = element_rect(color = "black", fill = NA))
panel_a <- ggplot(pwa %>% filter(!is.na(AQ), !is.na(PC1)), aes(AQ, PC1)) +
  geom_point(size = 3, alpha = 0.7, color = group_colors["PWA"]) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2, color = group_colors["PWA"],
              fill = group_colors["PWA"]) +
  labs(title = "A. PC1 vs. Aphasia Quotient", subtitle = annotation("PWA PC1~AQ"), x = "Aphasia Quotient", y = "PC1 Score") +
  fig_theme
mood_panel <- function(score, short, title, x_label) {
  ggplot(d %>% filter(!is.na(.data[[score]]), !is.na(PC1)), aes(.data[[score]], PC1, color = Group)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(fill = Group), alpha = 0.15) +
    scale_color_manual(values = group_colors) + scale_fill_manual(values = group_colors) +
    labs(title = title, x = x_label, y = "PC1 Score", color = "Group:", fill = "Group:",
         subtitle = paste0("HC: ", annotation(paste0("HC PC1~", short)), "\nPWA: ", annotation(paste0("PWA PC1~", short)))) +
    fig_theme
}
figure5 <- panel_a + mood_panel("PHQ_Total", "PHQ", "B. PC1 vs. PHQ Total Score", "PHQ Total Score") +
  mood_panel("GAD_Total", "GAD", "C. PC1 vs. GAD Total Score", "GAD Total Score") +
  plot_layout(guides = "collect", widths = c(1, 1, 1)) & theme(legend.position = "bottom")
ggsave(output_path("figure5.png"), figure5, width = 18, height = 6, dpi = 200, bg = "white")
ggsave(output_path("figure5.pdf"), figure5, width = 18, height = 6, bg = "white", device = grDevices::cairo_pdf)
for (i in seq_len(nrow(correlations))) with(correlations[i, ], cat(sprintf("%-12s r = %6.3f, p = %.3g, n = %d\n", label, r, p, n)))
