# Supplementary Figure S1: PC1 versus main concept distance of the spoken clip summaries (mean cosine distance of a
# participant's utterances to the centroid of the clip's main concept inventory, z-scored within clip over the analytic sample
# and averaged over clips). Missing distances are multiply imputed within group (predictive mean matching, 50 imputations) and
# correlations pooled on the Fisher z scale with Rubin's rules, including the PWA partial correlation controlling for AQ.
# Writes figure_s1.png, figure_s1_correlations.csv and figure_s1_imputed.csv (plotted values of the imputed participants).
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(mice)
})

ids <- read.csv(output_path("analytic_ids.csv"))
# rows stay in the order of the PC scores: the seeded imputation depends on row order
pcs <- read.csv(output_path("pc_scores.csv")) %>%
  select(Participant, PC1) %>%
  inner_join(ids, by = "Participant")
mc <- read_csv(input_file("main_concept_measures"), show_col_types = FALSE) %>%
  filter(data_type == "summary") %>%
  rename(Participant = participant_id, Clip = clip) %>%
  distinct(Participant, Clip, .keep_all = TRUE) %>%
  filter(Participant %in% ids$Participant) %>%
  group_by(Clip) %>%
  mutate(z = (avg_distance_to_centroid - mean(avg_distance_to_centroid)) / sd(avg_distance_to_centroid)) %>%
  group_by(Participant) %>%
  summarise(MC = mean(z), .groups = "drop")
aq <- read_csv(input_file("aphasia_quotient"), show_col_types = FALSE) %>% select(Participant, AQ)
stopifnot(!anyDuplicated(aq$Participant))
full <- pcs %>%
  left_join(mc, by = "Participant") %>%
  left_join(aq, by = "Participant") %>%
  mutate(Group = factor(Group, levels = c("HC", "PWA")))
stopifnot(nrow(full) == 100, sum(full$Group == "PWA") == 57, !anyNA(full$AQ[full$Group == "PWA"]))

# Imputation within group: PWA model PC1 + MC + AQ, HC model PC1 + MC. Pooled within variance 1/(n - 3), or 1/(n - 4) for the
# partial correlation.
M <- 50
fp <- full %>% filter(Group == "PWA")
fh <- full %>% filter(Group == "HC")
imp_p <- mice(as.data.frame(fp %>% select(PC1, MC, AQ)), m = M, method = "pmm", seed = 2026, printFlag = FALSE)
imp_h <- mice(as.data.frame(fh %>% select(PC1, MC)), m = M, method = "pmm", seed = 2026, printFlag = FALSE)
pool_r <- function(rs, n, k = 0) {
  z <- atanh(rs); W <- 1 / (n - 3 - k); B <- var(z); Tv <- W + (1 + 1 / length(z)) * B; zb <- mean(z)
  c(r = tanh(zb), lo = tanh(zb - qnorm(.975) * sqrt(Tv)), hi = tanh(zb + qnorm(.975) * sqrt(Tv)),
    p = 2 * pnorm(-abs(zb) / sqrt(Tv)), n = n)
}
partial_r <- function(y, x, z) cor(resid(lm(y ~ z)), resid(lm(x ~ z)))
cp <- lapply(1:M, function(i) complete(imp_p, i))
ch <- lapply(1:M, function(i) complete(imp_h, i))
mi_res <- rbind(
  PWA = pool_r(sapply(cp, function(x) cor(x$PC1, x$MC)), nrow(fp)),
  HC = pool_r(sapply(ch, function(x) cor(x$PC1, x$MC)), nrow(fh)),
  PWA_partial_AQ = pool_r(sapply(cp, function(x) partial_r(x$PC1, x$MC, x$AQ)), nrow(fp), k = 1),
  PWA_MC_AQ = pool_r(sapply(cp, function(x) cor(x$MC, x$AQ)), nrow(fp)),
  PWA_PC1_AQ = pool_r(sapply(cp, function(x) cor(x$PC1, x$AQ)), nrow(fp)),
  All = pool_r(sapply(1:M, function(i) cor(c(cp[[i]]$PC1, ch[[i]]$PC1), c(cp[[i]]$MC, ch[[i]]$MC))), nrow(full)))
write.csv(mi_res, output_path("figure_s1_correlations.csv"))
# An imputed participant is plotted at the mean of its 50 imputed distances
fp$MC_plot <- fp$MC
fp$MC_plot[is.na(fp$MC)] <- rowMeans(as.matrix(imp_p$imp$MC))
fh$MC_plot <- fh$MC
fh$MC_plot[is.na(fh$MC)] <- rowMeans(as.matrix(imp_h$imp$MC))
dplot <- bind_rows(fp, fh) %>% mutate(imputed = is.na(MC))
write.csv(dplot %>% filter(imputed) %>% select(Participant, Group, PC1, MC_plot), output_path("figure_s1_imputed.csv"),
          row.names = FALSE)

# Figure S2, in the style of Figure 5; r and p in the subtitle are the pooled values
group_colors <- c(HC = "#0072B2", PWA = "#CC0000")
manuscript_theme <- function() {
  theme_classic(base_size = 15) +
    theme(
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black", size = 14),
      axis.title = element_text(face = "bold", size = 16),
      legend.title = element_text(face = "plain", size = 15),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 14, color = "gray40", lineheight = 1.2),
      panel.border = element_rect(color = "black", fill = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}
fmt_r <- function(x) sprintf("%.3f", x)
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
subtitle <- sprintf("HC: r = %s, p = %s\nPWA: r = %s, p = %s", fmt_r(mi_res["HC", "r"]), fmt_p(mi_res["HC", "p"]),
                    fmt_r(mi_res["PWA", "r"]), fmt_p(mi_res["PWA", "p"]))
update_geom_defaults("point", list(shape = 19, colour = "black", fill = NA, size = 1.5, stroke = 0.5))
figure_s1 <- ggplot(dplot, aes(x = MC_plot, y = PC1, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(fill = Group), alpha = 0.15) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(title = "PC1 vs. Main Concept Distance", subtitle = subtitle, x = "Main Concept Distance (z-score)", y = "PC1 Score",
       color = "Group:", fill = "Group:") +
  manuscript_theme()
ggsave(output_path("figure_s1.png"), figure_s1, width = 6, height = 6, dpi = 300, bg = "white")

cat(sprintf("main concept distance imputed for %d PWA and %d HC (%d participants with transcribed summaries)\n",
            sum(is.na(fp$MC)), sum(is.na(fh$MC)), sum(!is.na(full$MC))))
for (k in rownames(mi_res)) {
  cat(sprintf("%-15s r = %.3f, 95%% CI [%.3f, %.3f], p = %.3g, n = %d\n", k, mi_res[k, "r"], mi_res[k, "lo"], mi_res[k, "hi"],
              mi_res[k, "p"], as.integer(mi_res[k, "n"])))
}
