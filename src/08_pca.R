# Principal component analysis of the seven movie features (Methods 4.6): participant means of medium-scale sample entropy,
# holistic and real-time MSZ, ISC mean and SD, and comprehension and antonym accuracy, standardized; missing feature values take
# the group mean. Main PCA (HC-referenced real-time MSZ): Figure 4 and the PC scores used later; writes pca_data.rds,
# pc_scores.csv, pca_variance.csv, pca_loadings.csv, pc_group_differences.csv (Welch t-tests), pca_results.rds and figure4.png.
# The variance, loadings and group differences reported in Results 2.4 come from the same PCA with real-time MSZ referenced to all
# participants; writes the same files with the suffix _pooled_msz.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})

entropy <- readRDS(output_path("sample_entropy.rds"))
holistic <- readRDS(output_path("holistic_msz.rds"))
isc <- readRDS(output_path("isc.rds"))
tasks <- readRDS(output_path("task_responses.rds"))

build_features <- function(realtime) entropy %>% filter(WindowBin == "Medium") %>% group_by(Participant, Group) %>%
  summarize(SampEn_Medium = mean(SampEn, na.rm = TRUE), .groups = "drop") %>%
  left_join(holistic %>% group_by(Participant, Group) %>% summarize(PANAS_MSZ_mean = mean(msz, na.rm = TRUE), .groups = "drop"),
            by = c("Participant", "Group")) %>%
  left_join(realtime %>% group_by(Participant, Group) %>% summarize(CARMA_MSZ_mean = mean(msz, na.rm = TRUE), .groups = "drop"),
            by = c("Participant", "Group")) %>%
  left_join(isc %>% group_by(Participant, Group) %>%
              summarize(ISC_mean = mean(isc_z, na.rm = TRUE), ISC_sd = sd(isc_z, na.rm = TRUE), .groups = "drop"),
            by = c("Participant", "Group")) %>%
  left_join(tasks %>% group_by(Participant, Group, Task) %>% summarize(accuracy = mean(Accuracy, na.rm = TRUE), .groups = "drop") %>%
              pivot_wider(names_from = Task, values_from = accuracy, names_glue = "{Task}_mean"),
            by = c("Participant", "Group")) %>%
  mutate(outcome = factor(Group, levels = c("HC", "PWA"))) %>%
  group_by(Group) %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>% ungroup() %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

features <- c("SampEn_Medium", "PANAS_MSZ_mean", "CARMA_MSZ_mean", "ISC_mean", "ISC_sd", "CompQ_mean", "AntGen_mean")
labels <- c(SampEn_Medium = "Sample Entropy", PANAS_MSZ_mean = "Holistic MSZ", CARMA_MSZ_mean = "Real-Time MSZ",
            ISC_mean = "Real-time ISC (Mean)", ISC_sd = "Real-time ISC (SD)", CompQ_mean = "Movie Comprehension",
            AntGen_mean = "Antonym Generation")
run_pca <- function(pca_data, suffix) {
  saveRDS(pca_data, output_path(paste0("pca_data", suffix, ".rds")))
  pca <- prcomp(scale(as.matrix(pca_data[, features])))
  importance <- summary(pca)$importance
  eigenvalues <- pca$sdev^2
  n_pc <- length(eigenvalues)
  variance <- data.frame(PC = paste0("PC", 1:n_pc), Eigenvalue = eigenvalues, Variance_Pct = importance[2, ] * 100,
                         Cumulative_Variance_Pct = importance[3, ] * 100)
  loadings <- data.frame(Feature = rownames(pca$rotation), pca$rotation) %>%
    mutate(Domain = ifelse(grepl("CompQ|AntGen", Feature), "Linguistic", "Emotional Reactivity"),
           Feature_Clean = unname(labels[Feature]))
  for (i in 1:n_pc) loadings[[paste0("PC", i, "_sq")]] <- loadings[[paste0("PC", i)]]^2
  scores <- data.frame(Participant = pca_data$Participant, Group = pca_data$Group, pca$x)
  differences <- map_dfr(1:n_pc, function(k) {
    x <- scores[[paste0("PC", k)]]; hc <- x[scores$Group == "HC"]; pwa <- x[scores$Group == "PWA"]
    t_result <- t.test(hc, pwa)
    pooled_sd <- sqrt(((length(hc) - 1) * sd(hc)^2 + (length(pwa) - 1) * sd(pwa)^2) / (length(hc) + length(pwa) - 2))
    data.frame(PC = paste0("PC", k), HC_Mean = mean(hc), HC_SD = sd(hc), PWA_Mean = mean(pwa), PWA_SD = sd(pwa),
               Difference = mean(hc) - mean(pwa), SE = unname(t_result$stderr),
               t_value = unname(t_result$statistic), df = unname(t_result$parameter), p_value = t_result$p.value,
               Cohen_d = (mean(hc) - mean(pwa)) / pooled_sd)
  })
  write.csv(variance, output_path(paste0("pca_variance", suffix, ".csv")), row.names = FALSE)
  write.csv(loadings, output_path(paste0("pca_loadings", suffix, ".csv")), row.names = FALSE)
  write.csv(scores, output_path(paste0("pc_scores", suffix, ".csv")), row.names = FALSE)
  write.csv(differences, output_path(paste0("pc_group_differences", suffix, ".csv")), row.names = FALSE)
  result <- list(pca_result = pca, pc_scores = scores, loadings = loadings, variance_explained = variance,
                 group_differences = differences, eigenvalues = eigenvalues)
  saveRDS(result, output_path(paste0("pca_results", suffix, ".rds")))
  cat(sprintf("PCA%s: variance PC1-PC3 %.2f / %.2f / %.2f%% (cumulative %.2f%%); eigenvalues >= 1: %d\n", suffix,
              variance$Variance_Pct[1], variance$Variance_Pct[2], variance$Variance_Pct[3], variance$Cumulative_Variance_Pct[3],
              sum(eigenvalues >= 1)))
  result
}
pca_data <- build_features(readRDS(output_path("realtime_msz.rds")))
stopifnot(nrow(pca_data) == 100)
main <- run_pca(pca_data, "")
invisible(run_pca(build_features(readRDS(output_path("realtime_msz_pooled.rds"))), "_pooled_msz"))
pca <- main$pca_result; scores <- main$pc_scores; loadings <- main$loadings; eigenvalues <- main$eigenvalues
importance <- summary(pca)$importance

# Figure 4
group_colors <- c(HC = "#0072B2", PWA = "#CC0000")
centroids <- scores %>% group_by(Group) %>% summarise(across(c(PC1, PC2, PC3), mean), .groups = "drop")
pc_theme <- theme_bw(base_size = 10) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 10), axis.text = element_text(size = 9))
pc_panel <- function(x, y, title) {
  ggplot(scores, aes(.data[[x]], .data[[y]])) +
    stat_ellipse(aes(color = Group), level = 0.68, linewidth = 1, alpha = 0.5) +
    stat_ellipse(aes(color = Group), level = 0.95, linewidth = 0.8, linetype = "dashed", alpha = 0.3) +
    geom_point(data = centroids, aes(fill = Group), shape = 23, size = 3, color = "black", stroke = 1) +
    geom_point(aes(color = Group), size = 1.5, alpha = 0.6) +
    scale_color_manual(values = group_colors, name = "Group") + scale_fill_manual(values = group_colors, name = "Group") +
    labs(title = title, x = sprintf("%s (%.1f%%)", x, importance[2, x] * 100), y = sprintf("%s (%.1f%%)", y, importance[2, y] * 100)) +
    pc_theme + coord_equal()
}
panel_c <- pc_panel("PC2", "PC3", "C: PC2 vs PC3") +
  theme(legend.position = "right", legend.title = element_text(size = 14), legend.text = element_text(size = 12))
n_show <- min(sum(eigenvalues > 1), 3)
loading_rows <- map_dfr(paste0("PC", 1:n_show), ~ loadings %>% select(Feature_Clean, Domain, Loading = all_of(.x)) %>% mutate(PC = .x)) %>%
  mutate(PC = factor(PC, levels = paste0("PC", 1:n_show)))
panel_d <- ggplot(loading_rows, aes(reorder(Feature_Clean, Loading), Loading, fill = Domain)) +
  theme_bw(base_size = 10) + geom_col(alpha = 0.8) + geom_hline(yintercept = 0, alpha = 0.5) +
  facet_wrap(~ PC, scales = "free_x", nrow = 1) + coord_flip() +
  scale_fill_manual(values = c(Linguistic = "#E69F00", `Emotional Reactivity` = "#56B4E9"), name = "Feature Domain") +
  labs(title = "D: Feature Loadings", x = NULL, y = "Loading Value") +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 10), strip.background = element_rect(fill = "gray90"),
        axis.text.y = element_text(size = 12, colour = "black"), axis.text.x = element_text(size = 9),
        axis.title.x = element_text(face = "bold", size = 10))
figure4 <- (pc_panel("PC1", "PC2", "A: PC1 vs PC2") | pc_panel("PC1", "PC3", "B: PC1 vs PC3") | panel_c) / panel_d
ggsave(output_path("figure4.png"), figure4, width = 16, height = 10, dpi = 300, bg = "white")
