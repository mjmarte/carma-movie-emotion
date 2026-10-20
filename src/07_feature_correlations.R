# Correlations among the 13 candidate movie features (entropy at three scales; holistic and real-time MSZ, ISC and task
# accuracy, each as participant mean and SD), used to select the seven PCA features (Supplementary text, Figure S2).
# Real-time MSZ referenced to all participants, as in the PCA of Results 2.4. Missing feature values take the group mean.
# Writes feature_correlations.csv and figure_s2.pdf (the published figure is this PDF exported at 200 dpi).
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(corrplot)
})

entropy <- readRDS(output_path("sample_entropy.rds"))
holistic <- readRDS(output_path("holistic_msz.rds"))
realtime <- readRDS(output_path("realtime_msz_pooled.rds"))
isc <- readRDS(output_path("isc.rds"))
tasks <- readRDS(output_path("task_responses.rds"))

features <- entropy %>% group_by(Participant, Group, WindowBin) %>%
  summarize(mean_SampEn = mean(SampEn, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = WindowBin, values_from = mean_SampEn, names_prefix = "SampEn_") %>%
  left_join(holistic %>% group_by(Participant, Group) %>%
              summarize(PANAS_MSZ_mean = mean(msz, na.rm = TRUE), PANAS_MSZ_sd = sd(msz, na.rm = TRUE), .groups = "drop"),
            by = c("Participant", "Group")) %>%
  left_join(realtime %>% group_by(Participant, Group) %>%
              summarize(CARMA_MSZ_mean = mean(msz, na.rm = TRUE), CARMA_MSZ_sd = sd(msz, na.rm = TRUE), .groups = "drop"),
            by = c("Participant", "Group")) %>%
  left_join(isc %>% group_by(Participant, Group) %>%
              summarize(ISC_mean = mean(isc_z, na.rm = TRUE), ISC_sd = sd(isc_z, na.rm = TRUE), .groups = "drop"),
            by = c("Participant", "Group")) %>%
  left_join(tasks %>% group_by(Participant, Group, Task) %>%
              summarize(accuracy_mean = mean(Accuracy, na.rm = TRUE), accuracy_sd = sd(Accuracy, na.rm = TRUE), .groups = "drop") %>%
              pivot_wider(names_from = Task, values_from = c(accuracy_mean, accuracy_sd), names_glue = "{Task}_{.value}"),
            by = c("Participant", "Group")) %>%
  group_by(Group) %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>% ungroup()
feature_cols <- names(features)[sapply(features, is.numeric)]
cor_matrix <- cor(features[feature_cols], use = "complete.obs")
write.csv(cor_matrix, output_path("feature_correlations.csv"))

pdf(output_path("figure_s2.pdf"), width = 12, height = 10)
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", tl.cex = 0.7, tl.col = "black",
         addCoef.col = "black", number.cex = 0.5, diag = FALSE, title = "Feature Correlation Matrix", mar = c(0, 0, 2, 0))
invisible(dev.off())

entropy_r <- cor_matrix[grep("SampEn", feature_cols, value = TRUE), grep("SampEn", feature_cols, value = TRUE)]
cat(sprintf(paste0("%d features, %d participants. Entropy scales: min r = %.3f; holistic MSZ mean vs SD r = %.3f; ",
                   "real-time MSZ mean vs SD r = %.3f; comprehension accuracy mean vs SD r = %.3f\n"),
            length(feature_cols), nrow(features), min(entropy_r[upper.tri(entropy_r)]), cor_matrix["PANAS_MSZ_mean", "PANAS_MSZ_sd"],
            cor_matrix["CARMA_MSZ_mean", "CARMA_MSZ_sd"], cor_matrix["CompQ_accuracy_mean", "CompQ_accuracy_sd"]))
