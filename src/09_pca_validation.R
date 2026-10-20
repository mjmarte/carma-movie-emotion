# Bootstrap stability of the seven-feature PCA reported in Results 2.4 (pca_data_pooled_msz.rds; Methods 4.6): 1,000
# participant resamples (seed 42); for each, the PCA is refit and the variance explained by PC1-PC3, the loadings, and the
# HC-PWA Welch t-test p and Cohen's d on each component are kept. Writes pca_bootstrap_variance.csv, pca_bootstrap_group_differences.csv, pca_bootstrap_cohen_d.csv and
# pca_bootstrap_loadings.csv.
source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

pca_data <- readRDS(output_path("pca_data_pooled_msz.rds"))
features <- c("SampEn_Medium", "PANAS_MSZ_mean", "CARMA_MSZ_mean", "ISC_mean", "ISC_sd", "CompQ_mean", "AntGen_mean")
n_boot <- 1000
n_pcs <- 3
set.seed(42)
variance <- matrix(NA, n_boot, n_pcs, dimnames = list(NULL, paste0("PC", 1:n_pcs)))
p_values <- matrix(NA, n_boot, n_pcs, dimnames = list(NULL, paste0("PC", 1:n_pcs, "_pvalue")))
cohen_d <- matrix(NA, n_boot, n_pcs, dimnames = list(NULL, paste0("PC", 1:n_pcs, "_d")))
loadings <- array(NA, c(length(features), n_pcs, n_boot), dimnames = list(features, paste0("PC", 1:n_pcs), NULL))
for (i in 1:n_boot) {
  tryCatch({
    boot_data <- pca_data[sample(nrow(pca_data), replace = TRUE), ]
    counts <- table(boot_data$Group)
    if (length(counts) < 2 || any(counts < 3)) next
    feature_matrix <- boot_data %>% select(all_of(features)) %>% as.matrix()
    column_var <- apply(feature_matrix, 2, var, na.rm = TRUE)
    if (any(column_var == 0 | is.na(column_var))) next
    boot_pca <- prcomp(scale(feature_matrix))
    variance[i, ] <- summary(boot_pca)$importance[2, 1:n_pcs]
    loadings[, , i] <- boot_pca$rotation[, 1:n_pcs, drop = FALSE]
    scores <- data.frame(Group = boot_data$Group, boot_pca$x[, 1:n_pcs, drop = FALSE])
    for (pc in 1:n_pcs) {
      hc <- scores[[paste0("PC", pc)]][scores$Group == "HC"]
      pwa <- scores[[paste0("PC", pc)]][scores$Group == "PWA"]
      if (length(hc) >= 2 && length(pwa) >= 2) {
        p_values[i, pc] <- t.test(hc, pwa)$p.value
        pooled_sd <- sqrt(((length(hc) - 1) * var(hc) + (length(pwa) - 1) * var(pwa)) / (length(hc) + length(pwa) - 2))
        cohen_d[i, pc] <- (mean(hc) - mean(pwa)) / pooled_sd
      }
    }
  }, error = function(e) NULL)
}

q <- function(x, p) apply(x, 2, quantile, p, na.rm = TRUE)
variance_stats <- data.frame(PC = paste0("PC", 1:n_pcs), Mean = colMeans(variance, na.rm = TRUE), SD = apply(variance, 2, sd, na.rm = TRUE),
                             CI_Lower = q(variance, 0.025), CI_Upper = q(variance, 0.975),
                             CV_Coefficient = apply(variance, 2, sd, na.rm = TRUE) / colMeans(variance, na.rm = TRUE))
difference_stats <- data.frame(PC = paste0("PC", 1:n_pcs), Prop_Significant_p05 = colMeans(p_values < 0.05, na.rm = TRUE),
                               Prop_Significant_p01 = colMeans(p_values < 0.01, na.rm = TRUE),
                               Prop_Significant_p001 = colMeans(p_values < 0.001, na.rm = TRUE),
                               Mean_p_value = colMeans(p_values, na.rm = TRUE), Median_p_value = apply(p_values, 2, median, na.rm = TRUE))
d_stats <- data.frame(PC = paste0("PC", 1:n_pcs), Mean_Cohen_d = colMeans(cohen_d, na.rm = TRUE), SD_Cohen_d = apply(cohen_d, 2, sd, na.rm = TRUE),
                      CI_Lower_d = q(cohen_d, 0.025), CI_Upper_d = q(cohen_d, 0.975),
                      Prop_Large_Effect = colMeans(abs(cohen_d) > 0.8, na.rm = TRUE))
loading_stats <- do.call(rbind, lapply(1:n_pcs, function(pc) {
  l <- loadings[, pc, ]
  data.frame(Feature = features, PC = paste0("PC", pc), Mean = apply(l, 1, mean, na.rm = TRUE), SD = apply(l, 1, sd, na.rm = TRUE),
             CI_Lower = apply(l, 1, quantile, 0.025, na.rm = TRUE), CI_Upper = apply(l, 1, quantile, 0.975, na.rm = TRUE),
             Stability_Score = 1 - apply(l, 1, sd, na.rm = TRUE) / abs(apply(l, 1, mean, na.rm = TRUE)))
}))
rownames(variance_stats) <- rownames(difference_stats) <- rownames(d_stats) <- NULL
write.csv(variance_stats, output_path("pca_bootstrap_variance.csv"), row.names = FALSE)
write.csv(difference_stats, output_path("pca_bootstrap_group_differences.csv"), row.names = FALSE)
write.csv(d_stats, output_path("pca_bootstrap_cohen_d.csv"), row.names = FALSE)
write.csv(loading_stats, output_path("pca_bootstrap_loadings.csv"), row.names = FALSE)
for (k in 1:n_pcs) cat(sprintf("PC%d: variance %.1f%% (95%% CI %.1f%% to %.1f%%, CV %.3f); significant at p < 0.01 in %.1f%% of resamples, mean p %.4f\n",
                               k, 100 * variance_stats$Mean[k], 100 * variance_stats$CI_Lower[k], 100 * variance_stats$CI_Upper[k],
                               variance_stats$CV_Coefficient[k], 100 * difference_stats$Prop_Significant_p01[k], difference_stats$Mean_p_value[k]))
