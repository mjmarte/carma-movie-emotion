# Figures 6 and 7 (Results 2.6). Figure 6: ROC curves of the three feature sets from the pooled out-of-fold predictions of
# repeat 1 with 95% participant-bootstrap bands, and SHAP values of the combined classifier. Figure 7: predicted versus actual
# Aphasia Quotient from the combined regression with its bootstrap statistics, and SHAP values of the combined regression.
# Reads the outputs of src/17_svm.R (svm_strict/) and src/18_svr.R; writes figure6.png and figure7.png.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggbeeswarm)
  library(patchwork)
  library(pROC)
  library(ggtext)
})
source("R/theme_marte.R"); marte_defaults()

SETCOL <- c(combined = "#400286", language = "#E87722", emotion = "#00A087")
SETLAB <- c(combined = "Combined (7 features)", language = "Language (2 features)", emotion = "Emotion (5 features)")
SETLT <- c(combined = "solid", language = "dashed", emotion = "dotted")

# Figure 6A: pooled out-of-fold ROC curves of repeat 1, with participant-bootstrap 95% bands
oof <- read_csv(output_path("svm_strict", "svm_oof_predictions_rep1.csv"), show_col_types = FALSE) %>% filter(set %in% names(SETCOL))
grid <- seq(0, 1, length.out = 201); set.seed(7)
roc_df <- bind_rows(lapply(names(SETCOL), function(s) {
  d <- oof %>% filter(set == s); y <- as.integer(d$actual == "PWA"); sc <- d$prob
  r <- roc(y, sc, quiet = TRUE); cc <- coords(r, "all", ret = c("specificity", "sensitivity")) %>% arrange(desc(specificity), sensitivity)
  bt <- replicate(1000, {
    ix <- sample(length(y), replace = TRUE); if (length(unique(y[ix])) < 2) return(rep(NA_real_, length(grid)))
    rb <- roc(y[ix], sc[ix], quiet = TRUE); cb <- coords(rb, "all", ret = c("specificity", "sensitivity"))
    approx(1 - cb$specificity, cb$sensitivity, grid, ties = max, rule = 2)$y
  })
  data.frame(set = s, fpr = 1 - cc$specificity, tpr = cc$sensitivity) %>% mutate(band = FALSE) %>%
    bind_rows(data.frame(set = s, fpr = grid, tpr = NA, lo = apply(bt, 1, quantile, .025, na.rm = TRUE),
                         hi = apply(bt, 1, quantile, .975, na.rm = TRUE), band = TRUE))
}))
curves <- roc_df %>% filter(!band); bands <- roc_df %>% filter(band)
p6a <- ggplot() + geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#999999") +
  geom_ribbon(data = bands, aes(x = fpr, ymin = lo, ymax = hi, fill = set), alpha = 0.13) +
  geom_path(data = curves, aes(fpr, tpr, colour = set, linetype = set), linewidth = 0.9) +
  scale_colour_manual(values = SETCOL, labels = SETLAB, name = NULL) + scale_fill_manual(values = SETCOL, labels = SETLAB, name = NULL) +
  scale_linetype_manual(values = SETLT, labels = SETLAB, name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + labs(x = "1 - Specificity", y = "Sensitivity") +
  panel_marte("A", "ROC curves by feature set") +
  theme(legend.position = c(0.68, 0.16), legend.background = element_rect(fill = "white", colour = NA),
        legend.key.width = unit(1.6, "lines"), legend.text = element_text(size = 9))

# Figure 6B: SHAP beeswarm of the combined classifier (repeat 1, all 100 held-out participants)
LAB6 <- c(AntGen_mean = "Linguistic:\nAntonym Generation", PANAS_MSZ_mean = "Emotional Reactivity:\nHolistic MSZ (Mean)",
          CompQ_mean = "Linguistic:\nMovie Comprehension", ISC_mean = "Emotional Reactivity:\nReal-time ISC (Mean)",
          SampEn_Medium = "Emotional Reactivity:\nSample Entropy (Medium)", ISC_sd = "Emotional Reactivity:\nReal-time ISC (SD)",
          CARMA_MSZ_mean = "Emotional Reactivity:\nReal-time MSZ (Mean)")
shap6 <- read_csv(output_path("svm_strict", "svm_shap_combined_rep1.csv"), show_col_types = FALSE) %>% group_by(feature) %>%
  mutate(z = as.numeric(scale(value)), z = pmax(pmin(z, 2.5), -2.5)) %>% ungroup()
ord6 <- shap6 %>% group_by(feature) %>% summarise(m = mean(abs(shap))) %>% arrange(m) %>% pull(feature)
shap6 <- shap6 %>% mutate(flab = factor(LAB6[feature], levels = LAB6[ord6]))
beeswarm_panel <- function(df, xlab, title) {
  ggplot(df, aes(x = shap, y = flab, fill = z)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "#999999", linewidth = 0.4) +
    geom_quasirandom(orientation = "y", shape = 21, colour = CAVIAR_EDGE, stroke = 0.2, size = 1.7, width = 0.35, alpha = 0.9) +
    scale_fill_viridis_c(option = "plasma", name = "Feature value\n(standardized)", limits = c(-2.5, 2.5), breaks = c(-2, 0, 2)) +
    labs(x = xlab, y = NULL) + panel_marte("B", title) +
    theme(legend.position = "right", legend.key.height = unit(1.4, "lines"),
          axis.text.y = element_text(size = 8.5, colour = "#222222", lineheight = 0.9), panel.grid.major.y = element_blank())
}
p6b <- beeswarm_panel(shap6, "SHAP value (impact on model output)", "SHAP values, combined classifier")
f6 <- p6a + p6b + plot_layout(widths = c(1, 1.35))
save_marte(f6, output_path("figure6.png"), width = 12, height = 5.5)
ggsave(output_path("figure6.pdf"), f6, width = 12, height = 5.5, device = grDevices::cairo_pdf)

# Figure 7A: predicted versus actual AQ, combined regression (nested leave-one-out), linear fit with 95% band, identity line
pr <- read_csv(output_path("svr_loo_predictions.csv"), show_col_types = FALSE) %>% filter(set == "combined")
met <- read_csv(output_path("svr_metrics_by_featureset.csv"), show_col_types = FALSE) %>% filter(set == "combined")
lab7 <- sprintf("**r %.3f** [%.3f, %.3f]<br>RMSE %.2f [%.2f, %.2f]<br>MAE %.2f [%.2f, %.2f]", met$r, met$r_lo, met$r_hi, met$RMSE,
                met$RMSE_lo, met$RMSE_hi, met$MAE, met$MAE_lo, met$MAE_hi)
p7a <- ggplot(pr, aes(actual, predicted)) + geom_abline(slope = 1, intercept = 0, colour = "#777777", linewidth = 0.5) +
  regline_marte(colour = "#CC0000") +
  geom_point(shape = 21, fill = "#CC0000", colour = CAVIAR_EDGE, stroke = 0.3, size = 2.4, alpha = 0.85) +
  coord_equal(xlim = c(20, 100), ylim = c(20, 100)) + scale_x_continuous(breaks = seq(20, 100, 20)) +
  scale_y_continuous(breaks = seq(20, 100, 20)) +
  labs(x = "Actual Aphasia Quotient", y = "Predicted Aphasia Quotient") + panel_marte("A", "Predicted versus actual Aphasia Quotient") +
  statbox_marte(lab7, loc = "ul", size = 3.2)

# Figure 7B: SHAP beeswarm of the combined regression (one instance per held-out participant)
LAB7 <- c(AntGen_mean = "Antonym Generation", CompQ_mean = "Movie Comprehension", SampEn_Medium = "Sample Entropy (Medium)",
          ISC_mean = "Real-time ISC (Mean)", PANAS_MSZ_mean = "Holistic MSZ (Mean)", ISC_sd = "Real-time ISC (SD)",
          CARMA_MSZ_mean = "Real-time MSZ (Mean)")
shap7 <- read_csv(output_path("svr_shap_combined.csv"), show_col_types = FALSE) %>% group_by(feature) %>%
  mutate(z = as.numeric(scale(value)), z = pmax(pmin(z, 2.5), -2.5)) %>% ungroup()
ord7 <- shap7 %>% group_by(feature) %>% summarise(m = mean(abs(shap))) %>% arrange(m) %>% pull(feature)
shap7 <- shap7 %>% mutate(flab = factor(LAB7[feature], levels = LAB7[ord7]))
p7b <- beeswarm_panel(shap7, "SHAP value (Aphasia Quotient points)", "SHAP values, combined regression") +
  theme(axis.text.y = element_text(size = 9))
f7 <- p7a + p7b + plot_layout(widths = c(1, 1.25))
save_marte(f7, output_path("figure7.png"), width = 12, height = 6)
ggsave(output_path("figure7.pdf"), f7, width = 12, height = 6, device = grDevices::cairo_pdf)

cat("Figure 6:", nrow(shap6) / 7, "held-out participants in the SHAP panel; feature order (top to bottom):",
    paste(gsub("\n", " ", rev(LAB6[ord6])), collapse = " | "), "\n")
cat("Figure 7: n =", nrow(pr), "; feature order (top to bottom):", paste(rev(LAB7[ord7]), collapse = " | "), "\n")
