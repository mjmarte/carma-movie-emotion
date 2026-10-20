# Prediction of the Aphasia Quotient within PWA by radial-kernel support vector regression (Results 2.6, Figure 7;
# Methods 4.8): nested leave-one-out cross-validation with 5-fold inner tuning, the holistic MSZ reference recomputed without the
# held-out participant, four feature sets, 2,000 participant bootstrap resamples of the predictions, SHAP for the combined model.
# Writes svr_loo_predictions.csv, svr_metrics_by_featureset.csv, svr_featureset_difference.csv, svr_hyperparameters.csv,
# svr_shap_combined.csv and svr_shap_summary.csv.
source("R/config.R")
suppressPackageStartupMessages({
  library(e1071)
  library(kernelshap)
})
source("R/foldwise_features.R")

tune_fit <- function(train, feats, seed) {
  set.seed(seed)
  tn <- tune(svm, as.formula(paste("AQ ~", paste(feats, collapse = " + "))), data = train,
             ranges = list(cost = c(0.1, 1, 10), epsilon = c(0.01, 0.1)), tunecontrol = tune.control(cross = 5))
  tn$best.model
}
preds <- list(); shap_rows <- list(); hp <- list(); t0 <- Sys.time()
for (i in seq_along(PWA_IDS)) {
  held <- PWA_IDS[i]; train_ids <- setdiff(ALL_IDS, held)
  f <- compute_features(train_ids) %>% inner_join(AQ_CORRECT, by = "Participant") %>% filter(Group == "PWA")
  stopifnot(nrow(f) == 57, !anyNA(f[, FEATS_ALL]))
  train <- f[f$Participant != held, ]; test <- f[f$Participant == held, ]
  for (s in names(FEATURE_SETS)) {
    m <- tune_fit(train, FEATURE_SETS[[s]], 123 + i)
    preds[[length(preds) + 1]] <- data.frame(set = s, Participant = held, actual = test$AQ, predicted = as.numeric(predict(m, test)))
    hp[[length(hp) + 1]] <- data.frame(set = s, fold = i, cost = m$cost, epsilon = m$epsilon)
    if (s == "combined") {
      set.seed(123 + i); bg <- train[sample(nrow(train), 10), FEATS_ALL]
      ks <- kernelshap(m, X = test[, FEATS_ALL], bg_X = bg, pred_fun = function(mod, nd) predict(mod, as.data.frame(nd)),
                       verbose = FALSE)
      shap_rows[[length(shap_rows) + 1]] <- data.frame(Participant = held, feature = FEATS_ALL, shap = as.numeric(ks$S[1, ]),
                                                       value = as.numeric(test[1, FEATS_ALL]))
    }
  }
  if (i %% 10 == 0) cat("leave-one-out fold", i, "of", length(PWA_IDS), "\n")
}
preds <- bind_rows(preds); shap <- bind_rows(shap_rows); hp <- bind_rows(hp)

metrics <- function(a, p) c(RMSE = sqrt(mean((a - p)^2)), MAE = mean(abs(a - p)), r = cor(a, p),
                            R2 = 1 - sum((a - p)^2) / sum((a - mean(a))^2))
set.seed(42); B <- 2000
wide <- preds %>% select(set, Participant, actual, predicted) %>% pivot_wider(names_from = set, values_from = predicted)
boot <- replicate(B, {
  ix <- sample(nrow(wide), replace = TRUE); w <- wide[ix, ]
  c(unlist(lapply(names(FEATURE_SETS), function(s) setNames(metrics(w$actual, w[[s]]), paste0(s, ".", c("RMSE", "MAE", "r", "R2"))))),
    d_r = cor(w$actual, w$combined) - cor(w$actual, w$language),
    d_RMSE = sqrt(mean((w$actual - w$combined)^2)) - sqrt(mean((w$actual - w$language)^2)))
})
ci <- function(v) c(quantile(v, .025), quantile(v, .975))
res <- bind_rows(lapply(names(FEATURE_SETS), function(s) {
  m <- metrics(wide$actual, wide[[s]])
  data.frame(set = s,
             RMSE = m["RMSE"], RMSE_lo = ci(boot[paste0(s, ".RMSE"), ])[1], RMSE_hi = ci(boot[paste0(s, ".RMSE"), ])[2],
             MAE = m["MAE"], MAE_lo = ci(boot[paste0(s, ".MAE"), ])[1], MAE_hi = ci(boot[paste0(s, ".MAE"), ])[2],
             r = m["r"], r_lo = ci(boot[paste0(s, ".r"), ])[1], r_hi = ci(boot[paste0(s, ".r"), ])[2],
             R2 = m["R2"], R2_lo = ci(boot[paste0(s, ".R2"), ])[1], R2_hi = ci(boot[paste0(s, ".R2"), ])[2])
}))
difference <- data.frame(d_r = cor(wide$actual, wide$combined) - cor(wide$actual, wide$language),
                         d_r_lo = ci(boot["d_r", ])[1], d_r_hi = ci(boot["d_r", ])[2],
                         d_RMSE = sqrt(mean((wide$actual - wide$combined)^2)) - sqrt(mean((wide$actual - wide$language)^2)),
                         d_RMSE_lo = ci(boot["d_RMSE", ])[1], d_RMSE_hi = ci(boot["d_RMSE", ])[2],
                         p_comb_gt_lang = mean(boot["d_r", ] > 0))
shap_sum <- shap %>% group_by(feature) %>% summarise(mean_abs_shap = mean(abs(shap)), .groups = "drop") %>% arrange(desc(mean_abs_shap))

write_csv(preds, output_path("svr_loo_predictions.csv"))
write_csv(shap, output_path("svr_shap_combined.csv"))
write_csv(hp, output_path("svr_hyperparameters.csv"))
write_csv(res, output_path("svr_metrics_by_featureset.csv"))
write_csv(difference, output_path("svr_featureset_difference.csv"))
write_csv(shap_sum, output_path("svr_shap_summary.csv"))

cat(sprintf("\nSVR, nested leave-one-out, n = %d PWA, %.1f min (95%% bootstrap CI)\n", nrow(wide),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
for (k in seq_len(nrow(res))) with(res[k, ], cat(sprintf(
  "  %-11s r %.3f [%.3f, %.3f]  RMSE %.2f [%.2f, %.2f]  MAE %.2f [%.2f, %.2f]  R2 %.3f\n",
  set, r, r_lo, r_hi, RMSE, RMSE_lo, RMSE_hi, MAE, MAE_lo, MAE_hi, R2)))
with(difference, cat(sprintf("  combined minus language: r %.3f [%.3f, %.3f]  RMSE %.2f [%.2f, %.2f]\n",
                             d_r, d_r_lo, d_r_hi, d_RMSE, d_RMSE_lo, d_RMSE_hi)))
cat("  mean |SHAP|, combined model:", paste(sprintf("%s %.2f", shap_sum$feature, shap_sum$mean_abs_shap), collapse = ", "), "\n")
sel <- hp %>% filter(set == "combined") %>% count(cost, epsilon)
cat("  combined model, selected (cost, epsilon):", paste(sprintf("(%g, %g) %d", sel$cost, sel$epsilon, sel$n), collapse = "; "), "\n")
