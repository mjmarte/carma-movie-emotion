# Classification of PWA versus HC by radial-kernel support vector machines (Results 2.6, Table 3, Figure 6; Methods 4.7):
# NREP repeats (default 30) of stratified 5-fold cross-validation with 3 inner tuning folds, reference-dependent features recomputed
# from each outer training set and, with STRICT=1 (default), each inner one; DeLong tests per repeat; SHAP values in repeat 1.
# Writes .csv files to svm_strict/ (svm_outer/ if STRICT=0): svm_summary, svm_per_repeat, svm_fold_metrics, svm_hyperparameters,
# svm_delong, svm_delong_summary, svm_featureset_differences, svm_oof_predictions_rep1, svm_shap_combined_rep1, svm_shap_summary_rep1.
source("R/config.R")
suppressPackageStartupMessages({
  library(e1071)
  library(pROC)
  library(caret)
  library(kernelshap)
})
source("R/foldwise_features.R")

STRICT <- Sys.getenv("STRICT", "1") == "1"
NREP <- as.integer(Sys.getenv("NREP", "30"))
OUTD <- if (STRICT) "svm_strict" else "svm_outer"
validate_features()
GRID <- expand.grid(cost = c(0.1, 1, 10), gamma = c(0.01, 0.1, 1))
addo <- function(f) f %>% mutate(outcome = factor(Group, levels = c("HC", "PWA")))
# Tunes cost and gamma by mean inner-fold AUC, refits on the outer training set and predicts the held-out fold
fit_fold <- function(train, test, feats, inner_sets) {
  res <- matrix(NA, nrow(GRID), length(inner_sets))
  for (j in seq_along(inner_sets)) {
    tr <- inner_sets[[j]]$tr; va <- inner_sets[[j]]$va
    for (p in 1:nrow(GRID)) {
      m <- svm(x = tr[, feats, drop = FALSE], y = tr$outcome, kernel = "radial", cost = GRID$cost[p], gamma = GRID$gamma[p],
               probability = TRUE)
      pr <- attr(predict(m, va[, feats, drop = FALSE], probability = TRUE), "probabilities")[, "PWA"]
      res[p, j] <- if (length(unique(va$outcome)) < 2) 0.5 else as.numeric(auc(roc(va$outcome == "PWA", pr, quiet = TRUE)))
    }
  }
  b <- which.max(rowMeans(res, na.rm = TRUE))
  m <- svm(x = train[, feats, drop = FALSE], y = train$outcome, kernel = "radial", cost = GRID$cost[b], gamma = GRID$gamma[b],
           probability = TRUE)
  list(model = m, prob = attr(predict(m, test[, feats, drop = FALSE], probability = TRUE), "probabilities")[, "PWA"],
       pred = as.character(predict(m, test[, feats, drop = FALSE])), cost = GRID$cost[b], gamma = GRID$gamma[b])
}
cm_metrics <- function(pred, actual) {
  tp <- sum(pred == "PWA" & actual == "PWA"); tn <- sum(pred == "HC" & actual == "HC")
  fp <- sum(pred == "PWA" & actual == "HC"); fn <- sum(pred == "HC" & actual == "PWA")
  sens <- tp / (tp + fn); spec <- tn / (tn + fp); prec <- if (tp + fp > 0) tp / (tp + fp) else NA
  f1 <- if (!is.na(prec) && prec + sens > 0) 2 * prec * sens / (prec + sens) else NA
  c(accuracy = (tp + tn) / length(pred), sensitivity = sens, specificity = spec, precision = prec, f1 = f1)
}
prob_pwa <- function(mod, nd) attr(predict(mod, as.data.frame(nd), probability = TRUE), "probabilities")[, "PWA"]

rows <- list(); delong <- list(); shap_rows <- list(); hp <- list(); oof1 <- NULL; t0 <- Sys.time()
for (r in 1:NREP) {
  set.seed(2026 + r); base <- addo(IDS); folds <- createFolds(base$outcome, k = 5); probs <- list(); fm <- list()
  for (i in 1:5) {
    rm(list = ls(.RT_CACHE), envir = .RT_CACHE)
    test_ids <- base$Participant[folds[[i]]]; train_ids <- setdiff(ALL_IDS, test_ids)
    f_outer <- addo(compute_features(train_ids)); stopifnot(!anyNA(f_outer[, FEATS_ALL]))
    train <- f_outer[f_outer$Participant %in% train_ids, ]; test <- f_outer[f_outer$Participant %in% test_ids, ]
    inner <- createFolds(train$outcome, k = 3)
    inner_sets <- lapply(1:3, function(j) {
      va_ids <- train$Participant[inner[[j]]]; tr_ids <- setdiff(train$Participant, va_ids)
      fj <- if (STRICT) addo(compute_features(tr_ids)) else f_outer
      list(tr = fj[fj$Participant %in% tr_ids, ], va = fj[fj$Participant %in% va_ids, ])
    })
    for (s in names(FEATURE_SETS)) {
      ff <- fit_fold(train, test, FEATURE_SETS[[s]], inner_sets)
      probs[[s]] <- rbind(probs[[s]], data.frame(Participant = test$Participant, actual = as.character(test$outcome),
                                                 prob = ff$prob, pred = ff$pred))
      fm[[length(fm) + 1]] <- data.frame(rep = r, fold = i, set = s,
                                         auc = as.numeric(auc(roc(test$outcome == "PWA", ff$prob, quiet = TRUE))),
                                         t(cm_metrics(ff$pred, as.character(test$outcome))))
      hp[[length(hp) + 1]] <- data.frame(rep = r, fold = i, set = s, cost = ff$cost, gamma = ff$gamma)
      if (s == "combined" && r == 1) {
        set.seed(100 + i); bg <- train[sample(nrow(train), 20), FEATS_ALL]
        ks <- kernelshap(ff$model, X = test[, FEATS_ALL], bg_X = bg, pred_fun = prob_pwa, verbose = FALSE)
        shap_rows[[length(shap_rows) + 1]] <- data.frame(fold = i, Participant = rep(test$Participant, each = length(FEATS_ALL)),
                                                         feature = rep(FEATS_ALL, nrow(test)), shap = as.numeric(t(ks$S)),
                                                         value = as.numeric(t(as.matrix(test[, FEATS_ALL]))))
      }
    }
  }
  fm <- bind_rows(fm); rows[[r]] <- fm
  rc <- roc(probs$combined$actual == "PWA", probs$combined$prob, quiet = TRUE)
  rl <- roc(probs$language$actual == "PWA", probs$language$prob, quiet = TRUE)
  re <- roc(probs$emotion$actual == "PWA", probs$emotion$prob, quiet = TRUE)
  delong[[r]] <- data.frame(rep = r, p_comb_vs_lang = roc.test(rc, rl, method = "delong", paired = TRUE)$p.value,
                            p_comb_vs_emo = roc.test(rc, re, method = "delong", paired = TRUE)$p.value,
                            pooled_auc_comb = as.numeric(auc(rc)), pooled_auc_lang = as.numeric(auc(rl)),
                            pooled_auc_emo = as.numeric(auc(re)))
  if (r == 1) oof1 <- bind_rows(lapply(names(probs), function(s) cbind(set = s, probs[[s]])))
  cat(sprintf("repeat %d/%d done (%.1f min elapsed)\n", r, NREP, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

fm <- bind_rows(rows); dl <- bind_rows(delong); hp <- bind_rows(hp); shap <- bind_rows(shap_rows)
metric_cols <- c("auc", "accuracy", "sensitivity", "specificity", "precision", "f1")
per_rep <- fm %>% group_by(rep, set) %>% summarise(across(all_of(metric_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop")
qs <- function(x) c(mean = mean(x), lo = unname(quantile(x, .025)), hi = unname(quantile(x, .975)))
summ <- per_rep %>% group_by(set) %>%
  summarise(across(all_of(metric_cols), list(mean = ~mean(.x), lo = ~unname(quantile(.x, .025)), hi = ~unname(quantile(.x, .975)))),
            .groups = "drop")
w <- per_rep %>% select(rep, set, auc) %>% pivot_wider(names_from = set, values_from = auc)
wa <- per_rep %>% select(rep, set, accuracy) %>% pivot_wider(names_from = set, values_from = accuracy)
diffs <- data.frame(contrast = c("combined_minus_language_auc", "combined_minus_emotion_auc", "language_minus_antgen_auc",
                                 "combined_minus_language_accuracy"),
                    rbind(qs(w$combined - w$language), qs(w$combined - w$emotion), qs(w$language - w$antgen_only),
                          qs(wa$combined - wa$language)),
                    p_gt0 = c(mean(w$combined - w$language > 0), mean(w$combined - w$emotion > 0),
                              mean(w$language - w$antgen_only > 0), mean(wa$combined - wa$language > 0)))
dls <- data.frame(median_p_comb_vs_lang = median(dl$p_comb_vs_lang), share_p05_comb_vs_lang = mean(dl$p_comb_vs_lang < .05),
                  median_p_comb_vs_emo = median(dl$p_comb_vs_emo),
                  t(c(pooled_comb = qs(dl$pooled_auc_comb), pooled_lang = qs(dl$pooled_auc_lang), pooled_emo = qs(dl$pooled_auc_emo))))
shap_sum <- shap %>% group_by(feature) %>% summarise(mean_abs_shap = mean(abs(shap)), .groups = "drop") %>% arrange(desc(mean_abs_shap))

write_csv(fm, output_path(OUTD, "svm_fold_metrics.csv"))
write_csv(per_rep, output_path(OUTD, "svm_per_repeat.csv"))
write_csv(summ, output_path(OUTD, "svm_summary.csv"))
write_csv(diffs, output_path(OUTD, "svm_featureset_differences.csv"))
write_csv(dl, output_path(OUTD, "svm_delong.csv"))
write_csv(dls, output_path(OUTD, "svm_delong_summary.csv"))
write_csv(hp, output_path(OUTD, "svm_hyperparameters.csv"))
write_csv(shap, output_path(OUTD, "svm_shap_combined_rep1.csv"))
write_csv(shap_sum, output_path(OUTD, "svm_shap_summary_rep1.csv"))
write_csv(oof1, output_path(OUTD, "svm_oof_predictions_rep1.csv"))

cat(sprintf("\nSVM, %d repeats, STRICT = %s, %.1f min (mean over folds and repeats [2.5th, 97.5th percentile across repeats])\n",
            NREP, STRICT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
for (k in seq_len(nrow(summ))) with(summ[k, ], cat(sprintf(
  "  %-11s AUC %.3f [%.3f, %.3f]  accuracy %.3f [%.3f, %.3f]  sensitivity %.3f  specificity %.3f  precision %.3f  F1 %.3f\n",
  set, auc_mean, auc_lo, auc_hi, accuracy_mean, accuracy_lo, accuracy_hi, sensitivity_mean, specificity_mean, precision_mean,
  f1_mean)))
for (k in seq_len(nrow(diffs))) cat(sprintf("  %-33s %.3f [%.3f, %.3f]  share of repeats > 0: %.3f\n", diffs$contrast[k],
                                            diffs$mean[k], diffs$lo[k], diffs$hi[k], diffs$p_gt0[k]))
cat(sprintf("  DeLong, combined vs language: median p %.2f (p < .05 in %.1f%% of repeats); combined vs emotion: median p %.2g\n",
            dls$median_p_comb_vs_lang, 100 * dls$share_p05_comb_vs_lang, dls$median_p_comb_vs_emo))
cat("  mean |SHAP|, combined model, repeat 1:",
    paste(sprintf("%s %.3f", shap_sum$feature, shap_sum$mean_abs_shap), collapse = ", "), "\n")
sel <- hp %>% filter(set == "combined") %>% count(cost, gamma) %>% arrange(desc(n))
cat("  combined model, selected (cost, gamma):", paste(sprintf("(%g, %g) %d", sel$cost, sel$gamma, sel$n), collapse = "; "), "\n")
