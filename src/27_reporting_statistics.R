# Reporting statistics: the degrees of freedom, 95% confidence intervals, effect sizes and exact p values reported in the
# manuscript and the Supplementary Tables, computed from the saved outputs of the earlier scripts (no value is typed by hand).
# Confidence intervals are t-based (estimate +/- qt(.975, df) * SE) for the mixed and linear models, Wald (estimate +/- 1.96 * SE)
# for the two binomial GLMM group terms, and Fisher-z for correlations (using n). Welch t tests and Cohen's d (pooled SD) come
# from the group means, SDs and ns; phi and Cramer's V come from the uncorrected chi-square. Writes reporting_statistics.csv.
source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

rows <- list()
add <- function(item, statistic, df = NA, estimate = NA, ci_lo = NA, ci_hi = NA, effect_size = NA, effect_size_type = NA, p = NA) {
  rows[[length(rows) + 1]] <<- data.frame(item = item, statistic = statistic, df = df, estimate = estimate, ci_lo = ci_lo,
                                          ci_hi = ci_hi, effect_size = effect_size, effect_size_type = effect_size_type, p = p)
}
t_ci <- function(est, se, df) est + c(-1, 1) * qt(.975, df) * se
r_ci <- function(r, n) tanh(atanh(r) + c(-1, 1) * 1.96 / sqrt(n - 3))
rd <- function(f) read.csv(output_path(f), check.names = FALSE)

# Group models and the entropy window terms
gm <- rd("group_models.csv")
lmm <- function(item, model, term) {
  r <- gm[gm$model == model & gm$term == term, ][1, ]
  ci <- t_ci(r$estimate, r$SE, r$df)
  add(item, "t", r$df, r$estimate, ci[1], ci[2], p = r$p)
}
lmm("holistic", "Holistic MSZ", "GroupPWA")
lmm("realtime", "Real-time MSZ, all-participant reference (Results 2.3)", "GroupPWA")
lmm("isc", "ISC", "GroupPWA")
lmm("ent_med", "Sample entropy", "WindowBinMedium:GroupPWA")
lmm("ent_long", "Sample entropy", "WindowBinLong:GroupPWA")
lmm("ent_age", "Sample entropy", "Age_centered")

# Comprehension and antonym accuracy: binomial GLMM group terms (log odds, Wald CI) and the odds ratio
ta <- rd("task_accuracy_results.csv")
for (k in list(c("CompQ", "compq"), c("AntGen", "antgen"))) {
  r <- ta[ta$result == "glmm" & ta$Task == k[1] & ta$term == "GroupPWA", ][1, ]
  lo <- r$estimate - 1.96 * r$se; hi <- r$estimate + 1.96 * r$se
  add(k[2], "z", NA, r$estimate, lo, hi, p = r$p)
  add(paste0(k[2], "_OR"), "z", NA, exp(r$estimate), exp(lo), exp(hi), exp(r$estimate), "odds ratio", r$p)
}

# Principal component group differences
pc <- rd("pc_group_differences_pooled_msz.csv")
for (i in which(pc$PC %in% c("PC1", "PC2", "PC3"))) {
  r <- pc[i, ]; ci <- t_ci(r$Difference, r$SE, r$df)
  add(r$PC, "t", r$df, r$Difference, ci[1], ci[2], r$Cohen_d, "Cohen's d", r$p_value)
}

# PC1 and AQ correlations (Figure 5)
f5 <- rd("figure5_correlations.csv")
for (i in seq_len(nrow(f5))) {
  r <- f5[i, ]; ci <- r_ci(r$r, r$n)
  add(paste0("cor_", r$label), "r", r$n - 2, r$r, ci[1], ci[2], r$r, "Pearson r", r$p)
}

# Table 1: Welch tests for age and, from the mood subsample, PHQ-8 and GAD-7
welch <- function(item, m1, s1, n1, m2, s2, n2) {
  se <- sqrt(s1^2 / n1 + s2^2 / n2); tt <- (m1 - m2) / se
  df <- (s1^2 / n1 + s2^2 / n2)^2 / ((s1^2 / n1)^2 / (n1 - 1) + (s2^2 / n2)^2 / (n2 - 1))
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)); ci <- t_ci(m1 - m2, se, df)
  add(item, "t", df, m1 - m2, ci[1], ci[2], (m1 - m2) / sp, "Cohen's d", 2 * pt(abs(tt), df, lower.tail = FALSE))
}
t1 <- rd("table1.csv")
a <- t1[t1$variable == "age" & t1$test == "Welch t-test", ][1, ]
welch("age", a$PWA_mean, a$PWA_sd, a$PWA_n, a$HC_mean, a$HC_sd, a$HC_n)
mo <- rd("mood_subsample.csv"); rownames(mo) <- mo$Group
welch("phq", mo["PWA", "PHQ_mean"], mo["PWA", "PHQ_sd"], mo["PWA", "n_PHQ"], mo["HC", "PHQ_mean"], mo["HC", "PHQ_sd"], mo["HC", "n_PHQ"])
welch("gad", mo["PWA", "GAD_mean"], mo["PWA", "GAD_sd"], mo["PWA", "n_GAD"], mo["HC", "GAD_mean"], mo["HC", "GAD_sd"], mo["HC", "n_GAD"])

# Sex and race: Yates chi-square, with phi from the uncorrected chi-square (N = 100)
sr <- rd("sex_race_tests.csv")
for (k in list(c("sex (self-reported)", "sex"), c("race white vs non-white (Table 1 check)", "race"))) {
  none <- sr[sr$variable == k[1] & sr$correction == "none", ][1, ]; yat <- sr[sr$variable == k[1] & sr$correction == "Yates", ][1, ]
  add(k[2], "chisq", 1, effect_size = sqrt(none$chisq / 100), effect_size_type = "phi", p = yat$p, estimate = yat$chisq)
}

# Ethnicity: Fisher exact odds ratio with Wald CI on the log odds ratio (unknown excluded), counts from Table 1
eth <- t1[t1$variable == "ethnicity" & t1$characteristic %in% c("Hispanic/Latino", "Not Hispanic/Latino"), ]
m <- rbind(c(eth$PWA_count[1], eth$PWA_count[2]), c(eth$HC_count[1], eth$HC_count[2]))
ft <- fisher.test(m)
or <- (m[1, 1] * m[2, 2]) / (m[1, 2] * m[2, 1]); w <- 1.96 * sqrt(sum(1 / m))
add("ethnicity", "Fisher exact", NA, or, exp(log(or) - w), exp(log(or) + w), or, "odds ratio", ft$p.value)

# Education: Cramer's V from the uncorrected chi-square on the four levels (not reported excluded)
ed <- t1[t1$variable == "education" & !t1$characteristic %in% c("Education, n (%)", "Not reported"), ]
cs <- suppressWarnings(chisq.test(cbind(PWA = ed$PWA_count, HC = ed$HC_count), correct = FALSE))
add("education", "chisq", unname(cs$parameter), unname(cs$statistic), effect_size = sqrt(unname(cs$statistic) / (sum(cs$observed) * 1)),
    effect_size_type = "Cramer's V", p = cs$p.value)

# Supplementary Table S1 (GroupPWA rows), S2 (correlations, BH-adjusted p) and S3 (Pearson chi-square, phi)
s1 <- rd("table_s1_sensitivity.csv"); s1 <- s1[s1$term == "GroupPWA", ]
for (i in seq_len(nrow(s1))) {
  r <- s1[i, ]; ci <- t_ci(r$estimate, r$SE, r$df)
  add(paste("S1", r$panel, r$feature, r$specification, sep = "|"), "t", r$df, r$estimate, ci[1], ci[2], p = r$p)
}
s2 <- rd("table_s2.csv")
for (i in seq_len(nrow(s2))) {
  r <- s2[i, ]; ci <- r_ci(r$r, r$n)
  add(paste("S2", r$feature, r$measure, sep = "|"), "r", r$n - 2, r$r, ci[1], ci[2], r$r, "Pearson r", r$p_bh)
}
s3 <- rd("table_s3.csv")
for (i in seq_len(nrow(s3))) {
  r <- s3[i, ]; N <- r$n_PWA + r$n_HC
  add(paste("S3", r$clip, sep = "|"), "chisq", 1, r$chisq_pearson, effect_size = sqrt(r$chisq_pearson / N), effect_size_type = "phi", p = r$p_pearson)
}
# Supplementary Figure S1 correlations (saved confidence limits)
fs <- read.csv(output_path("figure_s1_correlations.csv"), row.names = 1)
for (g in c("PWA", "HC")) add(paste0("S1fig_", g), "r", fs[g, "n"] - 2, fs[g, "r"], fs[g, "lo"], fs[g, "hi"], fs[g, "r"], "Pearson r", fs[g, "p"])

out <- do.call(rbind, rows)
write.csv(out, output_path("reporting_statistics.csv"), row.names = FALSE)
