# Age-adjusted group models (Results 2.3, Methods 4.5): linear mixed models DV ~ Group + Age + (1|Participant) + (1|Clip),
# Satterthwaite df, for holistic MSZ, real-time MSZ, window-level ISC, and sample entropy with Window Size and its interaction
# with Group; age centred on the sample mean. Results 2.3 reports real-time MSZ referenced to all participants; the HC-referenced
# version is the Section 2.7 baseline. PC1-PC3 group contrasts (Methods 4.6): Welch t-tests and PC ~ Group + Age, on the main PCA
# and on the PCA of Results 2.4. Writes group_models.csv and pc_contrasts.csv.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

ids <- read.csv(output_path("analytic_ids.csv"))
age <- ids %>%
  left_join(read.csv(input_file("demographics")) %>% select(Participant, Age), by = "Participant") %>%
  mutate(Age_centered = Age - mean(Age, na.rm = TRUE)) %>%
  select(Participant, Age, Age_centered)
stopifnot(nrow(age) == nrow(ids), !anyNA(age$Age))
add_age <- function(d) {
  out <- d %>% left_join(age, by = "Participant")
  stopifnot(nrow(out) == nrow(d), all(out$Participant %in% ids$Participant), !anyNA(out$Age))
  out
}

fit_lmm <- function(label, formula, d) {
  model <- lmer(formula, data = add_age(d))
  co <- summary(model)$coefficients
  data.frame(model = label, term = rownames(co), estimate = co[, "Estimate"], SE = co[, "Std. Error"], df = co[, "df"],
             t = co[, "t value"], p = co[, "Pr(>|t|)"], n_participants = unname(ngrps(model)["Participant"]),
             n_obs = nobs(model), row.names = NULL)
}
# lme4 flags the holistic MSZ fit (max|grad| 0.004) because its clip variance is near zero; other optimizers give the same
# estimates.
group_models <- bind_rows(
  fit_lmm("Holistic MSZ", msz ~ Group + Age_centered + (1 | Participant) + (1 | Clip), readRDS(output_path("holistic_msz.rds"))),
  fit_lmm("Real-time MSZ, all-participant reference (Results 2.3)", msz ~ Group + Age_centered + (1 | Participant) + (1 | Clip),
          readRDS(output_path("realtime_msz_pooled.rds"))),
  fit_lmm("Real-time MSZ, HC reference (Results 2.7)", msz ~ Group + Age_centered + (1 | Participant) + (1 | Clip),
          readRDS(output_path("realtime_msz.rds"))),
  fit_lmm("ISC", isc_z ~ Group + Age_centered + (1 | Participant) + (1 | Clip), readRDS(output_path("isc.rds"))),
  fit_lmm("Sample entropy", SampEn ~ WindowBin * Group + Age_centered + (1 | Participant) + (1 | Clip),
          readRDS(output_path("sample_entropy.rds"))))
write.csv(group_models, output_path("group_models.csv"), row.names = FALSE)

pc_group_contrasts <- function(file, pca) {
  scores <- add_age(read.csv(output_path(file)))
  bind_rows(lapply(paste0("PC", 1:3), function(pc) {
    hc <- scores[[pc]][scores$Group == "HC"]
    pwa <- scores[[pc]][scores$Group == "PWA"]
    welch <- t.test(hc, pwa)
    pooled_sd <- sqrt(((length(hc) - 1) * var(hc) + (length(pwa) - 1) * var(pwa)) / (length(hc) + length(pwa) - 2))
    model <- lm(as.formula(paste(pc, "~ Group + Age_centered")), data = scores)
    co <- summary(model)$coefficients
    bind_rows(
      data.frame(PC = pc, test = "Welch t-test", term = "HC - PWA", estimate = mean(hc) - mean(pwa),
                 SE = sqrt(var(hc) / length(hc) + var(pwa) / length(pwa)), df = unname(welch$parameter),
                 t = unname(welch$statistic), p = welch$p.value, cohens_d = (mean(hc) - mean(pwa)) / pooled_sd),
      data.frame(PC = pc, test = "PC ~ Group + Age", term = rownames(co), estimate = co[, "Estimate"], SE = co[, "Std. Error"],
                 df = df.residual(model), t = co[, "t value"], p = co[, "Pr(>|t|)"], cohens_d = NA_real_, row.names = NULL)) %>%
      mutate(n_HC = length(hc), n_PWA = length(pwa), pca = pca)
  }))
}
pc_contrasts <- bind_rows(pc_group_contrasts("pc_scores.csv", "main PCA (Figures 4 and 5)"),
                          pc_group_contrasts("pc_scores_pooled_msz.csv", "PCA of Results 2.4"))
write.csv(pc_contrasts, output_path("pc_contrasts.csv"), row.names = FALSE)

for (k in grep("Group|Age", group_models$term)) {
  g <- group_models[k, ]
  cat(sprintf("%s, %s: b = %.4f, SE = %.4f, t = %.3f, p = %.3g (%d participants, %d observations)\n", g$model, g$term,
              g$estimate, g$SE, g$t, g$p, g$n_participants, g$n_obs))
}
for (k in which(pc_contrasts$term %in% c("HC - PWA", "GroupPWA"))) {
  g <- pc_contrasts[k, ]
  cat(sprintf("%s, %s, %s, %s: b = %.3f, SE = %.3f, t = %.2f, p = %.3g\n", g$pca, g$PC, g$test, g$term, g$estimate, g$SE, g$t, g$p))
}
