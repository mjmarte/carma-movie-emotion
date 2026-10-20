# Supplementary and descriptive analyses: holistic MSZ with an HC-only reference (Supplementary text), the correlations of the
# movie features with the executive and attention measures within PWA (input to Supplementary Table S2), education, sex, race
# and mood by group (Results 2.1, Table 1), the group effects refit with education level as a covariate, and the joint models
# of PC1 on aphasia severity and a mood score within PWA (reported in the response to the reviewers). Writes
# holistic_msz_reference*.csv, cognitive_*.csv, education_*.csv, sex_*.csv, mood_*.csv and pc1_*.csv.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

ids <- read.csv(output_path("analytic_ids.csv"))
holistic <- readRDS(output_path("holistic_msz.rds"))
realtime <- readRDS(output_path("realtime_msz.rds"))
isc <- readRDS(output_path("isc.rds"))
entropy <- readRDS(output_path("sample_entropy.rds"))
tasks <- readRDS(output_path("task_responses.rds"))
stopifnot(all(unlist(lapply(list(holistic, realtime, isc, entropy, tasks), `[[`, "Group")) %in% c("HC", "PWA")))
panas_items <- read.csv(input_file("panas_items"))
dem <- read_csv(input_file("demographics"), show_col_types = FALSE)
aq <- read_csv(input_file("aphasia_quotient"), show_col_types = FALSE) %>% select(Participant, AQ)
pc <- read_csv(output_path("pc_scores.csv"), show_col_types = FALSE) %>% select(Participant, PC1)
questionnaire <- read_csv(output_path("demographics_questionnaire.csv"), show_col_types = FALSE)
clin <- dem %>% select(Participant, Age, PHQ_Total, GAD_Total, elevator_strings, tmtpartbtime, tmtratio) %>%
  mutate(tmtpartbtime = suppressWarnings(as.numeric(tmtpartbtime)))
d <- as.data.frame(readRDS(output_path("pca_data.rds"))) %>% select(-outcome) %>% left_join(pc, by = "Participant") %>%
  left_join(clin, by = "Participant") %>% left_join(aq, by = "Participant") %>%
  left_join(questionnaire %>% select(-Group), by = "Participant") %>% mutate(Group = factor(Group, levels = c("HC", "PWA")))
stopifnot(nrow(d) == 100, sum(d$Group == "PWA") == 57, !anyNA(d$PC1), !anyNA(d$AQ[d$Group == "PWA"]))
pwa <- d %>% filter(Group == "PWA")
ct <- function(x, y, m = "pearson") {
  ok <- complete.cases(x, y)
  t <- suppressWarnings(cor.test(x[ok], y[ok], method = m, exact = FALSE))
  data.frame(r = unname(t$estimate), p = t$p.value, n = sum(ok))
}
tidy_lm <- function(m, label) {
  cf <- summary(m)$coefficients
  data.frame(model = label, term = rownames(cf), estimate = cf[, 1], SE = cf[, 2], t = cf[, 3], p = cf[, 4], n = nobs(m),
             R2 = summary(m)$r.squared, row.names = NULL)
}

# Joint models of PC1 on AQ and a mood score within PWA; group-by-score interaction and simple slopes across groups
z <- function(x) as.numeric(scale(x))
pp <- pwa %>% filter(!is.na(PHQ_Total))
joint <- bind_rows(tidy_lm(lm(PC1 ~ AQ + PHQ_Total, pwa), "PC1 ~ AQ + PHQ-8"),
                   tidy_lm(lm(PC1 ~ AQ + PHQ_Total + Age, pwa), "PC1 ~ AQ + PHQ-8 + Age"),
                   tidy_lm(lm(PC1 ~ AQ + GAD_Total, pwa), "PC1 ~ AQ + GAD-7"),
                   tidy_lm(lm(PC1 ~ AQ + GAD_Total + Age, pwa), "PC1 ~ AQ + GAD-7 + Age"),
                   tidy_lm(lm(z(PC1) ~ z(AQ) + z(PHQ_Total) + z(Age), pp), "standardized: PC1 ~ AQ + PHQ-8 + Age"))
phq_data <- d %>% filter(!is.na(PHQ_Total))
gad_data <- d %>% filter(!is.na(GAD_Total))
m_phq <- lm(PC1 ~ PHQ_Total * Group + Age, phq_data)
m_gad <- lm(PC1 ~ GAD_Total * Group + Age, gad_data)
interaction <- bind_rows(tidy_lm(m_phq, "PC1 ~ PHQ-8 x Group + Age"),
                         tidy_lm(lm(PC1 ~ PHQ_Total * Group, phq_data), "PC1 ~ PHQ-8 x Group"),
                         tidy_lm(m_gad, "PC1 ~ GAD-7 x Group + Age"))
slopes <- bind_rows(cbind(score = "PHQ-8", as.data.frame(summary(emtrends(m_phq, ~Group, var = "PHQ_Total")))),
                    cbind(score = "GAD-7", as.data.frame(summary(emtrends(m_gad, ~Group, var = "GAD_Total")))))
write_csv(joint, output_path("pc1_joint_models.csv"))
write_csv(interaction, output_path("pc1_group_interaction.csv"))
write_csv(slopes, output_path("pc1_simple_slopes.csv"))

# Supplement: holistic MSZ with an HC-only reference (leave-one-out for HC), against the pooled reference of the main text
age <- dem %>% select(Participant, Age)
il <- panas_items %>% rename(Participant = participant_id, Clip = clip, Rating = rating, PANAS_Item = panas_item) %>%
  select(Participant, Clip, PANAS_Item, Rating) %>% inner_join(ids, by = "Participant")
pooled <- il %>% mutate(.by = c(Clip, PANAS_Item), z = (Rating - mean(Rating)) / sd(Rating)) %>%
  group_by(Participant, Group, Clip) %>% summarise(msz = mean(z^2), .groups = "drop")
hc_stats <- il %>% filter(Group == "HC") %>% group_by(Clip, PANAS_Item) %>%
  summarise(hs = sum(Rating), hss = sum(Rating^2), hn = n(), .groups = "drop")
hcref <- il %>% left_join(hc_stats, by = c("Clip", "PANAS_Item")) %>%
  mutate(self = ifelse(Group == "HC", Rating, 0), selfn = ifelse(Group == "HC", 1, 0), n_ref = hn - selfn,
         mean_ref = (hs - self) / n_ref, var_ref = ((hss - self^2) - n_ref * mean_ref^2) / (n_ref - 1),
         sd_ref = sqrt(pmax(var_ref, 0)))
hc_reference <- hcref %>% mutate(z = (Rating - mean_ref) / sd_ref) %>% filter(is.finite(z)) %>%
  group_by(Participant, Group, Clip) %>% summarise(msz = mean(z^2), .groups = "drop")
fit_grp <- function(df, label) {
  x <- df %>% left_join(age, by = "Participant") %>% mutate(Group = factor(Group, levels = c("HC", "PWA")))
  cf <- summary(lmer(msz ~ Group + Age + (1 | Participant) + (1 | Clip), x))$coefficients["GroupPWA", ]
  data.frame(reference = label, estimate = cf[1], SE = cf[2], df = cf[3], t = cf[4], p = cf[5], n_participant_clips = nrow(x),
             row.names = NULL)
}
reference <- bind_rows(fit_grp(pooled, "Pooled reference (all participants), as reported"),
                       fit_grp(hc_reference, "HC reference, leave-one-out for HC, undefined z-scores excluded"))
pl <- pooled %>% group_by(Participant) %>% summarise(pooled = mean(msz)) %>%
  inner_join(hc_reference %>% group_by(Participant) %>% summarise(hcref = mean(msz)), by = "Participant")
reference_info <- data.frame(item_observations = nrow(hcref), zero_sd_cells = sum(hcref$sd_ref == 0),
                             participant_level_r_pooled_vs_hcref = cor(pl$pooled, pl$hcref), n_participants = nrow(pl))
write_csv(reference, output_path("holistic_msz_reference.csv"))
write_csv(reference_info, output_path("holistic_msz_reference_info.csv"))

feats <- c(PANAS_MSZ_mean = "Holistic MSZ", CARMA_MSZ_mean = "Real-time MSZ", ISC_mean = "Real-time ISC (mean)",
           ISC_sd = "Real-time ISC (SD)", SampEn_Medium = "Sample entropy (medium)", CompQ_mean = "Movie comprehension accuracy",
           AntGen_mean = "Antonym generation accuracy")
# Supplementary Table S2 (built by 23_table_s2.R): the five emotion features, the two language accuracies and PC1 against the
# executive and attention measures within PWA (Pearson, Spearman, and partial correlations controlling for AQ)
cogs <- c(tmtratio = "TMT B/A ratio", tmtpartbtime = "TMT Part B time (s)", elevator_strings = "TEA Elevator Strings")
partial_aq <- function(f, cg) {
  sub <- pwa %>% select(all_of(c(f, cg, "AQ"))) %>% na.omit()
  names(sub) <- c("x", "y", "AQ")
  cor.test(resid(lm(x ~ AQ, sub)), resid(lm(y ~ AQ, sub)))
}
cognitive_load <- bind_rows(lapply(names(feats)[1:5], function(f) bind_rows(lapply(names(cogs), function(cg) {
  pe <- ct(pwa[[f]], pwa[[cg]])
  sp <- ct(pwa[[f]], pwa[[cg]], "spearman")
  pt <- partial_aq(f, cg)
  data.frame(feature = feats[[f]], measure = cogs[[cg]], n = pe$n, r = pe$r, p = pe$p, rho = sp$r, p_rho = sp$p,
             r_partial_AQ = unname(pt$estimate), p_partial = pt$p.value)
})))) %>% mutate(p_fdr = p.adjust(p, "BH"), p_partial_fdr = p.adjust(p_partial, "BH"))
language <- c(CompQ_mean = "Movie comprehension accuracy", AntGen_mean = "Antonym generation accuracy", PC1 = "PC1")
cognitive_language <- bind_rows(lapply(names(language), function(f) bind_rows(lapply(names(cogs), function(cg) {
  pt <- partial_aq(f, cg)
  cbind(feature = language[[f]], measure = cogs[[cg]], ct(pwa[[f]], pwa[[cg]]), r_partial_AQ = unname(pt$estimate),
        p_partial = pt$p.value)
}))))
write_csv(cognitive_load, output_path("cognitive_load_correlations.csv"))
write_csv(cognitive_language, output_path("cognitive_language_correlations.csv"))

# Table 1, Results 2.1: education by group, and the group effects refit with education level as a covariate
lev <- c("High School/GED", "Part College or 2-year degree", "4-year College degree", "Part or completed Graduate degree")
d <- d %>% mutate(edu_level = match(education, lev), edu_label = ifelse(is.na(edu_level), "Not reported", lev[edu_level]))
education_counts <- d %>% count(Group, edu_label) %>% pivot_wider(names_from = Group, values_from = n, values_fill = 0) %>%
  mutate(order = match(edu_label, c(lev, "Not reported"))) %>% arrange(order) %>% select(-order) %>%
  mutate(PWA_pct = 100 * PWA / sum(d$Group == "PWA"), HC_pct = 100 * HC / sum(d$Group == "HC"))
ft <- fisher.test(table(d$Group, d$edu_level))
wt <- suppressWarnings(wilcox.test(edu_level ~ Group, d))
fhs <- fisher.test(table(d$Group, d$edu_level == 1))
education_tests <- data.frame(test = c("Fisher exact, four levels", "Wilcoxon rank-sum on ordinal level",
                                       "Fisher exact, high school only vs beyond"),
                              statistic = c(NA, unname(wt$statistic), NA), p = c(ft$p.value, wt$p.value, fhs$p.value),
                              n_PWA = sum(!is.na(d$edu_level[d$Group == "PWA"])),
                              n_HC = sum(!is.na(d$edu_level[d$Group == "HC"])))
sub <- d %>% filter(!is.na(edu_level))
covs <- sub %>% select(Participant, Age, edu_level) %>% mutate(Age_c = as.numeric(scale(Age)))
grab <- function(m, model, spec) {
  cf <- summary(m)$coefficients
  r <- cf[grep("GroupPWA", rownames(cf)), , drop = FALSE]
  data.frame(model = model, specification = spec, term = rownames(r), estimate = r[, 1], SE = r[, 2], p = r[, ncol(r)],
             row.names = NULL)
}
lmm_pair <- function(df, dv, label) {
  df <- df %>% mutate(Group = factor(Group, levels = c("HC", "PWA"))) %>% inner_join(covs, by = "Participant")
  m0 <- lmer(as.formula(paste(dv, "~ Group + Age + (1|Participant) + (1|Clip)")), df)
  m1 <- lmer(as.formula(paste(dv, "~ Group + Age + edu_level + (1|Participant) + (1|Clip)")), df)
  rbind(grab(m0, label, "baseline"), grab(m1, label, "+ education")) %>% mutate(n_participants = length(unique(df$Participant)))
}
sm <- entropy %>%
  mutate(Group = factor(Group, levels = c("HC", "PWA")), WindowBin = factor(WindowBin, levels = c("Short", "Medium", "Long"))) %>%
  inner_join(covs, by = "Participant")
mt <- tasks %>% mutate(Group = factor(Group, levels = c("HC", "PWA"))) %>% inner_join(covs, by = "Participant") %>%
  filter(!is.na(Accuracy))
glmm_pair <- function(task, label) {
  dt <- mt %>% filter(Task == task)
  g0 <- glmer(Accuracy ~ Group + Age_c + (1 | Stimulus) + (1 | Clip) + (1 | Participant), dt, family = binomial,
              control = glmerControl(optimizer = "bobyqa"))
  g1 <- glmer(Accuracy ~ Group + Age_c + edu_level + (1 | Stimulus) + (1 | Clip) + (1 | Participant), dt, family = binomial,
              control = glmerControl(optimizer = "bobyqa"))
  rbind(grab(g0, label, "baseline"), grab(g1, label, "+ education")) %>% mutate(n_participants = length(unique(dt$Participant)))
}
education_sensitivity <- bind_rows(
  rbind(grab(lm(PC1 ~ Group + Age, sub), "PC1 ~ Group + Age", "baseline"),
        grab(lm(PC1 ~ Group + Age + edu_level, sub), "PC1 ~ Group + Age", "+ education")) %>% mutate(n_participants = nrow(sub)),
  lmm_pair(holistic, "msz", "Holistic MSZ"), lmm_pair(realtime, "msz", "Real-time MSZ"),
  lmm_pair(isc, "isc_z", "Time-resolved ISC"),
  rbind(grab(lmer(SampEn ~ WindowBin * Group + Age + (1 | Participant) + (1 | Clip), sm), "Sample entropy", "baseline"),
        grab(lmer(SampEn ~ WindowBin * Group + Age + edu_level + (1 | Participant) + (1 | Clip), sm), "Sample entropy",
             "+ education")) %>%
    mutate(n_participants = length(unique(sm$Participant))),
  glmm_pair("CompQ", "Movie comprehension accuracy (GLMM)"), glmm_pair("AntGen", "Antonym generation accuracy (GLMM)"))
write_csv(education_counts, output_path("education_by_group.csv"))
write_csv(education_tests, output_path("education_tests.csv"))
write_csv(education_sensitivity, output_path("education_sensitivity.csv"))

# Table 1: sex and race (white vs non-white) by group, chi-square without and with continuity correction
sex <- table(d$Group, d$sex)
race <- d %>% select(Participant, Group) %>% left_join(dem %>% select(Participant, Race), by = "Participant")
stopifnot(nrow(race) == 100, !anyNA(race$Race))
race_white <- table(factor(race$Group, levels = c("PWA", "HC")),
                    factor(ifelse(race$Race == "White", "white", "nonwhite"), levels = c("white", "nonwhite")))
chi <- list(chisq.test(sex, correct = FALSE), chisq.test(sex), chisq.test(race_white, correct = FALSE), chisq.test(race_white))
sex_race <- data.frame(variable = rep(c("sex (self-reported)", "race white vs non-white (Table 1 check)"), each = 2),
                       correction = rep(c("none", "Yates"), 2), chisq = sapply(chi, function(x) unname(x$statistic)),
                       p = sapply(chi, function(x) x$p.value))
write_csv(as.data.frame(sex), output_path("sex_by_group.csv"))
write_csv(sex_race, output_path("sex_race_tests.csv"))

# Results 2.1, Methods 4.1, Table 1: the mood subsample (PHQ-8, GAD-7) and Welch t-tests
mood <- d %>% group_by(Group) %>%
  summarise(n = n(), n_PHQ = sum(!is.na(PHQ_Total)), n_GAD = sum(!is.na(GAD_Total)),
            PHQ_mean = mean(PHQ_Total, na.rm = TRUE), PHQ_sd = sd(PHQ_Total, na.rm = TRUE),
            PHQ_min = min(PHQ_Total, na.rm = TRUE), PHQ_max = max(PHQ_Total, na.rm = TRUE),
            GAD_mean = mean(GAD_Total, na.rm = TRUE), GAD_sd = sd(GAD_Total, na.rm = TRUE),
            GAD_min = min(GAD_Total, na.rm = TRUE), GAD_max = max(GAD_Total, na.rm = TRUE))
tp <- t.test(PHQ_Total ~ Group, d)
tg <- t.test(GAD_Total ~ Group, d)
mood_tests <- data.frame(test = c("PHQ-8 Welch t", "GAD-7 Welch t"), t = c(tp$statistic, tg$statistic),
                         p = c(tp$p.value, tg$p.value))
write_csv(mood, output_path("mood_subsample.csv"))
write_csv(mood_tests, output_path("mood_tests.csv"))

j <- joint %>% filter(model == "PC1 ~ AQ + PHQ-8")
k <- interaction %>% filter(model == "PC1 ~ PHQ-8 x Group + Age", term == "PHQ_Total:GroupPWA")
cat(sprintf("PC1 ~ AQ + PHQ-8 (PWA, n = %d): AQ b = %.3f (SE %.3f), PHQ-8 b = %.3f (SE %.3f, p = %.3f)\n",
            j$n[1], j$estimate[2], j$SE[2], j$estimate[3], j$SE[3], j$p[3]))
cat(sprintf("PHQ-8 x Group (age adjusted): b = %.3f (SE %.3f, p = %.3f)\n", k$estimate, k$SE, k$p))
cat(sprintf("holistic MSZ group effect: pooled reference %.3f, HC reference %.3f (SE %.3f, t = %.2f); r = %.3f\n",
            reference$estimate[1], reference$estimate[2], reference$SE[2], reference$t[2],
            reference_info$participant_level_r_pooled_vs_hcref))
cat(sprintf("education Fisher p = %.3f (four levels), %.3f (high school only)\n", ft$p.value, fhs$p.value))
cat(sprintf("chi-square with continuity correction: sex %.3f (p = %.3f), race %.3f (p = %.3f)\n",
            sex_race$chisq[2], sex_race$p[2], sex_race$chisq[4], sex_race$p[4]))
cat(sprintf("mood subsample: PHQ-8 %d PWA, %d HC (t = %.3f); GAD-7 %d PWA, %d HC (t = %.3f)\n", mood$n_PHQ[mood$Group == "PWA"],
            mood$n_PHQ[mood$Group == "HC"], abs(tp$statistic), mood$n_GAD[mood$Group == "PWA"], mood$n_GAD[mood$Group == "HC"],
            abs(tg$statistic)))
