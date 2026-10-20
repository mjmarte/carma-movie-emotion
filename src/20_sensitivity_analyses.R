# Comprehension sensitivity analyses of the four emotion features (Results 2.7, Supplementary Table S1): the Group effect
# without and with participant-clip comprehension accuracy as a covariate (Panel A), in the PWA participant-clips with all
# three comprehension items correct plus all HC participant-clips (Panel B), and on Partly Cloudy, the clip without dialogue
# (Panel C). Writes table_s1_sensitivity.csv and table_s1_comprehension_subset.csv (the participant-clip counts of Panel B).
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

ids <- read.csv(output_path("analytic_ids.csv"))
# Age as in Table S1: PWA from the aphasia quotient table, HC from the demographics table
age <- bind_rows(
  read.csv(input_file("aphasia_quotient")) %>% select(Participant, Age) %>%
    filter(Participant %in% ids$Participant[ids$Group == "PWA"]),
  read.csv(input_file("demographics")) %>% select(Participant, Age) %>%
    filter(Participant %in% ids$Participant[ids$Group == "HC"]))
stopifnot(nrow(age) == nrow(ids), !anyDuplicated(age$Participant), !anyNA(age$Age))

# Comprehension accuracy per participant-clip: share of the three items answered correctly (unanswered items count as
# incorrect)
compq <- readRDS(output_path("task_responses.rds")) %>%
  filter(Task == "CompQ") %>%
  group_by(Participant, Group, Clip) %>%
  summarise(compq_n = n(), compq_corr = sum(Accuracy, na.rm = TRUE), .groups = "drop") %>%
  mutate(compq_acc = compq_corr / compq_n)
high_comp <- compq %>% filter(Group == "HC" | compq_acc == 1)

# Feature clip labels that differ from those of the task responses. As in Table S1, the holistic ratings' label for No
# Country for Old Men is not mapped, so those participant-clips have no comprehension score in Panels A and B.
task_clip <- c(AkeelahAndTheBee = "Akeelah", CatchMe = "CatchMeIfYouCan")

features <- list(
  "Holistic MSZ" = list(file = "holistic_msz.rds", dv = "msz", fixed = "Group"),
  "Real-time MSZ" = list(file = "realtime_msz.rds", dv = "msz", fixed = "Group"),
  "Time-resolved ISC" = list(file = "isc.rds", dv = "isc_z", fixed = "Group"),
  "Sample entropy" = list(file = "sample_entropy.rds", dv = "SampEn", fixed = "WindowBin * Group"))

group_terms <- function(m, terms) {
  cf <- coef(summary(m))
  rows <- grep(terms, rownames(cf))
  data.frame(term = rownames(cf)[rows], estimate = cf[rows, "Estimate"], SE = cf[rows, "Std. Error"],
             df = if ("df" %in% colnames(cf)) cf[rows, "df"] else m$df.residual,
             t = cf[rows, "t value"], p = cf[rows, "Pr(>|t|)"], row.names = NULL)
}
fit <- function(d, formula, feature, panel, specification, type = "LMM", terms = "^GroupPWA$", control = lmerControl()) {
  m <- if (type == "LM") lm(as.formula(formula), d) else lmer(as.formula(formula), d, control = control)
  stopifnot(nobs(m) == nrow(d))
  data.frame(panel = panel, feature = feature, specification = specification, model = type, group_terms(m, terms),
             n_obs = nrow(d), n_PWA = n_distinct(d$Participant[d$Group == "PWA"]),
             n_HC = n_distinct(d$Participant[d$Group == "HC"]))
}

results <- list()
for (feature in names(features)) {
  f <- features[[feature]]
  d <- readRDS(output_path(f$file)) %>%
    mutate(Group = factor(Group, levels = c("HC", "PWA")),
           task_clip = ifelse(Clip %in% names(task_clip), unname(task_clip[Clip]), Clip)) %>%
    left_join(age, by = "Participant")
  re <- "+ (1 | Participant) + (1 | Clip)"

  # A. Comprehension as a covariate; the baseline is refit on the participant-clips with a comprehension score
  da <- d %>%
    left_join(compq %>% select(Participant, task_clip = Clip, compq_acc), by = c("Participant", "task_clip")) %>%
    filter(!is.na(compq_acc))
  results[[length(results) + 1]] <- fit(da, paste(f$dv, "~", f$fixed, "+ Age", re), feature, "A", "Baseline (subsample)")
  results[[length(results) + 1]] <- fit(da, paste(f$dv, "~", f$fixed, "+ Age + compq_acc", re), feature, "A",
                                        "+ Comprehension")

  # B. PWA participant-clips with all comprehension items correct, all HC participant-clips
  db <- d %>% semi_join(high_comp %>% select(Participant, task_clip = Clip), by = c("Participant", "task_clip"))
  results[[length(results) + 1]] <- fit(db, paste(f$dv, "~", f$fixed, "+ Age", re), feature, "B",
                                        "PWA participant-clips with full comprehension")

  # C. Partly Cloudy only: linear models at the participant-clip level, participant random intercept at the window level
  dc <- d %>% filter(Clip == "PartlyCloudy")
  if (feature %in% c("Holistic MSZ", "Real-time MSZ")) {
    results[[length(results) + 1]] <- fit(dc, paste(f$dv, "~ Group + Age"), feature, "C", "Partly Cloudy only", "LM")
  } else {
    results[[length(results) + 1]] <- fit(dc, paste(f$dv, "~", f$fixed, "+ Age + (1 | Participant)"), feature, "C",
                                          "Partly Cloudy only", terms = "GroupPWA",
                                          control = lmerControl(optimizer = "bobyqa"))
  }
}
table_s1 <- bind_rows(results) %>% arrange(panel)

subset_counts <- compq %>%
  group_by(Group) %>%
  summarise(participant_clips = n(), participant_clips_retained = sum(Group == "HC" | compq_acc == 1),
            participants = n_distinct(Participant),
            participants_retained = n_distinct(Participant[Group == "HC" | compq_acc == 1]), .groups = "drop") %>%
  mutate(percent_retained = 100 * participant_clips_retained / participant_clips)

write.csv(table_s1, output_path("table_s1_sensitivity.csv"), row.names = FALSE)
write.csv(subset_counts, output_path("table_s1_comprehension_subset.csv"), row.names = FALSE)
print(table_s1 %>% mutate(across(c(estimate, SE, df, t), ~ round(., 3)), p = signif(p, 3)), row.names = FALSE)
print(subset_counts %>% mutate(percent_retained = round(percent_retained, 1)), row.names = FALSE)
