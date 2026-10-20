# Item-level movie comprehension (CompQ) and antonym generation (AntGen) accuracy of the analytic sample. Reproduces Results 2.3
# (accuracy by group, pooled over items) and the Methods 4.5 binomial GLMMs, Accuracy ~ Group + Age + (1|Stimulus) + (1|Clip) +
# (1|Participant), with age centred on the sample mean. Writes task_responses.rds (read by the later scripts) and
# task_accuracy_results.csv.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
})

ids <- read.csv(output_path("analytic_ids.csv"))
tasks <- read.csv(input_file("task_responses"))
for (r in RULES$relabel_participants) tasks$Participant[tasks$Participant %in% r$from] <- r$to
tasks <- tasks %>%
  filter(Participant %in% ids$Participant) %>%
  left_join(ids %>% rename(sample_group = Group), by = "Participant") %>%
  mutate(Group = sample_group) %>%
  select(-sample_group)
stopifnot(n_distinct(tasks$Participant) == nrow(ids))
saveRDS(tasks, output_path("task_responses.rds"))

accuracy <- tasks %>%
  group_by(Group, Task) %>%
  summarise(n_participants = n_distinct(Participant), n_items = sum(!is.na(Accuracy)),
            mean_accuracy = mean(Accuracy, na.rm = TRUE), sd_accuracy = sd(Accuracy, na.rm = TRUE), .groups = "drop")

age <- ids %>%
  left_join(read.csv(input_file("demographics")) %>% select(Participant, Age), by = "Participant") %>%
  mutate(Age_centered = Age - mean(Age, na.rm = TRUE))
stopifnot(!anyNA(age$Age))
fit_task <- function(task_name) {
  task_data <- tasks %>%
    left_join(age %>% select(Participant, Age, Age_centered), by = "Participant") %>%
    filter(Task == task_name)
  model <- glmer(Accuracy ~ Group + Age_centered + (1 | Stimulus) + (1 | Clip) + (1 | Participant),
                 data = task_data, family = binomial)
  co <- summary(model)$coefficients
  data.frame(Task = task_name, term = rownames(co), n_participants = unname(ngrps(model)["Participant"]),
             n_items = nobs(model), estimate = co[, "Estimate"], se = co[, "Std. Error"], z = co[, "z value"],
             p = co[, "Pr(>|z|)"], row.names = NULL)
}
glmms <- bind_rows(fit_task("CompQ"), fit_task("AntGen"))

results <- bind_rows(accuracy %>% mutate(result = "accuracy"), glmms %>% mutate(result = "glmm")) %>%
  select(result, Task, Group, term, n_participants, n_items, mean_accuracy, sd_accuracy, estimate, se, z, p)
write.csv(results, output_path("task_accuracy_results.csv"), row.names = FALSE)
cat("task responses:", nrow(tasks), "rows,", n_distinct(tasks$Participant), "participants\n")
for (k in c("CompQ", "AntGen")) {
  a <- accuracy %>% filter(Task == k)
  g <- glmms %>% filter(Task == k, term == "GroupPWA")
  cat(sprintf("%s: HC %.1f%%, PWA %.1f%%; Group (PWA) b = %.3f, SE = %.3f, z = %.3f, p = %.2g\n", k,
              100 * a$mean_accuracy[a$Group == "HC"], 100 * a$mean_accuracy[a$Group == "PWA"], g$estimate, g$se, g$z, g$p))
}
