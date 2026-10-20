# Robustness of the PC correlations to influential observations (Methods 4.6): for PC1-PC3 with AQ in PWA and PC1 with
# PHQ-8 and GAD-7 in each group, observations flagged by at least two of Cook's distance > 4/n, |studentized residual| > 2
# and leverage > 4/n (simple regression of the score on the PC) are removed and the Pearson correlation is recomputed.
# Writes outlier_checks.csv.
source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

ids <- read.csv(output_path("analytic_ids.csv"))
data <- read.csv(output_path("pc_scores.csv")) %>% select(Participant, PC1, PC2, PC3) %>%
  inner_join(ids, by = "Participant") %>%
  left_join(read.csv(input_file("aphasia_quotient")) %>% select(Participant, AQ), by = "Participant") %>%
  left_join(read.csv(input_file("demographics")) %>% select(Participant, PHQ_Total, GAD_Total), by = "Participant")
stopifnot(nrow(data) == 100)
check <- function(group, pc, score) {
  d <- data %>% filter(Group == group, !is.na(.data[[score]]))
  model <- lm(as.formula(paste(score, "~", pc)), data = d)
  n <- nrow(d)
  flags <- (cooks.distance(model) > 4 / n) + (abs(rstudent(model)) > 2) + (hatvalues(model) > 4 / n)
  keep <- flags < 2
  full <- cor.test(d[[pc]], d[[score]]); clean <- cor.test(d[[pc]][keep], d[[score]][keep])
  data.frame(group = group, pc = pc, score = score, n = n, n_flagged = sum(!keep), r = unname(full$estimate), p = full$p.value,
             r_without_flagged = unname(clean$estimate), p_without_flagged = clean$p.value)
}
results <- bind_rows(lapply(c("PC1", "PC2", "PC3"), function(pc) check("PWA", pc, "AQ")),
                     lapply(c("PWA", "HC"), function(g) bind_rows(check(g, "PC1", "PHQ_Total"), check(g, "PC1", "GAD_Total"))))
write.csv(results, output_path("outlier_checks.csv"), row.names = FALSE)
with(results, cat(sprintf("%-3s %s ~ %-9s n = %2d, flagged %d: r = %.3f (p = %.3g); without flagged r = %.3f (p = %.3g)\n",
                          group, pc, score, n, n_flagged, r, p, r_without_flagged, p_without_flagged), sep = ""))
