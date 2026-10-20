# Supplementary Table S2: within-PWA Pearson correlations of the five emotion features and the two language accuracies with
# the TMT B/A ratio and the TEA Elevator Strings scaled score, with Benjamini-Hochberg correction within the ten emotion tests
# and within the four language tests. Reads the correlations written by 16_supplementary.R. Writes table_s2.csv.
source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

measures <- c("TMT B/A ratio", "TEA Elevator Strings")
emotion <- read.csv(output_path("cognitive_load_correlations.csv")) %>%
  filter(measure %in% measures) %>%
  transmute(domain = "emotion", feature, measure, n, r, p)
language <- read.csv(output_path("cognitive_language_correlations.csv")) %>%
  filter(feature != "PC1", measure %in% measures) %>%
  transmute(domain = "language", feature, measure, n, r, p)
table_s2 <- bind_rows(emotion, language) %>%
  group_by(domain) %>%
  mutate(p_bh = p.adjust(p, "BH")) %>%
  ungroup() %>%
  mutate(cell = sprintf("%.3f (%s)", r, sub("^0", "", sprintf("%.3f", p_bh))))
stopifnot(nrow(table_s2) == 14)
write.csv(table_s2, output_path("table_s2.csv"), row.names = FALSE)
with(table_s2, cat(sprintf("%-30s %-21s %s (n = %d)\n", feature, measure, cell, n), sep = ""))
