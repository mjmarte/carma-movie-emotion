# Prior familiarity with the clips by group (Supplementary Table S3; Methods 4.5): percentage of PWA and HC who had seen each
# clip before, Pearson chi-square test without continuity correction (as reported), the continuity-corrected test and Fisher's
# exact test, each also Bonferroni-corrected over the eight clips. Writes table_s3.csv.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

ids <- read.csv(output_path("analytic_ids.csv"))
stopifnot(nrow(ids) == 100, sum(ids$Group == "PWA") == 57)
# Column of the export holding each clip's "seen this clip before" item (1 = had seen the clip, 2 = had not)
clips <- c(aatb_seen_before_v2 = "Akeelah and the Bee", cmiyc_seen_before_v2 = "Catch Me If You Can",
           gwh_seen_before_v2 = "Good Will Hunting", mi_seen_before_v2 = "Miracle", m_seen_before_v2 = "Moonlight",
           ncom_seen_before_v2 = "No Country for Old Men", dos_seen_before_v2 = "Partly Cloudy",
           pt_seen_clip_v2 = "The Parent Trap")
export <- input_file("familiarity")
header <- names(read_excel(export, sheet = "Data", n_max = 0))
seen <- read_excel(export, sheet = "Data", col_types = ifelse(header %in% c("participant_id", names(clips)), "guess", "skip")) %>%
  rename(Participant = participant_id) %>%
  filter(Participant %in% ids$Participant)
stopifnot(!anyDuplicated(seen$Participant))
d <- ids %>% left_join(seen, by = "Participant")
stopifnot(nrow(d) == nrow(ids))

familiarity <- bind_rows(lapply(names(clips), function(cc) {
  v <- d[[cc]]
  stopifnot(all(v %in% c(1, 2, NA)))
  k <- d %>%
    transmute(Group = factor(Group, c("PWA", "HC")), seen = factor(ifelse(v == 1, "seen", "not seen"), c("seen", "not seen"))) %>%
    filter(!is.na(seen))
  tb <- table(k$Group, k$seen)
  p0 <- suppressWarnings(chisq.test(tb, correct = FALSE))
  p1 <- suppressWarnings(chisq.test(tb, correct = TRUE))
  data.frame(clip = clips[[cc]], n_PWA = sum(tb["PWA", ]), seen_PWA = tb["PWA", "seen"], n_HC = sum(tb["HC", ]),
             seen_HC = tb["HC", "seen"], pct_PWA = 100 * tb["PWA", "seen"] / sum(tb["PWA", ]),
             pct_HC = 100 * tb["HC", "seen"] / sum(tb["HC", ]), min_expected = min(p0$expected),
             chisq_pearson = unname(p0$statistic), p_pearson = p0$p.value,
             chisq_yates = unname(p1$statistic), p_yates = p1$p.value,
             p_fisher = fisher.test(tb)$p.value)
})) %>%
  mutate(p_pearson_bonf = pmin(1, 8 * p_pearson), p_yates_bonf = pmin(1, 8 * p_yates), p_fisher_bonf = pmin(1, 8 * p_fisher))
write.csv(familiarity, output_path("table_s3.csv"), row.names = FALSE)

with(familiarity, cat(sprintf(paste0("%-23s PWA %d/%d (%.1f%%), HC %d/%d (%.1f%%); chi2 = %.2f, p = %.3f (Bonferroni %.3f); ",
                                     "Fisher p = %.3f\n"), clip, seen_PWA, n_PWA, pct_PWA, seen_HC, n_HC, pct_HC, chisq_pearson,
                              p_pearson, p_pearson_bonf, p_fisher), sep = ""))
