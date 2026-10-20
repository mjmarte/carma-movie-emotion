# Table 2 emotional response columns and Results 2.2, as published: per clip, mean (SD) over participant-clips of the mean
# PANAS-SF rating of the five positive and of the five negative items, over all participants with PANAS ratings in the
# item-level file, and mean (SD) of the 1 Hz real-time valence ratings over the analytic sample. Writes table2.csv.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

panas_abbrev <- c(Akeelah = "AKB", CatchMe = "CMIYC", GoodWillHunting = "GWH", Miracle = "MIR", Moonlight = "MOON",
                  `No Country for Old Men` = "NCOM", PartlyCloudy = "PC", ParentTrap = "PT")
panas <- read.csv(input_file("panas_items")) %>%
  group_by(participant_id, clip, valence) %>%
  summarise(m = mean(rating, na.rm = TRUE), .groups = "drop") %>%
  mutate(acr = panas_abbrev[clip]) %>%
  group_by(acr, valence) %>%
  summarise(mean = mean(m), sd = sd(m), n = n(), .groups = "drop")
realtime <- readRDS(output_path("realtime_ratings.rds")) %>%
  mutate(acr = as.character(clip_abbrev)) %>%
  group_by(acr) %>%
  summarise(mean = mean(Rating), sd = sd(Rating), n_part = n_distinct(Participant), .groups = "drop")
table2 <- panas %>%
  pivot_wider(names_from = valence, values_from = c(mean, sd, n)) %>%
  inner_join(realtime, by = "acr") %>%
  transmute(acr, positive = sprintf("%.2f (%.2f)", mean_positive, sd_positive),
            negative = sprintf("%.2f (%.2f)", mean_negative, sd_negative), valence = sprintf("%.2f (%.2f)", mean, sd),
            n_panas_participant_clips = n_positive, n_realtime_participants = n_part, mean_positive, mean_negative,
            mean_valence = mean)
stopifnot(nrow(table2) == 8, !anyNA(table2))
write.csv(table2, output_path("table2.csv"), row.names = FALSE)

with(table2, cat(sprintf("%-6s positive %s, negative %s, real-time valence %s\n", acr, positive, negative, valence), sep = ""))
cat("highest positive:", table2$acr[which.max(table2$mean_positive)], "; highest negative:",
    table2$acr[which.max(table2$mean_negative)], "; most positive valence:", table2$acr[which.max(table2$mean_valence)],
    "; most negative valence:", table2$acr[which.min(table2$mean_valence)], "\n")
