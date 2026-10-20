# Holistic MSZ (Methods 4.4): each post-clip PANAS-SF rating z-scored against the mean and SD of all participants for that
# clip and item (pooled reference), squared, and averaged over the 10 items of each participant-clip.
# Writes holistic_msz.rds and holistic_msz.csv.
source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

ids <- read.csv(output_path("analytic_ids.csv"))
items <- read.csv(input_file("panas_items"))
holistic_msz <- items %>%
  rename(Participant = participant_id, Clip = clip, Rating = rating, PANAS_Item = panas_item) %>%
  select(Participant, Clip, PANAS_Item, Rating) %>%
  inner_join(ids, by = "Participant") %>%
  select(Participant, Group, Clip, PANAS_Item, Rating) %>%
  mutate(.by = c(Clip, PANAS_Item),
         z_score = (Rating - mean(Rating, na.rm = TRUE)) / sd(Rating, na.rm = TRUE),
         z_squared = z_score^2) %>%
  group_by(Participant, Group, Clip) %>%
  summarise(msz = mean(z_squared, na.rm = TRUE), n_items = n(), .groups = "drop")
saveRDS(holistic_msz, output_path("holistic_msz.rds"))
write.csv(holistic_msz, output_path("holistic_msz.csv"), row.names = FALSE)
cat("holistic MSZ:", nrow(holistic_msz), "participant-clips\n")
