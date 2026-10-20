# Real-time valence ratings of the analytic sample (57 PWA, 43 HC; Methods 4.1 and 4.4). Applies the participant rules in
# config.json (groups and participants outside the analytic sample, corrected IDs) and orders the clips by mean valence.
# Writes realtime_ratings.rds (read by the later scripts) and analytic_ids.csv.
source("R/config.R")
suppressPackageStartupMessages(library(dplyr))

ratings <- read.csv(input_file("realtime_ratings")) %>%
  filter(!(Group %in% unlist(RULES$exclude_groups))) %>%
  filter(!(Participant %in% unlist(RULES$exclude_participants))) %>%
  mutate(Group = if_else(Group == "AP", "PWA", Group)) %>%
  filter(Group %in% c("HC", "PWA"))
for (r in RULES$relabel_participants) {
  hit <- ratings$Participant == r$from
  ratings$Participant[hit] <- r$to
  ratings$Group[ratings$Participant == r$to] <- r$group
}
abbrev <- c(AkeelahAndTheBee = "AKB", CatchMeIfYouCan = "CMIYC", NoCountryForOldMen = "NCOM", ParentTrap = "PT",
            PartlyCloudy = "PC", GoodWillHunting = "GWH", Miracle = "MIR", Moonlight = "MOON")
ratings <- ratings %>% mutate(clip_abbrev = unname(ifelse(Clip %in% names(abbrev), abbrev[Clip], Clip)))
ratings$Second <- as.numeric(ratings$Second)
ratings$Timepoint <- ratings$Second
valence_order <- ratings %>% group_by(clip_abbrev) %>% summarise(m = mean(Rating, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(clip_abbrev)
ratings <- ratings %>% mutate(clip_abbrev = factor(clip_abbrev, levels = valence_order))

ids <- ratings %>% distinct(Participant, Group) %>% filter(Group %in% c("HC", "PWA"))
stopifnot(nrow(ids) == 100, sum(ids$Group == "PWA") == 57, sum(ids$Group == "HC") == 43)
saveRDS(ratings, output_path("realtime_ratings.rds"))
write.csv(ids, output_path("analytic_ids.csv"), row.names = FALSE)
cat("real-time ratings:", nrow(ratings), "rows,", nrow(ids), "participants (", sum(ids$Group == "PWA"), "PWA,",
    sum(ids$Group == "HC"), "HC )\n")
