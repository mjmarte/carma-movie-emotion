# Real-time MSZ (Methods 4.4): each participant's valence rating at each second of a clip z-scored against the HC mean and SD at
# that second (the other HC for HC participants; SD floor 0.5), squared, and averaged over the seconds of each participant-clip
# (realtime_msz.rds: the PCA scores, Figures 4 and 5, Results 2.5 to 2.7, Tables S1 and S3, the SVM and SVR). The real-time MSZ
# group effect in Results 2.3, the PCA statistics in Results 2.4 and Supplementary Figure S2 use the same measure referenced to the
# mean and SD of all participants at each second, without an SD floor (realtime_msz_pooled.rds).
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

ratings <- readRDS(output_path("realtime_ratings.rds"))
MIN_SD <- 0.5

hc_stats <- ratings %>%
  filter(Group == "HC") %>%
  group_by(Clip, Second) %>%
  summarise(HC_mean = mean(Rating, na.rm = TRUE), HC_sd = sd(Rating, na.rm = TRUE), n_hc = n(), .groups = "drop")

clip_sds <- ratings %>%
  filter(Group == "HC") %>%
  group_by(Clip) %>%
  summarise(clip_sd = sd(Rating, na.rm = TRUE), .groups = "drop")

# Leave-one-out reference for each HC participant
hc_participants <- unique(ratings$Participant[ratings$Group == "HC"])
hc_loo_stats <- list()
for (hc_part in hc_participants) {
  hc_loo_stats[[hc_part]] <- ratings %>%
    filter(Group == "HC", Participant != hc_part) %>%
    group_by(Clip, Second) %>%
    summarise(HC_mean = mean(Rating, na.rm = TRUE), HC_sd = sd(Rating, na.rm = TRUE), n_hc = n(), .groups = "drop") %>%
    left_join(clip_sds, by = "Clip")
}

realtime_msz <- ratings %>%
  group_by(Participant, Group, Clip) %>%
  group_split() %>%
  map_dfr(function(part_data) {
    participant <- part_data$Participant[1]
    group <- part_data$Group[1]
    clip <- part_data$Clip[1]

    if (group == "HC") {
      reference <- hc_loo_stats[[participant]] %>%
        filter(Clip == clip)
    } else {
      reference <- hc_stats %>%
        filter(Clip == clip) %>%
        left_join(clip_sds, by = "Clip")
    }

    z_data <- part_data %>%
      left_join(reference, by = c("Clip", "Second")) %>%
      mutate(
        sd_to_use = pmax(HC_sd, MIN_SD, na.rm = TRUE),
        sd_to_use = ifelse(is.na(sd_to_use), clip_sd, sd_to_use),
        z_score = (Rating - HC_mean) / sd_to_use,
        z_squared = z_score^2,
        sd_type = case_when(
          is.na(HC_sd) ~ "clip",
          HC_sd < MIN_SD ~ "minimum",
          TRUE ~ "actual"
        )
      ) %>%
      filter(!is.na(z_score))

    if (nrow(z_data) > 0) {
      tibble(
        Participant = participant,
        Group = group,
        Clip = clip,
        msz = mean(z_data$z_squared, na.rm = TRUE),
        n_timepoints = nrow(z_data),
        n_excluded = nrow(part_data) - nrow(z_data),
        mean_z = mean(z_data$z_score, na.rm = TRUE),
        sd_z = sd(z_data$z_score, na.rm = TRUE),
        max_abs_z = max(abs(z_data$z_score), na.rm = TRUE),
        pct_min_sd = mean(z_data$sd_type != "actual") * 100
      )
    } else {
      tibble()
    }
  })

saveRDS(realtime_msz, output_path("realtime_msz.rds"))

pooled_reference <- ratings %>%
  group_by(Clip, Second) %>%
  summarise(pooled_mean = mean(Rating), pooled_sd = sd(Rating), .groups = "drop")
realtime_msz_pooled <- ratings %>%
  left_join(pooled_reference, by = c("Clip", "Second")) %>%
  mutate(z_score = (Rating - pooled_mean) / pooled_sd) %>%
  filter(is.finite(z_score)) %>%
  group_by(Participant, Group, Clip) %>%
  summarise(msz = mean(z_score^2), n_timepoints = n(), .groups = "drop")
saveRDS(realtime_msz_pooled, output_path("realtime_msz_pooled.rds"))
cat("real-time MSZ:", nrow(realtime_msz), "participant-clips (", sum(realtime_msz$Group == "PWA"), "PWA,",
    sum(realtime_msz$Group == "HC"), "HC )\n")
