# Time-resolved intersubject correlation of the real-time valence ratings (Methods 4.4): in 15 s windows moved in 5 s steps, each
# participant's ratings are correlated with the HC mean trajectory (the mean of the other HC for HC participants) and Fisher
# z-transformed; windows with fewer than 10 points or zero variance are skipped. Input to the ISC group effects (Results 2.3 and
# 2.7, Tables S1 and S3) and the ISC features of the PCA (Results 2.4). Writes isc.rds (one row per participant, clip and window).
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

compute_dynamic_isc <- function(data, window_size = 15, step = 5) {
  hc_mean <- data %>%
    filter(Group == "HC") %>%
    group_by(Clip, Second) %>%
    summarise(mean_rating = mean(Rating, na.rm = TRUE), .groups = "drop")

  data %>%
    group_by(Participant, Group, Clip) %>%
    group_modify(function(participant_data, keys) {
      clip_hc <- hc_mean %>%
        filter(Clip == keys$Clip) %>%
        arrange(Second)
      # Leave-one-out reference for HC participants
      if (keys$Group == "HC") {
        clip_hc <- data %>%
          filter(Group == "HC", Clip == keys$Clip, Participant != keys$Participant) %>%
          group_by(Second) %>%
          summarise(mean_rating = mean(Rating, na.rm = TRUE), .groups = "drop") %>%
          arrange(Second)
      }

      max_second <- max(participant_data$Second, na.rm = TRUE)
      windows <- seq(1, max_second - window_size + 1, by = step)

      map_dfr(windows, function(start) {
        end <- start + window_size - 1
        window_participant <- participant_data %>%
          filter(Second >= start & Second <= end) %>%
          arrange(Second)
        window_hc <- clip_hc %>%
          filter(Second >= start & Second <= end) %>%
          arrange(Second)
        if (nrow(window_participant) < 10 || nrow(window_hc) < 10) return(NULL)

        merged <- window_participant %>%
          inner_join(window_hc, by = "Second")
        if (nrow(merged) < 10) return(NULL)
        if (sd(merged$Rating) == 0 || sd(merged$mean_rating) == 0) return(NULL)

        isc_val <- cor(merged$Rating, merged$mean_rating, use = "complete.obs")
        tibble(
          window_start = start,
          window_center = start + window_size / 2,
          isc = isc_val,
          isc_z = atanh(isc_val),
          n_points = nrow(merged)
        )
      })
    }) %>%
    ungroup()
}

ratings <- readRDS(output_path("realtime_ratings.rds"))
isc <- compute_dynamic_isc(ratings, window_size = 15, step = 5)
saveRDS(isc, output_path("isc.rds"))
cat("time-resolved ISC:", nrow(isc), "windows from", nrow(distinct(isc, Participant, Clip)), "participant-clips\n")
