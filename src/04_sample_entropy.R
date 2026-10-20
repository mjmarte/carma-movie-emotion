# Sample entropy of the real-time valence ratings (Methods 4.4): pracma::sample_entropy with m = 2 and r = 0.15 rating
# units, in non-overlapping windows of 5 to 35 s per participant and clip, binned into short (5-15 s), medium (20-25 s) and
# long (30-35 s) scales. Feeds the entropy models of Results 2.3 and the medium-scale feature of the PCA, SVM and SVR.
# Writes sample_entropy.rds.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(pracma)
  library(foreach)
  library(doParallel)
})

sample_entropy_windows <- function(data, clip, window_sizes = c(5, 10, 15, 20, 25, 30, 35), r = 0.15) {
  clip_data <- data %>%
    filter(Clip == clip)
  participants <- unique(clip_data$Participant)
  max_second <- max(clip_data$Second)
  param_combinations <- expand.grid(
    participant = participants,
    window_length = window_sizes,
    stringsAsFactors = FALSE
  )
  foreach(i = 1:nrow(param_combinations), .combine = rbind, .packages = c("dplyr", "pracma")) %dopar% {
    participant <- param_combinations$participant[i]
    window_length <- param_combinations$window_length[i]
    participant_data <- clip_data %>%
      filter(Participant == participant) %>%
      arrange(Second)
    if (nrow(participant_data) < window_length) {
      return(NULL)
    }
    participant_results <- data.frame(
      Clip = character(),
      Participant = character(),
      Group = character(),
      StartTime = numeric(),
      WindowSize = numeric(),
      SampEn = numeric(),
      stringsAsFactors = FALSE
    )
    # Non-overlapping windows: the step equals the window length
    step_size <- window_length
    for (start_time in seq(1, max_second - window_length + 1, by = step_size)) {
      end_time <- start_time + window_length - 1
      window_data <- participant_data %>%
        filter(Second >= start_time & Second <= end_time) %>%
        pull(Rating)
      # Windows with ratings for fewer than 80% of their seconds are skipped; undefined or infinite values are dropped
      if (length(window_data) >= window_length * 0.8) {
        sampleen <- tryCatch({
          sample_entropy(window_data, edim = 2, r = r)
        }, error = function(e) {
          return(NA)
        })
        if (!is.na(sampleen) && !is.infinite(sampleen)) {
          participant_results <- rbind(participant_results, data.frame(
            Clip = clip,
            Participant = participant,
            Group = unique(participant_data$Group)[1],
            StartTime = start_time,
            WindowSize = window_length,
            SampEn = sampleen,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    return(participant_results)
  }
}

ratings <- readRDS(output_path("realtime_ratings.rds"))
window_sizes <- c(5, 10, 15, 20, 25, 30, 35)

# Tasks run in parallel and are combined in task order, so the result does not depend on the number of workers
cl <- makeCluster(max(1, min(8, parallel::detectCores() - 1), na.rm = TRUE))
registerDoParallel(cl)
sampleen_results <- data.frame()
for (clip in unique(ratings$Clip)) {
  clip_results <- sample_entropy_windows(ratings, clip, window_sizes = window_sizes, r = 0.15)
  sampleen_results <- rbind(sampleen_results, clip_results)
}
stopCluster(cl)

sampleen_binned <- sampleen_results %>%
  mutate(
    WindowBin = case_when(
      WindowSize %in% c(5, 10, 15) ~ "Short",
      WindowSize %in% c(20, 25) ~ "Medium",
      WindowSize %in% c(30, 35) ~ "Long"
    ),
    WindowBin = factor(WindowBin, levels = c("Short", "Medium", "Long"))
  ) %>%
  filter(!is.na(WindowBin))
saveRDS(sampleen_binned, output_path("sample_entropy.rds"))
cat("sample entropy:", nrow(sampleen_binned), "windows,", n_distinct(sampleen_binned$Participant), "participants,",
    n_distinct(sampleen_binned$Clip), "clips\n")
