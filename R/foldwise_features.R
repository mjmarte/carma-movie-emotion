# Fold-wise recomputation of the reference-dependent movie features for the classification and regression models (Methods 4.7
# and 4.8). Real-time MSZ (SD floor 0.5) and time-resolved ISC (15 s windows, 5 s step) are referenced to the training HC, and
# training HC are referenced leave-one-out. Holistic MSZ is referenced to all training participants. Sourced by src/17_svm.R
# and src/18_svr.R after R/config.R. With all participants in training, validate_features() reproduces the PCA features.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

IDS <- read.csv(output_path("analytic_ids.csv")) %>% arrange(Participant)
stopifnot(all(IDS$Group %in% c("HC", "PWA")))
ALL_IDS <- IDS$Participant
GROUP <- setNames(IDS$Group, IDS$Participant)
HC_IDS <- ALL_IDS[GROUP == "HC"]
PWA_IDS <- ALL_IDS[GROUP == "PWA"]
MIN_SD <- 0.5; WIN <- 15; STEP <- 5
carma <- readRDS(output_path("realtime_ratings.rds")) %>% filter(Participant %in% ALL_IDS) %>%
  select(Participant, Clip, Second, Rating)
CLIPS <- sort(unique(carma$Clip))
W <- lapply(CLIPS, function(cl) {
  d <- carma %>% filter(Clip == cl); nsec <- max(d$Second)
  m <- matrix(NA_real_, nrow = nsec, ncol = length(ALL_IDS), dimnames = list(seq_len(nsec), ALL_IDS))
  m[cbind(d$Second, match(d$Participant, ALL_IDS))] <- d$Rating
  m
})
names(W) <- CLIPS
PC_MAXSEC <- carma %>% group_by(Participant, Clip) %>% summarise(max_sec = max(Second), .groups = "drop")
PANAS <- read_csv(input_file("panas_items"), show_col_types = FALSE) %>%
  rename(Participant = participant_id, Clip = clip, Rating = rating, PANAS_Item = panas_item) %>%
  filter(Participant %in% ALL_IDS) %>% select(Participant, Clip, PANAS_Item, Rating)
SAMPEN_FEAT <- readRDS(output_path("sample_entropy.rds")) %>% filter(Participant %in% ALL_IDS, WindowBin == "Medium") %>%
  group_by(Participant) %>% summarise(SampEn_Medium = mean(SampEn, na.rm = TRUE), .groups = "drop")
TASK_FEAT <- readRDS(output_path("task_responses.rds")) %>% filter(Participant %in% ALL_IDS) %>%
  group_by(Participant, Task) %>% summarise(acc = mean(Accuracy, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Task, values_from = acc, names_glue = "{Task}_mean")
AQ_CORRECT <- read_csv(input_file("aphasia_quotient"), show_col_types = FALSE) %>% select(Participant, AQ)
FEATS_LANG <- c("AntGen_mean", "CompQ_mean")
FEATS_EMO <- c("SampEn_Medium", "PANAS_MSZ_mean", "CARMA_MSZ_mean", "ISC_mean", "ISC_sd")
FEATS_ALL <- c(FEATS_EMO, FEATS_LANG)
FEATURE_SETS <- list(language = FEATS_LANG, emotion = FEATS_EMO, combined = FEATS_ALL, antgen_only = "AntGen_mean")

# Running sums per timepoint, so that a reference without one participant is a subtraction
ref_stats <- function(m, cols) {
  sub <- m[, cols, drop = FALSE]
  list(s = rowSums(sub, na.rm = TRUE), ss = rowSums(sub^2, na.rm = TRUE), cnt = rowSums(!is.na(sub)))
}
loo_adjust <- function(st, x) {
  pr <- !is.na(x); xx <- ifelse(pr, x, 0)
  list(s = st$s - xx, ss = st$ss - xx^2, cnt = st$cnt - pr)
}
mean_sd_from <- function(st) {
  mu <- ifelse(st$cnt > 0, st$s / st$cnt, NA_real_)
  v <- ifelse(st$cnt > 1, (st$ss - st$cnt * mu^2) / (st$cnt - 1), NA_real_)
  v[!is.na(v) & v < 0] <- 0
  list(mean = mu, sd = sqrt(v))
}
rt_msz_one <- function(x, ref) {
  ok <- !is.na(x) & !is.na(ref$mean); sd_use <- pmax(ref$sd[ok], MIN_SD, na.rm = TRUE)
  z <- (x[ok] - ref$mean[ok]) / sd_use
  mean(z^2)
}
isc_one <- function(x, refmean, max_sec) {
  out <- numeric(0); if (max_sec < WIN) return(out)
  for (st in seq(1, max_sec - WIN + 1, by = STEP)) {
    idx <- st:(st + WIN - 1); xp <- x[idx]; rp <- refmean[idx]
    if (sum(!is.na(xp)) < 10 || sum(!is.na(rp)) < 10) next
    both <- !is.na(xp) & !is.na(rp); if (sum(both) < 10) next
    a <- xp[both]; b <- rp[both]; if (sd(a) == 0 || sd(b) == 0) next
    out <- c(out, atanh(cor(a, b)))
  }
  out
}
realtime_features <- function(train_hc) {
  rt <- vector("list", nrow(PC_MAXSEC)); isc <- vector("list", nrow(PC_MAXSEC)); k <- 0
  for (cl in CLIPS) {
    m <- W[[cl]]; st_full <- ref_stats(m, train_hc); full <- mean_sd_from(st_full)
    parts <- PC_MAXSEC[PC_MAXSEC$Clip == cl, ]
    for (j in seq_len(nrow(parts))) {
      p <- parts$Participant[j]; x <- m[, p]; k <- k + 1
      ref <- if (p %in% train_hc) mean_sd_from(loo_adjust(st_full, x)) else full
      rt[[k]] <- data.frame(Participant = p, Clip = cl, msz = rt_msz_one(x, ref))
      z <- isc_one(x, ref$mean, parts$max_sec[j])
      isc[[k]] <- data.frame(Participant = rep(p, length(z)), Clip = rep(cl, length(z)), isc_z = z)
    }
  }
  list(rt = bind_rows(rt), isc = bind_rows(isc))
}
holistic_features <- function(train_ids) {
  ref <- PANAS %>% filter(Participant %in% train_ids) %>% group_by(Clip, PANAS_Item) %>%
    summarise(mu = mean(Rating), s = sd(Rating), .groups = "drop")
  PANAS %>% left_join(ref, by = c("Clip", "PANAS_Item")) %>% mutate(z = (Rating - mu) / s, z = ifelse(is.finite(z), z, NA_real_)) %>%
    group_by(Participant, Clip) %>% summarise(msz = mean(z^2, na.rm = TRUE), .groups = "drop")
}
# The real-time features depend only on the training HC, so they are cached by that set
.RT_CACHE <- new.env()
compute_features <- function(train_ids) {
  train_hc <- sort(intersect(train_ids, HC_IDS)); key <- paste(train_hc, collapse = "|")
  if (!exists(key, envir = .RT_CACHE)) assign(key, realtime_features(train_hc), envir = .RT_CACHE)
  rt <- get(key, envir = .RT_CACHE)
  rt_feat <- rt$rt %>% group_by(Participant) %>% summarise(CARMA_MSZ_mean = mean(msz), .groups = "drop")
  isc_feat <- rt$isc %>% group_by(Participant) %>% summarise(ISC_mean = mean(isc_z), ISC_sd = sd(isc_z), .groups = "drop")
  hol_feat <- holistic_features(train_ids) %>% group_by(Participant) %>% summarise(PANAS_MSZ_mean = mean(msz), .groups = "drop")
  IDS %>% left_join(hol_feat, by = "Participant") %>% left_join(rt_feat, by = "Participant") %>%
    left_join(isc_feat, by = "Participant") %>% left_join(SAMPEN_FEAT, by = "Participant") %>%
    left_join(TASK_FEAT, by = "Participant") %>% as.data.frame()
}
validate_features <- function() {
  t0 <- Sys.time(); f <- compute_features(ALL_IDS); el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  p <- as.data.frame(readRDS(output_path("pca_data.rds"))) %>% arrange(Participant); f <- f %>% arrange(Participant)
  stopifnot(identical(p$Participant, f$Participant))
  mx <- sapply(FEATS_ALL, function(v) max(abs(p[[v]] - f[[v]])))
  stopifnot(all(mx < 1e-9))
  rt <- get(paste(HC_IDS, collapse = "|"), envir = .RT_CACHE)
  cat(sprintf("fold-wise features, all participants in training: %.1f s; max |difference| from pca_data.rds %.1e; ",
              el, max(mx)),
      sprintf("ISC windows %d (isc.rds %d)\n", nrow(rt$isc), nrow(readRDS(output_path("isc.rds")))), sep = "")
  invisible(f)
}
