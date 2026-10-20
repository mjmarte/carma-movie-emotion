# Response-method check within PWA (Supplementary Results S1, Methods 4.3): joystick (in person) versus keyboard (remote)
# as a fixed effect on each real-time feature, with aphasia severity (AQ) and age, and random intercepts for participant and
# clip: sample entropy at each scale and ISC at the window level, real-time MSZ per participant-clip. Writes input_method.csv.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

ids <- read.csv(output_path("analytic_ids.csv"))
method <- read.csv(input_file("aphasia_quotient")) %>%
  filter(Participant %in% ids$Participant[ids$Group == "PWA"]) %>%
  transmute(Participant, AQ, Method = factor(ifelse(Joystick == 1, "Joystick", "Keyboard"), levels = c("Joystick", "Keyboard"))) %>%
  left_join(read.csv(input_file("demographics")) %>% select(Participant, Age), by = "Participant")
stopifnot(nrow(method) == 57, !anyNA(method))
fit <- function(data, dv, label) {
  model <- lmer(as.formula(paste(dv, "~ Method + AQ + Age + (1 | Participant) + (1 | Clip)")), data = data %>% inner_join(method, by = "Participant"))
  co <- coef(summary(model))["MethodKeyboard", ]
  data.frame(measure = label, estimate = co[["Estimate"]], se = co[["Std. Error"]], df = co[["df"]], t = co[["t value"]],
             p = co[["Pr(>|t|)"]], n_obs = nobs(model))
}
entropy <- readRDS(output_path("sample_entropy.rds"))
results <- bind_rows(
  lapply(c("Short", "Medium", "Long"), function(b) fit(entropy %>% filter(WindowBin == b), "SampEn", paste("Sample entropy,", tolower(b)))),
  fit(readRDS(output_path("isc.rds")), "isc_z", "ISC"),
  fit(readRDS(output_path("realtime_msz.rds")), "msz", "Real-time MSZ"))
write.csv(results, output_path("input_method.csv"), row.names = FALSE)
cat(sprintf("Joystick %d, keyboard %d PWA\n", sum(method$Method == "Joystick"), sum(method$Method == "Keyboard")))
with(results, cat(sprintf("%-25s b = %.3f, SE = %.3f, t = %.2f, p = %.3f\n", measure, estimate, se, t, p), sep = ""))
