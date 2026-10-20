# Participant characteristics (Table 1, Results 2.1, Methods 4.1): highest education, self-reported sex and gender identity
# of the analytic sample from the questionnaire export, then every Table 1 value and test (Welch t-tests; chi-square tests
# with continuity correction for sex and for race, white vs non-white; Fisher's exact tests for ethnicity, Hispanic/Latino
# vs not, and for education, four levels). Writes demographics_questionnaire.csv (read by later scripts) and table1.csv.
source("R/config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

ids <- read.csv(output_path("analytic_ids.csv"))

# The raw export's header holds the variable names; the label export has the same columns in the same order, so both are
# read under those names. Only the four columns used here are read.
vars <- c(Participant = "participant_id", education = "education2", sex = "what_is_your_sex", gender_identity = "gender")
header <- function(file) names(read_csv(file, n_max = 0, col_types = cols(.default = col_character()), name_repair = "minimal"))
raw_names <- header(input_file("questionnaire_raw"))
stopifnot(all(vars %in% raw_names), length(header(input_file("questionnaire_labels"))) == length(raw_names))
read_answers <- function(file) {
  read_csv(file, col_names = raw_names, skip = 1, col_select = all_of(vars), col_types = cols(.default = col_character()),
           na = character(), progress = FALSE) %>%
    mutate(across(everything(), str_trim))
}
codes <- read_answers(input_file("questionnaire_raw"))
answers <- read_answers(input_file("questionnaire_labels"))
stopifnot(identical(codes$Participant, answers$Participant))
for (v in names(vars)[-1]) {
  pairs <- distinct(data.frame(code = codes[[v]], label = answers[[v]]))
  stopifnot(!anyDuplicated(pairs$code), !anyDuplicated(pairs$label), all((pairs$code == "") == (pairs$label == "")))
}
for (r in RULES$relabel_participants) answers$Participant[answers$Participant %in% r$from] <- r$to
stopifnot(all(ids$Participant %in% answers$Participant))
first_answer <- function(x) c(x[x != ""], "")[1]
questionnaire <- ids %>%
  left_join(answers %>% filter(Participant %in% ids$Participant) %>% group_by(Participant) %>%
              summarise(across(c(education, sex, gender_identity), first_answer)), by = "Participant") %>%
  arrange(Participant)
write_csv(questionnaire, output_path("demographics_questionnaire.csv"), na = "", eol = "\r\n")

# Table 1: age, race, ethnicity and clinical scores from the demographic table; time post onset and the WAB-R AQ from the
# AQ table; education and sex from the questionnaire. Percentages are of each group's sample size.
education_labels <- c("High School/GED" = "High school or GED", "Part College or 2-year degree" = "Some college or 2-year degree",
                      "4-year College degree" = "4-year college degree",
                      "Part or completed Graduate degree" = "Graduate degree (part or completed)")
d <- ids %>%
  left_join(read.csv(input_file("demographics")) %>%
              select(Participant, Age, Race, Ethnicity, MoCA_Score, sct_total_perc, elevator_strings, tmtratio, PHQ_Total, GAD_Total),
            by = "Participant") %>%
  left_join(read.csv(input_file("aphasia_quotient")) %>% select(Participant, MPO, AQ), by = "Participant") %>%
  left_join(questionnaire %>% select(Participant, education, sex), by = "Participant") %>%
  mutate(race = case_when(str_detect(Race, "Black|African") ~ "Black/African American", str_detect(Race, "White") ~ "White",
                          str_detect(Race, "Asian") ~ "Asian", str_detect(Race, "Hispanic|Latino") ~ "Hispanic/Latino",
                          TRUE ~ "Other/Not specified"),
         ethnicity = case_when(str_detect(Ethnicity, "Not Hispanic") ~ "Not Hispanic/Latino",
                               str_detect(Ethnicity, "Hispanic") ~ "Hispanic/Latino", TRUE ~ "Unknown"),
         education = ifelse(education %in% names(education_labels), education_labels[education], "Not reported"),
         sex = na_if(sex, ""),
         years_post_onset = round(MPO / 12, 1),
         navs_percent = 100 * sct_total_perc)
stopifnot(nrow(d) == nrow(ids), !anyNA(d$Age), !anyNA(d$AQ[d$Group == "PWA"]), !anyNA(d$MPO[d$Group == "PWA"]))

# symbols as printed in Table 1, built with intToUtf8() so that they are written as UTF-8 in any locale
dash <- intToUtf8(0x2014)
plus_minus <- intToUtf8(0xB1)
en_dash <- intToUtf8(0x2013)
chi_squared <- intToUtf8(c(0x3C7, 0xB2))
p_text <- function(p) if (p < 0.001) "p < 0.001" else if (round(p, 3) == 1) "p = 1" else sprintf("p = %.3f", p)
table_row <- function(variable, characteristic, PWA = dash, HC = dash, group_difference = dash, ...) {
  data.frame(variable, characteristic, PWA, HC, group_difference, ...)
}
continuous <- function(variable, label, x, groups = c("PWA", "HC"), digits = 1, unit = "") {
  mean_row <- table_row(variable, label)
  range_row <- table_row(variable, "Range")
  for (g in groups) {
    v <- x[d$Group == g & !is.na(x)]
    mean_row[[g]] <- sprintf("%.*f%s %s %.*f%s", digits, mean(v), unit, plus_minus, digits, sd(v), unit)
    mean_row[paste0(g, c("_n", "_mean", "_sd"))] <- list(length(v), mean(v), sd(v))
    range_row[[g]] <- paste0(round(min(v), digits), en_dash, round(max(v), digits), unit)
    range_row[paste0(g, c("_min", "_max"))] <- list(min(v), max(v))
  }
  if (length(groups) == 2) {
    tt <- t.test(x[d$Group == "PWA"], x[d$Group == "HC"])
    mean_row$group_difference <- sprintf("t = %.3f, %s", tt$statistic, p_text(tt$p.value))
    mean_row[c("test", "statistic", "df", "p")] <- list("Welch t-test", unname(tt$statistic), unname(tt$parameter), tt$p.value)
  }
  bind_rows(mean_row, range_row)
}
categorical <- function(variable, label, x, levels, test, test_name, tested) {
  chisq <- !is.null(test$statistic)
  difference <- if (chisq) sprintf("%s = %.3f, %s", chi_squared, test$statistic, p_text(test$p.value)) else
    paste("Fisher's exact", p_text(test$p.value))
  test_row <- table_row(variable, label, PWA = "", HC = "", group_difference = difference, test = test_name,
                        statistic = if (chisq) unname(test$statistic) else NA, df = if (chisq) unname(test$parameter) else NA,
                        p = test$p.value, PWA_n = sum(tested & d$Group == "PWA"), HC_n = sum(tested & d$Group == "HC"))
  level_rows <- lapply(levels, function(level) {
    out <- table_row(variable, level)
    for (g in c("PWA", "HC")) {
      n <- sum(d$Group == g & x %in% level)
      total <- sum(d$Group == g)
      out[[g]] <- paste0(n, " (", round(100 * n / total, 1), "%)")
      out[paste0(g, c("_n", "_count", "_percent"))] <- list(total, n, 100 * n / total)
    }
    out
  })
  bind_rows(test_row, level_rows)
}

known_ethnicity <- d$ethnicity != "Unknown"
known_education <- d$education != "Not reported"
group_n <- table(d$Group)
table1 <- bind_rows(
  table_row("n", "Sample size, n", PWA = as.character(group_n[["PWA"]]), HC = as.character(group_n[["HC"]]),
            PWA_n = group_n[["PWA"]], HC_n = group_n[["HC"]]),
  table_row("section", "Demographic Characteristics", PWA = "", HC = "", group_difference = ""),
  continuous("age", "Age (years), mean + SD", d$Age),
  categorical("sex", "Sex, n (%)", d$sex, sort(unique(na.omit(d$sex))), chisq.test(table(d$Group, d$sex)),
              "Chi-square test with continuity correction", !is.na(d$sex)),
  categorical("race", "Race, n (%)", d$race, sort(unique(d$race)), chisq.test(table(d$Group, d$race == "White")),
              "Chi-square test with continuity correction, white vs non-white", rep(TRUE, nrow(d))),
  categorical("ethnicity", "Ethnicity, n (%)", d$ethnicity, sort(unique(d$ethnicity)),
              fisher.test(table(d$Group[known_ethnicity], d$ethnicity[known_ethnicity])),
              "Fisher's exact test, Hispanic/Latino vs not Hispanic/Latino", known_ethnicity),
  categorical("education", "Education, n (%)", d$education, c(unname(education_labels), "Not reported"),
              fisher.test(table(d$Group[known_education], d$education[known_education])),
              "Fisher's exact test, four levels", known_education),
  table_row("section", "Clinical Characteristics", PWA = "", HC = "", group_difference = ""),
  continuous("years_post_onset", "Time post aphasia onset (years), mean + SD", d$years_post_onset, "PWA"),
  continuous("aq", "WAB-R Aphasia Quotient, mean + SD", d$AQ, "PWA"),
  continuous("moca", "MoCA Score, mean + SD", d$MoCA_Score, "HC"),
  continuous("navs", "NAVS Sentence Comprehension (%), mean + SD", d$navs_percent, "PWA", unit = "%"),
  continuous("tea_elevator", "TEA Elevator Strings (scaled), mean + SD", d$elevator_strings, "PWA"),
  continuous("tmt_ratio", "TMT-B/A Ratio, mean + SD", d$tmtratio, "PWA", digits = 2),
  continuous("phq8", "PHQ-8 Total, mean + SD", d$PHQ_Total),
  continuous("gad7", "GAD-7 Total, mean + SD", d$GAD_Total))
numbers <- c("n", "count", "percent", "mean", "sd", "min", "max")
table1 <- table1 %>% select(variable, characteristic, PWA, HC, group_difference, test, statistic, df, p,
                            all_of(paste0("PWA_", numbers)), all_of(paste0("HC_", numbers)))
write_csv(table1, output_path("table1.csv"), na = "")

cat(sprintf("questionnaire: %d participants; education reported by %d, sex by %d, gender identity by %d\n", nrow(questionnaire),
            sum(known_education), sum(!is.na(d$sex)), sum(questionnaire$gender_identity != "")))
cat(sprintf("Table 1: %d rows; PWA n = %d, HC n = %d\n", nrow(table1), group_n[["PWA"]], group_n[["HC"]]))
tests <- table1 %>% filter(!is.na(p))
cat(sprintf("  %s %s (n = %d PWA, %d HC)\n", tests$characteristic, tests$group_difference, tests$PWA_n, tests$HC_n), sep = "")
