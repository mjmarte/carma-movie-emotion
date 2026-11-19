# Figures 1 and 2 (Results 2.3): per clip and group, the participant-clip observations (jittered points), the group mean with
# its 95% bootstrap CI (diamond and whisker), and the HC and PWA kernel densities back to back on the clip's centre line.
# Figure 1: antonym generation and comprehension accuracy; Figure 2: mean PANAS-SF rating of the five positive and the five
# negative items. Writes figure1.png, figure2.png and the plotted data, figure1_data.csv and figure2_data.csv.
source("R/config.R")
source("R/theme_marte.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
marte_defaults()

ids <- read.csv(output_path("analytic_ids.csv"))
accuracy <- readRDS(output_path("task_responses.rds")) %>%
  group_by(Participant, Group, Clip, Task) %>%
  summarise(accuracy = mean(Accuracy, na.rm = TRUE), n_items = sum(!is.na(Accuracy)), .groups = "drop") %>%
  filter(!is.nan(accuracy))
panas <- read.csv(input_file("panas_items")) %>%
  select(-Group) %>%
  inner_join(ids, by = c("participant_id" = "Participant")) %>%
  group_by(Participant = participant_id, Group, Clip = clip, valence) %>%
  summarise(mean_rating = mean(rating), n_items = n(), .groups = "drop")
write.csv(accuracy, output_path("figure1_data.csv"), row.names = FALSE)
write.csv(panas, output_path("figure2_data.csv"), row.names = FALSE)

group_colors <- c(HC = "#0072B2", PWA = "#CC0000")
abbrev <- c(Akeelah = "AKB", CatchMe = "CMIYC", CatchMeIfYouCan = "CMIYC", GoodWillHunting = "GWH", Miracle = "MIR",
            Moonlight = "MOON", `No Country for Old Men` = "NCOM", NoCountryForOldMen = "NCOM", ParentTrap = "PT",
            PartlyCloudy = "PC")
clip_order <- c("AKB", "CMIYC", "NCOM", "PT", "PC", "GWH", "MIR", "MOON")
fig1 <- accuracy %>% mutate(clip = factor(abbrev[Clip], levels = clip_order))
fig2 <- panas %>% mutate(clip = factor(abbrev[Clip], levels = clip_order))
stopifnot(!anyNA(fig1$clip), !anyNA(fig2$clip), nrow(fig1) == 1412,
          n_distinct(fig1$Participant) == 100, n_distinct(fig1$Participant[fig1$Group == "PWA"]) == 57,
          n_distinct(fig2$Participant) == 100, n_distinct(fig2$Participant[fig2$Group == "PWA"]) == 57)

boot_ci <- function(v, B = 2000) {
  set.seed(42)
  m <- replicate(B, mean(sample(v, length(v), replace = TRUE)))
  c(lo = unname(quantile(m, .025)), hi = unname(quantile(m, .975)))
}
# Horizontal positions relative to the clip (one clip = 1 x unit): cloud baseline and direction (HC left, PWA right), mean
# diamond, centre and half-width of the band of points
cloud_width <- 0.19
positions <- data.frame(Group = factor(c("HC", "PWA"), levels = c("HC", "PWA")), cloud0 = c(0, 0), dir = c(-1, 1),
                        diam = c(-0.26, 0.26), rain = c(-0.375, 0.375), rain_hw = 0.045)
# Density (nrd0 bandwidth) trimmed to the observed range, scaled to the same maximum width in every cloud
cloud_poly <- function(v, x0, dir) {
  d <- density(v, n = 512, from = min(v), to = max(v))
  w <- cloud_width * d$y / max(d$y)
  data.frame(x = c(x0, x0 + dir * w, x0), y = c(d$x[1], d$x, d$x[512]))
}
rain_panel <- function(df, yvar, ylab, ylim, ybreaks, ylabels = waiver(), title, seed) {
  df <- df %>% mutate(Group = factor(Group, levels = c("HC", "PWA")), y = .data[[yvar]], ci = as.integer(clip)) %>%
    left_join(positions, by = "Group")
  set.seed(seed)
  df$xr <- df$ci + df$rain + runif(nrow(df), -1, 1) * df$rain_hw
  st <- df %>% group_by(clip, Group, ci, diam) %>%
    summarise(m = mean(y), lo = boot_ci(y)["lo"], hi = boot_ci(y)["hi"], .groups = "drop")
  cells <- df %>% group_by(clip, Group, ci, cloud0, dir) %>% summarise(k = n_distinct(y), v1 = first(y), .groups = "drop")
  polys <- df %>% semi_join(cells %>% filter(k > 1), by = c("clip", "Group")) %>% group_by(clip, Group) %>%
    reframe(cloud_poly(y, first(ci) + first(cloud0), first(dir))) %>% mutate(id = paste(clip, Group))
  # a group whose values are all identical gets a flat cloud: a thin bar at that value, 1.2% of the axis span high
  h <- 0.012 * diff(range(ybreaks))
  flat <- cells %>% filter(k == 1) %>%
    mutate(ymin = ifelse(v1 >= max(ybreaks), v1 - h, ifelse(v1 <= min(ybreaks), v1, v1 - h / 2)), ymax = ymin + h,
           x0 = ci + cloud0, x1 = ci + cloud0 + dir * cloud_width)
  ggplot() +
    geom_polygon(data = polys, aes(x, y, group = id, fill = Group), colour = CAVIAR_EDGE, linewidth = 0.3, alpha = 0.55,
                 show.legend = FALSE) +
    geom_rect(data = flat, aes(xmin = pmin(x0, x1), xmax = pmax(x0, x1), ymin = ymin, ymax = ymax, fill = Group),
              colour = CAVIAR_EDGE, linewidth = 0.3, alpha = 0.55, show.legend = FALSE) +
    geom_point(data = df, aes(xr, y, fill = Group), shape = 21, colour = CAVIAR_EDGE, stroke = 0.15, size = 1.0, alpha = 0.45,
               show.legend = FALSE) +
    geom_linerange(data = st, aes(x = ci + diam, ymin = lo, ymax = hi), colour = "#222222", linewidth = 0.8) +
    geom_point(data = st, aes(x = ci + diam, y = m, fill = Group), shape = 23, size = 3.4, colour = "#222222", stroke = 0.6) +
    scale_fill_manual(values = group_colors, name = NULL) +
    scale_x_continuous(breaks = seq_along(clip_order), labels = clip_order, limits = c(0.5, length(clip_order) + 0.5),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, breaks = ybreaks, labels = ylabels, expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = NULL, y = ylab, title = title) +
    theme(legend.position = "bottom", plot.title = element_text(size = 12),
          axis.text.x = element_text(size = 10, colour = "#222222"), panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
}

panels <- list(`Figure 1A` = fig1 %>% filter(Task == "AntGen"), `Figure 1B` = fig1 %>% filter(Task == "CompQ"),
               `Figure 2A` = fig2 %>% filter(valence == "positive"), `Figure 2B` = fig2 %>% filter(valence == "negative"))
panel_a <- rain_panel(panels$`Figure 1A`, "accuracy", "Accuracy", c(-0.02, 1.02), seq(0, 1, .25), scales::percent,
                      "A. Antonym generation", 11) + theme(legend.position = "none")
panel_b <- rain_panel(panels$`Figure 1B`, "accuracy", "Accuracy", c(-0.02, 1.02), seq(0, 1, .25), scales::percent,
                      "B. Comprehension questions", 12)
save_marte(panel_a / panel_b, output_path("figure1.png"), width = 10, height = 6.4)
ggsave(output_path("figure1.pdf"), panel_a / panel_b, width = 10, height = 6.4, device = grDevices::cairo_pdf)
panel_a <- rain_panel(panels$`Figure 2A`, "mean_rating", "Mean rating", c(0.95, 5.05), 1:5, waiver(), "A. Positive affect", 21) +
  theme(legend.position = "none")
panel_b <- rain_panel(panels$`Figure 2B`, "mean_rating", "Mean rating", c(0.95, 5.05), 1:5, waiver(), "B. Negative affect", 22)
save_marte(panel_a / panel_b, output_path("figure2.png"), width = 10, height = 6.4)
ggsave(output_path("figure2.pdf"), panel_a / panel_b, width = 10, height = 6.4, device = grDevices::cairo_pdf)

for (p in names(panels)) {
  d <- panels[[p]]
  cat(sprintf("%s: %d HC and %d PWA participants; %d HC and %d PWA participant-clip observations\n", p,
              n_distinct(d$Participant[d$Group == "HC"]), n_distinct(d$Participant[d$Group == "PWA"]),
              sum(d$Group == "HC"), sum(d$Group == "PWA")))
}
