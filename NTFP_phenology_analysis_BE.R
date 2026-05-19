####################################################################
###   Brazil nut phenology analysis      

#Set working directory
setwd()

# Libraries 
library(dplyr)
library(tidyverse)
library(readxl)
library(lubridate)
library(knitr)
library(gt)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(zoo)
library(ordinal)
library(circular)
library(lme4)
library(lmerTest)
library(patchwork)
library(MuMIn)


###############################
# Read and clean phenology data

BE_pheno <- read.csv("BE_phenology_data.csv")

# Rename columns and remove empty columns/rows
BE_pheno <- BE_pheno %>%
  rename(
    tree.no      = 1,
    flwrs.ground = 6,
    nuts.tree    = 7,
    nuts.ground  = 8,
    canopy.leaves = 9,
    green.leaves  = 10
  )
BE_pheno[, c(11, 12)] <- NULL
BE_pheno[BE_pheno == ""] <- NA
BE_pheno <- BE_pheno[rowSums(is.na(BE_pheno)) != ncol(BE_pheno), ]
rownames(BE_pheno) <- NULL

# Parse dates
BE_pheno$Date2 <- as.Date(parse_date_time(BE_pheno$Date, orders = c("dmy", "mdy")))

# Manually correct two batches of mis-entered dates
BE_pheno[c(4669:4692), ]$Date2 <- as.Date("2025-04-08")
BE_pheno[c(5125:5148), ]$Date2 <- as.Date("2026-01-08")

# Remove rows where date is still NA
BE_pheno <- BE_pheno[!is.na(BE_pheno$Date2), ]

# Standardize tree IDs
BE_pheno$tree.no <- ifelse(BE_pheno$tree.no %in% c("Z", "z"), 44, BE_pheno$tree.no)
BE_pheno <- BE_pheno %>%
  filter(!tree.no %in% c("Y", "X"))
BE_pheno$tree.no <- as.numeric(BE_pheno$tree.no)

# Add date-derived columns
BE_pheno <- BE_pheno %>%
  mutate(
    year   = year(Date2),
    month  = month.abb[month(Date2)],
    julian = as.numeric(format(Date2, "%j")),
    week   = isoweek(Date2)
  )

# Remove observations for trees after mortality
BE_pheno <- BE_pheno %>%
  filter(!(tree.no == 39 & Date2 >= as.Date("2021-12-13"))) # tree fell prior to this census
BE_pheno <- BE_pheno %>%
  filter(!(tree.no == 34 & Date2 >= as.Date("2019-11-14"))) # tree reported to show signs of lightning strike on this census
BE_pheno <- BE_pheno %>%
  filter(!(tree.no == 29 & Date2 >= as.Date("2023-11-06"))) # tree lost the majority of it crown prior to this census, shortly after lost the remainder of its crown
BE_pheno <- BE_pheno %>%
  filter(!(tree.no == 23 & Date2 >= as.Date("2022-06-07"))) # tree declined to below 50% leaves in the canopy, signaling mortality. within two months there were no leaves and tree was leafless for the remainder of monitoring
BE_pheno <- BE_pheno %>%
  filter(Date2 != as.Date("2018-07-26")) # no flower or fruit observations were made this census

# Read in tree metadata
BE_trees <- read.csv("BE_locations.csv")
BE_trees <- BE_trees %>%
  rename(tree.no = name)
BE_trees$habitat <- ifelse(BE_trees$habitat == "Forest", "Primary", "Secondary")

# Temporal summary per tree (for QC / visualization)
temporal_summary <- BE_pheno %>%
  group_by(tree.no) %>%
  summarise(
    first_obs      = min(Date2),
    last_obs       = max(Date2),
    n_observations = n(),
    .groups = "drop"
  ) %>%
  arrange(first_obs) %>%
  left_join(BE_trees, by = "tree.no")

BE_time_series <- ggplot(temporal_summary, aes(y = tree.no, shape = habitat, fill = habitat)) +
  geom_segment(aes(x = first_obs, xend = last_obs, yend = tree.no), linewidth = 1) +
  geom_point(aes(x = first_obs), color = "#238443", size = 3) +
  geom_point(aes(x = last_obs),  color = "black",    size = 3) +
  scale_shape_manual(values = c(19,24)) +
  scale_fill_manual(values = c("#238443", "white")) +
  labs(title = "Bertholletia excelsa", x = "", y = "") +
  theme_minimal() +
  guides(shape = guide_legend(title = "Forest type"),
         fill = guide_legend(title = "Forest type")) +
  theme(axis.text.y = element_text(size = 8))


##############################
# Read and clean climate data

weather <- read_excel(
  "ASA_FLP_weather_climate_data_MASTER.xlsx",
  col_types = c("date", "numeric", "numeric", "numeric",
                "text", "numeric", "numeric", "text", "text")
) %>%
  rename(
    Date        = 1,
    MinTemp     = 2,
    MaxTemp     = 3,
    CurrentTemp = 4,
    Time        = 5,
    Precip.mm   = 6
  ) %>%
  select(-c(7, 8)) %>%
  filter(Date <= "2026-04-01") %>%
  mutate(
    Date      = as.Date(Date),
    Precip.mm = ifelse(is.na(Precip.mm), 0, Precip.mm),
    Precip.mm.avg = rollmean(Precip.mm, k = 30, fill = NA),
    year   = year(Date),
    month  = month.abb[month(Date)],
    julian = as.numeric(format(Date, "%j"))
  )

# Weekly climate summaries (for joining to phenology data)
weather_weekly <- weather %>%
  mutate(week = isoweek(Date)) %>%
  group_by(year, week) %>%
  summarize(
    Precip.mm  = sum(Precip.mm,     na.rm = TRUE),
    MinTemp    = mean(MinTemp,      na.rm = TRUE),
    MaxTemp    = mean(MaxTemp,      na.rm = TRUE),
    Precip.avg = mean(Precip.mm.avg, na.rm = TRUE),
    .groups = "drop"
  )

# Monthly climate averages (for circular phenology plots)
weather_monthly <- weather %>%
  mutate(week = isoweek(Date)) %>%
  group_by(year, month) %>%
  summarize(
    precip.mm = sum(Precip.mm.avg, na.rm = TRUE),
    MinTemp   = mean(MinTemp,      na.rm = TRUE),
    MaxTemp   = mean(MaxTemp,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(month) %>%
  summarize(
    precip.mm = mean(precip.mm, na.rm = TRUE),
    MinTemp   = mean(MinTemp,   na.rm = TRUE),
    MaxTemp   = mean(MaxTemp,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(month.no = match(month, month.abb)) %>%
  arrange(month.no)


################################################
# Merge phenology and climate, clean some scores

# Join weekly climate to phenology observations
BE_pheno_weather <- full_join(BE_pheno, weather_weekly, by = c("year", "week")) %>%
  filter(!is.na(Date2))

# Canopy fruit score: numeric
BE_pheno_weather$nuts.tree <- as.numeric(BE_pheno_weather$nuts.tree)

# Fallen fruit score: recode messy text entries to 0-3 ordinal scale
recode_nuts_ground <- function(x) {
  x_clean <- str_trim(str_to_lower(x))
  case_when(
    x_clean %in% c("na", "none", "no", "0", ".", "-", "noen", "none ")   ~ 0,
    x_clean == ""                                                          ~ 0,
    str_detect(x_clean, "^few")                                           ~ 1,
    x_clean %in% c("1", "2", "3")                                         ~ 1,
    str_detect(x_clean, "^som")                                           ~ 2,
    x_clean == "one"                                                       ~ 2,
    x_clean %in% c("7", "14", "18")                                       ~ 2,
    x_clean == "several" | str_detect(x_clean, "^several")               ~ 2,
    str_detect(x_clean, "^man")                                           ~ 3,
    x_clean %in% c("32", "59")                                            ~ 3,
    TRUE                                                                   ~ NA_real_
  )
}

BE_pheno_weather <- BE_pheno_weather %>%
  mutate(
    nuts.ground    = recode_nuts_ground(nuts.ground),
    biweek_period  = ceiling(week / 2)
  )

# Check recoding results
table(BE_pheno_weather$nuts.ground, useNA = "always")

# Clean flower scores
BE_pheno_weather$flwrs.ground <- ifelse(BE_pheno_weather$month == "Aug", 0, BE_pheno_weather$flwrs.ground) #August flower scores >0 were actually for flowers from a nearby Eschweilera tree.
bad_flwr_scores <- c("-", ".", "- ")
BE_pheno_weather$flwrs.ground <- ifelse(BE_pheno_weather$flwrs.ground %in% bad_flwr_scores, 0, BE_pheno_weather$flwrs.ground)

# Convert julian day to circular radians for circular statistics
BE_pheno_weather$radians <- circular(
  BE_pheno_weather$julian * (2 * pi / 365),
  type = "angles", units = "radians",
  template = "geographics", modulo = "2pi"
)


##########################
# Circular phenology plots

# Summarize scores by biweekly period
flower_counts <- BE_pheno_weather %>%
  mutate(flwrs.ground = as.numeric(flwrs.ground)) %>%
  group_by(biweek_period, flwrs.ground) %>%
  summarise(n_observations = n(), .groups = "keep") %>%
  group_by(biweek_period) %>%
  mutate(
    total_period_obs = sum(n_observations),
    prop_score       = n_observations / total_period_obs
  ) %>%
  ungroup() %>%
  filter(!is.na(flwrs.ground))

nuts_tree_counts <- BE_pheno_weather %>%
  filter(!is.na(nuts.tree)) %>%
  group_by(biweek_period, nuts.tree) %>%
  summarise(n_observations = n(), .groups = "keep") %>%
  group_by(biweek_period) %>%
  mutate(
    total_period_obs = sum(n_observations),
    prop_score       = n_observations / total_period_obs
  ) %>%
  ungroup() %>%
  filter(!is.na(nuts.tree))

nuts_ground_counts <- BE_pheno_weather %>%
  group_by(biweek_period, nuts.ground) %>%
  summarise(n_observations = n(), .groups = "keep") %>%
  group_by(biweek_period) %>%
  mutate(
    total_period_obs = sum(n_observations),
    prop_score       = n_observations / total_period_obs
  ) %>%
  ungroup() %>%
  filter(!is.na(nuts.ground))

# Biweekly precipitation for polar plot overlay
precip_biweekly <- BE_pheno_weather %>%
  group_by(biweek_period) %>%
  summarise(mean_precip = mean(Precip.mm, na.rm = TRUE), .groups = "drop")

# Add plotting coordinates
add_plot_coords <- function(df) {
  df %>%
    mutate(week_equiv = (biweek_period - 1) * 2 + 1.5) %>%
    left_join(precip_biweekly, by = "biweek_period")
}

flower_counts      <- add_plot_coords(flower_counts)
nuts_tree_counts   <- add_plot_coords(nuts_tree_counts)
nuts_ground_counts <- add_plot_coords(nuts_ground_counts)

# Factor levels for scores
flower_counts$flwrs.ground   <- factor(flower_counts$flwrs.ground,   levels = c("0","1","2","3","4"))
nuts_tree_counts$nuts.tree   <- factor(nuts_tree_counts$nuts.tree,    levels = c("0","1","2","3","4"))
nuts_ground_counts$nuts.ground <- factor(nuts_ground_counts$nuts.ground, levels = c("0","1","2","3"))

# Color palettes
score_colors <- c("0" = "transparent", "1" = "#ffffcc", "2" = "#c2e699",
                  "3" = "#78c679", "4" = "#238443")
fruit_colors <- c("0" = "transparent", "1" = "#dde5b6", "2" = "#adc178", "3" = "#7f4f24")

# x-axis month breaks and labels
month_breaks <- seq(0.5, 51.5, 4.3)
month_labels <- c("Jan","Feb","Mar","Apr","May","Jun",
                  "Jul","Aug","Sep","Oct","Nov","Dec")

weather_monthly <- weather_monthly %>%
  mutate(
    week_equiv = c(month_breaks[1:11] + 2.15, 47.8 + 2.35),
    bar_width  = c(rep(4.3, 11), 4.75)
  )

# Shared polar plot theme
polar_theme <- theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 10),
    plot.title   = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "gray70", linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom"
  )

# Helper to build each polar plot
make_polar_plot <- function(data, score_col, title, colors) {
  ggplot(data[data[[score_col]] != 0, ],
         aes(x = week_equiv, y = prop_score, fill = as.factor(.data[[score_col]]))) +
    geom_col(data = weather_monthly,
             aes(y = precip.mm * 0.0025, fill = NULL, width = bar_width),
             fill = "lightblue", alpha = 1) +
    geom_hline(yintercept = c(0, 0.25, 0.5, 0.75, 1), color = "gray70", linewidth = 0.3) +
    geom_vline(xintercept = month_breaks,               color = "gray70", linewidth = 0.3) +
    geom_col(position = position_stack(reverse = TRUE), color = "black", width = 2, alpha = 1) +
    coord_polar(theta = "x", start = 0) +
    scale_fill_manual(values = colors, name = "Score") +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1),
      sec.axis = sec_axis(~ . / 0.0025, name = "")
    ) +
    scale_x_continuous(breaks = month_breaks, labels = month_labels) +
    labs(title = title, x = "", y = "") +
    polar_theme
}

BE_fl  <- make_polar_plot(flower_counts,      "flwrs.ground", "Fallen flowers", score_colors)
BE_fr  <- make_polar_plot(nuts_tree_counts,   "nuts.tree",    "Canopy fruit",   score_colors)
BE_fr2 <- make_polar_plot(nuts_ground_counts, "nuts.ground",  "Fallen fruit",   fruit_colors)

fig1a <- ggarrange(BE_fl, BE_fr, ncol = 1, nrow = 2,
                   common.legend = TRUE, labels = c("A", "B"))
fig1a

######################################
# Circular statistics

# Helper: run circular stats for one phenophase, return a results row
circular_stats <- function(data, radians_col, label) {
  
  dens       <- density.circular(data[[radians_col]], bw = 20)
  mode_idx   <- which.max(dens$y)
  mode_doy   <- (as.numeric(dens$x[mode_idx]) %% (2 * pi)) * (365 / (2 * pi))
  mode_date  <- as.Date(mode_doy, origin = "2023-12-31")
  circ_var   <- var.circular(data[[radians_col]])
  rayleigh   <- rayleigh.test(data[[radians_col]])
  
  data.frame(
    species          = "Brazil nut",
    phenophase       = label,
    n                = nrow(data),
    test_statistic   = round(as.numeric(rayleigh$statistic), 3),
    p_value          = round(rayleigh$p.value, 3),
    peak_doy         = round(mode_doy, 1),
    peak_date        = format(mode_date, "%B %d"),
    circular_variance = round(circ_var, 3)
  )
}

# Helper: compute CVp, xCVi, xPCC for one phenophase
masting_stats <- function(annual_df, score_col) {
  
  pop_means <- annual_df %>%
    group_by(year) %>%
    summarize(pop_mean = mean(.data[[score_col]], na.rm = TRUE), .groups = "drop")
  CVp <- sd(pop_means$pop_mean) / mean(pop_means$pop_mean)
  
  ind_CVs <- annual_df %>%
    group_by(tree.no) %>%
    summarize(
      m  = mean(.data[[score_col]], na.rm = TRUE),
      s  = sd(.data[[score_col]],   na.rm = TRUE),
      CV = s / m * 100,
      .groups = "drop"
    ) %>%
    filter(!is.na(CV) & is.finite(CV))
  xCVi <- mean(ind_CVs$CV, na.rm = TRUE) / 100
  
  wide <- annual_df %>%
    select(year, tree.no, all_of(score_col)) %>%
    pivot_wider(names_from = tree.no, values_from = all_of(score_col))
  cm   <- cor(wide %>% select(-year), use = "pairwise.complete.obs")
  xPCC <- mean(cm[upper.tri(cm)], na.rm = TRUE)
  
  c(CVp = CVp, xCVi = xCVi, xPCC = xPCC)
}

# Flowers
BE_fl_data <- BE_pheno_weather %>% filter(flwrs.ground > 0)
fl_circular <- circular_stats(BE_fl_data, "radians", "Flowers")

fl_annual <- BE_pheno_weather %>%
  mutate(score = as.numeric(flwrs.ground)) %>%
  group_by(year, tree.no) %>%
  summarize(mean_flower_score = mean(score, na.rm = TRUE), .groups = "drop")
fl_masting <- masting_stats(fl_annual, "mean_flower_score")

# Canopy fruit
BE_ft_data <- BE_pheno_weather %>% filter(nuts.tree > 0)
ft_circular <- circular_stats(BE_ft_data, "radians", "Canopy fruit")

ft_annual <- BE_pheno_weather %>%
  mutate(score = as.numeric(nuts.tree)) %>%
  group_by(year, tree.no) %>%
  summarize(mean_fruit_score = mean(score, na.rm = TRUE), .groups = "drop")
ft_masting <- masting_stats(ft_annual, "mean_fruit_score")

# Fallen fruit
BE_fg_data <- BE_pheno_weather %>% filter(nuts.ground > 0)
fg_circular <- circular_stats(BE_fg_data, "radians", "Fallen fruit")

fg_annual <- BE_pheno_weather %>%
  mutate(score = as.numeric(nuts.ground)) %>%
  group_by(year, tree.no) %>%
  summarize(mean_fruit_score = mean(score, na.rm = TRUE), .groups = "drop")
fg_masting <- masting_stats(fg_annual, "mean_fruit_score")

# Combine circular results table
circular_results <- bind_rows(fl_circular, ft_circular, fg_circular)
# write.csv(circular_results, "BE_circle.csv")

# Combine masting metrics table
masting_results_BE <- data.frame(
  Species   = "Brazil nut",
  Phenophase = c("Flowers", "Canopy fruit", "Fallen fruit"),
  CVp  = c(fl_masting["CVp"],  ft_masting["CVp"],  fg_masting["CVp"]),
  xCVi = c(fl_masting["xCVi"], ft_masting["xCVi"], fg_masting["xCVi"]),
  xPCC = c(fl_masting["xPCC"], ft_masting["xPCC"], fg_masting["xPCC"])
)


#######################################
# Seasonal climate variables with lags

dry_months <- c("May", "Jun", "Jul", "Aug", "Sep")

# Seasonal climate totals/means per water year (Oct–Sep)
season_precip <- weather %>%
  mutate(
    water_year = ifelse(month(Date) >= 10, year + 1, year),
    season     = ifelse(month %in% dry_months, "dry", "wet")
  ) %>%
  filter(Date >= as.Date("2017-10-01")) %>%
  group_by(water_year, season) %>%
  summarize(
    total_precip_mm = sum(Precip.mm,  na.rm = TRUE),
    mean_min_temp   = mean(MinTemp,   na.rm = TRUE),
    mean_max_temp   = mean(MaxTemp,   na.rm = TRUE),
    mean_precip     = mean(Precip.mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_temp = (mean_min_temp + mean_max_temp) / 2,
    season    = factor(season, levels = c("wet", "dry"))
  ) %>%
  arrange(water_year, season)

# Add 1- and 2-season lags to climate predictors
season_precip_lagged <- season_precip %>%
  arrange(water_year, season) %>%
  mutate(
    precip_prev_season      = lag(total_precip_mm, 1),
    temp_min_prev_season    = lag(mean_min_temp,   1),
    temp_max_prev_season    = lag(mean_max_temp,   1),
    precip_2prev_season     = lag(total_precip_mm, 2),
    temp_min_2prev_season   = lag(mean_min_temp,   2),
    temp_max_2prev_season   = lag(mean_max_temp,   2)
  ) %>%
  ungroup()


#########################################################
# Lagged canopy fruit scores (for fallen fruit analysis)

# Seasonal mean canopy fruit score per individual tree
fruit_tree_seasonal <- BE_pheno_weather %>%
  mutate(
    fruit_score = as.numeric(nuts.tree),
    Date2       = as.Date(Date2),
    water_year  = ifelse(month(Date2) >= 10, year + 1, year),
    season      = ifelse(month %in% dry_months, "dry", "wet"),
    season      = factor(season, levels = c("wet", "dry"))
  ) %>%
  filter(Date2 >= as.Date("2017-10-01")) %>%
  group_by(water_year, season, tree.no) %>%
  summarize(
    mean_canopy_fruit = mean(fruit_score, na.rm = TRUE),
    n_obs             = n(),
    .groups           = "drop"
  ) %>%
  arrange(tree.no, water_year, season)

# Lag canopy fruit score within each tree (1 and 2 seasons back)
fruit_tree_seasonal <- fruit_tree_seasonal %>%
  group_by(tree.no) %>%
  mutate(
    canopy_fruit_prev_season  = lag(mean_canopy_fruit, 1),
    canopy_fruit_2prev_season = lag(mean_canopy_fruit, 2)
  ) %>%
  ungroup()


######################################################
# Calculate phenological year for canopy fruit scores

# Identifies the calendar start day that minimizes within-year
# variance of the score distribution — i.e., the optimal
# phenological year boundary for canopy fruit
calculate_phenological_year <- function(dates, scores) {
  doy        <- as.numeric(format(dates, "%j")) - 1
  valid_idx  <- scores > 0 & !is.na(scores) & !is.na(doy)
  doy        <- doy[valid_idx]
  scores     <- scores[valid_idx]
  agg_data   <- aggregate(scores, by = list(doy = doy), FUN = sum)
  names(agg_data) <- c("day", "quantity")
  day <- agg_data$day
  q   <- agg_data$quantity
  n   <- sum(q)
  
  variances <- numeric(366)
  for (i in 0:365) {
    x          <- ifelse(day > i, day - i, 365 + day - i)
    variances[i + 1] <- sum(x^2 * q) - sum(x * q)^2 / n
  }
  
  min_var  <- min(variances, na.rm = TRUE)
  min_days <- which(variances == min_var) - 1
  
  if (length(min_days) > 1) {
    rads      <- circular(2 * pi * min_days / 366, units = "radians", modulo = "2pi")
    med_rad   <- as.numeric(median.circular(rads))
    rad_vals  <- as.numeric(rads)
    closest   <- which(rad_vals >= med_rad)[1]
    if (is.na(closest)) closest <- 1
    optimal_shift <- min_days[closest]
  } else {
    optimal_shift <- min_days[1]
  }
  
  start_doy  <- optimal_shift + 1
  if (start_doy > 365) start_doy <- 1
  start_date <- as.Date(start_doy - 1, origin = "2020-01-01")
  
  list(
    optimal_shift = optimal_shift,
    start_doy     = start_doy,
    start_date    = format(start_date, "%B %d"),
    min_variance  = min_var,
    all_variances = variances,
    n_min_days    = length(min_days)
  )
}

fruit_pheno_year <- calculate_phenological_year(
  dates  = BE_pheno_weather$Date2,
  scores = as.numeric(BE_pheno_weather$nuts.tree)
)

# Assign a phenological year label to each observation date
assign_pheno_year <- function(date, start_doy) {
  doy           <- as.numeric(format(date, "%j"))
  calendar_year <- as.numeric(format(date, "%Y"))
  ifelse(doy >= start_doy, calendar_year, calendar_year - 1)
}

# Climate aggregated by phenological year (not currently used in analysis)
fruit_pheno_climate <- weather %>%
  mutate(fruit_pheno_year = assign_pheno_year(Date, fruit_pheno_year$start_doy)) %>%
  group_by(fruit_pheno_year) %>%
  summarize(
    total_precip_mm = sum(Precip.mm,  na.rm = TRUE),
    mean_min_temp   = mean(MinTemp,   na.rm = TRUE),
    mean_max_temp   = mean(MaxTemp,   na.rm = TRUE),
    mean_temp       = (mean_min_temp + mean_max_temp) / 2,
    n_days          = n(),
    .groups         = "drop"
  ) %>%
  rename(pheno_year = fruit_pheno_year)


##########################
# Build modeling datasets

# Observation-level dataset for CLMM 
# Joins: per-tree lagged canopy fruit + seasonal lagged climate
# Used for models of canopy fruit and fallen fruit

model_data_obs <- BE_pheno_weather %>%
  mutate(
    Date2      = as.Date(Date2),
    water_year = ifelse(month(Date2) >= 10, year + 1, year),
    season     = factor(
      ifelse(month %in% dry_months, "dry", "wet"),
      levels = c("wet", "dry")),
    nuts.ground = as.ordered(nuts.ground),
    nuts.tree = as.ordered(nuts.tree),
    flwrs.ground = as.ordered(flwrs.ground)) %>%
  filter(Date2 >= as.Date("2017-10-01"), !is.na(nuts.ground)) %>%
  # Per-tree lagged canopy fruit score
  left_join(
    fruit_tree_seasonal %>%
      select(water_year, season, tree.no,
             canopy_fruit_prev_season, canopy_fruit_2prev_season),
    by = c("water_year", "season", "tree.no")
  ) %>%
  # Seasonal lagged climate (population-level)
  left_join(
    season_precip_lagged %>%
      select(water_year, season,
             total_precip_mm,
             mean_min_temp, mean_max_temp,
             precip_prev_season,
             temp_min_prev_season, temp_max_prev_season,
             precip_2prev_season,
             temp_min_2prev_season, temp_max_2prev_season),
    by = c("water_year", "season")
  ) %>%
  # Phenological-year climate (for canopy fruit model)
  mutate(fruit_pheno_year = assign_pheno_year(Date2, fruit_pheno_year$start_doy)) %>%
  left_join(
    fruit_pheno_climate %>%
      select(
        pheno_year,
        pheno_year_total_precip = total_precip_mm,
        pheno_year_mean_temp    = mean_temp,
        pheno_year_max_temp     = mean_max_temp,
        pheno_year_min_temp     = mean_min_temp
      ),
    by = c("fruit_pheno_year" = "pheno_year")
  )

# Standardize predictors (z-scores) 
climate_vars_to_scale <- c(
  "total_precip_mm",
  "mean_min_temp", "mean_max_temp",
  "precip_prev_season",
  "temp_min_prev_season", "temp_max_prev_season",
  "precip_2prev_season",
  "temp_min_2prev_season", "temp_max_2prev_season",
  "canopy_fruit_prev_season", "canopy_fruit_2prev_season",
  "pheno_year_total_precip", "pheno_year_min_temp", "pheno_year_max_temp"
)

scale_vars <- function(df, vars) {
  df %>%
    mutate(across(
      all_of(intersect(vars, names(df))),
      ~ scale(.)[, 1],
      .names = "{.col}_z"
    ))
}

model_data_obs <- scale_vars(model_data_obs, climate_vars_to_scale)

##################################
# Ordinal regressions using CLMMs

# Canopy fruit score 

clmm_canopy <- clmm(
  nuts.tree ~ total_precip_mm_z + precip_prev_season_z +
    mean_min_temp_z + mean_max_temp_z +
    temp_min_prev_season_z + temp_max_prev_season_z +
    (1 | tree.no),
  data = filter(model_data_obs, season == "dry"),
  link = "logit"
)
summary(clmm_canopy)

clmm_canopy_summ <- summary(clmm_canopy)
ci_canopy        <- confint(clmm_canopy)
n_thresholds_canopy <- length(clmm_canopy$alpha)   # number of threshold parameters

fr_canopy <- data.frame(
  Model    = "Canopy fruit score",
  Variable = c("Concurrent season precipitation",
               "Previous season precipitation",
               "Concurrent season minimum temperature",
               "Concurrent season maximum temperature",
               "Prior season minimum temperature",
               "Prior season maximum temperature"),
  Coefficient = coef(clmm_canopy)[-seq_len(n_thresholds_canopy)],
  Std_Error   = clmm_canopy_summ$coefficients[-seq_len(n_thresholds_canopy), "Std. Error"],
  Z_value     = clmm_canopy_summ$coefficients[-seq_len(n_thresholds_canopy), "z value"],
  P_value     = clmm_canopy_summ$coefficients[-seq_len(n_thresholds_canopy), "Pr(>|z|)"]
) %>%
  mutate(
    Odds_Ratio   = exp(Coefficient),
    OR_Lower_CI  = exp(ci_canopy[-seq_len(n_thresholds_canopy), 1]),
    OR_Upper_CI  = exp(ci_canopy[-seq_len(n_thresholds_canopy), 2])
  )

# Flower score (wet season only, seasonal climate)
clmm_flower <- clmm(
  flwrs.ground ~ total_precip_mm_z + precip_prev_season_z +
    mean_min_temp_z + mean_max_temp_z +
    temp_min_prev_season_z + temp_max_prev_season_z +
    (1 | tree.no),
  data = filter(model_data_obs, season == "wet"),
  link = "logit"
)

summary(clmm_flower)

clmm_flower_summ  <- summary(clmm_flower)
ci_flower         <- confint(clmm_flower)
n_thresholds_fl   <- length(clmm_flower$alpha)

fl_score <- data.frame(
  Model    = "Flower score",
  Variable = c(
    "Concurrent season precipitation",
    "Previous season precipitation",
    "Concurrent season minimum temperature",
    "Concurrent season maximum temperature",
    "Prior season minimum temperature",
    "Prior season maximum temperature"
  ),
  Coefficient = coef(clmm_flower)[-seq_len(n_thresholds_fl)],
  Std_Error   = clmm_flower_summ$coefficients[-seq_len(n_thresholds_fl), "Std. Error"],
  Z_value     = clmm_flower_summ$coefficients[-seq_len(n_thresholds_fl), "z value"],
  P_value     = clmm_flower_summ$coefficients[-seq_len(n_thresholds_fl), "Pr(>|z|)"]
) %>%
  mutate(
    Odds_Ratio  = exp(Coefficient),
    OR_Lower_CI = exp(ci_flower[-seq_len(n_thresholds_fl), 1]),
    OR_Upper_CI = exp(ci_flower[-seq_len(n_thresholds_fl), 2])
  )

# Fallen fruit score — climate + lagged canopy fruit control 
# Control for prior canopy fruit abundance per tree (accounting for
# the fact that fallen fruit directly reflects what was in the canopy)
clmm_fallen_null <- clmm(
  nuts.ground ~ canopy_fruit_prev_season_z +
    (1 | tree.no),
  data = model_data_obs,
  link = "logit"
)
summary(clmm_fallen_null)

clmm_fallen <- clmm(
  nuts.ground ~ canopy_fruit_prev_season_z +           # control: prior canopy fruit
    precip_prev_season_z +
    temp_min_prev_season_z + temp_max_prev_season_z +  # 1-season lag climate
    precip_2prev_season_z +
    temp_min_2prev_season_z + temp_max_2prev_season_z + # 2-season lag climate
    (1 | tree.no),
  data = model_data_obs,
  link = "logit"
)
summary(clmm_fallen)

# Model comparison
AIC(clmm_fallen_null, clmm_fallen)

clmm_fallen_summ  <- summary(clmm_fallen)
ci_fallen         <- confint(clmm_fallen)
n_thresholds_fg   <- length(clmm_fallen$alpha)

fr_ground_score <- data.frame(
  Model    = "Fallen fruit score",
  Variable = c(
    "Prior canopy fruit (control)",
    "Previous season precipitation",
    "Previous season minimum temperature",
    "Previous season maximum temperature",
    "2-season prior precipitation",
    "2-season prior minimum temperature",
    "2-season prior maximum temperature"
  ),
  Coefficient = coef(clmm_fallen)[-seq_len(n_thresholds_fg)],
  Std_Error   = clmm_fallen_summ$coefficients[-seq_len(n_thresholds_fg), "Std. Error"],
  Z_value     = clmm_fallen_summ$coefficients[-seq_len(n_thresholds_fg), "z value"],
  P_value     = clmm_fallen_summ$coefficients[-seq_len(n_thresholds_fg), "Pr(>|z|)"]
) %>%
  mutate(
    Odds_Ratio  = exp(Coefficient),
    OR_Lower_CI = exp(ci_fallen[-seq_len(n_thresholds_fg), 1]),
    OR_Upper_CI = exp(ci_fallen[-seq_len(n_thresholds_fg), 2])
  )

# Combine and round model output table
ordinal_summary <- bind_rows(fl_score, fr_canopy, fr_ground_score) %>%  
  mutate(across(where(is.numeric), ~ round(., 3)),
         P_value = round(P_value, 4))

#write.csv(ordinal_summary, "supp_table3.csv")



##########################
# clmm to test difference in fruit production between forest and edge/pasture trees

#combine BE_trees with model_data_obs
habitat_data_obs <- left_join(model_data_obs, BE_trees, by = c("tree.no"))
habitat_data_obs <- habitat_data_obs %>%
  filter(!is.na(tree.no))

### model scores based on habitat
# flower scores
m <- clmm(flwrs.ground ~ habitat + (1 | tree.no), data = habitat_data_obs)
summary(m)
n_thresholds_hb <- length(m$alpha)

# Make summary table
clmm_hb_summ  <- summary(m)
ci_hb        <- confint(m)
n_thresholds_hb   <- length(m$alpha)

habitat_table <- data.frame(
  Model = "Flower score",
  Variable = "Habitat",
  Coefficient = coef(m)[-seq_len(n_thresholds_hb)],
  Std_Error   = clmm_hb_summ$coefficients[-seq_len(n_thresholds_hb), "Std. Error"],
  Z_value     = clmm_hb_summ$coefficients[-seq_len(n_thresholds_hb), "z value"],
  P_value     = clmm_hb_summ$coefficients[-seq_len(n_thresholds_hb), "Pr(>|z|)"]
) %>%
  mutate(
    Odds_Ratio  = exp(Coefficient),
    OR_Lower_CI = exp(ci_hb[-seq_len(n_thresholds_hb), 1]),
    OR_Upper_CI = exp(ci_hb[-seq_len(n_thresholds_hb), 2])
  )

# Canopy fruit
m <- clmm(nuts.tree ~ habitat + (1 | tree.no), data = habitat_data_obs)
summary(m)
n_thresholds_hb <- length(m$alpha)

# Make summary table
clmm_hb_summ  <- summary(m)
ci_hb        <- confint(m)
n_thresholds_hb   <- length(m$alpha)

habitat_table <- rbind(habitat_table, data.frame(
  Model = "Canopy fruit",
  Variable = c("Habitat"),
  Coefficient = coef(m)[-seq_len(n_thresholds_hb)],
  Std_Error   = clmm_hb_summ$coefficients[-seq_len(n_thresholds_hb), "Std. Error"],
  Z_value     = clmm_hb_summ$coefficients[-seq_len(n_thresholds_hb), "z value"],
  P_value     = clmm_hb_summ$coefficients[-seq_len(n_thresholds_hb), "Pr(>|z|)"]
) %>%
  mutate(
    Odds_Ratio  = exp(Coefficient),
    OR_Lower_CI = exp(ci_hb[-seq_len(n_thresholds_hb), 1]),
    OR_Upper_CI = exp(ci_hb[-seq_len(n_thresholds_hb), 2])
  ))

# Fallen fruit
m <- clmm(nuts.ground ~ habitat + (1 | tree.no), data = habitat_data_obs)
summary(m)
n_thresholds_hb <- length(m$alpha)

# Make summary table
clmm_hb_summ  <- summary(m)
ci_hb        <- confint(m)
n_thresholds_hb   <- length(m$alpha)

habitat_table <- rbind(habitat_table, data.frame(
  Model = "Fallen fruit",
  Variable = c("Habitat"),
  Coefficient = coef(m)[-seq_len(n_thresholds_hb)],
  Std_Error   = clmm_hb_summ$coefficients[-seq_len(n_thresholds_hb), "Std. Error"],
  Z_value     = clmm_hb_summ$coefficients[-seq_len(n_thresholds_hb), "z value"],
  P_value     = clmm_hb_summ$coefficients[-seq_len(n_thresholds_hb), "Pr(>|z|)"]
) %>%
  mutate(
    Odds_Ratio  = exp(Coefficient),
    OR_Lower_CI = exp(ci_hb[-seq_len(n_thresholds_hb), 1]),
    OR_Upper_CI = exp(ci_hb[-seq_len(n_thresholds_hb), 2])
  ))

habitat_table <- habitat_table %>%  
  mutate(across(where(is.numeric), ~ round(., 3)))

#write.csv(habitat_table, file = "supp_table5.csv")

# Figure to show the time series of Brazil nut flower and fruit scores through time
# with seasonal precipitation
BE_time_series <- habitat_data_obs %>%
  mutate(
    flwrs.ground = as.numeric(levels(flwrs.ground))[flwrs.ground],
    nuts.tree  = as.numeric(levels(nuts.tree))[nuts.tree],
    nuts.ground  = as.numeric(levels(nuts.ground))[nuts.ground]
  ) %>%
  group_by(water_year, season, habitat) %>%
  summarize(
    flwrs.ground.mean = mean(flwrs.ground, na.rm = TRUE),
    flwrs.sd          = sd(flwrs.ground, na.rm = TRUE),
    nuts.tree.mean = mean(nuts.tree, na.rm = TRUE),
    nuts.tree.sd          = sd(nuts.tree, na.rm = TRUE),
    nuts.ground.mean       = mean(nuts.ground, na.rm = TRUE),
    nuts.ground.sd    = sd(nuts.ground, na.rm = TRUE),
    .groups = "drop"
  )

BE_time_series <- BE_time_series %>%
  mutate(flwrs.ground.mean = ifelse(season == "dry", NA, flwrs.ground.mean),
         nuts.ground.mean       = ifelse(season == "dry", NA, nuts.ground.mean))

BE_time_series <- BE_time_series %>%
  left_join(season_precip) %>%
  mutate(season = factor(season, levels = c("wet", "dry"))) %>%
  arrange(water_year, season)

BE_time_series$time <- interaction(BE_time_series$water_year, BE_time_series$season, sep = "-")
BE_time_series$time <- factor(BE_time_series$time, levels = unique(BE_time_series$time))

# Scale factors 
y_max        <- max(BE_time_series$nuts.ground.mean, na.rm = TRUE)
precip_scale <- y_max / max(BE_time_series$total_precip_mm, na.rm = TRUE)
temp_scale <- y_max / mean(BE_time_series$mean_min_temp, na.rm = TRUE)

fl_time_plot <- ggplot(BE_time_series, aes(x = time)) +
  geom_col(aes(y = total_precip_mm * precip_scale),fill = "lightblue") +
  #geom_errorbar(aes(ymin = nuts.ground.mean - nuts.ground.sd,
  #                  ymax = nuts.ground.mean + nuts.ground.sd),
  #              width = 0.2, color = "black", na.rm = TRUE) +
  geom_point(aes(y = flwrs.ground.mean, shape = habitat), fill = "#FFEB33", size = 5, na.rm = TRUE) +
  scale_shape_manual(values = c(21,24)) +
  scale_y_continuous(
    name = "Mean flower score",
    sec.axis = sec_axis(~ . / precip_scale / 2, name = "Precipitation (mm)")
  ) +
  labs(x = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold")) +
  guides(shape = guide_legend(title = "Forest type"),
         fill = "black") 

fr_time_plot <- ggplot(BE_time_series, aes(x = time)) +
  geom_col(aes(y = total_precip_mm * precip_scale),fill = "lightblue") +
  #geom_errorbar(aes(ymin = nuts.ground.mean - nuts.ground.sd,
  #                  ymax = nuts.ground.mean + nuts.ground.sd),
  #              width = 0.2, color = "black", na.rm = TRUE) +
  geom_point(aes(y = nuts.ground.mean, shape = habitat), fill = "#238443", size = 5, na.rm = TRUE) +
  scale_shape_manual(values = c(21,24)) +
  scale_y_continuous(
    name = "Mean fallen fruit score",
    sec.axis = sec_axis(~ . / precip_scale / 2, name = "Precipitation (mm)")) +
  labs(x = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold")) +
  guides(color = "none",
         shape = guide_legend(title = "Forest type"),
         fill = "none") 

fig4 <- ggarrange(fl_time_plot, fr_time_plot, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")
annotate_figure(fig4, bottom = text_grob("Year-season", size = 14))

combined_plot_data <- pivot_longer(
  data = BE_time_series,
  cols = c(flwrs.ground.mean, nuts.ground.mean), # Columns to reshape
  names_to = "Phenophase",     # New column for old column names
  values_to = "Score"       # New column for the values
)

#BE_score_calendar <- 
ggplot(combined_plot_data, aes(x = time, y = Score, fill = Phenophase, shape = habitat)) +
  geom_col(aes(x = time, y = total_precip_mm * precip_scale / 2), fill = "lightblue", inherit.aes = FALSE) +
  #geom_errorbar(aes(ymin = nuts.ground.mean - nuts.ground.sd,
  #                  ymax = nuts.ground.mean + nuts.ground.sd),
  #              width = 0.2, color = "black", na.rm = TRUE) +
  geom_point(size = 5, na.rm = TRUE) +
  coord_cartesian(ylim = c(0, 3)) +
  scale_fill_manual(values = c("#FFEB33", "#238443"), labels = c("Flowers", "Fallen fruit")) +
  scale_shape_manual(values = c(21,24)) +
  scale_y_continuous(
    name = "Mean score",
    sec.axis = sec_axis(~ . / precip_scale / 2, name = "Precipitation (mm)")) +
  labs(title = "", x = "Year-season") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title  = element_text(face = "italic"),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12)) +
  guides(shape = guide_legend(title = "Forest type"),
         fill = guide_legend(title = "Phase", override.aes = list(shape = 21))) 
