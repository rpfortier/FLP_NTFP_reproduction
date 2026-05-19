### Data analysis for NTFP phenology at FLP ###
### This script requires running the Brazil nut phenology script first

#### Read in the phenology data and clean dates ####

#### Mauritia flexuosa
#read in data
MF_pheno <- read.csv("MF_phenology_data.csv")

MF_pheno <- MF_pheno %>%
  rename(tree.no = 1,
         n.old.bracts = 6,
         n.old.fruit = 7,
         color.old.fruit = 8,
         n.new.bracts = 9,
         perc.fruit = 10,
         color.new.fruit = 11,
         flwrs.ground = 12)
MF_pheno[,c(13:26)] <- NULL
#change blank cells to NA and remove rows with all NA
MF_pheno[MF_pheno == ""] <- NA
MF_pheno <- MF_pheno[rowSums(is.na(MF_pheno)) != ncol(MF_pheno), ]

#reset rownames
rownames(MF_pheno) <- NULL

#change Date column to date format. 
MF_pheno$Date2 <- as.Date(parse_date_time(MF_pheno$Date, orders = c("dmy", "mdy")))

#manually change some weeks of dates that were entered incorrectly
MF_pheno[c(4282:4301),]$Date2 <- as.Date("2026-01-07")
MF_pheno[c(4342:4361),]$Date2 <- as.Date("2026-02-10")
MF_pheno[c(4382:4401),]$Date2 <- as.Date("2026-03-10")

#find dates that are NA
x <- MF_pheno[is.na(MF_pheno$Date2),]
#remove the NA rows
MF_pheno <- MF_pheno[!is.na(MF_pheno$Date2),]

#Change some dates
a <- ymd(MF_pheno$Date2)
MF_pheno$year <- year(a)
MF_pheno$month <- month(a)
MF_pheno$month <- month.abb[MF_pheno$month]
MF_pheno$julian <- format(MF_pheno$Date2, "%j")
MF_pheno$julian <- as.numeric(MF_pheno$julian)
MF_pheno$week <- isoweek(MF_pheno$Date2)

# Temporal summary per tree (for QC / visualization)
temporal_summary_MF <- MF_pheno %>%
  group_by(tree.no) %>%
  summarise(
    first_obs      = min(Date2),
    last_obs       = max(Date2),
    n_observations = n(),
    .groups = "drop"
  ) %>%
  arrange(first_obs)
temporal_summary_MF$sex <- ifelse(temporal_summary_MF$tree.no > 9, "Female", "Male")

MF_time_series <- ggplot(temporal_summary_MF, aes(y = tree.no, shape = sex, fill = sex)) +
  geom_segment(aes(x = first_obs, xend = last_obs, yend = tree.no), linewidth = 1) +
  geom_point(aes(x = first_obs), color = "#bd0026", size = 3) +
  geom_point(aes(x = last_obs),  color = "black",    size = 3) +
  scale_shape_manual(values = c(19,25)) +
  scale_fill_manual(values = c("#bd0026", "white")) +
  labs(title = "Mauritia flexuosa", x = "", y = "") +
  theme_minimal() +
  guides(shape = guide_legend(title = "Sex"),
         fill = guide_legend(title = "Sex")) +
  theme(axis.text.y = element_text(size = 8))

figS2 <- ggarrange(BE_time_series, MF_time_series, ncol = 1, nrow = 2)
annotate_figure(figS2, 
                left = text_grob("Tree ID", rot = 90, size = 14),
                bottom = text_grob("Date", size = 14))


###########################
#### Mauritia flexuosa ####

#combine weather data with phenology data
MF_pheno_weather <- full_join(MF_pheno, weather_weekly, by = c("year","week"))

#remove rows with no Date2
MF_pheno_weather <- MF_pheno_weather[!is.na(MF_pheno_weather$Date2),]

#make perc.fruit numeric
MF_pheno_weather$perc.fruit <- as.numeric(MF_pheno_weather$perc.fruit)
#change NAs to 0
MF_pheno_weather$perc.fruit <- ifelse(is.na(MF_pheno_weather$perc.fruit), 0, MF_pheno_weather$perc.fruit)

#clean flower data
x <- c("Muchos" , "muchos", "4")
MF_pheno_weather$flwrs.ground <- ifelse(MF_pheno_weather$flwrs.ground %in% x, 3, MF_pheno_weather$flwrs.ground)
x <- c("-" , "o", "?")
MF_pheno_weather$flwrs.ground <- ifelse(MF_pheno_weather$flwrs.ground %in% x, 0, MF_pheno_weather$flwrs.ground)
MF_pheno_weather$flwrs.ground <- ifelse(MF_pheno_weather$flwrs.ground == "Si", 1, MF_pheno_weather$flwrs.ground)
MF_pheno_weather$flwrs.ground <- as.numeric(MF_pheno_weather$flwrs.ground)
MF_pheno_weather$flwrs.ground <- ifelse(MF_pheno_weather$flwrs.ground > 3, 1, MF_pheno_weather$flwrs.ground)

#convert perc.fruit to an ordinal score
MF_pheno_weather <- MF_pheno_weather %>%
  mutate(fruit.score = case_when(
    perc.fruit == 0            ~ 0,
    perc.fruit > 0 & perc.fruit <= 200   ~ 1,
    perc.fruit > 200 & perc.fruit <= 400 ~ 2,
    perc.fruit > 400           ~ 3,
    TRUE ~ NA_real_
  ))

#clean fruit color
x <- c("No", "-")
MF_pheno_weather$color.new.fruit <- ifelse(MF_pheno_weather$tree.no <= 9, NA, 
                                           ifelse(MF_pheno_weather$color.new.fruit %in% x | is.na(MF_pheno_weather$color.new.fruit), 0, MF_pheno_weather$color.new.fruit))


######################################
# Circular statistics

#first convert each week to biweeks
MF_pheno_weather <- MF_pheno_weather %>%
  mutate(biweek_period = ceiling(week / 2))

#visualize with circular statistics
flower_counts <- MF_pheno_weather %>%
  mutate(biweek_period = ceiling(week / 2)) %>%
  filter(!is.na(flwrs.ground)) %>%
  filter(Date2 >= "2018-07-01") %>% #skip first year (no flower scores)
  group_by(biweek_period, flwrs.ground) %>%
  summarise(n_observations = n(), .groups = "keep") %>%
  group_by(biweek_period) %>%
  mutate(total_period_obs = sum(n_observations),
         prop_score = n_observations / total_period_obs) %>%
  ungroup()

#filter out male trees for the fruit scores
aguaje_female <- subset(MF_pheno_weather, tree.no > 9) 

aguaje_color <- aguaje_female %>%
  mutate(biweek_period = ceiling(week / 2),
         color.new.fruit = suppressWarnings(as.numeric(color.new.fruit))) %>%
  filter(!is.na(color.new.fruit)) %>%
  filter(color.new.fruit < 4) %>%
  group_by(biweek_period, color.new.fruit) %>%
  summarise(n_observations = n(), .groups = "keep") %>%
  group_by(biweek_period) %>%
  mutate(total_period_obs = sum(n_observations),
         prop_score = n_observations / total_period_obs) %>%
  ungroup()

aguaje_counts <- aguaje_female %>%
  mutate(biweek_period = ceiling(week / 2)) %>%
  filter(!is.na(fruit.score)) %>%
  group_by(biweek_period, fruit.score) %>%
  summarise(n_observations = n(), .groups = "keep") %>%
  group_by(biweek_period) %>%
  mutate(total_period_obs = sum(n_observations),
         prop_score = n_observations / total_period_obs) %>%
  ungroup()

###quantify CHANGE in color to find a seasonal peak in fruit maturation since mature fruit can stay on the plant a long time.
#ensure each tree's record is in chronological order
aguaje_female <- aguaje_female %>%
  arrange(Date2, tree.no)

aguaje_female <- aguaje_female %>%
  group_by(tree.no) %>%
  mutate(color.new.fruit = as.numeric(color.new.fruit),
         color.new.fruit = ifelse(color.new.fruit > 3, 3, color.new.fruit),
         color.new.fruit = ifelse(is.na(color.new.fruit), 0, color.new.fruit),
         maturation_score = ifelse(color.new.fruit > lag(color.new.fruit, 1), 1, 0),
         maturation_score = ifelse(is.na(maturation_score), 0, maturation_score)) %>%
  ungroup()

aguaje_maturation <- aguaje_female %>%
  mutate(biweek_period = ceiling(week / 2)) %>%
  group_by(biweek_period, maturation_score) %>%
  summarise(n_observations = n(), .groups = "keep") %>%
  group_by(biweek_period) %>%
  mutate(total_period_obs = sum(n_observations),
         prop_score = n_observations / total_period_obs) %>%
  ungroup()

# For plotting, convert biweek_period back to a week equivalent
# This uses the middle week of each bi-weekly period for positioning
flower_counts <- flower_counts %>%
  mutate(week_equiv = (biweek_period - 1) * 2 + 1.5) %>%  # Centers each period
  left_join(precip_biweekly, by = "biweek_period")
aguaje_color <- aguaje_color %>%
  mutate(week_equiv = (biweek_period - 1) * 2 + 1.5) %>%   # Centers each period
  left_join(precip_biweekly, by = "biweek_period")
aguaje_maturation <- aguaje_maturation %>%
  mutate(week_equiv = (biweek_period - 1) * 2 + 1.5) %>%   # Centers each period
  left_join(precip_biweekly, by = "biweek_period")
aguaje_counts <- aguaje_counts %>%
  mutate(week_equiv = (biweek_period - 1) * 2 + 1.5) %>%   # Centers each period
  left_join(precip_biweekly, by = "biweek_period")

#assign levels to  score
flower_counts$flwrs.ground <- factor(flower_counts$flwrs.ground, levels = c("0","1","2","3"))
aguaje_color$color.new.fruit <- factor(aguaje_color$color.new.fruit, levels = c("0","1","2","3"))
aguaje_counts$fruit.score <- factor(aguaje_counts$fruit.score, levels = c("0","1","2","3"))

#month breaks for plotting
month_breaks <- seq(0.5, 51, 4.3) 
month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
#add same month breaks to monthly weather
weather_monthly <- weather_monthly %>%
  mutate(week_equiv = c(month_breaks[1:11] + 2.15, 47.8 +2.35),
         bar_width = c(4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.75))

score_colors <- c("0" = "transparent", "1" = "#ffffb2", "2" = "#fd8d3c", "3" = "#bd0026")
maturation_colors <- c("0" = "transparent", "1" = "#ffffb2", "2" = "#fd8d3c", "3" = "#bd0026")

#Flower plot
MF_fl <- ggplot(subset(flower_counts, flwrs.ground != 0), 
             aes(x = week_equiv, y = prop_score, fill = as.factor(flwrs.ground))) +
  # Add precipitation as a background layer (scaled to fit)
  geom_col(data = weather_monthly, aes(y = precip.mm * 0.00125, fill = NULL, width = bar_width), 
           fill = "lightblue", alpha = 1) +
  geom_hline(yintercept = c(0, 0.25, 0.5), 
             color = "gray70", size = 0.3) +
  geom_vline(xintercept = month_breaks, 
             color = "gray70", size = 0.3) +
  geom_col(position = position_stack(reverse = TRUE), 
           color = "black", width = 2, alpha = 1) +
  coord_polar(theta = "x", start = 0) +
  scale_fill_manual(values = score_colors, name = "Score") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5), 
                     limits = c(0, 0.5), 
                     sec.axis = sec_axis(~ . / 0.00125, name = "")) +
  scale_x_continuous(breaks = month_breaks, labels = month_labels) +
  labs(
    title = "Fallen flowers",
    x = "",
    y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "gray70", size = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
MF_fl

#Fruit color plot
MF_fr <- ggplot(subset(aguaje_color, color.new.fruit != 0), 
       aes(x = week_equiv, y = prop_score, fill = as.factor(color.new.fruit))) +
  # Add precipitation as a background layer (scaled to fit)
  geom_col(data = weather_monthly, aes(y = precip.mm * 0.0025, fill = NULL, width = bar_width), 
           fill = "lightblue", alpha = 1) +
  geom_hline(yintercept = c(0, 0.25, 0.5, 0.75, 1), 
             color = "gray70", size = 0.3) +
  geom_vline(xintercept = month_breaks, 
             color = "gray70", size = 0.3) +
  geom_col(position = position_stack(reverse = TRUE), 
           color = "black", width = 2, alpha = 1) +
  coord_polar(theta = "x", start = 0) +
  scale_fill_manual(values = maturation_colors, name = "Color score") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     limits = c(0, 1), 
                     sec.axis = sec_axis(~ . / 0.0025, name = "")) +
  scale_x_continuous(breaks = month_breaks, labels = month_labels) +
  labs(
    title = "Fruit color",
    x = "",
    y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "gray70", size = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
MF_fr

#fruit maturation 
MF_fr2 <- ggplot(subset(aguaje_maturation, maturation_score != 0), 
              aes(x = week_equiv, y = prop_score)) +
  # Add precipitation as a background layer (scaled to fit)
  geom_col(data = weather_monthly, aes(y = precip.mm * 0.00125, fill = NULL, width = bar_width), 
           fill = "lightblue", alpha = 1) +
  geom_hline(yintercept = c(0, 0.125, 0.25, 0.375, 0.5), 
             color = "gray70", size = 0.3) +
  geom_vline(xintercept = month_breaks, 
             color = "gray70", size = 0.3) +
  geom_col(position = position_stack(reverse = TRUE), 
           color = "black", fill = "#fd8d3c", width = 2, alpha = 1) +
  coord_polar(theta = "x", start = 0) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5), 
                     limits = c(0, 0.5), 
                     sec.axis = sec_axis(~ . / 0.00125, name = "")) +
  scale_x_continuous(breaks = month_breaks, labels = month_labels) +
  labs(
    title = "Fruit maturation",
    x = "",
    y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "gray70", size = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
MF_fr2

#fruit score
ggplot(subset(aguaje_counts, fruit.score != 0), 
                aes(x = week_equiv, y = prop_score, fill = as.factor(fruit.score))) +
  # Add precipitation as a background layer (scaled to fit)
  geom_col(data = weather_monthly, aes(y = precip.mm * 0.0025, fill = NULL, width = bar_width), 
           fill = "lightblue", alpha = 1) +
  geom_hline(yintercept = c(0, 0.25, 0.5, 0.75, 1), 
             color = "gray70", size = 0.3) +
  geom_vline(xintercept = month_breaks, 
             color = "gray70", size = 0.3) +
  geom_col(position = position_stack(reverse = TRUE), 
           color = "black", width = 2, alpha = 1) +
  coord_polar(theta = "x", start = 0) +
  scale_fill_manual(values = score_colors, name = "Score") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     limits = c(0, 1), 
                     sec.axis = sec_axis(~ . / 0.0025, name = "")) +
  scale_x_continuous(breaks = month_breaks, labels = month_labels) +
  labs(
    title = "Fruiting structures",
    x = "",
    y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "gray70", size = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )


#combine graphs together
fig1b <- ggarrange(MF_fl, MF_fr, ncol = 1, nrow = 2, align = "v", common.legend = TRUE, labels = c("C", "D"))
fig1b
fig1 <- ggarrange(fig1a, fig1b, ncol = 2, nrow = 1)
annotate_figure(fig1, 
                top = text_grob("Brazil nut                                            Aguaje", face = "bold", size = 16),
                left = text_grob("Proportion of trees with score", rot = 90, size = 14),
                right = text_grob("Precipitation (mm)", rot = 270, size = 14))

#ggsave(filename = "fig1.png", width = 8.5, height = 8.5, units = "in", dpi = 300)

fig2 <- ggarrange(BE_fr2, MF_fr2, ncol = 2, nrow = 1, labels = c("A", "B"))
annotate_figure(fig2,
                top = text_grob("Brazil nut                                                       Aguaje", face = "bold", size = 16),
                left = text_grob("Proportion of trees", rot = 90, size = 14),
                right = text_grob("Precipitation (mm)", rot = 270, size = 14))


######################################
# Circular statistics

#convert julian day to radians on the MF_pheno_weather dataframe
MF_pheno_weather$radians <- circular(MF_pheno_weather$julian * (2 * pi / 365), 
                                     type = "angles",
                                     units = "radians",
                                     template = "geographics",
                                     modulo = "2pi") 
aguaje_female$radians <- circular(aguaje_female$julian * (2 * pi / 365), 
                                  type = "angles",
                                  units = "radians",
                                  template = "geographics",
                                  modulo = "2pi") 

# Flowers
MF_fl_data <- MF_pheno_weather %>% 
  filter(flwrs.ground > 0) %>%
  filter(Date2 >= "2018-07-01") #skip first year (no data)
fl_circular <- circular_stats(MF_fl_data, "radians", "Flowers")

fl_annual <- MF_pheno_weather %>%
  mutate(score = as.numeric(flwrs.ground)) %>%
  group_by(year, tree.no) %>%
  summarize(mean_flower_score = mean(score, na.rm = TRUE), .groups = "drop")
fl_masting <- masting_stats(fl_annual, "mean_flower_score")

# Fruit color
MF_fr_data <- aguaje_female %>% 
  filter(fruit.score > 0) %>%
  filter(Date2 >= "2018-07-01") #skip first year
fr_circular <- circular_stats(MF_fr_data, "radians", "Fruit color")

fr_annual <- aguaje_female %>%
  mutate(score = as.numeric(fruit.score)) %>%
  group_by(year, tree.no) %>%
  summarize(mean_fr_score = mean(score, na.rm = TRUE), .groups = "drop")
fr_masting <- masting_stats(fr_annual, "mean_fr_score")

# Fruit maturation
MF_fm_data <- aguaje_female %>% 
  filter(maturation_score > 0) %>%
  filter(Date2 >= "2018-07-01") #skip first year
fm_circular <- circular_stats(MF_fm_data, "radians", "Fruit maturation")

fm_annual <- aguaje_female %>%
  mutate(score = as.numeric(maturation_score)) %>%
  group_by(year, tree.no) %>%
  summarize(mean_fm_score = mean(score, na.rm = TRUE), .groups = "drop")
fm_masting <- masting_stats(fm_annual, "mean_fm_score")

# Combine circular results table
MF_circular_results <- bind_rows(fl_circular, fr_circular, fm_circular)
MF_circular_results$species <- "Aguaje"
circular_results <- rbind(circular_results, MF_circular_results)
#write.csv(circular_results, file = "supp_table1.csv")

# Combine masting metrics table
masting_results_MF <- data.frame(
  Species   = "Aguaje",
  Phenophase = c("Flowers", "Fruit score", "Fruit maturation"),
  CVp  = c(fl_masting["CVp"],  fr_masting["CVp"],  fm_masting["CVp"]),
  xCVi = c(fl_masting["xCVi"], fr_masting["xCVi"], fm_masting["xCVi"]),
  xPCC = c(fl_masting["xPCC"], fr_masting["xPCC"], fm_masting["xPCC"])
)
masting_results <- rbind(masting_results_BE, masting_results_MF)
#write.csv(masting_results, file = "supp_table2.csv")


#time series plot
ggplot(subset(MF_pheno_weather, tree.no > 9), aes(x = Date2)) +
  geom_col(data = weather, aes(x = Date, y = Precip.mm.avg * 0.25), alpha = 1, width = 1.5, na.rm = TRUE, fill = "skyblue") +
  # Individual tree flower production lines
  geom_line(aes(y = as.numeric(fruit.score), group = tree.no, color = as.factor(tree.no)), 
            alpha = 0.3, size = 0.8) +
  #geom_smooth(aes(y = as.numeric(flwrs.ground)), 
  #            method = "loess", span = 0.3, 
  #            color = "black", size = 1.5, alpha = 0.8) +
  scale_y_continuous(
    name = "Score",
    limits = c(0, 3),
    sec.axis = sec_axis(~ . / 0.25, name = "Precipitation (mm/biweek)")
  ) +
  # Color palette for individual trees
  #scale_color_viridis_d(name = "Tree ID", guide = "none") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Date", 
       title = "") +
  theme_minimal() 

ggplot(aguaje_female, aes(x = Date2)) +
  geom_col(data = weather, aes(x = Date, y = Precip.mm.avg * 0.1), alpha = 1, width = 1.5, na.rm = TRUE, fill = "skyblue") +
  # Individual tree flower production lines
  geom_line(aes(y = as.numeric(flwrs.ground), group = tree.no, color = as.factor(tree.no)), 
            alpha = 0.3, size = 0.8) +
  #geom_smooth(aes(y = as.numeric(flwrs.ground)), 
  #            method = "loess", span = 0.3, 
  #            color = "black", size = 1.5, alpha = 0.8) +
  scale_y_continuous(
    name = "Score",
    sec.axis = sec_axis(~ . / 0.1, name = "Precipitation (mm/biweek)")
  ) +
  # Color palette for individual trees
  #scale_color_viridis_d(name = "Tree ID", guide = "none") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Date", 
       title = "") +
  theme_minimal() 


######################################################
# Calculate phenological year for  fruit scores

fruit_pheno_year <- calculate_phenological_year(
  dates = aguaje_female$Date2,
  scores = as.numeric(aguaje_female$fruit.score)
)

# Calculate phenological years (not currently used in analysis)
weather_pheno <- weather %>%
  mutate(fruit_pheno_year = assign_pheno_year(Date, fruit_pheno_year$start_doy))

fruit_pheno_climate <- weather_pheno %>%
  group_by(fruit_pheno_year) %>%
  summarize(
    total_precip_mm = sum(Precip.mm, na.rm = TRUE),
    mean_min_temp = mean(MinTemp, na.rm = TRUE),
    mean_max_temp = mean(MaxTemp, na.rm = TRUE),
    mean_temp = (mean_min_temp + mean_max_temp) / 2,
    n_days = n(),
    .groups = "drop"
  ) %>%
  rename(pheno_year = fruit_pheno_year)

##########################
# Build modeling datasets

#match climate with MF phenology data and create new df for modeling
model_data_MF <- MF_pheno_weather %>%
  mutate(water_year = ifelse(month(Date2) >= 10, year(Date2) + 1, year(Date2)),
         season = ifelse(month %in% dry_months, "dry", "wet")) %>%
  left_join(season_precip_lagged, by = c("water_year", "season")) 

#model_data_MF <- model_data_MF %>%
#  left_join(
#    fruit_pheno_climate %>%
#      select(pheno_year,
#             fruit_total_precip = total_precip_mm,
#             fruit_mean_temp = mean_temp,
#             fruit_max_temp = mean_min_temp,
#             fruit_min_temp = mean_max_temp),
#    by = c("fruit_pheno_year" = "pheno_year"))


model_data_MF <- model_data_MF %>%
  mutate(
    # Ensure scores are correct format
    flwrs.ground = factor(flwrs.ground, 
                       levels = 0:3, 
                       ordered = TRUE),
    fruit.score = factor(fruit.score,
                         levels = 0:3,
                         ordered = TRUE),
    color.new.fruit = as.numeric(color.new.fruit),
    color.new.fruit = ifelse(is.na(color.new.fruit), 0, color.new.fruit),
    color.new.fruit = ifelse(color.new.fruit >3, 3, color.new.fruit),
    color.new.fruit = factor(color.new.fruit,
                         levels = 0:3,
                         ordered = TRUE),
    # Tree ID as factor
    tree_id = as.factor(tree.no),
    # Standardize climate variables for better model convergence
    precip_current_z = scale(total_precip_mm)[,1],
    precip_prev_season_z = scale(precip_prev_season)[,1],
    temp_current_z = scale(mean_temp)[,1],
    temp_min_current_z = scale(MinTemp)[,1],
    temp_max_current_z = scale(MaxTemp)[,1],
    temp_min_prev_season_z = scale(temp_min_prev_season)[,1],
    temp_max_prev_season_z = scale(temp_max_prev_season)[,1],
#    fruit_precip_z = scale(fruit_total_precip)[,1],
#    fruit_temp_z = scale(fruit_mean_temp)[,1],
#    fruit_max_temp_z = scale(fruit_max_temp)[,1],
#    fruit_min_temp_z = scale(fruit_min_temp)[,1]
  ) %>%
  filter(!is.na(flwrs.ground) & !is.na(tree.no))


##################################
# Ordinal regressions using CLMMs

# Flower score (highest during wet season)
clmm_flower <- clmm(flwrs.ground ~ precip_current_z + precip_prev_season_z + temp_min_current_z + temp_max_current_z + temp_min_prev_season_z + temp_max_prev_season_z +
                        (1 | tree.no),
                      data = filter(model_data_MF, season == "wet"),
                      link = "logit")

summary(clmm_flower)
summary <- summary(clmm_flower)

fl <- data.frame(
  Model = "Flower score",
  Variable = c("Precipitation current season","Precipitation prior season", "Minimum temperature current season", "Maximum temperature current season", "Minimum temperature prior season", "Maximum temperature prior season"),
  Coefficient = coef(clmm_flower)[-c(1:3)],
  Std_Error = summary$coefficients[-c(1:3), "Std. Error"],
  Z_value = summary$coefficients[-c(1:3), "z value"],
  P_value = summary$coefficients[-c(1:3), "Pr(>|z|)"]
)

# Add odds ratios and confidence intervals
ci <- confint(clmm_flower)
fl$Odds_Ratio <- exp(fl$Coefficient)
fl$OR_Lower_CI <- exp(ci[-c(1:3), 1])
fl$OR_Upper_CI <- exp(ci[-c(1:3), 2])

# Fruit score
clmm_fruit <- clmm(fruit.score ~ precip_current_z + precip_prev_season_z + temp_min_current_z + temp_max_current_z + temp_min_prev_season_z + temp_max_prev_season_z +
                        (1 | tree.no),
  data = filter(model_data_MF, tree.no > 9 & season == "wet"),
  link = "logit")

summary(clmm_fruit)
summary <- summary(clmm_fruit)

fr <- data.frame(
  Model = "Fruit score",
  Variable = c("Precipitation current season","Precipitation prior season", "Minimum temperature current season", "Maximum temperature current season", "Minimum temperature prior season", "Maximum temperature prior season"),
  Coefficient = coef(clmm_fruit)[-c(1:3)],
  Std_Error = summary$coefficients[-c(1:3), "Std. Error"],
  Z_value = summary$coefficients[-c(1:3), "z value"],
  P_value = summary$coefficients[-c(1:3), "Pr(>|z|)"]
)

# Add odds ratios and confidence intervals
ci <- confint(clmm_fruit)
fr$Odds_Ratio <- exp(fr$Coefficient)
fr$OR_Lower_CI <- exp(ci[-c(1:3), 1])
fr$OR_Upper_CI <- exp(ci[-c(1:3), 2])

#combine fl and fr
ordinal_summary_MF <- rbind(fl, fr)

#round
ordinal_summary_MF <- ordinal_summary_MF %>%
  mutate(Coefficient = round(Coefficient, 3),
         Std_Error = round(Std_Error, 3),
         Z_value = round(Z_value, 3),
         P_value = round(P_value, 3),
         Odds_Ratio = round(Odds_Ratio, 3),
         OR_Lower_CI = round(OR_Lower_CI, 3),
         OR_Upper_CI = round(OR_Upper_CI, 3))

#write.csv(ordinal_summary_MF, file = "supp_table4.csv")


# Figure to show time series of aguaje flower and fruit scores through time
# with seasonal precipitation
MF_time_series <- model_data_MF %>%
  mutate(
    flwrs.ground = as.numeric(levels(flwrs.ground))[flwrs.ground],
    fruit.score  = as.numeric(levels(fruit.score))[fruit.score]
  ) %>%
  group_by(water_year, season) %>%
  summarize(
    flwrs.ground.mean = mean(flwrs.ground, na.rm = TRUE),
    flwrs.sd          = sd(flwrs.ground, na.rm = TRUE),
    fruit.score.mean = mean(fruit.score, na.rm = TRUE),
    fruit.score.sd          = sd(fruit.score, na.rm = TRUE),
    .groups = "drop"
  )

MF_time_series <- MF_time_series %>%
  mutate(flwrs.ground.mean = ifelse(season == "dry", NA, flwrs.ground.mean))

MF_time_series <- MF_time_series %>%
  left_join(season_precip) %>%
  mutate(season = factor(season, levels = c("wet", "dry"))) %>%
  arrange(water_year, season)

MF_time_series$time <- interaction(MF_time_series$water_year, MF_time_series$season, sep = "-")
MF_time_series$time <- factor(MF_time_series$time, levels = unique(MF_time_series$time))

# Scale factors 
y_max        <- max(MF_time_series$fruit.score.mean, na.rm = TRUE)
precip_scale <- y_max / max(MF_time_series$total_precip_mm, na.rm = TRUE)
temp_scale <- y_max / mean(MF_time_series$mean_min_temp, na.rm = TRUE)

combined_plot_data_MF <- pivot_longer(
  data = MF_time_series,
  cols = c(flwrs.ground.mean, fruit.score.mean), # Columns to reshape
  names_to = "Phenophase",     # New column for old column names
  values_to = "Score"       # New column for the values
  ) %>%
  filter(water_year > 2017)

#MF_score_calendar <- 
ggplot(combined_plot_data_MF, aes(x = time, y = Score, fill = Phenophase)) +
  geom_col(aes(x = time, y = total_precip_mm * precip_scale / 2), fill = "steelblue", alpha = 0.4, inherit.aes = FALSE) +
  #geom_errorbar(aes(ymin = nuts.ground.mean - nuts.ground.sd,
  #                  ymax = nuts.ground.mean + nuts.ground.sd),
  #              width = 0.2, color = "black", na.rm = TRUE) +
  geom_point(size = 5, shape = 21, na.rm = TRUE) +
  #coord_cartesian(ylim = c(0, 3)) +
  scale_fill_manual(values = c("#ffffb2", "#bd0026"), labels = c("Flowers", "Fruit")) +
  scale_y_continuous(
    name = "Mean score",
    sec.axis = sec_axis(~ . / precip_scale / 2, name = "Precipitation (mm)")) +
  labs(title = "Mauritia flexuosa", x = "Year-season") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "italic")) +
  guides(fill = guide_legend(title = "M. flexuosa score", override.aes = list(shape = 21))) 
