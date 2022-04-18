#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(nlme)
library(broom)
library(janitor)
library(patchwork)
library(pals)

# figure settings
source("code/settings/figure_settings.R")

# functions
source("code/generic-functions/continuous_time_interval.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")


#### edit data ####

# select hydrilla
# percent difference
inv_plant2 <- inv_plant %>% filter(CommonName == "Hydrilla") %>%
  mutate(PercCovered = PropCovered * 100,
         PercDiffCovered = PercCovered - InitPercCovered)

# first management year
min(inv_ctrl$GSYear)

# add management data
hydr_dat <- inv_plant2 %>%
  inner_join(inv_ctrl)

# add hydrilla cover from year before management
hydr_dat2 <- hydr_dat %>%
  filter(GSYear == min(inv_ctrl$GSYear) & !is.na(InitPercCovered)) %>%
  select(PermanentID) %>% # select waterbodies with previous year data
  mutate(GSYear = min(inv_ctrl$GSYear) - 1) %>%
  inner_join(inv_plant2 %>% # select data
               mutate(PercDiffCovered = NA_real_)) %>% # remove %change from year before
  full_join(hydr_dat) # combine

# histogram of hydrilla differences
hydr_dat2 %>%
  filter(PercDiffCovered < 0) %>%
  ggplot(aes(x = PercDiffCovered)) +
  geom_histogram(binwidth = 1)

# waterbodies with crashes >= 30%
hydr_dat3 <- hydr_dat2 %>%
  filter(PercDiffCovered <= -30) %>%
  select(PermanentID) %>%
  unique() %>%
  inner_join(hydr_dat2)

# complete time intervals
# surveys were not conducted every year on every lake for every species
# control data is implicitly complete -- missing interpreted as no control
hydr_time_cont <- hydr_dat3 %>%
  filter(PercDiffCovered <= -30) %>%
  select(PermanentID, GSYear) %>%
  unique() %>%
  mutate(out = pmap(., function(GSYear, PermanentID) 
    time_cont_fun(yearT = GSYear, permID = PermanentID, dat_in = hydr_dat3))) %>%
  unnest(cols = out)

# multiple time series for different crashes
hydr_time_cont %>%
  select(PermanentID, Year, YearDir) %>%
  unique() %>%
  get_dupes(PermanentID, YearDir) %>%
  left_join(hydr_time_cont) %>%
  left_join(hydr_dat3 %>%
              select(PermanentID, AreaName) %>%
              unique()) %>%
  arrange(AreaName, GSYear, Year)
# two waterbodies have two time series

# make year cut-offs wide
hydr_time_cont2 <- hydr_time_cont %>%
  select(PermanentID, Year, YearDir) %>%
  unique() %>%
  group_by(PermanentID, YearDir) %>%
  mutate(ID = 1:n()) %>%
  ungroup() %>%
  pivot_wider(names_from = YearDir,
              values_from = Year)

# check for duplicates
filter(hydr_time_cont2, ID > 1)

# select continuous time
# change prop treated to percent
# make lake name/county variable
hydr_dat4 <- hydr_dat3 %>%
  left_join(hydr_time_cont2) %>%
  filter(GSYear > before & GSYear < after ) %>%
  mutate(across(ends_with("PropTreated"), ~ .x * 100),
         AreaCounty = paste0(AreaName, " (", str_to_title(County), ")"),
         AreaCountyID = if_else(ID == 1, AreaCounty, paste(AreaCounty, ID))) %>%
  rename_with(~ str_replace(.x, "Prop", "Perc"), ends_with("PropTreated"))

# waterbodies
n_distinct(hydr_dat4$PermanentID)

perm_ids <- hydr_dat4 %>%
  select(PermanentID, AreaCounty) %>%
  unique() %>%
  arrange(AreaCounty)

# visualize
pdf("output/fwc_hydrilla_crash_exploratory_times_series.pdf")
for(i in 1:nrow(perm_ids)){
  
  subdat <- hydr_dat4 %>%
    filter(PermanentID == perm_ids$PermanentID[i])
  
  print(ggplot(subdat, aes(x = GSYear, y = PercCovered, color = AreaCountyID)) +
          geom_point() +
          geom_line() +
          labs(x = "Year", y = "Hydrilla PAC", title = perm_ids$AreaCounty[i]) +
          def_theme_paper +
          theme(legend.position = "none"))
  
}
dev.off()

# remove years without previous data
# arrange data by year
hydr_dat5 <- hydr_dat4 %>%
  filter(!is.na(InitPercCovered) & !is.na(PercDiffCovered)) %>%
  group_by(AreaCountyID) %>%
  arrange(GSYear) %>%
  ungroup()


#### initial visualizations ####

ggplot(hydr_dat5, aes(x = Lag1PercTreated, y = PercDiffCovered,
                      color = PermanentID)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ AreaCountyID, scales = "free") +
  def_theme_paper +
  theme(legend.position = "none")


#### fit models ####

# fit model to each lake
crash_mods <- hydr_dat5 %>%
  select(PermanentID, AreaCountyID, GSYear, PercDiffCovered, Lag1PercTreated) %>%
  nest(data = c(GSYear, PercDiffCovered, Lag1PercTreated)) %>%
  mutate(fit = map(data, ~lm(PercDiffCovered ~ Lag1PercTreated, data = .)))

# plot autocorrelations
pdf("output/fwc_hydrilla_crash_autocorrelations.pdf")
for(i in 1:nrow(crash_mods)){
  
  print(acf(residuals(crash_mods$fit[i][[1]]),
      main = crash_mods$AreaCountyID[i]))
  
  print(acf(residuals(crash_mods$fit[i][[1]]),
            type = "partial",
            main = crash_mods$AreaCountyID[i]))
  
}
dev.off()

#### start here ####


# all coefficients
crash_coef <- crash_mods %>%
  mutate(tidied1 = map(fit1, tidy),
         tidied2 = map(fit2, tidy),
         tidied3 = map(fit3, tidy)) %>%
  select(PermanentID, AreaCounty, tidied1, tidied2, tidied3) %>%
  unnest(c(tidied1, tidied2, tidied3),
         names_sep = "_") %>%
  filter(tidied1_term != "(Intercept)") %>% # remove intercepts
  pivot_longer(cols = tidied1_term:tidied3_p.value,
               names_to = c("lag", ".value"),
               names_pattern = "tidied(.)_(.*)")

# visualize coefficients
ggplot(crash_coef, aes(x = estimate, color = lag)) +
  geom_density()

# extreme values
filter(crash_coef, estimate < -5000 | estimate > 2000)
# unrealistic estimates from Rodman reservoir
# probably because it has very few treatments
# and they're very small

# non-finite values
filter(crash_coef, estimate == Inf | estimate == -Inf | is.na(estimate))
# Little Lake Sawgrass had no treatment

# remove unrealistic estimates
crash_coef2 <- crash_coef %>%
  filter(!(PermanentID %in% c("120053782", "164160423")))

# visualize coefficients
ggplot(crash_coef2, aes(x = estimate, color = lag)) +
  geom_density()

crash_coef2 %>%
  group_by(lag) %>%
  summarize(est = mean(estimate),
            est_sd = sd(estimate),
            p_sig = sum(p.value < 0.1 & estimate < 0))
# estimates become more negative with longer lag time
# fewer relationships are statistically significant

# remove lakes that couldn't be used to fit models properly
hydr_dat4 <- hydr_dat3 %>%
  filter(!(PermanentID %in% c("120053782", "164160423")))


#### visualization ####

# P values for figure
crash_coef3 <- crash_coef2 %>%
  filter(lag == 1) %>%
  add_column(hydr_dat4 %>%
              transmute(maxTreat = max(Lag1PercTreated),
                        minTreat = min(Lag1PercTreated),
                        maxDiff = max(PercDiffCovered),
                        minDiff = min(PercDiffCovered)) %>%
               unique()) %>%
  mutate(p.round = round_half_up(p.value, 3),
         p.text = paste0("italic('P')~`=`~", p.round),
         est.round = round_half_up(estimate, 3),
         est.text = paste0("est. = ", est.round))

# subset for figure
area_names <- sort(crash_coef3$AreaCounty)
crash_coef3_a <- crash_coef3 %>% filter(AreaCounty %in% area_names[1:32])
crash_coef3_b <- crash_coef3 %>% filter(AreaCounty %in% area_names[33:63])
hydr_dat4_a <- hydr_dat4 %>% filter(AreaCounty %in% area_names[1:32])
hydr_dat4_b <- hydr_dat4 %>% filter(AreaCounty %in% area_names[33:63])

# time series fig
ts_fig_a <- ggplot(hydr_dat4_a, aes(x = GSYear, y = PercCovered,
                      color = AreaCounty)) +
  geom_line(size = 0.5) +
  facet_grid(rows = vars(AreaCounty)) +
  labs(x = "Year", title = "PAC") +
  scale_color_manual(values = rep(kelly()[-c(1, 22)], 2)) +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

ts_fig_b <- ggplot(hydr_dat4_b, aes(x = GSYear, y = PercCovered,
                                    color = AreaCounty)) +
  geom_line(size = 0.5) +
  facet_grid(rows = vars(AreaCounty)) +
  labs(x = "Year", title = "PAC") +
  scale_color_manual(values = rep(kelly()[-c(1, 22)], 2)) +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

# model fit fig
mf_fig_a <- ggplot(hydr_dat4_a, aes(x = Lag1PercTreated, y = PercDiffCovered,
                      color = AreaCounty)) +
  geom_hline(yintercept = 0, size = 0.1, linetype = "dashed") +
  geom_smooth(method = "lm", formula = "y ~ x", size = 0.5) +
  geom_text(data = crash_coef3_a, 
            aes(x = maxTreat, y = maxDiff, label = p.round),
            hjust = 1, vjust = 1, size = 2, color = "black", fontface = "italic") +
  geom_text(data = crash_coef3_a, 
            aes(x = minTreat, y = maxDiff, label = est.round),
            hjust = 0, vjust = 1, size = 2, color = "black") +
  facet_grid(rows = vars(AreaCounty)) +
  labs(x = "Percent area managed", title = "Annual difference in PAC") +
  scale_color_manual(values = rep(kelly()[-c(1, 22)], 2)) +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text.y = element_text(size = 6, color="black",
                                    angle = 0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

mf_fig_b <- ggplot(hydr_dat4_b, aes(x = Lag1PercTreated, y = PercDiffCovered,
                                    color = AreaCounty)) +
  geom_hline(yintercept = 0, size = 0.1, linetype = "dashed") +
  geom_smooth(method = "lm", formula = "y ~ x", size = 0.5) +
  geom_text(data = crash_coef3_b, 
            aes(x = maxTreat, y = maxDiff, label = p.round),
            hjust = 1, vjust = 1, size = 2, color = "black", fontface = "italic") +
  geom_text(data = crash_coef3_b, 
            aes(x = minTreat, y = maxDiff, label = est.round),
            hjust = 0, vjust = 1, size = 2, color = "black") +
  facet_grid(rows = vars(AreaCounty)) +
  labs(x = "Percent area managed", title = "Annual difference in PAC") +
  scale_color_manual(values = rep(kelly()[-c(1, 22)], 2)) +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text.y = element_text(size = 6, color="black",
                                    angle = 0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

# combine and save
comb_fig <- ts_fig_a + mf_fig_a + ts_fig_b + mf_fig_b +
  plot_layout(nrow = 1)

ggsave("output/fwc_hydrilla_crash_time_series_prediction.png", comb_fig,
       device = "png", width = 7, height = 9, units = "in")


#### values for text ####

crash_sum <- crash_coef3 %>%
  mutate(summ = case_when(estimate < 0 & p.value < 0.1 ~ "sig neg",
                          estimate > 0 & p.value < 0.1 ~ "sig pos",
                          estimate < 0 & p.value >= 0.1 ~ "n.s. neg",
                          estimate > 0 & p.value >= 0.1 ~ "n.s. pos")) %>%
  count(summ)

write_csv(crash_sum, "output/fwc_hydrilla_crash_prediction_summary.csv")
