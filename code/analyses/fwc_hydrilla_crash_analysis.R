#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(nlme)
library(broom.mixed)
library(janitor)
library(patchwork)
library(pals)
library(car)

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
# (including first year added on for visualization)
# arrange data by year
hydr_dat5 <- hydr_dat4 %>%
  filter(!is.na(InitPercCovered) & !is.na(PercDiffCovered)) %>%
  group_by(AreaCountyID) %>%
  arrange(GSYear) %>%
  ungroup()


#### initial visualizations ####

ggplot(hydr_dat5, aes(x = PercTreated, y = PercDiffCovered,
                      color = PermanentID)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ AreaCountyID, scales = "free") +
  def_theme_paper +
  theme(legend.position = "none")


#### fit models ####

# fit model to each lake
crash_mods <- hydr_dat5 %>%
  select(PermanentID, AreaCountyID, GSYear, PercDiffCovered, 
         PercTreated, Lag2PercTreated, Lag3PercTreated) %>%
  nest(data = c(GSYear, PercDiffCovered, 
                PercTreated, Lag2PercTreated, Lag3PercTreated)) %>%
  mutate(fit1 = map(data, ~lm(PercDiffCovered ~ PercTreated, data = .)),
         fit2 = map(data, ~lm(PercDiffCovered ~ Lag2PercTreated, data = .)),
         fit3 = map(data, ~lm(PercDiffCovered ~ Lag3PercTreated, data = .))) %>%
  arrange(AreaCountyID)

# coefficients
crash_coef <- crash_mods %>%
  mutate(tidied1 = map(fit1, tidy),
         tidied2 = map(fit2, tidy),
         tidied3 = map(fit3, tidy)) %>%
  select(PermanentID, AreaCountyID, tidied1, tidied2, tidied3) %>%
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

crash_mods2 <- crash_mods %>%
  filter(!(PermanentID %in% c("120053782", "164160423")))

# visualize coefficients
ggplot(crash_coef2, aes(x = estimate, color = lag)) +
  geom_density()

crash_coef2 %>%
  group_by(lag) %>%
  summarize(est = mean(estimate),
            est_sd = sd(estimate),
            p_sig = sum(p.value < 0.1 & estimate < 0))
# estimates similar for lag1 and lag3
# fewer relationships are statistically significant

# plot autocorrelations
pdf("output/fwc_hydrilla_crash_autocorrelations.pdf")
for(i in 1:nrow(crash_mods2)){
  
  temp_resid <- residuals(crash_mods2$fit1[i][[1]])
  temp_name <- crash_mods2$AreaCountyID[i]
  
  print(plot(x = 1:length(temp_resid),
             y = temp_resid,
             type = "l",
             main = temp_name))
  
  print(acf(temp_resid,
            main = temp_name))
  
  print(acf(temp_resid,
            type = "partial",
            main = temp_name))
  
}
dev.off()
# most have no significant autocorrelation
# some have sig negative, indicating more fluctuations
# around zero than expected from random
# none show wave-like patterns

# test for autocorrelation
crash_auto <- crash_mods2 %>%
  mutate(autocor = map(fit1, ~durbinWatsonTest(.x, max.lag = 2)$r), # extract coefficient
         pval = map(fit1, ~durbinWatsonTest(.x, max.lag = 2)$p)) %>% # and p-value
  select(PermanentID, AreaCountyID, autocor, pval) %>%
  unnest(c(autocor, pval)) %>%
  mutate(lag = rep(1:2, n_distinct(AreaCountyID))) # add lag number

# check for sig autocorrelations
(crash_auto_sig <- crash_auto %>%
    filter(pval < 0.05))
# most are negative

get_dupes(crash_auto_sig, AreaCountyID) # only 1 lag is ever sig

# refit models with autocorrelated errors
crash_mods3 <- crash_mods2 %>%
  mutate(fit1b = map(data, ~gls(PercDiffCovered ~ PercTreated, data = .,
                               correlation = NULL, method = "ML")),
         fit1c = map(data, ~gls(PercDiffCovered ~ PercTreated, data = .,
                               correlation = corAR1(), method = "ML")),
         CorComp = map2(fit1b, fit1c, ~anova(.x, .y)),
         DeltaAIC = map(CorComp, ~.x$AIC[2] - .x$AIC[1]),
         CorCompP = map(CorComp, ~.x$p[2])) %>%
  unnest(c(DeltaAIC, CorCompP)) %>%
  left_join(crash_auto_sig) %>% # visually compared to anova -- generally, but not completely consistent
  mutate(fitFin = case_when(CorCompP >= 0.05 ~ fit1b, # use simpler model when models are the same
                            CorCompP < 0.05 & DeltaAIC < 0 ~ fit1c, # use autocorrelation when significantly better
                            CorCompP < 0.05 & DeltaAIC > 0 ~ fit1b), # use random correlation when significantly better
         corAR = case_when(CorCompP >= 0.05 ~ "",
                           CorCompP < 0.05 & DeltaAIC < 0 ~ "1",
                           CorCompP < 0.05 & DeltaAIC > 0 ~ ""))

# coefficients after correcting for autocorrelation
crash_coef3 <- crash_mods3 %>%
  mutate(tidied = map(fitFin, tidy)) %>%
  select(PermanentID, AreaCountyID, tidied, corAR) %>%
  unnest(tidied) %>%
  filter(term != "(Intercept)") # remove intercepts


#### visualization ####

# edit raw data
hydr_dat6 <- hydr_dat4 %>% # use 4 to include initial year, this is for visualization only
  filter(!(PermanentID %in% c("120053782", "164160423")) & !is.na(PercCovered))

# P values for figure
# facet labels to indicate autocorrelation
crash_coef4 <- crash_coef3 %>%
  mutate(AreaCountyID2 = factor(AreaCountyID,
                               labels = paste0(str_replace_all(AreaCountyID, 
                                                               c("\\, Lake" = "", 
                                                                 " " = "~", 
                                                                 "\\'" = "")), 
                                               "^{", corAR, "}"))) %>%
  add_column(hydr_dat6 %>%
              transmute(maxTreat = max(PercTreated, na.rm = T),
                        minTreat = min(PercTreated, na.rm = T),
                        maxDiff = max(PercDiffCovered, na.rm = T),
                        minDiff = min(PercDiffCovered, na.rm = T)) %>%
               unique()) %>%
  mutate(p.round = round_half_up(p.value, 3),
         est.round = round_half_up(estimate, 3))

# facet labels to indicate autocorrelation
hydr_dat7 <- hydr_dat6 %>%
  left_join(crash_coef4 %>%
              select(AreaCountyID, AreaCountyID2))

# subset for figure
area_names <- sort(crash_coef4$AreaCountyID)
crash_coef4_a <- crash_coef4 %>% filter(AreaCountyID %in% area_names[1:33])
crash_coef4_b <- crash_coef4 %>% filter(AreaCountyID %in% area_names[34:65])
ts_dat_1 <- hydr_dat7 %>% filter(AreaCountyID %in% area_names[1:33])
ts_dat_2 <- hydr_dat7 %>% filter(AreaCountyID %in% area_names[34:65])
ts_dat_3 <- hydr_dat7 %>% filter(AreaCountyID %in% area_names[1:33] & !is.na(PercDiffCovered))
ts_dat_4 <- hydr_dat7 %>% filter(AreaCountyID %in% area_names[34:65] & !is.na(PercDiffCovered))

# time series fig
ts_fig_a <- ggplot(ts_dat_1, aes(x = GSYear, y = PercCovered,
                      color = AreaCountyID)) +
  geom_line(size = 0.5) +
  facet_grid(rows = vars(AreaCountyID2)) +
  labs(x = "Year", title = "PAC") +
  scale_color_manual(values = rep(kelly()[-c(1, 22)], 2)) +
  scale_y_continuous(breaks = c(50, 100)) +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

ts_fig_b <- ggplot(ts_dat_2, aes(x = GSYear, y = PercCovered,
                                    color = AreaCountyID)) +
  geom_line(size = 0.5) +
  facet_grid(rows = vars(AreaCountyID2)) +
  labs(x = "Year", title = "PAC") +
  scale_color_manual(values = rep(kelly()[-c(1, 22)], 2)) +
  scale_y_continuous(breaks = c(50, 100)) +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

# model fit fig
mf_fig_a <- ggplot(ts_dat_3, aes(x = PercTreated, y = PercDiffCovered,
                      color = AreaCountyID)) +
  geom_hline(yintercept = 0, size = 0.1, linetype = "dashed") +
  geom_smooth(method = "lm", formula = "y ~ x", size = 0.5) +
  geom_text(data = crash_coef4_a, 
            aes(x = maxTreat, y = maxDiff, label = p.round),
            hjust = 1, vjust = 1, size = 2, color = "black", fontface = "italic") +
  geom_text(data = crash_coef4_a, 
            aes(x = minTreat, y = maxDiff, label = est.round),
            hjust = 0, vjust = 1, size = 2, color = "black") +
  facet_grid(rows = vars(AreaCountyID2), labeller = label_parsed) +
  labs(x = "Percent area managed", title = "Annual difference in PAC") +
  scale_color_manual(values = rep(kelly()[-c(1, 22)], 2)) +
  scale_y_continuous(breaks = 0) +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text.y = element_text(size = 6, color="black",
                                    angle = 0, hjust = 0),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

mf_fig_b <- ggplot(ts_dat_4, aes(x = PercTreated, y = PercDiffCovered,
                                    color = AreaCountyID)) +
  geom_hline(yintercept = 0, size = 0.1, linetype = "dashed") +
  geom_smooth(method = "lm", formula = "y ~ x", size = 0.5) +
  geom_text(data = crash_coef4_b, 
            aes(x = maxTreat, y = maxDiff, label = p.round),
            hjust = 1, vjust = 1, size = 2, color = "black", fontface = "italic") +
  geom_text(data = crash_coef4_b, 
            aes(x = minTreat, y = maxDiff, label = est.round),
            hjust = 0, vjust = 1, size = 2, color = "black") +
  facet_grid(rows = vars(AreaCountyID2), labeller = label_parsed) +
  labs(x = "Percent area managed", title = "Annual difference in PAC") +
  scale_color_manual(values = rep(kelly()[-c(1, 22)], 2)) +
  scale_y_continuous(breaks = 0) +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text.y = element_text(size = 6, color="black",
                                    angle = 0, hjust = 0),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

# combine and save
comb_fig <- ts_fig_a + mf_fig_a + ts_fig_b + mf_fig_b +
  plot_layout(nrow = 1)

ggsave("output/fwc_hydrilla_crash_time_series_prediction.png", comb_fig,
       device = "png", width = 7, height = 9, units = "in")


#### values for text ####

crash_sum <- crash_coef3 %>%
  mutate(summ = case_when(estimate < 0 & p.value < 0.05 ~ "sig neg",
                          estimate > 0 & p.value < 0.05 ~ "sig pos",
                          estimate < 0 & p.value >= 0.05 ~ "n.s. neg",
                          estimate > 0 & p.value >= 0.05 ~ "n.s. pos")) %>%
  count(summ)

write_csv(crash_sum, "output/fwc_hydrilla_crash_prediction_summary.csv")
