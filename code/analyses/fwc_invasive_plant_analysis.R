#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(GGally)
library(fixest) # FE models
library(modelsummary)
library(cowplot)

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv",
                      col_types = list(PrevPropCovered = col_double(),
                                       PrevAreaCoveredRaw_ha = col_double(),
                                       SurveyDays = col_double(),
                                       RatioCovered = col_double(),
                                       LogRatioCovered = col_double(),
                                       LogitPrevPropCovered = col_double(),
                                       LogRatioCovered = col_double(),
                                       LogitPrevPropCovered = col_double(),
                                       LogRatioCovered = col_double()))
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")
lw_qual <- read_csv("intermediate-data/LW_quality_formatted.csv")
wa_qual <- read_csv("intermediate-data/water_atlas_quality_formatted.csv")


#### edit data ####

# combine datasets
inv_dat <- inv_plant %>%
  filter(!is.na(MinSurveyorExperience)) %>%
  inner_join(inv_ctrl) %>%
  filter(!is.na(PrevPropCovered)) %>% # need initial pop size
  mutate(PrevPercCovered = PrevPropCovered * 100)

# checked before removing NA's
filter(inv_dat, is.na(MinSurveyorExperience)) # 18 rows

# water quality
qual <- lw_qual %>%
  full_join(wa_qual) %>%
  filter(QualityMetric == "Secchi_ft") %>%
  group_by(PermanentID, GSYear) %>%
  summarize(QualityValue = mean(QualityValue),
            MonthsSampled = mean(MonthsSampled))

# add quality data
inv_qual_dat <- inv_dat %>%
  inner_join(qual %>%
               select(PermanentID, GSYear, QualityValue, MonthsSampled))

# split by species
hydr_dat <- filter(inv_dat, CommonName == "Hydrilla")
wale_dat <- filter(inv_dat, CommonName == "Water lettuce")
wahy_dat <- filter(inv_dat, CommonName == "Water hyacinth")
hydr_qual_dat <- filter(inv_qual_dat, CommonName == "Hydrilla")
wale_qual_dat <- filter(inv_qual_dat, CommonName == "Water lettuce")
wahy_qual_dat <- filter(inv_qual_dat, CommonName == "Water hyacinth")


#### initial visualizations ####

# all and species-specific control correlated?
cor.test(~ Lag0Treated + Lag0AllTreated, data = hydr_dat)
cor.test(~ Lag0Treated + Lag0AllTreated, data = wale_dat)
cor.test(~ Lag0Treated + Lag0AllTreated, data = wahy_dat)  
cor.test(~ Lag5Treated + Lag5AllTreated, data = hydr_dat)
cor.test(~ Lag5Treated + Lag5AllTreated, data = wale_dat)
cor.test(~ Lag5Treated + Lag5AllTreated, data = wahy_dat) 
# yes

# covariate correlations
hydr_dat %>%
  select(Lag0Treated, Lag1Treated, PrevPercCovered, MinSurveyorExperience) %>%
  ggpairs()

wahy_dat %>%
  select(Lag0Treated, Lag1Treated, PrevPercCovered, MinSurveyorExperience) %>%
  ggpairs()

wale_dat %>%
  select(Lag0Treated, Lag1Treated, PrevPercCovered, MinSurveyorExperience) %>%
  ggpairs()

# log ratio prop covered
ggplot(hydr_dat, aes(x = PrevPercCovered, y = LogRatioCovered, 
                     color = as.factor(Lag0Treated))) +
  geom_point()

ggplot(wale_dat, aes(x = PrevPercCovered, y = LogRatioCovered, 
                     color = as.factor(Lag0Treated))) +
  geom_point()

ggplot(wahy_dat, aes(x = PrevPercCovered, y = LogRatioCovered, 
                     color = as.factor(Lag0Treated))) +
  geom_point()

# treated and change in prop
ggplot(hydr_dat, aes(x = Lag5Treated, y = LogRatioCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

ggplot(wale_dat, aes(x = Lag5Treated, y = LogRatioCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

ggplot(wahy_dat, aes(x = Lag5Treated, y = LogRatioCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# data availability for lags
hydr_dat %>%
  select(PermanentID, GSYear,
         Lag0Treated, Lag1Treated, Lag2Treated, Lag3Treated, Lag4Treated, Lag5Treated) %>%
  pivot_longer(cols = starts_with("Lag"),
               names_to = "Lag",
               values_to = "Treated") %>%
  filter(!is.na(Treated)) %>%
  ggplot(aes(x = Lag)) +
  geom_bar()


#### model-fitting functions ####

# functions to fit models
mod_fit <- function(dat_in){
  
  # subset data
  dat_mod <- dat_in %>%
    filter(!is.na(Lag0Treated) & !is.na(Lag1Treated) & !is.na(Lag2Treated) & !is.na(Lag3Treated) & !is.na(Lag4Treated) & !is.na(Lag5Treated)) %>%
    mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
           PrevPercCovered_c = PrevPercCovered - mean(PrevPercCovered))
  
  # fit models
  mod0 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag0Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  mod1 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag1Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  mod2 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag2Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  mod3 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag3Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  mod4 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag4Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  mod5 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag5Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  
  # output
  return(list(mod0, mod1, mod2, mod3, mod4, mod5))

}

mod_qual_fit <- function(dat_in){
  
  # subset data
  dat_mod <- dat_in %>%
    filter(!is.na(Lag0Treated) & !is.na(Lag1Treated) & !is.na(Lag2Treated) & !is.na(Lag3Treated) & !is.na(Lag4Treated) & !is.na(Lag5Treated)) %>%
    mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
           PrevPercCovered_c = PrevPercCovered - mean(PrevPercCovered),
           Clarity_s = (QualityValue - mean(QualityValue)) / sd(QualityValue))
  
  # fit models
  mod0 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag0Treated + SurveyorExperience_s + Clarity_s | PermanentID + GSYear, data = dat_mod)
  mod1 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag1Treated + SurveyorExperience_s + Clarity_s| PermanentID + GSYear, data = dat_mod)
  mod2 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag2Treated + SurveyorExperience_s + Clarity_s| PermanentID + GSYear, data = dat_mod)
  mod3 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag3Treated + SurveyorExperience_s + Clarity_s| PermanentID + GSYear, data = dat_mod)
  mod4 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag4Treated + SurveyorExperience_s + Clarity_s| PermanentID + GSYear, data = dat_mod)
  mod5 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag5Treated + SurveyorExperience_s + Clarity_s| PermanentID + GSYear, data = dat_mod)
  
  # output
  return(list(mod0, mod1, mod2, mod3, mod4, mod5))
  
}

# function to compare models
mod_lag_comp <- function(mod_list){
  
  # compare models
  aic_mod <- tibble(Lag = c(0, 1, 2, 3, 4, 5),
                    AIC = AIC(mod_list[[1]], mod_list[[2]], mod_list[[3]],
                              mod_list[[4]], mod_list[[5]], mod_list[[6]])) %>%
    mutate(deltaAIC = AIC - min(AIC))
  
  # output
  return(aic_mod)
  
}


#### treatment, surveyor model figure ####

# fit models
hydr_mods <- mod_fit(hydr_dat)
wahy_mods <- mod_fit(wahy_dat)
wale_mods <- mod_fit(wale_dat)

# name models
names(hydr_mods) <- c("1", "2", "3", "4", "5", "6")
names(wahy_mods) <- c("1", "2", "3", "4", "5", "6")
names(wale_mods) <- c("1", "2", "3", "4", "5", "6")

# rename coefficients
coef_names <- c("SurveyorExperience_s" = "Surveyor experience",
                "PrevPercCovered_c:Lag0Treated" = "Treatment x abundance",
                "PrevPercCovered_c:Lag1Treated" = "Treatment x abundance",
                "PrevPercCovered_c:Lag2Treated" = "Treatment x abundance",
                "PrevPercCovered_c:Lag3Treated" = "Treatment x abundance",
                "PrevPercCovered_c:Lag4Treated" = "Treatment x abundance",
                "PrevPercCovered_c:Lag5Treated" = "Treatment x abundance",
                "PrevPercCovered_c" = "Initial abundance (%)",
                "Lag0Treated" = "Treatment frequency",
                "Lag1Treated" = "Treatment frequency",
                "Lag2Treated" = "Treatment frequency",
                "Lag3Treated" = "Treatment frequency",
                "Lag4Treated" = "Treatment frequency",
                "Lag5Treated" = "Treatment frequency")

# panels
hydr_fig <- modelplot(hydr_mods,
          coef_map = coef_names,
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1) +
  labs(x = "",
       title = "(A) hydrilla") +
  def_theme_paper +
  theme(legend.position = "none") +
  guides(color = guide_legend(reverse = TRUE))

wahy_fig <- modelplot(wahy_mods,
          coef_map = coef_names,
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       title = "(B) water hyacinth") +
  def_theme_paper +
  theme(legend.position = "none",
        axis.text.y = element_blank())

wale_fig <- modelplot(wale_mods,
          coef_map = coef_names,
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(name = "Years\nof\ndata", direction = -1) +
  labs(x = "",
       title = "(C) water lettuce") +
  def_theme_paper +
  theme(axis.text.y = element_blank(),
        legend.box.margin = margin(-10, 0, -10, -10))

# combine figures
pdf("output/fwc_invasive_plant_treatment_model.pdf", width = 6.5, height = 4)
plot_grid(hydr_fig, wahy_fig, wale_fig,
          nrow = 1,
          rel_widths = c(1, 0.57, 0.75))
dev.off()


#### treatment, surveyor, clarity model figure ####

# fit models
hydr_qual_mods <- mod_qual_fit(hydr_qual_dat)
wahy_qual_mods <- mod_qual_fit(wahy_qual_dat)
wale_qual_mods <- mod_qual_fit(wale_qual_dat)

# name models
names(hydr_qual_mods) <- c("1", "2", "3", "4", "5", "6")
names(wahy_qual_mods) <- c("1", "2", "3", "4", "5", "6")
names(wale_qual_mods) <- c("1", "2", "3", "4", "5", "6")

# rename coefficients
coef_qual_names <- c("Clarity_s" = "Water clarity", 
                     "SurveyorExperience_s" = "Surveyor experience",
                     "PrevPercCovered_c:Lag0Treated" = "Treatment x abundance",
                     "PrevPercCovered_c:Lag1Treated" = "Treatment x abundance",
                     "PrevPercCovered_c:Lag2Treated" = "Treatment x abundance",
                     "PrevPercCovered_c:Lag3Treated" = "Treatment x abundance",
                     "PrevPercCovered_c:Lag4Treated" = "Treatment x abundance",
                     "PrevPercCovered_c:Lag5Treated" = "Treatment x abundance",
                     "PrevPercCovered_c" = "Initial abundance (%)",
                     "Lag0Treated" = "Treatment frequency",
                     "Lag1Treated" = "Treatment frequency",
                     "Lag2Treated" = "Treatment frequency",
                     "Lag3Treated" = "Treatment frequency",
                     "Lag4Treated" = "Treatment frequency",
                     "Lag5Treated" = "Treatment frequency")

# panels
hydr_qual_fig <- modelplot(hydr_qual_mods,
                      coef_map = coef_qual_names,
                      background = list(geom_vline(xintercept = 0, color = "black",
                                                   size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       title = "(A) hydrilla") +
  def_theme_paper +
  theme(legend.position = "none") +
  guides(color = guide_legend(reverse = TRUE))

wahy_qual_fig <- modelplot(wahy_qual_mods,
                      coef_map = coef_qual_names,
                      background = list(geom_vline(xintercept = 0, color = "black",
                                                   size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       title = "(B) water hyacinth") +
  def_theme_paper +
  theme(legend.position = "none",
        axis.text.y = element_blank()) +
  guides(color = guide_legend(reverse = TRUE))

wale_qual_fig <- modelplot(wale_qual_mods,
                      coef_map = coef_qual_names,
                      background = list(geom_vline(xintercept = 0, color = "black",
                                                   size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(name = "Years\nof\ndata", direction = -1) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       title = "(C) water lettuce") +
  def_theme_paper +
  theme(axis.text.y = element_blank(),
        legend.box.margin = margin(-10, 0, -10, -10)) +
  guides(color = guide_legend(reverse = TRUE))

# combine figures
pdf("output/fwc_invasive_plant_treatment_clarity_model.pdf", width = 6.5, height = 4)
plot_grid(hydr_qual_fig, wahy_qual_fig, wale_qual_fig,
          nrow = 1,
          rel_widths = c(1, 0.57, 0.75))
dev.off()


#### model prediction figure ####

# function to format data
pred_fig <- function(dat_in, mod_list){
  
  # subset data
  raw_dat <- dat_in %>%
    filter(!is.na(Lag0Treated) & !is.na(Lag1Treated) & !is.na(Lag2Treated) & !is.na(Lag3Treated) & !is.na(Lag4Treated) & !is.na(Lag5Treated) & !is.na(SurveyorExperience)) %>%
    mutate(PrevPercCovered_c = PrevPercCovered - mean(PrevPercCovered))
  
  # prediction dataset
  pred_dat <- tibble(PrevPercCovered = c(0, 15, 30)) %>%
    mutate(PrevPercCovered_c = PrevPercCovered - mean(raw_dat$PrevPercCovered),
           Clarity_s = 0) %>%
    expand_grid(tibble(Lag5Treated = c(0, 1/6, 2/6, 1/2, 4/6, 5/6, 1))) %>%
    expand_grid(raw_dat %>%
                  select(PermanentID, GSYear) %>%
                  unique()) %>%
    mutate(SurveyorExperience_s = 0,
           PrevPercCovered_f = as.factor(PrevPercCovered) %>%
             fct_rev()) %>%
    mutate(Pred = predict(mod_list[[6]], newdata = .))
  
  fig_out <- ggplot(pred_dat, aes(x = Lag5Treated, y = Pred)) +
    geom_point(data = filter(raw_dat, PrevPercCovered <= 1), 
               alpha = 0.2, position = position_jitter(width = 0.01),
               shape = 21, color = "transparent",
               aes(y = LogRatioCovered, fill = PrevPercCovered)) +
    geom_point(data = filter(raw_dat, PrevPercCovered > 1), 
               alpha = 0.5, position = position_jitter(width = 0.01),
               shape = 21, color = "transparent",
               aes(y = LogRatioCovered, fill = PrevPercCovered)) +
    stat_summary(geom = "line", size = 1.3, fun = "median", aes(color = PrevPercCovered_f)) +
    scale_color_viridis_d(name = "Modeled\ninitial\nabundance (%)", begin = 0.7) +
    scale_fill_viridis_c(name = "Observed\ninitial\nabundance (%)", direction = -1) +
    labs(x = "", y = "Change in abundance") +
    def_theme_paper
  
  return(fig_out)
  
}

# figures
hydr_pred_fig <- pred_fig(hydr_dat, hydr_mods) +
  labs(title = "(A) hydrilla") +
  theme(legend.position = "none")
wahy_pred_fig <- pred_fig(wahy_dat, wahy_mods) +
  labs(title = "(B) water hyacinth", x = "Treatment frequency") +
  theme(legend.position = "none",
        axis.title.y = element_blank())
wale_pred_fig <- pred_fig(wale_dat, wale_mods) +
  labs(title = "(C) water lettuce") +
  theme(axis.title.y = element_blank())

hydr_qual_pred_fig <- pred_fig(hydr_qual_dat, hydr_qual_mods) +
  labs(title = "(A) hydrilla") +
  theme(legend.position = "none")
wahy_qual_pred_fig <- pred_fig(wahy_qual_dat, wahy_qual_mods) +
  labs(title = "(B) water hyacinth", x = "Treatment frequency") +
  theme(legend.position = "none",
        axis.title.y = element_blank())
wale_qual_pred_fig <- pred_fig(wale_qual_dat, wale_qual_mods) +
  labs(title = "(C) water lettuce") +
  theme(axis.title.y = element_blank())

# combine figures
pdf("output/fwc_invasive_plant_treatment_predictions.pdf", width = 6.5, height = 3)
plot_grid(hydr_pred_fig, wahy_pred_fig, wale_pred_fig,
          nrow = 1,
          rel_widths = c(1, 0.9, 1.25))
dev.off()

pdf("output/fwc_invasive_plant_treatment_clarity_predictions.pdf", width = 6.5, height = 3)
plot_grid(hydr_qual_pred_fig, wahy_qual_pred_fig, wale_qual_pred_fig,
          nrow = 1,
          rel_widths = c(1, 0.9, 1.25))
dev.off()



### model comparison ####

# compare models (use delta 4 as cut-off)
mod_lag_comp(hydr_mods) # 4, 5 similar
mod_lag_comp(wahy_mods) # all except 1 similar
mod_lag_comp(wale_mods) # 5 best

mod_lag_comp(hydr_qual_mods) # all except 1 similar
mod_lag_comp(wahy_qual_mods) # 3, 4, 5 similar
mod_lag_comp(wale_qual_mods) # 4, 5 similar

# output AIC values
aic_out <- tibble(Lag = mod_lag_comp(hydr_mods)$Lag,
                  Hydrilla = mod_lag_comp(hydr_mods)$deltaAIC,
                  Waterhyacinth = mod_lag_comp(wahy_mods)$deltaAIC,
                  Waterlettuce = mod_lag_comp(wale_mods)$deltaAIC)
write_csv(aic_out, "output/fwc_invasive_plant_delta_aic.csv")

aic_qual_out <- tibble(Lag = mod_lag_comp(hydr_qual_mods)$Lag,
                       Hydrilla = mod_lag_comp(hydr_qual_mods)$deltaAIC,
                       Waterhyacinth = mod_lag_comp(wahy_qual_mods)$deltaAIC,
                       Waterlettuce = mod_lag_comp(wale_qual_mods)$deltaAIC)
write_csv(aic_qual_out, "output/fwc_invasive_plant_clarity_delta_aic.csv")
