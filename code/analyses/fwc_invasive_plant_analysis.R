#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(GGally)
library(fixest) # FE models
library(modelsummary)

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
qual <- read_csv("intermediate-data/LW_quality_formatted.csv")


#### edit data ####

# combine datasets
inv_dat <- inv_plant %>%
  inner_join(inv_ctrl) %>%
  filter(!is.na(PrevPropCovered)) %>% # need initial pop size
  mutate(PrevPercCovered = PrevPropCovered * 100)

# add quality data
inv_qual_dat <- inv_dat %>%
  inner_join(qual) %>%
  filter(!is.na(Secchi))

filter(inv_dat, is.na(SurveyorExperience))
# only three entries are missing surveyor experience

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
  select(Lag0Treated, Lag1Treated, PrevPercCovered, SurveyorExperience) %>%
  ggpairs()

wahy_dat %>%
  select(Lag0Treated, Lag1Treated, PrevPercCovered, SurveyorExperience) %>%
  ggpairs()

wale_dat %>%
  select(Lag0Treated, Lag1Treated, PrevPercCovered, SurveyorExperience) %>%
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


#### fit models ####

# function to compare models
mod_lag_comp <- function(dat_in){
  
  # subset data
  dat_mod <- dat_in %>%
    filter(!is.na(Lag0Treated) & !is.na(Lag1Treated) & !is.na(Lag2Treated) & !is.na(Lag3Treated) & !is.na(Lag4Treated) & !is.na(Lag5Treated) & !is.na(SurveyorExperience)) %>%
    mutate(SurveyorExperience_s = (SurveyorExperience - mean(SurveyorExperience)) / sd(SurveyorExperience),
           PrevPercCovered_c = PrevPercCovered - mean(PrevPercCovered))
  
  # fit models
  mod0 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag0Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  mod1 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag1Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  mod2 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag2Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  mod3 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag3Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  mod4 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag4Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  mod5 <- feols(LogRatioCovered ~ PrevPercCovered_c * Lag5Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod)
  
  # compare models
  aic_mod <- tibble(Lag = c(0, 1, 2, 3, 4, 5),
                    AIC = AIC(mod0, mod1, mod2, mod3, mod4, mod5)) %>%
    mutate(deltaAIC = AIC - min(AIC))
  
  # output
  return(aic_mod)
  
}

# compare models
mod_lag_comp(hydr_dat) # 0, 4, 5 similar
mod_lag_comp(wahy_dat) # 0, 4, 5 similar
mod_lag_comp(wale_dat) # 5 best
# use lag 5

# output AIC values
aic_out <- tibble(Lag = mod_lag_comp(hydr_dat)$Lag,
                  Hydrilla = mod_lag_comp(hydr_dat)$deltaAIC,
                  Waterhyacinth = mod_lag_comp(wahy_dat)$deltaAIC,
                  Waterlettuce = mod_lag_comp(wale_dat)$deltaAIC)
write_csv(aic_out, "output/fwc_invasive_plant_delta_aic.csv")

# function to fit models
mod_lag_fit <- function(dat_in, dat_qual_in){
  
  # subset data
  dat_mod <- dat_in %>%
    filter(!is.na(Lag5Treated) & !is.na(SurveyorExperience)) %>%
    mutate(SurveyorExperience_s = (SurveyorExperience - mean(SurveyorExperience)) / sd(SurveyorExperience),
           PrevPercCovered_c = PrevPercCovered - mean(PrevPercCovered))
  
  dat_qual_mod <- dat_qual_in %>%
    filter(!is.na(Lag5Treated) & !is.na(SurveyorExperience)) %>%
    mutate(SurveyorExperience_s = (SurveyorExperience - mean(SurveyorExperience)) / sd(SurveyorExperience),
           PrevPercCovered_c = PrevPercCovered - mean(PrevPercCovered),
           Turbidity_s = (Secchi - mean(Secchi)) / sd(Secchi))
  
  # fit models
  mod <- feols(LogRatioCovered ~ Lag5Treated * PrevPercCovered_c + SurveyorExperience_s  | PermanentID + GSYear, data = dat_mod)
  mod_qual <- feols(LogRatioCovered  ~ Lag5Treated * PrevPercCovered_c + SurveyorExperience_s + Turbidity_s | PermanentID + GSYear, data = dat_qual_mod)
  mod_simp <- feols(LogRatioCovered ~ Lag5Treated * PrevPercCovered_c  | PermanentID + GSYear, data = dat_mod)
  
  # output
  return(list(mod, mod_qual, mod_simp))
  
}

# fit models
hydr_mod <- mod_lag_fit(hydr_dat, hydr_qual_dat)[[1]]
hydr_qual_mod <- mod_lag_fit(hydr_dat, hydr_qual_dat)[[2]]
wahy_mod <- mod_lag_fit(wahy_dat, wahy_qual_dat)[[1]]
wahy_qual_mod <- mod_lag_fit(wahy_dat, wahy_qual_dat)[[2]]
wale_mod <- mod_lag_fit(wale_dat, wale_qual_dat)[[1]]
wale_qual_mod <- mod_lag_fit(wale_dat, wale_qual_dat)[[2]]

# model summaries
summary(hydr_mod)
summary(wahy_mod)
summary(wale_mod)
# surveyor experience not significant for any

summary(hydr_qual_mod)
summary(wahy_qual_mod)
summary(wale_qual_mod)
# turbidity increases hydrilla growth (?)

# simplify models (remove surveyor experience)
hydr_simp_mod <- mod_lag_fit(hydr_dat, hydr_qual_dat)[[3]]
wahy_simp_mod <- mod_lag_fit(wahy_dat, wahy_qual_dat)[[3]]
wale_simp_mod <- mod_lag_fit(wale_dat, wale_qual_dat)[[3]]

summary(hydr_simp_mod)
summary(wahy_simp_mod)
summary(wale_simp_mod)


#### figures ####

# compile models
fig_mods <- list()
fig_mods[['water lettuce']] <- wale_simp_mod
fig_mods[['water hyacinth']] <- wahy_simp_mod
fig_mods[['hydrilla']] <- hydr_simp_mod

pdf("output/fwc_invasive_plant_treatment_model.pdf", width = 5, height = 5)
modelplot(fig_mods,
          coef_map = c("Lag5Treated:PrevPercCovered_c" = "Treatment x abundance",
                       "PrevPercCovered_c" = "Initial abundance",
                       "Lag5Treated" = "6-year treatment\nfrequency"),
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed"))) +
  scale_color_brewer(type = "qual", palette = "Dark2", direction = -1) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = ""))) +
  def_theme_paper +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.1)) +
  guides(color = guide_legend(reverse = TRUE))
dev.off()
