#### set-up ####

# this script is a condensed version of herbicides_analysis.R

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(magrittr) # abundance_dataset
library(lubridate) # abundance_dataset
library(glmmTMB)
library(lme4)
library(effects)

# library(GGally)
# library(lme4)
# library(ggeffects)

# import data
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
ctrl_old <- read_csv("intermediate-data/FWC_control_old_formatted.csv")
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
qual <- read_csv("intermediate-data/LW_quality_formatted.csv",
                 col_types = list(PermanentID = col_character(),
                                  JoinNotes = col_character()))

# 
# herb_type <- read_csv("intermediate-data/herbicide_types.csv")
# plant_detect <- read_csv("intermediate-data/FWC_plant_survey_first_detection_manual.csv")

# source functions
source("code/figure_settings.R")
source("code/beta_transformations.R")
source("code/abundance_duplicates.R")
source("code/cumulative_herbicides.R")
source("code/cumulative_abundance.R")
source("code/invasion_interval.R")
source("code/native_interval.R")
source("code/okeechobee_growth.R")
source("code/surveyor_experience.R")
source("code/plant_abundance_formatting.R")
source("code/herbicide_formatting.R")
source("code/abundance_herbicide_formatting.R")
source("code/growth_model.R")
source("code/cumulative_quality.R")
source("code/abundance_quality_formatting.R")
source("code/growth_quality_model.R")


#### invasive plant data formatting ####

inv_taxa <- tibble(TaxonName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"),
                   CommonName = c("Hydrilla", "Water lettuce", "Water hyacinth"))

inv_fwc2 <- abundance_dataset(plant_fwc, inv_taxa) # keep 2 for continuity with other scripts, change all at once

# remove missing data
inv_fwc3 <- inv_fwc2 %>%
  filter(!is.na(PropCovered))

# compare with previous version
inv_fwc3_check <- read_csv("intermediate-data/FWC_hydrilla_pistia_eichhornia_survey_formatted.csv")
all_equal(inv_fwc3 %>% mutate(SurveyorExperienceB = as.character(SurveyorExperienceB),
                              SurveyorExperience = if_else(is.na(SurveyorExperience), NA_real_, SurveyorExperience)), 
          inv_fwc3_check)


#### herbicide data formatting ####

herb_taxa <- tibble(Species = c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)"))

ctrl <- herbicide_dataset(ctrl_old, ctrl_new, herb_taxa)

# compare with previous version
ctrl_check <- read_csv("intermediate-data/FWC_hydrilla_pistia_eichhornia_herbicide_formatted.csv")
all_equal(ctrl, ctrl_check)


#### combine invasive plant and ctrl data ####

herb_inv_taxa <- tibble(Species = c("Hydrilla verticillata", rep("Floating Plants (Eichhornia and Pistia)", 2)), # double each row that has floating plants
                        TaxonName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"))

inv_ctrl <- herbicide_abundance_dataset(ctrl, inv_fwc3, herb_inv_taxa)

# compare with previous version
inv_ctrl_check <- read_csv("intermediate-data/FWC_hydrilla_pistia_eichhornia_survey_herbicide_formatted.csv") %>%
  rename(TaxonName = SpeciesName)
all_equal(inv_ctrl %>% mutate(SurveyorExperienceB = as.character(SurveyorExperienceB),
                              SurveyorExperience = if_else(is.na(SurveyorExperience), NA_real_, SurveyorExperience)), 
          inv_ctrl_check)


#### combine invasive plant and quality data ####

inv_qual <- quality_abundance_dataset(qual, inv_fwc3)


#### subset data ####

# subset for hydrilla
hydr_ctrl <- inv_ctrl %>% filter(CommonName == "Hydrilla")
hydr_qual <- inv_qual %>% filter(CommonName == "Hydrilla")


#### hydrilla full treatment models ####

# subset for all herbicide lags
hydr_ctrl2 <- hydr_ctrl %>%
  filter(!is.na(Lag0PropTreated) & 
           !is.na(Lag1PropTreated) & 
           !is.na(Lag2PropTreated) & 
           !is.na(Lag3PropTreated) & 
           !is.na(Lag4PropTreated) & 
           !is.na(Lag5PropTreated) & 
           !is.na(SurveyorExperienceB)) %>%
  mutate(PrevPropCoveredAdjCS = (PrevPropCoveredAdj - mean(PrevPropCoveredAdj)) / sd(PrevPropCoveredAdj))

# figures
hydr_ctrl2 %>% ggplot(aes(LogPropCovered)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(PrevPropCoveredAdjCS)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(SurveyorExperienceB)) + geom_bar() # should these be re-binned with each data subset?
hydr_ctrl2 %>% ggplot(aes(Lag0PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag1PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag2PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag3PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag4PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag5PropTreated)) + geom_histogram(binwidth = 0.1)

# random effects
length(unique(hydr_ctrl2$PermanentID))
length(unique(hydr_ctrl2$GSYear))

# tried beta models, but it was difficult to interpret the coefficients

# models
hydr_prop0_mod <- growth_model(hydr_ctrl2, Lag0PropTreated, "hydrilla_lag0_prop_treated")
hydr_prop1_mod <- growth_model(hydr_ctrl2, Lag1PropTreated, "hydrilla_lag1_prop_treated")
hydr_prop2_mod <- growth_model(hydr_ctrl2, Lag2PropTreated, "hydrilla_lag2_prop_treated")
hydr_prop3_mod <- growth_model(hydr_ctrl2, Lag3PropTreated, "hydrilla_lag3_prop_treated")
hydr_prop4_mod <- growth_model(hydr_ctrl2, Lag4PropTreated, "hydrilla_lag4_prop_treated")
hydr_prop5_mod <- growth_model(hydr_ctrl2, Lag5PropTreated, "hydrilla_lag5_prop_treated")
hydr_trtd0_mod <- growth_model(hydr_ctrl2, Lag0Treated, "hydrilla_lag0_bin_treated")
hydr_trtd1_mod <- growth_model(hydr_ctrl2, Lag1Treated, "hydrilla_lag1_bin_treated")
hydr_trtd2_mod <- growth_model(hydr_ctrl2, Lag2Treated, "hydrilla_lag2_bin_treated")
hydr_trtd3_mod <- growth_model(hydr_ctrl2, Lag3Treated, "hydrilla_lag3_bin_treated")
hydr_trtd4_mod <- growth_model(hydr_ctrl2, Lag4Treated, "hydrilla_lag4_bin_treated")
hydr_trtd5_mod <- growth_model(hydr_ctrl2, Lag5Treated, "hydrilla_lag5_bin_treated")

# compare models
AIC(hydr_prop0_mod, hydr_prop1_mod, hydr_prop2_mod, hydr_prop3_mod, hydr_prop4_mod, hydr_prop5_mod,
    hydr_trtd0_mod, hydr_trtd1_mod, hydr_trtd2_mod, hydr_trtd3_mod, hydr_trtd4_mod, hydr_trtd5_mod) %>%
  arrange(AIC) %>%
  mutate(deltaAIC = AIC - min(AIC))

# best model
summary(hydr_prop0_mod)

# use largest possible dataset
hydr_ctrl3 <- hydr_ctrl %>%
  filter(!is.na(Lag0PropTreated) & !is.na(SurveyorExperienceB)) %>%
  mutate(PrevPropCoveredAdjCS = (PrevPropCoveredAdj - mean(PrevPropCoveredAdj)) / sd(PrevPropCoveredAdj))

# figures
hydr_ctrl3 %>% ggplot(aes(LogPropCovered)) + geom_histogram(binwidth = 0.1)
hydr_ctrl3 %>% ggplot(aes(PrevPropCoveredAdjCS)) + geom_histogram(binwidth = 0.1)
hydr_ctrl3 %>% ggplot(aes(SurveyorExperience)) + geom_histogram(binwidth = 10)
hydr_ctrl3 %>% ggplot(aes(SurveyorExperienceB)) + geom_bar()
hydr_ctrl3 %>% ggplot(aes(Lag0PropTreated)) + geom_histogram(binwidth = 0.1)

# refit model
hydr_mod <- growth_model(hydr_ctrl3, Lag0PropTreated, "hydrilla_treated")
summary(hydr_mod)


#### hydrilla new treatment models ####

# earliest year for new dataset
min_ctrl_new <- min(filter(ctrl, CtrlSet == "new")$GSYear)

# subset for all herbicide lags
hydr_ctrl_new <- hydr_ctrl %>%
  filter(!is.na(Lag0PropTreated) & 
           !is.na(Lag1PropTreated) & 
           !is.na(Lag2PropTreated) & 
           !is.na(Lag3PropTreated) & 
           !is.na(Lag4PropTreated) & 
           !is.na(Lag5PropTreated) & 
           !is.na(SurveyorExperienceB) &
           GSYear > min_ctrl_new + 5) %>%
  mutate(PrevPropCoveredAdjCS = (PrevPropCoveredAdj - mean(PrevPropCoveredAdj)) / sd(PrevPropCoveredAdj))

# figures
hydr_ctrl_new %>% ggplot(aes(LogPropCovered)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new %>% ggplot(aes(PrevPropCoveredAdjCS)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new %>% ggplot(aes(SurveyorExperienceB)) + geom_bar() # should these be re-binned with each data subset?
hydr_ctrl_new %>% ggplot(aes(Lag0PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new %>% ggplot(aes(Lag1PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new %>% ggplot(aes(Lag2PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new %>% ggplot(aes(Lag3PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new %>% ggplot(aes(Lag4PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new %>% ggplot(aes(Lag5PropTreated)) + geom_histogram(binwidth = 0.1)

# random effects
length(unique(hydr_ctrl_new$PermanentID))
length(unique(hydr_ctrl_new$GSYear))

# tried beta models, but it was difficult to interpret the coefficients

# models
hydr_new_prop0_mod <- growth_model(hydr_ctrl_new, Lag0PropTreated, "hydrilla_new_herb_lag0_prop_treated")
hydr_new_prop1_mod <- growth_model(hydr_ctrl_new, Lag1PropTreated, "hydrilla_new_herb_lag1_prop_treated")
hydr_new_prop2_mod <- growth_model(hydr_ctrl_new, Lag2PropTreated, "hydrilla_new_herb_lag2_prop_treated")
hydr_new_prop3_mod <- growth_model(hydr_ctrl_new, Lag3PropTreated, "hydrilla_new_herb_lag3_prop_treated")
hydr_new_prop4_mod <- growth_model(hydr_ctrl_new, Lag4PropTreated, "hydrilla_new_herb_lag4_prop_treated")
hydr_new_prop5_mod <- growth_model(hydr_ctrl_new, Lag5PropTreated, "hydrilla_new_herb_lag5_prop_treated")
hydr_new_trtd0_mod <- growth_model(hydr_ctrl_new, Lag0Treated, "hydrilla_new_herb_lag0_bin_treated")
hydr_new_trtd1_mod <- growth_model(hydr_ctrl_new, Lag1Treated, "hydrilla_new_herb_lag1_bin_treated")
hydr_new_trtd2_mod <- growth_model(hydr_ctrl_new, Lag2Treated, "hydrilla_new_herb_lag2_bin_treated")
hydr_new_trtd3_mod <- growth_model(hydr_ctrl_new, Lag3Treated, "hydrilla_new_herb_lag3_bin_treated")
hydr_new_trtd4_mod <- growth_model(hydr_ctrl_new, Lag4Treated, "hydrilla_new_herb_lag4_bin_treated")
hydr_new_trtd5_mod <- growth_model(hydr_ctrl_new, Lag5Treated, "hydrilla_new_herb_lag5_bin_treated")

# compare models
AIC(hydr_new_prop0_mod, hydr_new_prop1_mod, hydr_new_prop2_mod, hydr_new_prop3_mod, hydr_new_prop4_mod, hydr_new_prop5_mod,
    hydr_new_trtd0_mod, hydr_new_trtd1_mod, hydr_new_trtd2_mod, hydr_new_trtd3_mod, hydr_new_trtd4_mod, hydr_new_trtd5_mod) %>%
  arrange(AIC) %>%
  mutate(deltaAIC = AIC - min(AIC))

# best models
summary(hydr_new_prop1_mod)
summary(hydr_new_prop0_mod)

# use largest possible dataset
hydr_ctrl_new2 <- hydr_ctrl %>%
  filter(!is.na(Lag0PropTreated) & 
           !is.na(SurveyorExperienceB) &
           GSYear > min_ctrl_new) %>%
  mutate(PrevPropCoveredAdjCS = (PrevPropCoveredAdj - mean(PrevPropCoveredAdj)) / sd(PrevPropCoveredAdj))

# figures
hydr_ctrl_new2 %>% ggplot(aes(LogPropCovered)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new2 %>% ggplot(aes(PrevPropCoveredAdjCS)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new2 %>% ggplot(aes(SurveyorExperienceB)) + geom_bar()
hydr_ctrl_new2 %>% ggplot(aes(Lag0PropTreated)) + geom_histogram(binwidth = 0.1)

# refit model
hydr_new_mod <- growth_model(hydr_ctrl_new2, Lag0PropTreated, "hydrilla_new_treated")
summary(hydr_new_mod)


#### hydrilla quality models ####

# subset for all quality lags
hydr_qual2 <- hydr_qual %>%
  select(c("PermanentID", "GSYear", "LogPropCovered", "PrevPropCoveredAdj"), starts_with("Lag")) %>%
  drop_na() %>%
  mutate(PrevPropCoveredAdjCS = (PrevPropCoveredAdj - mean(PrevPropCoveredAdj)) / sd(PrevPropCoveredAdj))

#### start here ####

# figures
hydr_qual2 %>% ggplot(aes(LogPropCovered)) + geom_histogram(binwidth = 0.1)
hydr_qual2 %>% ggplot(aes(PrevPropCoveredAdjCS)) + geom_histogram(binwidth = 0.1)
hydr_qual2 %>% ggplot(aes(Lag0PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_qual2 %>% ggplot(aes(Lag1PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_qual2 %>% ggplot(aes(Lag2PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_qual2 %>% ggplot(aes(Lag3PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_qual2 %>% ggplot(aes(Lag4PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_qual2 %>% ggplot(aes(Lag5PropTreated)) + geom_histogram(binwidth = 0.1)

# random effects
length(unique(hydr_ctrl_new$PermanentID))
length(unique(hydr_ctrl_new$GSYear))

# tried beta models, but it was difficult to interpret the coefficients

# models
hydr_new_prop0_mod <- growth_model(hydr_ctrl_new, Lag0PropTreated, "hydrilla_new_herb_lag0_prop_treated")
hydr_new_prop1_mod <- growth_model(hydr_ctrl_new, Lag1PropTreated, "hydrilla_new_herb_lag1_prop_treated")
hydr_new_prop2_mod <- growth_model(hydr_ctrl_new, Lag2PropTreated, "hydrilla_new_herb_lag2_prop_treated")
hydr_new_prop3_mod <- growth_model(hydr_ctrl_new, Lag3PropTreated, "hydrilla_new_herb_lag3_prop_treated")
hydr_new_prop4_mod <- growth_model(hydr_ctrl_new, Lag4PropTreated, "hydrilla_new_herb_lag4_prop_treated")
hydr_new_prop5_mod <- growth_model(hydr_ctrl_new, Lag5PropTreated, "hydrilla_new_herb_lag5_prop_treated")
hydr_new_trtd0_mod <- growth_model(hydr_ctrl_new, Lag0Treated, "hydrilla_new_herb_lag0_bin_treated")
hydr_new_trtd1_mod <- growth_model(hydr_ctrl_new, Lag1Treated, "hydrilla_new_herb_lag1_bin_treated")
hydr_new_trtd2_mod <- growth_model(hydr_ctrl_new, Lag2Treated, "hydrilla_new_herb_lag2_bin_treated")
hydr_new_trtd3_mod <- growth_model(hydr_ctrl_new, Lag3Treated, "hydrilla_new_herb_lag3_bin_treated")
hydr_new_trtd4_mod <- growth_model(hydr_ctrl_new, Lag4Treated, "hydrilla_new_herb_lag4_bin_treated")
hydr_new_trtd5_mod <- growth_model(hydr_ctrl_new, Lag5Treated, "hydrilla_new_herb_lag5_bin_treated")

# compare models
AIC(hydr_new_prop0_mod, hydr_new_prop1_mod, hydr_new_prop2_mod, hydr_new_prop3_mod, hydr_new_prop4_mod, hydr_new_prop5_mod,
    hydr_new_trtd0_mod, hydr_new_trtd1_mod, hydr_new_trtd2_mod, hydr_new_trtd3_mod, hydr_new_trtd4_mod, hydr_new_trtd5_mod) %>%
  arrange(AIC) %>%
  mutate(deltaAIC = AIC - min(AIC))

# best models
summary(hydr_new_prop1_mod)
summary(hydr_new_prop0_mod)

# use largest possible dataset
hydr_ctrl_new2 <- hydr_ctrl %>%
  filter(!is.na(Lag0PropTreated) & 
           !is.na(SurveyorExperienceB) &
           GSYear > min_ctrl_new) %>%
  mutate(PrevPropCoveredAdjCS = (PrevPropCoveredAdj - mean(PrevPropCoveredAdj)) / sd(PrevPropCoveredAdj))

# figures
hydr_ctrl_new2 %>% ggplot(aes(LogPropCovered)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new2 %>% ggplot(aes(PrevPropCoveredAdjCS)) + geom_histogram(binwidth = 0.1)
hydr_ctrl_new2 %>% ggplot(aes(SurveyorExperienceB)) + geom_bar()
hydr_ctrl_new2 %>% ggplot(aes(Lag0PropTreated)) + geom_histogram(binwidth = 0.1)

# refit model
hydr_new_mod <- growth_model(hydr_ctrl_new2, Lag0PropTreated, "hydrilla_new_treated")
summary(hydr_new_mod)

#### add in native species code from herbicides_analysis ####

