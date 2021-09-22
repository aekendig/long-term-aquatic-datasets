#### set-up ####

# this script is a condensed version of herbicides_analysis.R
# left off condensing sections and fitting models with water quality

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(magrittr) # abundance_dataset
library(lubridate) # abundance_dataset
library(glmmTMB)
library(effects)
library(car) # for logit
library(fixest) # FE models
# library(lfe) # FE models
# library(alpaca) # FE models
# library(lme4) # use glmmTMB unless it's too slow

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
source("code/proportion_transformations.R")
source("code/abundance_duplicates.R")
source("code/cumulative_herbicides.R")
source("code/cumulative_abundance.R")
source("code/invasion_interval.R")
source("code/native_interval.R")
source("code/okeechobee_growth.R")
source("code/surveyor_experience.R")
source("code/plant_abundance_formatting.R")
source("code/herbicide_formatting.R")
source("code/herbicide_old_formatting.R")
source("code/herbicide_new_formatting.R")
source("code/abundance_herbicide_formatting.R")
source("code/growth_model.R")
source("code/cumulative_quality.R")
source("code/abundance_quality_formatting.R")
source("code/growth_quality_model.R")

# prevent scientific notation
options(scipen = 999)


#### invasive plant data formatting ####

inv_taxa <- tibble(TaxonName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"),
                   CommonName = c("Hydrilla", "Water lettuce", "Water hyacinth"))

inv_fwc2 <- abundance_dataset(plant_fwc, inv_taxa) # keep 2 for continuity with other scripts, change all at once

# remove missing data
inv_fwc3 <- inv_fwc2 %>%
  filter(!is.na(EstAreaCoveredRaw_ha)) # missing surveys for year-lake combos

# compare with previous version
# inv_fwc3_check <- read_csv("intermediate-data/FWC_hydrilla_pistia_eichhornia_survey_formatted.csv")
# all_equal(inv_fwc3 %>% mutate(SurveyorExperienceB = as.character(SurveyorExperienceB),
#                               SurveyorExperience = if_else(is.na(SurveyorExperience), NA_real_, SurveyorExperience)), 
#           inv_fwc3_check)


#### herbicide data formatting ####

herb_taxa <- tibble(Species = c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)"))

# ctrl <- herbicide_dataset(ctrl_old, ctrl_new, herb_taxa) # new and old control datasets combined
herb_old <- herbicide_old_dataset(ctrl_old, herb_taxa) # note that this includes herbicide and non-herbicide methods
herb_new <- herbicide_new_dataset(ctrl_new, herb_taxa)

# overlap
herb_old %>%
  filter(AreaTreated_ha > 0) %>%
  select(PermanentID, Species, GSYear) %>%
  unique() %>%
  inner_join(herb_new %>%
               filter(AreaTreated_ha > 0) %>%
               select(PermanentID, Species, GSYear) %>%
               unique())
# yes, 2010 is repeated in both, 154 cases

# remove overlapping data
herb_new2 <- herb_new %>%
  anti_join(herb_old %>%
              filter(AreaTreated_ha > 0) %>%
              select(PermanentID, Species, GSYear) %>%
              unique()) # unclear what the total herbicide amount was for this growing season

# compare with previous version
# ctrl_check <- read_csv("intermediate-data/FWC_hydrilla_pistia_eichhornia_herbicide_formatted.csv")
# all_equal(ctrl, ctrl_check)


#### combine invasive plant and ctrl data ####

herb_inv_taxa <- tibble(Species = c("Hydrilla verticillata", rep("Floating Plants (Eichhornia and Pistia)", 2)), # double each row that has floating plants
                        TaxonName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"))

# inv_ctrl <- herbicide_abundance_dataset(ctrl, inv_fwc3, herb_inv_taxa)
inv_herb_new <- herbicide_abundance_dataset(herb_new2, inv_fwc3, herb_inv_taxa)

# compare with previous version
# inv_ctrl_check <- read_csv("intermediate-data/FWC_hydrilla_pistia_eichhornia_survey_herbicide_formatted.csv") %>%
#   rename(TaxonName = SpeciesName)
# all_equal(inv_ctrl %>% mutate(SurveyorExperienceB = as.character(SurveyorExperienceB),
#                               SurveyorExperience = if_else(is.na(SurveyorExperience), NA_real_, SurveyorExperience)), 
#           inv_ctrl_check)


#### combine invasive plant and quality data ####

# inv_qual <- quality_abundance_dataset(qual, inv_fwc3)


#### how does invasion affect treatment? ####

# subset for lakes with invasions
inv_treat_new <- inv_herb_new %>%
  filter(!is.na(Lag0PropTreated) & PrevSpeciesPresent == 1) %>%
  mutate(Treated = fct_recode(as.character(Lag0Treated), "Untreated" = "0", "Treated" = "1"),
         AreaTreatedAdj_ha = case_when(Lag0AreaTreated_ha > Area_ha ~ Area_ha, # if area treated in a year exceeds lake size, use lake size
                                       TRUE ~ Lag0AreaTreated_ha),
         AreaTreatedRounded_ha = round(AreaTreatedAdj_ha), # round for glm
         AreaUntreatedRounded_ha = round(Area_ha) - AreaTreatedRounded_ha,
         PropTreatedAdj = AreaTreatedAdj_ha/Area_ha,
         Log10Area_ha = log10(Area_ha),
         LogitPropTreated = logit(PropTreatedAdj, adjust = prop_adjust),
         Log10PrevAreaCovered_ha = log10(PrevAreaCoveredRaw_ha),
         LogitPrevPropCovered = logit(PrevPropCovered, adjust = prop_adjust),
         PrevAreaCoveredCS_ha = (PrevAreaCoveredRaw_ha - mean(PrevAreaCoveredRaw_ha)) / sd(PrevAreaCoveredRaw_ha),
         PrevAreaCoveredRounded_ha = case_when(PrevAreaCoveredRaw_ha > Area_ha ~ round(Area_ha), # if area covered exceed lake size, use lake size
                                               TRUE ~ round(PrevAreaCoveredRaw_ha)), # round for glm
         PrevAreaUncoveredRounded_ha = round(Area_ha) - PrevAreaCoveredRounded_ha)

# divide lake areas into groups
inv_treat_area <- inv_treat_new %>%
  select(PermanentID, Area_ha) %>%
  unique() %>%
  mutate(AreaGroup = cut_number(Area_ha, n = 3))

# rename lake groups
levels(inv_treat_area$AreaGroup) <- c("small", "medium", "large")

inv_treat_new2 <- inv_treat_new %>%
  left_join(inv_treat_area)

inv_treat_new_sum <- inv_treat_new2 %>%
  group_by(Treated) %>%
  summarise(Log10PrevAreaCovered_ha = mean(Log10PrevAreaCovered_ha),
            LogitPrevPropCovered = mean(LogitPrevPropCovered)) %>%
  ungroup()

# t-tests
inv_treat_mod1 <- t.test(Log10PrevAreaCovered_ha ~ Treated, data = inv_treat_new2)
inv_treat_mod1
# treated invasions have larger total area than untreated invasions

inv_treat_mod2 <- t.test(LogitPrevPropCovered ~ Treated, data = inv_treat_new2)
inv_treat_mod2
# treated invasions have larger proportion of lake covered than untreated invasions
# magnitude of difference seems small

# histogram to illustrate
ggplot(inv_treat_new2, aes(x = Log10PrevAreaCovered_ha)) +
  geom_histogram() +
  geom_vline(data = inv_treat_new_sum, aes(xintercept = Log10PrevAreaCovered_ha)) +
  facet_wrap(~ Treated) +
  def_theme_facet

ggplot(inv_treat_new2, aes(x = LogitPrevPropCovered)) +
  geom_histogram() +
  geom_vline(data = inv_treat_new_sum, aes(xintercept = LogitPrevPropCovered)) +
  facet_wrap(~ Treated) +
  def_theme_facet

# check proportion conversion
ggplot(inv_treat_new2, aes(x = PropTreatedAdj, y = LogitPropTreated)) +
  geom_point()

# invasion/treatment relationships
ggplot(inv_treat_new2, aes(x = PrevPropCovered, y = PropTreatedAdj)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new2, aes(x = Log10PrevAreaCovered_ha, y = PropTreatedAdj)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new2, aes(x = PrevAreaCoveredCS_ha, y = PropTreatedAdj)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new2, aes(x = PrevPropCovered, y = LogitPropTreated)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new2, aes(x = Log10PrevAreaCovered_ha, y = LogitPropTreated)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new2, aes(x = LogitPrevPropCovered, y = LogitPropTreated)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

# area models
inv_area_mod1 <- glmmTMB(Log10PrevAreaCovered_ha ~ Log10Area_ha + (1|PermanentID) + (1|GSYear), data = inv_treat_new2)
summary(inv_area_mod1)
# larger lakes have larger invasions

inv_area_mod2 <- glmmTMB(cbind(PrevAreaCoveredRounded_ha, PrevAreaUncoveredRounded_ha) ~ Log10Area_ha + (1|PermanentID) + (1|GSYear), data = inv_treat_new2, family = "binomial")
summary(inv_area_mod2)
# larger lakes have smaller proportions covered

# mixed-effects herbicide models
inv_herb_mod1 <- glmmTMB(cbind(AreaTreatedRounded_ha, AreaUntreatedRounded_ha) ~ Log10PrevAreaCovered_ha * CommonName + Log10Area_ha + (1|PermanentID) + (1|GSYear), data = inv_treat_new2, family= "binomial")
summary(inv_herb_mod1)
# larger invasions have more proportion (lake) treated
# larger lakes have less proportion treated

inv_herb_mod2 <- glmmTMB(cbind(AreaTreatedRounded_ha, AreaUntreatedRounded_ha) ~ PrevPropCovered * CommonName + Log10Area_ha + (1|PermanentID) + (1|GSYear), data = inv_treat_new2, family= "binomial")
summary(inv_herb_mod2)
# larger invasions have more proportion (lake) treated
# larger lakes have less proportion treated

inv_herb_mod3 <- glmmTMB(cbind(AreaTreatedRounded_ha, AreaUntreatedRounded_ha) ~ LogitPrevPropCovered * CommonName + Log10Area_ha + (1|PermanentID) + (1|GSYear), data = inv_treat_new2, family= "binomial")
summary(inv_herb_mod3)

# fixed-effects herbicide models
# no fe methods available for cbind(success, failure) or beta distributions
inv_herb_mod4 <- feols(LogitPropTreated ~ Log10PrevAreaCovered_ha * CommonName | PermanentID + GSYear, data = inv_treat_new2)
summary(inv_herb_mod4)
# Log10Area_ha and AreaGroup are collinear with PermanentID
# invasion size increases treatment

inv_herb_mod5 <- feols(LogitPropTreated ~ PrevPropCovered * CommonName | PermanentID + GSYear, data = inv_treat_new2)
summary(inv_herb_mod5)
# invasion size increases treatment

inv_herb_mod6 <- feols(LogitPropTreated ~ LogitPrevPropCovered * CommonName | PermanentID + GSYear, data = inv_treat_new2)
summary(inv_herb_mod6)

# predicted data
inv_treat_pred <- inv_treat_new2 %>%
  select(CommonName, PrevPropCovered, LogitPrevPropCovered, Log10PrevAreaCovered_ha, PrevAreaCoveredRaw_ha) %>%
  unique() %>%
  mutate(PermanentID = inv_treat_new2 %>%
           group_by(PermanentID) %>%
           count() %>%
           arrange(desc(n)) %>%
           select(PermanentID) %>%
           head(n = 1) %>% 
           pull(), # select most sampled lake
         GSYear = "2018") %>%
  mutate(pred4 = predict(inv_herb_mod4, newdata = .),
         pred6 = predict(inv_herb_mod6, newdata = .),
         pred4_prop = inv_logit_adjust(pred4, a = prop_adjust),
         pred6_prop = inv_logit_adjust(pred6, a = prop_adjust))

# check proportion conversion
ggplot(inv_treat_pred, aes(x = pred4_prop, y = pred4)) +
  geom_point()

# change in prop for each species
inv_treat_change <- inv_treat_pred %>%
  group_by(CommonName) %>%
  summarise(Log10PrevAreaCovered_ha = max(Log10PrevAreaCovered_ha),
            PrevAreaCoveredRaw_ha = max(PrevAreaCoveredRaw_ha),
            LogitPrevPropCovered = max(LogitPrevPropCovered),
            MinProp4 = min(pred4_prop),
            MaxProp4 = max(pred4_prop),
            MinProp6 = min(pred6_prop),
            MaxProp6 = max(pred6_prop)) %>%
  ungroup() %>%
  mutate(ChangeProp4 = paste("Delta == ", round(MaxProp4 - MinProp4, 2), sep = ""),
         ChangeProp6 = paste("Delta == ", round(MaxProp6 - MinProp6, 2), sep = ""))

# figure
pdf("output/new_herbicide_invasion_area.pdf", width = 11, height = 4)
inv_treat_pred %>%
  ggplot(aes(x = PrevAreaCoveredRaw_ha, y = pred4)) +
  geom_point(data = inv_treat_new2, alpha = 0.3, aes(y = LogitPropTreated)) +
  geom_line(color = "#009193", size = 1.5) +
  geom_text(data = inv_treat_change, y = 6.5, aes(label = ChangeProp4), 
            parse = T, check_overlap = T, hjust = 1, vjust = 1, size = 4) +
  facet_wrap(~ CommonName, scales = "free_x") +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000),
                labels = c(0.1, 1, 10, 100, 1000, 10000)) +
  xlab("Area covered by species (ha)") +
  ylab("Proportion of waterbody treated (log-odds)") +
  def_theme_facet +
  theme(axis.text.x = element_text(size = 10, color="black"))
dev.off()

pdf("output/new_herbicide_invasion_proportion.pdf", width = 11, height = 4)
inv_treat_pred %>%
  ggplot(aes(x = LogitPrevPropCovered, y = pred6)) +
  geom_point(data = inv_treat_new2, alpha = 0.3, aes(y = LogitPropTreated)) +
  geom_line(color = "#009193", size = 1.5) +
  geom_text(data = inv_treat_change, y = 6.5, aes(label = ChangeProp6), 
            parse = T, check_overlap = T, hjust = 1, vjust = 1, size = 4,
            nudge_x = -0.7) +
  facet_wrap(~ CommonName, scales = "free_x") +
  xlab("Proportion of waterbody covered by species (log-odds)") +
  ylab("Proportion of waterbody treated (log-odds)") +
  def_theme_facet
dev.off()


#### start here ####
#### how does treatment affect invasion? ####

# treatment variables
ggplot(inv_treat_new2, aes(x = LogitPropTreated, y = Lag0TreatmentDays)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)

# fit model
inv_chg_mod1 <- feols(LogRatioCovered ~ (LogitPropTreated + Log10PrevAreaCovered_ha) * CommonName | PermanentID + GSYear, data = inv_treat_new2)
summary(inv_chg_mod1)

# predicted data
inv_chg_pred <- inv_treat_new2 %>%
  select(CommonName, LogitPropTreated, Log10PrevAreaCovered_ha) %>%
  unique() %>%
  mutate(PermanentID = inv_treat_new2 %>%
           group_by(PermanentID) %>%
           count() %>%
           arrange(desc(n)) %>%
           select(PermanentID) %>%
           head(n = 1) %>% 
           pull(), # select most sampled lake
         GSYear = "2018", 
         Log10PrevAreaCovered_ha = mean(Log10PrevAreaCovered_ha)) %>%
  mutate(pred1 = predict(inv_chg_mod1, newdata = .))

# figure
inv_chg_pred %>%
  ggplot(aes(x = LogitPropTreated, y = pred1)) +
  geom_point(data = inv_treat_new2, alpha = 0.3, aes(y = LogRatioCovered)) +
  geom_line(color = "#009193", size = 1.5) +
  # geom_text(data = inv_treat_change, y = 6.5, aes(label = ChangeProp4), 
  #           parse = T, check_overlap = T, hjust = 1, vjust = 1, size = 4) +
  facet_wrap(~ CommonName, scales = "free_x") +
  xlab("Proportion of waterbody treated (log-odds)") +
  ylab("Change in area covered by species (log-ratio)") +
  def_theme_facet +
  theme(axis.text.x = element_text(size = 10, color="black"))


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

