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
library(GGally)
library(vegan)
library(ggvegan)
# library(lfe) # FE models
# library(alpaca) # FE models
# library(lme4) # use glmmTMB unless it's too slow
# library(ggeffects)

# import data
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
ctrl_old <- read_csv("intermediate-data/FWC_control_old_formatted.csv")
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
herb_type <- read_csv("intermediate-data/herbicide_types.csv")
qual <- read_csv("intermediate-data/LW_quality_formatted.csv",
                 col_types = list(PermanentID = col_character(),
                                  JoinNotes = col_character()))

# 
# herb_type <- read_csv("intermediate-data/herbicide_types.csv")
# plant_detect <- read_csv("intermediate-data/FWC_plant_survey_first_detection_manual.csv")

# source functions
source("code/figure_settings.R")
source("code/mode_function.R")
source("code/proportion_transformations.R")
source("code/abundance_duplicates.R")
source("code/cumulative_herbicides.R")
source("code/cumulative_abundance.R")
source("code/invasion_interval.R")
source("code/native_interval.R")
source("code/okeechobee_growth.R")
source("code/surveyor_experience.R")
source("code/plant_abundance_formatting.R")
# source("code/herbicide_formatting.R")
source("code/control_old_formatting.R")
# source("code/control_new_formatting.R")
source("code/herbicide_new_formatting.R")
source("code/non_herbicide_new_formatting.R")
source("code/abundance_old_control_formatting.R")
source("code/abundance_herbicide_formatting.R")
source("code/abundance_non_herbicide_formatting.R")
source("code/growth_model.R")
# source("code/cumulative_quality.R")
source("code/quality_formatting.R")
source("code/abundance_quality_formatting.R")
source("code/growth_quality_model.R")

# prevent scientific notation
options(scipen = 999)


#### survey timing ####

# # extract info
# survey_time <- plant_fwc %>%
#   select(AreaOfInterestID, SurveyDate) %>%
#   unique() %>%
#   mutate(SurveyMonth = month(SurveyDate)) %>%
#   group_by(AreaOfInterestID) %>%
#   mutate(MostMonth = mean(Modes(SurveyMonth)),
#          Surveys = n()) %>% # average multiple modes
#   ungroup() %>%
#   mutate(MonthDiff = as.numeric(SurveyMonth) - as.numeric(MostMonth),
#          MonthDiff = case_when(MonthDiff < 0 ~ ceiling(MonthDiff), # reduce differences for multiple modes
#                                MonthDiff > 0 ~ floor(MonthDiff), # same as above
#                                TRUE ~ MonthDiff),
#          MonthDiff = case_when(MonthDiff > 6 ~ MonthDiff - 12, # correct for late/early surveys
#                                MonthDiff < -6 ~ MonthDiff + 12,
#                                TRUE ~ MonthDiff)) %>%
#   filter(Surveys >= 3)
# 
# # most month
# survey_time %>%
#   select(AreaOfInterestID, MostMonth) %>%
#   unique() %>%
#   ggplot(aes(x = MostMonth)) +
#   geom_histogram(binwidth = 1)
# 
# # figure
# pdf("output/survey_timing.pdf", width = 3.5, height = 3.5)
# ggplot(survey_time, aes(x = MonthDiff)) +
#   geom_histogram(binwidth = 1) +
#   labs(x = "Difference from most frequent month", y = "Surveys") +
#   def_theme_paper
# dev.off()
# 
# # non-zero months
# survey_time %>%
#   mutate(OnTime = if_else(MonthDiff == 0, 1, 0)) %>%
#   group_by(OnTime) %>%
#   count()


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


#### control data formatting ####

herb_taxa <- tibble(Species = c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)"))

# ctrl <- herbicide_dataset(ctrl_old, ctrl_new, herb_taxa) # new and old control datasets combined
ctrl_old2 <- ctrl_old_dataset(ctrl_old, herb_taxa) # note that this includes herbicide and non-herbicide methods
herb_new <- herbicide_new_dataset(ctrl_new, herb_taxa)
non_herb_new <- non_herb_new_dataset(ctrl_new, herb_taxa)

# overlap
ctrl_old2 %>%
  filter(AreaTreated_ha > 0) %>%
  select(PermanentID, Species, GSYear) %>%
  unique() %>%
  inner_join(herb_new %>%
               filter(AreaTreated_ha > 0) %>%
               select(PermanentID, Species, GSYear) %>%
               unique())
# yes, 2010 is repeated in both, 154 cases
# remove before joining with abundance

# make sure appropriate years are covered
range(herb_new$GSYear)
range(non_herb_new$GSYear)


#### combine invasive plant and ctrl data ####

herb_inv_taxa <- tibble(Species = c("Hydrilla verticillata", rep("Floating Plants (Eichhornia and Pistia)", 2)), # double each row that has floating plants
                        TaxonName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"))

# adjust gs year for old ctrl
inv_ctrl_old <- old_ctrl_abundance_dataset(ctrl_old2, inv_fwc3, herb_inv_taxa)

# format herbicide and non-herbicide datasets
inv_herb_new <- herbicide_abundance_dataset(ctrl_new, inv_ctrl_old, herb_new, inv_fwc3, herb_inv_taxa)
inv_ctrl_new <- non_herb_abundance_dataset(ctrl_new, inv_ctrl_old, non_herb_new, inv_herb_new, herb_inv_taxa)
# note that the above contains missing control data


#### combine invasive plant and quality data ####

# format quality data
qual2 <- quality_dataset(qual, 9)

# select variables for PCA 
qual_vals <- qual2 %>%
  select(TP_ug_L, TN_ug_L, CHL_ug_L, Secchi) %>%
  drop_na %>%
  mutate(across(.fns = log))

# visualize variables
ggpairs(qual_vals)

# PCA
qual_rda <- rda(qual_vals)
summary(qual_rda)
summary(qual_rda)$species

# format PCA scores
qual_scores <- data.frame(summary(qual_rda)$sites)

# output dataset
qual3 <- qual2 %>%
  select(PermanentID, GSYear, TP_ug_L, TN_ug_L, CHL_ug_L, Secchi) %>%
  drop_na %>%
  bind_cols(qual_vals %>%
              rename_with(~ paste0("log_", .x)) %>%
              bind_cols(select(qual_scores, PC1)))


# combine data
inv_ctrl_new_qual <- quality_abundance_dataset(qual3, inv_ctrl_new)

# sample size with quality info
inv_ctrl_new_qual %>%
  filter(!is.na(TN_ug_L) & !is.na(AreaCovered_ha) & !is.na(PropHerbTreated)) %>%
  summarise(waterbodies = n_distinct(PermanentID),
            years = n_distinct(GSYear),
            obs = n())


#### how does invasion affect herbicide treatment? ####

# subset for lakes with invasions and available herbicide data
inv_treat_new <- inv_ctrl_new_qual %>%
  filter(!is.na(PropHerbTreated) & PrevSpeciesPresent == 1) %>%
  mutate(Treated = fct_recode(as.character(HerbTreated), "Untreated" = "0", "Treated" = "1"),
         AreaHerbTreatedAdj_ha = case_when(AreaHerbTreated_ha > Area_ha ~ Area_ha, # if area treated in a year exceeds lake size, use lake size
                                           TRUE ~ AreaHerbTreated_ha),
         AreaHerbTreatedRounded_ha = round(AreaHerbTreatedAdj_ha), # round for glm
         AreaHerbUntreatedRounded_ha = round(Area_ha) - AreaHerbTreatedRounded_ha,
         PropHerbTreatedAdj = AreaHerbTreatedAdj_ha/Area_ha,
         Log10Area_ha = log10(Area_ha),
         LogitPropHerbTreated = logit(PropHerbTreatedAdj, adjust = prop_adjust),
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
ggplot(inv_treat_new2, aes(x = PropHerbTreatedAdj, y = LogitPropHerbTreated)) +
  geom_point()

# invasion/treatment relationships
ggplot(inv_treat_new2, aes(x = PrevPropCovered, y = PropHerbTreatedAdj)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new2, aes(x = Log10PrevAreaCovered_ha, y = PropHerbTreatedAdj)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new2, aes(x = PrevAreaCoveredCS_ha, y = PropHerbTreatedAdj)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new2, aes(x = PrevPropCovered, y = LogitPropHerbTreated)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new2, aes(x = Log10PrevAreaCovered_ha, y = LogitPropHerbTreated)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new2, aes(x = LogitPrevPropCovered, y = LogitPropHerbTreated)) +
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
inv_herb_mod1 <- glmmTMB(cbind(AreaHerbTreatedRounded_ha, AreaHerbUntreatedRounded_ha) ~ Log10PrevAreaCovered_ha * CommonName + Log10Area_ha + (1|PermanentID) + (1|GSYear), data = inv_treat_new2, family= "binomial")
summary(inv_herb_mod1)
# larger invasions have more proportion (lake) treated
# larger lakes have less proportion treated

inv_herb_mod2 <- glmmTMB(cbind(AreaHerbTreatedRounded_ha, AreaHerbUntreatedRounded_ha) ~ PrevPropCovered * CommonName + Log10Area_ha + (1|PermanentID) + (1|GSYear), data = inv_treat_new2, family= "binomial")
summary(inv_herb_mod2)
# larger invasions have more proportion (lake) treated
# larger lakes have less proportion treated

inv_herb_mod3 <- glmmTMB(cbind(AreaHerbTreatedRounded_ha, AreaHerbUntreatedRounded_ha) ~ LogitPrevPropCovered * CommonName + Log10Area_ha + (1|PermanentID) + (1|GSYear), data = inv_treat_new2, family= "binomial")
summary(inv_herb_mod3)
# same as above

# fixed-effects herbicide models
# no fe methods available for cbind(success, failure) or beta distributions
inv_herb_mod4 <- feols(LogitPropHerbTreated ~ Log10PrevAreaCovered_ha * CommonName | PermanentID + GSYear, data = inv_treat_new2)
summary(inv_herb_mod4)
# Log10Area_ha and AreaGroup are collinear with PermanentID
# invasion size increases treatment

inv_herb_mod5 <- feols(LogitPropHerbTreated ~ PrevPropCovered * CommonName | PermanentID + GSYear, data = inv_treat_new2)
summary(inv_herb_mod5)
# invasion size increases treatment

inv_herb_mod6 <- feols(LogitPropHerbTreated ~ LogitPrevPropCovered * CommonName | PermanentID + GSYear, data = inv_treat_new2)
summary(inv_herb_mod6)
# same as above

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
  geom_point(data = inv_treat_new2, alpha = 0.3, aes(y = LogitPropHerbTreated)) +
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
  geom_point(data = inv_treat_new2, alpha = 0.3, aes(y = LogitPropHerbTreated)) +
  geom_line(color = "#009193", size = 1.5) +
  geom_text(data = inv_treat_change, y = 6.5, aes(label = ChangeProp6), 
            parse = T, check_overlap = T, hjust = 1, vjust = 1, size = 4,
            nudge_x = -0.7) +
  facet_wrap(~ CommonName, scales = "free_x") +
  xlab("Proportion of waterbody covered by species (log-odds)") +
  ylab("Proportion of waterbody treated (log-odds)") +
  def_theme_facet
dev.off()


#### how does herbicide treatment affect invasion? ####

# previous prop groups
prev_prop_group <- inv_treat_new2 %>%
  select(CommonName, LogitPrevPropCovered) %>%
  unique() %>%
  group_by(CommonName) %>%
  mutate(PrevPropGroup = cut_interval(LogitPrevPropCovered, n = 3)) %>%
  ungroup()

levels(prev_prop_group$PrevPropGroup) <- rep(c("small", "medium", "large"), 3)

# subset for lakes with non-herbicide info
inv_treat_new3 <- inv_treat_new2 %>%
  filter(!is.na(PropNonHerbTreated)) %>% # should be same number of rows as inv_treat_new2
  mutate(AreaNonHerbTreatedAdj_ha = case_when(AreaNonHerbTreated_ha > Area_ha ~ Area_ha, # if area treated in a year exceeds lake size, use lake size
                                              TRUE ~ AreaNonHerbTreated_ha),
         PropNonHerbTreatedAdj = AreaNonHerbTreatedAdj_ha/Area_ha,
         LogitPropNonHerbTreated = logit(PropNonHerbTreatedAdj, adjust = prop_adjust),
         PerHaCoverChange = (EstAreaCoveredRaw_ha - PrevAreaCoveredRaw_ha) / PrevAreaCoveredRaw_ha,
         LogitPropCovered = logit(PropCovered, adjust = prop_adjust),
         LogRatioGrowth = if_else(LogRatioCovered > 0, 1, 0)) %>%
  left_join(prev_prop_group)

# dataset for invasive plants
inv_chg <- inv_treat_new3 %>%
  filter(CommonName == "Water hyacinth") %>%
  select(PermanentID, GSYear, LogRatioCovered) %>%
  rename("WaterHyacinth" = "LogRatioCovered") %>%
  full_join(inv_treat_new3 %>%
              filter(CommonName == "Water lettuce") %>%
              select(PermanentID, GSYear, LogRatioCovered) %>%
              rename("WaterLettuce" = "LogRatioCovered")) %>%
  full_join(inv_treat_new3 %>%
              filter(CommonName == "Hydrilla") %>%
              select(PermanentID, GSYear, LogRatioCovered) %>%
              rename("Hydrilla" = "LogRatioCovered"))

# plant relationships
ggpairs(inv_chg %>% select(-c(PermanentID, GSYear)))
# missing values from years where initial abundance of one is zero
# growth is positively correlated

# treatment variables
ggplot(inv_treat_new3, aes(x = LogitPropHerbTreated, y = HerbTreatmentDays)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)
# positively correlated

inv_treat_new3 %>%
  filter(PerHaCoverChange < 100) %>%
  ggplot(aes(x = PerHaCoverChange)) +
  geom_histogram(binwidth = 1)

ggplot(inv_treat_new3, aes(x = LogitPropHerbTreated, y = SurveyTreatDays)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)
# negatively correlated because untreated lakes necessary have long time windows (time between surveys, ~ 1 year)

# regression relationships
ggplot(inv_treat_new3, aes(x = LogitPrevPropCovered, y = LogRatioCovered, color = as.factor(HerbTreated))) +
  # geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new3, aes(x = LogitPrevPropCovered, y = PerHaCoverChange, color = as.factor(HerbTreated))) +
  # geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")
# looks similar to above

ggplot(inv_treat_new3, aes(x = LogitPrevPropCovered, y = LogitPropCovered, color = as.factor(HerbTreated))) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new3, aes(x = PrevPropGroup, y = LogRatioCovered)) +
  geom_boxplot() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new3, aes(x = LogitPropHerbTreated, y = LogRatioCovered)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_treat_new3, aes(x = LogitPropHerbTreated, y = PerHaCoverChange)) +
  # geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")
# water lettuce becomes negative

ggplot(inv_treat_new3, aes(x = LogitPropHerbTreated, y = LogRatioCovered)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(CommonName ~ PrevPropGroup, scales = "free")
# herbicide has positive effects even after taking initial abundance into consideration

ggplot(inv_treat_new3, aes(x = LogitPropHerbTreated, y = LogRatioGrowth)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "glm", formula = y ~ x) +
  facet_grid(CommonName ~ PrevPropGroup, scales = "free")
# still positive effects of herbicides for floating plants

ggplot(inv_treat_new3, aes(x = LogitPropNonHerbTreated, y = LogRatioCovered)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ CommonName, scales = "free")
# non-herb methods not that useful for water hyacinth and lettuce (very few)

# divide data by species
# without intercept in FE model, not sure how to interpret species coefficients
hyd_treat_new <- inv_treat_new3 %>%
  filter(CommonName == "Hydrilla") %>%
  mutate(SurveyTreatDaysCS = (SurveyTreatDays - mean(SurveyTreatDays)) / sd(SurveyTreatDays))
why_treat_new <- inv_treat_new3 %>%
  filter(CommonName == "Water hyacinth") %>%
  mutate(SurveyTreatDaysCS = (SurveyTreatDays - mean(SurveyTreatDays)) / sd(SurveyTreatDays))
wle_treat_new <- inv_treat_new3 %>%
  filter(CommonName == "Water lettuce") %>%
  mutate(SurveyTreatDaysCS = (SurveyTreatDays - mean(SurveyTreatDays)) / sd(SurveyTreatDays))

# treatment days and herbicide used -- correlated?
cor.test(~ SurveyTreatDaysCS + LogitPropHerbTreated, data = hyd_treat_new)
cor.test(~ SurveyTreatDaysCS + LogitPropHerbTreated, data = why_treat_new)
cor.test(~ SurveyTreatDaysCS + LogitPropHerbTreated, data = wle_treat_new)

# dictionary for coefficient plot
dict <- c("LogitPropHerbTreated" = "Herbicide",
          "LogitPrevPropCovered" = "Initial invasion",
          "SurveyTreatDaysCS" = "Days post treatment",
          "SurveyorExperienceBmedium" = "Med. surveyor experience",
          "SurveyorExperienceBlow" = "Low surveyor experience",
          "LogitPropHerbTreated:LogitPrevPropCovered" = "Herbicide:Initial invasion",
          "LogitPropHerbTreated:SurveyTreatDaysCS" = "Herbicide:Days post treatment",
          "LogitPropHerbTreated:SurveyorExperienceBmedium" = "Herbicide:Med. experience",
          "LogitPropHerbTreated:SurveyorExperienceBlow" = "Herbicide:Low experience",
          "PC1" = "Water quality",
          "LogitPropHerbTreated:PC1" = "Herbicide:Water quality")
setFixest_coefplot(dict = dict)

# fixed effect models
inv_chg_mod1 <- feols(LogRatioCovered ~ (LogitPropHerbTreated + LogitPrevPropCovered) * CommonName | PermanentID + GSYear, data = inv_treat_new3)
summary(inv_chg_mod1)

hyd_chg_mod1 <- feols(LogRatioCovered ~ LogitPropHerbTreated * LogitPrevPropCovered | PermanentID + GSYear, data = hyd_treat_new)
summary(hyd_chg_mod1)

hyd_chg_mod1b <- feols(LogRatioCovered ~ LogitPrevPropCovered * (LogitPropHerbTreated + LogitPropNonHerbTreated) | PermanentID + GSYear, data = hyd_treat_new)
summary(hyd_chg_mod1b)
# estimate for herbicide effect doesn't change much when non-herbicide treatments are accounted for

why_chg_mod1 <- feols(LogRatioCovered ~ LogitPropHerbTreated * LogitPrevPropCovered | PermanentID + GSYear, data = why_treat_new)
summary(why_chg_mod1)

wle_chg_mod1 <- feols(LogRatioCovered ~ LogitPropHerbTreated * LogitPrevPropCovered | PermanentID + GSYear, data = wle_treat_new)
summary(wle_chg_mod1)

pdf("output/fixed_effect_change_models_simple_coefficient_plots.pdf")
coefplot(hyd_chg_mod1, horiz = T, main = "Hydrilla annual growth")
coefplot(why_chg_mod1, horiz = T, main = "Water hyacinth annual growth")
coefplot(wle_chg_mod1, horiz = T, main = "Water lettuce annual growth")
dev.off()

hyd_chg_mod2 <- feols(LogitPropCovered ~ LogitPropHerbTreated * LogitPrevPropCovered | PermanentID + GSYear, data = hyd_treat_new)
summary(hyd_chg_mod2)

why_chg_mod2 <- feols(LogitPropCovered ~ LogitPropHerbTreated * LogitPrevPropCovered | PermanentID + GSYear, data = why_treat_new)
summary(why_chg_mod2)

wle_chg_mod2 <- feols(LogitPropCovered ~ LogitPropHerbTreated * LogitPrevPropCovered | PermanentID + GSYear, data = wle_treat_new)
summary(wle_chg_mod2)

pdf("output/fixed_effect_cover_models_simple_coefficient_plots.pdf")
coefplot(hyd_chg_mod2, horiz = T, main = "Hydrilla annual growth")
coefplot(why_chg_mod2, horiz = T, main = "Water hyacinth annual growth")
coefplot(wle_chg_mod2, horiz = T, main = "Water lettuce annual growth")
dev.off()


# mixed effect models
# hyd_chg_mod2 <- glmmTMB(LogRatioCovered ~ LogitPropHerbTreated + LogitPrevPropCovered + (1|PermanentID) + (1|GSYear), data = hyd_treat_new)
# summary(hyd_chg_mod2)
# 
# why_chg_mod2 <- glmmTMB(LogRatioCovered ~ LogitPropHerbTreated + LogitPrevPropCovered + (1|PermanentID) + (1|GSYear), data = why_treat_new)
# summary(why_chg_mod2)
# 
# wle_chg_mod2 <- glmmTMB(LogRatioCovered ~ LogitPropHerbTreated + LogitPrevPropCovered + (1|PermanentID) + (1|GSYear), data = wle_treat_new)
# summary(wle_chg_mod2)

# predicted data
inv_chg_pred <- inv_treat_new2 %>%
  select(CommonName, LogitPropHerbTreated) %>%
  unique() %>%
  left_join(inv_treat_new2 %>%
              group_by(CommonName) %>%
              summarise(LogitPrevPropCovered = mean(LogitPrevPropCovered)) %>%
              ungroup()) %>%
  mutate(PermanentID = inv_treat_new2 %>%
           group_by(PermanentID) %>%
           count() %>%
           arrange(desc(n)) %>%
           select(PermanentID) %>%
           head(n = 1) %>% 
           pull(), # select most sampled lake
         GSYear = "2018") %>%
  mutate(pred1 = predict(inv_chg_mod1, newdata = .),
         pValue = case_when(CommonName == "Hydrilla" ~ "P = 0.12",
                            CommonName == "Water hyacinth" ~ "P = 0.09",
                            CommonName == "Water lettuce" ~ "P < 0.001"))

# figure
pdf("output/new_invasion_herbicide.pdf", width = 11, height = 4)
inv_chg_pred %>%
  ggplot(aes(x = LogitPropHerbTreated, y = pred1)) +
  geom_hline(yintercept = 0) +
  geom_point(data = inv_treat_new2, alpha = 0.3, aes(y = LogRatioCovered)) +
  geom_line(color = "#009193", size = 1.5) +
  geom_text(aes(label = pValue, x = max(LogitPropHerbTreated)), y = 7.5, check_overlap = T, 
            hjust = 1, vjust = 1, size = 4) +
  facet_wrap(~ CommonName, scales = "free_x") +
  xlab("Proportion of waterbody treated (log-odds)") +
  ylab("Change in area covered by sp. (log-ratio)") +
  def_theme_facet +
  theme(axis.text.x = element_text(size = 10, color="black"))
dev.off()

# add time between treatment and survey
# only use treated lakes

# subset data
hyd_treat_new2 <- hyd_treat_new %>% filter(HerbTreated == 1) %>%
  mutate(SurveyTreatDaysCS = (SurveyTreatDays - mean(SurveyTreatDays)) / sd(SurveyTreatDays))
why_treat_new2 <- why_treat_new %>% filter(HerbTreated == 1) %>%
  mutate(SurveyTreatDaysCS = (SurveyTreatDays - mean(SurveyTreatDays)) / sd(SurveyTreatDays))
wle_treat_new2 <- wle_treat_new %>% filter(HerbTreated == 1) %>%
  mutate(SurveyTreatDaysCS = (SurveyTreatDays - mean(SurveyTreatDays)) / sd(SurveyTreatDays))

# figures
ggplot(hyd_treat_new2, aes(x = SurveyTreatDaysCS, y = LogitPropHerbTreated)) +
  geom_point() + geom_smooth(method = "lm")
cor.test(~ SurveyTreatDaysCS + LogitPropHerbTreated, data = hyd_treat_new2)
ggplot(why_treat_new2, aes(x = SurveyTreatDaysCS, y = LogitPropHerbTreated)) +
  geom_point() + geom_smooth(method = "lm")
cor.test(~ SurveyTreatDaysCS + LogitPropHerbTreated, data = why_treat_new2)
ggplot(wle_treat_new2, aes(x = SurveyTreatDaysCS, y = LogitPropHerbTreated)) +
  geom_point() + geom_smooth(method = "lm")
cor.test(~ SurveyTreatDaysCS + LogitPropHerbTreated, data = wle_treat_new2)
# all are negatively correlated, but correlation magnitudes are small (< -0.3)

ggplot(hyd_treat_new2, aes(x = SurveyTreatDaysCS, y = LogRatioCovered)) +
  geom_point() + geom_smooth(method = "lm")
ggplot(why_treat_new2, aes(x = SurveyTreatDaysCS, y = LogRatioCovered)) +
  geom_point() + geom_smooth(method = "lm")
ggplot(wle_treat_new2, aes(x = SurveyTreatDaysCS, y = LogRatioCovered)) +
  geom_point() + geom_smooth(method = "lm")

# fixed effect models
hyd_chg_mod3 <- feols(LogRatioCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB) | PermanentID + GSYear, data = hyd_treat_new2)
summary(hyd_chg_mod3)

why_chg_mod3 <- feols(LogRatioCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB) | PermanentID + GSYear, data = why_treat_new2)
summary(why_chg_mod3)

wle_chg_mod3 <- feols(LogRatioCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB) | PermanentID + GSYear, data = wle_treat_new2)
summary(wle_chg_mod3)

pdf("output/fixed_effect_change_models_coefficient_plots.pdf")
coefplot(hyd_chg_mod3, horiz = T, main = "Hydrilla annual growth")
coefplot(why_chg_mod3, horiz = T, main = "Water hyacinth annual growth")
coefplot(wle_chg_mod3, horiz = T, main = "Water lettuce annual growth")
dev.off()

# fixed effect models with prop covered as response
hyd_cov_mod <- feols(LogitPropCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB) | PermanentID + GSYear, data = hyd_treat_new2)
summary(hyd_cov_mod)

why_cov_mod <- feols(LogitPropCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB) | PermanentID + GSYear, data = why_treat_new2)
summary(why_cov_mod)

wle_cov_mod <- feols(LogitPropCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB) | PermanentID + GSYear, data = wle_treat_new2)
summary(wle_cov_mod)

pdf("output/fixed_effect_cover_models_coefficient_plots.pdf")
coefplot(hyd_cov_mod, horiz = T, main = "Hydrilla abundance")
coefplot(why_cov_mod, horiz = T, main = "Water hyacinth abundance")
coefplot(wle_cov_mod, horiz = T, main = "Water lettuce abundance")
dev.off()


#### how does water quality mediate the effects of herbicides? ####

# subset data
hyd_treat_new3 <- hyd_treat_new2 %>% 
  filter(!is.na(PC1))
why_treat_new3 <- why_treat_new2 %>% 
  filter(!is.na(PC1))
wle_treat_new3 <- wle_treat_new2 %>% 
  filter(!is.na(PC1))

# correlations
hyd_treat_new3 %>%
  select(TN_ug_L, TP_ug_L, CHL_ug_L, Secchi) %>%
  ggpairs()
why_treat_new3 %>%
  select(TN_ug_L, TP_ug_L, CHL_ug_L, Secchi) %>%
  ggpairs()
wle_treat_new3 %>%
  select(TN_ug_L, TP_ug_L, CHL_ug_L, Secchi) %>%
  ggpairs()
# many are highly correlated, use PC1

# fixed effect models
hyd_chg_mod4 <- feols(LogRatioCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB + PC1) | PermanentID + GSYear, data = hyd_treat_new3)
summary(hyd_chg_mod4)

why_chg_mod4 <- feols(LogRatioCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB + PC1) | PermanentID + GSYear, data = why_treat_new3)
summary(why_chg_mod4)

wle_chg_mod4 <- feols(LogRatioCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB + PC1) | PermanentID + GSYear, data = wle_treat_new3)
summary(wle_chg_mod4)

pdf("output/fixed_effect_change_quality_models_coefficient_plots.pdf")
coefplot(hyd_chg_mod4, horiz = T, main = "Hydrilla annual growth")
coefplot(why_chg_mod4, horiz = T, main = "Water hyacinth annual growth")
coefplot(wle_chg_mod4, horiz = T, main = "Water lettuce annual growth")
dev.off()


# fixed effect models with prop covered as response
hyd_cov_mod2 <- feols(LogitPropCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB + PC1) | PermanentID + GSYear, data = hyd_treat_new3)
summary(hyd_cov_mod2)

why_cov_mod2 <- feols(LogitPropCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB + PC1) | PermanentID + GSYear, data = why_treat_new3)
summary(why_cov_mod2)

wle_cov_mod2 <- feols(LogitPropCovered ~ LogitPropHerbTreated * (LogitPrevPropCovered + SurveyTreatDaysCS + SurveyorExperienceB + PC1) | PermanentID + GSYear, data = wle_treat_new3)
summary(wle_cov_mod2)

pdf("output/fixed_effect_cover_quality_models_coefficient_plots.pdf")
coefplot(hyd_cov_mod2, horiz = T, main = "Hydrilla abundance")
coefplot(why_cov_mod2, horiz = T, main = "Water hyacinth abundance")
coefplot(wle_cov_mod2, horiz = T, main = "Water lettuce abundance")
dev.off()


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

