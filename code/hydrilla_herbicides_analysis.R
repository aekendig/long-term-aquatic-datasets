#### info ####

# goal: evaluate effects of herbicides on hydrilla


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)

# import data
ctrl_old <- read_csv("intermediate-data/FWC_control_old_formatted.csv")
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
plant_lw <- read_csv("intermediate-data/LW_plant_formatted.csv")


#### edit data ####

# select data associated with Hydrilla
# convert all to hectares
ctrl_old_hyd <- ctrl_old %>%
  filter(Species == "Hydrilla verticillata") %>%
  mutate(AreaTreated_ha= TotalAcres * 0.405,
         Area_ha = ShapeArea * 100)

ctrl_new_hyd <- ctrl_new %>%
  filter(Species == "Hydrilla verticillata") %>%
  mutate(AreaTreated_ha = TotalAcres * 0.405,
         Area_ha = ShapeArea * 100)

plant_fwc_hyd <- plant_fwc %>%
  filter(SpeciesName == "Hydrilla verticillata") %>%
  mutate(AreaCovered_ha = SpeciesAcres * 0.405,
         Area_ha = ShapeArea * 100,
         SpeciesFrequency_ha = AreaCovered_ha / Area_ha)

plant_lw_hyd <- plant_lw %>%
  filter(GenusSpecies == "Hydrilla verticillata") %>%
  mutate(Area_ha = ShapeArea * 100)

# see if species frequency depends on area source in FWC
ggplot(plant_fwc_hyd, aes(x = SpeciesFrequency, y = SpeciesFrequency_ha)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point()
# a lot more are greater than 1
# cluster around zero: probably total shape area includes larger waterbody connected to smaller one

plant_fwc_hyd %>%
  filter(SpeciesFrequency_ha > 1) %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, AreaCovered_ha, SpeciesFrequency_ha) %>%
  unique() %>%
  data.frame()
# may be inaccurate area estimates - round down to 1

# see if treated area exceeds total area
ctrl_old_hyd %>%
  filter(AreaTreated_ha > Area_ha) %>%
  select(AreaOfInterest, County_FWC, AreaOfInterestID, PermanentID, AreaTreated_ha, Area_ha) %>%
  unique() %>%
  data.frame()
# checked a few of the large differences and the areas exceed reported areas on websites - round down to 1

# how many lw plant surveys are relevant?
plant_lw_hyd %>%
  select(PermanentID, Date) %>%
  unique() %>%
  inner_join(ctrl_old_hyd %>%
               select(PermanentID) %>%
               unique() %>%
               full_join(ctrl_new_hyd %>%
                           select(PermanentID) %>%
                           unique())) %>% # 184 surveys
  select(PermanentID) %>%
  unique() %>% # 51 lakes
  anti_join(plant_fwc_hyd %>%
               select(PermanentID) %>%
               unique) # 1 lake was not sampled by FWC, the other 50 were

# how consistent are LW and FWC plant surveys?
plant_lw_hyd %>%
  select(Year, PermanentID, SpeciesFrequency) %>%
  inner_join(plant_fwc_hyd %>%
               select(SurveyYear, PermanentID, SpeciesFrequency_ha) %>%
               rename(Year = SurveyYear)) %>%
  ggplot(aes(x = SpeciesFrequency, y = SpeciesFrequency_ha)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  xlab("LW frequency") +
  ylab("FWC frequency")
# generally much lower estimates from FWC than LW

# check plant dataset for unique ID-name combos
plant_fwc_hyd %>%
  select(AreaOfInterest, AreaOfInterestID) %>%
  unique() %>%
  mutate(dup = duplicated(AreaOfInterestID)) %>%
  filter(dup == T)

# FWC average plant change
plant_change_fwc_hyd <- plant_fwc_hyd %>%
  mutate(AreaCovered_ha = case_when(AreaCovered_ha > Area_ha ~ Area_ha,
                                      TRUE ~ AreaCovered_ha),
         SpeciesFrequency_ha = AreaCovered_ha/Area_ha) %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha) %>%
  summarise(AreaInitial_ha = AreaCovered_ha[SurveyDate == min(SurveyDate)],
            AreaChange_ha = coef(lm(AreaCovered_ha ~ SurveyDate))[2] * 365,
            FreqChange = coef(lm(SpeciesFrequency_ha ~ SurveyDate))[2] * 365,
            DaysChange = max(SurveyDate) - min(SurveyDate),
            LastSurvey = max(SurveyDate)) %>%
  ungroup()

# FWC survey-to-survey plant change
survey_change_fwc_hyd <- plant_fwc_hyd %>%
  mutate(AreaCovered_ha = case_when(AreaCovered_ha > Area_ha ~ Area_ha,
                                       TRUE ~ AreaCovered_ha)) %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, SurveyDate, AreaCovered_ha) %>%
  arrange(AreaOfInterestID, SurveyDate) %>%
  group_by(AreaOfInterestID, PermanentID) %>%
  mutate(DaysChange = lead(SurveyDate) - SurveyDate,
         AreaChange_ha = (lead(AreaCovered_ha) - AreaCovered_ha)/as.numeric(DaysChange) * 365,
         SpeciesFrequency_ha = AreaCovered_ha/Area_ha,
         FrequencyChange = (lead(SpeciesFrequency_ha) - SpeciesFrequency_ha)/as.numeric(DaysChange) * 365) %>%
  ungroup()

# check ctrl dataset for unique ID-name combos
ctrl_old_hyd %>%
  select(AreaOfInterest, AreaOfInterestID) %>%
  unique() %>%
  mutate(dup = duplicated(AreaOfInterestID)) %>%
  filter(dup == T)

ctrl_new_hyd %>%
  select(AreaOfInterest, AreaOfInterestID) %>%
  unique() %>%
  mutate(dup = duplicated(AreaOfInterestID)) %>%
  filter(dup == T)

# non-herbicide methods (from herbicide_initial_visualizations)
non_herb <- c("Mechanical Harvester", 
              "Snagging (tree removal)", 
              "Aquatic Dye (for shading)", 
              "Grass Carp", "Hand Removal", 
              "Mechanical (Other)", 
              "Mechanical Shredder", 
              "Prescribed Fire")

# FWC average herbicide and plants
plant_herb_fwc_hyd <- plant_change_fwc_hyd %>%
  inner_join(ctrl_old_hyd %>%
               filter(AreaTreated_ha > 0) %>%
               group_by(AreaOfInterestID, PermanentID, Year, Area_ha) %>%
               summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
               ungroup() %>%
               full_join(ctrl_new_hyd %>%
                           filter(AreaTreated_ha > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
                           group_by(AreaOfInterestID, PermanentID, Year, Area_ha) %>%
                           summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
                           ungroup()) %>%
               mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha,
                                                 TRUE ~ AreaTreated_ha),
                      YearTreated = paste0(Year, "-12-31") %>% as.Date("%Y-%m-%d"),
                      MinYear = min(Year)) %>%
               select(-c(Area_ha, Year))) %>%
  filter(YearTreated < LastSurvey) %>%
  mutate(YearsTotal = year(LastSurvey) - MinYear - 1) %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, AreaInitial_ha, AreaChange_ha, FreqChange, DaysChange) %>%
  summarise(TreatmentYears = unique(YearsTotal),
            TreatmentFrequency = n()/unique(YearsTotal),
            TreatmentIntensity = mean(AreaTreated_ha/Area_ha)) %>%
  ungroup()


#### start here: add treatment data to survey-to-survey data ####


#### plant_change_fwc_hyd exploratory figures ####

ggplot(plant_change_fwc_hyd, aes(DaysChange, AreaChange_ha)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(AreaInitial_ha/Area_ha, AreaChange_ha)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(Area_ha, AreaChange_ha)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(DaysChange, FreqChange_ha)) +
  geom_point()


#### survey_change_fwc_hyd exploratory figures ####

ggplot(survey_change_fwc_hyd, aes(x = AreaChange_ha)) +
  geom_histogram()

ggplot(survey_change_fwc_hyd, aes(x = FrequencyChange)) +
  geom_histogram()

ggplot(survey_change_fwc_hyd, aes(Area_ha, FrequencyChange)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(survey_change_fwc_hyd, aes(log(Area_ha), FrequencyChange)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(survey_change_fwc_hyd, aes(SpeciesFrequency_ha, FrequencyChange)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(survey_change_fwc_hyd, aes(SpeciesFrequency_ha, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(survey_change_fwc_hyd, aes(SpeciesFrequency_ha, AreaChange_ha, color = PermanentID)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.3) +
  theme(legend.position = "none")

ggplot(survey_change_fwc_hyd, aes(SpeciesFrequency_ha, AreaChange_ha, color = PermanentID)) +
  geom_smooth(method = "lm", se = F, size = 0.3) +
  theme(legend.position = "none")

ggplot(survey_change_fwc_hyd, aes(SurveyDate, AreaChange_ha, color = PermanentID)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.3) +
  theme(legend.position = "none")

ggplot(survey_change_fwc_hyd, aes(SurveyDate, AreaChange_ha, color = PermanentID)) +
  geom_smooth(method = "lm", se = F, size = 0.3) +
  theme(legend.position = "none")

ggplot(survey_change_fwc_hyd, aes(log(as.numeric(DaysChange)), AreaChange_ha, color = PermanentID)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.3) +
  theme(legend.position = "none")


#### plant_herb_fwc_hyd  exploratory figures ####

ggplot(plant_herb_fwc_hyd, aes(TreatmentFrequency, AreaChange_ha)) +
  geom_point()

ggplot(plant_herb_fwc_hyd, aes(TreatmentIntensity, AreaChange_ha)) +
  geom_point()

ggplot(plant_herb_fwc_hyd, aes(TreatmentFrequency, FreqChange)) +
  geom_point()

ggplot(plant_herb_fwc_hyd, aes(TreatmentIntensity, FreqChange)) +
  geom_point()
