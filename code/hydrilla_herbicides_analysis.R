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

# assumptions
MinHerbLag = 30
MaxHerbLag = 365 * 2


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

# see if species cover exceeds total area
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
            FirstSurvey = min(SurveyDate),
            LastSurvey = max(SurveyDate),
            DaysChange = max(SurveyDate) - min(SurveyDate)) %>%
  ungroup()

# FWC survey-to-survey plant change
survey_change_fwc_hyd <- plant_fwc_hyd %>%
  mutate(AreaCovered_ha = case_when(AreaCovered_ha > Area_ha ~ Area_ha,
                                       TRUE ~ AreaCovered_ha)) %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, SurveyDate, AreaCovered_ha) %>%
  arrange(AreaOfInterestID, SurveyDate) %>%
  group_by(AreaOfInterestID, PermanentID) %>%
  mutate(DaysChange = SurveyDate - lag(SurveyDate),
         OrigAreaCovered_ha = lag(AreaCovered_ha),
         AreaChange_ha = (AreaCovered_ha - lag(AreaCovered_ha))/as.numeric(DaysChange) * 365,
         SpeciesFrequency_ha = AreaCovered_ha/Area_ha,
         OrigFreq = lag(SpeciesFrequency_ha),
         FreqChange = (SpeciesFrequency_ha - lag(SpeciesFrequency_ha)),
         FreqChangePerYear = (SpeciesFrequency_ha - lag(SpeciesFrequency_ha))/as.numeric(DaysChange) * 365,
         FreqChangePreFreqPerYear = (SpeciesFrequency_ha - lag(SpeciesFrequency_ha))/(lag(SpeciesFrequency_ha) * as.numeric(DaysChange)) * 365) %>%
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
# assign 12/31 to the year herbicide was applied for old dataset (no date info)
# the earliest last survey is in April
# only include herbicide treatments before the last survey (4 month+ lag)
# for known dates, required MinHerbLag days before survey
# herbicide treatments more than MaxHerbLag days before the first survey are not included
plant_herb_fwc_hyd <- plant_change_fwc_hyd %>%
  left_join(ctrl_old_hyd %>%
              filter(AreaTreated_ha > 0) %>%
              group_by(AreaOfInterestID, PermanentID, Year, Area_ha) %>%
              summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
              ungroup() %>%
              mutate(TreatmentDate = paste0(Year, "-12-31") %>% as.Date("%Y-%m-%d")) %>%
              full_join(ctrl_new_hyd %>%
                          filter(AreaTreated_ha > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
                          group_by(AreaOfInterestID, PermanentID, BeginDate, Year, Area_ha) %>%
                          summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
                          ungroup() %>%
                          rename(TreatmentDate = BeginDate)) %>%
              mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha,
                                                TRUE ~ AreaTreated_ha),
                     MinYear = min(Year)) %>%
              select(-c(Area_ha, Year))) %>%
  mutate(AreaTreated_ha = case_when(FirstSurvey - TreatmentDate > MaxHerbLag ~ 0,
                                    LastSurvey - TreatmentDate <= MinHerbLag ~ 0,
                                    is.na(TreatmentDate) ~ 0,
                                    TRUE ~ AreaTreated_ha),
         TreatmentDate = case_when(AreaTreated_ha == 0 ~ NA_Date_,
                                 TRUE ~ TreatmentDate),
         MinYear = case_when(AreaTreated_ha == 0 ~ NA_real_,
                             TRUE ~ MinYear),
         YearsTotal = case_when(AreaTreated_ha > 0 ~ year(LastSurvey) - MinYear,
                                TRUE ~ 1)) %>%
  unique() %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID) %>%
  mutate(rows = n()) %>%
  ungroup() %>%
  filter(!(AreaTreated_ha == 0 & rows > 1)) %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, AreaInitial_ha, AreaChange_ha, FreqChange, DaysChange) %>%
  summarise(TreatmentYears = case_when(sum(AreaTreated_ha) > 0 ~ unique(YearsTotal),
                                       TRUE ~ 0),
            TreatmentFrequency = sum(AreaTreated_ha > 0)/unique(YearsTotal),
            TreatmentIntensity = mean(AreaTreated_ha/AreaInitial_ha),
            TreatmentIntensityLake = mean(AreaTreated_ha/Area_ha)) %>%
  ungroup()

# FWC survey-to-survey herbicide and plants
# herbicide treatments that occurred between the two surveys are included
# unless it is within MinHerbLag days of the last survey
survey_herb_fwc_hyd <- survey_change_fwc_hyd %>%
  left_join(ctrl_old_hyd %>%
              filter(AreaTreated_ha > 0) %>%
              group_by(AreaOfInterestID, PermanentID, Year, Area_ha) %>%
              summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
              ungroup() %>%
              mutate(TreatmentDate = paste0(Year, "-12-31") %>% as.Date("%Y-%m-%d")) %>%
              full_join(ctrl_new_hyd %>%
                          filter(AreaTreated_ha > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
                          group_by(AreaOfInterestID, PermanentID, BeginDate, Area_ha) %>%
                          summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
                          ungroup() %>%
                          rename(TreatmentDate = BeginDate)) %>%
              mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha,
                                                TRUE ~ AreaTreated_ha)) %>%
              select(-c(Area_ha, Year))) %>%
  mutate(SurvTreatDays = SurveyDate - TreatmentDate,
         AreaTreated_ha = case_when(SurvTreatDays > MinHerbLag & SurvTreatDays <= DaysChange ~ AreaTreated_ha,
                                    TRUE ~ 0),
         SurvTreatDays = case_when(AreaTreated_ha == 0 ~ NA_real_,
                                   TRUE ~ as.numeric(SurvTreatDays)),
         TreatmentDate = case_when(AreaTreated_ha == 0 ~ NA_Date_,
                                 TRUE ~ TreatmentDate)) %>%
  unique() %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, SurveyDate) %>%
  mutate(rows = n()) %>%
  ungroup() %>%
  filter(!(AreaTreated_ha == 0 & rows > 1)) %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, SurveyDate, AreaCovered_ha, DaysChange, AreaChange_ha, SpeciesFrequency_ha, FreqChange) %>%
  summarise(TreatmentYears = case_when(sum(AreaTreated_ha) > 0 ~ unique((as.numeric(DaysChange) - MinHerbLag)/365),
                                       TRUE ~ 0),
            TreatmentFrequency = sum(AreaTreated_ha > 0)/unique((as.numeric(DaysChange) - MinHerbLag) * 365),
            TreatmentIntensity = mean(AreaTreated_ha/OrigAreaCovered_ha),
            TreatmentIntensityLake = mean(AreaTreated_ha/Area_ha)) %>%
  ungroup()

# FWC survey-to-survey cumulative herbicide and plants
cltv_herb_fwc_hyd <- survey_change_fwc_hyd %>%
  left_join(ctrl_old_hyd %>%
              filter(AreaTreated_ha > 0) %>%
              group_by(AreaOfInterestID, PermanentID, Year, Area_ha) %>%
              summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
              ungroup() %>%
              mutate(TreatmentDate = paste0(Year, "-12-31") %>% as.Date("%Y-%m-%d")) %>%
              full_join(ctrl_new_hyd %>%
                          filter(AreaTreated_ha > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
                          group_by(AreaOfInterestID, PermanentID, BeginDate, Year, Area_ha) %>%
                          summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
                          ungroup() %>%
                          rename(TreatmentDate = BeginDate)) %>%
              mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha,
                                                TRUE ~ AreaTreated_ha),
                     MinYear = min(Year)) %>%
              select(-c(Area_ha, Year))) %>%
  mutate(SurvTreatDays = SurveyDate - TreatmentDate,
         AreaTreated_ha = case_when(SurvTreatDays > MinHerbLag & SurvTreatDays < MaxHerbLag & !is.na(DaysChange) ~ AreaTreated_ha,
                                    TRUE ~ 0),
         SurvTreatDays = case_when(AreaTreated_ha == 0 ~ NA_real_,
                                   TRUE ~ as.numeric(SurvTreatDays)),
         TreatmentDate = case_when(AreaTreated_ha == 0 ~ NA_Date_,
                                 TRUE ~ TreatmentDate),
         MinYear = case_when(AreaTreated_ha == 0 ~ NA_real_,
                             TRUE ~ MinYear)) %>%
  unique() %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, SurveyDate) %>%
  mutate(rows = n()) %>%
  ungroup() %>%
  filter(!(AreaTreated_ha == 0 & rows > 1)) %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, SurveyDate, AreaCovered_ha, DaysChange, AreaChange_ha, SpeciesFrequency_ha, FreqChange) %>%
  summarise(TreatmentFrequency = sum(AreaTreated_ha > 0)/MaxHerbLag,
            TreatmentIntensity = mean(AreaTreated_ha/OrigAreaCovered_ha),
            TreatmentIntensityLake = mean(AreaTreated_ha/Area_ha)) %>%
  ungroup()

# check LW data for freq > 1
plant_lw_hyd %>%
  filter(SpeciesFrequency > 1)

# LW survey-to-survey plant change
survey_change_lw_hyd <- plant_lw_hyd %>%
  select(Lake, PermanentID, County_LW, Area_ha, Date, SpeciesFrequency) %>%
  arrange(PermanentID, Lake, County_LW, Date) %>%
  group_by(PermanentID, Lake, County_LW, ) %>%
  mutate(DaysChange = Date - lag(Date),
         OrigFreq = lag(SpeciesFrequency),
         FreqChange = (SpeciesFrequency - lag(SpeciesFrequency)),
         FreqChangePerYear = (SpeciesFrequency - lag(SpeciesFrequency))/as.numeric(DaysChange) * 365,
         FreqChangePreFreqPerYear = (SpeciesFrequency - lag(SpeciesFrequency))/(lag(SpeciesFrequency) * as.numeric(DaysChange)) * 365) %>%
  ungroup() %>%
  mutate(uniqueID = paste(Lake, County_LW, PermanentID, sep = "_"))
  

#### plant_fwc_hyd exploratory figures ####

# small invasions over time
plant_fwc_hyd %>%
  filter(SpeciesAcres < 10) %>%
  ggplot(aes(x = SurveyDate, y = SpeciesAcres, color = as.factor(AreaOfInterestID))) +
  geom_line(alpha = 0.5) +
  theme(legend.position = "none")

# acres by surveyor
plant_fwc_hyd %>%
  ggplot(aes(x = SpeciesAcres)) +
  geom_histogram() +
  facet_wrap(~ Surveyor, scales = "free")

# acres by surveyor when Hydrilla present
plant_fwc_hyd %>%
  filter(SpeciesAcres > 0) %>%
  ggplot(aes(x = SpeciesAcres)) +
  geom_histogram() +
  facet_wrap(~ Surveyor, scales = "free")

# distribution of low values
plant_fwc_hyd %>%
  filter(SpeciesAcres < 20) %>%
  ggplot(aes(x = SpeciesAcres)) +
  geom_histogram(binwidth = 0.1)

# distribution of medium values
plant_fwc_hyd %>%
  filter(SpeciesAcres > 20 & SpeciesAcres < 500) %>%
  ggplot(aes(x = SpeciesAcres)) +
  geom_histogram(binwidth = 10)

# distribution of high values
plant_fwc_hyd %>%
  filter(SpeciesAcres > 500) %>%
  ggplot(aes(x = SpeciesAcres)) +
  geom_histogram(binwidth = 100)

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

#### plant_change_fwc_hyd exploratory figures ####

ggplot(plant_change_fwc_hyd, aes(DaysChange, AreaChange_ha)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(AreaInitial_ha/Area_ha, AreaChange_ha)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(Area_ha, AreaChange_ha)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(AreaInitial_ha/Area_ha, FreqChange)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(DaysChange, FreqChange)) +
  geom_point()


#### survey_change_fwc_hyd exploratory figures ####

ggplot(survey_change_fwc_hyd, aes(DaysChange, AreaChange_ha)) +
  geom_vline(xintercept = 365, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.5)

survey_change_fwc_hyd %>%
  filter(!is.na(AreaChange_ha)) %>%
  ggplot(aes(DaysChange, log(AreaCovered_ha))) +
  geom_hline(yintercept = log(0.0405), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 365, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.5)

ggplot(survey_change_fwc_hyd, aes(x = FreqChange)) +
  geom_histogram()

ggplot(survey_change_fwc_hyd, aes(x = FreqChangePerYear)) +
  geom_histogram()

ggplot(survey_change_fwc_hyd, aes(Area_ha, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(survey_change_fwc_hyd, aes(log(Area_ha), FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(survey_change_fwc_hyd, aes(OrigFreq, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm")

survey_change_fwc_hyd %>%
  mutate(OrigFreqBin = cut_interval(OrigFreq, n = 10)) %>%
  filter(!is.na(OrigFreqBin))  %>%
  ggplot(aes(x = FreqChange)) +
  geom_histogram() +
  facet_wrap(~ OrigFreqBin, scales = "free_y")

ggplot(survey_change_fwc_hyd, aes(DaysChange, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm")


#### plant_herb_fwc_hyd  exploratory figures ####

ggplot(plant_herb_fwc_hyd, aes(TreatmentFrequency, AreaChange_ha)) +
  geom_point()

ggplot(plant_herb_fwc_hyd, aes(TreatmentIntensity, AreaChange_ha)) +
  geom_point()

ggplot(plant_herb_fwc_hyd, aes(TreatmentFrequency, FreqChange)) +
  geom_point()

ggplot(plant_herb_fwc_hyd, aes(TreatmentIntensity, FreqChange)) +
  geom_point()

plant_herb_fwc_hyd %>%
  filter(!is.na(FreqChange)) %>%
  ggplot(aes(TreatmentFrequency, TreatmentIntensity)) +
  geom_point()

plant_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentFrequency * TreatmentIntensity, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm")


#### survey_herb_fwc_hyd exploratory figures ####

# hydrilla change over time
pdf("output/hydrilla_survey_herb_over_time_by_lake.pdf")
for(i in sort(unique(survey_herb_fwc_hyd$AreaOfInterestID))){
  subDat <- filter(survey_herb_fwc_hyd, AreaOfInterestID == i)
  lakeName <- unique(subDat$AreaOfInterest)
  print(ggplot(subDat, aes(SurveyDate, SpeciesFrequency_ha)) +
          geom_line() +
          geom_point(aes(color = TreatmentIntensity*TreatmentFrequency), size = 2) +
          ggtitle(lakeName) +
          scale_colour_viridis_c(name = "Herbicide\n(prop. treated/year)"))
}
dev.off()

ggplot(survey_herb_fwc_hyd, aes(TreatmentFrequency, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(survey_herb_fwc_hyd, aes(TreatmentIntensity, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(survey_herb_fwc_hyd, aes(TreatmentFrequency, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

survey_herb_fwc_hyd %>%
  filter(TreatmentFrequency > 0) %>%
  ggplot(aes(TreatmentFrequency, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(survey_herb_fwc_hyd, aes(TreatmentIntensity, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

survey_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

survey_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity * TreatmentFrequency, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

survey_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentFrequency, TreatmentIntensity)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)


#### cltv_herb_fwc_hyd exploratory figures ####

# hydrilla change over time
pdf("output/hydrilla_cumulative_herb_over_time_by_lake.pdf")
for(i in sort(unique(cltv_herb_fwc_hyd$AreaOfInterestID))){
  subDat <- filter(cltv_herb_fwc_hyd, AreaOfInterestID == i)
  lakeName <- unique(subDat$AreaOfInterest)
  print(ggplot(subDat, aes(SurveyDate, SpeciesFrequency_ha)) +
          geom_line() +
          geom_point(aes(color = TreatmentIntensity*TreatmentFrequency), size = 2) +
          ggtitle(lakeName) +
          scale_colour_viridis_c(name = "Herbicide\n(prop. treated/year)"))
}
dev.off()

ggplot(cltv_herb_fwc_hyd, aes(TreatmentFrequency, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(cltv_herb_fwc_hyd, aes(TreatmentIntensity, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(cltv_herb_fwc_hyd, aes(TreatmentFrequency, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

cltv_herb_fwc_hyd %>%
  filter(TreatmentFrequency > 0) %>%
  ggplot(aes(TreatmentFrequency, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(cltv_herb_fwc_hyd, aes(TreatmentIntensity, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

cltv_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

cltv_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity * TreatmentFrequency, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

cltv_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentFrequency, TreatmentIntensity)) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  geom_point() +
  geom_smooth(method = "lm", se = F)

cltv_herb_fwc_hyd %>%
  filter(TreatmentFrequency <= 1 & TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

#### survey_change_lw_hyd exploratory figures ####

pdf("output/hydrilla_lw_over_time_by_lake.pdf")
for(i in sort(unique(survey_change_lw_hyd$uniqueID))){
  subDat <- filter(survey_change_lw_hyd, uniqueID == i)
  lakeName <- unique(subDat$Lake)
  print(ggplot(subDat, aes(Date, SpeciesFrequency)) +
          geom_line() +
          ggtitle(lakeName))
}
dev.off()