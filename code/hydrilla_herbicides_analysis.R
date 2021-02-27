#### info ####

# goal: evaluate effects of herbicides on hydrilla


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)
library(pracma)
library(zoo)
library(ggfortify)

# stan settings
source("code/stan_settings.R")

# import data
ctrl_old <- read_csv("intermediate-data/FWC_control_old_formatted.csv")
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
plant_lw <- read_csv("intermediate-data/LW_plant_formatted.csv")

# assumptions
MinTime = 5 # minimum number of time points needed

# figure settings
def_theme <- theme_bw() +
  theme(axis.text = element_text(size = 8, color="black"),
        axis.title = element_text(size = 10, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.text = element_text(size = 10, color="black"),
        strip.background = element_blank())


#### edit FWC plant data ####

# Perm IDs per AOI
plant_fwc %>%
  group_by(AreaOfInterestID) %>%
  summarise(IDs = length(unique(PermanentID))) %>%
  filter(IDs > 1)
# one Perm ID per AOI

plant_fwc_hyd <- plant_fwc %>% # start with all surveys (no hydrilla means abundance = 0)
  filter(SpeciesName != "Hydrilla verticillata") %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate) %>%
  unique() %>% # don't need a row for each species
  left_join(plant_fwc %>% # add hydrilla information
              filter(SpeciesName == "Hydrilla verticillata") %>%
              select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, SpeciesAcres)) %>%
  unique() %>% # remove duplicate reports
  group_by(AreaOfInterestID, SurveyDate) %>%  # remove duplicate reports from combining two datasets
  mutate(NumSurveysPerDate = n()) %>%
  ungroup() %>%
  filter(!(NumSurveysPerDate > 1 & is.na(SpeciesAcres))) %>%
  mutate(SpeciesAcres = replace_na(SpeciesAcres, 0), # add columns
         Area_ha = ShapeArea * 100,
         AreaCovered_ha = SpeciesAcres * 0.405,
         AreaCovered_ha = case_when(AreaCovered_ha > Area_ha ~ Area_ha,
                                    TRUE ~ AreaCovered_ha),
         log_AreaCovered_ha = log(AreaCovered_ha + 0.001),
         SurveyYear = year(SurveyDate),
         SurveyMonth = month(SurveyDate)) %>%
  group_by(AreaOfInterestID, SurveyYear) %>% # remove reports of zero when there is another report that year
  mutate(AreaCoveredAnnAvg_ha = mean(AreaCovered_ha),
         FirstSurveyPerYear = min(SurveyDate)) %>%
  ungroup() %>%
  filter(!(AreaCoveredAnnAvg_ha > 0 & AreaCovered_ha == 0) & # other report is non-zero
           !(AreaCovered_ha == 0 & SurveyDate != FirstSurveyPerYear)) %>% # other report is zero
  group_by(AreaOfInterestID) %>%
  arrange(SurveyDate) %>%
  mutate(LastSurvey_days = SurveyDate - lag(SurveyDate),
         log_AreaCoveredChange_ha = log_AreaCovered_ha - lag(log_AreaCovered_ha),
         MostFreqMonth =  Mode(SurveyMonth)) %>%
  ungroup() %>%
  mutate(MonthDisplacement = SurveyMonth - MostFreqMonth, # months away from most frequent month
         SurveyYearAdj = case_when(MonthDisplacement > 6 ~ SurveyYear + 1, # move into next year if more than 6 months away
                                   MonthDisplacement < -6 ~ SurveyYear - 1, # move into previous year if more than 6 months away
                                   TRUE ~ SurveyYear),
         MonthDisplacementAdj = case_when(MonthDisplacement > 6 ~ 12 - SurveyMonth + MostFreqMonth,
                                          MonthDisplacement < -6 ~ 12 - MostFreqMonth + SurveyMonth,
                                          TRUE ~ MonthDisplacement)) %>%
  group_by(AreaOfInterestID, SurveyYearAdj) %>% # three cases of duplicate surveys (displacements: -2 and 3, -6 and 5, 0 and 4)
  mutate(MinDisp = min(MonthDisplacementAdj)) %>% # choose survey with lower displacement
  ungroup() %>%
  filter(MonthDisplacementAdj == MinDisp) %>%
  select(-c(NumSurveysPerDate, AreaCoveredAnnAvg_ha, FirstSurveyPerYear, MinDisp))


#### edit control data ####

# non-herbicide methods (from herbicide_initial_visualizations)
non_herb <- c("Mechanical Harvester", 
              "Snagging (tree removal)", 
              "Aquatic Dye (for shading)", 
              "Grass Carp", "Hand Removal", 
              "Mechanical (Other)", 
              "Mechanical Shredder", 
              "Prescribed Fire")

# old herbicide data
ctrl_old_hyd <- ctrl_old %>%
  filter(Species == "Hydrilla verticillata" & TotalAcres > 0) %>%
  mutate(AreaTreated_ha= TotalAcres * 0.405,
         Area_ha = ShapeArea * 100) %>%
  group_by(AreaOfInterestID, PermanentID, Year, Area_ha) %>%
  summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
  ungroup() %>%
  mutate(TreatmentMonth = 12) %>%
  rename(TreatmentYear = Year)

# new herbicide data
ctrl_new_hyd <- ctrl_new %>%
  filter(Species == "Hydrilla verticillata" & TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
  mutate(AreaTreated_ha = TotalAcres * 0.405,
         Area_ha = ShapeArea * 100) %>%
  group_by(AreaOfInterestID, PermanentID, BeginDate, Area_ha) %>%
  summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
  ungroup() %>%
  mutate(TreatmentYear = year(BeginDate),
         TreatmentMonth = month(BeginDate)) %>%
  rename(TreatmentDate = BeginDate)

# combine herbicide data
ctrl_hyd <- ctrl_old_hyd %>%
  full_join(ctrl_new_hyd) %>%
  mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha,
                                    TRUE ~ AreaTreated_ha))
  

#### combine FWC plant and ctrl data ####

# ctrl data, adjusted for plant data
ctrl_hyd_adj <- ctrl_hyd %>%
  left_join(plant_fwc_hyd %>%
              select(AreaOfInterestID, PermanentID, MostFreqMonth) %>% # add most frequent month for waterbody
              unique()) %>%
  mutate(TreatmentYearAdj = case_when(TreatmentMonth < MostFreqMonth ~ TreatmentYear - 1, # move treatment into previous year if it occurred before survey
                                      TRUE ~ TreatmentYear)) %>%
  group_by(AreaOfInterestID, PermanentID, TreatmentYearAdj) %>%
  summarise(TreatmentFrequency = n(), # treatments per year
            TreatmentIntensity = mean(AreaTreated_ha/Area_ha), # average area treated
            Treatment = TreatmentFrequency * TreatmentIntensity) %>%
  ungroup()

# add ctrl data for years with plant surveys
plant_fwc_hyd_ctrl <- plant_fwc_hyd %>%
  rename(YearAdj = SurveyYearAdj) %>%
  left_join(ctrl_hyd_adj %>%
              rename(YearAdj = TreatmentYearAdj) %>%
              mutate(YearAdj = YearAdj + 1)) %>% # add year to control data so that it's in the same row as the survey that followed it
  mutate(TreatmentFrequency = replace_na(TreatmentFrequency, 0),
         TreatmentIntensity = replace_na(TreatmentIntensity, 0),
         Treatment = replace_na(Treatment, 0))

# split data for before/after control data available
plant_fwc_hyd_pre <- plant_fwc_hyd %>%
  filter(SurveyYearAdj < (min(ctrl_hyd_adj$TreatmentYearAdj) + 1)) %>%
  group_by(AreaOfInterestID) %>%
  mutate(NumSurveys = n()) %>%
  ungroup() %>%
  filter(NumSurveys >= MinTime) # remove lakes with too few surveys

plant_fwc_hyd_post <- plant_fwc_hyd_ctrl %>%
  filter(YearAdj > min(ctrl_hyd_adj$TreatmentYearAdj)) %>%
  group_by(AreaOfInterestID) %>%
  mutate(NumSurveys = n()) %>%
  ungroup() %>%
  filter(NumSurveys >= MinTime) # remove lakes with too few surveys


#### visualizations ####

# original month displacement
ggplot(plant_fwc_hyd, aes(x = MonthDisplacement)) +
  geom_histogram(binwidth = 1)

# adjusted month displacement
pdf("output/adjusted_month_displacement_hydrilla_fwc.pdf", width = 4, height = 4)
ggplot(plant_fwc_hyd, aes(x = MonthDisplacementAdj)) +
  geom_histogram(binwidth = 1) +
  xlab("Month displacement (J)") +
  ylab("Surveys") +
  def_theme
dev.off()

# treatment
pdf("output/herbicide_distribution_hydrilla.pdf", width = 4, height = 4)
ggplot(plant_fwc_hyd_post, aes(x = Treatment)) +
  geom_histogram(binwidth = 0.01) +
  xlab("Herbicide (H, proportion lake x applications)") +
  ylab("Surveys") +
  def_theme
dev.off()

# initial treatment visualization
ggplot(plant_fwc_hyd_post, aes(x = Treatment, y = log_AreaCoveredChange_ha)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (H, proportion lake x applications)") +
  ylab("Change in hydrilla cover (log(final/initial))") +
  def_theme
# only 28 missing values because the others have surveys pre-herbicide


#### start here: make pre-herbicide model ####


#### edit data ####

# select data associated with Hydrilla
# convert all to hectares




plant_lw_hyd <- plant_lw %>%
  filter(GenusSpecies == "Hydrilla verticillata") %>%
  mutate(Area_ha = ShapeArea * 100)

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

# FWC average plant change
plant_change_fwc_hyd <- plant_fwc_hyd %>%
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
  select(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, SurveyDate, AreaCovered_ha, log_AreaCovered_ha) %>%
  arrange(AreaOfInterestID, SurveyDate) %>%
  group_by(AreaOfInterestID, PermanentID) %>%
  mutate(AreaNotCovered_ha = Area_ha - AreaCovered_ha,
         DaysChange = SurveyDate - lag(SurveyDate),
         OrigAreaCovered_ha = lag(AreaCovered_ha),
         Lag2AreaCovered_ha = lag(AreaCovered_ha, n = 2),
         AreaChange_haPerYear = (AreaCovered_ha - lag(AreaCovered_ha))/as.numeric(DaysChange) * 365,
         RelAreaChange = (AreaCovered_ha - lag(AreaCovered_ha))/lag(AreaCovered_ha),
         SpeciesFrequency_ha = AreaCovered_ha/Area_ha,
         OrigFreq = lag(SpeciesFrequency_ha),
         FreqChange = (SpeciesFrequency_ha - lag(SpeciesFrequency_ha)),
         FreqChangePerYear = (SpeciesFrequency_ha - lag(SpeciesFrequency_ha))/as.numeric(DaysChange) * 365,
         AbsFreqChangePerYear = abs(FreqChangePerYear),
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
         AreaTreated_ha = case_when(SurvTreatDays > MinHerbLag & SurvTreatDays <= (DaysChange + MinHerbLag) ~ AreaTreated_ha,
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
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, SurveyDate, AreaCovered_ha, OrigAreaCovered_ha, log_AreaCovered_ha, DaysChange, AreaChange_haPerYear, SpeciesFrequency_ha, RelAreaChange, FreqChange, FreqChangePerYear, AbsFreqChangePerYear) %>%
  summarise(TreatmentYears = case_when(sum(AreaTreated_ha) > 0 ~ unique((as.numeric(DaysChange) - MinHerbLag)/365),
                                       TRUE ~ 0),
            TreatmentFrequency = sum(AreaTreated_ha > 0)/unique((as.numeric(DaysChange) - MinHerbLag) * 365),
            TreatmentIntensity = mean(AreaTreated_ha/OrigAreaCovered_ha),
            TreatmentIntensityLake = mean(AreaTreated_ha/Area_ha)) %>%
  ungroup() %>%
  mutate(Herbicide = TreatmentFrequency * TreatmentIntensityLake)

# subset with minimum time points
survey_ts_fwc_hyd <- survey_herb_fwc_hyd %>%
  group_by(AreaOfInterestID, PermanentID) %>%
  mutate(TimePoints = n()) %>%
  ungroup() %>%
  filter(TimePoints >= 5)

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
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, SurveyDate, AreaCovered_ha, OrigAreaCovered_ha, DaysChange, AreaChange_haPerYear, SpeciesFrequency_ha, FreqChangePerYear) %>%
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
  arrange(PermanentID, Lake, Date) %>%
  group_by(PermanentID, Lake, County_LW, Area_ha) %>%
  mutate(DaysChange = Date - lag(Date),
         OrigFreq = lag(SpeciesFrequency),
         FreqChange = (SpeciesFrequency - lag(SpeciesFrequency)),
         FreqChangePerYear = (SpeciesFrequency - lag(SpeciesFrequency))/as.numeric(DaysChange) * 365,
         FreqChangePreFreqPerYear = (SpeciesFrequency - lag(SpeciesFrequency))/(lag(SpeciesFrequency) * as.numeric(DaysChange)) * 365) %>%
  ungroup() %>%
  mutate(uniqueID = paste(Lake, County_LW, PermanentID, sep = "_"))

# LW survey-to-survey herbicide and plants
# herbicide treatments that occurred between the two surveys are included
# unless it is within MinHerbLag days of the last survey
survey_herb_lw_hyd <- survey_change_lw_hyd %>%
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
  mutate(SurvTreatDays = Date - TreatmentDate,
         AreaTreated_ha = case_when(SurvTreatDays > MinHerbLag & SurvTreatDays <= DaysChange ~ AreaTreated_ha,
                                    TRUE ~ 0),
         SurvTreatDays = case_when(AreaTreated_ha == 0 ~ NA_real_,
                                   TRUE ~ as.numeric(SurvTreatDays)),
         TreatmentDate = case_when(AreaTreated_ha == 0 ~ NA_Date_,
                                   TRUE ~ TreatmentDate)) %>%
  unique() %>%
  group_by(uniqueID, Date) %>%
  mutate(rows = n()) %>%
  ungroup() %>%
  filter(!(AreaTreated_ha == 0 & rows > 1)) %>%
  group_by(uniqueID, PermanentID, Lake, County_LW, Area_ha, Date, DaysChange, SpeciesFrequency, FreqChangePerYear) %>%
  summarise(TreatmentYears = case_when(sum(AreaTreated_ha) > 0 ~ unique((as.numeric(DaysChange) - MinHerbLag)/365),
                                       TRUE ~ 0),
            TreatmentFrequency = sum(AreaTreated_ha > 0)/unique((as.numeric(DaysChange) - MinHerbLag) * 365),
            TreatmentIntensityLake = mean(AreaTreated_ha/Area_ha)) %>%
  ungroup()

  

#### plant_fwc_hyd exploratory figures ####

# small invasions over time
plant_fwc_hyd %>%
  filter(SpeciesAcres < 10) %>%
  ggplot(aes(x = SurveyDate, y = SpeciesFrequency_ha, color = as.factor(AreaOfInterestID))) +
  geom_line(alpha = 0.5) +
  theme(legend.position = "none")
# lots of low values
# non-linear changes over time

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

# see if species frequency depends on area source
ggplot(plant_fwc_hyd, aes(x = SpeciesFrequency, y = SpeciesFrequency_ha)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point()
# a lot more are greater than 1
# cluster around zero: probably total shape area includes larger waterbody connected to smaller one


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

# lagged population regulation
pdf("output/hydrilla_survey_lag_by_lake.pdf")
for(i in sort(unique(survey_change_fwc_hyd$AreaOfInterestID))){
  subDat <- filter(survey_change_fwc_hyd, AreaOfInterestID == i)
  lakeName <- unique(subDat$AreaOfInterest)
  print(ggplot(subDat, aes(OrigAreaCovered_ha, AreaCovered_ha)) +
          geom_point() +
          geom_smooth(method = lm) +
          xlab("Area covered by hydrilla at t-1 (log-ha)") +
          ylab("Area covered by hydrilla at t (log-ha)") +
          ggtitle(lakeName) +
          theme_bw())
}
dev.off()

ggplot(survey_change_fwc_hyd, aes(DaysChange, AreaChange_ha)) +
  geom_vline(xintercept = 365, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.5)
# smaller changes with longer times between surveys

survey_change_fwc_hyd %>%
  filter(!is.na(AreaChange_ha)) %>%
  ggplot(aes(DaysChange, log(AreaCovered_ha))) +
  geom_hline(yintercept = log(0.0405), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 365, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.5)
# not all small changes are 0.01 acres (common entry)
# many surveys taken 1 year apart

ggplot(survey_change_fwc_hyd, aes(x = FreqChangePerYear)) +
  geom_histogram()

ggplot(survey_change_fwc_hyd, aes(x = AreaChange_haPerYear)) +
  geom_histogram()

ggplot(survey_change_fwc_hyd, aes(log(Area_ha), FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm")
# no strong relationship between lake size and freq change

ggplot(survey_change_fwc_hyd, aes(log(Area_ha), abs(AreaChange_haPerYear))) +
  geom_point() +
  geom_smooth(method = "lm")
# the absolute value in area change increases with area (expected, constrained variable)

ggplot(survey_change_fwc_hyd, aes(OrigFreq, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm")
# frequency change is constrained by initial frequency

ggplot(survey_change_fwc_hyd, aes(OrigFreq, SpeciesFrequency_ha)) +
  geom_point() +
  geom_smooth(method = "lm")
# trend of temporal autocorrelation, but tons of variation

survey_change_fwc_hyd %>%
  mutate(OrigFreqBin = cut_interval(OrigFreq, n = 10)) %>%
  filter(!is.na(OrigFreqBin))  %>%
  ggplot(aes(x = FreqChangePerYear)) +
  geom_histogram() +
  facet_wrap(~ OrigFreqBin, scales = "free_y")
# the distribution of frequency changes shifts from positive to negative with increasing original frequency
# this is partially because of the constraint noted above

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
          geom_point(aes(color = Herbicide), size = 2) +
          xlab("Date") +
          ylab(expression(paste(italic("Hydrilla"), " abundance (prop. lake)", sep = ""))) +
          ggtitle(lakeName) +
          scale_colour_viridis_c(name = "Herbicide\n(prop. lake\ntreated/year)") +
          coord_cartesian(ylim = c(0, 1)) +
          theme_bw())
}
dev.off()

ggplot(survey_herb_fwc_hyd, aes(x = Herbicide)) +
  geom_histogram()

ggplot(survey_herb_fwc_hyd, aes(x = TreatmentFrequency)) +
  geom_histogram()

ggplot(survey_herb_fwc_hyd, aes(x = TreatmentIntensityLake)) +
  geom_histogram()

ggplot(survey_herb_fwc_hyd, aes(Herbicide, AbsFreqChangePerYear, color = as.factor(AreaOfInterestID))) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme(legend.position = "none")

ggplot(survey_herb_fwc_hyd, aes(TreatmentFrequency, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(survey_herb_fwc_hyd, aes(TreatmentIntensity, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(survey_herb_fwc_hyd, aes(TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

survey_herb_fwc_hyd %>%
  filter(TreatmentFrequency > 0) %>%
  ggplot(aes(TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(survey_herb_fwc_hyd, aes(TreatmentIntensity, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

survey_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

pdf("output/hydrilla_change_survey_herbicide.pdf")
survey_herb_fwc_hyd %>%
  #filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()

survey_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()
dev.off()

survey_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentFrequency, TreatmentIntensity)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)


#### survey_ts_fwc_hyd exploratory figures ####

# longest time series?
max_ts <- survey_ts_fwc_hyd %>%
  filter(TimePoints == max(TimePoints))
# max = 38

pdf("output/longest_time_series_lakes.pdf")
survey_ts_fwc_hyd %>%
  filter(AreaOfInterestID %in% max_ts$AreaOfInterestID) %>%
  mutate(Lake_size = paste(str_replace(AreaOfInterest, ", Lake", ""), " (", round(Area_ha), " ha)", sep = "")) %>%
  ggplot(aes(SurveyDate, log_AreaCovered_ha)) +
  geom_vline(xintercept = as.Date("1998-12-31", "%Y-%m-%d"), linetype = "dashed") +
  geom_line() +
  geom_point(aes(color = Herbicide)) +
  facet_wrap(~ Lake_size, scales = "free_x") +
  scale_color_viridis_c(name = "Herbicide\n(prop. lake\ntreated/year)") +
  xlab("Date") +
  ylab("Area covered by hydrilla (log-ha)") +
  def_theme
dev.off()

# minimum time interval
min(survey_ts_fwc_hyd$DaysChange, na.rm = T)
# 162 days
162/30 # ~5.5 months

# time interval distribution
ggplot(survey_ts_fwc_hyd, aes(x = DaysChange)) +
  geom_histogram()

# month samples
survey_ts_fwc_hyd %>%
  mutate(Month = month(SurveyDate)) %>%
  ggplot(aes(x = Month)) +
  geom_histogram(bins = 12)
# clustered between May and October

# samples per year
survey_ts_fwc_hyd %>%
  mutate(Year = year(SurveyDate)) %>%
  group_by(AreaOfInterestID, PermanentID, Year) %>%
  count() %>%
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 1)
# most are once per year

survey_ts_fwc_hyd %>%
  mutate(Year = year(SurveyDate)) %>%
  group_by(AreaOfInterestID, PermanentID, Year) %>%
  count() %>%
  filter(n > 1)
# one was sampled twice in one year - maybe use average


#### cltv_herb_fwc_hyd exploratory figures ####

# hydrilla change over time
pdf("output/hydrilla_cumulative_herb_over_time_by_lake.pdf")
for(i in sort(unique(cltv_herb_fwc_hyd$AreaOfInterestID))){
  subDat <- filter(cltv_herb_fwc_hyd, AreaOfInterestID == i)
  lakeName <- unique(subDat$AreaOfInterest)
  print(ggplot(subDat, aes(SurveyDate, SpeciesFrequency_ha)) +
          geom_line() +
          geom_point(aes(color = TreatmentIntensityLake*TreatmentFrequency), size = 2) +
          xlab("Date") +
          ylab(expression(paste(italic("Hydrilla"), " abundance (prop. lake)", sep = ""))) +
          ggtitle(lakeName) +
          scale_colour_viridis_c(name = "Herbicide\n(prop. lake\ntreated/year)") +
          coord_cartesian(ylim = c(0, 1)) +
          theme_bw())
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
  ggplot(aes(TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(cltv_herb_fwc_hyd, aes(TreatmentIntensity, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

cltv_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

pdf("output/hydrilla_change_cumulative_herbicide.pdf")
cltv_herb_fwc_hyd %>%
  #filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()

cltv_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()
dev.off()

cltv_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentFrequency, TreatmentIntensity)) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  geom_point() +
  geom_smooth(method = "lm", se = F)

cltv_herb_fwc_hyd %>%
  filter(TreatmentFrequency <= 1 & TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

#### survey_herb_lw_hyd exploratory figures ####

# hydrilla change over time
pdf("output/hydrilla_lw_survey_herb_over_time_by_lake.pdf")
for(i in sort(unique(survey_herb_lw_hyd$uniqueID))){
  subDat <- filter(survey_herb_lw_hyd, uniqueID == i)
  lakeName <- unique(subDat$Lake)
  print(ggplot(subDat, aes(Date, SpeciesFrequency)) +
          geom_line() +
          geom_point(aes(color = TreatmentIntensityLake*TreatmentFrequency), size = 2) +
          xlab("Date") +
          ylab(expression(paste(italic("Hydrilla"), " abundance (prop. lake)", sep = ""))) +
          ggtitle(lakeName) +
          scale_colour_viridis_c(name = "Herbicide\n(prop. lake\ntreated/year)") +
          coord_cartesian(ylim = c(0, 1)) +
          theme_bw())
}
dev.off()

pdf("output/hydrilla_change_lw_survey_herbicide.pdf")
survey_herb_lw_hyd %>%
  #filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()

survey_herb_lw_hyd %>%
  filter(TreatmentIntensityLake > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()
dev.off()


#### test model ####

# longest time series?
max_ts <- plant_fwc_hyd %>%
  group_by(AreaOfInterestID) %>%
  count() %>%
  ungroup() %>%
  filter(n == max(n))

# visualize
plant_fwc_hyd %>%
  filter(AreaOfInterestID %in% max_ts$AreaOfInterestID) %>%
  ggplot(aes(SurveyDate, SpeciesFrequency_ha)) + 
  geom_line() +
  facet_wrap(~AreaOfInterestID)

# test time series
test_dat <- filter(plant_fwc_hyd, AreaOfInterestID == 206) %>%
  mutate(MinDate = min(SurveyDate),
         DateNum = as.numeric(SurveyDate - MinDate))
ggplot(test_dat, aes(DateNum, SpeciesFrequency_ha)) + geom_line()
ggplot(test_dat, aes(DateNum, log_AreaCovered_ha)) + geom_line()
test_ts <- zoo(test_dat$log_AreaCovered_ha, test_dat$DateNum)
plot(test_ts)

# stan data
stan_data <- within(list(), {
  y <- as.vector(test_ts)
  n <- length(test_ts)
})

# herbicide time series
test_dat2 <- filter(survey_ts_fwc_hyd, AreaOfInterestID == 206) %>%
  mutate(MinDate = min(SurveyDate),
         DateNum = as.numeric(SurveyDate - MinDate)) %>%
  filter(!is.na(FreqChangePerYear))
ggplot(test_dat2, aes(DateNum, FreqChangePerYear)) + geom_line()
ggplot(test_dat2, aes(DateNum, Herbicide)) + geom_line()
ggplot(test_dat2, aes(DateNum, log_AreaCovered_ha)) + 
  geom_line() +
  geom_point(aes(color = Herbicide)) +
  scale_color_viridis_c()
ggplot(test_dat2, aes(DateNum, abs(FreqChangePerYear))) + geom_line()
ggplot(test_dat2, aes(Herbicide, abs(FreqChangePerYear))) + geom_point()
ggplot(test_dat2, aes(Herbicide, log_AreaCovered_ha)) + geom_point()
test_ts2 <- zoo(test_dat2$log_AreaCovered_ha, test_dat2$DateNum)
test_x2 <- zoo(test_dat2$Herbicide, test_dat2$DateNum)

# stan data
stan_data2 <- within(list(), {
  y <- as.vector(test_ts2)
  x <- as.vector(test_x2)
  n <- length(test_ts2)
})

# simulated data
test_dat3 <- tibble(days = 0:30,
                    yhat = 10,
                    y = yhat + rnorm(1, 0, 1))

for(i in 2:nrow(test_dat3)){
  test_dat3$yhat[i] <- test_dat3$yhat[i-1] * exp(1.1 - 0.055 * test_dat3$yhat[i-1])
  test_dat3$y[i] = test_dat3$yhat[i] + rnorm(1, 0, 1)
}

ggplot(test_dat3, aes(x = days, y = yhat)) +
  geom_line() +
  geom_point() + 
  geom_line(aes(y = y))

test_ts3 <- zoo(test_dat3$y, test_dat3$days)

# stan data
stan_data3 <- within(list(), {
  y <- as.vector(test_ts3)
  n <- length(test_ts3)
  mu0 <- log(test_dat3$y[1])
})


# model 1: level only, observation error only

# model file
model_file1 <- 'models/test_level_only_obs_error_only.stan'
cat(paste(readLines(model_file1)), sep = '\n')

# fit model
test_fit1 <- stan(file = model_file1, data = stan_data,
            iter = 2000, chains = 4)

mu1 <- get_posterior_mean(test_fit1, par = 'mu')[, 'mean-all chains']

# figures
autoplot(test_ts) + 
  geom_hline(yintercept = mu1, linetype = "dashed") +
  ggtitle("level only, observation error only")

autoplot(test_ts - mu1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggtitle("observation error")

# model 2: level only, observation and process error

# model file
model_file2 <- 'models/test_level_only_obs_error_proc_error.stan'
cat(paste(readLines(model_file2)), sep = '\n')

# fit model
test_fit2 <- stan(file = model_file2, data = stan_data,
                 iter = 2000, chains = 4)

mu2 <- get_posterior_mean(test_fit2, par = 'mu')[, 'mean-all chains']
plot_dat2 <- tibble(model = mu2, 
               Days = test_dat$DateNum,
               observation = test_dat$log_AreaCovered_ha) %>%
  mutate(ObsError = observation - model,
         ProcError = model - lag(model)) %>%
  pivot_longer(-Days, names_to = "yType", values_to = "y")

# figures
plot_dat2 %>%
  filter(!(yType %in% c("ObsError", "ProcError"))) %>%
  ggplot(aes(Days, y, color = yType)) +
  geom_line() +
  ylab("log(Area Covered (ha))") +
  theme_bw() +
  theme(legend.title = element_blank())

plot_dat2 %>%
  filter(yType == "ObsError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Observation error") +
  theme_bw()

plot_dat2 %>%
  filter(yType == "ProcError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Process error") +
  theme_bw()

# model 3: level and slope, observation and process error

# model file
model_file3 <- 'models/test_level_slope_obs_error_proc_error.stan'
cat(paste(readLines(model_file3)), sep = '\n')

# fit model
test_fit3 <- stan(file = model_file3, data = stan_data,
                  iter = 2000, chains = 4)

# extract values
mu3 <- get_posterior_mean(test_fit3, par = 'mu')[, 'mean-all chains']
v3 <- get_posterior_mean(test_fit3, par = 'v')[, 'mean-all chains']

plot_dat3 <- tibble(model = mu3, 
                    slope = c(v3, NA),
                    Days = test_dat$DateNum,
                    observation = test_dat$log_AreaCovered_ha) %>%
  mutate(ObsError = observation - model,
         SlopeProcError = slope - lag(slope),
         LevelProcError = model - lag(model) - lag(slope)) %>%
  pivot_longer(-Days, names_to = "yType", values_to = "y")

# figures
plot_dat3 %>%
  filter(yType %in% c("model", "observation")) %>%
  ggplot(aes(Days, y, color = yType)) +
  geom_line() +
  ylab("log(Area Covered (ha))") +
  theme_bw() +
  theme(legend.title = element_blank())

plot_dat3 %>%
  filter(yType == "ObsError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Observation error") +
  theme_bw()

plot_dat3 %>%
  filter(yType %in% c("SlopeProcError", "LevelProcError")) %>%
  ggplot(aes(Days, y, color = yType)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Process error") +
  theme_bw()

# model 4: deterministic model explanatory variable (regression)

# model file
model_file4 <- 'models/test_level_explanatory_obs_error_only.stan'
cat(paste(readLines(model_file4)), sep = '\n')

# fit model
test_fit4 <- stan(file = model_file4, data = stan_data2,
                  iter = 2000, chains = 4)

# extract values
yhat4 <- get_posterior_mean(test_fit4, par = 'yhat')[, 'mean-all chains']

plot_dat4 <- tibble(model = yhat4, 
                    Days = test_dat2$DateNum,
                    Herbicide = test_dat2$Herbicide,
                    observation = test_dat2$FreqChangePerYear) %>%
  mutate(ObsError = observation - model) %>%
  pivot_longer(-c(Days, Herbicide), names_to = "yType", values_to = "y")

# figures
plot_dat4 %>%
  filter(yType %in% c("model", "observation")) %>%
  ggplot(aes(Days, y, color = yType)) +
  geom_line() +
  ylab("Frequency change per year") +
  theme_bw() +
  theme(legend.title = element_blank())

plot_dat4 %>%
  filter(yType == "ObsError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Observation error") +
  theme_bw()


# model 5: stochastic model with explanatory variable

# model file
model_file5 <- 'models/test_level_explanatory_obs_error_proc_error.stan'
cat(paste(readLines(model_file5)), sep = '\n')

# fit model
test_fit5 <- stan(file = model_file5, data = stan_data2,
                  iter = 2000, chains = 4)

# extract values
yhat5 <- get_posterior_mean(test_fit5, par = 'yhat')[, 'mean-all chains']
mu5 <- get_posterior_mean(test_fit5, par = 'mu')[, 'mean-all chains']

plot_dat5 <- tibble(mu = mu5,
                    model = yhat5, 
                    Days = test_dat2$DateNum,
                    Herbicide = test_dat2$Herbicide,
                    observation = test_dat2$FreqChangePerYear) %>%
  mutate(ObsError = observation - model,
         ProcError = mu - lag(mu)) %>%
  pivot_longer(-c(Days, Herbicide), names_to = "yType", values_to = "y")

# figures
plot_dat5 %>%
  filter(yType %in% c("model", "observation")) %>%
  ggplot(aes(Days, y, color = yType)) +
  geom_line() +
  ylab("Frequency change per year") +
  theme_bw() +
  theme(legend.title = element_blank())

plot_dat5 %>%
  filter(yType == "ObsError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Observation error") +
  theme_bw()

plot_dat5 %>%
  filter(yType == "ProcError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Process error") +
  theme_bw()

# model 6: Ricker model, observation error only

# model file
model_file6 <- 'models/test_ricker_obs_error_only.stan'
cat(paste(readLines(model_file6)), sep = '\n')

# fit model
test_fit6 <- stan(file = model_file6, data = stan_data3,
                  iter = 2000, chains = 4)

print(test_fit6, max = 50)
# generally bad at estimating parameters, could give stricter priors...

# try simple model
test_fit6b <- stan(file = model_file1, data = stan_data3,
                  iter = 2000, chains = 4)

mu6b <- get_posterior_mean(test_fit6b, par = 'mu')[, 'mean-all chains']

# figures
autoplot(test_ts3) + 
  geom_hline(yintercept = mu6b, linetype = "dashed") +
  ggtitle("level only, observation error only")

autoplot(test_ts3 - mu6b) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggtitle("observation error")
