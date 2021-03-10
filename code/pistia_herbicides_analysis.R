#### info ####

# goal: evaluate effects of herbicides on pistia


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)
library(pracma) # for Mode
library(reshape2) # for melt

# stan settings
source("code/stan_settings.R")

# import data
ctrl_old <- read_csv("intermediate-data/FWC_control_old_formatted.csv")
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
plant_lw <- read_csv("intermediate-data/LW_plant_formatted.csv")

# assumptions
MinTime = 5 # minimum number of time points needed per waterbody
SurveyDist = 30 # number of days between comparable LW and FWC surveys

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

plant_fwc_pis <- plant_fwc %>% # start with all surveys (no pistia means abundance = 0)
  filter(SpeciesName != "Pistia stratiotes") %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate) %>%
  unique() %>% # don't need a row for each species
  left_join(plant_fwc %>% # add pistia information
              filter(SpeciesName == "Pistia stratiotes") %>%
              select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, SpeciesAcres)) %>%
  unique() %>% # remove duplicate reports
  group_by(AreaOfInterestID, SurveyDate) %>%
  mutate(NumSurveysPerDate = n()) %>%
  ungroup() %>%
  filter(!(NumSurveysPerDate > 1 & is.na(SpeciesAcres))) %>%  # remove duplicate reports from combining two datasets
  mutate(SpeciesAcres = replace_na(SpeciesAcres, 0), # pistia cover 0 when it wasn't in a survey
         Area_ha = ShapeArea * 100, # convert lake area from km-squared to hectares
         AreaCovered_ha = SpeciesAcres * 0.405, # convert plant cover from acres to hectares
         AreaCovered_ha = case_when(AreaCovered_ha > Area_ha ~ Area_ha,
                                    TRUE ~ AreaCovered_ha), # make plant cover the size of the lake area if it exceeds it
         log_AreaCovered_ha = log(AreaCovered_ha + 0.001), # log-transform (min AreaCovered_ha = 0.00405)
         SurveyYear = year(SurveyDate),
         SurveyMonth = month(SurveyDate)) %>%
  group_by(AreaOfInterestID, SurveyYear) %>% # remove reports of zero when there is another report that year
  mutate(AreaCoveredAnnAvg_ha = mean(AreaCovered_ha),
         FirstSurveyPerYear = min(SurveyDate)) %>%
  ungroup() %>%
  filter(!(AreaCoveredAnnAvg_ha > 0 & AreaCovered_ha == 0) & # other report is non-zero
           !(AreaCovered_ha == 0 & SurveyDate != FirstSurveyPerYear)) %>% # other report is zero
  group_by(AreaOfInterestID) %>%
  mutate(MostFreqMonth =  Mode(SurveyMonth)) %>% # most frequently sampled month
  ungroup() %>%
  mutate(MonthDisplacement = SurveyMonth - MostFreqMonth, # months away from most frequent month
         SurveyYearAdj = case_when(MonthDisplacement > 6 ~ SurveyYear + 1, # move into next year if more than 6 months away
                                   MonthDisplacement < -6 ~ SurveyYear - 1, # move into previous year if more than 6 months away
                                   TRUE ~ SurveyYear),
         MonthDisplacementAdj = case_when(MonthDisplacement > 6 ~ 12 - SurveyMonth + MostFreqMonth, # adjust month displacement in new year
                                          MonthDisplacement < -6 ~ 12 - MostFreqMonth + SurveyMonth,
                                          TRUE ~ MonthDisplacement)) %>%
  group_by(AreaOfInterestID, SurveyYearAdj) %>% # three cases of duplicate surveys (displacements: -2 and 3, -6 and 5, 0 and 4)
  mutate(MinDisp = min(MonthDisplacementAdj)) %>% # choose survey with lower displacement
  ungroup() %>%
  filter(MonthDisplacementAdj == MinDisp) %>%
  select(-c(NumSurveysPerDate, AreaCoveredAnnAvg_ha, FirstSurveyPerYear, MinDisp))


#### pre-herbicide data FWC surveys ####

# select years
plant_fwc_pis_pre <- plant_fwc_pis %>%
  filter(SurveyYearAdj <= (min(ctrl_old$Year))) %>% # goes up to year of first application
  group_by(AreaOfInterestID) %>%
  mutate(NumSurveys = n()) %>%
  ungroup() %>%
  filter(NumSurveys >= MinTime) # remove lakes with too few surveys

# missing data?
plant_fwc_pis %>%
  group_by(AreaOfInterestID) %>%
  summarise(YearRange = max(SurveyYearAdj) - min(SurveyYearAdj) + 1,
            SurveyYears = n()) %>%
  ungroup() %>%
  filter(YearRange != SurveyYears)
# 229 lakes

# lowest cover
min(plant_fwc_pis$log_AreaCovered_ha)
# -6.9

# add -99 for missing years
plant_fwc_pis_pre2 <- plant_fwc_pis_pre %>%
  group_by(AreaOfInterestID, PermanentID, AreaOfInterest, Area_ha) %>%
  summarise(MinYear = min(SurveyYearAdj),
            MaxYear = max(SurveyYearAdj)) %>% # year range for a waterbody
  ungroup() %>%
  mutate(SurveyYearAdj = map2(MinYear, MaxYear, seq, by = 1)) %>%
  unnest(cols = SurveyYearAdj) %>% # populate all years for a waterbody
  select(-c(MinYear, MaxYear)) %>%
  full_join(plant_fwc_pis_pre) %>% # missing years will have NA for area metrics
  group_by(AreaOfInterestID) %>%
  arrange(SurveyYearAdj) %>%
  mutate(log_AreaCoveredChange_ha = log_AreaCovered_ha - lag(log_AreaCovered_ha), # change in cover since last survey (NA if missing data)
         log_AreaCoveredMis_ha = replace_na(log_AreaCovered_ha, -99),
         MonthDisplacementMis = replace_na(abs(MonthDisplacementAdj), 0), # make month displacement positive
         SigmaJ = case_when(MonthDisplacementMis != 0 ~ log(1 + MonthDisplacementMis * 0.1), # 10% error by month
                            MonthDisplacementMis == 0 ~ 0.001)) %>% # small error if in survey month
  ungroup() %>%
  arrange(AreaOfInterestID, SurveyYearAdj)


#### post-herbicide data FWC surveys ####

# select data
plant_fwc_pis_post <-  plant_fwc_pis %>%
  filter(SurveyYearAdj > min(ctrl_old$Year)) %>% # starts with year after first application
  group_by(AreaOfInterestID) %>%
  mutate(NumSurveys = n()) %>%
  ungroup() %>%
  filter(NumSurveys >= MinTime) # remove lakes with too few surveys

# missing data?
plant_fwc_pis_post %>%
  group_by(AreaOfInterestID) %>%
  summarise(YearRange = max(SurveyYearAdj) - min(SurveyYearAdj) + 1,
            SurveyYears = n()) %>%
  ungroup() %>%
  filter(YearRange != SurveyYears)
# 109 lakes

# add -99 for missing years
plant_fwc_pis_post2 <- plant_fwc_pis_post %>%
  group_by(AreaOfInterestID, PermanentID, AreaOfInterest, Area_ha) %>%
  summarise(MinYear = min(SurveyYearAdj),
            MaxYear = max(SurveyYearAdj)) %>% # year range for a waterbody
  ungroup() %>%
  mutate(SurveyYearAdj = map2(MinYear, MaxYear, seq, by = 1)) %>%
  unnest(cols = SurveyYearAdj) %>% # populate all years for a waterbody
  select(-c(MinYear, MaxYear)) %>%
  full_join(plant_fwc_pis_post) %>% # missing years will have NA for area metrics
  group_by(AreaOfInterestID) %>%
  arrange(SurveyYearAdj) %>%
  mutate(log_AreaCoveredChange_ha = log_AreaCovered_ha - lag(log_AreaCovered_ha), # change in cover since last survey (NA if missing data)
         log_AreaCoveredMis_ha = replace_na(log_AreaCovered_ha, -99),
         MonthDisplacementMis = replace_na(abs(MonthDisplacementAdj), 0), # make month displacement positive
         SigmaJ = case_when(MonthDisplacementMis != 0 ~ log(1 + MonthDisplacementMis * 0.1), # 10% error by month
                            MonthDisplacementMis == 0 ~ 0.001)) %>% # small error if in survey month
  ungroup() %>%
  arrange(AreaOfInterestID, SurveyYearAdj)


#### edit control data ####

# non-herbicide methods (from herbicide_initial_visualizations)
non_herb <- c("Mechanical Harvester", 
              "Snagging (tree removal)", 
              "Aquatic Dye (for shading)", 
              "Grass Carp", "Hand Removal", 
              "Mechanical (Other)", 
              "Mechanical Shredder", 
              "Prescribed Fire")

# duplicates per year/date?
ctrl_old %>%
  filter(str_detect(Species, "Pistia") == T & TotalAcres > 0) %>%
  group_by(AreaOfInterestID, Year) %>%
  summarise(apps = n()) %>%
  filter(apps > 1)
# yes, 72

ctrl_new %>%
  filter(str_detect(Species, "Pistia") == T & TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
  group_by(AreaOfInterestID, BeginDate) %>%
  summarise(apps = n()) %>%
  filter(apps > 1)
# yes, 4928
# for treatments involving more than one herbicide, separate lines for each herbicide type, but they share a treatment ID and acres

# old herbicide data
ctrl_old_pis <- ctrl_old %>%
  filter(str_detect(Species, "Pistia") == T & TotalAcres > 0) %>%
  mutate(AreaTreated_ha= TotalAcres * 0.405,
         Area_ha = ShapeArea * 100,
         TreatmentMonth = 12) %>% # no date given for these surveys
  select(AreaOfInterestID, PermanentID, Year, Area_ha, AreaTreated_ha, TreatmentMonth) %>%
  rename(TreatmentYear = Year)

# new herbicide data
ctrl_new_pis <- ctrl_new %>%
  filter(str_detect(Species, "Pistia") == T & TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>% # herbicide control only
  select(AreaOfInterestID, PermanentID, BeginDate, TreatmentID, TotalAcres, ShapeArea) %>%
  unique() %>%
  mutate(AreaTreated_ha = TotalAcres * 0.405,
         Area_ha = ShapeArea * 100,
         TreatmentYear = year(BeginDate),
         TreatmentMonth = month(BeginDate)) %>%
  rename(TreatmentDate = BeginDate)

# combine herbicide data
ctrl_pis <- ctrl_old_pis %>%
  full_join(ctrl_new_pis) %>%
  mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha, # make treatment area size of lake if it exceeds it
                                    TRUE ~ AreaTreated_ha))

# add zero for missing years (assume no application)
ctrl_pis2 <- ctrl_old %>%
  select(AreaOfInterestID, PermanentID) %>%
  unique() %>%
  full_join(ctrl_new %>%
              select(AreaOfInterestID, PermanentID) %>%
              unique()) %>%
  expand_grid(tibble(TreatmentYear = min(ctrl_old_pis$TreatmentYear):max(ctrl_new_pis$TreatmentYear))) %>% # one row per year
  full_join(ctrl_pis) %>%
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0)) # make area treated zero if the year wasn't included


#### combine FWC plant and ctrl data ####

# ctrl data, adjusted for plant data timing
ctrl_pis_adj <- ctrl_pis %>% # use data without zeros
  left_join(plant_fwc_pis %>%
              select(AreaOfInterestID, PermanentID, MostFreqMonth) %>% # add most frequent month for waterbody
              unique()) %>%
  mutate(TreatmentYearAdj = case_when(TreatmentMonth < MostFreqMonth ~ TreatmentYear - 1, # move treatment into previous year if it occurred before survey
                                      TRUE ~ TreatmentYear)) %>%
  group_by(AreaOfInterestID, PermanentID, Area_ha, TreatmentYearAdj) %>%
  summarise(TreatmentFrequency = n(), # treatments per year
            TreatmentIntensity = mean(AreaTreated_ha/Area_ha), # average area treated
            Treatment = TreatmentFrequency * TreatmentIntensity) %>%
  ungroup() %>%
  full_join(ctrl_old %>% # all waterbodies from the dataset
              select(AreaOfInterestID, PermanentID, ShapeArea) %>%
              unique() %>%
              mutate(Area_ha = ShapeArea * 100) %>%
              select(-ShapeArea) %>%
              full_join(ctrl_new %>%
                          select(AreaOfInterestID, PermanentID, ShapeArea) %>%
                          unique() %>%
                          mutate(Area_ha = ShapeArea * 100) %>%
                          select(-ShapeArea)) %>%
              expand_grid(tibble(TreatmentYearAdj = min(ctrl_old_pis$TreatmentYear):max(ctrl_new_pis$TreatmentYear)))) %>%
  mutate(TreatmentFrequency = replace_na(TreatmentFrequency, 0), # make values 0 if no herbicide data available
         TreatmentIntensity = replace_na(TreatmentIntensity, 0),
         Treatment = replace_na(Treatment, 0)) %>%
  filter(TreatmentYearAdj < max(plant_fwc_pis_post2$SurveyYearAdj)) # remove data the year of or after plant surveys

# post herbicide surveys with control data
plant_fwc_ctrl_pis_post <-  plant_fwc_pis_post2 %>% # plant data has NA for missing values
  rename(YearAdj = SurveyYearAdj) %>%
  left_join(ctrl_pis_adj %>%
              rename(YearAdj = TreatmentYearAdj) %>%
              mutate(YearAdj = YearAdj + 1)) %>% # add year to control data so that it's in the same row as the survey that followed it
  mutate(TreatmentFrequency = replace_na(TreatmentFrequency, 0), # make values 0 if no control data exists for a waterbody (revisit)
         TreatmentIntensity = replace_na(TreatmentIntensity, 0),
         Treatment = replace_na(Treatment, 0))


#### edit LW plant data ####

# select pistia data
# create date range to merge with FWC
plant_lw_pis <- plant_lw %>%
  filter(GenusSpecies == "Pistia stratiotes") %>% # add pistia information
  select(County_LW, Lake, PermanentID, Date, SpeciesFrequency) %>%
  mutate(SurveyMin = Date - SurveyDist,
         SurveyMax = Date + SurveyDist) %>%
  rename(SurveyDate_LW = Date,
         PropCover_LW = SpeciesFrequency)


#### combine FWC and LW plant data ####

# merge based on permanentID
# merge if fwc survey is within SurveyDist days of lw survey
plant_lw_fwc_pis <- plant_fwc_pis %>%
  rename(SurveyDate_FWC = SurveyDate) %>%
  inner_join(plant_lw_pis) %>%
  filter(SurveyDate_FWC >= SurveyMin & SurveyDate_FWC <= SurveyMax) %>%
  mutate(PropCover_FWC = AreaCovered_ha / Area_ha,
         SurveyDiff = SurveyDate_LW - SurveyDate_FWC)

