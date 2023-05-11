#### set-up ####

# note: previous version: invasive_plant_data_processing.R (has FWRI processing and surveyor data too)

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)
library(magrittr)
library(lubridate)

# import data
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")

# Lake O data and model values
source("code/data-processing/okeechobee_growth.R")


#### check data ####

# consistent waterbody area?
plant_fwc %>%
  distinct(AreaOfInterestID, WaterbodyAcres) %>%
  get_dupes(AreaOfInterestID)
# should return zero

# consistent waterbody name?
plant_fwc %>%
  distinct(AreaOfInterestID, AreaOfInterest) %>%
  get_dupes(AreaOfInterestID)
# 465

# duplicates?
get_dupes(plant_fwc)
# none


#### edit data ####

# specify taxa of interest
inv_taxa <- tibble(TaxonName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes", "Panicum repens", "Urochloa mutica", "Cyperus blepharoleptos"),
                   CommonName = c("Hydrilla", "Water lettuce", "Water hyacinth", "Torpedograss", "Para grass", "Cuban bulrush"),
                   Code = c("HYDR", "WALE", "WAHY", "TORP", "PAGR", "BUSE"))

# multiple surveys in one year?
(dup_inv <- plant_fwc %>%
  filter(TaxonName %in% inv_taxa$TaxonName) %>%
  group_by(TaxonName, AreaOfInterestID, GSYear) %>%
  mutate(surveys = length(unique(SurveyDate))) %>%
  ungroup() %>%
  filter(surveys > 1) %>%
  select(TaxonName, AreaOfInterest, GSYear, SurveyDate, SpeciesAcres) %>%
  arrange(TaxonName, AreaOfInterest, GSYear, SurveyDate))
# 428

# plant abundance dataset
inv_fwc <- plant_fwc %>%
  filter(!(AreaOfInterestID == 476 & SurveyYear == 2017) &  # incomplete survey
           TaxonName %in% inv_taxa$TaxonName) %>% # select desired taxa
  mutate(AreaOfInterest = if_else(AreaOfInterest == "Watermelon Pond" & AreaOfInterestID == 465,
                                  "Watermellon Pond", # correct dual-spelled name
                                  AreaOfInterest),
         Area_ha = ShapeArea * 100, # convert PermID area from km-squared to hectares
         Waterbody_ha = WaterbodyAcres * 0.405, # convert FWC waterbody size to hectares
         AreaCovered_ha = SpeciesAcres * 0.405, # convert plant cover from acres to hectares,
         MonthDay = case_when(SurveyMonth >= 4 ~ as.Date(paste("2020", SurveyMonth, SurveyDay, sep = "-")), # start "year" in April (2020/2021 are arbitrary)
                              SurveyMonth < 4 ~ as.Date(paste("2021", SurveyMonth, SurveyDay, sep = "-")))) %>% # this is for joining dayDat
  left_join(dayDat) %>% # add standardized days (proportion between April 1 and March 31)
  mutate(AreaChangeSD = lakeO_area_beta1 * (lakeO_area_days-Days) + lakeO_area_beta2 * (lakeO_area_days^2 - Days^2)) %>% # calculate the number of sd's to change to get est. max abundance
  group_by(AreaOfInterestID, TaxonName) %>% # take standard deviation by survey area and species
  mutate(SDAreaCovered = sd(AreaCovered_ha, na.rm = T)) %>%
  ungroup() %>%
  mutate(EstAreaCoveredRaw_ha = case_when(SpeciesAcres > 0.01 ~ AreaCovered_ha + AreaChangeSD * SDAreaCovered, # adjust by sd
                                          TRUE ~ AreaCovered_ha),
         EstAreaCoveredAdj_ha = if_else(EstAreaCoveredRaw_ha > Waterbody_ha, Waterbody_ha, EstAreaCoveredRaw_ha)) %>%
  group_by(AreaOfInterestID, GSYear, TaxonName) %>% # account for multiple surveys per GSyear
  mutate(MaxSpeciesArea = if_else(sum(!is.na(EstAreaCoveredAdj_ha)) > 0,
                                  max(EstAreaCoveredAdj_ha, na.rm = T),
                                  NA_real_),
         TempSurveyDate = if_else(EstAreaCoveredAdj_ha == MaxSpeciesArea |
                                   (is.na(EstAreaCoveredAdj_ha) & is.na(MaxSpeciesArea)),
                                 SurveyDate,
                                 NA_Date_),
         MaxSurveyDate = max(TempSurveyDate, na.rm = T),
         nMaxArea = if_else(!is.na(MaxSpeciesArea),
                             sum(EstAreaCoveredAdj_ha == MaxSpeciesArea),
                            sum(is.na(EstAreaCoveredAdj_ha)))) %>%
  ungroup() %>%
  filter((nMaxArea == 1 & EstAreaCoveredAdj_ha == MaxSpeciesArea) | # choose max area or max date
           (nMaxArea > 1 & SurveyDate == MaxSurveyDate)) %>%
  select(AreaOfInterest, AreaOfInterestID, County, PermanentID, Area_ha, Waterbody_ha, 
         SpeciesAcres, AreaCovered_ha, EstAreaCoveredRaw_ha, EstAreaCoveredAdj_ha,
         SurveyDate, Surveyor, SurveyMonth, SurveyDay, SurveyYear, GSYear,
         TaxonName) %>%
  left_join(inv_taxa) %>% # add common name
  group_by(AreaOfInterestID, TaxonName) %>%
  arrange(GSYear) %>%
  mutate(GSYearDiff = GSYear - lag(GSYear),
         PercCovered = 100 * EstAreaCoveredAdj_ha / Waterbody_ha,
         PrevPercCovered = if_else(GSYearDiff == 1, lag(PercCovered), NA_real_),
         PercCovDiff = PercCovered - PrevPercCovered,
         PercCovRatio = (PercCovered + 1) / (PrevPercCovered + 1),
         PercCovLogRatio = log(PercCovRatio)) %>%
  ungroup()

# check for raw percent > 100
inv_fwc %>%
  filter(EstAreaCoveredRaw_ha > Waterbody_ha) %>%
  distinct(TaxonName, AreaOfInterest, AreaCovered_ha, EstAreaCoveredRaw_ha, Waterbody_ha) %>%
  mutate(EstDiff = EstAreaCoveredRaw_ha - Waterbody_ha) %>%
  data.frame()
# all hydrilla, pretty close to total waterbody

# missing data?
inv_fwc %>%
  filter(is.na(PercCovered) & !is.na(SpeciesAcres))
# none

# duplicates?
inv_fwc %>%
  group_by(TaxonName, AreaOfInterestID, GSYear) %>%
  mutate(surveys = n_distinct(SurveyDate)) %>%
  ungroup() %>%
  filter(surveys > 1)
# none

get_dupes(inv_fwc, TaxonName, AreaOfInterestID, GSYear)

# remove missing years
inv_fwc2 <- inv_fwc %>%
  filter(!is.na(PercCovered))

# waterbody-level rep
inv_fwc2 %>%
  group_by(AreaOfInterestID, TaxonName) %>%
  summarize(surveys = n_distinct(SurveyDate)) %>%
  ungroup() %>%
  ggplot(aes(x = surveys)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ TaxonName)

# save
write_csv(inv_fwc2, "intermediate-data/FWC_only_invasive_plant_formatted.csv")


#### edit data for water quality ####