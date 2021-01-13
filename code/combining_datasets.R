#### info ####

# goal: combine datasets
# naming structure:
# camel case for column names
# suffix _FWC or _LW to indicate column source
# units after underscores


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)

# import data
ctrl_old <- read_csv("original-data/PrePMARS_IPMData.csv")
ctrl <- read_csv("original-data/FWC_Herbicide_Treatments_mid2010-Oct2020.csv")
qual <- read_csv("original-data/Lakewatch_1986-2019_TP_TN_Chl_Secchi_Color_Cond.csv")
qual_base <- read_csv("original-data/Lakewatch_Base_File_for_Amy_2020.csv")
lw_plant <- read_csv("original-data/Lakewatch_Plant_Surveys.csv")
fwc_plant <- read_csv("original-data/FWC Plant Surveys.csv")
fwc_plant_new <- read_csv("original-data/FWCSurvey_2018_2019.csv")
gnis <- read_csv("original-data/LW_matching_Herbicide_lakes_with_GNIS.csv")
depth <- read_csv("original-data/Mean_Depth_Data.csv")
vol <- read_csv("original-data/Volume_Calculation.csv")
hydro <- read_csv("original-data/Florida_Lakes_National_Hydrography_Dataset.csv")


#### old control data ####

# see old_herbicide_initial_visualizations.R

# remove entries missing AreaOfInterest (don't include county either)
# capitalize lake and county
# rename columns
ctrl_old2 <- ctrl_old %>%
  filter(!is.na(AreaOfInterest)) %>%
  mutate(County_FWC = toupper(County) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")) %>%
           str_replace("SAINTJOHNS", "SAINT JOHNS"),
         AreaOfInterest = ifelse(AreaOfInterest == "Watermellon Pond", "Watermelon Pond", AreaOfInterest),
         Lake_FWC = str_replace_all(AreaOfInterest, c(", Lake" = "", "Lake " = "", " Lake" = "","'" = "")) %>%
           toupper(),
         Lake_FWC = case_when(AreaOfInterest == "Ella, Lake" & County_FWC == "LAKE" ~ "ELLA 2",
                              Lake_FWC == "HALF MOON" & County_FWC == "MARION" ~ "HALFMOON",
                              Lake_FWC == "LITTLE RED WATER" & County_FWC == "HIGHLANDS" ~ "LITTLE REDWATER",
                              TRUE ~ Lake_FWC)) %>%
  select(-County) %>%
  rename("Longitude_FWC" = "longitude",
         "Latitude_FWC" = "latitude",
         "SpeciesOrig" = "Species_orig",
         "TotalContFWC" = "Total_cont_fwc") %>%
  unique()

# duplicate name/county combinations
ctrl_old2 %>%
  mutate(Lake_County = paste(Lake_FWC, County_FWC, sep = "_")) %>%
  select(AreaOfInterest, AreaOfInterestID, Lake_County) %>%
  unique() %>%
  mutate(dups = duplicated(Lake_County)) %>%
  filter(dups == T)
# no


#### control data ####

# see herbicide_initial_visualizations.R

# capitalize lake and county
# format date
# rename columns
ctrl2 <- ctrl %>%
  mutate(County_FWC = toupper(County) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")),
         Lake_FWC = str_replace_all(AreaOfInterest, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper(),
         Lake_FWC = case_when(AreaOfInterest == "Ella, Lake" & County_FWC == "LAKE" ~ "ELLA 2",
                              Lake_FWC == "HALF MOON" & County_FWC == "MARION" ~ "HALFMOON",
                              Lake_FWC == "LITTLE RED WATER" & County_FWC == "HIGHLANDS" ~ "LITTLE REDWATER",
                              TRUE ~ Lake_FWC),
         BeginDate = as.Date(BeginDate, "%m/%d/%y")) %>%
  select(-County) %>%
  rename("Year" = "year",
         "Month" = "month",
         "Species" = "species",
         "Herbicide" = "herbicide",
         "TotalHerbicideUsed" = "totalherbicideused",
         "Longitude_FWC" = "longitude",
         "Latitude_FWC" = "latitude") %>%
  unique()

# duplicate name/county combinations
ctrl2 %>%
  mutate(Lake_County = paste(Lake_FWC, County_FWC, sep = "_")) %>%
  select(AreaOfInterest, AreaOfInterestID, Lake_County) %>%
  unique() %>%
  mutate(dups = duplicated(Lake_County)) %>%
  filter(dups == T)


#### water quality data ####

# see quality_initial_visualizations.R

# capitalize lake and county
# correct lake spellingS
# format date
# make secchi disk data numeric
qual2 <- qual %>%
  mutate(County_LW = toupper(County) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")),
         County_LW = ifelse(County_LW == "SAINTLUCIE", "SAINT LUCIE", County_LW),
         Lake_LW = str_replace_all(Lake, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper(),
         Lake_LW = case_when(Lake_LW == "LITTLE RED FISH" & County_LW == "WALTON" ~ "LITTLE REDFISH",
                             Lake_LW == "REDWATER" & County_LW == "HIGHLANDS" ~ "RED WATER",
                             Lake_LW == "BIVANS ARM" & County_LW == "ALACHUA" ~ "BIVENS ARM",
                             Lake_LW == "ALLIGATOR SOUTH" & County_LW == "COLUMBIA" ~ "ALLIGATOR",
                             TRUE ~ Lake_LW),
         Date = as.Date(Date, "%m/%d/%y"),
         SecchiCombined = ifelse(is.na(SECCHI_ft), SECCHI_2, SECCHI_ft) %>%
           tolower(),
         SecchiBottom = ifelse(str_detect(SecchiCombined, "bottom") == T, 1, 0),
         SecchiWeeds = ifelse(str_detect(SecchiCombined, "weeds") == T, 1, 0),
         Secchi = case_when(SecchiCombined == "." ~ NA_character_, 
                            SecchiCombined == "weeds" ~ NA_character_,
                            SecchiCombined == "weeds (surface)" ~ NA_character_,
                            SecchiCombined == "bottom" ~ NA_character_,
                            TRUE ~ SecchiCombined) %>%
           parse_number()) %>%
  select(-c(County, Lake)) %>%
  rename("Secchi1_ft" = "SECCHI_ft",
         "Secchi2_ft" = "SECCHI_2") %>%
  unique()


#### lakewatch plant survey data ####

# see lakewatch_plant_initial_visualizations.R

# capitalize lake and county
# format date
# combine biomass values
# use the lower value of Frequency_percent and N_rows/Stations for species occurrance
# rename columns
# use max value for species with multiple rows per sample (n = 1)
lw_plant2 <- lw_plant %>%
  mutate(County_LW = toupper(County) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")),
         County_LW = ifelse(County_LW == "SAINTLUCIE", "SAINT LUCIE", County_LW),
         Lake_LW = str_replace_all(Lake, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper(),
         Lake_LW = case_when(Lake_LW == "LITTLE RED FISH" & County_LW == "WALTON" ~ "LITTLE REDFISH",
                             Lake_LW == "REDWATER" & County_LW == "HIGHLANDS" ~ "RED WATER",
                             Lake_LW == "MILLDAM" & County_LW == "MARION" ~ "MILL DAM",
                             TRUE ~ Lake_LW),
         Date = as.Date(paste(Month, Day, Year, sep = "-"), "%m-%d-%Y"),
         TotalBiomass_kg_m2 = rowSums(.[c("Em_biomass_kg_m2", "Fl_biomass_kg_m2", "Sub_biomass_kg_m2")]),
         RowsStations = N_rows / Stations * 100,
         RowsStations = case_when(RowsStations - round(RowsStations) == 0.5 ~ RowsStations + 0.5,
                                   TRUE ~ round(RowsStations)),
         SpeciesFrequency = ifelse(RowsStations < round(Frequency_percent),
                                    RowsStations/100,
                                    Frequency_percent/100)) %>%
  select(-c(County, Lake)) %>%
  rename("EmFlZoneWidth_ft" = "Em_fl_zone_width_ft",
         "EmBiomass_kg_m2" = "Em_biomass_kg_m2",
         "FlBiomass_kg_m2" = "Fl_biomass_kg_m2",
         "SubBiomass_kg_m2" = "Sub_biomass_kg_m2",
         "LakeDepth_m" = "Lake_depth_m",
         "PercentAreaCovered" = "Percent_area_covered",
         "PercentVolumeInhabited" = "Percent_volume_inhabited",
         "StationsPresent" = "N_rows",
         "FrequencyOrig" = "Frequency_percent",
         "CommonName" = "Common_name",
         "GenusSpecies" = "Genus_species") %>%
  group_by(Year, Month, Day, Date, Stations, EmFlZoneWidth_ft, EmBiomass_kg_m2, FlBiomass_kg_m2, SubBiomass_kg_m2, LakeDepth_m, PercentAreaCovered, PercentVolumeInhabited, CommonName, Genus, Species, GenusSpecies, County_LW, Lake_LW, TotalBiomass_kg_m2) %>%
  summarise(RowsStations = max(RowsStations),
            SpeciesFrequency = max(SpeciesFrequency),
            StationsPresent = max(StationsPresent),
            FrequencyOrig = max(FrequencyOrig)) %>%
  ungroup() %>%
  unique()


#### FWC plant survey data ####

# see FWC_plant_initial_visualizations.R

# format dates
# correct a mispelling
# combine dataframes (no overlapping rows
# remove rows with missing species names
# sum cover for duplicate non-species
# use max cover for duplicate species
# add columns
fwc_plant2 <- fwc_plant %>%
  mutate(SurveyDate = as.Date(SurveyDate, "%m/%d/%Y"),
         WaterbodyName = ifelse(WaterbodyName == "Watermellon Pond", "Watermelon Pond", WaterbodyName)) %>%
  full_join(fwc_plant_new %>%
              mutate(SurveyDate = as.Date(SurveyDate, "%m/%d/%y"))) %>%
  filter(!is.na(SpeciesName)) %>%
  group_by(WaterbodyName, WaterbodyAcres, WaterbodyType, County, WMD, SurveyYear, SurveyDate, IsUnableToSurvey, Surveyor, IsDetected, IsAcreageRequire, SpeciesName, Origin, Eppc, Habitat, HabitatShortName, AreaOfInterestID) %>%
  mutate(County_FWC = toupper(County) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")),
         Lake_FWC = str_replace_all(WaterbodyName, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper(),
         Lake_FWC = case_when(Lake_FWC == "HALF MOON" & County_FWC == "MARION" ~ "HALFMOON",
                              Lake_FWC == "LITTLE RED WATER" & County_FWC == "HIGHLANDS" ~ "LITTLE REDWATER",
                              TRUE ~ Lake_FWC),
         SpeciesAcres = case_when(str_detect(SpeciesName, "Filamentous algae|other|spp.|/") == T ~ sum(SpeciesAcres, na.rm = T),
                                     TRUE ~ max(SpeciesAcres, na.rm = T))) %>%
  ungroup() %>%
  unique() %>%
  select(-County) %>%
  mutate(SpeciesAcres = ifelse(SpeciesAcres == -Inf, NA_real_, SpeciesAcres),
         SpeciesFrequency = SpeciesAcres / WaterbodyAcres)
# will give warnings, but these are addressed by making -Inf into NAs

# duplicate name/county combinations
# AOI ID only in new dataset
fwc_plant2 %>%
  mutate(Lake_County = paste(Lake_FWC, County_FWC, sep = "_")) %>%
  select(WaterbodyName, Lake_County) %>%
  unique() %>%
  mutate(dups = duplicated(Lake_County)) %>%
  filter(dups == T)


#### ID data ####

# make all county information uppercase
# remove "lake" from names
gnis2 <- gnis %>%
  mutate(County_FWC = toupper(County_FWC) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")),
         County_LW = toupper(County_LW) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")),
         Lake_LW = str_replace_all(Lake_LW, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper(),
         Lake_FWC = str_replace_all(AreaOfInterest, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper(),
         Lake_FWC = case_when(Lake_FWC == "LITTLE RED WATER" & County_FWC == "HIGHLANDS" ~ "LITTLE REDWATER",
                              TRUE ~ Lake_FWC)) %>%
  unique() %>%
  select(-AreaOfInterest)

# county overlap
gnis2 %>%
  filter(County_FWC != County_LW)
# Apopka: in both counties
# Isabell: in Polk, but on the edge of Highlands
# Monroe: in both counties

# lake overlap
gnis2 %>%
  filter(Lake_LW != Lake_FWC) %>%
  arrange(County_LW, Lake_LW) %>%
  data.frame()
# multiple GNIS values assigned to a single FWC lake
  # Little Red Water in Highlands
  # Blue in Polk County
  # Rodman Reservoir in Putnam County
# others are different spellings

# check for missing data
sum(is.na(gnis2$Lake_LW))
sum(is.na(gnis2$Lake_FWC))
sum(is.na(gnis2$County_LW))
sum(is.na(gnis2$County_FWC))
# none


#### combine datasets 1 ####

# check for missing data data
sum(is.na(ctrl_old2$Year))
sum(is.na(ctrl2$BeginDate))
sum(is.na(qual2$Date))
sum(is.na(lw_plant2$Date))
sum(is.na(fwc_plant2$SurveyDate))

# check for missing AOI data
sum(is.na(ctrl_old2$AreaOfInterestID))
sum(is.na(ctrl_old2$AreaOfInterest))
sum(is.na(ctrl2$AreaOfInterestID))
sum(is.na(ctrl2$AreaOfInterest))
sum(is.na(fwc_plant2$AreaOfInterestID)) # a lot
sum(is.na(fwc_plant2$AreaOfInterest)) # no column

# number of times each lake was sampled for each variable
# start and end dates
# leave out AreaOfInterestID from plant data (missing for a lot)
# leave out AreaOfInterest (some mismatches)
# note that ADMINISTRATION refers to "General lakes"
lakes <- gnis2 %>%
  full_join(ctrl_old2  %>%
              mutate(CtrlOldStart = min(Year),
                     CtrlOldEnd = max(Year)) %>%
              group_by(County_FWC, Lake_FWC, CtrlOldStart, CtrlOldEnd) %>%
              summarise(CtrlOld = length(unique(Year))) %>%
              ungroup()) %>%
  full_join(ctrl2 %>%
              mutate(CtrlStart = min(BeginDate),
                     CtrlEnd = max(BeginDate)) %>%
              group_by(County_FWC, Lake_FWC, CtrlStart, CtrlEnd) %>%
              summarise(Ctrl = length(unique(BeginDate))) %>%
              ungroup()) %>%
  full_join(qual2 %>%
              group_by(County_LW, Lake_LW) %>%
              summarise(Qual = length(unique(Date)),
                        QualStart = min(Date),
                        QualEnd = max(Date)) %>%
              ungroup()) %>%
  full_join(lw_plant2 %>%
              group_by(County_LW, Lake_LW) %>%
              summarise(LWPlant = length(unique(Date)),
                        LWPlantStart = min(Date),
                        LWPlantEnd = max(Date)) %>%
              ungroup()) %>%
  full_join(fwc_plant2 %>%
              group_by(County_FWC, Lake_FWC) %>%
              summarise(FWCPlant = length(unique(SurveyDate)),
                        FWCPlantStart = min(SurveyDate),
                        FWCPlantEnd = max(SurveyDate)) %>%
              ungroup())


#### combine water bodies ####

# data table of names
lake_names <- lakes %>%
  select(Lake_LW, Lake_FWC, County_LW, County_FWC, GNIS_ID) %>%
  pivot_longer(-GNIS_ID,
               names_to = c(".value", "source"),
               names_pattern = "(.+)_(.+)") %>%
  mutate(Lake_County = paste(Lake, County, sep = "_")) %>%
  unique()

# export data to use text clustering algorithms in OpenRefine
write_csv(lake_names, "output/lake_names_for_clustering.csv")
# used metaphone3 algorithm
# checked lakes on Google Maps
# made intermediate data file waterbody_split_combine_121620 to list

# different spellings
lakes %>%
  filter(County_FWC == "SAINT JOHNS") %>%
  data.frame()

lakes %>%
  filter(Lake_LW == "LITTLE RED FISH") %>%
  data.frame()

lakes %>%
  filter(Lake_LW == "REDWATER" | Lake_FWC == "RED WATER") %>%
  data.frame()

lakes %>%
  filter(Lake_LW == "ELLA 2" | Lake_LW == "ELLA" | Lake_FWC == "ELLA") %>%
  data.frame()
# two lake Ella's in one county 
# Ella 2 in LW is probably "Ella, Lake" in FWC

# import list of lakes to split/combine
lake_combos <- read_csv("intermediate-data/waterbody_split_combine_121620.csv")

# edit data
lake_combos2 <- lake_combos %>%
  mutate(common_lake_name = toupper(common_lake_name),
         county = toupper(county)) %>%
  filter(waterbody_guess == "combine") %>%
  select(common_lake_name, county)

# do in two steps so that I can check that the matched names are correct

# gnis
gnis3 <- gnis2 %>%
  left_join(lake_combos2 %>%
              rename("County_LW" = "county")) %>%
  left_join(lake_combos2 %>%
              rename("County_FWC" = "county")) %>%
  mutate(Lake_FWC2 = case_when(str_detect(Lake_FWC, common_lake_name) == T ~ common_lake_name,
                               TRUE ~ NA_character_),
         Lake_LW2 = case_when(str_detect(Lake_LW, common_lake_name) == T ~ common_lake_name, 
                              TRUE ~ NA_character_)) %>%
  select(-common_lake_name) %>%
  filter(!is.na(Lake_FWC2) | !(is.na(Lake_LW2))) %>%
  full_join(gnis2) %>%
  mutate(Lake_LW = ifelse(!is.na(Lake_LW2), Lake_LW2, Lake_LW),
         Lake_FWC = ifelse(!is.na(Lake_FWC2), Lake_FWC2, Lake_FWC)) %>%
  select(-c(GNIS_ID, Latitude, Longitude, Lake_LW2, Lake_FWC2)) %>%
  unique()

# ctrl_old
ctrl_old3 <- ctrl_old2 %>%
  left_join(lake_combos2 %>%
              rename("County_FWC" = "county")) %>%
  mutate(Lake_FWC2 = case_when(str_detect(Lake_FWC, common_lake_name) == T ~ common_lake_name,
                               TRUE ~ NA_character_),
         Lake_FWC2 = ifelse(Lake_FWC == "ROCK SPRING RUN", NA_character_, Lake_FWC2)) %>%
  select(-common_lake_name) %>%
  filter(!is.na(Lake_FWC2) & Lake_FWC2 != Lake_FWC) %>%
  full_join(ctrl_old2) %>%
  mutate(Lake_FWC = ifelse(!is.na(Lake_FWC2), Lake_FWC2, Lake_FWC)) %>%
  select(-Lake_FWC2) %>%
  unique()

# ctrl
ctrl3 <- ctrl2 %>%
  left_join(lake_combos2 %>%
              rename("County_FWC" = "county")) %>%
  mutate(Lake_FWC2 = case_when(str_detect(Lake_FWC, common_lake_name) == T ~ common_lake_name,
                               TRUE ~ NA_character_)) %>%
  select(-common_lake_name) %>%
  filter(!is.na(Lake_FWC2) & Lake_FWC2 != Lake_FWC) %>%
  full_join(ctrl2) %>%
  mutate(Lake_FWC = ifelse(!is.na(Lake_FWC2), Lake_FWC2, Lake_FWC)) %>%
  select(-Lake_FWC2) %>%
  unique()

# fwc_plant
fwc_plant3 <- fwc_plant2 %>%
  left_join(lake_combos2 %>%
              rename("County_FWC" = "county")) %>%
  mutate(Lake_FWC2 = case_when(str_detect(Lake_FWC, common_lake_name) == T ~ common_lake_name,
                               TRUE ~ NA_character_),
         Lake_FWC2 = ifelse(Lake_FWC %in% c("ROCK SPRING RUN", "LITTLE HARRIS"), 
                            NA_character_, 
                            Lake_FWC2)) %>%
  select(-common_lake_name) %>%
  filter(!is.na(Lake_FWC2) & Lake_FWC2 != Lake_FWC) %>%
  full_join(fwc_plant2) %>%
  mutate(Lake_FWC = ifelse(!is.na(Lake_FWC2), Lake_FWC2, Lake_FWC)) %>%
  select(-Lake_FWC2) %>%
  unique()

# water quality
qual3 <- qual2 %>%
  left_join(lake_combos2 %>%
              rename("County_LW" = "county")) %>%
  mutate(Lake_LW2 = case_when(str_detect(Lake_LW, common_lake_name) == T ~ common_lake_name,
                               TRUE ~ NA_character_),
         Lake_LW2 = ifelse(Lake_LW %in% c("ELDORADO", "LITTLE HARRIS", "MOUNT DORA", "MOORE POND", "CRYSTAL BOWL", "MIDDLE BEAR"),
                           NA_character_,
                           Lake_LW2)) %>%
  select(-common_lake_name) %>%
  filter(!is.na(Lake_LW2) & Lake_LW2 != Lake_LW) %>%
  full_join(qual2) %>%
  mutate(Lake_LW = ifelse(!is.na(Lake_LW2), Lake_LW2, Lake_LW)) %>%
  select(-Lake_LW2) %>%
  unique()
# removed some duplicates

# lw_plant
lw_plant3 <- lw_plant2 %>%
  left_join(lake_combos2 %>%
              rename("County_LW" = "county")) %>%
  mutate(Lake_LW2 = case_when(str_detect(Lake_LW, common_lake_name) == T ~ common_lake_name,
                              TRUE ~ NA_character_),
         Lake_LW2 = ifelse(Lake_LW %in% c("ELDORADO", "LITTLE HARRIS", "MOUNT DORA", "MOORE POND", "CRYSTAL BOWL", "MIDDLE BEAR"),
                           NA_character_,
                           Lake_LW2)) %>%
  select(-common_lake_name) %>%
  filter(!is.na(Lake_LW2) & Lake_LW2 != Lake_LW) %>%
  full_join(lw_plant2) %>%
  mutate(Lake_LW = ifelse(!is.na(Lake_LW2), Lake_LW2, Lake_LW)) %>%
  select(-Lake_LW2) %>%
  unique()


#### combine datasets 2 ####

# FWC data
lakes_FWC <- ctrl_old3  %>%
  mutate(CtrlOldStart = min(Year),
         CtrlOldEnd = max(Year)) %>%
  group_by(County_FWC, Lake_FWC, CtrlOldStart, CtrlOldEnd) %>%
  summarise(CtrlOld = length(unique(Year))) %>%
  ungroup() %>%
  full_join(ctrl3 %>%
              mutate(CtrlStart = min(BeginDate),
                     CtrlEnd = max(BeginDate)) %>%
              group_by(County_FWC, Lake_FWC, CtrlStart, CtrlEnd) %>%
              summarise(Ctrl = length(unique(BeginDate))) %>%
              ungroup()) %>%
    full_join(fwc_plant3 %>%
                group_by(County_FWC, Lake_FWC) %>%
                summarise(FWCPlant = length(unique(SurveyDate)),
                          FWCPlantStart = min(SurveyDate),
                          FWCPlantEnd = max(SurveyDate)) %>%
                ungroup()) %>%
  left_join(gnis3) %>%
  mutate(Lake_LW = ifelse(is.na(Lake_LW), Lake_FWC, Lake_LW),
         County_LW = ifelse(is.na(County_LW), County_FWC, County_LW))

# LW data
lakes_LW <- qual3 %>%
  group_by(County_LW, Lake_LW) %>%
  summarise(Qual = length(unique(Date)),
            QualStart = min(Date),
            QualEnd = max(Date)) %>%
  ungroup() %>%
  full_join(lw_plant3 %>%
              group_by(County_LW, Lake_LW) %>%
              summarise(LWPlant = length(unique(Date)),
                        LWPlantStart = min(Date),
                        LWPlantEnd = max(Date)) %>%
              ungroup()) %>%
  left_join(gnis3) %>%
  mutate(Lake_FWC = ifelse(is.na(Lake_FWC), Lake_LW, Lake_FWC),
         County_FWC = ifelse(is.na(County_FWC), County_LW, County_FWC))

# combine data
lakes2 <- lakes_FWC %>%
  full_join(lakes_LW) %>%
  mutate(CtrlOld = replace_na(CtrlOld, 0),
         Ctrl = replace_na(Ctrl, 0),
         Qual = replace_na(Qual, 0),
         LWPlant = replace_na(LWPlant, 0),
         FWCPlant = replace_na(FWCPlant, 0),
         AnyCtrl = CtrlOld + Ctrl,
         AnyPlant = LWPlant + FWCPlant,
         DataType = case_when(AnyCtrl > 0 & Qual > 0 & AnyPlant > 0 ~ "all",
                              AnyCtrl > 0 & Qual > 0 & AnyPlant == 0 ~ "ctrl + quality",
                              AnyCtrl > 0 & Qual == 0 & AnyPlant > 0 ~ "ctrl + plant",
                              AnyCtrl == 0 & Qual > 0 & AnyPlant > 0 ~ "quality + plant",
                              AnyCtrl > 0 & Qual == 0 & AnyPlant == 0 ~ "ctrl only",
                              AnyCtrl == 0 & Qual > 0 & AnyPlant == 0 ~ "quality only",
                              AnyCtrl == 0 & Qual == 0 & AnyPlant > 0 ~ "plant only",
                              AnyCtrl == 0 & Qual == 0 & AnyPlant == 0 ~ "none") %>%
           fct_relevel("all", "ctrl + quality", "ctrl + plant", "quality + plant", "ctrl only", "quality only", "plant only"))

# counts per data type
pdf("output/combined_data_types.pdf")
lakes2 %>%
  group_by(DataType) %>%
  count() %>%
  ggplot(aes(x = DataType, y = n)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust = -0.25, size = 3) +
    xlab("Data available") +
    ylab("Waterbodies") +
    theme_bw()
dev.off()

# time lines
unique(lakes2$CtrlOldStart)

lake_times <- lakes2 %>%
  rowwise() %>%
  mutate(AnyCtrlStart = ifelse(AnyCtrl > 0,"1/1/1998", NA_character_) %>%
           as.Date("%m/%d/%Y"),
         AnyCtrlEnd = CtrlEnd,
         AnyPlantStart = ifelse(AnyPlant > 0, as.character(min(c(FWCPlantStart, LWPlantStart), na.rm = T)), NA_character_) %>%
           as.Date("%Y-%m-%d"),
         AnyPlantEnd = ifelse(AnyPlant > 0, as.character(max(c(FWCPlantEnd, LWPlantEnd), na.rm = T)), NA_character_) %>%
           as.Date("%Y-%m-%d")) %>%
  ungroup() %>%
  select(Lake_LW, County_LW, AnyCtrlStart, AnyCtrlEnd, AnyPlantStart, AnyPlantEnd, QualStart, QualEnd) %>%
  rename(start_ctrl = AnyCtrlStart,
         end_ctrl = AnyCtrlEnd,
         start_plant = AnyPlantStart,
         end_plant = AnyPlantEnd,
         start_quality = QualStart,
         end_quality = QualEnd) %>%
  pivot_longer(cols = -c(Lake_LW, County_LW),
               names_to = c(".value", "DataType"),
               names_pattern = "(.+)_(.+)") %>%
  pivot_longer(cols = c(start, end),
               names_to = "TimePoints",
               values_to = "Date") %>%
  mutate(LakeCounty = paste(Lake_LW, County_LW, sep = "_") %>%
           as.factor(),
         LakeGroup = cut(as.numeric(LakeCounty), breaks = 5),
         LakeCounty = fct_rev(LakeCounty))

lake_times_grp = sort(unique(lake_times$LakeGroup))

pdf("output/combined_data_time_series.pdf")
for(i in 1:length(lake_times_grp)){
  
  lake_times_sub <- lake_times %>%
    filter(LakeGroup == lake_times_grp[i])
  
  print(ggplot(data = lake_times_sub,
               aes(x = Date,
                   y = LakeCounty,
                   color = DataType,
                   linetype = DataType,
                   shape = DataType)) +
          geom_line(size = 0.3, alpha = 0.5) +
          geom_point(size = 0.3) +
          xlab("Date") +
          ylab("Waterbody") +
          theme_bw() +
          theme(axis.text.y = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()))
  
}
dev.off()



#### check combined data ####

# missing plant data (n = 6)
lakes2 %>%
  filter(DataType == "ctrl + quality") %>%
  data.frame()

lakes2 %>%
  filter(DataType == "plant only" &
           (str_detect(Lake_LW, "EM|MARSH|CONS|GRADY|WOOD|TUCKER|BONNET|BLUE") == T |
              str_detect(Lake_FWC, "EM|MARSH|CONS|GRADY|WOOD|TUCKER|BONNET|BLUE") == T)) %>%
  data.frame()
# probably don't have plant data 

# missing quality data (n = 154)
lakes2 %>%
  filter(DataType == "ctrl + plant") %>%
  filter(!is.na(LWPlantStart)) %>%
  select(Lake_LW, County_LW)
# 3 lakes have LW plant data, should have quality data

lakes2 %>%
  filter(!is.na(QualStart) & str_detect(Lake_LW, "ALLIGATOR|TOWNSEND|THOMAS")) %>%
  select(Lake_LW, County_LW, DataType)
# called "Alligator South" in quality dataset (corrected above)

# control data missing FWC plant data
fwc_lake_match <- lakes2 %>%
  filter((AnyCtrl > 0 & is.na(FWCPlantStart) | (AnyCtrl == 0 & !is.na(FWCPlantStart)))) %>%
  select(Lake_FWC, County_FWC, AnyCtrl, FWCPlantStart)

# export to use clustering in OpenRefine
write_csv(fwc_lake_match, "output/fwc_lake_names_for_clustering.csv")
# checked all clustering algorithms and didn't find any to group


#### combine datasets hydrilla ####

# volume data
acres <- vol %>%
  mutate(County_LW = toupper(County) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")),
         County_LW = ifelse(County_LW == "SAINTLUCIE", "SAINT LUCIE", County_LW),
         Lake_LW = str_replace_all(Lake, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper(),
         Lake_LW = case_when(Lake_LW == "LITTLE RED FISH" & County_LW == "WALTON" ~ "LITTLE REDFISH",
                             Lake_LW == "REDWATER" & County_LW == "HIGHLANDS" ~ "RED WATER",
                             Lake_LW == "BIVANS ARM" & County_LW == "ALACHUA" ~ "BIVENS ARM",
                             Lake_LW == "ALLIGATOR SOUTH" & County_LW == "COLUMBIA" ~ "ALLIGATOR",
                             Lake_LW == "MILLDAM" & County_LW == "MARION" ~ "MILL DAM",
                             TRUE ~ Lake_LW)) %>%
  group_by(County_LW, Lake_LW) %>%
  summarise(SurfaceAreaAcres = mean(Surface_area_acres)) %>%
  ungroup()

# FWC data
hydrilla_FWC <- ctrl_old3  %>%
  filter(Species == "Hydrilla verticillata") %>%
  mutate(CtrlOldStart = min(Year),
         CtrlOldEnd = max(Year),
         CtrlOldYears = CtrlOldEnd - CtrlOldStart) %>%
  group_by(County_FWC, Lake_FWC, CtrlOldStart, CtrlOldEnd, CtrlOldYears) %>%
  summarise(CtrlOld = length(unique(Year)),
            CtrlOldAcres = mean(TotalAcres, na.rm = T),
            CtrlOldFirst = paste("01/01/", min(Year), sep = "") %>%
              as.Date("%m/%d/%Y"),
            CtrlOldLast = paste("01/01/", max(Year), sep = "") %>%
              as.Date("%m/%d/%Y")) %>%
  ungroup() %>%
  mutate(CtrlOldFreq = CtrlOld/CtrlOldYears) %>%
  full_join(ctrl3 %>%
              filter(Species == "Hydrilla verticillata") %>%
              mutate(CtrlStart = min(BeginDate),
                     CtrlEnd = max(BeginDate),
                     CtrlYears = as.numeric((CtrlEnd - CtrlStart)/365)) %>%
              group_by(County_FWC, Lake_FWC, CtrlStart, CtrlEnd, CtrlYears) %>%
              summarise(Ctrl = length(unique(BeginDate)),
                        CtrlAcres = mean(TotalAcres, na.rm = T),
                        CtrlFirst = min(BeginDate),
                        CtrlLast = max(BeginDate)) %>%
              ungroup() %>%
              mutate(CtrlFreq = Ctrl/CtrlYears)) %>%
  left_join(fwc_plant3 %>%
              filter(SpeciesFrequency <= 1) %>%
              group_by(County_FWC, Lake_FWC) %>%
              summarise(FWCPlant = length(unique(SurveyDate)),
                        FWCPlantStart = min(SurveyDate),
                        FWCPlantEnd = max(SurveyDate)) %>%
              ungroup()) %>%
  left_join(gnis3) %>%
  mutate(Lake_LW = ifelse(is.na(Lake_LW), Lake_FWC, Lake_LW),
         County_LW = ifelse(is.na(County_LW), County_FWC, County_LW))

# combine data
hydrilla_lakes <- hydrilla_FWC %>%
  left_join(lakes_LW) %>%
  left_join(acres) %>%
  mutate(CtrlOld = replace_na(CtrlOld, 0),
         Ctrl = replace_na(Ctrl, 0),
         LWPlant = replace_na(LWPlant, 0),
         FWCPlant = replace_na(FWCPlant, 0),
         AnyCtrl = CtrlOld + Ctrl,
         AnyPlant = LWPlant + FWCPlant,
         DataType = case_when(AnyCtrl > 0 & AnyPlant > 0 ~ "ctrl + plant",
                              AnyCtrl > 0 & AnyPlant == 0 ~ "ctrl only",
                              AnyPlant > 0 ~ "plant only",
                              AnyCtrl == 0 & AnyPlant == 0 ~ "none") %>%
           fct_relevel("ctrl + plant", "ctrl only", "plant only"),
         AnyCtrlFreq = case_when(CtrlOld > 0 & Ctrl > 0 ~ (Ctrl + CtrlOld)/(CtrlYears + CtrlOldYears),
                                 CtrlOld > 0 & Ctrl == 0 ~ CtrlOld/(unique(hydrilla_FWC$CtrlYears)[1] + CtrlOldYears),
                                 CtrlOld == 0 & Ctrl > 0 ~ Ctrl/(CtrlYears + unique(hydrilla_FWC$CtrlOldYears)[1]),
                                 TRUE ~ 0)) %>%
  rowwise() %>%
  mutate(AnyCtrlAcres = ifelse(AnyCtrl > 0, mean(c(CtrlAcres, CtrlOldAcres), na.rm = T)/SurfaceAreaAcres, 0),
         AnyCtrlFirst = min(c(CtrlOldFirst, CtrlFirst), na.rm = T),
         AnyCtrlLast = max(c(CtrlOldLast, CtrlLast), na.rm = T),
         AnyPlantFirst = min(c(FWCPlantStart, LWPlantStart), na.rm = T),
         AnyPlantLast = max(c(FWCPlantEnd, LWPlantEnd), na.rm = T)) %>%
  ungroup() %>%
  mutate(AnyCtrlAcres = ifelse(AnyCtrlAcres > 1, 1, AnyCtrlAcres))

# save figures
pdf("output/hydrilla_data_types.pdf")

# data types
hydrilla_lakes %>%
  group_by(DataType) %>%
  count() %>%
  ggplot(aes(x = DataType, y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.25, size = 3) +
  xlab("Hydrilla data available") +
  ylab("Waterbodies") +
  theme_bw()

# frequency histogram
hydrilla_lakes %>%
  mutate(FreqGroup = ifelse(AnyCtrlFreq <= 1, "<= 1", "> 1")) %>%
  ggplot(aes(x = AnyCtrlFreq)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ FreqGroup, scales = "free_x") +
  xlab(expression(paste("Hydrilla treatment frequency (", year^-1, ")", sep = ""))) +
  ylab("Waterbodies") +
  theme_bw() +
  theme(strip.background = element_blank())

# intensity histogram
hydrilla_lakes %>%
  filter(!is.na(AnyCtrlAcres)) %>%
  ggplot(aes(x = AnyCtrlAcres)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Proportion of lake treated for Hydrilla") +
  ylab("Waterbodies") +
  theme_bw() +
  theme(strip.background = element_blank())

# frequency by intensity
hydrilla_lakes %>%
  filter(!is.na(AnyCtrlAcres)) %>%
  ggplot(aes(x = AnyCtrlAcres, y = AnyCtrlFreq)) +
  geom_point() +
  xlab("Proportion of lake treated for Hydrilla") +
  ylab(expression(paste("Hydrilla treatment frequency (", year^-1, ")", sep = ""))) +
  theme_bw()

# close pdf
dev.off()

# timeline
hydrilla_lakes %>%
  filter(AnyPlantFirst < AnyCtrlFirst) # surveys before herbicide: 311

hydrilla_lakes %>%
  filter(AnyPlantLast > AnyCtrlFirst) # surveys after herbicide: 318

hydrilla_lakes %>%
  filter(AnyPlantFirst < AnyCtrlFirst & AnyPlantLast > AnyCtrlFirst) # surveys before and after herbicide: 303

hydrilla_lakes %>%
  filter(QualStart < AnyCtrlFirst) # quality before herbicide: 194

hydrilla_lakes %>%
  filter(QualEnd > AnyCtrlFirst) # quality after herbicide: 200

hydrilla_lakes %>%
  filter(QualStart < AnyCtrlFirst & QualEnd > AnyCtrlFirst) # quality before and after herbicide: 162

hydrilla_lakes %>%
  filter(AnyPlantFirst < AnyCtrlFirst & AnyPlantLast > AnyCtrlFirst &
           QualStart < AnyCtrlFirst & QualEnd > AnyCtrlFirst) # everything before and after: 158


#### combine datasets floating ####

# FWC data
floating_FWC <- ctrl_old3  %>%
  filter(str_detect(Species, "Eichhornia") == T | str_detect(Species, "Pistia") == T) %>%
  mutate(CtrlOldStart = min(Year),
         CtrlOldEnd = max(Year),
         CtrlOldYears = CtrlOldEnd - CtrlOldStart) %>%
  group_by(County_FWC, Lake_FWC, CtrlOldStart, CtrlOldEnd, CtrlOldYears) %>%
  summarise(CtrlOld = length(unique(Year)),
            CtrlOldAcres = mean(TotalAcres, na.rm = T),
            CtrlOldFirst = paste("01/01/", min(Year), sep = "") %>%
              as.Date("%m/%d/%Y"),
            CtrlOldLast = paste("01/01/", max(Year), sep = "") %>%
              as.Date("%m/%d/%Y")) %>%
  ungroup() %>%
  mutate(CtrlOldFreq = CtrlOld/CtrlOldYears) %>%
  full_join(ctrl3 %>%
              filter(str_detect(Species, "Eichhornia") == T | str_detect(Species, "Pistia") == T) %>%
              mutate(CtrlStart = min(BeginDate),
                     CtrlEnd = max(BeginDate),
                     CtrlYears = as.numeric((CtrlEnd - CtrlStart)/365)) %>%
              group_by(County_FWC, Lake_FWC, CtrlStart, CtrlEnd, CtrlYears) %>%
              summarise(Ctrl = length(unique(BeginDate)),
                        CtrlAcres = mean(TotalAcres, na.rm = T),
                        CtrlFirst = min(BeginDate),
                        CtrlLast = max(BeginDate)) %>%
              ungroup() %>%
              mutate(CtrlFreq = Ctrl/CtrlYears)) %>%
  left_join(fwc_plant3 %>%
              filter(SpeciesName %in% c("Eichhornia crassipes", "Pistia stratiotes") & SpeciesFrequency <= 1) %>%
              group_by(County_FWC, Lake_FWC) %>%
              summarise(FWCPlant = length(unique(SurveyDate)),
                        FWCPlantStart = min(SurveyDate),
                        FWCPlantEnd = max(SurveyDate)) %>%
              ungroup()) %>%
  left_join(gnis3) %>%
  mutate(Lake_LW = ifelse(is.na(Lake_LW), Lake_FWC, Lake_LW),
         County_LW = ifelse(is.na(County_LW), County_FWC, County_LW))

# combine data
floating_lakes <- floating_FWC%>%
  left_join(lakes_LW) %>%
  left_join(acres) %>%
  mutate(CtrlOld = replace_na(CtrlOld, 0),
         Ctrl = replace_na(Ctrl, 0),
         LWPlant = replace_na(LWPlant, 0),
         FWCPlant = replace_na(FWCPlant, 0),
         AnyCtrl = CtrlOld + Ctrl,
         AnyPlant = LWPlant + FWCPlant,
         DataType = case_when(AnyCtrl > 0 & AnyPlant > 0 ~ "ctrl + plant",
                              AnyCtrl > 0 & AnyPlant == 0 ~ "ctrl only",
                              AnyPlant > 0 ~ "plant only",
                              AnyCtrl == 0 & AnyPlant == 0 ~ "none") %>%
           fct_relevel("ctrl + plant", "ctrl only", "plant only"),
         AnyCtrlFreq = case_when(CtrlOld > 0 & Ctrl > 0 ~ (Ctrl + CtrlOld)/(CtrlYears + CtrlOldYears),
                                 CtrlOld > 0 & Ctrl == 0 ~ CtrlOld/(unique(hydrilla_FWC$CtrlYears)[1] + CtrlOldYears),
                                 CtrlOld == 0 & Ctrl > 0 ~ Ctrl/(CtrlYears + unique(hydrilla_FWC$CtrlOldYears)[1]),
                                 TRUE ~ 0)) %>%
  rowwise() %>%
  mutate(AnyCtrlAcres = ifelse(AnyCtrl > 0, mean(c(CtrlAcres, CtrlOldAcres), na.rm = T)/SurfaceAreaAcres, 0),
         AnyCtrlFirst = min(c(CtrlOldFirst, CtrlFirst), na.rm = T),
         AnyCtrlLast = max(c(CtrlOldLast, CtrlLast), na.rm = T),
         AnyPlantFirst = min(c(FWCPlantStart, LWPlantStart), na.rm = T),
         AnyPlantLast = max(c(FWCPlantEnd, LWPlantEnd), na.rm = T)) %>%
  ungroup() %>%
  mutate(AnyCtrlAcres = ifelse(AnyCtrlAcres > 1, 1, AnyCtrlAcres))

# save figures
pdf("output/floating_plant_data_types.pdf")

# data types
floating_lakes %>%
  group_by(DataType) %>%
  count() %>%
  ggplot(aes(x = DataType, y = n)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.25, size = 3) +
  xlab("Floating plant data available") +
  ylab("Waterbodies") +
  theme_bw()

# frequency histogram
floating_lakes %>%
  mutate(FreqGroup = ifelse(AnyCtrlFreq <= 1, "<= 1", "> 1")) %>%
  ggplot(aes(x = AnyCtrlFreq)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ FreqGroup, scales = "free_x") +
  xlab(expression(paste("Floating plant treatment frequency (", year^-1, ")", sep = ""))) +
  ylab("Waterbodies") +
  theme_bw() +
  theme(strip.background = element_blank())

# intensity histogram
floating_lakes %>%
  filter(!is.na(AnyCtrlAcres)) %>%
  ggplot(aes(x = AnyCtrlAcres)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Proportion of lake treated for floating plants") +
  ylab("Waterbodies") +
  theme_bw() +
  theme(strip.background = element_blank())

# frequency by intensity
floating_lakes %>%
  filter(!is.na(AnyCtrlAcres)) %>%
  ggplot(aes(x = AnyCtrlAcres, y = AnyCtrlFreq)) +
  geom_point() +
  xlab("Proportion of lake treated for floating plants") +
  ylab(expression(paste("Floating plant treatment frequency (", year^-1, ")", sep = ""))) +
  theme_bw()

# close pdf
dev.off()

# timeline
floating_lakes %>%
  filter(AnyPlantFirst < AnyCtrlFirst) # surveys before herbicide: 324

floating_lakes %>%
  filter(AnyPlantLast > AnyCtrlFirst) # surveys after herbicide: 336

floating_lakes %>%
  filter(AnyPlantFirst < AnyCtrlFirst & AnyPlantLast > AnyCtrlFirst) # surveys before and after herbicide: 305

floating_lakes %>%
  filter(QualStart < AnyCtrlFirst) # quality before herbicide: 182

floating_lakes %>%
  filter(QualEnd > AnyCtrlFirst) # quality after herbicide: 205

floating_lakes %>%
  filter(QualStart < AnyCtrlFirst & QualEnd > AnyCtrlFirst) # quality before and after herbicide: 150

floating_lakes %>%
  filter(AnyPlantFirst < AnyCtrlFirst & AnyPlantLast > AnyCtrlFirst &
           QualStart < AnyCtrlFirst & QualEnd > AnyCtrlFirst) # everything before and after: 135
