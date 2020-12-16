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
lw_plant <- read_csv("original-data/Lakewatch_Plant_Surveys.csv")
fwc_plant <- read_csv("original-data/FWC Plant Surveys.csv")
fwc_plant_new <- read_csv("original-data/FWCSurvey_2018_2019.csv")
gnis <- read_csv("original-data/LW_matching_Herbicide_lakes_with_GNIS.csv")
depth <- read_csv("original-data/Mean_Depth_Data.csv")
vol <- read_csv("original-data/Volume_Calculation.csv")


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
         Lake_FWC = str_replace_all(AreaOfInterest, c(", Lake" = "", "Lake " = "", " Lake" = "","'" = "")) %>%
           toupper(),
         Lake_FWC = case_when(AreaOfInterest == "Ella, Lake" & County_FWC == "LAKE" ~ "ELLA 2",
                              TRUE ~ Lake_FWC)) %>%
  select(-County) %>%
  rename("Longitude_FWC" = "longitude",
         "Latitude_FWC" = "latitude",
         "SpeciesOrig" = "Species_orig",
         "TotalContFWC" = "Total_cont_fwc") %>%
  unique()


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


#### water quality data ####

# see quality_initial_visualizations.R

# capitalize lake and county
# correct lake spellingS
# format date
# make secchi disk data numeric
qual2 <- qual %>%
  mutate(County_LW = toupper(County) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")),
         Lake_LW = str_replace_all(Lake, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper(),
         Lake_LW = case_when(Lake_LW == "LITTLE RED FISH" & County_LW == "WALTON" ~ "LITTLE REDFISH",
                             Lake_LW == "REDWATER" & County_LW == "HIGHLANDS" ~ "RED WATER",
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
         Lake_LW = str_replace_all(Lake, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper(),
         Lake_LW = case_when(Lake_LW == "LITTLE RED FISH" & County_LW == "WALTON" ~ "LITTLE REDFISH",
                             Lake_LW == "REDWATER" & County_LW == "HIGHLANDS" ~ "RED WATER",
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
         Lake_FWC = case_when(AreaOfInterest == "Ella, Lake" & County_FWC == "LAKE" ~ "ELLA 2",
                              TRUE ~ Lake_FWC),
         SpeciesAcres = case_when(str_detect(SpeciesName, "Filamentous algae|other|spp.|/") == T ~ sum(SpeciesAcres, na.rm = T),
                                     TRUE ~ max(SpeciesAcres, na.rm = T))) %>%
  ungroup() %>%
  unique() %>%
  select(-County) %>%
  mutate(SpeciesAcres = ifelse(SpeciesAcres == -Inf, NA_real_, SpeciesAcres),
         SpeciesFrequency = SpeciesAcres / WaterbodyAcres)
# will give warnings, but these are addressed by making -Inf into NAs


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
           toupper()) %>%
  unique()

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


#### combine datasets ####

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
lakes <- gnis2 %>%
  full_join(ctrl_old2  %>%
              mutate(CtrlOldStart = min(Year),
                     CtrlOldEnd = max(Year)) %>%
              group_by(AreaOfInterest, Longitude_FWC, Latitude_FWC, County_FWC, Lake_FWC, CtrlOldStart, CtrlOldEnd) %>%
              summarise(CtrlOld = length(unique(Year))) %>%
              ungroup()) %>%
  full_join(ctrl2 %>%
              mutate(CtrlStart = min(BeginDate),
                     CtrlEnd = max(BeginDate)) %>%
              group_by(AreaOfInterest, Longitude_FWC, Latitude_FWC, County_FWC, Lake_FWC, CtrlStart, CtrlEnd) %>%
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
              ungroup()) %>%
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



#### match data ####

# missing plant data (n = 1)
lakes %>%
  filter(DataType == "ctrl + quality") %>%
  data.frame()

lakes %>%
  filter(DataType == "plant only" &
           (str_detect(Lake_LW, "EM|MARSH|CONS") == T |
           str_detect(Lake_FWC, "EM|MARSH|CONS") == T)) %>%
  data.frame()
# probably don't have plant data

# data table of names
lake_names <- lakes %>%
  filter(is.na(GNIS_ID)) %>%
  select(Lake_LW, Lake_FWC, County_LW, County_FWC) %>%
  pivot_longer(everything(),
               names_to = c(".value", "source"),
               names_pattern = "(.+)_(.+)") %>%
  mutate(Lake_County = paste(Lake, County, sep = "_")) %>%
  unique()

# export data to use text clustering algorithms in OpenRefine
write_csv(lake_names, "output/lake_names_for_clustering.csv")
# used metaphone3 algorithm
# checked lakes on Google Maps

# different spellings
lakes %>%
  filter(County_FWC == "SAINT JOHNS") %>%
  data.frame()

lakes %>%
  filter(Lake_LW %in% c("JANE", "JEAN") & County_LW == "LEON") %>%
  data.frame()
# ask Mark about these

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
# Ella 2 in LW is probably Ella, Lake in FWC

lakes %>%
  filter(Lake_LW == "CLEARWATER 2" | Lake_LW == "CLEARWATER" | Lake_FWC == "CLEARWATER") %>%
  filter(County_LW == "PUTNAM" | County_FWC == "PUTNAM") %>%
  data.frame()

ctrl2 %>%
  filter(AreaOfInterest == "Clearwater Lake" & County_FWC == "PUTNAM" & Lake_FWC == "CLEARWATER")
# didn't have an AOI ID earlier, so it's not in the GNIS dataset
#### START HERE: RESOLVE ABOVE ISSUE ####
# In OpenRefine, check the similarly names lakes that are in LW


#### visualize ####

# counts per data type
lakes %>%
  group_by(DataType) %>%
  count() %>%
  ggplot(aes(x = DataType, y = n)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust = -0.25, size = 3) +
    xlab("Data available") +
    ylab("Water bodies") +
    theme_bw()




