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
library(readxl)

# import data
ctrl_old <- read_csv("original-data/PrePMARS_IPMData.csv")
ctrl <- read_csv("original-data/FWC_Herbicide_Treatments_mid2010-Oct2020.csv")
fwc_plant <- read_csv("original-data/FWC Plant Surveys.csv")
fwc_plant_new <- read_csv("original-data/FWCSurvey_2018_2019.csv")
fwc_ID <- read_csv("gis/data/FWC_Replaced_Coordinates.csv")

qual <- read_csv("original-data/Lakewatch_1986-2019_TP_TN_Chl_Secchi_Color_Cond.csv")
lw_plant <- read_csv("original-data/Lakewatch_Plant_Surveys.csv")
qual_base <- read_csv("original-data/Lakewatch_Base_File_for_Amy_2020.csv")
lw_plant_names <- read_csv("intermediate-data/Lakewatch_Plant_Base_Lake_Names.csv")
# made this manually based on intermediate-data/Lakewatch_Plant_Missing_Coordinates.csv
lw_missing_gis <- read_csv("intermediate-data/Lakewatch_Missing_Coordinates.csv")

gis_ed <- read_csv("gis/intermediate-data/Lakewatch_FWC_Waterbody_merge_join_edited.csv")
gis_man <- read_csv("gis/intermediate-data/Lakewatch_FWC_Waterbody_merge_join_manual.csv")

# gnis <- read_csv("original-data/LW_matching_Herbicide_lakes_with_GNIS.csv")
# depth <- read_csv("original-data/Mean_Depth_Data.csv")
# vol <- read_csv("original-data/Volume_Calculation.csv")


#### gis data ####

# look at notes
unique(gis_ed$JoinNotes)
unique(gis_man$JoinNotes)

# remove lakes without shapes
# add manually matched lakes back in
# convert shape area to km squared
gis <- gis_ed %>%
  filter(!is.na(PermanentID)) %>%
  full_join(gis_man) %>%
  mutate(ShapeArea = ifelse(ShapeSource == "FDEP", ShapeArea*10^-6, ShapeArea*10^4),
         County = toupper(County)) %>%
  filter(JoinNotes != "outer boundary" | is.na(JoinNotes))


#### old control data ####

# rename columns
# add permanent ID based on AreaOfInterestID
# remove missing ID
# remove duplicate rows
ctrl_old2 <- ctrl_old %>%
  rename("Longitude_FWC" = "longitude",
         "Latitude_FWC" = "latitude",
         "SpeciesOrig" = "Species_orig",
         "TotalContFWC" = "Total_cont_fwc",
         "County_FWC" = "County") %>%
  mutate(County_FWC = toupper(County_FWC)) %>%
  left_join(gis %>%
              filter(CoordSource %in% c("FWC", "FWC_2")) %>%
              select(AreaOfInterestID, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique()) %>%
  filter(!is.na(PermanentID)) %>%
  unique()


#### control data ####

# format date
# rename columns
# add permanent ID based on AreaOfInterestID
# remove missing ID
# remove duplicate rows
ctrl2 <- ctrl %>%
  mutate(BeginDate = as.Date(BeginDate, "%m/%d/%y"),
         County = toupper(County)) %>%
  rename("Year" = "year",
         "Month" = "month",
         "Species" = "species",
         "Herbicide" = "herbicide",
         "TotalHerbicideUsed" = "totalherbicideused",
         "Longitude_FWC" = "longitude",
         "Latitude_FWC" = "latitude",
         "County_FWC" = "County")  %>%
  left_join(gis %>%
              filter(CoordSource %in% c("FWC", "FWC_2")) %>%
              select(AreaOfInterestID, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique()) %>%
  filter(!is.na(PermanentID)) %>%
  unique()


#### new FWC plant survey data ####

# format date
# rename columns
# add permanent ID based on AreaOfInterestID
# remove missing ID
# remove duplicate rows
fwc_plant_new2 <- fwc_plant_new %>%
  mutate(SurveyDate = as.Date(SurveyDate, "%m/%d/%y"),
         County = toupper(County)) %>%
  rename("AreaOfInterest" = "WaterbodyName") %>%
  left_join(gis %>%
              filter(CoordSource %in% c("FWC", "FWC_2")) %>%
              select(AreaOfInterestID, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique()) %>%
  filter(!is.na(PermanentID)) %>%
  unique()


#### FWC plant survey data ####

# this dataset is missing AOI ID's
# get ID's from other FWC datasets
fwc_lakes <- ctrl_old %>%
  select(AreaOfInterest, County, AreaOfInterestID) %>%
  unique() %>%
  mutate(County = toupper(County)) %>%
  full_join(ctrl %>%
              select(AreaOfInterest, County, AreaOfInterestID) %>%
              unique() %>%
              mutate(County = toupper(County))) %>%
  full_join(fwc_plant_new %>%
              rename("AreaOfInterest" = "WaterbodyName") %>%
              select(AreaOfInterest, County, AreaOfInterestID) %>%
              unique() %>%
              mutate(County = toupper(County))) %>%
  full_join(fwc_ID %>%
              select(AreaOfInterest, County, AreaOfInterestID) %>%
              unique() %>%
              mutate(County = toupper(County)) %>%
              filter(!(AreaOfInterestID %in% c(714, 44, 872, 218)))) %>%
  filter()
# the four removed are duplicate AreaOfInterest names with different IDs
# Alex Dew advised on which should be removed for survey data
# previous to 6/24/21, I had removed 219 instead of 218 and 45 instead of 44, which was incorrect
# 218 is another lake with the same name as 219 in a different county
# 44 is another lake with the same name as 45 in a different county
# so analyses prior to 6/24/21 combined these two lakes

# duplicate Lake/County combos
fwc_lakes %>%
  group_by(AreaOfInterest, County) %>%
  mutate(IDs = length(unique(AreaOfInterestID))) %>%
  ungroup() %>%
  filter(IDs > 1)
# none

# missing IDs
fwc_lakes %>%
  filter(is.na(AreaOfInterestID))
# none

# matches in plant surveys
fwc_plant %>%
  rename("AreaOfInterest" = "WaterbodyName") %>%
  select(AreaOfInterest, County) %>%
  unique() %>%
  mutate(County = toupper(County)) %>%
  anti_join(fwc_lakes)
# all are matched

# add ID's
# format dates
# add permanent ID based on AreaOfInterestID
# remove missing ID
# remove duplicate rows
fwc_plant2 <- fwc_plant %>%
  mutate(SurveyDate = as.Date(SurveyDate, "%m/%d/%Y"),
         AreaOfInterest = WaterbodyName,
         County = toupper(County)) %>%
  left_join(fwc_lakes) %>%
  left_join(gis %>%
              filter(CoordSource %in% c("FWC", "FWC_2")) %>%
              select(AreaOfInterestID, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique()) %>%
  filter(!is.na(PermanentID)) %>%
  unique()


# combine plant surveys (no overlapping rows)
# remove rows with missing species names
# sum cover for duplicate non-species
# use max cover for duplicate species
# add columns
fwc_plant3 <- fwc_plant2 %>%
  full_join(fwc_plant_new2) %>%
  filter(!is.na(SpeciesName)) %>%
  group_by(AreaOfInterestID, WaterbodyAcres, WaterbodyType, County, WMD, SurveyYear, SurveyDate, IsUnableToSurvey, Surveyor, IsDetected, IsAcreageRequire, SpeciesName, Origin, Eppc, Habitat, HabitatShortName) %>%
  mutate(SpeciesAcres = case_when(str_detect(SpeciesName, "Filamentous algae|other|spp.|/") == T ~ sum(SpeciesAcres, na.rm = T),
                                  TRUE ~ max(SpeciesAcres, na.rm = T))) %>%
  ungroup() %>%
  unique() %>%
  mutate(SpeciesAcres = ifelse(SpeciesAcres == -Inf, NA_real_, SpeciesAcres),
         SpeciesFrequency = SpeciesAcres / WaterbodyAcres) %>%
  rename(County_FWC = County)
# will give warnings, but these are addressed by making -Inf into NAs

# check area from FWC against shapes
fwc_plant3 %>%
  select(AreaOfInterest, AreaOfInterestID, WaterbodyAcres, ShapeArea) %>%
  unique() %>%
  mutate(Waterbody_sqKm = WaterbodyAcres * 0.004) %>%
  ggplot(aes(x = Waterbody_sqKm, y = ShapeArea)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.5)

fwc_plant3 %>%
  select(AreaOfInterest, County_FWC, AreaOfInterestID, ShapeArea, WaterbodyAcres) %>%
  unique() %>%
  mutate(Waterbody_sqKm = WaterbodyAcres * 0.004,
         AreaDev = abs(ShapeArea/Waterbody_sqKm)) %>%
  filter(AreaDev > 1.5) %>%
  arrange(County_FWC, AreaOfInterest) %>%
  data.frame()
# many of these are because the named waterbody is attached to a larger one

filter(fwc_plant3, WaterbodyAcres == max(WaterbodyAcres))
# Okeechobee


#### water quality data ####

# format date
# make secchi disk data numeric
# add permanent ID based on Lake name and county
# remove missing ID (all have been checked, include rivers, creeks, etc.)
# remove duplicate rows
qual2 <- qual %>%
  mutate(Date = as.Date(Date, "%m/%d/%y"),
         SecchiCombined = ifelse(is.na(SECCHI_ft), SECCHI_2, SECCHI_ft) %>%
           tolower(),
         SecchiBottom = ifelse(str_detect(SecchiCombined, "bottom") == T, 1, 0),
         SecchiWeeds = ifelse(str_detect(SecchiCombined, "weeds") == T, 1, 0),
         Secchi = case_when(SecchiCombined == "." ~ NA_character_, 
                            SecchiCombined == "weeds" ~ NA_character_,
                            SecchiCombined == "weeds (surface)" ~ NA_character_,
                            SecchiCombined == "bottom" ~ NA_character_,
                            TRUE ~ SecchiCombined) %>%
           parse_number(),
         AreaOfInterest = Lake,
         County = toupper(County)) %>%
  rename("Secchi1_ft" = "SECCHI_ft",
         "Secchi2_ft" = "SECCHI_2") %>%
  left_join(gis %>%
              filter(CoordSource == "Lakewatch") %>%
              select(AreaOfInterest, County, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique()) %>%
  filter(!is.na(PermanentID)) %>%
  unique() %>%
  select(-AreaOfInterest) %>%
  rename(County_LW = County)


#### lakewatch plant survey data ####

# mismatched plant lake names
# multiple permanent IDs?
lw_plant_names %>%
  rename("AreaOfInterest" = "Lake_base") %>%
  left_join(gis %>%
              filter(CoordSource == "Lakewatch") %>%
              select(AreaOfInterest, County, PermanentID) %>%
              unique()) %>%
  data.frame()
# one part of Ivanhoe was assigned a different ID, this is a different lake
# use first option for each

# format date
# combine biomass values
# use the lower value of Frequency_percent and N_rows/Stations for species occurrence
# rename columns
# use max value for species with multiple rows per sample (n = 1)
# add permanent ID based on Lake name and county
# remove missing ID (all have been checked, include rivers, creeks, etc.)
# remove duplicate rows
lw_plant2 <- lw_plant %>%
  mutate(Date = as.Date(paste(Month, Day, Year, sep = "-"), "%m-%d-%Y"),
         TotalBiomass_kg_m2 = rowSums(.[c("Em_biomass_kg_m2", "Fl_biomass_kg_m2", "Sub_biomass_kg_m2")]),
         RowsStations = N_rows / Stations * 100,
         RowsStations = case_when(RowsStations - round(RowsStations) == 0.5 ~ RowsStations + 0.5,
                                   TRUE ~ round(RowsStations)),
         SpeciesFrequency = ifelse(RowsStations < round(Frequency_percent),
                                    RowsStations/100,
                                    Frequency_percent/100),
         County = toupper(County)) %>%
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
  group_by(County, Lake, Year, Month, Day, Date, Stations, EmFlZoneWidth_ft, EmBiomass_kg_m2, FlBiomass_kg_m2, SubBiomass_kg_m2, LakeDepth_m, PercentAreaCovered, PercentVolumeInhabited, CommonName, Genus, Species, GenusSpecies, TotalBiomass_kg_m2) %>%
  summarise(RowsStations = max(RowsStations),
            SpeciesFrequency = max(SpeciesFrequency),
            StationsPresent = max(StationsPresent),
            FrequencyOrig = max(FrequencyOrig)) %>%
  ungroup() %>%
  left_join(lw_plant_names %>%
              filter(!is.na(Lake_base)) %>%
              group_by(County, Lake_plant) %>%
              summarise(Lake_base = unique(Lake_base)[1]) %>%
              ungroup() %>%
              mutate(County = toupper(County)) %>%
              rename(Lake = Lake_plant)) %>%
  mutate(AreaOfInterest = ifelse(is.na(Lake_base), Lake, Lake_base)) %>%
  left_join(gis %>%
              filter(CoordSource == "Lakewatch") %>%
              select(AreaOfInterest, County, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique()) %>%
  filter(!is.na(PermanentID)) %>%
  unique() %>%
  rename(County_LW = County)


#### combine datasets ####

# check for unique permanent IDs
ctrl_old2 %>%
  group_by(AreaOfInterestID) %>%
  summarise(IDs = length(unique(PermanentID))) %>%
  filter(IDs > 1)

ctrl2 %>%
  group_by(AreaOfInterestID) %>%
  summarise(IDs = length(unique(PermanentID))) %>%
  filter(IDs > 1)

fwc_plant3 %>%
  group_by(AreaOfInterestID) %>%
  summarise(IDs = length(unique(PermanentID))) %>%
  filter(IDs > 1)

qual2 %>%
  group_by(Lake, County_LW) %>%
  summarise(IDs = length(unique(PermanentID))) %>%
  filter(IDs > 1)
filter(gis, AreaOfInterest == "Burrell" & County == "HILLSBOROUGH") %>% select(JoinNotes) # stations in separate water bodies
filter(gis, AreaOfInterest == "Eel" & County == "LEON") %>% select(JoinNotes) # stations in separate water bodies
filter(gis, AreaOfInterest == "Hubbert" & County == "ORANGE") %>% select(JoinNotes) # stations in separate water bodies
filter(gis, AreaOfInterest == "Seneca" & County == "BROWARD") %>% select(JoinNotes) # stations in separate water bodies

lw_plant2 %>%
  group_by(Lake, County_LW) %>%
  summarise(IDs = length(unique(PermanentID))) %>%
  filter(IDs > 1)

# see if LW missing gis are needed
# add alternative lake names (without cardinal directions)
lw_missing_gis2 <- lw_missing_gis %>%
  mutate(County_LW = toupper(County),
         Lake = toupper(Lake)) %>%
  full_join(tibble(County_LW = c("COLUMBIA", "FLAGLER", "SEMINOLE"),
                   Lake = c("ALLIGATOR", "BELLE AIRE", "JESUP"))) %>%
  left_join(lw_plant2 %>%
              select(County_LW, Lake) %>%
              unique() %>%
              mutate(Lake = toupper(Lake),
                     Plant = "yes")) %>%
  left_join(qual2 %>%
              select(County_LW, Lake) %>%
              unique() %>%
              mutate(Lake = toupper(Lake),
                     Quality = "yes"))

lw_missing_gis2 %>%
  inner_join(fwc_lakes %>%
               mutate(County_LW = County,
                      Lake = toupper(AreaOfInterest) %>%
                        str_replace(", LAKE", "") %>%
                        str_replace(" LAKE", "") %>%
                        str_replace("LAKE", "")) %>%
               select(-County))
# one lake matches, but it's the alternative name

lw_missing_gis2 %>%
  filter(Lake == "ALLIGATOR SOUTH")
# it doesn't have quality or plant data - must be in base set
# can ignore these missing coordinate lakes

# check for missing date data
sum(is.na(ctrl_old2$Year))
sum(is.na(ctrl2$BeginDate))
sum(is.na(fwc_plant3$SurveyDate))
sum(is.na(qual2$Date))
sum(is.na(lw_plant2$Date))

# number of times each lake was sampled for each variable
# merge by permanent
# start and end dates
lakes <- ctrl_old2  %>%
  mutate(CtrlOldStart = min(Year),
         CtrlOldEnd = max(Year)) %>%
  group_by(PermanentID, ShapeSource, ShapeArea, CtrlOldStart, CtrlOldEnd) %>%
  summarise(CtrlOld = length(unique(Year))) %>%
  ungroup() %>%
  full_join(ctrl2 %>%
              mutate(CtrlStart = min(BeginDate),
                     CtrlEnd = max(BeginDate)) %>%
              group_by(PermanentID, ShapeSource, ShapeArea, CtrlStart, CtrlEnd) %>%
              summarise(Ctrl = length(unique(BeginDate))) %>%
              ungroup()) %>%
  full_join(fwc_plant3 %>%
              group_by(PermanentID, ShapeSource, ShapeArea) %>%
              summarise(FWCPlant = length(unique(SurveyDate)),
                        FWCPlantStart = min(SurveyDate),
                        FWCPlantEnd = max(SurveyDate)) %>%
              ungroup()) %>%
  full_join(qual2 %>%
              group_by(PermanentID, ShapeSource, ShapeArea) %>%
              summarise(Qual = length(unique(Date)),
                        QualStart = min(Date),
                        QualEnd = max(Date)) %>%
              ungroup()) %>%
  full_join(lw_plant2 %>%
              group_by(PermanentID, ShapeSource, ShapeArea) %>%
              summarise(LWPlant = length(unique(Date)),
                        LWPlantStart = min(Date),
                        LWPlantEnd = max(Date)) %>%
              ungroup()) %>%
  mutate(CtrlOld = replace_na(CtrlOld, 0),
         Ctrl = replace_na(Ctrl, 0),
         FWCPlant = replace_na(FWCPlant, 0),
         Qual = replace_na(Qual, 0),
         LWPlant = replace_na(LWPlant, 0),
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
lakes %>%
  group_by(DataType) %>%
  count() %>%
  ggplot(aes(x = DataType, y = n)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust = -0.25, size = 3) +
    xlab("Data available") +
    ylab("Waterbodies") +
    theme_bw()
dev.off()

# size distribution per data type
pdf("output/combined_data_types_size_dist.pdf")
lakes %>%
  ggplot(aes(x = ShapeArea)) +
  geom_histogram() +
  geom_vline(xintercept = 1e-02, alpha = 0.7, linetype = "dashed") +
  geom_vline(xintercept = 1, alpha = 0.7, linetype = "dashed") +
  geom_vline(xintercept = 1e+02, alpha = 0.7, linetype = "dashed") +
  scale_x_log10() +
  xlab("Square km") +
  ylab("Waterbodies") +
  facet_wrap(~ DataType, scales = "free_y") +
  theme_bw()
dev.off()

# time lines
unique(lakes$CtrlOldStart)

lake_times <- lakes %>%
  rowwise() %>%
  mutate(AnyCtrlStart = ifelse(AnyCtrl > 0,"1/1/1998", NA_character_) %>%
           as.Date("%m/%d/%Y"),
         AnyCtrlEnd = CtrlEnd,
         AnyPlantStart = ifelse(AnyPlant > 0, as.character(min(c(FWCPlantStart, LWPlantStart), na.rm = T)), NA_character_) %>%
           as.Date("%Y-%m-%d"),
         AnyPlantEnd = ifelse(AnyPlant > 0, as.character(max(c(FWCPlantEnd, LWPlantEnd), na.rm = T)), NA_character_) %>%
           as.Date("%Y-%m-%d")) %>%
  ungroup() %>%
  select(PermanentID, AnyCtrlStart, AnyCtrlEnd, AnyPlantStart, AnyPlantEnd, QualStart, QualEnd) %>%
  rename(start_ctrl = AnyCtrlStart,
         end_ctrl = AnyCtrlEnd,
         start_plant = AnyPlantStart,
         end_plant = AnyPlantEnd,
         start_quality = QualStart,
         end_quality = QualEnd) %>%
  pivot_longer(cols = -PermanentID,
               names_to = c(".value", "DataType"),
               names_pattern = "(.+)_(.+)") %>%
  pivot_longer(cols = c(start, end),
               names_to = "TimePoints",
               values_to = "Date") %>%
  mutate(PermanentF = as.factor(PermanentID),
         LakeGroup = cut(as.numeric(PermanentF), breaks = 5),
         PermanentF = fct_rev(PermanentF))

lake_times_grp = sort(unique(lake_times$LakeGroup))

pdf("output/combined_data_time_series.pdf")
for(i in 1:length(lake_times_grp)){
  
  lake_times_sub <- lake_times %>%
    filter(LakeGroup == lake_times_grp[i])
  
  print(ggplot(data = lake_times_sub,
               aes(x = Date,
                   y = PermanentF,
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


#### output data ####

# FWC-only lakes
fwc_lakes2 <- ctrl_old2 %>%
  select(AreaOfInterest, AreaOfInterestID, GNISID, PermanentID, ShapeSource) %>%
  unique() %>%
  mutate(ctrl_old = 1) %>%
  full_join(ctrl2 %>%
              select(AreaOfInterest, AreaOfInterestID, GNISID, PermanentID, ShapeSource) %>%
              unique() %>%
              mutate(ctrl_new = 1)) %>%
  full_join(fwc_plant3 %>%
              select(AreaOfInterest, AreaOfInterestID, GNISID, PermanentID, ShapeSource) %>%
              unique() %>%
              mutate(plants = 1)) %>%
  mutate(ctrl_old = replace_na(ctrl_old, 0),
         ctrl_new = replace_na(ctrl_new, 0),
         plants = replace_na(plants, 0))

write_csv(ctrl_old2, "intermediate-data/FWC_control_old_formatted.csv")
write_csv(ctrl2, "intermediate-data/FWC_control_new_formatted.csv")
write_csv(fwc_plant3, "intermediate-data/FWC_plant_formatted.csv")
write_csv(qual2, "intermediate-data/LW_quality_formatted.csv")
write_csv(lw_plant2, "intermediate-data/LW_plant_formatted.csv")
write_csv(lakes, "intermediate-data/FWC_LW_combined.csv")
write_csv(fwc_lakes2, "intermediate-data/FWC_control_plants_lakes.csv")
write_csv(fwc_lakes2 %>%
            filter(ShapeSource == "FDEP"), "gis/intermediate-data/FWC_control_plants_lakes_FDEP.csv")
write_csv(fwc_lakes2 %>%
            filter(ShapeSource == "NHD"), "gis/intermediate-data/FWC_control_plants_lakes_NHD.csv")



#### old code: combine datasets hydrilla ####

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


#### old code: combine datasets floating ####

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
