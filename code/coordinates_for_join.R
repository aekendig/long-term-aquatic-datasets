#### info ####

# goal: visualize management activity


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(fst)

# import data
ctrl <- read_csv("original-data/FWC_Herbicide_Treatments_mid2010-Oct2020.csv")
ctrl_old <- read_fst("original-data/PrePMARS_IPMData.fst")
fwc_plant <- read_csv("original-data/FWC Plant Surveys.csv")
fwc_plant_new <- read_csv("original-data/FWCSurvey_2018_2019.csv")
lw_base <- read_csv("original-data/Lakewatch_Base_File_for_Amy_2020.csv")
lw_plant <- read_csv("original-data/Lakewatch_Plant_Surveys.csv")
qual <- read_csv("original-data/Lakewatch_1986-2019_TP_TN_Chl_Secchi_Color_Cond.csv")


#### edit FWC data ####

# select columns from control data
ctrl_gis <- ctrl %>%
  select(AreaOfInterestID, longitude, latitude, AreaOfInterest) %>%
  unique()

ctrl_old_gis <- ctrl_old %>%
  select(AreaOfInterestID, longitude, latitude, AreaOfInterest) %>%
  unique()

# combine old and new control data
fwc_gis <- full_join(ctrl_gis, ctrl_old_gis)

# missing coordinates
fwc_gis %>%
  filter(is.na(latitude) | is.na(longitude))
# admin areas - not needed

# check for duplicate IDs
dupIDs <- fwc_gis[duplicated(fwc_gis$AreaOfInterestID) == T, "AreaOfInterestID"]
fwc_gis %>%
  filter(AreaOfInterestID %in% dupIDs$AreaOfInterestID) %>%
  arrange(AreaOfInterestID)
# some have different coordinates, some have different names
# make sure they are assigned the same waterbody

# no lat/long for plant data - needed?
fwc_missing_gis <- fwc_plant %>%
  rename(AreaOfInterest = WaterbodyName) %>%
  select(AreaOfInterest, County) %>%
  unique() %>%
  full_join(fwc_plant_new %>%
              rename(AreaOfInterest = WaterbodyName) %>%
              select(AreaOfInterestID, AreaOfInterest, County) %>%
              unique()) %>%
  anti_join(fwc_gis)
# 98 out of 497 missing coordinates

fwc_missing_gis %>%
  filter(is.na(AreaOfInterestID))
# 32 missing IDs


#### edit LW data ####

# see how many coords per LW lake
lw_base %>%
  group_by(County, Lake) %>%
  summarise(lats = length(unique(round(Latitude, 1))),
            longs = length(unique(round(Longitude, 1)))) %>%
  ungroup() %>%
  filter(lats > 1 | longs > 1)
# lots

# missing coordinates
lw_missing_gis <- lw_base %>%
  filter(is.na(Latitude) | is.na(Longitude)) %>%
  select(County, Lake, GNIS_ID) %>%
  unique()
# some have GNIS_ID - use to merge

lw_base %>%
  filter(!is.numeric(Latitude) | !is.numeric(Longitude))

# summarise
lw_gis <- lw_base %>%
  group_by(County, Lake, GNIS_ID, Surface_area_hectares) %>%
  summarise(avg_station_latitude = mean(Latitude, na.rm = T),
            avg_station_longitude = mean(Longitude, na.rm = T))

# duplicate lake/county combos?
lw_gis %>%
  group_by(County, Lake) %>%
  summarise(IDs = length(unique(GNIS_ID)),
           areas = length(unique(Surface_area_hectares))) %>%
  filter(IDs > 1 | areas > 1) %>%
  inner_join(lw_base) %>%
  arrange(County, Lake) %>%
  select(County, Lake, GNIS_ID, Station, Station_ID, Latitude, 
         Longitude, Surface_area_hectares, Surface_area_citation) %>%
  data.frame()
# 7 lake/county combos
# Bro-Seneca-1 and 2 are in a different section of the lake than 3 (different GNIS ID)
# one set of stations missing GNIS ID or surface area data

# missing from qual/plant datasets?
qual %>%
  select(County, Lake) %>%
  unique() %>%
  left_join(lw_gis %>%
              mutate(gis = "yes")) %>%
  filter(is.na(gis))
# none missing

lw_plant_missing <- lw_plant %>%
  select(County, Lake) %>%
  unique() %>%
  left_join(lw_gis %>%
              mutate(gis = "yes")) %>%
  filter(is.na(gis))
# 20 lakes - some spelling differences, some missing


#### new FWC coordinates ####

# import coordinates
fwc_coord <- read_csv("original-data/FWC_coordinates_full.csv")

# add to the missing dataset
fwc_replaced_gis <- fwc_missing_gis %>%
  filter(!is.na(AreaOfInterestID)) %>%
  left_join(fwc_coord %>%
              rename(AreaOfInterest = Name)) %>%
  full_join(fwc_missing_gis %>%
              filter(is.na(AreaOfInterestID)) %>%
              select(-AreaOfInterestID) %>%
              left_join(fwc_coord %>%
                          rename(AreaOfInterest = Name)))

fwc_replaced_gis %>%
  filter(is.na(Longitude))
# these are in the dataset with coordinates, just mispellings


#### check coordinates (gis mismatches) ####

# Little Sawyer Lake
lw_base %>%
  filter(Lake == "Little Sawyer" & County == "Orange") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()
# type-o in the second longitude

lw_base %>%
  filter(Lake == "Little Sawyer" & County == "Orange" & Station != 2) %>%
  summarise(latitude = mean(Latitude),
            longitude = mean(Longitude)) %>%
  data.frame()

# Riley Lake
lw_base %>%
  filter(Lake == "Riley" & County == "Putnam") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()
# station 3 longitude is off

lw_base %>%
  filter(Lake == "Riley" & County == "Putnam" & Station != 3) %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame() %>%
  summarise(latitude = mean(Latitude),
            longitude = mean(Longitude)) %>%
  data.frame()

# Tojo Lake
lw_base %>%
  filter(Lake == "Tojo" & County == "Manatee") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()
# station 2 latitude is off

lw_base %>%
  filter(Lake == "Tojo" & County == "Manatee" & Station != 2) %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame() %>%
  summarise(latitude = mean(Latitude),
            longitude = mean(Longitude)) %>%
  data.frame()

# Banyan
lw_base %>%
  filter(Lake == "Banyan" & County == "Palm Beach") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Bear2
lw_base %>%
  filter(Lake == "Bear 2" & County == "Lake") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Bream
lw_base %>%
  filter(Lake == "Bream" & County == "Washington") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Eel
lw_base %>%
  filter(Lake == "Eel" & County == "Leon") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# English
lw_base %>%
  filter(Lake == "English" & County == "Putnam") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Estancia
lw_base %>%
  filter(Lake == "Estancia" & County == "Broward") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Fisher
lw_base %>%
  filter(Lake == "Fisher" & County == "Leon") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

lw_base %>%
  filter(Lake == "Fisher" & County == "Leon" & Station != 1) %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame() %>%
  summarise(latitude = mean(Latitude),
            longitude = mean(Longitude)) %>%
  data.frame()

# Gannett Pond
lw_base %>%
  filter(Lake == "Gannett Pond" & County == "Leon") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Glade
lw_base %>%
  filter(Lake == "Glade" & County == "Miami-Dade") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Griffin West
lw_base %>%
  filter(Lake == "Griffin West" & County == "Lake") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Hammock
lw_base %>%
  filter(Lake == "Hammock" & County == "Hillsborough") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()
# station 2 seems off

lw_base %>%
  filter(Lake == "Hammock" & County == "Hillsborough" & Station != 2) %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame() %>%
  summarise(latitude = mean(Latitude),
            longitude = mean(Longitude)) %>%
  data.frame()

# Harmony Estates-Retention Pond
lw_base %>%
  filter(Lake == "Harmony Estates-Retention Pond") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()
# station 2 coordinates are not a pond

# Helen
lw_base %>%
  filter(Lake == "Helen" & County == "Broward") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Heron
lw_base %>%
  filter(Lake == "Heron" & County == "Lake") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()
# coordinates off a little

# Holly
lw_base %>%
  filter(Lake == "Holly" & County == "Hillsborough") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Hubbert
lw_base %>%
  filter(Lake == "Hubbert" & County == "Orange") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Ironwood
lw_base %>%
  filter(Lake == "Ironwood B") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Johnson Pond
lw_base %>%
  filter(Lake == "Johnson Pond" & County == "Alachua") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Lee
lw_base %>%
  filter(Lake == "Lee" & County == "Orange") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()
# station 3 coordinate is not a lake

lw_base %>%
  filter(Lake == "Lee" & County == "Orange" & Station != 3) %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame() %>%
  summarise(latitude = mean(Latitude),
            longitude = mean(Longitude)) %>%
  data.frame()

# Lost
lw_base %>%
  filter(Lake == "Lost" & County == "Polk") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Lucy
lw_base %>%
  filter(Lake == "Lucy" & County == "Lake") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()
# third coordinate is off

lw_base %>%
  filter(Lake == "Lucy" & County == "Lake" & Station != 3) %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame() %>%
  summarise(latitude = mean(Latitude),
            longitude = mean(Longitude)) %>%
  data.frame()

# Mill Pond
lw_base %>%
  filter(Lake == "Mill Pond" & County == "Marion") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Offspring
lw_base %>%
  filter(Lake == "Offspring" & County == "Volusia") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()

# Parkway Meadows
lw_base %>%
  filter(Lake == "Parkway Meadows" & County == "Brevard") %>%
  select(Station_ID, Latitude, Longitude) %>%
  data.frame()


#### output ####
write_csv(fwc_gis, "gis/data/FWC_Herbicide_Coordinates.csv")
write_csv(fwc_missing_gis, "intermediate-data/FWC_Missing_Coordinates.csv")
write_csv(fwc_replaced_gis %>% filter(!is.na(Longitude)), "gis/data/FWC_Replaced_Coordinates.csv")

write_csv(lw_gis, "gis/data/Lakewatch_Coordinates.csv")
write_csv(lw_missing_gis, "intermediate-data/Lakewatch_Missing_Coordinates.csv")
write_csv(lw_plant_missing, "intermediate-data/Lakewatch_Plant_Missing_Coordinates.csv")
