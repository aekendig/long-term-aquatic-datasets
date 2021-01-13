#### info ####

# goal: visualize surface area data


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
fwc_plant <- read_csv("original-data/FWC Plant Surveys.csv")
fwc_plant_new <- read_csv("original-data/FWCSurvey_2018_2019.csv")
vol <- read_csv("original-data/Volume_Calculation.csv")
hydro <- read_csv("original-data/Florida_Lakes_National_Hydrography_Dataset.csv") # not sure what units this is in, likely m^2
gnis <- read_csv("original-data/LW_matching_Herbicide_lakes_with_GNIS.csv")


#### edit data ####

# gnis data
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

# edit plant data (modified from combining_datasets.R)
fwc_plant2 <- fwc_plant %>%
  mutate(WaterbodyName = ifelse(WaterbodyName == "Watermellon Pond", "Watermelon Pond", WaterbodyName)) %>%
  full_join(fwc_plant_new) %>%
  filter(!is.na(SpeciesName)) %>%
  mutate(County_FWC = toupper(County) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")),
         Lake_FWC = str_replace_all(WaterbodyName, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper(),
         Lake_FWC = case_when(Lake_FWC == "HALF MOON" & County_FWC == "MARION" ~ "HALFMOON",
                              Lake_FWC == "LITTLE RED WATER" & County_FWC == "HIGHLANDS" ~ "LITTLE REDWATER",
                              TRUE ~ Lake_FWC),
         WaterbodyMeters = WaterbodyAcres*4046.86) %>%
  group_by(Lake_FWC, County_FWC) %>%
  summarise(WaterbodyAcres = mean(WaterbodyAcres),
            WaterbodyMeters = mean(WaterbodyMeters)) %>%
  left_join(gnis2) %>%
  mutate(Lake_LW = ifelse(is.na(Lake_LW), Lake_FWC, Lake_LW),
         County_LW = ifelse(is.na(County_LW), County_FWC, County_LW))

# volume data
vol2 <- vol %>%
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
                             TRUE ~ Lake_LW),
         Surface_area_meters = Surface_area_acres*4046.86) %>%
  group_by(County_LW, Lake_LW) %>%
  summarise(SurfaceAreaAcres = mean(Surface_area_acres),
            SurfaceAreaMeters = mean(Surface_area_meters)) %>%
  ungroup()

# area data
hydro2 <- hydro %>%
  mutate(County_LW = toupper(COUNTY) %>%
           str_replace_all(c("ST." = "SAINT", "ST"= "SAINT")),
         Lake_LW = str_replace_all(NAME, c(", Lake" = "", "Lake " = "", " Lake" = "", "'" = "")) %>%
           toupper())

# check area data for duplicates
hydro2 %>%
  group_by(County_LW, Lake_LW) %>%
  summarise(lakes = n(),
            names = paste(NAME, collapse = ", "),
            counties = paste(COUNTY, collapse = ", "),
            IDs = paste(OBJECTID, collapse = ", "),
            areas = paste(SHAPEAREA, collapse = ", ")) %>%
  ungroup() %>%
  filter(lakes > 1)
# 300 cases, the ones shown have the same exact name
#### start here: emailed FDEP ####

# combine data
area <- full_join(hydro2, fwc_plant2) %>%
  full_join(vol2)


#### compare area estimates ####

# hydro data and FWC (n = 371)
hydro_fwc <- area %>%
  filter(!is.na(SHAPEAREA) & !is.na(WaterbodyMeters)) %>%
  select(County_LW, Lake_LW, SHAPEAREA, WaterbodyMeters, WaterbodyAcres) %>%
  unique() %>%
  mutate(meters_diff = (SHAPEAREA - WaterbodyMeters)/WaterbodyMeters,
         acres_diff = (SHAPEAREA - WaterbodyAcres)/WaterbodyAcres)

hydro_fwc %>%
  ggplot(aes(x = meters_diff)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10()

hydro_fwc %>%
  ggplot(aes(x = acres_diff)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10()

median(hydro_fwc$meters_diff)
median(hydro_fwc$acres_diff)
# more are around zero for meters

# hydro data and volume data (n = 351)
hydro_vol <- area %>%
  filter(!is.na(SHAPEAREA) & !is.na(SurfaceAreaMeters)) %>%
  select(County_LW, Lake_LW, SHAPEAREA, SurfaceAreaMeters, SurfaceAreaAcres) %>%
  unique() %>%
  mutate(meters_diff = (SHAPEAREA - SurfaceAreaMeters)/SurfaceAreaMeters,
         acres_diff = (SHAPEAREA - SurfaceAreaAcres)/SurfaceAreaAcres)

hydro_vol %>%
  ggplot(aes(x = meters_diff)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10()

hydro_vol %>%
  ggplot(aes(x = acres_diff)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10()

median(hydro_vol$meters_diff)
median(hydro_vol$acres_diff)
# more are similar for meters than acres

# FWC and volume data (163)
fwc_vol <- area %>%
  filter(!is.na(WaterbodyMeters) & !is.na(SurfaceAreaMeters)) %>%
  select(County_LW, Lake_LW, SHAPEAREA, WaterbodyMeters, SurfaceAreaMeters, WaterbodyAcres, SurfaceAreaAcres) %>%
  unique() %>%
  mutate(meters_diff = (WaterbodyMeters - SurfaceAreaMeters)/SurfaceAreaMeters,
         acres_diff = (WaterbodyAcres - SurfaceAreaAcres)/SurfaceAreaAcres)

fwc_vol %>%
  ggplot(aes(x = meters_diff)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10()

fwc_vol %>%
  ggplot(aes(x = acres_diff)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_log10()

median(fwc_vol$meters_diff)

# combine comparisons
comp <- hydro_fwc %>%
  select(-c(acres_diff, WaterbodyAcres)) %>%
  rename("hydro_fwc" = "meters_diff") %>%
  full_join(hydro_vol %>%
              select(-c(acres_diff, SurfaceAreaAcres)) %>%
              rename("hydro_vol" = "meters_diff")) %>%
  full_join(fwc_vol %>%
              select(-c(acres_diff, WaterbodyAcres, SurfaceAreaAcres)) %>%
              rename("fwc_vol" = "meters_diff"))
  

# multiple comparisons?
comp %>%
  filter(!is.na(hydro_fwc) & !is.na(hydro_vol)) # 149

comp %>%
  filter(!is.na(hydro_fwc) & !is.na(fwc_vol)) # 149

comp %>%
  filter(!is.na(hydro_vol) & !is.na(fwc_vol)) # 149

# hydro_fwc
comp %>%
  filter(abs(hydro_fwc) > 0.1 & (abs(hydro_vol) > 0.1 | abs(fwc_vol) > 0.1)) %>%
  data.frame()
# there are repeat lakes with different surface areas according to the hydro data