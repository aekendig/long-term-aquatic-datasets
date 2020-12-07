#### info ####

# goal: combine datasets


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
ctrl_old <- read_csv("original-data/PrePMARS_IPMData.csv")
ctrl <- read_csv("original-data/FWC_Herbicide_Treatments_mid2010-Oct2020.csv")
gnis <- read_csv("original-data/LW_matching_Herbicide_lakes_with_GNIS.csv")
qual <- read_csv("original-data/Lakewatch_1986-2019_TP_TN_Chl_Secchi_Color_Cond.csv")
lw_plant <- read_csv("original-data/Lakewatch_Plant_Surveys.csv")


#### ID data ####

# make all county information uppercase
# remove "lake" from names
gnis2 <- gnis %>%
  mutate(County_FWC = toupper(County_FWC) %>%
           str_replace("ST", "SAINT"),
         County_LW = toupper(County_LW) %>%
           str_replace("ST", "SAINT"),
         Lake_LW = str_replace_all(Lake_LW, c(", Lake" = "", "Lake " = "", " Lake" = "")) %>%
           toupper(),
         Lake_FWC = str_replace_all(AreaOfInterest, c(", Lake" = "", "Lake " = "", " Lake" = "")) %>%
           toupper()) %>%
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


#### control data ####

# county-lake combinations are specific
# see herbicide_initial_visualizations.R
# see old_herbicide_initial_visualizations.R

# add ID to old data
ctrl_old2 <- ctrl_old %>%
  filter(!is.na(AreaOfInterest)) %>%
  mutate(County_FWC = toupper(County) %>%
           str_replace("ST", "SAINT"),
         Lake_FWC = str_replace_all(AreaOfInterest, c(", Lake" = "", "Lake " = "", " Lake" = "")) %>%
           toupper()) %>%
  left_join(gnis2 %>%
              select(GNIS_ID, Lake_FWC, County_FWC)) %>%
  mutate(merge_ID = ifelse(!is.na(GNIS_ID), 
                           GNIS_ID, 
                           paste(Lake_FWC, County_FWC, sep = "_")))

# add ID to new data
ctrl2 <- ctrl %>%
  mutate(County_FWC = toupper(County) %>%
           str_replace("ST", "SAINT"),
         Lake_FWC = str_replace_all(AreaOfInterest, c(", Lake" = "", "Lake " = "", " Lake" = "")) %>%
           toupper()) %>%
  left_join(gnis2 %>%
              select(GNIS_ID, Lake_FWC, County_FWC)) %>%
  mutate(merge_ID = ifelse(!is.na(GNIS_ID), 
                           GNIS_ID, 
                           paste(Lake_FWC, County_FWC, sep = "_")))


#### water quality data ####

# no missing location data
# see quality_initial_visualizations.R

# add ID to data
qual2 <- qual %>%
  mutate(County_LW = toupper(County) %>%
           str_replace("ST", "SAINT"),
         Lake_LW = str_replace_all(Lake, c(", Lake" = "", "Lake " = "", " Lake" = "")) %>%
           toupper()) %>%
  left_join(gnis2 %>%
              select(GNIS_ID, Lake_LW, County_LW)) %>%
  mutate(merge_ID = ifelse(!is.na(GNIS_ID), 
                           GNIS_ID, 
                           paste(Lake_LW, County_LW, sep = "_")))


#### lakewatch plant survey data ####

# no missing location data
# see lakewatch_plant_initial_visualizations.R

# add ID to data
lw_plant2 <- lw_plant %>%
  mutate(County_LW = toupper(County) %>%
           str_replace("ST", "SAINT"),
         Lake_LW = str_replace_all(Lake, c(", Lake" = "", "Lake " = "", " Lake" = "")) %>%
           toupper()) %>%
  left_join(gnis2 %>%
              select(GNIS_ID, Lake_LW, County_LW)) %>%
  mutate(merge_ID = ifelse(!is.na(GNIS_ID), 
                           GNIS_ID, 
                           paste(Lake_LW, County_LW, sep = "_")))

#### start here ####
# fwc plant survey data (openrefine) and initial visualizations