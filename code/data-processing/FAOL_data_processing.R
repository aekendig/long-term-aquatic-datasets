#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
gis <- read_csv("gis/intermediate-data/FAOL_Lakewatch_FWC.csv")
chnep <- read_csv("original-data/water_atlas_CHNEP_121321.csv")
hillsborough <- read_csv("original-data/water_atlas_Hillsborough_121321.csv")
lake <- read_csv("original-data/water_atlas_Lake_121321.csv")
manatee <- read_csv("original-data/water_atlas_Manatee_121321.csv")
orange1 <- read_csv("original-data/water_atlas_Orange1_121321.csv")
orange2 <- read_csv("original-data/water_atlas_Orange2_121321.csv")
pinellas <- read_csv("original-data/water_atlas_Pinellas_121321.csv")
polk <- read_csv("original-data/water_atlas_Polk_121321.csv")
sarasota <- read_csv("original-data/water_atlas_Sarasota_121321.csv")
seminole <- read_csv("original-data/water_atlas_Seminole_121321.csv")
tampa <- read_csv("original-data/water_atlas_Tampa_Bay_Estuary_121321.csv")


#### edit data ####

# rename GIS columns
gis2 <- gis %>%
  rename(WBodyID = FAOL_WBODY,
         PermanentID = Permanent_) %>%
  select(WBodyID, PermanentID)

# CHNEP: check characters that should be numbers
CHNEP %>%
  mutate(WBodyID2 = as.numeric(WBodyID)) %>%
  filter(is.na(WBodyID2)) %>%
  select(WBodyID, WBodyID2) %>%
  unique()
# WBodyID column was overwritten with notes
# the data seems generally messed up -- don't include

CHNEP %>%
  mutate(Original_Result_Value2 = as.numeric(Original_Result_Value)) %>%
  filter(is.na(Original_Result_Value2)) %>%
  select(Original_Result_Value) %>%
  unique()
# these are NA values anyway

# CHNEP: turn characters to numeric
CHNEP2 <- CHNEP %>%
  mutate(WBodyID = as.numeric(WBodyID),
         Original_Result_Value = as.numeric(Original_Result_Value)) %>%
  filter(!is.na(WBodyID))

# Manatee: check characters that should be numbers
Manatee %>%
  mutate(Original_Result_Value2 = as.numeric(Original_Result_Value)) %>%
  filter(is.na(Original_Result_Value2)) %>%
  select(Original_Result_Value) %>%
  unique()
# these are NA values anyway

# Manatee: turn characters to numeric
Manatee2 <- Manatee %>%
  mutate(Original_Result_Value = as.numeric(Original_Result_Value))

# combine data
atlas <- CHNEP2 %>%
  full_join(Hillsborough) %>%
  full_join(Lake) %>%
  full_join(Manatee2) %>%
  inner_join(gis2)

# lakes included
length(unique(atlas$PermanentID))
# 502
