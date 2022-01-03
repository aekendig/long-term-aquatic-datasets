#### set-up ####

#### start here ####
# QA codes for chnep. Not sure how to resolve these: different sources, multiple pasted together

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)

# import data
gis <- read_csv("gis/intermediate-data/FAOL_Lakewatch_FWC.csv")
chnep <- read_csv("original-data/water_atlas_CHNEP_121321.csv")
hillsborough <- read_csv("original-data/water_atlas_Hillsborough_121321.csv")
lake <- read_csv("original-data/water_atlas_Lake_121321.csv")
manatee <- read_csv("original-data/water_atlas_Manatee_121321.csv")
orange1 <- read_csv("original-data/water_atlas_Orange1_121321.csv")
orange2 <- read_csv("original-data/water_atlas_Orange2_121321.csv")
pinellas <- read_csv("original-data/water_atlas_Pinellas_121321.csv") # manually cleaned to remove running rows
polk <- read_csv("original-data/water_atlas_Polk_121321.csv")
sarasota <- read_csv("original-data/water_atlas_Sarasota_121321.csv")
seminole <- read_csv("original-data/water_atlas_Seminole_121321.csv")
tampa <- read_csv("original-data/water_atlas_Tampa_Bay_Estuary_121321.csv") # manually cleaned to remove running rows


#### edit GIS ####

# rename GIS columns
gis2 <- gis %>%
  rename(WBodyID = FAOL_WBODY,
         PermanentID = Permanent_) %>%
  select(WBodyID, PermanentID)


#### edit CHNEP ####

# CHNEP: check characters that should be numbers
chnep %>%
  mutate(WBodyID2 = as.numeric(WBodyID)) %>%
  filter(is.na(WBodyID2)) %>%
  select(WBodyID, WBodyID2) %>%
  unique()
# WBodyID column was overwritten with notes
# the data seems generally messed up -- don't include

chnep %>%
  mutate(Original_Result_Value2 = as.numeric(Original_Result_Value)) %>%
  filter(is.na(Original_Result_Value2)) %>%
  select(Original_Result_Value) %>%
  unique()
# these are NA values anyway

# QA codes
unique(chnep$QACode)

# CHNEP: turn characters to numeric
# edit date
chnep2 <- chnep %>%
  mutate(WBodyID = as.numeric(WBodyID),
         Original_Result_Value = as.numeric(Original_Result_Value),
         DateTime = as_datetime(SampleDate, format = "%m/%d/%y %H:%M", tz = "America/New_York"),
         Year = if_else(DateTime > as_date("2021-12-13"), # remove future dates
                        year(DateTime) - 100,
                        year(DateTime)),
         Date = as_date(paste0(Year, "-", month(DateTime), "-", day(DateTime)), format = "%Y-%m-%d")) %>%
  select(-DateTime) %>%
  filter(!is.na(WBodyID))

# check times
ggplot(chnep2, aes(x = Year)) +
  geom_histogram(binwidth = 1)


#### edit hillsborough ####

hillsborough




# Manatee: check characters that should be numbers
manatee %>%
  mutate(Original_Result_Value2 = as.numeric(Original_Result_Value)) %>%
  filter(is.na(Original_Result_Value2)) %>%
  select(Original_Result_Value) %>%
  unique()
# these are NA values anyway

# Manatee: turn characters to numeric
manatee2 <- manatee %>%
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
