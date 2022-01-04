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
qa <- read_csv("original-data/water_atlas_qa_codes.csv")


#### edit GIS ####

# rename GIS columns
gis2 <- gis %>%
  rename(WBodyID = FAOL_WBODY,
         PermanentID = Permanent_) %>%
  select(WBodyID, PermanentID)


#### edit CHNEP ####

chnep

# check characters that should be numbers
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

# turn characters to numeric
chnep2 <- chnep %>%
  mutate(WBodyID = as.numeric(WBodyID),
         Original_Result_Value = as.numeric(Original_Result_Value)) %>%
  filter(!is.na(WBodyID))


#### edit hillsborough ####

hillsborough


#### edit lake ####

lake


#### edit manatee ####

# check characters that should be numbers
manatee %>%
  mutate(Original_Result_Value2 = as.numeric(Original_Result_Value)) %>%
  filter(is.na(Original_Result_Value2)) %>%
  select(Original_Result_Value) %>%
  unique()
# these are NA values anyway

# turn characters to numeric
manatee2 <- manatee %>%
  mutate(Original_Result_Value = as.numeric(Original_Result_Value))


#### edit orange ####

orange1
orange2


#### edit pinellas ####

pinellas


#### edit polk ####

polk


#### edit sarasota ####

sarasota


#### edit seminole ####

seminole

# check characters that should be numbers
seminole %>%
  mutate(Original_Result_Value2 = as.numeric(Original_Result_Value)) %>%
  filter(is.na(Original_Result_Value2)) %>%
  select(Original_Result_Value) %>%
  unique()
# these are NA values anyway

# turn characters to numeric
seminole2 <- seminole %>%
  mutate(Original_Result_Value = as.numeric(Original_Result_Value))


#### edit tampa ####

tampa


#### combine ####

# Tampa, Pinellas, Manatee, and Sarasota should be included in Hillsborough and CHNEP
# https://wateratlas.usf.edu/about/
counties <- tampa %>%
  full_join(pinellas) %>%
  full_join(manatee2) %>%
  full_join(sarasota)

regions <- chnep2 %>%
  full_join(hillsborough)

counties %>%
  anti_join(regions)
# many rows not included

# combine data
atlas <- chnep2 %>%
  full_join(hillsborough) %>%
  full_join(lake) %>%
  full_join(manatee2) %>%
  full_join(orange1) %>%
  full_join(orange2) %>%
  full_join(pinellas) %>%
  full_join(polk) %>%
  full_join(sarasota) %>%
  full_join(seminole2) %>%
  full_join(tampa) %>%
  inner_join(gis2)

# lakes included
length(unique(atlas$PermanentID))
# 849


#### edit combined data ####

# edit date
# remove missing data
atlas2 <- atlas %>%
  mutate(DateTime = as_datetime(SampleDate, format = "%m/%d/%y %H:%M", tz = "America/New_York"),
         Year = if_else(DateTime > as_date("2021-12-13"), # revise future dates
                        year(DateTime) - 100,
                        year(DateTime)),
         Date = as_date(paste(Year, month(DateTime), day(DateTime), sep = "-"), format = "%Y-%m-%d")) %>%
  select(-DateTime) %>%
  filter(!is.na(Result_Value))

# check years
ggplot(atlas2, aes(x = Year)) +
  geom_histogram(binwidth = 1)

# qa codes with/out datasource
qa_ds <- qa %>%
  filter(!is.na(DataSource))

qa_nods <- qa %>%
  filter(is.na(DataSource)) %>%
  select(-DataSource)

# QA codes
atlas_qa <- atlas2 %>%
  filter(!is.na(QACode) & QACode != "NULL") %>%
  select(DataSource, QACode) %>%
  unique() %>%
  left_join(qa_ds) %>% # match codes with defined data sources
  rename(DSFull = QAMeaning) %>%
  mutate(Chars = nchar(QACode), # count characters
         QACode1 = str_sub(QACode, 1, 1),
         QACode2 = str_sub(QACode, 2, 2),
         QACode3 = str_sub(QACode, 3, 3),
         QACode4 = str_sub(QACode, 4, 4),
         QACode2 = if_else(QACode1 == "J" & str_detect(QACode2, "[0-9]") == T, 
                           paste0("J", QACode2), QACode2),
         QACode3 = if_else(QACode2 == "J" & str_detect(QACode3, "[0-9]") == T, 
                           paste0("J", QACode3), QACode3),
         QACode4 = if_else(QACode3 == "J" & str_detect(QACode4, "[0-9]") == T, 
                           paste0("J", QACode4), QACode4),
         QACode1 = if_else(QACode1 == "J" & QACode2 %in% c("J1", "J2", "J3", "J4", "J5"),
                           "", QACode1),
         QACode2 = if_else(QACode2 == "J" & QACode3 %in% c("J1", "J2", "J3", "J4", "J5"),
                           "", QACode2),
         QACode3 = if_else(QACode3 == "J" & QACode4 %in% c("J1", "J2", "J3", "J4", "J5"),
                           "", QACode3)) %>%
  left_join(qa_ds %>%
              rename(QACode1 = QACode,
                     DS1 = QAMeaning)) %>%
  left_join(qa_ds %>%
              rename(QACode2 = QACode,
                     DS2 = QAMeaning)) %>%
  left_join(qa_ds %>%
              rename(QACode3 = QACode,
                     DS3 = QAMeaning)) %>%
  left_join(qa_ds %>%
              rename(QACode4 = QACode,
                     DS4 = QAMeaning)) %>%
  left_join(qa_nods %>%
              rename(QACode1 = QACode,
                     NoDS1 = QAMeaning)) %>%
  left_join(qa_nods %>%
              rename(QACode2 = QACode,
                     NoDS2 = QAMeaning)) %>%
  left_join(qa_nods %>%
              rename(QACode3 = QACode,
                     NoDS3 = QAMeaning)) %>%
  left_join(qa_nods %>%
              rename(QACode4 = QACode,
                     NoDS4 = QAMeaning)) %>%
  mutate(QAM1 = if_else(!is.na(DS1), DS1, NoDS1),
         QAM2 = if_else(!is.na(DS2), DS2, NoDS2),
         QAM3 = if_else(!is.na(DS3), DS3, NoDS3),
         QAM4 = if_else(!is.na(DS4), DS4, NoDS4)) %>%
  rowwise() %>%
  mutate(QAM = paste(QAM1, QAM2, QAM3, QAM4)) %>%
  ungroup() %>%
  mutate(QAM = str_replace_all(QAM, "NA", ""),
         QAMeaning = if_else(!is.na(DSFull), DSFull, QAM))

#### start here: resolve below, if possible ####
# check missing values
atlas_qa %>%
  filter((is.na(QAM1) & QACode1 != "")) %>%
  select(DataSource, QACode, QACode1, QAM1) %>%
  data.frame()

atlas_qa %>%
  filter((is.na(QAM2) & QACode2 != "")) %>%
  select(DataSource, QACode, QACode2, QAM2) %>%
  data.frame()

atlas_qa %>%
  filter((is.na(QAM3) & QACode3 != "")) %>%
  select(DataSource, QACode, QACode3, QAM3) %>%
  data.frame()

atlas_qa %>%
  filter((is.na(QAM4) & QACode4 != "")) %>%
  select(DataSource, QACode, QACode4, QAM4) %>%
  data.frame()
