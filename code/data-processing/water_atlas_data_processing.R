#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)
library(janitor)

# import data
gis <- read_csv("gis/intermediate-data/FAOL_Lakewatch_FWC_v2.csv") 
# v2 replaced previous (no _v included) because some waterbodies were accidentally omitted from previous
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
# remove duplicates
chnep2 <- chnep %>%
  mutate(WBodyID = as.numeric(WBodyID),
         Original_Result_Value = as.numeric(Original_Result_Value)) %>%
  filter(!is.na(WBodyID)) %>%
  unique()


#### edit hillsborough ####

hillsborough

# remove duplicates
hillsborough2 <- hillsborough %>%
  unique()


#### edit lake ####

lake

# remove duplicates
lake2 <- lake %>%
  unique()


#### edit manatee ####

# check characters that should be numbers
manatee %>%
  mutate(Original_Result_Value2 = as.numeric(Original_Result_Value)) %>%
  filter(is.na(Original_Result_Value2)) %>%
  select(Original_Result_Value) %>%
  unique()
# these are NA values anyway

# turn characters to numeric
# remove duplicates
manatee2 <- manatee %>%
  mutate(Original_Result_Value = as.numeric(Original_Result_Value)) %>%
  unique()


#### edit orange ####

orange1
orange2

# combine and remove duplicates
orange <- orange1 %>%
  full_join(orange2) %>%
  unique()


#### edit pinellas ####

pinellas

# remove duplicates
pinellas2 <- pinellas %>%
  unique()

#### edit polk ####

polk

# remove duplicates
polk2 <- polk %>%
  unique()


#### edit sarasota ####

sarasota

# remove duplicates
sarasota2 <- sarasota %>%
  unique()


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
# remove duplicates
seminole2 <- seminole %>%
  mutate(Original_Result_Value = as.numeric(Original_Result_Value)) %>%
  unique()


#### edit tampa ####

tampa

# remove duplicates
tampa2 <- tampa %>%
  unique()


#### combine ####

# Tampa, Pinellas, Manatee, and Sarasota should be included in Hillsborough and CHNEP
# https://wateratlas.usf.edu/about/
counties <- tampa2 %>%
  full_join(pinellas2) %>%
  full_join(manatee2) %>%
  full_join(sarasota2)

regions <- chnep2 %>%
  full_join(hillsborough2)

counties %>%
  anti_join(regions)
# many rows not included

# combine data
atlas <- chnep2 %>%
  full_join(hillsborough2) %>%
  full_join(lake2) %>%
  full_join(manatee2) %>%
  full_join(orange) %>%
  full_join(pinellas2) %>%
  full_join(polk2) %>%
  full_join(sarasota2) %>%
  full_join(seminole2) %>%
  full_join(tampa2) %>%
  inner_join(gis2) # selects for lakes with Permanent IDs

# before joining with gis2
# lakes included
n_distinct(atlas$WBodyID) # 1569

# stations per waterbody
stat_sum <- atlas %>%
  group_by(WBodyID) %>%
  summarize(nStat = n_distinct(StationID)) 

max(stat_sum$nStat)

stat_sum %>%
  ggplot(aes(x = nStat)) +
  geom_histogram(binwidth = 1)

# lakes included after joining with gis2
n_distinct(atlas$PermanentID)
# 864


#### edit dates ####

# date errors
date_err <- atlas %>%
  mutate(DateTime = as_datetime(SampleDate, format = "%m/%d/%y %H:%M", tz = "America/New_York")) %>%
  filter(is.na(DateTime)) %>%
  select(SampleDate) %>%
  unique()
# I think these are all daylight savings time errors
# switch to 2 PM

# edit date
# remove missing data
atlas2 <- atlas %>%
  mutate(SampleDate = if_else(SampleDate %in% date_err$SampleDate, 
                              str_replace(SampleDate, "2:", "14:"),
                              SampleDate),
         DateTime = as_datetime(SampleDate, format = "%m/%d/%y %H:%M", tz = "America/New_York"),
         Year = if_else(DateTime > as_date("2021-12-13"), # revise future dates
                        year(DateTime) - 100,
                        year(DateTime)),
         Date = as_date(paste(Year, month(DateTime), day(DateTime), sep = "-"), format = "%Y-%m-%d"),
         Month = month(Date)) %>%
  select(-DateTime) %>%
  filter(!is.na(Result_Value)) # only 4 rows

# year range
range(atlas2$Year)

# check years
ggplot(atlas2, aes(x = Year)) +
  geom_histogram(binwidth = 1)

# check months
ggplot(atlas2, aes(x = Month)) +
  geom_histogram(binwidth = 1)
# slightly more in spring/summer

# samples per year
date_sum <- atlas2 %>%
  group_by(WBodyID, Year) %>%
  summarize(nSamps = n_distinct(Date)) %>%
  ungroup

range(date_sum$nSamps)

date_sum %>%
  filter(nSamps <= 25) %>%
  ggplot(aes(x = nSamps)) +
  geom_histogram(binwidth = 1)


#### edit QA codes ####

# qa codes with/out datasource
qa_ds <- qa %>%
  filter(!is.na(DataSource)) %>%
  select(-Notes)

qa_nods <- qa %>%
  filter(is.na(DataSource)) %>%
  select(-c(DataSource, Notes))

# QA codes
# J# codes are multiple characters long, others are one
atlas_qa <- atlas2 %>%
  filter(!is.na(QACode) & QACode != "NULL") %>%
  select(DataSource, QACode) %>%
  unique() %>%
  left_join(qa_ds) %>% # match codes with defined data sources
  rename(DSFull = QAMeaning) %>% # full QA code text
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
                           "", QACode3),
         QACode1 = if_else(QACode1 %in% c("j", "x", "u", "l"), toupper(QACode1), QACode1)) %>%
  left_join(qa_ds %>% # match concatenated codes
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
  mutate(QAM1 = if_else(!is.na(DS1), DS1, NoDS1), # use no data source code if data source code is missing
         QAM2 = if_else(!is.na(DS2), DS2, NoDS2),
         QAM3 = if_else(!is.na(DS3), DS3, NoDS3),
         QAM4 = if_else(!is.na(DS4), DS4, NoDS4)) %>%
  rowwise() %>%
  mutate(QAM = paste(QAM1, QAM2, QAM3, QAM4)) %>% # paste meanings
  ungroup() %>%
  mutate(QAM = str_replace_all(QAM, "NA", ""),
         QAMeaning = if_else(!is.na(DSFull), DSFull, QAM)) # use full or concatenated meaning

# max characters
max(atlas_qa$Chars) # used to write code above

# check missing values
atlas_qa %>%
  filter((is.na(QAM4) & QACode4 != "")) %>%
  select(DataSource, QACode, QACode4, QAM4) %>%
  data.frame()
# none

atlas_qa %>%
  filter((is.na(QAM3) & QACode3 != "")) %>%
  select(DataSource, QACode, QACode3, QAM3) %>%
  data.frame()
# space

atlas_qa %>%
  filter((is.na(QAM2) & QACode2 != "")) %>%
  select(DataSource, QACode, QACode2, QAM2) %>%
  data.frame()
# commas and spaces
# 0.5 is assigned "unknown code" because of the 5

atlas_qa %>%
  filter((is.na(QAM1) & QACode1 != "")) %>%
  select(DataSource, QACode, QACode1, QAM1) %>%
  data.frame()

# simplify
atlas_qa2 <- atlas_qa %>%
  select(DataSource, QACode, QAMeaning) %>%
  arrange(QAMeaning, QACode, DataSource)

# save to manuall edit in Excel
write_csv(atlas_qa2, "intermediate-data/water_atlas_qa_codes_combined.csv")
# remove all error codes except (review with an expert):
# actual value is known to be greater/less than
# analysis from unpreservered/improperly preserved sample
# composite sample from multiple stations/mean of 2 or more determinations
# deviation from historical range
# detected value/good record/good quality
# diluted
# field measurement
# analyzed, but not detected
# held beyond accepted time
# secchi disk hits bottom
# significant rain in past 48 hours
# note that some combinations of these or these in combination with another warning were removed

# import with remove column
atlas_qa3 <- read_csv("intermediate-data/water_atlas_qa_codes_combined_manual.csv")

# add QA codes
atlas3 <- atlas2 %>%
  left_join(atlas_qa3) %>%
  mutate(Remove = replace_na(Remove, 0)) %>%
  filter(Remove == 0) %>%
  select(-Remove)

# before replacing NA above, check that all codes are accounted for
# filter(atlas3, !is.na(QACode) & QACode != "NULL" & is.na(Remove))
# should return 0


#### evaluate duplicates ####

# check for consistency in units
atlas3 %>%
  select(Parameter, Result_Unit) %>%
  unique() %>%
  arrange(Parameter)
# capitalization issues
# some are missing, but the units are in the name

# standardize parameters
atlas4 <- atlas3 %>%
  mutate(Parameter = tolower(Parameter))

# look at activity depth
atlas4 %>%
  select(Parameter, DepthUnits) %>%
  unique() %>%
  arrange(Parameter)
# all are in meters

atlas4 %>%
  ggplot(aes(x = ActivityDepth)) +
  geom_histogram() +
  facet_wrap(~ Parameter, scales = "free")
# no unusual values

# data sources per lake
atlas4 %>%
  group_by(WBodyID) %>%
  summarize(sources = n_distinct(DataSource)) %>%
  filter(sources > 1)
# several
  
# check for duplicates with same result value
atlas_dup <- get_dupes(atlas4, PermanentID, WBodyID, WaterBodyName, 
                       StationID, StationName, Actual_StationID,
                       Date, DataSource, Parameter, Result_Value)
# vary in activity depth, sample fraction, result_comment, but all have same value

# check for multiple result values for one sample
# included code for resolving multiples with different QA values, but decided to average them
atlas_dup2 <- get_dupes(atlas4, PermanentID, WBodyID, WaterBodyName, 
                        StationID, StationName, Actual_StationID, 
                        Date, DataSource, Parameter) %>%
  group_by(PermanentID, WBodyID, WaterBodyName, 
           StationID, StationName, Actual_StationID, 
           Date, DataSource, Parameter) %>%
  mutate(Result_vals = n_distinct(Result_Value), # multiple result values?
         # KeepCodes = case_when(sum(is.na(QACode)) > 0 ~ 1, # Positive codes
         #                       sum(QACode %in% c("D", "H", "A", "S", "E", "2", "x")) > 0 ~ 1,
         #                       TRUE ~ 0),
         QACodes = n_distinct(QACode) - as.numeric(sum(is.na(QACode)) > 0)) %>% # number of non-NA QA codes
  ungroup() %>%
  filter(Result_vals > 1) # %>%
  # mutate(KeepCode = case_when(is.na(QACode) ~ 1, # Positive codes
  #                             QACode %in% c("D", "H", "A", "S", "E", "2", "x") ~ 1,
  #                             TRUE ~ 0))

# QA codes - used to create list of acceptable codes above
# atlas_dup_qa <- atlas_dup2 %>%
#   select(QACode, QAMeaning) %>%
#   unique() %>% data.frame()

# multiple QA codes?
filter(atlas_dup2, QACodes > 1) %>% data.frame()
# resolve manually

# resolve multiple result values
# included code for resolving multiples with different QA values, but decided to average them
# atlas_dup3 <- atlas_dup2 %>%
#   filter(QACodes > 1) %>%
#   mutate(Remove = case_when(QACode == "U" ~ 1,
#                             QACode == "<" ~ 1,
#                             TRUE ~ 0)) %>%
#   full_join(atlas_dup2 %>%
#               filter(QACodes <= 1) %>%
#               mutate(Remove = case_when(QACodes == 1 & KeepCodes > 0 & KeepCode == 0 ~ 1)))

# how often are different data sources reporting the same value?
# how often are different data sources reporting different values?
atlas_source <- atlas4 %>%
  group_by(PermanentID, WBodyID, WaterBodyName, 
           StationID, StationName, Actual_StationID, 
           Date, Parameter) %>%
  summarize(TotSources = n_distinct(DataSource),
            TotResValues = n_distinct(Result_Value)) %>%
  ungroup() %>%
  filter(TotSources > 1) %>%
  left_join(atlas4 %>%
              group_by(PermanentID, WBodyID, WaterBodyName, 
                       StationID, StationName, Actual_StationID, 
                       Date, Parameter, Result_Value) %>%
              summarize(SameSources = n_distinct(DataSource)) %>%
              ungroup() %>%
              select(-Result_Value) %>%
              unique()) %>%
  mutate(SameValue = case_when(SameSources == 1 ~ 0,
                               TotSources == SameSources ~ 1))

# used to figure out SameValue
# atlas_source %>%
#   select(TotSources, TotResValues, SameSources) %>%
#   unique()

atlas_source %>%
  group_by(SameValue) %>%
  count() %>%
  ungroup() %>%
  mutate(tot = sum(n),
         perc = n/tot)
# 31% of duplicate sources have different values
# 69% of duplicate sources have same values


#### summarize by permanent ID and year ####

# check result comments
unique(atlas4$Result_Comment)
# most are uninterpretable
# others should have been accounted for with QA

# check for NA's
sum(is.na(atlas4$Result_Value))

# select relevant columns
# remove duplicate values (within and across data sources)
# summarize by permanent ID and month
atlas5 <- atlas4 %>%
  group_by(PermanentID, WBodyID, WaterBodyName, 
         StationID, StationName, Actual_StationID, 
         Year, Month, Date, Parameter, Result_Value) %>%  # remove duplicate result values (same value)
  summarize(QACode = paste(unique(QACode), collapse = "; "),
            QAMeaning = paste(unique(QAMeaning), collapse = "; ")) %>%
  ungroup() %>%
  group_by(PermanentID, WBodyID, WaterBodyName, 
           StationID, StationName, Actual_StationID, 
           Year, Month, Date, Parameter) %>%
  summarize(Result_Value = mean(Result_Value), # average multiple result values
            QACode = paste(unique(QACode), collapse = "; "),
            QAMeaning = paste(unique(QAMeaning), collapse = "; ")) %>%
  ungroup() %>%
  group_by(PermanentID, Year, Month, Date, Parameter) %>%
  summarize(Result_Value = mean(Result_Value), # average across stations
            QACode = paste(unique(QACode), collapse = "; "),
            QAMeaning = paste(unique(QAMeaning), collapse = "; "),
            WaterBodyName = paste(unique(WaterBodyName, collapse = "/")),
            StationsPerDate = n()) %>%
  ungroup() %>%
  group_by(PermanentID, WaterBodyName, Year, Month, Parameter) %>%
  summarize(Result_Value = mean(Result_Value), # average across dates within a month
            QACode = paste(unique(QACode), collapse = "; "),
            QAMeaning = paste(unique(QAMeaning), collapse = "; "),
            AvgStationsPerDate = mean(StationsPerDate),
            DatesSampled = n()) %>%
  ungroup()

# visulaize month distribution
ggplot(atlas5, aes(x = Month)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ Parameter)
# grouping samples by quarters includes at least one highly sampled month

# year quarters
quarts <- tibble(Month = 1:12,
                 Quarter = rep(c(4, 1:3), each = 3))

# summarize by permanent ID and quarter
atlas6 <- atlas5 %>%
  left_join(quarts) %>%
  group_by(PermanentID, WaterBodyName, Year, Quarter, Parameter) %>%
  summarize(Result_Value = mean(Result_Value), # average across months within a GS year
            QACode = paste(unique(QACode), collapse = "; "),
            QAMeaning = paste(unique(QAMeaning), collapse = "; "),
            AvgStationsPerDate = mean(AvgStationsPerDate),
            AvgDatesPerMonth = mean(DatesSampled),
            MonthsSampled = n_distinct(Month)) %>%
  ungroup() %>%
  mutate(Parameter = fct_recode(Parameter,
                                "TN_ug_L" = "tn_ugl",
                                "Secchi_ft" = "secchi_ft",
                                "CHL_ug_L" = "chla_ugl",
                                "TP_ug_L" = "tp_ugl"),
         QACode = if_else(QACode == "NA", NA_character_, QACode),
         QAMeaning = if_else(QAMeaning == "NA", NA_character_, QAMeaning)) %>%
  rename("QualityMetric" = "Parameter",
         "QualityValue" = "Result_Value")

# save data
write_csv(atlas6, "intermediate-data/water_atlas_quality_formatted.csv")
