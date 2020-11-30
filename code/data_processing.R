#### info ####

# goal: combine datasets


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
ctrl <- read_csv("original-data/FWC_Herbicide_Treatments_mid2010-Oct2020.csv")
gnis <- read_csv("original-data/LW_matching_Herbicide_lakes_with_GNIS.csv")
qual <- read_csv("original-data/Lakewatch_1986-2019_TP_TN_Chl_Secchi_Color_Cond.csv")


#### control data ####

# unique AreaOfInterest
unique(ctrl$AreaOfInterest) %>%
  length() # 446

# unique AreaOfInterest ID
unique(ctrl$AreaOfInterestID) %>%
  length() # 454

# extract location info
ctrl_loc <- ctrl %>%
  select(AreaOfInterest,
         AreaOfInterestID,
         County,
         WMDNAME,
         longitude,
         latitude) %>%
  unique() # 454

# check for unique ID
duplicated(ctrl_loc$AreaOfInterestID) %>%
  sum() # all unique

duplicated(ctrl_loc$AreaOfInterest) %>%
  sum() # 8 duplicates

# add GNIS ID's to control data
# matching variable: AreaOfInterest
ctrl_loc2 <- ctrl_loc %>%
  left_join(gnis)

# check county names and remove if okay
ctrl_loc2 %>%
  filter(County != County_FWC) %>%
  select(County, County_FWC, AreaOfInterest, Lake_LW) %>%
  unique()
# 14 entries, 6 lake names

# check missing GNIS
ctrl_loc2 %>%
  filter(is.na(GNIS_ID)) %>%
  nrow()
# 218 entries

# duplicated GNIS
ctrl_loc3 <- ctrl_loc2 %>%
  mutate(dup_GNIS_ID = duplicated(GNIS_ID),
         dup_GNIS_ID = ifelse(is.na(GNIS_ID), NA, dup_GNIS_ID))

dup_GID <- ctrl_loc3 %>%
  filter(dup_GNIS_ID == T) %>%
  select(GNIS_ID) # 27 duplicated

ctrl_loc4 <- ctrl_loc3 %>%
  mutate(dup_GNIS_ID = ifelse(GNIS_ID %in% dup_GID$GNIS_ID, T, dup_GNIS_ID))

sum(ctrl_loc4$dup_GNIS_ID, na.rm = T)

# save problem data for consultation
write_csv(ctrl_loc4 %>%
            filter(dup_GNIS_ID == T |
                     County != County_FWC),
          "intermediate-data/FWC_herbicides_GNIS_issue_20201118.csv")

write_csv(ctrl_loc4 %>%
            filter(is.na(GNIS_ID)),
          "intermediate-data/FWC_herbicides_GNIS_missing_20201118.csv")
