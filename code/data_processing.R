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
         AreaOfInterest = str_replace_all(AreaOfInterest, c(", Lake" = "", "Lake " = "", " Lake" = "")) %>%
           toupper())

# county overlap
gnis2 %>%
  filter(County_FWC != County_LW)
# Apopka: in both counties
# Isabell: in Polk, but on the edge of Highlands
# Monroe: in both counties

# lake overlap
gnis2 %>%
  filter(Lake_LW != AreaOfInterest) %>%
  arrange(County_LW, Lake_LW) %>%
  data.frame()
# multiple GNIS values assigned to a single FWC lake
  # Little Red Water in Highlands
  # Blue in Polk County
  # Rodman Reservoir in Putnam County
# others are different spellings

#### start here ####


#### control data ####

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

# make county names all uppercase
# add GNIS ID's to control data
# matching variable: AreaOfInterest
# for duplicate AreaOfInterest, check County
ctrl_loc2 <- ctrl_loc %>%
  mutate(County = toupper(County)) %>%
  left_join(gnis %>%
              mutate(County = toupper(County_FWC)))

# check if duplicate AreaOfInterests were matched correctly
ctrl_loc2 %>%
  filter(!is.na(GNIS_ID) & AreaOfInterest %in% ctrl_AOI$AreaOfInterest)

# check missing GNIS
ctrl_loc2 %>%
  filter(is.na(GNIS_ID)) %>%
  nrow()
# 220 entries

# duplicated GNIS
ctrl_loc3 <- ctrl_loc2 %>%
  mutate(dup_GNIS_ID = duplicated(GNIS_ID),
         dup_GNIS_ID = ifelse(is.na(GNIS_ID), NA, dup_GNIS_ID))

dup_GID <- ctrl_loc3 %>%
  filter(dup_GNIS_ID == T) %>%
  select(GNIS_ID) # 8 duplicated

ctrl_loc4 <- ctrl_loc3 %>%
  mutate(dup_GNIS_ID = ifelse(GNIS_ID %in% dup_GID$GNIS_ID, T, dup_GNIS_ID))

sum(ctrl_loc4$dup_GNIS_ID, na.rm = T) # 15

# save problem data for consultation
write_csv(ctrl_loc4 %>%
            filter(dup_GNIS_ID == T |
                     County != County_FWC),
          "intermediate-data/FWC_herbicides_GNIS_issue_20201118.csv")

write_csv(ctrl_loc4 %>%
            filter(is.na(GNIS_ID)),
          "intermediate-data/FWC_herbicides_GNIS_missing_20201118.csv")
