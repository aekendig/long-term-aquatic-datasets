#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)
library(lubridate)

# import data
lw_qual <- read_csv("intermediate-data/LW_quality_formatted.csv")
wa_qual <- read_csv("intermediate-data/water_atlas_quality_formatted.csv")
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
qual_ctrl <- read_csv("intermediate-data/FWC_quality_control_formatted.csv")


#### iniital visualizations ####
lw_qual %>%
  filter(QualityMetric %in% c("CHL_ug_L", "TN_ug_L", "TP_ug_L", "Secchi_ft", "Color_Pt_Co_Units", "Cond_uS_cm")) %>%
ggplot(aes(x = Year, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_grid(QualityMetric ~ Quarter, scales = "free_y") +
  theme(legend.position = "none")
# color and conductivity don't start until 2000-2006
# some weirdly high P values
# one very high chlorophyll value in quarter 2
# assume missing is missing, not zero values

ggplot(wa_qual, aes(x = Year, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_grid(QualityMetric ~ Quarter, scales = "free_y") +
  theme(legend.position = "none")
# a handful of high TN values, but not huge

# high value in final dataset (fwc_nitrogen_analysis)
wa_qual %>%
  filter(PermanentID == "16798757" & Year == 2012 & QualityMetric == "TN_ug_L")
# no QA codes, this extreme value shows up in time series
# https://www.polk.wateratlas.usf.edu/waterbodies/lakes/160778/

wa_qual %>%
  filter(PermanentID == "c3aa4be1-a57b-47a6-8389-edde0c50422f" & Year == 2016 & QualityMetric == "TN_ug_L") %>%
  data.frame()


#### edit data ####

# add row for every year for each site/species combo (NA's for missing surveys)
lw_qual2 <- lw_qual %>%
  filter(QualityMetric == "TN_ug_L") %>%
  full_join(lw_qual %>%
              select(PermanentID) %>%
              unique() %>%
              expand_grid(Year = min(lw_qual$Year):max(lw_qual$Year)) %>%
              expand_grid(Quarter = 1:4)) %>%
  group_by(PermanentID, Quarter) %>%
  arrange(Year) %>% 
  mutate(PrevValue = lag(QualityValue)) %>% # previous year's value
  ungroup() %>%
  filter(!is.na(QualityValue))

# repeat for Water Atlas
wa_qual2 <- wa_qual %>%
  filter(QualityMetric == "TN_ug_L") %>%
  full_join(lw_qual %>%
              select(PermanentID) %>%
              unique() %>%
              expand_grid(Year = min(lw_qual$Year):max(lw_qual$Year)) %>%
              expand_grid(Quarter = 1:4)) %>%
  group_by(PermanentID, Quarter) %>%
  arrange(Year) %>% 
  mutate(PrevValue = lag(QualityValue)) %>% # previous year's value
  ungroup() %>%
  filter(!is.na(QualityValue))

# combine data
qual <- lw_qual %>%
  full_join(wa_qual) %>%
  filter(QualityMetric == "TN_ug_L") %>%
  group_by(PermanentID, Year, Quarter) %>%
  summarize(QualityValue = mean(QualityValue),
            MonthsSampled = mean(MonthsSampled)) %>%
  ungroup()

qual2 <-  qual %>%
  full_join(qual %>%
              select(PermanentID) %>%
              unique() %>%
              expand_grid(Year = min(qual$Year):max(qual$Year)) %>%
              expand_grid(Quarter = 1:4)) %>%
  group_by(PermanentID, Quarter) %>%
  arrange(Year) %>% 
  mutate(PrevValue = lag(QualityValue)) %>% # previous year's value
  ungroup() %>%
  filter(!is.na(QualityValue))


#### save water quality datasets #####

write_csv(lw_qual2, "intermediate-data/LW_nitrogen_formatted.csv")
write_csv(wa_qual2, "intermediate-data/water_atlas_nitrogen_formatted.csv")
write_csv(qual2, "intermediate-data/LW_water_atlas_nitrogen_formatted.csv")

#### add invasive plant and management info ####

lwwa_nit
