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
  filter(QualityMetric %in% c("CHL_ug_L", "CHL_ug_L", "TP_ug_L", "Secchi_ft", "Color_Pt_Co_Units", "Cond_uS_cm")) %>%
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
# some super high chlorophyll values (don't exceed 400 in LW data)
# some very high P values

wa_qual %>%
  filter(QualityMetric == "CHL_ug_L" & QualityValue > 1e4) %>%
  data.frame()
# most are pre-1980 -- won't be going back that far
# Noreast Lake in 2006
# neither has a QACode

wa_qual %>%
  filter(Year >= 2000 & !(QualityMetric == "CHL_ug_L" & QualityValue > 1e4)) %>%
  ggplot(aes(x = Year, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_grid(QualityMetric ~ Quarter, scales = "free_y") +
  theme(legend.position = "none")


#### edit data ####

# add row for every year for each site/species combo (NA's for missing surveys)
lw_qual2 <- lw_qual %>%
  filter(QualityMetric == "CHL_ug_L") %>%
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
  filter(QualityMetric == "CHL_ug_L") %>%
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
  filter(QualityMetric == "CHL_ug_L") %>%
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

write_csv(lw_qual2, "intermediate-data/LW_chlorophyll_formatted.csv")
write_csv(wa_qual2, "intermediate-data/water_atlas_chlorophyll_formatted.csv")
write_csv(qual2, "intermediate-data/LW_water_atlas_chlorophyll_formatted.csv")

#### add invasive plant and management info ####