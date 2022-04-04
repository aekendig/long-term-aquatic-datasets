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


#### iniital visualizations ####
lw_qual %>%
  filter(QualityMetric %in% c("CHL_ug_L", "TN_ug_L", "TP_ug_L", "Secchi_ft", "Color_Pt_Co_Units", "Cond_uS_cm")) %>%
  ggplot(aes(x = GSYear, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_grid(QualityMetric ~ Quarter, scales = "free_y") +
  theme(legend.position = "none")
# color and conductivity don't start until 2000-2006
# some weirdly high P values
# one very high chlorophyll value in quarter 2
# assume missing is missing, not zero values

ggplot(wa_qual, aes(x = GSYear, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_grid(QualityMetric ~ Quarter, scales = "free_y") +
  theme(legend.position = "none")
# no extreme secchi values


#### edit data ####

# add row for every year for each site/species combo (NA's for missing surveys)
lw_sec <- lw_qual %>%
  filter(QualityMetric == "Secchi_ft") %>%
  full_join(lw_qual %>%
              select(PermanentID) %>%
              unique() %>%
              expand_grid(GSYear = min(lw_qual$GSYear):max(lw_qual$GSYear)) %>%
              expand_grid(Quarter = 1:4)) %>%
  group_by(PermanentID, Quarter) %>%
  arrange(GSYear) %>% 
  mutate(PrevValue = lag(QualityValue)) %>% # previous year's value
  ungroup() %>%
  filter(!is.na(QualityValue))

# repeat for Water Atlas
wa_sec <- wa_qual %>%
  filter(QualityMetric == "Secchi_ft") %>%
  full_join(lw_qual %>%
              select(PermanentID) %>%
              unique() %>%
              expand_grid(GSYear = min(lw_qual$GSYear):max(lw_qual$GSYear)) %>%
              expand_grid(Quarter = 1:4)) %>%
  group_by(PermanentID, Quarter) %>%
  arrange(GSYear) %>% 
  mutate(PrevValue = lag(QualityValue)) %>% # previous year's value
  ungroup() %>%
  filter(!is.na(QualityValue))

# combine data
qual <- lw_qual %>%
  full_join(wa_qual) %>%
  filter(QualityMetric == "Secchi_ft") %>%
  group_by(PermanentID, GSYear, Quarter) %>%
  summarize(QualityValue = mean(QualityValue),
            MonthsSampled = mean(MonthsSampled)) %>%
  ungroup()

lwwa_sec <-  qual %>%
  full_join(qual %>%
              select(PermanentID) %>%
              unique() %>%
              expand_grid(GSYear = min(qual$GSYear):max(qual$GSYear)) %>%
              expand_grid(Quarter = 1:4)) %>%
  group_by(PermanentID, Quarter) %>%
  arrange(GSYear) %>% 
  mutate(PrevValue = lag(QualityValue)) %>% # previous year's value
  ungroup() %>%
  filter(!is.na(QualityValue))


#### output #####
write_csv(lw_sec, "intermediate-data/LW_secchi_formatted.csv")
write_csv(wa_sec, "intermediate-data/water_atlas_secchi_formatted.csv")
write_csv(lwwa_sec, "intermediate-data/LW_water_atlas_secchi_formatted.csv")
