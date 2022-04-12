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
# some weirdly high P values
# probably won't stay in dataset because time series are incomplete

ggplot(wa_qual, aes(x = GSYear, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_grid(QualityMetric ~ Quarter, scales = "free_y") +
  theme(legend.position = "none")
# some very high P values, but a lot are early in time series

# high values found in analysis
wa_qual %>%
  filter((PermanentID == "16798757" & GSYear == 2012) | (PermanentID == "78690681" & GSYear == 2016)) %>%
  filter(QualityMetric == "TP_ug_L")
# no QA codes
# Banana lake is just high that year (saw in N data)
# high P in time series for Garfield in 2016: https://polk.wateratlas.usf.edu/waterbodies/lakes/161045/lake-garfield

#### edit data ####

# add row for every year for each site/species combo (NA's for missing surveys)
lw_pho <- lw_qual %>%
  filter(QualityMetric == "TP_ug_L") %>%
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
wa_pho <- wa_qual %>%
  filter(QualityMetric == "TP_ug_L") %>%
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
  filter(QualityMetric == "TP_ug_L") %>%
  group_by(PermanentID, GSYear, Quarter) %>%
  summarize(QualityValue = mean(QualityValue),
            MonthsSampled = mean(MonthsSampled)) %>%
  ungroup()

lwwa_pho <-  qual %>%
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
write_csv(lw_pho, "intermediate-data/LW_phosphorus_formatted.csv")
write_csv(wa_pho, "intermediate-data/water_atlas_phosphorus_formatted.csv")
write_csv(lwwa_pho, "intermediate-data/LW_water_atlas_phosphorus_formatted.csv")
