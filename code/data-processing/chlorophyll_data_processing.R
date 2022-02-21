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

ggplot(lw_qual, aes(x = GSYear, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ QualityMetric, scales = "free_y") +
  theme(legend.position = "none")
# color and conductivity don't start until 2000-2006
# some weirdly high P values
# assume missing is missing, not zero values

ggplot(wa_qual, aes(x = GSYear, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ QualityMetric, scales = "free_y") +
  theme(legend.position = "none")
# one super high chlorophyll value (don't exceed 400 in LW data)
# some very high P values

wa_qual %>%
  filter(QualityMetric == "CHL_ug_L" & QualityValue > 2e4) %>%
  data.frame()
# Lawne Lake in 1971 -- won't be going back that far
# Noreast Lake in 2006
# neither has a QACode

wa_qual %>%
  filter(GSYear >= 2000 & !(QualityMetric == "CHL_ug_L" & QualityValue > 2e4)) %>%
  ggplot(aes(x = GSYear, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ QualityMetric, scales = "free_y") +
  theme(legend.position = "none")


#### edit data ####

# add row for every year for each site/species combo (NA's for missing surveys)
lw_chl <- lw_qual %>%
  filter(QualityMetric == "CHL_ug_L") %>%
  full_join(lw_qual %>%
              select(PermanentID) %>%
              unique() %>%
              expand_grid(GSYear = min(lw_qual$GSYear):max(lw_qual$GSYear))) %>%
  group_by(PermanentID) %>%
  arrange(GSYear) %>% 
  mutate(PrevValue = lag(QualityValue)) %>% # previous year's value
  ungroup() %>%
  mutate(RatioQual = QualityValue / PrevValue,
         LogRatioQual = log(RatioQual)) %>%
  filter(!is.na(QualityValue) & !is.na(PrevValue) & !(LogRatioQual %in% c(Inf, -Inf)))
# one lake has a chlorophyll value of 0 in 2002, resulting in -Inf

# check values
ggplot(lw_chl, aes(x = GSYear, y = LogRatioQual, color = PermanentID)) +
  geom_line() +
  theme(legend.position = "none")

filter(lw_chl, LogRatioQual < -3) %>% data.frame()

# repeat for Water Atlas
wa_chl <- wa_qual %>%
  filter(QualityMetric == "CHL_ug_L") %>%
  full_join(lw_qual %>%
              select(PermanentID) %>%
              unique() %>%
              expand_grid(GSYear = min(lw_qual$GSYear):max(lw_qual$GSYear))) %>%
  group_by(PermanentID) %>%
  arrange(GSYear) %>% 
  mutate(PrevValue = lag(QualityValue)) %>% # previous year's value
  ungroup() %>%
  mutate(RatioQual = QualityValue / PrevValue,
         LogRatioQual = log(RatioQual)) %>%
  filter(!is.na(QualityValue) & !is.na(PrevValue) & !(LogRatioQual %in% c(Inf, -Inf)))
# one lake has a zero value in one year

# check values
ggplot(wa_chl, aes(x = GSYear, y = LogRatioQual, color = PermanentID)) +
  geom_line() +
  theme(legend.position = "none")
# much larger values than LW

# combine data
qual <- lw_qual %>%
  full_join(wa_qual) %>%
  filter(QualityMetric == "CHL_ug_L") %>%
  group_by(PermanentID, GSYear) %>%
  summarize(QualityValue = mean(QualityValue),
            MonthsSampled = mean(MonthsSampled)) %>%
  ungroup()

lwwa_chl <-  qual %>%
  full_join(qual %>%
              select(PermanentID) %>%
              unique() %>%
              expand_grid(GSYear = min(qual$GSYear):max(qual$GSYear))) %>%
  group_by(PermanentID) %>%
  arrange(GSYear) %>% 
  mutate(PrevValue = lag(QualityValue)) %>% # previous year's value
  ungroup() %>%
  mutate(RatioQual = QualityValue / PrevValue,
         LogRatioQual = log(RatioQual)) %>%
  filter(!is.na(QualityValue) & !is.na(PrevValue) & !(LogRatioQual %in% c(Inf, -Inf)))

# check values
ggplot(lwwa_chl, aes(x = GSYear, y = LogRatioQual, color = PermanentID)) +
  geom_line() +
  theme(legend.position = "none")


#### output #####
write_csv(lw_chl, "intermediate-data/LW_chlorophyll_formatted.csv")
write_csv(wa_chl, "intermediate-data/water_atlas_chlorophyll_formatted.csv")
write_csv(lwwa_chl, "intermediate-data/LW_water_atlas_chlorophyll_formatted.csv")
