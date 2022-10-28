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
  filter(QualityMetric %in% c("CHL_ug_L", "TP_ug_L", "TP_ug_L", "Secchi_ft", "Color_Pt_Co_Units", "Cond_uS_cm")) %>%
  ggplot(aes(x = Year, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_grid(QualityMetric ~ Quarter, scales = "free_y") +
  theme(legend.position = "none")
# some weirdly high P values
# probably won't stay in dataset because time series are incomplete

ggplot(wa_qual, aes(x = Year, y = QualityValue, color = PermanentID)) +
  geom_line() +
  facet_grid(QualityMetric ~ Quarter, scales = "free_y") +
  theme(legend.position = "none")
# some very high P values, but a lot are early in time series

# high values found in analysis
wa_qual %>%
  filter((PermanentID == "16798757" & Year == 2012) | (PermanentID == "78690681" & Year == 2016)) %>%
  filter(QualityMetric == "TP_ug_L")
# no QA codes
# Banana lake is just high that year (saw in N data)
# high P in time series for Garfield in 2016: https://polk.wateratlas.usf.edu/waterbodies/lakes/161045/lake-garfield

#### edit data ####

# add row for every year for each site/species combo (NA's for missing surveys)
lw_qual2 <- lw_qual %>%
  filter(QualityMetric == "TP_ug_L") %>%
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
  filter(QualityMetric == "TP_ug_L") %>%
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
  filter(QualityMetric == "TP_ug_L") %>%
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

write_csv(lw_qual2, "intermediate-data/LW_phosphorus_formatted.csv")
write_csv(wa_qual2, "intermediate-data/water_atlas_phosphorus_formatted.csv")
write_csv(qual2, "intermediate-data/LW_water_atlas_phosphorus_formatted.csv")


#### add invasive plant and management info ####

# source script for formatting invasive plant and water quality datasets
source("code/data-processing/invasive_plant_control_data_processing_for_quality_datasets.R")
source("code/generic-functions/continuous_time_interval.R")

# add invasive plant 
qual3 <- qual2 %>%
  mutate(GSYear = if_else(Quarter == 4, Year - 1, Year), # match quality in Jan-Apr with surveys from prev year (through April)
         TreatmentYear = Year - 1) %>% # match with treatment from previous year
  inner_join(inv_plant2) %>% # select waterbodies and years in both datasets
  inner_join(qual_ctrl2) %>%
  inner_join(perm_plant_ctrl) # select waterbodies that have had species and at least one year of management

# sample sizes
qual3 %>%
  group_by(CommonName, PermanentID) %>%
  summarize(Years = n_distinct(Year),
            Quarters = n_distinct(Quarter)) %>%
  ungroup() %>%
  ggplot(aes(x = Years)) +
  geom_histogram(binwidth = 1) +
  facet_grid(CommonName ~ Quarters, scales = "free")
# most have four quarters sampled
  
# complete time intervals
# surveys were not conducted every year on every lake for every species
# water quality was not measured every year on every lake for every metric
# control data is implicitly complete -- missing interpreted as no control
qual_time_int <- qual3 %>%
  distinct(GSYear, CommonName) %>% 
  mutate(out = pmap(., function(GSYear, CommonName) 
    time_int_qual_fun(year1 = GSYear, taxon = CommonName, dat_in = qual3))) %>%
  unnest(cols = out)

qual_time_int %>%
  distinct(CommonName, years_out, waterbodies) %>%
  ggplot(aes(x = years_out, y = waterbodies)) +
  geom_point() +
  facet_wrap(~ CommonName, scales = "free")
# low sample sizes

# select largest number of datapoints
qual_time_int2 <- qual_time_int %>%
  group_by(CommonName) %>%
  mutate(max_data_points = max(data_points)) %>%
  ungroup() %>%
  filter(data_points == max_data_points)

# filter dataset for lakes with complete surveys
qual4 <- qual3 %>%
  inner_join(qual_time_int2 %>%
               mutate(MinGSYear = GSYear,
                      MaxGSYear = GSYear + years_out - 1) %>% # 1 year added to years_out to count initial year
               select(CommonName, PermanentID, MinGSYear, MaxGSYear)) %>%
  filter(GSYear >= MinGSYear & GSYear <= MaxGSYear)

ggplot(qual4, aes(x = GSYear)) +
  geom_histogram() +
  facet_grid(Quarter ~ CommonName, scales = "free")

qual4 %>% 
  select(PermanentID, Year, CommonName) %>%
  unique() %>%
  ggplot(aes(x = Year)) +
  geom_histogram() +
  facet_grid(~ CommonName, scales = "free")

# select waterbodies sampled throughout
qual5 <- qual4 %>%
  group_by(CommonName) %>%
  mutate(maxYears = n_distinct(GSYear)) %>% # number of years per species
  ungroup() %>%
  group_by(CommonName, PermanentID, Quarter) %>%
  mutate(nYears = n_distinct(GSYear)) %>% # years per waterbody
  ungroup() %>%
  filter(nYears == maxYears) %>%
  mutate(ValueDiff = QualityValue - PrevValue,  # change over time
         across(ends_with("AvgPropCovered"), ~ .x * 100),
         logQual = log(QualityValue),
         Lag3APCsq = Lag3AvgPropCovered^2) %>% # square perc covered
  rename_with(str_replace, pattern = "AvgPropCovered", replacement = "AvgPercCovered")
