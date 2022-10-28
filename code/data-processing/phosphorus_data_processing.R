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


#### select complete times and waterbodies ####

# sample sizes
qual3 %>%
  group_by(CommonName, PermanentID) %>%
  summarize(Years = n_distinct(Year),
            Quarters = n_distinct(Quarter)) %>%
  ungroup() %>%
  ggplot(aes(x = Years)) +
  geom_histogram(binwidth = 1) +
  facet_grid(CommonName ~ Quarters, scales = "free")
# most have four quarters sampled, but these may not be continuous
  
# complete time intervals
# surveys were not conducted every year on every lake for every species
# water quality was not measured every year on every lake for every metric
# control data is implicitly complete -- missing interpreted as no control
qual_time_int <- qual3 %>%
  distinct(GSYear, CommonName, Quarter) %>% 
  mutate(out = pmap(., function(GSYear, CommonName, Quarter) 
    time_int_qual_fun(year1 = GSYear, taxon = CommonName, quarter = Quarter,
                      dat_in = qual3))) %>%
  unnest(cols = out)

qual_time_int %>%
  distinct(CommonName, Quarter, years_out, waterbodies) %>%
  ggplot(aes(x = years_out, y = waterbodies)) +
  geom_point() +
  facet_grid(CommonName ~ Quarter, scales = "free")

# select largest number of datapoints
qual_time_int2 <- qual_time_int %>%
  group_by(CommonName, Quarter) %>%
  mutate(max_data_points = max(data_points)) %>%
  ungroup() %>%
  filter(data_points == max_data_points) %>%
  group_by(CommonName, Quarter) %>%
  mutate(max_waterbodies = max(waterbodies)) %>%
  ungroup() %>%
  filter(waterbodies == max_waterbodies)

# summary
qual_time_int2 %>%
  distinct(CommonName, Quarter, GSYear, years_out, waterbodies, data_points) %>%
  arrange(CommonName, Quarter) %>%
  data.frame()
# most species have the most data points in the first quarter

# filter dataset for lakes with complete surveys
qual_time_int3 <- qual3 %>%
  inner_join(qual_time_int2 %>%
               mutate(MinGSYear = GSYear,
                      MaxGSYear = GSYear + years_out - 1) %>% # 1 year added to years_out to count initial year
               select(CommonName, Quarter, PermanentID, MinGSYear, MaxGSYear)) %>%
  filter(GSYear >= MinGSYear & GSYear <= MaxGSYear)

# check for irregularities
ggplot(qual_time_int3, aes(x = GSYear)) +
  geom_histogram(binwidth = 1) +
  facet_grid(Quarter ~ CommonName, scales = "free")

# modify columns
qual_time_int4 <- qual_time_int3 %>%
  mutate(ValueDiff = QualityValue - PrevValue,  # change over time
         logQual = log(QualityValue)) # square perc covered

# save
write_csv(qual_time_int4, "intermediate-data/phosphorus_inv_ctrl_formatted.csv")


#### select complete dataset for recent treatment ####

# complete time intervals
qual_rec_int <- qual3 %>%
  distinct(GSYear, CommonName, Quarter) %>% 
  mutate(out = pmap(., function(GSYear, CommonName, Quarter) 
    time_int_qual_recent_fun(year1 = GSYear, taxon = CommonName, quarter = Quarter,
                      dat_in = qual3))) %>%
  unnest(cols = out)

qual_rec_int %>%
  distinct(CommonName, Quarter, years_out, waterbodies) %>%
  ggplot(aes(x = years_out, y = waterbodies)) +
  geom_point() +
  facet_grid(CommonName ~ Quarter, scales = "free")

# select largest number of datapoints
qual_rec_int2 <- qual_rec_int %>%
  group_by(CommonName, Quarter) %>%
  mutate(max_data_points = max(data_points)) %>%
  ungroup() %>%
  filter(data_points == max_data_points) %>%
  group_by(CommonName, Quarter) %>%
  mutate(max_waterbodies = max(waterbodies)) %>%
  ungroup() %>%
  filter(waterbodies == max_waterbodies)

# summary
qual_rec_int2 %>%
  distinct(CommonName, Quarter, GSYear, years_out, waterbodies, data_points) %>%
  arrange(CommonName, Quarter) %>%
  data.frame()
# most species have the most data points in the first quarter

# filter dataset for lakes with complete surveys
qual_rec_int3 <- qual3 %>%
  inner_join(qual_rec_int2 %>%
               mutate(MinGSYear = GSYear,
                      MaxGSYear = GSYear + years_out - 1) %>% # 1 year added to years_out to count initial year
               select(CommonName, Quarter, PermanentID, MinGSYear, MaxGSYear)) %>%
  filter(GSYear >= MinGSYear & GSYear <= MaxGSYear)

# check for irregularities
ggplot(qual_rec_int3, aes(x = GSYear)) +
  geom_histogram(binwidth = 1) +
  facet_grid(Quarter ~ CommonName, scales = "free")

# modify columns
qual_rec_int4 <- qual_rec_int3 %>%
  mutate(ValueDiff = QualityValue - PrevValue,  # change over time
         logQual = log(QualityValue)) # square perc covered

# save
write_csv(qual_rec_int4, "intermediate-data/phosphorus_inv_ctrl_recent_formatted.csv")


#### uninvaded waterbodies ####

# import data
uninv <- read_csv("output/fwc_uninvaded_permID.csv") # lakes with no recorded invasion

# floating uninv
uninv_float <- uninv %>%
  filter(CommonName %in% c("Water hyacinth", "Water lettuce")) %>%
  group_by(PermanentID, Treatments) %>%
  summarize(nUninv = n()) %>%
  ungroup() %>%
  filter(nUninv == 2) %>%
  select(-nUninv) %>%
  mutate(CommonName = "floating plants",
         Established = 0)

# add water quality to uninvaded dataset
# select years to match invasion dataset
uninv2 <- qual2 %>%
  mutate(GSYear = if_else(Quarter == 4, Year - 1, Year)) %>% # match quality in Jan-Apr with surveys from prev year (through April)
  inner_join(uninv %>%
               full_join(uninv_float) %>%
               filter(!(CommonName %in% c("Water hyacinth", "Water lettuce")))) %>%
  left_join(qual5 %>%
              distinct(CommonName, Quarter, MinGSYear, MaxGSYear)) %>%
  filter(GSYear >= MinGSYear & GSYear <= MaxGSYear) %>%
  mutate(logQual = log(QualityValue))

# save
write_csv(uninv2, "intermediate-data/phosphorus_fwc_uninvaded_permID.csv")