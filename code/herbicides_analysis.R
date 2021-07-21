
#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)
library(glmmTMB)
library(effects)

# import data
ctrl_old <- read_csv("intermediate-data/FWC_control_old_formatted.csv")
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
dayDat <- read_csv("intermediate-data/LakeO_day_data_for_model.csv")

# load models
load("output/lakeO_intraannual_veg_change_mod.rda")

# figure settings
def_theme <- theme_bw() +
  theme(axis.text = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.text = element_text(size = 14, color="black"),
        strip.background = element_blank())

# functions to transform data to account for 0's and 1's
transform01 <- function(x) {
  n <- sum(!is.na(x))
  (x * (n - 1) + 0.5) / n
}

backtransform01 <- function(x) {
  n <- sum(!is.na(x))
  (x * n - 0.5) / (n - 1)
}  

# function to remove duplicates from invasive plant data
rem_dups_fun <- function(dat) {
  
  if (nrow(dat) == 1){
    dat2 <- dat
  } else if (var(dat$AreaCovered_ha) > 0){
    dat2 <- dat %>%
      filter(AreaCovered_ha == max(AreaCovered_ha)) # choose survey with maximum area
  } else if (var(dat$AreaCovered_ha) == 0){
    dat2 <- dat %>%
      mutate(DaysDiff = abs(lakeO_days - Days)) %>%
      filter(DaysDiff == min(DaysDiff)) %>% # choose survey closes to max. abundance time
      select(-DaysDiff)
  }
  
  return(dat2)
}

# function for cumulative herbicide
ctrl_lag_fun <- function(GSYear, Lag, DSet){
  
  # parameters
  ModYear <- GSYear
  
  # filter dataset
  if(DSet == "inv"){
    subdat <- ctrl_inv %>%
      filter(GSYear <= ModYear & GSYear >= (ModYear - Lag))
  }else if(DSet == "nat_init_imm"){
    subdat <- ctrl_nat_init_imm %>%
      filter(GSYear <= ModYear & GSYear >= (ModYear - Lag))
  }else if(DSet == "nat_ext"){
    subdat <- ctrl_nat_ext %>%
      filter(GSYear <= ModYear & GSYear >= (ModYear - Lag))
  }else if(DSet == "nat_rept_imm"){
    subdat <- ctrl_nat_rept_imm %>%
      filter(GSYear <= ModYear & GSYear >= (ModYear - Lag))
  }
  
  # summarize
  outdat <- subdat %>%
    group_by(PermanentID, SpeciesName, Species, GSYear) %>% # annual summary
    summarise(PropTreated = sum(PropTreated)) %>% # add proportion lake treated for all treatments in a year (can exceed 1)
    ungroup() %>%
    group_by(PermanentID, SpeciesName, Species) %>% # summarize over lag years
    summarise(PropTreated = mean(PropTreated), # average proportion treated per year 
              Treated = as.numeric(PropTreated > 0)) %>%
    ungroup() %>%
    mutate(GSYear = ModYear,
           Lag = Lag)
  
  # return
  return(outdat)
}

# function for cumulative invasive abundance
inv_lag_fun <- function(GSYear, Lag){
  
  # parameters
  ModYear <- GSYear
  
  # filter dataset
  subdat <- inv_fwc2 %>% # use version with NA values
    filter(GSYear <= ModYear & GSYear >= (ModYear - Lag))
  
  # summarize
  outdat <- subdat %>%
    group_by(PermanentID, SpeciesName) %>% # summarize over lag
    summarise(PropCovered = mean(PropCovered)) %>%
    ungroup() %>%
    pivot_wider(names_from = SpeciesName,
                values_from = PropCovered) %>%
    rename(Hydrilla = "Hydrilla verticillata",
           WaterLettuce = "Pistia stratiotes",
           WaterHyacinth = "Eichhornia crassipes") %>%
    mutate(GSYear = ModYear,
           Lag = Lag)
  
  # return
  return(outdat)
}

# # function to add invasive plants and herbicides to native plant data
# inv_ctrl_fun <- function(FirstGS, GSYear, PermanentID){
#   
#   # parameters
#   window_start <- FirstGS
#   window_end <- GSYear
#   perm_ID <- PermanentID
#   
#   # average invasive plant cover
#   inv_dat <- inv_fwc2 %>%
#     filter(GSYear <= window_end & GSYear >= window_start & PermanentID == perm_ID) %>%
#     mutate(GSYear = window_end) %>%
#     group_by(GSYear, PermanentID, CommonName) %>% # one value for each invasive species
#     summarise(AvgPropCovered = mean(PropCovered)) %>%
#     ungroup() %>%
#     pivot_wider(names_from = CommonName,
#                 values_from = AvgPropCovered) %>%
#     rename(WaterLettuce = "Water lettuce",
#            WaterHyacinth = "Water hyacinth")
#   
#   # output
#   return(inv_dat)
# }


#### lake O parameters ####

# extract coefficients
lakeO_beta1 <- coef(summary(lake0_mod))[2, "Estimate"] # days
lakeO_beta2 <- coef(summary(lake0_mod))[3, "Estimate"] # days^2

# date of max abundance (-b/2a)
lakeO_days <- -lakeO_beta1 / (2 * lakeO_beta2)


#### surveyor experience ####

# surveys per surveyor
surveyor <- plant_fwc %>%
  select(Surveyor, AreaOfInterest, AreaOfInterestID, PermanentID, SurveyDate) %>%
  unique() %>%
  mutate(Survey = 1) %>%
  arrange(SurveyDate, AreaOfInterestID) %>%
  group_by(Surveyor) %>%
  mutate(SurveyorExperience = cumsum(Survey)) %>%
  ungroup() %>%
  select(-Survey)

# surveys without surveyor
filter(surveyor, is.na(Surveyor))

# remove NAs
surveyor2 <- surveyor %>%
  filter(!is.na(Surveyor)) %>%
  mutate(SurveyorExperienceCS = (SurveyorExperience - mean(SurveyorExperience)) / sd(SurveyorExperience))


#### edit invasive plant data ####

# Perm IDs per AOI
plant_fwc %>%
  group_by(AreaOfInterestID) %>%
  summarise(IDs = length(unique(PermanentID))) %>%
  filter(IDs > 1)
# one Perm ID per AOI

# min area
plant_fwc %>%
  filter(SpeciesAcres > 0) %>%
  summarise(min = min(SpeciesAcres)) # 0.01

# invasive plant dataset
inv_fwc <- plant_fwc %>% # start with all surveys
  select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, Surveyor) %>%
  unique() %>% # one row per survey
  expand_grid(tibble(SpeciesName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"))) %>% # one row per species per survey
  full_join(plant_fwc %>% # add invasive plant information
              filter(SpeciesName %in% c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes")) %>%
              select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, Surveyor, SpeciesName, SpeciesAcres)) %>%
  mutate(SpeciesAcres = replace_na(SpeciesAcres, 0), # cover 0 when it wasn't in a survey
         Area_ha = ShapeArea * 100, # convert lake area from km-squared to hectares
         AreaCovered_ha = SpeciesAcres * 0.405, # convert plant cover from acres to hectares
         SurveyMonth = month(SurveyDate),
         SurveyDay = day(SurveyDate),
         SurveyYear = year(SurveyDate),
         GSYear = case_when(SurveyMonth >= 4 ~ SurveyYear,
                            SurveyMonth < 4 ~ SurveyYear - 1), # assume growing season starts in April
         MonthDay = case_when(SurveyMonth >= 4 ~ as.Date(paste("2020", SurveyMonth, SurveyDay, sep = "-")), # start "year" in April (2020/2021 are arbitrary)
                              SurveyMonth < 4 ~ as.Date(paste("2021", SurveyMonth, SurveyDay, sep = "-"))), # this is for joining dayDat
         CommonName = case_when(SpeciesName == "Eichhornia crassipes" ~ "Water hyacinth", 
                                SpeciesName == "Hydrilla verticillata" ~ "Hydrilla", 
                                SpeciesName == "Pistia stratiotes" ~ "Water lettuce")) %>%
  left_join(dayDat) %>% # add standardized days (proportion between April 1 and March 31)
  mutate(AreaChangeSD = lakeO_beta1 * (lakeO_days-Days) + lakeO_beta2 * (lakeO_days^2 - Days^2)) %>% # calculate the number of sd's to change to get est. max abundance
  group_by(AreaOfInterestID, SpeciesName) %>% # take standard deviation by survey area and species
  mutate(EstAreaCoveredRaw_ha = AreaCovered_ha + AreaChangeSD * sd(AreaCovered_ha)) %>% # calculate est. max abundance, NA if only one value is available
  ungroup() %>%
  left_join(surveyor2) # surveyor experience

# remove duplicates
# summarize by waterbody
inv_fwc2 <- inv_fwc %>%
  nest(data = c(SurveyDate, Surveyor, SurveyorExperience, SurveyorExperienceCS, SpeciesAcres, AreaCovered_ha, SurveyMonth, SurveyDay, SurveyYear, MonthDay, Days, AreaChangeSD, EstAreaCoveredRaw_ha)) %>% # find multiple surveys within area of interest, growing season year, and species
  mutate(newdata = map(data, ~rem_dups_fun(.))) %>% # remove duplicates
  select(-data) %>% # removes 123 rows of data
  unnest(newdata) %>%
  group_by(PermanentID, Area_ha, GSYear, SpeciesName, CommonName) %>% # summarize for multiple AOIs in one PermanentID (i.e., waterbody)
  summarise(AreaName = paste(AreaOfInterest, collapse = "/"),
            SurveyDate = max(SurveyDate),
            SurveyorExperience = mean(SurveyorExperience, na.rm = T),
            SurveyorExperienceCS = mean(SurveyorExperienceCS, na.rm = T),
            SpeciesAcres = sum(SpeciesAcres),
            AreaCovered_ha = sum(AreaCovered_ha),
            EstAreaCoveredRaw_ha = sum(EstAreaCoveredRaw_ha)) %>%
  ungroup() %>%
  mutate(EstAreaCovered_ha = case_when(EstAreaCoveredRaw_ha > Area_ha ~ Area_ha, # reduce areas covered to total area
                                       TRUE ~ EstAreaCoveredRaw_ha),
         PropCovered = EstAreaCovered_ha / Area_ha,
         PropCoveredAdj = case_when(PropCovered < 1e-3 ~ 1e-3, # avoid super small values that skew ratios
                           TRUE ~ PropCovered),
         SpeciesPresent = case_when(SpeciesAcres > 0 ~ 1,
                                    SpeciesAcres == 0 ~ 0)) %>%
  full_join(inv_fwc %>% # add row for every year for each site/species combo (NA's for missing surveys)
              select(PermanentID, SpeciesName) %>%
              unique() %>%
              expand_grid(GSYear = min(inv_fwc$GSYear):max(inv_fwc$GSYear))) %>%
  group_by(PermanentID, SpeciesName) %>%
  arrange(GSYear) %>% 
  mutate(PrevPropCovered = lag(PropCovered)) %>% # previous year's PropCovered
  ungroup()

# check that there are no duplicates in same year
inv_fwc2 %>%
  group_by(PermanentID, GSYear, SpeciesName) %>%
  count() %>%
  filter(n > 1)

# remove missing data
inv_fwc3 <- inv_fwc2 %>%
  filter(!is.na(PropCovered))

# save data
write_csv(inv_fwc3, "intermediate-data/FWC_hydrilla_pistia_eichhornia_survey_formatted.csv")


#### edit control data ####

# non-herbicide methods (from herbicide_initial_visualizations)
non_herb <- c("Mechanical Harvester", 
              "Snagging (tree removal)", 
              "Aquatic Dye (for shading)", 
              "Grass Carp", "Hand Removal", 
              "Mechanical (Other)", 
              "Mechanical Shredder", 
              "Prescribed Fire")

# duplicates per year/date?
ctrl_old %>%
  filter(Species %in% c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)") & TotalAcres > 0) %>%
  group_by(AreaOfInterestID, Year, Species) %>%
  summarise(apps = n(),
            acres = paste(TotalAcres, collapse = ", ")) %>%
  filter(apps > 1)
# yes, 106
# different acres in same year, probably different events

ctrl_new %>%
  filter(Species %in% c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)") & 
           TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
  group_by(AreaOfInterestID, BeginDate, Species) %>%
  summarise(apps = n(),
            acres = paste(TotalAcres, collapse = ", "),
            method = paste(ControlMethod, collapse = ", "),
            id = paste(TreatmentID, collapse = ", ")) %>%
  filter(apps > 1)
# yes, 5866
# for treatments involving more than one herbicide, separate lines for each herbicide type, but they share a treatment ID and acres
ctrl_new %>%
  group_by(TreatmentID) %>%
  summarise(acres = length(unique(TotalAcres)),
            species = length(unique(Species)),
            dates = length(unique(BeginDate))) %>%
  filter(acres > 1 | species > 1 | dates > 1)
# each treatment ID is a unique application for each species, area, and date
sum(is.na(ctrl_new$TreatmentID))

# overlap in datasets
ctrl_old %>%
  select(AreaOfInterestID, Year) %>%
  unique() %>%
  inner_join(ctrl_new %>%
               select(AreaOfInterestID, Year) %>%
               unique())
# 163 lakes
# potentially different events (some have different acres)

# average treatment month
ctrl_new  %>%
  filter(Species %in% c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)") & 
           TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>% # use same filters as ctrl_new2
  group_by(Species) %>%
  summarise(month = mean(month(BeginDate)))
# 7 for both

# old herbicide data
ctrl_old2 <- ctrl_old %>%
  filter(Species %in% c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)") & TotalAcres > 0) %>%
  mutate(AreaTreated_ha = TotalAcres * 0.405,
         Area_ha = ShapeArea * 100,
         AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha, # make full area if it exceeds it
                                    TRUE ~ AreaTreated_ha),
         PropTreated = AreaTreated_ha / Area_ha,
         TreatmentMethod = "unknown",
         TreatmentID = paste("old", Year, substr(Species, 1, 1), TotalAcres, sep = "_"),
         TreatmentYear = Year,
         CtrlSet = "old",
         TreatmentDate = as.Date(paste0(TreatmentYear, "-07-01")),
         TreatmentMonth = 7,
         GSYear = case_when(TreatmentMonth >= 4 ~ TreatmentYear,
                            TreatmentMonth < 4 ~ TreatmentYear - 1)) %>%
  select(AreaOfInterestID, PermanentID, TreatmentYear, Species, Area_ha, AreaTreated_ha, PropTreated, TreatmentMethod, TreatmentMonth, TreatmentDate, TreatmentID, CtrlSet, GSYear)

# new herbicide data
ctrl_new2 <- ctrl_new %>%
  filter(Species %in% c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)") & 
           TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>% # herbicide control only
  group_by(AreaOfInterestID, PermanentID, Species, BeginDate, TreatmentID, TotalAcres, ShapeArea) %>%  # captures area treated for an event without duplication due to multiple herbicides
  mutate(AreaTreated_ha = TotalAcres * 0.405,
         Area_ha = ShapeArea * 100,
         AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha, # make full area if it exceeds it
                                    TRUE ~ AreaTreated_ha),
         PropTreated = AreaTreated_ha / Area_ha,
         TreatmentMethod = paste(unique(ControlMethod), collapse = ", "),
         TreatmentYear = year(BeginDate),
         TreatmentMonth = month(BeginDate),
         TreatmentID = as.character(TreatmentID),
         CtrlSet = "new",
         GSYear = case_when(TreatmentMonth >= 4 ~ TreatmentYear,
                            TreatmentMonth < 4 ~ TreatmentYear - 1)) %>%
  ungroup() %>%
  select(AreaOfInterestID, PermanentID, TreatmentYear, Species, Area_ha, AreaTreated_ha, PropTreated, TreatmentMethod, TreatmentMonth, BeginDate, TreatmentID, CtrlSet, GSYear) %>%
  rename(TreatmentDate = BeginDate) %>%
  unique() # some duplication due to different TotalHerbicideUsed (but nothing else) for the same TreatmentID

# combine herbicide data
ctrl <- ctrl_old2 %>%
  full_join(ctrl_new2) %>%
  full_join(ctrl_old2 %>% # all possible measurements
              select(PermanentID) %>%
              unique() %>%
              expand_grid(GSYear = min(ctrl_old2$TreatmentYear):max(ctrl_old2$TreatmentYear)) %>% 
              full_join(ctrl_new2 %>%
                          select(PermanentID) %>%
                          unique() %>%
                          expand_grid(GSYear = min(ctrl_new2$TreatmentYear):max(ctrl_new2$TreatmentYear))) %>%
              expand_grid(Species = unique(ctrl_new2$Species))) %>% # one row for each species
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0),
         PropTreated = replace_na(PropTreated, 0)) %>% # no herbicide applied that year
  full_join(ctrl_old2 %>% # all ID's across all years
              select(PermanentID) %>%
              unique() %>%
              full_join(ctrl_new2 %>%
                          select(PermanentID) %>%
                          unique()) %>%
              expand_grid(GSYear = min(ctrl_old2$TreatmentYear):max(ctrl_new2$TreatmentYear)) %>%
              expand_grid(Species = unique(ctrl_new2$Species))) 
# leave AreaTreated_ha as NA for new rows from last join (lakes in old dataset not in new and vice versa)

# save data
write_csv(ctrl, "intermediate-data/FWC_hydrilla_pistia_eichhornia_herbicide_formatted.csv")


#### combine invasive plant and ctrl data ####

# modify control data to join invasion data
# if plant survey occurred before treatment, move treatment to following year
ctrl_inv <- ctrl %>%
  full_join(tibble(Species = c("Hydrilla verticillata", rep("Floating Plants (Eichhornia and Pistia)", 2)), # double each row that has floating plants
                   SpeciesName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"))) %>%
  left_join(inv_fwc3 %>%
              select(PermanentID, SpeciesName, SurveyDate, GSYear) %>% 
              unique()) %>% # add survey dates for each lake and year
  mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
         GSYear = case_when(is.na(SurveyTreatDays) | SurveyTreatDays >= 14 ~ GSYear, # no survey/survey after treatment -> keep year
                            SurveyTreatDays < 14 ~ GSYear + 1)) # survey before or very soon after treatment -> move treatment to next year

# treatments by year and lag
ctrl_inv2 <- ctrl_inv %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 0:5)) %>% # remove repeat row for each species
  mutate(DSet = "inv") %>%
  pmap(ctrl_lag_fun) %>% # summarizes ctrl_inv for each GS, Lag, PermID, and Sp
  bind_rows() %>%
  pivot_wider(names_from = Lag,
              values_from = c(PropTreated, Treated),
              names_glue = "Lag{Lag}{.value}") # make treatments wide by lag

# combine treatment and invasion datasets
inv_ctrl <- inv_fwc3 %>%
  inner_join(ctrl_inv2) %>% # only include data from both datasets
  group_by(SpeciesName) %>%
  mutate(PropCoveredBeta = transform01(PropCovered)) %>% # uses sample size within species, leaving out NA's
  ungroup() %>%
  filter(!is.na(PrevPropCovered)) #  missing initial pop abundance

# check for missing data
inv_fwc3 %>%
  select(PermanentID, SpeciesName, GSYear) %>%
  unique() %>%
  filter(GSYear >= 1998) %>%
  anti_join(inv_ctrl %>%
              select(PermanentID, SpeciesName, GSYear) %>%
              unique()) %>%
  select(PermanentID) %>%
  unique() %>% # 210 total
  inner_join(ctrl %>%
              select(PermanentID) %>%
              unique()) # 140
# treatment data missing for some years for 140 lakes
# 70 lakes not in ctrl dataset with overlapping time

# save data
write_csv(inv_ctrl, "intermediate-data/FWC_hydrilla_pistia_eichhornia_survey_herbicide_formatted.csv")


#### edit native plant data ####

# one origin per species?
plant_fwc %>%
  group_by(SpeciesName) %>%
  summarise(orig = length(unique(Origin))) %>%
  ungroup() %>%
  filter(orig != 1)
# yes

# Eppc and origin
plant_fwc %>%
  select(Eppc, Origin) %>%
  unique()

# first year detected
nat_first_detect <- plant_fwc %>%
  filter(Origin == "Native") %>%
  group_by(PermanentID, SpeciesName) %>%
  summarise(FirstDetect = min(SurveyDate)) %>%
  ungroup()

# select native species
# used clustering algorithm in Open Refine to look for mispelled names (none)
nat_fwc <- plant_fwc %>%
  filter(Origin == "Native") %>%
  select(SpeciesName, Habitat, HabitatShortName) %>%
  unique() %>% # full native species list
  expand_grid(plant_fwc %>%
                select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, Surveyor) %>%
                unique()) %>% # full survey list (row for every species in every survey)
  full_join(plant_fwc %>%
              filter(Origin == "Native") %>%
              select(SpeciesName, Habitat, HabitatShortName, 
                     AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, Surveyor, 
                     IsDetected)) %>% # add detection data (only "Yes")
  full_join(nat_first_detect) %>% # add first detection date
  left_join(surveyor2) %>% # add surveyor experience
  mutate(IsDetected = replace_na(IsDetected, "No"),
         Detected = case_when(IsDetected == "Yes" ~ 1,
                              IsDetected == "No" ~ 0),
         FirstDetect = replace_na(FirstDetect, as.Date("2021-06-24")),  # use today's date if not detected
         SurveyMonth = month(SurveyDate),
         GSYear = case_when(SurveyMonth >= 4 ~ year(SurveyDate),
                            SurveyMonth < 4 ~ year(SurveyDate) - 1),
         Area_ha = ShapeArea * 100) # convert lake area from km-squared to hectares

# duplicate surveys in a year
nat_fwc %>%
  group_by(AreaOfInterestID, GSYear) %>%
  summarise(surveys = length(unique(SurveyDate))) %>%
  ungroup() %>%
  filter(surveys > 1) 
# 41 AOIs have multiple surveys in a year

# check that FirstDetect is accurate
nat_fwc %>%
  filter(SurveyDate < FirstDetect) %>%
  select(IsDetected) %>%
  unique()
# all are not detected

# variable mean/sd for centering/scaling
nat_interval <- nat_fwc %>%
  group_by(PermanentID, GSYear) %>% # remove redundancy from each species in survey
  summarise(SurveyDate = max(SurveyDate)) %>%
  ungroup() %>%
  group_by(PermanentID) %>%
  arrange(SurveyDate) %>%
  mutate(SurveyInterval = as.numeric(SurveyDate - lag(SurveyDate))) %>%
  ungroup() %>%
  summarise(MeanInterval = mean(SurveyInterval, na.rm = T)) %>%
  pull(MeanInterval)

nat_area <- nat_fwc %>%
  select(PermanentID, Area_ha) %>%
  unique() %>% # remove redundancy from each survey and species in each lake
  summarise(MeanArea = mean(Area_ha),
            SDArea = sd(Area_ha))

# initial immigration dataset
nat_init_imm <- nat_fwc %>%
  filter(SurveyDate <= FirstDetect) %>% # remove surveys after first detection
  group_by(PermanentID, Area_ha, GSYear, SpeciesName, Habitat, HabitatShortName) %>% # summarize over multiple surveys per year and lake
  summarise(AreaName = paste(AreaOfInterest, collapse = "/"),
            Detected = as.numeric(sum(Detected) > 0), # was species detected that growing season?
            SurveyDate = max(SurveyDate), # last survey each growing season
            SurveyorExperience = mean(SurveyorExperience, na.rm = T),
            SurveyorExperienceCS = mean(SurveyorExperienceCS, na.rm = T)) %>%
  ungroup() %>%
  group_by(PermanentID, SpeciesName) %>%
  arrange(SurveyDate) %>%
  mutate(SurveyInterval = as.numeric(SurveyDate - lag(SurveyDate)), # days between surveys
         FirstGS = min(GSYear)) %>%
  ungroup() %>%
  filter(GSYear > FirstGS) # remove first year (no interval)

# ext/repeat imm dataset
nat_post_imm <- nat_fwc %>%
  filter(SurveyDate >= FirstDetect) %>% # remove surveys before first detection
  group_by(PermanentID, Area_ha, GSYear, SpeciesName, Habitat, HabitatShortName) %>% # summarize over multiple surveys per year and lake
  summarise(AreaName = paste(AreaOfInterest, collapse = "/"),
            Detected = as.numeric(sum(Detected) > 0), # was species detected that growing season?
            SurveyDate = max(SurveyDate), # last survey each growing season
            SurveyorExperience = mean(SurveyorExperience, na.rm = T),
            SurveyorExperienceCS = mean(SurveyorExperienceCS, na.rm = T)) %>%
  ungroup() %>%
  group_by(PermanentID, SpeciesName) %>%
  arrange(GSYear) %>%
  mutate(FirstGS = min(GSYear),
         Extinct = case_when(lag(Detected) == 1 & Detected == 1 ~ 0,
                             lag(Detected) == 1 & Detected == 0 ~ 1,
                             lag(Detected) == 0 ~ NA_real_,
                             is.na(lag(Detected)) ~ NA_real_),
         RepeatImm = case_when(lag(Detected) == 0 & Detected == 0 ~ 0,
                               lag(Detected) == 0 & Detected == 1 ~ 1,
                               lag(Detected) == 1 ~ NA_real_,
                               is.na(lag(Detected)) ~ NA_real_),
         Switch = case_when(lag(Detected) == 0 & Detected == 1 ~ 1,
                            lag(Detected) == 1 & Detected == 0 ~ 1,
                            GSYear == FirstGS ~ 1, # first detection
                            TRUE ~ 0),
         LagSwitch = lag(Switch),
         LagSwitch = replace_na(LagSwitch, 0),
         Window = cumsum(LagSwitch),
         SurveyInterval = as.numeric(SurveyDate - lag(SurveyDate))) %>%
  ungroup()

# extinction dataset
nat_ext <- nat_post_imm %>%
  filter(!is.na(Extinct)) %>%
  select(-c(Switch, LagSwitch))

# check that it worked
unique(nat_ext$Window)
# all odd numbers

# repeat imm dataset
nat_rept_imm <- nat_post_imm %>%
  filter(!is.na(RepeatImm)) %>%
  select(-c(Switch, LagSwitch))

# check that it worked
unique(nat_rept_imm$Window)
# all even numbers


#### combine native, invasive, ctrl data ####

# modify control data to join native datasets
# if plant survey occurred before treatment, move treatment to following year
ctrl_nat_init_imm <- ctrl %>%
  left_join(nat_init_imm %>%
              select(PermanentID, SpeciesName, GSYear, SurveyDate) %>%
              unique()) %>% # add survey dates for each lake and year
  mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
         GSYear = case_when(SurveyTreatDays >= 14 | is.na(SurveyTreatDays) ~ GSYear, # survey after treatment/no survey -> keep year
                            SurveyTreatDays < 14 ~ GSYear + 1)) # survey before or very soon after treatment -> move treatment to next year

ctrl_nat_ext <- ctrl %>%
  left_join(nat_ext %>%
              select(PermanentID, SpeciesName, GSYear, SurveyDate) %>%
              unique()) %>% # add survey dates for each lake and year
  mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
         GSYear = case_when(SurveyTreatDays >= 14 | is.na(SurveyTreatDays) ~ GSYear, # survey after treatment/no survey -> keep year
                            SurveyTreatDays < 14 ~ GSYear + 1)) # survey before or very soon after treatment -> move treatment to next year

ctrl_nat_rept_imm <- ctrl %>%
  left_join(nat_rept_imm %>%
              select(PermanentID, SpeciesName, GSYear, SurveyDate) %>%
              unique()) %>% # add survey dates for each lake and year
  mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
         GSYear = case_when(SurveyTreatDays >= 14 | is.na(SurveyTreatDays) ~ GSYear, # survey after treatment/no survey -> keep year
                            SurveyTreatDays < 14 ~ GSYear + 1)) # survey before or very soon after treatment -> move treatment to next year

# treatments by year and lag
ctrl_nat_init_imm2 <- ctrl_nat_init_imm %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 0:5)) %>% # remove repeat row for each species
  mutate(DSet = "nat_init_imm") %>%
  pmap(ctrl_lag_fun) %>% # summarizes ctrl_nat_init_imm for each GS, Lag, PermID, and Sp
  bind_rows() %>%
  mutate(Species = fct_recode(Species, "F" = "Floating Plants (Eichhornia and Pistia)",
                              "H" = "Hydrilla verticillata")) %>%
  pivot_wider(names_from = c(Lag, Species),
              values_from = c(PropTreated, Treated),
              names_glue = "Lag{Lag}{.value}{Species}") # make treatments wide by lag

ctrl_nat_ext2 <- ctrl_nat_ext %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 0:5)) %>% # remove repeat row for each species
  mutate(DSet = "nat_ext") %>%
  pmap(ctrl_lag_fun) %>% # summarizes ctrl_nat_ext for each GS, Lag, PermID, and Sp
  bind_rows() %>%
  mutate(Species = fct_recode(Species, "F" = "Floating Plants (Eichhornia and Pistia)",
                              "H" = "Hydrilla verticillata")) %>%
  pivot_wider(names_from = c(Lag, Species),
              values_from = c(PropTreated, Treated),
              names_glue = "Lag{Lag}{.value}{Species}") # make treatments wide by lag

ctrl_nat_rept_imm2 <- ctrl_nat_rept_imm %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 0:5)) %>% # remove repeat row for each species
  mutate(DSet = "nat_rept_imm") %>%
  pmap(ctrl_lag_fun) %>% # summarizes ctrl_nat_init_imm for each GS, Lag, PermID, and Sp
  bind_rows() %>%
  mutate(Species = fct_recode(Species, "F" = "Floating Plants (Eichhornia and Pistia)",
                              "H" = "Hydrilla verticillata")) %>%
  pivot_wider(names_from = c(Lag, Species),
              values_from = c(PropTreated, Treated),
              names_glue = "Lag{Lag}{.value}{Species}") # make treatments wide by lag

# invasive abundance by year and lag
inv_nat_init_imm <- nat_init_imm %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 0:5)) %>%
  pmap(inv_lag_fun) %>% # summarizes inv_fwc2 for each GS, Lag, PermID, and Sp
  bind_rows() %>%
  pivot_wider(names_from = Lag,
              values_from = c(Hydrilla, WaterLettuce, WaterHyacinth),
              names_glue = "Lag{Lag}{.value}") # make treatments wide by lag

inv_nat_ext <- nat_ext %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 0:5)) %>%
  pmap(inv_lag_fun) %>% # summarizes inv_fwc2 for each GS, Lag, PermID, and Sp
  bind_rows() %>%
  pivot_wider(names_from = Lag,
              values_from = c(Hydrilla, WaterLettuce, WaterHyacinth),
              names_glue = "Lag{Lag}{.value}") # make treatments wide by lag

inv_nat_rept_imm <- nat_rept_imm %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 0:5)) %>%
  pmap(inv_lag_fun) %>% # summarizes inv_fwc2 for each GS, Lag, PermID, and Sp
  bind_rows() %>%
  pivot_wider(names_from = Lag,
              values_from = c(Hydrilla, WaterLettuce, WaterHyacinth),
              names_glue = "Lag{Lag}{.value}") # make treatments wide by lag

# combine dataset for initial immigration
nat_init_imm2 <- nat_init_imm %>%
  inner_join(ctrl_nat_init_imm2) %>%
  left_join(inv_nat_init_imm) %>%
  mutate(SurveyIntervalC = SurveyInterval - nat_interval, # center/scale variables
         Area_haCS = (Area_ha - nat_area$MeanArea) / nat_area$SDArea)

nat_ext2 <- nat_ext %>%
  inner_join(ctrl_nat_ext2) %>%
  left_join(inv_nat_ext) %>%
  mutate(SurveyIntervalC = SurveyInterval - nat_interval, # center/scale variables
         Area_haCS = (Area_ha - nat_area$MeanArea) / nat_area$SDArea)

nat_rept_imm2 <- nat_rept_imm %>%
  inner_join(ctrl_nat_rept_imm2) %>%
  left_join(inv_nat_rept_imm) %>%
  mutate(SurveyIntervalC = SurveyInterval - nat_interval, # center/scale variables
         Area_haCS = (Area_ha - nat_area$MeanArea) / nat_area$SDArea)

# save data
write_csv(nat_init_imm2, "intermediate-data/FWC_native_initial_immigration_invasion_herbicide_formatted.csv")
write_csv(nat_ext2, "intermediate-data/FWC_native_extinction_invasion_herbicide_formatted.csv")
write_csv(nat_rept_imm2, "intermediate-data/FWC_native_repeat_immigration_invasion_herbicide_formatted.csv")


#### hydrilla models ####

# subset for hydrilla
hydr_ctrl <- inv_ctrl %>%
  filter(CommonName == "Hydrilla")

# subset for hydrilla and all herbicide lags
hydr_ctrl2 <- hydr_ctrl %>%
  filter(!is.na(Lag0PropTreated) & !is.na(Lag1PropTreated) & !is.na(Lag2PropTreated) & !is.na(Lag3PropTreated) & !is.na(Lag4PropTreated) & !is.na(Lag5PropTreated) & !is.na(SurveyorExperienceCS))

# figures
hydr_ctrl2 %>% ggplot(aes(PropCoveredBeta)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(PrevPropCovered)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(SurveyorExperienceCS)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag0PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag1PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag2PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag3PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag4PropTreated)) + geom_histogram(binwidth = 0.1)
hydr_ctrl2 %>% ggplot(aes(Lag5PropTreated)) + geom_histogram(binwidth = 0.1)

# random effects
length(unique(hydr_ctrl2$PermanentID))
length(unique(hydr_ctrl2$GSYear))

# models couldn't converge with surveyor experience - add in after treatment form is chosen

# models
hydr_prop0_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag0PropTreated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag0PropTreated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_prop1_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag1PropTreated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag1PropTreated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_prop2_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag2PropTreated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag2PropTreated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_prop3_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag3PropTreated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag3PropTreated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_prop4_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag4PropTreated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag4PropTreated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_prop5_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag5PropTreated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag5PropTreated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_trtd0_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag0Treated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag0Treated:SurveyorExperienceCS + (1|PermanentID) + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_trtd1_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag1Treated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag1Treated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_trtd2_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag2Treated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag2Treated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_trtd3_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag3Treated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag3Treated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_trtd4_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag4Treated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag4Treated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)
hydr_trtd5_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag5Treated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag5Treated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl2)

# compare models
AIC(hydr_prop0_mod, hydr_prop1_mod, hydr_prop2_mod, hydr_prop3_mod, hydr_prop4_mod, hydr_prop5_mod,
    hydr_trtd0_mod, hydr_trtd1_mod, hydr_trtd2_mod, hydr_trtd3_mod, hydr_trtd4_mod, hydr_trtd5_mod) %>%
  arrange(AIC) %>%
  mutate(deltaAIC = AIC - min(AIC))

# best model
summary(hydr_prop5_mod)
summary(hydr_prop1_mod) # surveyor exp. still sig

# use largest possible dataset
hydr_ctrl3 <- hydr_ctrl %>%
  filter(!is.na(Lag5PropTreated) & !is.na(SurveyorExperienceCS))

# figures
hydr_ctrl3 %>% ggplot(aes(PropCoveredBeta)) + geom_histogram(binwidth = 0.1)
hydr_ctrl3 %>% ggplot(aes(PrevPropCovered)) + geom_histogram(binwidth = 0.1)
hydr_ctrl3 %>% ggplot(aes(SurveyorExperienceCS)) + geom_histogram(binwidth = 0.1)
hydr_ctrl3 %>% ggplot(aes(Lag5PropTreated)) + geom_histogram(binwidth = 0.1)

# refit model
hydr_mod <- glmmTMB(PropCoveredBeta ~ PrevPropCovered + PrevPropCovered:Lag5PropTreated + PrevPropCovered:SurveyorExperienceCS + PrevPropCovered:Lag5PropTreated:SurveyorExperienceCS + (1|PermanentID) + (1|GSYear), family = beta_family(), data = hydr_ctrl)
summary(hydr_mod)

# figures
plot(predictorEffect("PrevPropCovered", hydr_mod))
# avg surveyor and not treatment:
  # increasing treatment decreases
#### start here ####
# I think it's problematic that surveyor experience is negative



#### older analyses below ####

#### summary stats ####

# overall invasion dataset
inv_ctrl %>%
  summarise(Lakes = length(unique(PermanentID)),
            Years = length(unique(GSYear)))

# detected lakes
lake_sum <- plant_ctrl_2yearAlt %>%
  group_by(CommonName) %>%
  summarise(Lakes = length(unique(AreaOfInterestID)),
            YearLakes = n()) %>%
  mutate(NameLakes = paste(CommonName, "\n(", Lakes, " lakes, ", YearLakes, " data points)", sep = ""))

# year ranges
pdf("output/year_ranges_herbicide_analysis.pdf", width = 11.5, height = 3)
plant_ctrl_2yearAlt %>%
  group_by(CommonName, AreaOfInterestID) %>%
  summarise(Years = max(SurveyYear) - min(SurveyYear) + 1) %>%
  ungroup() %>%
  left_join(lake_sum) %>%
  ggplot(aes(x = Years)) +
  geom_bar() +
  facet_wrap(~ NameLakes) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  labs(x = "Years of data", y = "Number of lakes") +
  def_theme
dev.off()

# population sizes
prop_sum <- plant_ctrl_2yearAlt %>%
  group_by(CommonName) %>%
  summarise(mean = mean(PropCoveredAdj))

pdf("output/pop_sizes_herbicide_analysis.pdf", width = 11.5, height = 2.5)
plant_ctrl_2yearAlt %>%
  ggplot(aes(PropCoveredAdj)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(data = prop_sum, aes(xintercept = mean), color = "blue", linetype = "dashed") +
  geom_text(data = prop_sum, aes(x =  mean, label = paste("mean = ", as.character(round(mean, 3)), sep = "")), 
            y = 6150, hjust = -0.05, size = 4, color = "blue") +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  labs(x = "Proportion of lake inhabited", y = "Data points") +
  ylim(0, 6200) +
  def_theme +
  theme(strip.text = element_blank())
dev.off()

# survey times
pdf("output/survey_times_1year_herbicide_analysis.pdf", width = 11.5, height = 2.75)
plant_ctrl_1yearAlt %>%
  ggplot(aes(x = Next1SurveyDays)) + 
  geom_histogram(binwidth = 10) +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  labs(x = expression(paste("Time between surveys for ", tau, " = 1 [days]", sep = "")),
       y = "Data points") +
  def_theme
dev.off()

pdf("output/survey_times_2years_herbicide_analysis.pdf", width = 11.5, height = 2.5)
plant_ctrl_2yearAlt %>%
  ggplot(aes(x = Next2SurveyDays)) + 
  geom_histogram(binwidth = 10) +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  labs(x = expression(paste("Time between surveys for ", tau, " = 2 [days]", sep = "")),
       y = "Data points") +
  def_theme +
  theme(strip.text = element_blank())
dev.off()

pdf("output/treatment_survey_times_herbicide_analysis.pdf", width = 11.5, height = 2.5)
plant_ctrl_2yearAlt %>%
  filter(!is.na(TreatmentDate)) %>%
  ggplot(aes(x = TreatSurveyDays)) + 
  geom_histogram(binwidth = 10) +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  labs(x = "Time from survey to treatment [days]",
       y = "Data points") +
  def_theme +
  theme(strip.text = element_blank())
dev.off()

pdf("output/pop_changes_1year_herbicide_analysis.pdf", width = 11.5, height = 2.75)
plant_ctrl_1yearAlt %>%
  ggplot(aes(log_Diff1PropCovered)) +
  geom_histogram() +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Change in population size [ln(", N[t+1], "/", N[t], ")]", sep = "")),
       y = "Data points") +
  def_theme
dev.off()

pdf("output/pop_changes_2year_herbicide_analysis.pdf", width = 11.5, height = 2.75)
plant_ctrl_2yearAlt %>%
  ggplot(aes(log_Diff2PropCovered)) +
  geom_histogram() +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Change in population size [ln(", N[t+2], "/", N[t], ")]", sep = "")),
       y = "Data points") +
  def_theme
dev.off()
# similar distribution, less data

pdf("output/initial_pop_herbicide_analysis.pdf", width = 11.5, height = 2.5)
plant_ctrl_2yearAlt %>%
  ggplot(aes(log_PropCoveredC)) +
  geom_histogram() +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Initial population size [ln(", N[t], "/A) - ", mu, "]", sep = "")),
       y = "Data points") +
  def_theme +
  theme(strip.text = element_blank())
dev.off()

pdf("output/treatment_herbicide_analysis.pdf", width = 11.5, height = 2.5)
plant_ctrl_2yearAlt %>%
  ggplot(aes(x = TreatedF)) +
  geom_bar() +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = "Treatment applied in initial year",
       y = "Data points") +
  def_theme +
  theme(strip.text = element_blank())
dev.off()

pdf("output/future_pop_1year_herbicide_analysis.pdf", width = 11.5, height = 2.75)
plant_ctrl_1yearAlt %>%
  ggplot(aes(x = Next1PropCoveredBeta)) +
  geom_histogram() +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Future population size [", N[t+1], "/A]", sep = "")),
       y = "Data points") +
  def_theme
dev.off()

plant_ctrl_2yearAlt %>%
  ggplot(aes(x = Next2PropCoveredBeta)) +
  geom_histogram() +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Future population size [", N[t+1], "/A]", sep = "")),
       y = "Data points") +
  def_theme

pdf("output/initial_pop_beta_herbicide_analysis.pdf", width = 11.5, height = 2.5)
plant_ctrl_2yearAlt %>%
  ggplot(aes(PropCoveredBetaC)) +
  geom_histogram() +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Initial population size [", N[t], "/A - ", mu, "]", sep = "")),
       y = "Data points") +
  def_theme +
  theme(strip.text = element_blank())
dev.off()


#### divide data ####

# hydrilla 1-year growth
hyd_dat_1year <- plant_ctrl_1year %>%
  filter(CommonName == "Hydrilla")

# water hyacinth 1-year growth
why_dat_1year <- plant_ctrl_1year %>%
  filter(CommonName == "Water hyacinth")

# water lettuce 1-year growth
wle_dat_1year <- plant_ctrl_1year %>%
  filter(CommonName == "Water lettuce")

# hydrilla 2-year growth
hyd_dat_2year <- plant_ctrl_2year %>%
  filter(CommonName == "Hydrilla")

# water hyacinth 2-year growth
why_dat_2year <- plant_ctrl_2year %>%
  filter(CommonName == "Water hyacinth")

# water lettuce 2-year growth
wle_dat_2year <- plant_ctrl_2year %>%
  filter(CommonName == "Water lettuce")

# alt hydrilla 1-year growth
hyd_dat_1yearAlt <- plant_ctrl_1yearAlt %>%
  filter(CommonName == "Hydrilla")

# alt water hyacinth 1-year growth
why_dat_1yearAlt <- plant_ctrl_1yearAlt %>%
  filter(CommonName == "Water hyacinth")

# alt water lettuce 1-year growth
wle_dat_1yearAlt <- plant_ctrl_1yearAlt %>%
  filter(CommonName == "Water lettuce")

# alt hydrilla 2-year growth
hyd_dat_2yearAlt <- plant_ctrl_2yearAlt %>%
  filter(CommonName == "Hydrilla")

# alt water hyacinth 2-year growth
why_dat_2yearAlt <- plant_ctrl_2yearAlt %>%
  filter(CommonName == "Water hyacinth")

# alt water lettuce 2-year growth
wle_dat_2yearAlt <- plant_ctrl_2yearAlt %>%
  filter(CommonName == "Water lettuce")


#### fit hydrilla models ####

# 1-year growth
hyd_mod_1year <- glmmTMB(log_Diff1PropCovered ~ log_PropCoveredC * Treated +
                           (1|AreaOfInterestID) + (1|SurveyYear), 
                         data = hyd_dat_1year)
summary(hyd_mod_1year)

# 2-year growth
hyd_mod_2year <- glmmTMB(log_Diff2PropCovered ~ log_PropCoveredC * Treated +
                           (1|AreaOfInterestID) + (1|SurveyYear), 
                         data = hyd_dat_2year)
summary(hyd_mod_2year)

# alt 1-year growth
hyd_mod_1yearAlt <- update(hyd_mod_1year, data = hyd_dat_1yearAlt)
summary(hyd_mod_1yearAlt)

# alt 2-year growth
hyd_mod_2yearAlt <- update(hyd_mod_2year, data = hyd_dat_2yearAlt)
summary(hyd_mod_2yearAlt)
# only model with marginal interaction and all had apositive treatment effects


#### fit water hyacinth models ####

# 1-year growth
why_mod_1year <- glmmTMB(log_Diff1PropCovered ~ log_PropCoveredC * Treated +
                           (1|AreaOfInterestID) + (1|SurveyYear), 
                         data = why_dat_1year)
summary(why_mod_1year)

# 2-year growth
why_mod_2year <- glmmTMB(log_Diff2PropCovered ~ log_PropCoveredC * Treated +
                           (1|AreaOfInterestID) + (1|SurveyYear), 
                         data = why_dat_2year)
summary(why_mod_2year)

# alt 1-year growth
why_mod_1yearAlt <- update(why_mod_1year, data = why_dat_1yearAlt)
summary(why_mod_1yearAlt)

# alt 2-year growth
why_mod_2yearAlt <- update(why_mod_2year, data = why_dat_2yearAlt)
summary(why_mod_2yearAlt)

# no models had sig interactions and all had positive treatment effects


#### fit water lettuce models ####

# 1-year growth
wle_mod_1year <- glmmTMB(log_Diff1PropCovered ~ log_PropCoveredC * Treated +
                           (1|AreaOfInterestID) + (1|SurveyYear), 
                         data = wle_dat_1year)
summary(wle_mod_1year)

# 2-year growth
wle_mod_2year <- glmmTMB(log_Diff2PropCovered ~ log_PropCoveredC * Treated +
                           (1|AreaOfInterestID) + (1|SurveyYear), 
                         data = wle_dat_2year)
summary(wle_mod_2year)

# alt 1-year growth
wle_mod_1yearAlt <- update(wle_mod_1year, data = wle_dat_1yearAlt)
summary(wle_mod_1yearAlt)

# alt 2-year growth
wle_mod_2yearAlt <- update(wle_mod_2year, data = wle_dat_2yearAlt)
summary(wle_mod_2yearAlt)

# no models had sig interactions and all had positive treatment effects


#### visualize 1-year alt models ####

# add predicted values
plant_ctrl_pred1 <- hyd_dat_1yearAlt %>%
  mutate(pred = predict(hyd_mod_1yearAlt, newdata = ., re.form = NA),
         pred_se = predict(hyd_mod_1yearAlt, newdata = ., re.form = NA, se.fit = T)$se.fit) %>%
  full_join(why_dat_1yearAlt %>%
              mutate(pred = predict(why_mod_1yearAlt, newdata = ., re.form = NA),
                     pred_se = predict(why_mod_1yearAlt, newdata = ., re.form = NA, se.fit = T)$se.fit)) %>%
  full_join(wle_dat_1yearAlt %>%
              mutate(pred = predict(wle_mod_1yearAlt, newdata = ., re.form = NA),
                     pred_se = predict(wle_mod_1yearAlt, newdata = ., re.form = NA, se.fit = T)$se.fit))

# figure
pdf("output/model_fit_1year_herbicide_analysis.pdf", width = 11.5, height = 4)
ggplot(plant_ctrl_pred1, aes(x = log_PropCoveredC, y = log_Diff1PropCovered)) +
  geom_point(alpha = 0.2, aes(color = TreatedF)) +
  geom_ribbon(aes(y = pred, ymin = pred - pred_se, ymax = pred + pred_se, fill = TreatedF), 
              color = NA, alpha = 0.5) +
  geom_line(aes(y = pred, color = TreatedF)) +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Initial population size [ln(", N[t], "/A) - ", mu, "]", sep = "")),
       y = expression(paste("Change in population size [ln(", N[t+1], "/", N[t], ")]", sep = ""))) +
  scale_color_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  scale_fill_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  def_theme +
  theme(legend.position = c(0.97, 0.87))
dev.off()

pdf("output/model_fit_1year_untransformed_herbicide_analysis.pdf", width = 11.5, height = 4)
ggplot(plant_ctrl_pred1, aes(x = PropCoveredAdj, y = log_Diff1PropCovered)) +
  geom_point(alpha = 0.2, aes(color = TreatedF)) +
  geom_ribbon(aes(y = pred, ymin = pred - pred_se, ymax = pred + pred_se, fill = TreatedF), 
              color = NA, alpha = 0.5) +
  geom_line(aes(y = pred, color = TreatedF)) +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Initial population size [", N[t], "/A]", sep = "")),
       y = expression(paste("Change in population size [ln(", N[t+1], "/", N[t], ")]", sep = ""))) +
  scale_color_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  scale_fill_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  def_theme +
  theme(legend.position = c(0.97, 0.87))
dev.off()


#### visualize 2-year alt models ####

# add predicted values
plant_ctrl_pred2 <- hyd_dat_2yearAlt %>%
  mutate(pred = predict(hyd_mod_2yearAlt, newdata = ., re.form = NA),
         pred_se = predict(hyd_mod_2yearAlt, newdata = ., re.form = NA, se.fit = T)$se.fit) %>%
  full_join(why_dat_2yearAlt %>%
              mutate(pred = predict(why_mod_2yearAlt, newdata = ., re.form = NA),
                     pred_se = predict(why_mod_2yearAlt, newdata = ., re.form = NA, se.fit = T)$se.fit)) %>%
  full_join(wle_dat_2yearAlt %>%
              mutate(pred = predict(wle_mod_2yearAlt, newdata = ., re.form = NA),
                     pred_se = predict(wle_mod_2yearAlt, newdata = ., re.form = NA, se.fit = T)$se.fit))

# figure
pdf("output/model_fit_2year_herbicide_analysis.pdf", width = 11.5, height = 3.5)
ggplot(plant_ctrl_pred2, aes(x = log_PropCoveredC, y = log_Diff2PropCovered)) +
  geom_point(alpha = 0.2, aes(color = TreatedF)) +
  geom_ribbon(aes(y = pred, ymin = pred - pred_se, ymax = pred + pred_se, fill = TreatedF), 
              color = NA, alpha = 0.5) +
  geom_line(aes(y = pred, color = TreatedF)) +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Initial population size [ln(", N[t], "/A) - ", mu, "]", sep = "")),
       y = expression(paste("Change in population size [ln(", N[t+2], "/", N[t], ")]", sep = ""))) +
  scale_color_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  scale_fill_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  def_theme +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.title.y = element_text(size = 14, color="black", hjust = 0.9))
dev.off()
  
pdf("output/model_fit_2year_untransformed_herbicide_analysis.pdf", width = 11.5, height = 3.5)
ggplot(plant_ctrl_pred2, aes(x = PropCoveredAdj, y = log_Diff2PropCovered)) +
  geom_point(alpha = 0.2, aes(color = TreatedF)) +
  geom_ribbon(aes(y = pred, ymin = pred - pred_se, ymax = pred + pred_se, fill = TreatedF), 
              color = NA, alpha = 0.5) +
  geom_line(aes(y = pred, color = TreatedF)) +
  #geom_hline(yintercept = 0) +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Initial population size [", N[t], "/A]", sep = "")),
       y = expression(paste("Change in population size [ln(", N[t+2], "/", N[t], ")]", sep = ""))) +
  scale_color_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  scale_fill_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  def_theme +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.title.y = element_text(size = 14, color="black", hjust = 0.9))
dev.off()


#### fit hydrilla beta models ####

# 1-year growth
hyd_beta_mod_1year <- glmmTMB(Next1PropCoveredBeta ~ PropCoveredBetaC * Treated +
                           (1|AreaOfInterestID) + (1|SurveyYear), 
                         family = beta_family(),
                         data = hyd_dat_1year)
summary(hyd_beta_mod_1year) # sig negative interaction

# 2-year growth
hyd_beta_mod_2year <- glmmTMB(Next2PropCoveredBeta ~ PropCoveredBetaC * Treated +
                           (1|AreaOfInterestID) + (1|SurveyYear), 
                         data = hyd_dat_2year)
summary(hyd_beta_mod_2year)

# alt 1-year growth
hyd_beta_mod_1yearAlt <- update(hyd_beta_mod_1year, data = hyd_dat_1yearAlt)
summary(hyd_beta_mod_1yearAlt) # sig negative interaction

# alt 2-year growth
hyd_beta_mod_2yearAlt <- update(hyd_beta_mod_2year, data = hyd_dat_2yearAlt)
summary(hyd_beta_mod_2yearAlt)


#### fit water hyacinth beta models ####

# 1-year growth
why_beta_mod_1year <- glmmTMB(Next1PropCoveredBeta ~ PropCoveredBetaC * Treated +
                                (1|AreaOfInterestID) + (1|SurveyYear), 
                              family = beta_family(),
                              data = why_dat_1year)
summary(why_beta_mod_1year) # sig negative interaction

# 2-year growth
why_beta_mod_2year <- glmmTMB(Next2PropCoveredBeta ~ PropCoveredBetaC * Treated +
                                (1|AreaOfInterestID) + (1|SurveyYear), 
                              data = why_dat_2year)
summary(why_beta_mod_2year)

# alt 1-year growth
why_beta_mod_1yearAlt <- update(why_beta_mod_1year, data = why_dat_1yearAlt)
summary(why_beta_mod_1yearAlt) # sig negative interaction

# alt 2-year growth
why_beta_mod_2yearAlt <- update(why_beta_mod_2year, data = why_dat_2yearAlt)
summary(why_beta_mod_2yearAlt)


#### fit water lettuce beta models ####

# 1-year growth
wle_beta_mod_1year <- glmmTMB(Next1PropCoveredBeta ~ PropCoveredBetaC * Treated +
                                (1|AreaOfInterestID) + (1|SurveyYear), 
                              family = beta_family(),
                              data = wle_dat_1year)
summary(wle_beta_mod_1year) # sig positive interaction

# 2-year growth
wle_beta_mod_2year <- glmmTMB(Next2PropCoveredBeta ~ PropCoveredBetaC * Treated +
                                (1|AreaOfInterestID) + (1|SurveyYear), 
                              data = wle_dat_2year)
summary(wle_beta_mod_2year) # sig positive interaction

# alt 1-year growth
wle_beta_mod_1yearAlt <- update(wle_beta_mod_1year, data = wle_dat_1yearAlt)
summary(wle_beta_mod_1yearAlt) # no interaction

# alt 2-year growth
wle_beta_mod_2yearAlt <- update(wle_beta_mod_2year, data = wle_dat_2yearAlt)
summary(wle_beta_mod_2yearAlt) # sig negative interaction


#### visualize 1-year alt beta models ####

# add predicted values
plant_ctrl_beta_pred1 <- hyd_dat_1yearAlt %>%
  mutate(pred = predict(hyd_beta_mod_1yearAlt, newdata = ., re.form = NA, type = "response"),
         pred_se = predict(hyd_beta_mod_1yearAlt, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit) %>%
  full_join(why_dat_1yearAlt %>%
              mutate(pred = predict(why_beta_mod_1yearAlt, newdata = ., re.form = NA, type = "response"),
                     pred_se = predict(why_beta_mod_1yearAlt, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)) %>%
  full_join(wle_dat_1yearAlt %>%
              mutate(pred = predict(wle_beta_mod_1yearAlt, newdata = ., re.form = NA, type = "response"),
                     pred_se = predict(wle_beta_mod_1yearAlt, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit))

# figure
pdf("output/beta_model_fit_1year_herbicide_analysis.pdf", width = 11.5, height = 4)
ggplot(plant_ctrl_beta_pred1, aes(x = PropCoveredBeta, y = Next1PropCoveredBeta)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.2, aes(color = TreatedF)) +
  geom_ribbon(aes(y = pred, ymin = pred - pred_se, ymax = pred + pred_se, fill = TreatedF), 
              color = NA, alpha = 0.5) +
  geom_line(aes(y = pred, color = TreatedF)) +
  facet_wrap(~ CommonName, scales = "free") +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Initial population size [", N[t], "/A]", sep = "")),
       y = expression(paste("Future population size [", N[t+1], "/A]", sep = ""))) +
  scale_color_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  scale_fill_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  def_theme +
  theme(legend.position = c(0.45, 0.77))
dev.off()

pdf("output/beta_model_only_1year_herbicide_analysis.pdf", width = 11.5, height = 3.5)
ggplot(plant_ctrl_beta_pred1, aes(x = PropCoveredBeta, y = Next1PropCoveredBeta)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_ribbon(aes(y = pred, ymin = pred - pred_se, ymax = pred + pred_se, fill = TreatedF), 
              color = NA, alpha = 0.5) +
  geom_line(aes(y = pred, color = TreatedF)) +
  facet_wrap(~ CommonName, scales = "free") +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Initial population size [", N[t], "/A]", sep = "")),
       y = expression(paste("Future population size [", N[t+1], "/A]", sep = ""))) +
  scale_color_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  scale_fill_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  def_theme +
  theme(legend.position = c(0.77, 0.77),
        strip.text = element_blank())
dev.off()


#### visualize 2-year alt beta models ####

# add predicted values
plant_ctrl_beta_pred2 <- hyd_dat_2yearAlt %>%
  mutate(pred = predict(hyd_beta_mod_2yearAlt, newdata = ., re.form = NA, type = "response"),
         pred_se = predict(hyd_beta_mod_2yearAlt, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit) %>%
  full_join(why_dat_2yearAlt %>%
              mutate(pred = predict(why_beta_mod_2yearAlt, newdata = ., re.form = NA, type = "response"),
                     pred_se = predict(why_beta_mod_2yearAlt, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit)) %>%
  full_join(wle_dat_2yearAlt %>%
              mutate(pred = predict(wle_beta_mod_2yearAlt, newdata = ., re.form = NA, type = "response"),
                     pred_se = predict(wle_beta_mod_2yearAlt, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit))

# figure
pdf("output/beta_model_fit_2year_herbicide_analysis.pdf", width = 11.5, height = 4)
ggplot(plant_ctrl_beta_pred2, aes(x = PropCoveredBeta, y = Next2PropCoveredBeta)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.2, aes(color = TreatedF)) +
  geom_ribbon(aes(y = pred, ymin = pred - pred_se, ymax = pred + pred_se, fill = TreatedF), 
              color = NA, alpha = 0.5) +
  geom_line(aes(y = pred, color = TreatedF)) +
  facet_wrap(~ CommonName, scales = "free") +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Initial population size [", N[t], "/A]", sep = "")),
       y = expression(paste("Future population size [", N[t+2], "/A]", sep = ""))) +
  scale_color_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  scale_fill_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  def_theme +
  theme(legend.position = c(0.45, 0.77))
dev.off()

pdf("output/beta_model_only_2year_herbicide_analysis.pdf", width = 11.5, height = 3.5)
ggplot(plant_ctrl_beta_pred2, aes(x = PropCoveredBeta, y = Next2PropCoveredBeta)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_ribbon(aes(y = pred, ymin = pred - pred_se, ymax = pred + pred_se, fill = TreatedF), 
              color = NA, alpha = 0.5) +
  geom_line(aes(y = pred, color = TreatedF)) +
  facet_wrap(~ CommonName, scales = "free") +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("Initial population size [", N[t], "/A]", sep = "")),
       y = expression(paste("Future population size [", N[t+2], "/A]", sep = ""))) +
  scale_color_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  scale_fill_manual(values = c("black", "blue"), name = "Treated\nin initial\nyear") +
  def_theme +
  theme(legend.position = c(0.43, 0.77),
        strip.text = element_blank())
dev.off()
