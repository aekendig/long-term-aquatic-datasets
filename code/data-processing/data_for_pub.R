# load libraries
library(tidyverse)

# import data
plants <- read_csv("intermediate-data/FWC_plant_formatted.csv")
mgmt <- read_csv("intermediate-data/FWC_management_formatted.csv")
loc_dat <- read_csv("intermediate-data/FWC_lakes_formatted.csv")
fwri <- read_csv("intermediate-data/fwri_double_sampling_data.csv")
fwc <- read_csv("intermediate-data/fwc_double_sampling_data.csv")

# select columns
plants2 <- plants %>% 
  filter(WaterbodyType == "Lake") %>% 
  select(AreaOfInterestID, AreaOfInterest, WaterbodyAcres,
         SurveyYear, SurveyDate, TaxonName, IsDetected, Origin, Habitat,
         SpeciesAcres, Incomplete = Outlier)

mgmt2 <- mgmt %>% 
  select(AreaOfInterestID, AreaOfInterest, TreatmentDate, TreatmentMonth,
         TreatmentYear, CtrlSet, Species, TaxonName, ControlMethod, 
         MethodHerbicide, Contact, TotalAcres) %>% 
  rename(TreatmentSet = CtrlSet)

fwc2 <- fwc %>% 
  rename_with(~str_remove(.x, "fwc_")) %>% 
  select(Year, Date, PermanentID, AreaOfInterest = AOIs, TaxonName, Habitat, 
         IsDetected, PAC, WaterbodyAcres)

fwri2 <- fwri %>% 
  rename_with(~str_remove(.x, "fwri_")) %>% 
  select(Year, Date, PermanentID, AreaOfInterest = AOIs, TaxonName,
         IsDetected, PAC1, PAC2, PAC3, TotalPoints)

# add permanent ID and shape source to location data
loc_dat2 <- loc_dat %>% 
  left_join(plants %>% 
              distinct(AreaOfInterestID, PermanentID, ShapeSource))

# check
filter(loc_dat2, is.na(PermanentID))
filter(loc_dat2, is.na(ShapeSource))

# export data
write_csv(plants2, "../aquatic-management-effects/data/plant_data.csv")
write_csv(mgmt2, "../aquatic-management-effects/data/management_data.csv")
write_csv(fwc2, "../aquatic-management-effects/data/meander_data.csv")
write_csv(fwri2, "../aquatic-management-effects/data/point_intercept_data.csv")
write_csv(loc_dat2, "../aquatic-management-effects/data/waterbody_data.csv")
