#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)

# figure settings
source("code/settings/figure_settings.R")

# import data
plants <- read_csv("intermediate-data/FWC_plant_formatted.csv")
mgmt <- read_csv("intermediate-data/FWC_management_formatted.csv")


#### format data ####

# types of waterbodies
unique(plants$WaterbodyType)

# look at rivers and canals
plants %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID, WaterbodyType) %>%
  filter(WaterbodyType %in% c("River", "Canal"))
# remove Boggy Creek, Josephine Creek, Silver Glen Springs, Swift Creek, Tarpon Lake Outfall Canal
# takes care of the lake-add-ons identified in fwc_plant_data_processing and double_sampling_data_processing
# check Mud Lake (categorized as River) - slightly wider part of a river

# features types of "Lake"
# https://www.usgs.gov/ngp-standards-and-specifications/national-hydrography-dataset-nhd-data-dictionary-feature-classes
plants %>%
  filter(WaterbodyType == "Lake") %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID, FType) %>%
  count(FType)
# 390: lake/pond
# 436: reservoir
# 466: swamp/marsh

# look more closely at swamp/marsh
plants%>%
  filter(FType == 466) %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID)
# Lafayette Lake looks like a wetland, but it borders an open waterbody names Lafayette
# open Lafayette permanent ID: 164324612
filter(plants, PermanentID == "164324612")
# used the wetland/swamp polygon because this is where the FWC herbicide coordinates landed
# the open water part looks like it's caused by a dam/structure
# okay to call this lake/reservoir because that what FWC calls it

# filter plant dataset for lakes
plants2 <- plants %>%
  filter(WaterbodyType == "Lake")

# waterbodies that were surveyed
waterbodies <- plants2 %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID)

# can I use AOI ID?
n_distinct(waterbodies$AreaOfInterestID)
n_distinct(waterbodies$AreaOfInterestID) == nrow(waterbodies)
# yes

# plant surveys
plant_surv <- plants2 %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID, SurveyYear, SurveyDate)

# select waterbodies from management data with surveys
mgmt2 <- mgmt %>%
  inner_join(waterbodies)

# were any waterbodies not treated?
(plant_no_mgmt <- plant_surv %>%
  anti_join(mgmt %>%
              distinct(PermanentID, AreaOfInterest, AreaOfInterestID)) %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID) %>%
    arrange(AreaOfInterest))
# 55

# any waterbodies treated, but not surveyed?
(mgmt_no_plant <- mgmt %>%
    distinct(PermanentID, AreaOfInterest, AreaOfInterestID) %>%
    anti_join(plant_surv %>%
                distinct(PermanentID, AreaOfInterest, AreaOfInterestID)) %>%
    arrange(AreaOfInterest))
# 37

# do any have same permanent ID? (i.e., were actually treated)
plant_no_mgmt %>%
  distinct(PermanentID) %>%
  inner_join(mgmt_no_plant %>%
               distinct(PermanentID))
# no, likely all different waterbodies
# manually checked names for mispellings/errors

# plot to explore
ggplot(plant_surv, aes(x = SurveyDate, y = fct_rev(AreaOfInterest))) +
  geom_hline(aes(yintercept = fct_rev(AreaOfInterest)), linewidth = 0.2) +
  geom_point(size = 0.4) +
  geom_point(data = mgmt2, aes(x = TreatmentDate, y = fct_rev(AreaOfInterest),
                               color = CtrlSet), size = 0.4, shape = 1)
# no mgmt data before 1998
# all the old treatments are piled on the same date (1998-2010)
# some of the management surveys go beyond the plant surveys


#### summarize whether or not waterbody was treated each year ####

# this is for annual changes in response variables
# between each plant survey we need to know:
# number of days between surveys
# whether it was treated

# plant survey timing
plant_surv2 <- plant_surv %>%
  arrange(AreaOfInterestID, SurveyDate) %>%
  group_by(AreaOfInterestID) %>%
  mutate(LastSurveyDate = lag(SurveyDate)) %>%
  ungroup() %>%
  mutate(BetweenSurveyDays = SurveyDate - LastSurveyDate)

# there should be 391 NA's
filter(plant_surv2, is.na(BetweenSurveyDays)) # yes
ggplot(plant_surv2, aes(x = BetweenSurveyDays)) +
  geom_density()
# some are super long -- will need to cut-off probably, but leave in all for now

#### start here ####
# add management dates to survey dates
# how to deal with July 1st control dates?


#### summarize new treatments per waterbody per year ####

# this is for annual changes in response variables
# between each plant survey we need to know:
# number of days between surveys
# management type (herbicide, other categories)
# timing (month) - maybe just for herbicide
# proportion of waterbody managed - maybe just for herbicide



### older code - potentially useful, check before deleting ####

# identify max survey date
plant_surv_max <- plant_surv %>%
  group_by(PermanentID, AreaOfInterest, AreaOfInterestID) %>%
  summarize(MaxSurveyDate = max(SurveyDate),
            .groups = "drop")

# duplicate lake names? need to specify county
mgmt_dup_lakes <- mgmt2 %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID) %>%
  get_dupes(AreaOfInterest)

# restrict mgmt to before last survey
# identify treatment targets
mgmt3 <- mgmt2 %>%
  left_join(plant_surv_max) %>%
  filter(TreatmentDate <= MaxSurveyDate) %>%
  mutate(TreatmentTarget = case_when(Species == "Hydrilla verticillata" ~ "hydrilla", 
                                     Species == "Eichhornia crassipes" ~ "hyacinth", 
                                     Species == "Pistia stratiotes" ~ "lettuce",
                                     Species == "Floating Plants (Eichhornia and Pistia)" ~ "floating",
                                     TRUE ~ "other"),
         Waterbody = if_else(AreaOfInterest %in% mgmt_dup_lakes$AreaOfInterest,
                             paste0(AreaOfInterest, " (", str_to_title(County_FWC), ")"),
                             AreaOfInterest))

# for the new dataset, how many treatments does each waterbody have until the
# last survey, and how do these breakdown by target?
mgmt_new_trts_sum <- mgmt3 %>%
  filter(CtrlSet == "new") %>%
  group_by(PermanentID, AreaOfInterest, AreaOfInterestID, Waterbody, TreatmentTarget) %>%
  summarize(Treatments = n_distinct(TreatmentID),
            .groups = "drop") %>%
  full_join(mgmt3 %>% # combine all floating
              filter(CtrlSet == "new" & TreatmentTarget %in% c("hyacinth", "lettuce", "floating")) %>%
              mutate(TreatmentTarget = "floating_combined") %>%
              group_by(PermanentID, AreaOfInterest, AreaOfInterestID, Waterbody, TreatmentTarget) %>%
              summarize(Treatments = n_distinct(TreatmentID),
                        .groups = "drop")) %>%
  full_join(mgmt3 %>% # combine all treatments
              filter(CtrlSet == "new") %>%
              mutate(TreatmentTarget = "total") %>%
              group_by(PermanentID, AreaOfInterest, AreaOfInterestID, Waterbody, TreatmentTarget) %>%
              summarize(Treatments = n_distinct(TreatmentID),
                        .groups = "drop")) %>%
  full_join(mgmt3 %>% # all combinations to add zeros
              filter(CtrlSet == "new") %>%
              distinct(PermanentID, AreaOfInterest, AreaOfInterestID, Waterbody) %>%
              expand_grid(TreatmentTarget = c("hydrilla", "hyacinth", "lettuce", 
                                              "floating", "other", "floating_combined", "total"))) %>%
  mutate(Treatments = replace_na(Treatments, 0))

# how many treatment years per waterbody by target?
mgmt_yrs_sum <- mgmt3 %>%
  group_by(PermanentID, AreaOfInterest, AreaOfInterestID, Waterbody, TreatmentTarget) %>%
  summarize(TreatmentYears = n_distinct(TreatmentYear),
            .groups = "drop") %>%
  full_join(mgmt3 %>% # combine all floating
              filter(TreatmentTarget %in% c("hyacinth", "lettuce", "floating")) %>%
              mutate(TreatmentTarget = "floating_combined") %>%
              group_by(PermanentID, AreaOfInterest, AreaOfInterestID, Waterbody, TreatmentTarget) %>%
              summarize(TreatmentYears = n_distinct(TreatmentYear),
                        .groups = "drop")) %>%
  full_join(mgmt3 %>% # combine all treatments
              mutate(TreatmentTarget = "total") %>%
              group_by(PermanentID, AreaOfInterest, AreaOfInterestID, Waterbody, TreatmentTarget) %>%
              summarize(TreatmentYears = n_distinct(TreatmentYear),
                        .groups = "drop")) %>%
  full_join(mgmt3 %>% # all combinations to add zeros
              distinct(PermanentID, AreaOfInterest, AreaOfInterestID, Waterbody) %>%
              expand_grid(TreatmentTarget = c("hydrilla", "hyacinth", "lettuce", 
                                              "floating", "other", "floating_combined", "total"))) %>%
  mutate(TreatmentYears = replace_na(TreatmentYears, 0))

# visualize
ggplot(mgmt_new_trts_sum, aes(x = Treatments)) +
  geom_density() +
  facet_wrap(~ TreatmentTarget, scales = "free") +
  def_theme

ggplot(mgmt_yrs_sum, aes(x = TreatmentYears)) +
  geom_density() +
  facet_wrap(~ TreatmentTarget, scales = "free") +
  def_theme

# look at some of the high values
(temp_high_floating <- filter(mgmt_new_trts_sum, 
                              TreatmentTarget == "floating" & Treatments > 2000) %>%
    inner_join(mgmt3))
# multiple treatments within the same day using different amounts of total herbicide or total acres treated,
# probably treating different areas of the lake
# need a more specific definition of a treatment than a unique treatment ID