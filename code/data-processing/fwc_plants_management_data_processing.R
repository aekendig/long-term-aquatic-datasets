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
  inner_join(waterbodies) %>%
  mutate(TreatmentTarget = case_when(Species == "Hydrilla verticillata" ~ "hydrilla", 
                                     Species %in% c("Eichhornia crassipes", "Floating Plants (Eichhornia and Pistia)") ~ "hyacinth", 
                                     Species %in% c("Pistia stratiotes", "Floating Plants (Eichhornia and Pistia)") ~ "lettuce",
                                     TRUE ~ "other"))

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
# all the old treatments are piled on the same date within a year (1998-2010)
# some of the management surveys go beyond the plant surveys

# zoom in on new management
ggplot(filter(plant_surv, SurveyYear >= 2010), 
       aes(x = SurveyDate, y = fct_rev(AreaOfInterest))) +
  geom_point(size = 0.4) +
  geom_point(data = filter(mgmt2, CtrlSet == "new"), 
             aes(x = TreatmentDate, y = fct_rev(AreaOfInterest)),
                 size = 0.4, shape = 1, color = "pink")

# indicate maximum survey date
waterbodies2 <- waterbodies %>%
  left_join(plant_surv %>%
              group_by(AreaOfInterestID) %>%
              summarize(MaxSurveyDate = max(SurveyDate),
                        .groups = "drop")) %>%
  mutate(MaxSurveyYear = year(MaxSurveyDate))


#### summarize all management data ####

# whether or not management occurred each year
mgmt_year <- mgmt2 %>%
  left_join(waterbodies2) %>%
  filter(TreatmentDate <= MaxSurveyDate) %>% # remove dates older than the max survey date for each waterbody (old management, date = 1/1)
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID, TreatmentYear, TreatmentTarget) %>% # isolate waterbodies and years in management dataset
  mutate(Treatment = 1) %>%
  full_join(plant_surv %>%
              distinct(PermanentID, AreaOfInterest, AreaOfInterestID) %>% # add all waterbodies from the plant survey dataset
              expand_grid(TreatmentYear = min(mgmt2$TreatmentYear):max(mgmt2$TreatmentYear)) %>% # expanded by all years in the management data
              expand_grid(TreatmentTarget = c("hydrilla", "hyacinth", "lettuce", "other"))) %>% # and all targets
  mutate(Treatment = replace_na(Treatment, 0)) %>%
  left_join(waterbodies2 %>%
              select(PermanentID, AreaOfInterest, AreaOfInterestID, MaxSurveyYear)) %>%
  filter(TreatmentYear <= MaxSurveyYear) # remove years older than max survey year

# waterbodies missing
anti_join(waterbodies2, mgmt_year) 
# all of these had surveys that ended before management started

# update waterbodies list
waterbodies3 <- waterbodies2 %>%
  filter(AreaOfInterestID %in% mgmt_year$AreaOfInterestID)

# for each waterbody, number of management years by target
mgmt_year_sum <- mgmt_year %>%
  group_by(PermanentID, AreaOfInterest, AreaOfInterestID, TreatmentTarget) %>%
  summarize(TreatmentYears = sum(Treatment),
            TotalYears = n_distinct(TreatmentYear),
            .groups = "drop")


#### summarize new management data by survey-span ####

# this is for changes in response variables between surveys
# between each plant survey we need to know:
# number of days between surveys
# management type (herbicide, other categories)
# timing (month) - maybe just for herbicide
# proportion of waterbody managed - maybe just for herbicide

# calculate plant survey timing
plant_surv2 <- plant_surv %>%
  arrange(AreaOfInterestID, SurveyDate) %>%
  group_by(AreaOfInterestID) %>%
  mutate(LastSurveyDate = lag(SurveyDate)) %>%
  ungroup() %>%
  mutate(BetweenSurveyDays = SurveyDate - LastSurveyDate) %>%
  filter(SurveyYear >= min(mgmt2$TreatmentYear) & SurveyYear <= max(mgmt2$TreatmentYear)) # select waterbodies/years surveyed within management data

# NA's for surveys that didn't start until later
filter(plant_surv2, is.na(BetweenSurveyDays)) # 33

# distribution of survey spans
ggplot(plant_surv2, aes(x = BetweenSurveyDays)) +
  geom_density()
# some are super long -- will need to cut-off probably, but leave in all for now

# new management data
mgmt2_new <- mgmt2 %>%
  filter(CtrlSet == "new")

# select surveys we have new management data for
plant_surv2_new <- plant_surv2 %>%
  filter(LastSurveyDate >= min(mgmt2_new$TreatmentDate) &
           SurveyDate <= max(mgmt2_new$TreatmentDate))

# combine new management and surveys
mgmt_between <- left_join(mgmt2_new, plant_surv2_new, relationship = "many-to-many") %>%
  filter(TreatmentDate >= LastSurveyDate & TreatmentDate <= SurveyDate) %>% # select management between surveys
  full_join(plant_surv2_new %>% # add all surveys (management will be NA)
              expand_grid(TreatmentTarget = c("hydrilla", "hyacinth", "lettuce", "other"))) # expanded by all targets

# for each survey, what management methods were used?
mgmt_between_method <- mgmt_between %>%
  mutate(MethodHerbicide = replace_na(MethodHerbicide, "none") %>%
           fct_recode(yes = "unknown")) %>%
  count(PermanentID, AreaOfInterest, AreaOfInterestID, SurveyDate, LastSurveyDate, BetweenSurveyDays, 
           TreatmentTarget, MethodHerbicide) %>% # 7 cases, seems more likely they'd be herbicide
  pivot_wider(names_from = MethodHerbicide,
              values_from = n) %>%
  mutate(across(.cols = c(none, yes, no), .fns = ~replace_na(.x, 0)),
         none = if_else(none == 0, 1, 0)) %>%
  rename(Treatment = none,
         Herbicide = yes,
         NonHerbicide = no)

#### start here: herbicide timing and area ####


### older code - potentially useful, check before deleting ####


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