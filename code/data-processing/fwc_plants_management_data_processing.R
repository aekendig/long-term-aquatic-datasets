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
key_all_pres <- read_csv("original-data/FWC_plant_survey_key_all_presence.csv")
key_exotic_acre <- read_csv("original-data/FWC_plant_survey_key_exotic_acreage.csv")


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
# remove problematic survey
plants2 <- plants %>%
  filter(WaterbodyType == "Lake" & Outlier == 0)

# waterbodies that were surveyed
waterbodies <- plants2 %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID, County)

# plant surveys
plant_surv <- plants2 %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID, County, WaterbodyAcres, 
           SurveyYear, SurveyDate)

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

# select plant survey data within management timeframe
plants3 <- plants2 %>%
  cross_join(mgmt %>%
               summarize(MinTreatmentYear = min(TreatmentYear),
                         MaxTreatmentYear = max(TreatmentYear))) %>%
  filter(SurveyYear >= MinTreatmentYear & SurveyYear <= MaxTreatmentYear) %>%
  select(-c(MinTreatmentYear, MaxTreatmentYear))

# waterbodies now included
# range of survey years
waterbodies2 <- plants3 %>%
  group_by(PermanentID, AreaOfInterest, AreaOfInterestID, County, 
           FType, WaterbodyAcres) %>%
  summarize(MinSurveyYear = min(SurveyYear),
            MaxSurveyYear = max(SurveyYear),
            .groups = "drop")

# mgmt record with missing info
filter(mgmt, is.na(Species)) %>%
  data.frame()
# 10 records have acres, but no other info

# select waterbodies from management data with surveys
mgmt2 <- mgmt %>%
  inner_join(waterbodies2) %>%
  filter(TreatmentYear >= MinSurveyYear & TreatmentYear <= MaxSurveyYear) %>%
  mutate(Species = replace_na(Species, "unknown"))

# check for duplicates
mgmt2 %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County, WaterbodyAcres) %>%
  get_dupes(AreaOfInterestID, County, WaterbodyAcres)

# final waterbodies not treated
waterbodies2 %>%
  anti_join(mgmt2 %>%
              distinct(PermanentID, AreaOfInterest, AreaOfInterestID)) %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID)


#### summarize plant data ####

# native vs. non-native
plants3 %>%
  distinct(TaxonName, Origin) %>%
  count(Origin)

# multiple surveys in a year
plants3 %>%
  distinct(AreaOfInterestID, SurveyYear, SurveyDate) %>%
  get_dupes(AreaOfInterestID, SurveyYear)

# richness
plant_sum <- plants3 %>%
  inner_join(key_all_pres) %>%
  filter(IsDetected == "Yes" & !(TaxonName %in% c("Hydrilla verticillata", "Eichhornia crassipes", "Pistia stratiotes"))) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, FType, SurveyYear, SurveyDate) %>%
  summarize(NativeRichness = sum(Origin == "Native"),
            NonNativeRichness = sum(Origin == "Exotic"),
            .groups = "drop") %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, FType, SurveyYear) %>%
  summarize(NativeRichness = max(NativeRichness),
            NonNativeRichness = max(NonNativeRichness),
            .groups = "drop") %>%
  group_by(AreaOfInterestID) %>%
  mutate(SurveyYears = n_distinct(SurveyYear)) %>%
  ungroup() %>%
  filter(SurveyYears >= 3)

# visualize
ggplot(plant_sum, aes(x = SurveyYear, y = NativeRichness, 
                      group = as.character(AreaOfInterestID), color = as.character(AreaOfInterestID))) +
  geom_smooth(formula = y~x, method = "lm", se = F, linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none")

ggplot(plant_sum, aes(x = SurveyYear, y = NonNativeRichness, 
                      group = as.character(AreaOfInterestID), color = as.character(AreaOfInterestID))) +
  geom_smooth(formula = y~x, method = "lm", se = F, linewidth = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none")

# mean invasive PAC
inv_sum <- plants3 %>%
  inner_join(key_exotic_acre) %>%
  filter(TaxonName %in% c("Hydrilla verticillata", "Eichhornia crassipes", "Pistia stratiotes")) %>%
  mutate(Invasive = case_when(TaxonName == "Hydrilla verticillata" ~ "Hydr", 
                                     TaxonName == "Eichhornia crassipes" ~ "Hyac", 
                                     TaxonName == "Pistia stratiotes" ~ "Lett")) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, Invasive,
           WaterbodyAcres, SurveyYear) %>%
  summarize(SpeciesAcres = max(SpeciesAcres), # for multiple surveys within a year
            .groups = "drop") %>%
  pivot_wider(names_from = Invasive,
              values_from = SpeciesAcres,
              names_glue = "{Invasive}Acres") %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County) %>%
  summarize(HydrPAC = mean(100 * HydrAcres / WaterbodyAcres),
            FloatPAC = mean(100 * (HyacAcres + LettAcres) / WaterbodyAcres), # these two are correlated
            .groups = "drop")

# distributions
ggplot(inv_sum, aes(x = HydrPAC)) +
  geom_histogram(binwidth = 1)

ggplot(inv_sum, aes(x = FloatPAC)) +
  geom_histogram(binwidth = 0.1)


#### summarize management data ####

# management targets
mgmt_targets <- mgmt2 %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County, 
           TreatmentYear, Species)

# most managed
mgmt_targets %>%
  count(Species) %>%
  arrange(desc(n))
# floating plants and hydrilla

# summarize across other species
mgmt_targets_sum <- mgmt2 %>%
  mutate(Target = case_when(Species %in% c("Eichhornia crassipes", "Pistia stratiotes", 
                                           "Floating Plants (Eichhornia and Pistia)") ~ 
                              "water hyacinth and/or water lettuce",
                            Species == "Hydrilla verticillata" ~ "hydrilla",
                            TRUE ~ "other") %>%
           fct_relevel("hydrilla", "water hyacinth and/or water lettuce")) %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County, 
           TreatmentYear, Target) %>%
  count(TreatmentYear, Target)

# visualize
ggplot(mgmt_targets_sum, aes(x = TreatmentYear, y = n, color = Target)) +
  geom_line() + 
  geom_point() +
  theme_bw()

# expand on other species
mgmt_targets_other <- mgmt_targets %>%
  filter(!(Species %in% c("Eichhornia crassipes", "Pistia stratiotes", 
                          "Floating Plants (Eichhornia and Pistia)",
                          "Hydrilla verticillata"))) %>%
  count(Species) %>%
  mutate(Species = fct_reorder(Species, -n))

# visualize
ggplot(mgmt_targets_other, aes(x = Species, y = n)) +
  geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# need to clean up synonyms

# summarize by waterbody and Target
mgmt_sum <- mgmt2 %>%
  mutate(SurveyYears = MaxSurveyYear - MinSurveyYear + 1,
         Target = case_when(Species %in% c("Eichhornia crassipes", "Pistia stratiotes", 
                                           "Floating Plants (Eichhornia and Pistia)") ~ 
                              "Float",
                            Species == "Hydrilla verticillata" ~ "Hydr",
                            TRUE ~ "Other"),
         PropTreated = if_else(TotalAcres <= WaterbodyAcres, # Some exceed waterbody acres. Not sure why.
                               TotalAcres / WaterbodyAcres,
                               1)) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, Target) %>%
  summarize(TrtFreq = 100 * n_distinct(TreatmentYear) / unique(SurveyYears),
            TrtArea = mean(100 * PropTreated),
            .groups = "drop") %>%
  full_join(waterbodies2 %>% # add in waterbodies without management
              select(PermanentID, AreaOfInterestID, AreaOfInterest, County) %>%
              expand_grid(Target = c("Float", "Hydr", "Other"))) %>%
  mutate(TrtFreq = replace_na(TrtFreq, 0),
         TrtArea = replace_na(TrtArea, 0)) %>%
  pivot_wider(names_from = Target,
              values_from = c(TrtFreq, TrtArea),
              names_glue = "{Target}{.value}")

# distribution
ggplot(mgmt_sum, aes(x = HydrTrtFreq)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_sum, aes(x = HydrTrtArea)) +
  geom_histogram(binwidth = 1)

ggplot(mgmt_sum, aes(x = FloatTrtFreq)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_sum, aes(x = FloatTrtArea)) +
  geom_histogram(binwidth = 1)

ggplot(mgmt_sum, aes(x = OtherTrtFreq)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_sum, aes(x = OtherTrtArea)) +
  geom_histogram(binwidth = 1)

# management methods
mgmt_methods_sum <- mgmt2 %>%
  filter(CtrlSet == "new") %>%
  mutate(Method = case_when(!is.na(MechanismOfAction) ~ str_to_lower(MechanismOfAction),
                            !is.na(ControlMethod) ~ str_to_lower(ControlMethod), # format these names
                            TRUE ~ "unknown")) %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County, 
           TreatmentYear, Method) %>%
  count(TreatmentYear, Method) %>%
  group_by(TreatmentYear) %>%
  mutate(Prop = n / sum(n)) %>%
  ungroup()

# visualize
ggplot(mgmt_methods_sum, aes(x = as.character(TreatmentYear), 
                             y = Prop, fill = Method)) +
  geom_col() +
  theme_bw()
# add names as text on last bar for clarity


#### to do: management metrics for lakes with management ####

# average month
# maybe something related to method


#### combine data ####

# combine data
# create a time variable
rich_dat <- plant_sum %>%
  inner_join(inv_sum) %>%
  inner_join(mgmt_sum) %>%
  mutate(Time = SurveyYear - min(SurveyYear))

# any missing rich_data?
filter(rich_dat, is.na(NativeRichness))
filter(rich_dat, is.na(NonNativeRichness))
filter(rich_dat, is.na(HydrPAC))
filter(rich_dat, is.na(FloatPAC))
filter(rich_dat, is.na(HydrTrtFreq))
filter(rich_dat, is.na(HydrTrtArea))
filter(rich_dat, is.na(FloatTrtFreq))
filter(rich_dat, is.na(FloatTrtArea))
filter(rich_dat, is.na(OtherTrtFreq))
filter(rich_dat, is.na(OtherTrtArea))

# total waterbodies
AOIs <- sort(unique(rich_dat$AreaOfInterestID))

# types of waterbodies
rich_dat %>%
  distinct(AreaOfInterestID, FType) %>%
  count(FType)

# waterbodies without management
rich_dat %>%
  filter(HydrTrtFreq == 0 & FloatTrtFreq == 0 & OtherTrtFreq == 0) %>%
  pull(AreaOfInterestID) %>%
  n_distinct()

# plot richness
pdf("output/analysis_data_richness_time_series.pdf")
for(i in AOIs){
  
  # filter rich_data
  rich_dat_sub <- filter(rich_dat, AreaOfInterestID == i) %>%
    pivot_longer(cols = c(NativeRichness, NonNativeRichness),
                 names_to = "Origin",
                 values_to = "Richness") %>%
    mutate(Origin = if_else(Origin == "NativeRichness", "native", "non-native"))
  
  # get name
  rich_dat_name <- unique(rich_dat_sub$AreaOfInterest)
  
  # figure
  print(ggplot(rich_dat_sub, aes(x = SurveyYear, y = Richness, color = Origin)) +
          geom_point() + 
          geom_line() +
          ggtitle(rich_dat_name) +
          theme_bw())
}
dev.off()



#### save data ####

write_csv(rich_dat, "intermediate-data/FWC_plant_management_richness_analysis_formatted.csv")


#### older code ####

# duplicate rows with floating plants
# add treatment targets
mgmt3 <- mgmt2 %>%
  mutate(Species = if_else(Species == "Floating Plants (Eichhornia and Pistia)", 
                           "Floating (Eichhornia)",
                           Species)) %>%
  full_join(mgmt2 %>%
              mutate(Species = if_else(Species == "Floating Plants (Eichhornia and Pistia)", 
                                       "Floating (Pistia)",
                                       Species))) %>%
  mutate(TreatmentTarget = case_when(Species == "Hydrilla verticillata" ~ "hydrilla", 
                                     Species %in% c("Eichhornia crassipes", "Floating (Eichhornia)") ~ "hyacinth", 
                                     Species %in% c("Pistia stratiotes", "Floating (Pistia)") ~ "lettuce",
                                     TRUE ~ "other"))

# check that the row duplication worked correctly
filter(mgmt2, Species == "Floating Plants (Eichhornia and Pistia)") %>% nrow() == nrow(mgmt3) - nrow(mgmt2)



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



#### summarize all management data ####

# whether or not management occurred each year
mgmt_year <- mgmt3 %>%
  left_join(waterbodies2) %>%
  filter(TreatmentDate <= MaxSurveyDate) %>% # remove dates older than the max survey date for each waterbody (old management, date = 1/1)
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID, TreatmentYear, TreatmentTarget) %>% # isolate waterbodies and years in management dataset
  mutate(Treatment = 1) %>%
  full_join(plant_surv %>%
              distinct(PermanentID, AreaOfInterest, AreaOfInterestID, County) %>% # add all waterbodies from the plant survey dataset
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

# duplicate lake names? need to specify county
wb3_dup_lakes <- waterbodies3 %>%
  distinct(PermanentID, AreaOfInterest, AreaOfInterestID, County) %>%
  get_dupes(AreaOfInterest)

# for each waterbody, number of management years by target
mgmt_year_sum <- mgmt_year %>%
  group_by(PermanentID, AreaOfInterest, AreaOfInterestID, County, TreatmentTarget) %>%
  summarize(TreatmentYears = sum(Treatment),
            TotalTreatmentYears = n_distinct(TreatmentYear),
            .groups = "drop") %>%
  full_join(mgmt_year %>%
              filter(TreatmentTarget %in% c("hyacinth", "lettuce")) %>%
              group_by(PermanentID, AreaOfInterest, AreaOfInterestID, County, TreatmentYear) %>% # summarize by waterbody and year
              summarize(Treatment = if_else(sum(Treatment) > 0, 1, 0), # if either species was treated, count
                        .groups = "drop") %>%
              group_by(PermanentID, AreaOfInterest, AreaOfInterestID, County) %>% # summarize by waterbody
              summarize(TreatmentYears = sum(Treatment),
                        TotalTreatmentYears = n_distinct(TreatmentYear),
                        .groups = "drop") %>%
              mutate(TreatmentTarget = "floating")) %>%
  mutate(WaterbodyName = if_else(AreaOfInterest %in% wb3_dup_lakes$AreaOfInterest,
                                 paste0(AreaOfInterest, " (", str_to_sentence(County), ")"),
                                 AreaOfInterest),
         PropTreatmentYears = TreatmentYears / TotalTreatmentYears)

# waterbodies
n_distinct(mgmt_year_sum$AreaOfInterestID) # 386
# targets
n_distinct(mgmt_year_sum$TreatmentTarget) # 5
# check rows
nrow(mgmt_year_sum) == n_distinct(mgmt_year_sum$AreaOfInterestID) * n_distinct(mgmt_year_sum$TreatmentTarget)

# distribution of values
ggplot(mgmt_year_sum, aes(x = PropTreatmentYears)) +
  geom_histogram(binwidth = 0.1)


#### summarize all invasive plant data ####

# select lakes and years in management data time span
plants3 <- plants2 %>%
  filter(AreaOfInterestID %in% waterbodies3$AreaOfInterestID &
           SurveyYear >= min(mgmt_year$TreatmentYear) &
           SurveyYear <= max(mgmt_year$TreatmentYear))

# select target species
# calculate proportion waterbody covered
inv_plants <- plants3 %>%
  filter(TaxonName %in% c("Hydrilla verticillata", "Eichhornia crassipes", "Pistia stratiotes")) %>%
  mutate(PropSpeciesAcres = SpeciesAcres / WaterbodyAcres,
         TreatmentTarget = case_when(TaxonName == "Hydrilla verticillata" ~ "hydrilla", 
                                     TaxonName == "Eichhornia crassipes" ~ "hyacinth", 
                                     TaxonName == "Pistia stratiotes" ~ "lettuce"))

# waterbodies
n_distinct(inv_plants$AreaOfInterestID) # 386
# species
n_distinct(inv_plants$TaxonName) # 3

# combine with management
inv_plants_mgmt <- left_join(inv_plants, mgmt_year_sum) %>%
  mutate(PropTreatmentYearsF = as.factor(round_half_up(PropTreatmentYears, 1)))

# visualize
ggplot(inv_plants_mgmt, aes(x = SurveyYear, y = PropSpeciesAcres, group = AreaOfInterestID)) +
  geom_smooth(method = "lm", se = F, linewidth = 0.5) +
  facet_grid(TreatmentTarget ~ PropTreatmentYearsF, scales = "free")

# save
write_csv(inv_plants_mgmt, "intermediate-data/FWC_invasive_plant_management_long_term.csv")


#### summarize all native plant data ####

# select native plants
native_plants <- plants3 %>%
  filter(Origin == "Native" & str_detect(TaxonName, "spp.") == F &
           TaxonName != "Filamentous algae")

# see if any are very rare
native_plants_freq <- native_plants %>%
  filter(IsDetected == "Yes") %>%
  group_by(TaxonName) %>%
  summarize(Waterbodies = n_distinct(AreaOfInterestID),
            Years = n_distinct(SurveyYear),
            Occurrences = n_distinct(paste(AreaOfInterestID, SurveyYear, sep = "_")),
            .groups = "drop")

ggplot(native_plants_freq, aes(x = Waterbodies)) +
  geom_histogram(binwidth = 1)
min(native_plants_freq$Waterbodies)

ggplot(native_plants_freq, aes(x = Years)) +
  geom_histogram(binwidth = 1)
min(native_plants_freq$Years)

ggplot(native_plants_freq, aes(x = Occurrences)) +
  geom_histogram(binwidth = 10)
min(native_plants_freq$Occurrences)

ggplot(native_plants_freq, aes(x = Waterbodies, y = Occurrences)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 5) # number of years present
# there are a handful with waterbodies < 50
# look at these more closely



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
              expand_grid(TreatmentTarget = c("hydrilla", "hyacinth", "lettuce", "other"))) %>% # expanded by all targets
  mutate(WaterbodyName = if_else(AreaOfInterest %in% wb3_dup_lakes$AreaOfInterest,
                                 paste0(AreaOfInterest, " (", str_to_sentence(County), ")"),
                                 AreaOfInterest))

# for each survey, what management methods were used?
mgmt_between_method <- mgmt_between %>%
  mutate(MethodHerbicide = replace_na(MethodHerbicide, "none") %>%
           fct_recode(yes = "unknown")) %>% # 7 cases, seems more likely they'd be herbicide
  count(PermanentID, AreaOfInterest, AreaOfInterestID, County, WaterbodyName, SurveyDate, LastSurveyDate, 
        BetweenSurveyDays, TreatmentTarget, MethodHerbicide) %>%
  pivot_wider(names_from = MethodHerbicide,
              values_from = n) %>%
  mutate(across(.cols = c(none, yes, no), .fns = ~replace_na(.x, 0)),
         none = if_else(none == 0, 1, 0)) %>%
  rename(Treatment = none,
         Herbicide = yes,
         NonHerbicide = no)

# for each survey with herbicide used in between, when was it applied?
mgmt_between_time <- mgmt_between %>%
  filter(MethodHerbicide == "yes") %>%
  count(PermanentID, AreaOfInterest, AreaOfInterestID, County, WaterbodyName, SurveyDate, LastSurveyDate, 
        BetweenSurveyDays, TreatmentTarget, TreatmentMonth) %>% # number of treatments per month
  mutate(TreatmentQuarter = case_when(TreatmentMonth %in% 1:3 ~ "Q1",
                                      TreatmentMonth %in% 4:6 ~ "Q2",
                                      TreatmentMonth %in% 7:9 ~ "Q3",
                                      TreatmentMonth %in% 10:12 ~ "Q4")) %>%
  group_by(PermanentID, AreaOfInterest, AreaOfInterestID, County, SurveyDate, LastSurveyDate, 
           BetweenSurveyDays, TreatmentTarget, TreatmentQuarter) %>%
  summarize(n = sum(n), # number of treatments per quarter
            .groups = "drop") %>%
  pivot_wider(names_from = TreatmentQuarter,
              values_from = n) %>%
  mutate(across(.cols = starts_with("Q"), .fns = ~ replace_na(.x, 0)))

ggplot(mgmt_between_time) +
  geom_density(aes(x = log(Q1 + 1))) +
  geom_density(aes(x = log(Q2 + 1)), color = "blue") +
  geom_density(aes(x = log(Q3 + 1)), color = "green") +
  geom_density(aes(x = log(Q4 + 1)), color = "yellow") +
  facet_wrap(~ TreatmentTarget)

# for each survey with herbicide used in between how much was treated?
mgmt_between_space <- mgmt_between %>%
  filter(MethodHerbicide == "yes") %>%
  mutate(PropTreated = TotalAcres / WaterbodyAcres) %>%
  group_by(PermanentID, AreaOfInterest, AreaOfInterestID, County, WaterbodyName, SurveyDate, LastSurveyDate, 
           BetweenSurveyDays, TreatmentTarget) %>%
  summarize(PropTreated = mean(PropTreated),
            .groups = "drop")

ggplot(mgmt_between_space, aes(x = log(PropTreated))) +
  geom_density() +
  facet_wrap(~ TreatmentTarget)

# combine space and time data
mgmt_between_herb <- full_join(mgmt_between_space, mgmt_between_time) %>%
  mutate(LogPropTreated = log(PropTreated),
         LogQ1 = log(Q1 + 1),
         LogQ2 = log(Q2 + 1),
         LogQ3 = log(Q3 + 1),
         LogQ4 = log(Q4 + 1))

#### start here: summarize plant responses ####





