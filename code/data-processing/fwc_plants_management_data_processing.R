#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)

# import data
plants <- read_csv("intermediate-data/FWC_plant_formatted.csv")
mgmt <- read_csv("intermediate-data/FWC_management_formatted.csv")
key_all_pres <- read_csv("original-data/FWC_plant_survey_key_all_presence.csv") # surveys of all species + acreage of three focal invasives

# function to divide data into categories
div_fun <- function(dat, orig_dat, colName){
  
  # get median value
  med <- orig_dat %>%
    filter(!!sym(colName) > 0) %>%
    pull(!!sym(colName)) %>%
    median()
  
  # increase to 1% if under
  if(med < 1){
    med <- 1
  }
  
  # divide data based on median
  dat_out <- dat %>%
    mutate(tempNew = case_when(!!sym(colName) == 0 ~"none",
                               !!sym(colName) > 0 & !!sym(colName) <= med ~ "low",
                               !!sym(colName) > med ~"high"))
  
  # return data with new column
  return(dat_out)
}


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
# select years where all species were surveyd
# select waterbodies with at least three years of data
plants3 <- plants2 %>%
  cross_join(mgmt %>%
               summarize(MinTreatmentYear = min(TreatmentYear),
                         MaxTreatmentYear = max(TreatmentYear))) %>%
  filter(SurveyYear >= MinTreatmentYear & SurveyYear <= MaxTreatmentYear) %>%
  select(-c(MinTreatmentYear, MaxTreatmentYear)) %>%
  inner_join(key_all_pres) %>%
  group_by(AreaOfInterestID) %>%
  mutate(SurveyYears = n_distinct(SurveyYear)) %>%
  ungroup() %>%
  filter(SurveyYears >= 3)

# waterbodies surveyed
# range of survey years
waterbodies <- plants3 %>%
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
# select years during survey span
mgmt2 <- mgmt %>%
  inner_join(waterbodies) %>%
  filter(TreatmentYear >= MinSurveyYear & TreatmentYear <= MaxSurveyYear) %>%
  mutate(Species = replace_na(Species, "unknown"))

# check for duplicates
mgmt2 %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County, WaterbodyAcres) %>%
  get_dupes(AreaOfInterestID, County, WaterbodyAcres)


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
  filter(IsDetected == "Yes" & !(TaxonName %in% c("Hydrilla verticillata", "Eichhornia crassipes", "Pistia stratiotes"))) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, FType, SurveyYear, SurveyDate) %>%
  summarize(NativeRichness = sum(Origin == "Native"),
            NonNativeRichness = sum(Origin == "Exotic"),
            .groups = "drop") %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, FType, SurveyYear) %>%
  summarize(NativeRichness = max(NativeRichness),
            NonNativeRichness = max(NonNativeRichness),
            .groups = "drop")

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

# divide into groups
inv_sum2 <- inv_sum %>%
  div_fun(inv_sum, "HydrPAC") %>%
  rename(HydrPACF = tempNew) %>%
  div_fun(inv_sum, "FloatPAC") %>%
  rename(FloatPACF = tempNew)

# distributions
ggplot(inv_sum2, aes(x = HydrPAC, fill = HydrPACF)) +
  geom_histogram(binwidth = 1)

ggplot(inv_sum2, aes(x = FloatPAC, fill = FloatPACF)) +
  geom_histogram(binwidth = 0.1)


#### summarize management data ####

# management targets
mgmt_targets <- mgmt2 %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County, 
           TreatmentYear, Species)

# rename targets
mgmt3 <- mgmt2 %>%
  mutate(Target = case_when(Species %in% c("Eichhornia crassipes", "Pistia stratiotes", 
                                           "Floating Plants (Eichhornia and Pistia)") ~ 
                              "Float",
                            Species == "Hydrilla verticillata" ~ "Hydr",
                            TRUE ~ "Other")) 

# summarize across groups
mgmt_targets_sum <- mgmt3 %>%
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
# need to clean up synonyms if using this

# summarize by waterbody and Target
mgmt_sum <- mgmt3 %>%
  mutate(SurveyYears = MaxSurveyYear - MinSurveyYear + 1,
         PropTreated = if_else(TotalAcres <= WaterbodyAcres, # Some exceed waterbody acres. Not sure why.
                               TotalAcres / WaterbodyAcres,
                               1)) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, Target) %>%
  summarize(TrtFreq = 100 * n_distinct(TreatmentYear) / unique(SurveyYears),
            TrtArea = mean(100 * PropTreated),
            .groups = "drop") %>%
  full_join(waterbodies %>% # add in waterbodies without management
              select(PermanentID, AreaOfInterestID, AreaOfInterest, County) %>%
              expand_grid(Target = c("Float", "Hydr", "Other"))) %>%
  mutate(TrtFreq = replace_na(TrtFreq, 0),
         TrtArea = replace_na(TrtArea, 0)) %>%
  pivot_wider(names_from = Target,
              values_from = c(TrtFreq, TrtArea),
              names_glue = "{Target}{.value}")

# divide into low, medium, high
mgmt_sum2 <- mgmt_sum %>%
  div_fun(mgmt_sum, "HydrTrtFreq") %>%
  rename(HydrTrtFreqF = tempNew) %>%
  div_fun(mgmt_sum, "HydrTrtArea") %>%
  rename(HydrTrtAreaF = tempNew) %>%
  div_fun(mgmt_sum, "FloatTrtFreq") %>%
  rename(FloatTrtFreqF = tempNew) %>%
  div_fun(mgmt_sum, "FloatTrtArea") %>%
  rename(FloatTrtAreaF = tempNew) %>%
  div_fun(mgmt_sum, "OtherTrtFreq") %>%
  rename(OtherTrtFreqF = tempNew) %>%
  div_fun(mgmt_sum, "OtherTrtArea") %>%
  rename(OtherTrtAreaF = tempNew) 

# distribution
ggplot(mgmt_sum2, aes(x = HydrTrtFreq, fill = HydrTrtFreqF)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_sum2, aes(x = HydrTrtArea, , fill = HydrTrtAreaF)) +
  geom_histogram(binwidth = 1)

ggplot(mgmt_sum2, aes(x = FloatTrtFreq, fill = FloatTrtFreqF)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_sum2, aes(x = FloatTrtArea, fill = FloatTrtAreaF)) +
  geom_histogram(binwidth = 1)

ggplot(mgmt_sum2, aes(x = OtherTrtFreq, fill = OtherTrtFreqF)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_sum2, aes(x = OtherTrtArea, fill = OtherTrtAreaF)) +
  geom_histogram(binwidth = 1)


#### management metrics for lakes with management ####

# average month
# maybe something related to method

# management methods
mgmt_methods_sum <- mgmt3 %>%
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

# make fewer groups
mgmt_methods <- mgmt3 %>%
  filter(CtrlSet == "new") %>%
  mutate(Method = case_when(MethodHerbicide == "yes" & Contact == 1 ~ "Con", # contact
                            MethodHerbicide == "yes" & Contact == 0 ~ "Sys", # systemtic
                            MethodHerbicide == "no" ~ "Non", # non-herbicide
                            TRUE ~ "Unk"), # unknown
         SurveyYears = MaxSurveyYear - MinSurveyYear + 1,
         PropTreated = if_else(TotalAcres <= WaterbodyAcres, # Some exceed waterbody acres. Not sure why.
                               TotalAcres / WaterbodyAcres,
                               1))

mgmt_methods_sum2 <- mgmt_methods %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County, 
           TreatmentYear, Method, Target) %>%
  count(TreatmentYear, Method, Target) %>%
  group_by(TreatmentYear) %>%
  mutate(Prop = n / sum(n)) %>%
  ungroup() %>%
  group_by(TreatmentYear, Target) %>%
  mutate(TargetProp = n / sum(n)) %>%
  ungroup()

ggplot(mgmt_methods_sum2, aes(x = as.character(TreatmentYear), 
                             y = Prop, fill = Method)) +
  geom_col() +
  theme_bw()

ggplot(mgmt_methods_sum2, aes(x = as.character(TreatmentYear), 
                              y = TargetProp, fill = Method)) +
  geom_col() +
  theme_bw() +
  facet_wrap(~ Target)

# have any lakes with hydrilla management been only non-herbicide?
mgmt_methods %>%
  filter(Target == "Hydr") %>%
  distinct(AreaOfInterestID, Method) %>%
  group_by(AreaOfInterestID) %>%
  mutate(methods = n_distinct(Method)) %>%
  ungroup() %>%
  filter(methods == 1 & Method == "Non")

filter(mgmt_methods, Target == "Hydr" & AreaOfInterestID == 290)
# just one, and it was only treated in one year in the new data

# summed area treated for each method
mgmt_methods_sum3 <- mgmt_methods %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, Method) %>%
  summarize(TrtFreq = 100 * n_distinct(TreatmentYear) / unique(SurveyYears),
            TrtArea = mean(100 * PropTreated),
            .groups = "drop") %>%
  pivot_wider(names_from = Method,
              values_from = c(TrtArea, TrtFreq),
              names_glue = "{.value}{Method}") %>%
  mutate(across(.cols = starts_with("Trt"), 
                .fns = ~replace_na(.x, 0))) 

# distributions
ggplot(mgmt_methods_sum3, aes(x = TrtAreaCon)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_methods_sum3, aes(x = TrtAreaSys)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_methods_sum3, aes(x = TrtAreaNon)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_methods_sum3, aes(x = TrtAreaUnk)) +
  geom_histogram(binwidth = 1) # don't use this one

ggplot(mgmt_methods_sum3, aes(x = TrtFreqCon)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_methods_sum3, aes(x = TrtFreqSys)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_methods_sum3, aes(x = TrtFreqNon)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_methods_sum3, aes(x = TrtFreqUnk)) +
  geom_histogram(binwidth = 1) # don't use this one

# management timing
mgmt_timing_sum <- mgmt_methods %>%
  mutate(TreatmentQuarter = case_when(TreatmentMonth %in% 1:3 ~ "Q1",
                                      TreatmentMonth %in% 4:6 ~ "Q2",
                                      TreatmentMonth %in% 7:9 ~ "Q3",
                                      TreatmentMonth %in% 10:12 ~ "Q4")) %>%
  group_by(PermanentID, AreaOfInterestID, AreaOfInterest, County, TreatmentQuarter) %>%
  summarize(TrtFreq = 100 * n_distinct(TreatmentYear) / unique(SurveyYears),
            TrtArea = mean(100 * PropTreated),
            .groups = "drop") %>%
  pivot_wider(names_from = TreatmentQuarter,
              values_from = c(TrtArea, TrtFreq),
              names_glue = "{.value}{TreatmentQuarter}") %>%
  mutate(across(.cols = starts_with("Trt"), 
                .fns = ~replace_na(.x, 0))) 

# distributions
ggplot(mgmt_timing_sum, aes(x = TrtAreaQ1)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_timing_sum, aes(x = TrtAreaQ2)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_timing_sum, aes(x = TrtAreaQ3)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_timing_sum, aes(x = TrtAreaQ4)) +
  geom_histogram(binwidth = 1)

ggplot(mgmt_timing_sum, aes(x = TrtFreqQ1)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_timing_sum, aes(x = TrtFreqQ2)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_timing_sum, aes(x = TrtFreqQ3)) +
  geom_histogram(binwidth = 1)
ggplot(mgmt_timing_sum, aes(x = TrtFreqQ4)) +
  geom_histogram(binwidth = 1) # don't use this one


#### combine data ####

# combine data
# create a time variable
rich_dat <- plant_sum %>%
  full_join(inv_sum2) %>%
  full_join(mgmt_sum2) %>%
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

# for individual species
taxa_dat <- plants3 %>%
  mutate(Detected = if_else(IsDetected == "Yes", 1, 0)) %>%
  full_join(inv_sum2) %>%
  full_join(mgmt_sum2)
  
# for specific management methods
methods_dat <- plant_sum %>%
  inner_join(mgmt_methods_sum3 %>% 
               full_join(mgmt_timing_sum))


#### save data ####

write_csv(rich_dat, "intermediate-data/FWC_plant_management_richness_analysis_formatted.csv")
write_csv(taxa_dat, "intermediate-data/FWC_plant_management_taxa_analysis_formatted.csv")
write_csv(methods_dat, "intermediate-data/FWC_plant_management_methods_analysis_formatted.csv")
