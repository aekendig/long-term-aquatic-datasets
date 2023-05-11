#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)
library(lubridate)

# import data
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
# plant_fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv",
#                        col_types = list(Depth_ft = col_double(),
#                                         PermanentID = col_character(),
#                                         YearF = col_character()))
key_all_acre <- read_csv("original-data/FWC_plant_survey_key_all_acreage.csv")
key_all_pres <- read_csv("original-data/FWC_plant_survey_key_all_presence.csv")
inv_plant <- read_csv("intermediate-data/FWC_only_invasive_plant_formatted.csv")
ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv") 

#### edit plant community data ####

# one origin per species?
plant_fwc %>%
  group_by(TaxonName) %>%
  summarise(orig = length(unique(Origin))) %>%
  ungroup() %>%
  filter(orig != 1)
# none with multiple origins

# Eppc and origin
plant_fwc %>%
  select(Eppc, Origin) %>%
  unique()

# years when all taxa were surveyed
surv_years <- key_all_acre %>% full_join(key_all_pres)

# select for years with all species surveyed based on keys
# Tohopekaliga survey in 2017 seems incomplete in other data exploration
plant_fwc2 <- plant_fwc %>%
  inner_join(surv_years) %>%
  filter(!(AreaOfInterestID == 476 & SurveyYear == 2017)) %>%
  mutate(Detected = case_when(IsDetected == "Yes" ~ 1,
                              IsDetected == "No" ~ 0),
         PreCtrl = if_else(GSYear < 1998, "pre ctrl data", "post ctrl data"),
         Area_ha = ShapeArea * 100)

# richness over time
pdf("output/plant_richness_raw_time_series.pdf", width = 6, height = 3.5)
plant_fwc2 %>%
  filter(IsDetected == "Yes") %>%
  group_by(AreaOfInterestID, SurveyYear, Origin) %>%
  summarize(Richness = n_distinct(TaxonName)) %>%
  ungroup() %>%
  ggplot(aes(x = SurveyYear, y = Richness, color = as.factor(AreaOfInterestID))) +
  geom_line() +
  facet_wrap(~ Origin, scales = "free") +
  theme(legend.position = "none")

plant_fwc2 %>%
  filter(IsDetected == "Yes") %>%
  group_by(AreaOfInterestID, SurveyYear, Origin) %>%
  summarize(Richness = n_distinct(TaxonName)) %>%
  ungroup() %>%
  ggplot(aes(x = SurveyYear)) +
  geom_bar() +
  facet_wrap(~ Origin, scales = "free") +
  labs(x = "Year", y = "Number of surveys with richness > 0")
dev.off()

# duplicate surveys in a year
plant_fwc2 %>%
  group_by(AreaOfInterestID, GSYear) %>%
  summarise(surveys = length(unique(SurveyDate))) %>%
  ungroup() %>%
  filter(surveys > 1) 
# 33 AOIs have multiple surveys in a year

# summarize to remove duplicates
plant_fwc3 <- plant_fwc2 %>%
  group_by(AreaOfInterestID, AreaOfInterest, PermanentID, Area_ha, GSYear, PreCtrl, TaxonName, Habitat, Origin) %>%
  summarize(SurveyDate = max(SurveyDate),
            Detected = as.numeric(sum(Detected) > 0)) %>%
  ungroup()

plant_fwc3 %>%
  group_by(AreaOfInterestID, GSYear, TaxonName) %>%
  summarise(surveys = length(unique(SurveyDate))) %>%
  ungroup() %>%
  filter(surveys > 1) 
# none

# add row for every year for each site/species combo (NA's for missing surveys)
plant_fwc4 <- plant_fwc3 %>%
  group_by(AreaOfInterestID, TaxonName) %>%
  arrange(GSYear) %>% 
  mutate(GSYearDiff = GSYear - lag(GSYear),
         PrevDetected = if_else(GSYearDiff == 1, lag(Detected), NA_real_)) %>% # previous year's detection
  ungroup()

# species richness over time
pdf("output/plant_richness_processed_time_series.pdf", width = 6, height = 3.5)
plant_fwc4 %>%
  filter(Detected == 1) %>%
  group_by(PermanentID, GSYear, Origin) %>%
  summarize(Richness = n_distinct(TaxonName)) %>%
  ungroup() %>%
  ggplot(aes(x = GSYear, y = Richness, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ Origin, scales = "free") +
  theme(legend.position = "none")

plant_fwc4 %>%
  filter(Detected == 1) %>%
  group_by(PermanentID, GSYear, Origin) %>%
  summarize(Richness = n_distinct(TaxonName)) %>%
  ungroup() %>%
  ggplot(aes(x = GSYear)) +
  geom_bar() +
  facet_wrap(~ Origin, scales = "free") +
  labs(x = "Year", y = "Number of surveys with richness > 0")
dev.off()

# make sure missing Detected/PrevDetected applies to all taxa
plant_fwc4 %>%
  group_by(AreaOfInterestID, GSYear, Origin) %>%
  summarize(NADetected = sum(is.na(Detected)),
            TotDetected = length(Detected),
            NAPrevDetected = sum(is.na(PrevDetected)),
            TotPrevDetected = length(PrevDetected)) %>%
  ungroup() %>%
  filter((NADetected > 0 & NADetected != TotDetected) | 
           (NAPrevDetected > 0 & NAPrevDetected != TotPrevDetected))
# yes

# zero detected
plant_fwc4 %>%
  group_by(AreaOfInterestID, GSYear, SurveyDate) %>%
  summarize(TotDetected = sum(Detected)) %>%
  ungroup() %>%
  filter(TotDetected == 0)
# none

# save
write_csv(plant_fwc4, "intermediate-data/FWC_plant_community_formatted.csv")


#### most common native taxa ####

# remove non-native taxa
nat_fwc <- plant_fwc4 %>%
  filter(Origin == "Native" & !is.na(Detected))

# total waterbody-year combos
TotWatYear <- nat_fwc %>%
  select(PreCtrl, AreaOfInterestID, GSYear) %>%
  unique() %>%
  group_by(PreCtrl) %>%
  summarize(TotWatYear = n(),
            TotWaterbody = n_distinct(AreaOfInterestID)) %>%
  ungroup()

# total habitat types
TotHabitat <- nat_fwc %>%
  group_by(Habitat) %>%
  summarize(TotTaxa = n_distinct(TaxonName)) %>%
  ungroup() %>%
  mutate(Taxa40 = round(TotTaxa * 0.4))
# only 2 floating taxa

# summarize waterbody-year occurrences
common_fwc <- nat_fwc %>%
  filter(Detected == 1) %>%
  group_by(TaxonName, Habitat, Origin, PreCtrl) %>%
  summarize(OccWatYear = n(),
            OccWaterbody = n_distinct(AreaOfInterestID)) %>%
  ungroup() %>%
  group_by(PreCtrl) %>%
  mutate(RankWatYear = rank(-OccWatYear),
         RankWaterbody = rank(-OccWaterbody)) %>%
  ungroup() %>%
  group_by(PreCtrl, Habitat) %>%
  mutate(RankWatYearHab = rank(-OccWatYear),
         RankWaterbodyHab = rank(-OccWaterbody)) %>%
  ungroup() %>%
  left_join(TotWatYear) %>%
  mutate(RatioWatYear = OccWatYear / TotWatYear,
         RatioWaterbody = OccWaterbody / TotWaterbody,
         Over20Waterbody = if_else(RatioWaterbody > 0.2, 1, 0), # thresholds based on post ctrl data distributions (below)
         Over7WatYear = if_else(RatioWatYear > 0.07, 1, 0)) # these are already split by PreCtrl

# distribution of occurrences
ggplot(common_fwc, aes(x = RankWatYear, y = RatioWatYear)) +
  geom_hline(yintercept = 0.07, linetype = "dashed") +
  geom_line() +
  geom_point(aes(color = Over20Waterbody)) +
  facet_wrap(~ PreCtrl)

 ggplot(common_fwc, aes(x = RankWaterbody, y = RatioWaterbody)) +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  geom_line() +
  geom_point(aes(color = Over7WatYear)) +
  facet_wrap(~ PreCtrl)

ggplot(common_fwc, aes(x = RankWatYearHab, y = RatioWatYear)) +
  geom_hline(yintercept = 0.07, linetype = "dashed") +
  geom_line() +
  geom_point(aes(color = Over20Waterbody)) +
  facet_grid(Habitat ~ PreCtrl)

ggplot(common_fwc, aes(x = RankWaterbodyHab, y = RatioWaterbody)) +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  geom_line() +
  geom_point(aes(color = Over7WatYear)) +
  facet_grid(Habitat ~ PreCtrl)

# common species in post control data
# taxa must be present in over 7% of water-year occurrences 
# and over 20% of lakes (this requirement is met with above)
common_fwc2 <- common_fwc %>%
  filter(PreCtrl == "post ctrl data" & Over7WatYear == 1 & Over20Waterbody == 1) %>%
  select(TaxonName) %>%
  inner_join(common_fwc)

common_fwc2 %>%
  group_by(PreCtrl, Habitat) %>%
  summarize(Taxa = n_distinct(TaxonName),
            MinWatYear = min(OccWatYear),
            MinWaterbody = min(OccWaterbody)) %>%
  left_join(TotHabitat) %>%
  mutate(RatioTaxa = Taxa / TotTaxa)

n_distinct(common_fwc2$TaxonName) # 64 taxa total

ggplot(common_fwc2, aes(x = RankWatYearHab, y = RatioWatYear)) +
  geom_line() +
  geom_point(aes(color = Over20Waterbody)) +
  facet_grid(Habitat ~ PreCtrl)

# select taxa
nat_fwc2 <- nat_fwc %>%
  filter(TaxonName %in% common_fwc2$TaxonName)

# save
write_csv(nat_fwc2, "intermediate-data/FWC_common_native_plants_formatted.csv")


#### START HERE: invasion and management data ####

# asked candice about absence of control data
# may need to update control dataset

# combine
inv_ctrl <- inner_join(inv_plant, ctrl %>%
                         rename(TaxonName = Species)) %>%
  group_by(AreaOfInterestID, TaxonName) %>%
  mutate(InvPresent = if_else(sum(EstAreaCoveredRaw_ha > 0) > 0, 
                              "yes", "no"),
         TreatPresent = if_else(sum(Treated > 0) > 0,
                                "yes", "no")) %>%
  ungroup()

# values
inv_ctrl %>%
  distinct(AreaOfInterestID, TaxonName, InvPresent, TreatPresent) %>%
  count(TaxonName, InvPresent, TreatPresent)

# waterbodies that have the species present and been managed at least once
perm_plant_ctrl <- inner_join(perm_plant, perm_ctrl) %>%
  select(PermanentID, CommonName)


#### combine data and select waterbodies/years ####

# add invasive plant and control data
nat_fwc3 <- nat_fwc2 %>%
  inner_join(inv_plant2) %>% # select waterbodies and years in both datasets
  inner_join(ctrl2) %>%
  left_join(perm_plant) %>%
  mutate(RecentTreatment = replace_na(RecentTreatment, 0),
         Established = replace_na(Established, 0))

# summarize for richness and add invasive plant and control data
rich <- nat_fwc2 %>%
  group_by(PermanentID, GSYear, Area_ha) %>%
  summarize(Richness = sum(Detected)) %>%
  ungroup() %>%
  mutate(LogRichness = log(Richness + 1)) %>%
  inner_join(inv_plant2) %>% # select waterbodies and years in both datasets
  inner_join(ctrl2) %>%
  left_join(perm_plant) %>%
  mutate(RecentTreatment = replace_na(RecentTreatment, 0),
         Established = replace_na(Established, 0))

# select waterbodies that have had species and at least one year of management
nat_inv <- nat_fwc3 %>%
  inner_join(perm_plant_ctrl)

rich_inv <- rich %>%
  inner_join(perm_plant_ctrl)

# save data
write_csv(nat_fwc3, "intermediate-data/FWC_common_native_plants_invasive_species_data_formatted.csv")
write_csv(nat_inv, "intermediate-data/FWC_common_native_plants_invaded_data_formatted.csv")

write_csv(rich, "intermediate-data/FWC_common_native_richness_invasive_species_data_formatted.csv")
write_csv(rich_inv, "intermediate-data/FWC_common_native_richness_invaded_data_formatted.csv")


#### waterbody counts ####
nat_fwc3 %>%
  group_by(CommonName, Established) %>%
  summarize(waterbodies = n_distinct(PermanentID)) %>%
  ungroup() %>%
  mutate(Established = fct_recode(as.factor(Established),
                                  "Invaded" = "1",
                                  "Uninvaded" = "0"),
         CommonName = tolower(CommonName)) %>%
  pivot_wider(names_from = Established,
              values_from = waterbodies) %>%
  write_csv("intermediate-data/native_plants_invaded_uninvaded_counts.csv")


#### FWRI data ####

# species names
common_taxa <- common_fwc2 %>%
  select(TaxonName) %>%
  unique() %>%
  rename(ScientificName = TaxonName) %>%
  mutate(ScientificName = str_replace_all(ScientificName, "spp.", "sp."))
  
# subset for species
nat_fwri <- plant_fwri %>%
  inner_join(common_taxa)

# species matched
n_distinct(nat_fwri$ScientificName)

# unmatched species
unmatch_fwc <- common_taxa %>%
  anti_join(nat_fwri) %>%
  arrange(ScientificName)

unmatch_fwri <- plant_fwri %>%
  select(ScientificName) %>%
  unique() %>%
  anti_join(nat_fwri) %>%
  arrange(ScientificName)

# revise names to match
common_taxa2 <- common_taxa %>%
  mutate(ScientificName = fct_recode(ScientificName,
                                     "Cladium mariscus jamaicense" = "Cladium jamaicense",
                                     "Lemna sp." = "Lemna/Spirodela sp.",
                                     "Luziola fluitans (syn. Hydrochloa caroliniensis)" = "Luziola fluitans",
                                     "Nuphar luteum" = "Nuphar advena",
                                     "Paspalum fluitans" = "Paspalum repens",
                                     "Utricularia inflata or Utricularia floridana" = "Utricularia floridana"))
# Nuphar leteum is the old name of Nuphar advena
# Paspalum fluitans is a synonym of Paspalum repens
# both Utricularia are on the FWC list

# update species list
nat_fwri2 <- plant_fwri %>%
  inner_join(common_taxa2)

# species left out
common_taxa2 %>%
  anti_join(nat_fwri2) %>%
  arrange(ScientificName)

n_distinct(nat_fwri2$ScientificName) # 50 taxa

# save
write_csv(nat_fwri2, "intermediate-data/FWRI_common_native_plants_formatted.csv")
