#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)
library(lubridate)

# import data
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted_temporal_coverage.csv",
                      col_types = list(JoinNotes = col_character(),
                                       PermanentID = col_character()))
plant_fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv",
                       col_types = list(Depth_ft = col_double(),
                                        PermanentID = col_character(),
                                        YearF = col_character()))
key_all_acre <- read_csv("original-data/FWC_plant_survey_key_all_acreage.csv")
key_all_pres <- read_csv("original-data/FWC_plant_survey_key_all_presence.csv")
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
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

# summarize by permanentID to remove duplicates
plant_fwc3 <- plant_fwc2 %>%
  group_by(PermanentID, Area_ha, GSYear, PreCtrl, TaxonName, Habitat, Origin) %>%
  summarize(AreaName = paste(sort(unique(AreaOfInterest)), collapse = "/"),
            SurveyDate = max(SurveyDate),
            Detected = as.numeric(sum(Detected) > 0)) %>%
  ungroup()

plant_fwc3 %>%
  group_by(PermanentID, GSYear, TaxonName) %>%
  summarise(surveys = length(unique(SurveyDate))) %>%
  ungroup() %>%
  filter(surveys > 1) 
# none

# add row for every year for each site/species combo (NA's for missing surveys)
plant_fwc4 <- plant_fwc3 %>%
  full_join(plant_fwc3 %>%
              select(PermanentID, TaxonName, Origin) %>%
              unique() %>%
              expand_grid(GSYear = min(plant_fwc3$GSYear):max(plant_fwc3$GSYear))) %>%
  group_by(PermanentID, TaxonName) %>%
  arrange(GSYear) %>% 
  mutate(PrevDetected = lag(Detected), # previous year's detection
         NextDetected = lead(Detected)) %>%
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
  group_by(PermanentID, GSYear, Origin) %>%
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
  group_by(PermanentID, AreaName, GSYear, SurveyDate) %>%
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
  select(PreCtrl, PermanentID, GSYear) %>%
  unique() %>%
  group_by(PreCtrl) %>%
  summarize(TotWatYear = n(),
            TotWaterbody = n_distinct(PermanentID)) %>%
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
            OccWaterbody = n_distinct(PermanentID)) %>%
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


#### invasive plant data ####

# combine water hyacinth and lettuce percent covered
floating_cover <- inv_plant %>%
  filter(CommonName %in% c("Water hyacinth", "Water lettuce")) %>%
  group_by(PermanentID, GSYear) %>%
  summarize(across(.cols = ends_with ("PropCovered"), sum),
            InitPercCovered = sum(InitPercCovered),
            SpeciesAcres = sum(SpeciesAcres),
            EstAreaCoveredRaw_ha = sum(EstAreaCoveredRaw_ha)) %>%
  ungroup() %>%
  mutate(across(.cols = ends_with ("PropCovered"), ~ if_else(.x > 1, 1, .x)),
         InitPercCovered = if_else(InitPercCovered > 1, 1, InitPercCovered),
         CommonName = "floating plants")

# add floating cover
inv_plant2 <- inv_plant %>%
  full_join(floating_cover) %>%
  mutate(across(ends_with("PropCovered"), ~ .x * 100), # change proportion to percent
         Lag2APCsq = Lag2AvgPropCovered^2) %>% # square perc covered
  rename_with(str_replace, pattern = "PropCovered", replacement = "PercCovered") %>%
  select(CommonName, PermanentID, GSYear, ends_with("PercCovered"), Lag2APCsq, SpeciesAcres, EstAreaCoveredRaw_ha)

# Species names
inv_taxa <- tibble(Species = c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)", "Panicum repens", "Urochloa mutica", "Cyperus blepharoleptos"),
                   CommonName = c("Hydrilla", "floating plants", "Torpedograss", "Para grass", "Cuban bulrush"))


#### management data ####

# use "Species" and "unique" to get floating plants combined
# unlike water quality, use match control to native plants by GSYear (rather than GSYear - 1)
# becaue GSYear has already been calibrated to precede the plant survey
ctrl2 <- ctrl %>%
  select(Species, PermanentID, GSYear, LastTreatment, RecentTreatment, ends_with("Treated")) %>%
  unique() %>%
  left_join(inv_taxa) 


#### identify waterbodies to use ####

# has the plant ever been established?
# by using EstAreaCoveredRaw_ha, there had to be more than one year per permanentID
perm_plant <- inv_plant2 %>%
  group_by(PermanentID, CommonName) %>%
  summarize(Established = as.numeric(sum(EstAreaCoveredRaw_ha) > 0)) %>%
  ungroup() %>%
  filter(Established > 0)

# has the plant ever been treated?
perm_ctrl <- ctrl2 %>%
  group_by(PermanentID, Species) %>%
  summarize(Treatments = as.numeric(sum(Treated, na.rm = T) > 0)) %>%
  ungroup() %>%
  filter(Treatments > 0) %>%
  left_join(inv_taxa) %>%
  select(-Species)

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

# select waterbodies that have had species and at least one year of management
nat_inv <- nat_fwc3 %>%
  inner_join(perm_plant_ctrl)

# save data
write_csv(nat_fwc3, "intermediate-data/FWC_common_native_plants_invasive_species_data_formatted.csv")
write_csv(nat_inv, "intermediate-data/FWC_common_native_plants_invaded_data_formatted.csv")


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
