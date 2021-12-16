#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)
library(lubridate)

# import data
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv",
                      col_types = list(JoinNotes = col_character(),
                                       PermanentID = col_character()))
plant_fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv",
                       col_types = list(Depth_ft = col_double(),
                                        PermanentID = col_character(),
                                        YearF = col_character()))
plant_detect <- read_csv("intermediate-data/FWC_plant_survey_first_detection_manual.csv")


#### edit plant community data ####

# one origin per species?
plant_fwc %>%
  group_by(TaxonName) %>%
  summarise(orig = length(unique(Origin))) %>%
  ungroup() %>%
  filter(orig != 1)

plant_fwc %>%
  filter(TaxonName == "Ludwigia grandiflora/hexapetala") %>%
  select(SurveyYear, Origin) %>%
  unique()
# early records have native origin
# early records called this species Ludwigia grandifolia (see fwc_plant_data_processing),
# which I can't confirm is an actual species
# remove from native species analysis

# Eppc and origin
plant_fwc %>%
  select(Eppc, Origin) %>%
  unique()

# first detection across all lakes
first_detect <- plant_fwc %>%
  group_by(TaxonName) %>%
  summarise(FirstYear = min(SurveyYear)) %>%
  ungroup()

first_detect %>%
  ggplot(aes(x = FirstYear)) +
  geom_bar()

# use species from first two years
sort(unique(first_detect$FirstYear))[1:2]

# save data to add continuous detection
# write_csv(first_detect, "intermediate-data/FWC_plant_survey_first_detection.csv")
# manual version imported above

# native species data missing 2000-20001
plant_fwc %>%
  filter(Origin == "Native") %>% 
  filter(SurveyYear %in% c(2000, 20001) & IsDetected == "Yes")

# species sampled continuously 
# remove focal invasive species
# remove confused origin species (not sampled continuously anyway)
plant_cont <- plant_detect %>%
  filter(FirstDetect %in% c(1982, 1983) & Survey2020 == 1 & 
           (str_detect(Notes, "meaning of this changes over time") == F | is.na(Notes)) &
           !(TaxonName %in% c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes",
                              "Ludwigia grandiflora/hexapetala", "Ludwigia octovalvis/peruviana",
                              "Filamentous algae")))
# Ludwigia octovalvis/peruviana removed because it combines native and non-native and Ludwigia are difficult to ID
# filamentous algae removed -- use chlorophyll to analyze algae
# 95 species

# format data
plant_fwc2 <- plant_fwc %>%
  filter(TaxonName %in% plant_cont$TaxonName) %>%
  select(TaxonName, Habitat, Origin) %>%
  unique() %>% # full species list
  expand_grid(plant_fwc %>%
                select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate) %>%
                unique()) %>% # full survey list (row for every species in every survey)
  full_join(plant_fwc %>%
              filter(TaxonName %in% plant_cont$TaxonName) %>%
              select(TaxonName, Habitat, Origin,
                     AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, 
                     IsDetected)) %>% # add detection data (only "Yes")
  mutate(IsDetected = replace_na(IsDetected, "No"),
         IsDetected = case_when(year(SurveyDate) %in% c(2000, 2001) ~ NA_character_, # no native species these years
                                TRUE ~ IsDetected),
         Detected = case_when(IsDetected == "Yes" ~ 1,
                              IsDetected == "No" ~ 0),
         SurveyMonth = month(SurveyDate),
         SurveyYear = year(SurveyDate),
         GSYear = case_when(SurveyMonth >= 4 ~ SurveyYear,
                            SurveyMonth < 4 ~ SurveyYear - 1),
         PreCtrl = if_else(GSYear < 1998, "pre ctrl data", "post ctrl data"),
         Area_ha = ShapeArea * 100) %>% # convert lake area from km-squared to hectares
  filter(!is.na(IsDetected)) # remove 2000-2001 surveys

# duplicate surveys in a year
plant_fwc2 %>%
  group_by(AreaOfInterestID, GSYear) %>%
  summarise(surveys = length(unique(SurveyDate))) %>%
  ungroup() %>%
  filter(surveys > 1) 
# 40 AOIs have multiple surveys in a year

# summarize by permanentID to remove duplicates
plant_fwc3 <- plant_fwc2 %>%
  group_by(PermanentID, Area_ha, GSYear, PreCtrl, TaxonName, Habitat, Origin) %>%
  summarize(AreaName = paste(AreaOfInterest, collapse = "/"),
            SurveyDate = max(SurveyDate),
            Detected = if_else(sum(Detected) > 0, 1, 0)) %>%
  ungroup()

# save
write_csv(plant_fwc3, "intermediate-data/FWC_plant_community_formatted.csv")


#### most common native taxa ####

# remove non-native taxa
nat_fwc <- plant_fwc3 %>%
  filter(Origin == "Native")

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
# only 4 floating taxa

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
# selected top 40% of each group, but the submersed were much rarer than the others
# taxa must be present in over 25% of water-year occurrences 
# and over 50% of lakes (this requirement is met with above)
common_fwc2 <- common_fwc %>%
  # left_join(TotHabitat) %>%
  # filter(PreCtrl == "post ctrl data" & RankWatYearHab <= Taxa40) %>%
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

n_distinct(common_fwc2$TaxonName) # 62 taxa total

ggplot(common_fwc2, aes(x = RankWatYearHab, y = RatioWatYear)) +
  geom_line() +
  geom_point(aes(color = Over20Waterbody)) +
  facet_grid(Habitat ~ PreCtrl)

# select taxa
nat_fwc2 <- nat_fwc %>%
  filter(TaxonName %in% common_fwc2$TaxonName)

# save
write_csv(nat_fwc2, "intermediate-data/FWC_common_native_plants_formatted.csv")


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
