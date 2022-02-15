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

# species sampled somewhat continuously 
# remove focal invasive species
# remove confused origin species (not sampled continuously anyway)
plant_cont <- plant_detect %>%
  filter(FirstDetect %in% c(1982, 1983) & Survey2020 == 1 & 
           (str_detect(Notes, "meaning of this changes over time") == F | is.na(Notes)) &
           !(TaxonName %in% c("Ludwigia grandiflora/hexapetala", "Ludwigia octovalvis/peruviana",
                              "Filamentous algae")))
# Ludwigia octovalvis/peruviana removed because it combines native and non-native and Ludwigia are difficult to ID
# filamentous algae removed -- use chlorophyll to analyze algae
# 95 species

# format data
plant_fwc2 <- plant_fwc %>%
  filter(TaxonName %in% plant_cont$TaxonName) %>% # select taxa monitored continuously
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
         Detected = case_when(IsDetected == "Yes" ~ 1,
                              IsDetected == "No" ~ 0),
         SurveyMonth = month(SurveyDate),
         SurveyYear = year(SurveyDate),
         GSYear = case_when(SurveyMonth >= 4 ~ SurveyYear,
                            SurveyMonth < 4 ~ SurveyYear - 1),
         PreCtrl = if_else(GSYear < 1998, "pre ctrl data", "post ctrl data"),
         Area_ha = ShapeArea * 100) # convert lake area from km-squared to hectares

# duplicate surveys in a year
plant_fwc2 %>%
  group_by(AreaOfInterestID, GSYear) %>%
  summarise(surveys = length(unique(SurveyDate))) %>%
  ungroup() %>%
  filter(surveys > 1) 
# 41 AOIs have multiple surveys in a year

# summarize by permanentID to remove duplicates
plant_fwc3 <- plant_fwc2 %>%
  filter(!(AreaOfInterest == "Tohopekaliga, Lake" & SurveyYear == 2017)) %>% # survey seems incomplete
  group_by(PermanentID, Area_ha, GSYear, PreCtrl, TaxonName, Habitat, Origin) %>%
  summarize(AreaName = paste(sort(unique(AreaOfInterest)), collapse = "/"),
            SurveyDate = max(SurveyDate),
            Detected = if_else(sum(Detected) > 0, 1, 0)) %>%
  ungroup() %>%
  full_join(plant_fwc2 %>% # add row for every year for each site/species combo (NA's for missing surveys)
              select(PermanentID, TaxonName) %>%
              unique() %>%
              expand_grid(GSYear = min(plant_fwc2$GSYear):max(plant_fwc2$GSYear))) %>%
  group_by(PermanentID, TaxonName) %>%
  arrange(GSYear) %>% 
  mutate(PrevDetected = lag(Detected), # previous year's detection
         NextDetected = lead(Detected)) %>%
  ungroup() %>%
  mutate(Surveyed = case_when(Detected == 1 ~ 1,
                              is.na(Detected) ~ 0,
                              Detected == 0 & GSYear %in% c(1985, 1987, 1989, 1991, 1993, 2000, 2001) ~ 0, # clear dips in species richness, will correct if others in origin group were detected (see below)
                              TRUE ~ 1)) # errs on the side of counting as surveyed

# make sure missing Detected/PrevDetected applies to all taxa
plant_fwc3 %>%
  group_by(PermanentID, GSYear) %>%
  summarize(NADetected = sum(is.na(Detected)),
            TotDetected = length(Detected),
            NAPrevDetected = sum(is.na(PrevDetected)),
            TotPrevDetected = length(PrevDetected)) %>%
  ungroup() %>%
  filter((NADetected > 0 & NADetected != TotDetected) | 
           (NAPrevDetected > 0 & NAPrevDetected != TotPrevDetected))
# yes

# zero detected
plant_fwc3 %>%
  group_by(PermanentID, AreaName, GSYear, SurveyDate) %>%
  summarize(TotDetected = sum(Detected)) %>%
  ungroup() %>%
  filter(TotDetected == 0) %>%
  rename(AreaOfInterest = AreaName) %>%
  inner_join(plant_fwc) %>%
  select(AreaOfInterest, SurveyDate, TaxonName, IsDetected)
# only have one taxon and they're not in continuously monitored group

# assign NA's when a survey probably wasn't conducted
plant_origin_count <- plant_fwc3 %>%
  filter(!is.na(Origin)) %>%
  group_by(PermanentID, GSYear, Origin) %>%
  summarize(Taxa = sum(Detected),
            Surveyed = sum(Surveyed)) %>%
  ungroup() %>%
  full_join(plant_fwc3 %>% # add every combo of area and survey
              select(PermanentID, GSYear) %>%
              unique() %>%
              expand_grid(tibble(Origin = c("Native", "Exotic")))) %>%
  mutate(Taxa = replace_na(Taxa, 0),
         Surveyed = replace_na(Surveyed, 0)) %>%
  pivot_wider(names_from = Origin, 
              values_from = c(Taxa, Surveyed),
              names_glue = "{Origin}_{.value}") 

plant_origin_count %>%
  filter(Native_Taxa == 0 & Exotic_Taxa != 0) %>% # select surveys with no native species recorded
  group_by(Exotic_Taxa) %>%
  summarize(Surveys = n(),
            MaxYear = max(GSYear)) %>% # number of surveys per exotic species count
  ungroup() %>%
  mutate(TotSurveys = sum(Surveys))
# 0 exotic -- all should be NA (no survey)
# 1 - 7 exotic species, 1714 surveys
# max year = 2005

plant_origin_count %>%
  filter(Exotic_Taxa == 0 & Native_Taxa != 0) %>% # select surveys with no exotic species recorded
  group_by(Native_Taxa) %>%
  summarize(Surveys = n(),
            MaxYear = max(GSYear)) %>% # number of surveys per native species count
  ungroup() %>%
  mutate(TotSurveys = sum(Surveys)) %>%
  data.frame()
# 1 - 28 native species, 279 surveys
# max year = 2019

ggplot(plant_origin_count, aes(x = GSYear, y = Native_Taxa, color = PermanentID)) +
  geom_line() +
  theme(legend.position = "none")
# obvious years where native species weren't sampled
# remove if value is zero and value in previous and next years are not zero
# remove if value is zero and year is 2000, 2001
# used to make Surveyed column in plant_fwc3

ggplot(plant_origin_count, aes(x = GSYear, y = Exotic_Taxa, color = PermanentID)) +
  geom_line() +
  theme(legend.position = "none")
# sampled some invasive species in those times

ggplot(plant_origin_count, aes(x = Native_Surveyed, y = Native_Taxa)) +
  geom_point(alpha = 0.1)

ggplot(plant_origin_count, aes(x = Exotic_Surveyed, y = Exotic_Taxa)) +
  geom_point(alpha = 0.1)

(max_native <- max(plant_origin_count$Native_Surveyed)) # 81
(max_exotic <- max(plant_origin_count$Exotic_Surveyed)) # 17

plant_origin_count %>%
  filter(Native_Surveyed == max_native & Native_Taxa == 0) %>%
  group_by(GSYear) %>%
  count()
# 8 are in the first survey year

plant_origin_count %>%
  filter(GSYear == 1982 & Native_Surveyed == max_native & Native_Taxa == 0) %>%
  select(PermanentID) %>%
  inner_join(plant_origin_count %>%
               filter(GSYear == 1983))
# all have at least 4 taxa the following year

plant_origin_count %>%
  filter(Exotic_Surveyed == max_exotic & Exotic_Taxa == 0) %>%
  group_by(GSYear) %>%
  count()

plant_origin_count %>%
  filter(GSYear == 1982 & Exotic_Surveyed == max_exotic & Exotic_Taxa == 0) %>%
  select(PermanentID) %>%
  inner_join(plant_origin_count %>%
               filter(GSYear == 1983)) %>%
  group_by(Exotic_Taxa) %>%
  count()

# check NA's
plant_fwc3 %>%
  filter(is.na(Origin)) %>%
  select(Detected) %>%
  unique()
# should all be NA's

# remove rows for no surveys
plant_fwc4 <- plant_fwc3 %>%
  left_join(plant_origin_count) %>%
  mutate(Surveyed = case_when(Native_Taxa > 0 & Origin == "Native" ~ 1, # if any taxa of that origin were detected, assume survey occurred
                              Exotic_Taxa > 0 & Origin == "Exotic" ~ 1,
                              GSYear == 1982 & Native_Surveyed == max_native & Native_Taxa == 0 & Origin == "Native" ~ 0, # remove first year surveys with no detects
                              GSYear == 1982 & Exotic_Surveyed == max_exotic & Exotic_Taxa == 0 & Origin == "Exotic" ~ 0,
                              TRUE ~ Surveyed),
         Detected = if_else(Surveyed == 0, NA_real_, Detected)) %>% # make detections NA if no survey
  group_by(PermanentID, TaxonName) %>%
  arrange(GSYear) %>% 
  mutate(PrevDetected = lag(Detected), # redo
         NextDetected = lead(Detected)) %>%
  ungroup() %>%
  select(-c(Native_Taxa, Native_Surveyed, Exotic_Taxa, Exotic_Surveyed))

# check surveyed/detected
plant_origin_check <- plant_fwc4 %>%
  filter(!is.na(Origin)) %>%
  group_by(PermanentID, GSYear, Origin) %>%
  summarize(Surveyed = sum(Surveyed),
            Detected = sum(Detected, na.rm = T)) %>%
  ungroup()

ggplot(plant_origin_check, aes(x = Surveyed, y = Detected)) +
  geom_point() +
  facet_wrap(~ Origin, scales = "free")
# x-axis should be all 0 or max

# remake richness figures
plant_fwc4 %>%
  filter(!is.na(Origin)) %>%
  group_by(Origin, PermanentID, GSYear) %>%
  summarize(Taxa = sum(Detected)) %>%
  ungroup() %>%
  ggplot(aes(x = GSYear, y = Taxa, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ Origin, scales = "free") +
  theme(legend.position = "none")

# save
write_csv(plant_fwc4, "intermediate-data/FWC_plant_community_formatted.csv")


#### most common native taxa ####

# remove non-native taxa
nat_fwc <- plant_fwc4 %>%
  filter(Origin == "Native" & Surveyed == 1)

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
# taxa must be present in over 7% of water-year occurrences 
# and over 20% of lakes (this requirement is met with above)
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
