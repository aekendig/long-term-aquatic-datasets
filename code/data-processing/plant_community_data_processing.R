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
# plant_fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv",
#                        col_types = list(Depth_ft = col_double(),
#                                         PermanentID = col_character(),
#                                         YearF = col_character())) # add in later, need full list over time
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
                              "Ludwigia grandiflora/hexapetala")))
# 97 species

# select all species except invasive
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


#### most common taxa ####

# total waterbody-year combos
TotWatYear <- plant_fwc3 %>%
  select(PreCtrl, PermanentID, GSYear) %>%
  unique() %>%
  group_by(PreCtrl) %>%
  summarize(TotWatYear = n(),
            TotWaterbody = n_distinct(PermanentID)) %>%
  ungroup()

# total habitat types
TotHabitat <- plant_fwc3 %>%
  group_by(Habitat) %>%
  summarize(TotTaxa = n_distinct(TaxonName)) %>%
  ungroup() %>%
  mutate(Taxa40 = round(TotTaxa * 0.4))

# summarize waterbody-year occurrences
common_fwc <- plant_fwc3 %>%
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
         Over50Waterbody = if_else(RatioWaterbody > 0.5, 1, 0), # thresholds based on post ctrl data distributions (below)
         Over25WatYear = if_else(RatioWatYear > 0.25, 1, 0)) # these are already split by PreCtrl

# distribution of occurrences
ggplot(common_fwc, aes(x = RankWatYear, y = RatioWatYear)) +
  geom_line() +
  geom_point(aes(color = Over50Waterbody)) +
  facet_wrap(~ PreCtrl)

ggplot(common_fwc, aes(x = RankWaterbody, y = RatioWaterbody)) +
  geom_line() +
  geom_point(aes(color = Over25WatYear)) +
  facet_wrap(~ PreCtrl)

ggplot(common_fwc, aes(x = RankWatYearHab, y = RatioWatYear)) +
  geom_line() +
  geom_point(aes(color = Over50Waterbody)) +
  facet_grid(Habitat ~ PreCtrl)

ggplot(common_fwc, aes(x = RankWaterbodyHab, y = RatioWaterbody)) +
  geom_line() +
  geom_point(aes(color = Over25WatYear)) +
  facet_grid(Habitat ~ PreCtrl)

# common species in post control data
# selected top 40% of each group, but the submersed were much rarer than the others
# taxa must be present in over 25% of water-year occurrences 
# and over 50% of lakes (this requirement is met with above)
common_fwc2 <- common_fwc %>%
  left_join(TotHabitat) %>%
  # filter(PreCtrl == "post ctrl data" & RankWatYearHab <= Taxa40) %>%
  filter(PreCtrl == "post ctrl data" & Over25WatYear == 1 & Over50Waterbody) %>%
  select(TaxonName) %>%
  inner_join(common_fwc)

common_fwc2 %>%
  group_by(Habitat) %>%
  summarize(Taxa = n_distinct(TaxonName),
            MinWatYear = min(OccWatYear),
            MinWaterbody = min(OccWaterbody)) %>%
  left_join(TotHabitat) %>%
  mutate(RatioTaxa = Taxa / TotTaxa)

ggplot(common_fwc2, aes(x = RankWatYearHab, y = RatioWatYear)) +
  geom_line() +
  geom_point(aes(color = Over50Waterbody)) +
  facet_grid(Habitat ~ PreCtrl)

# select taxa
plant_fwc4 <- plant_fwc3 %>%
  filter(TaxonName %in% common_fwc2$TaxonName)

# save
write_csv(plant_fwc4, "intermediate-data/FWC_common_plant_community_formatted.csv")
