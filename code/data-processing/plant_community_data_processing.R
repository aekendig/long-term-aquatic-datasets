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


#### edit native plant data ####

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

# species samples continuously 
# remove focal invasive species
# remove confused origin species
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
