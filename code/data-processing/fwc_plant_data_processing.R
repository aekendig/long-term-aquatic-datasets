#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(data.table)
library(tidyverse)
library(janitor)
library(taxize)
library(lubridate)

# load data
fwc_plant <- read_csv("original-data/FWC Plant Surveys.csv")
fwc_plant_new <- read_csv("original-data/FWCSurvey_2018_2019.csv")
gis <- read_csv("intermediate-data/gis_fwc_lakewatch_fwri.csv",
                col_types = list(wkt_geom = col_character(),
                                 AOI = col_character(),
                                 Lake = col_character()))
ctrl_old <- read_csv("original-data/PrePMARS_IPMData.csv")
ctrl <- read_csv("original-data/FWC_Herbicide_Treatments_mid2010-Oct2020.csv")
fwc_ID <- read_csv("gis/data/FWC_Replaced_Coordinates.csv")
eppc_1999 <- read_csv("original-data/EPPC_1999.csv")
key_all_acre <- read_csv("original-data/FWC_plant_survey_key_all_acreage.csv")
key_exotic_acre <- read_csv("original-data/FWC_plant_survey_key_exotic_acreage.csv")
key_all_pres <- read_csv("original-data/FWC_plant_survey_key_all_presence.csv")


#### resolve species names ####

# list of taxa
taxa_list <- fwc_plant %>%
  select(SpeciesName) %>%
  unique() %>%
  full_join(fwc_plant_new %>%
              select(SpeciesName) %>%
              unique()) %>%
  unique() %>%
  mutate(Taxon = str_replace_all(SpeciesName, " spp.| spp| sp.| sp", ""), # for genus-level
         Taxon = str_replace_all(Taxon, ", sub| \\(exotic\\)| \\(other natives\\)| \\(other\\)|, natives|, emersed", ""), # for origin/growth type
         Taxon = str_replace(Taxon, "\\/.*", ""), # anything after /
         Taxon = str_replace(Taxon, "sub\\.", "ssp."),
         words = str_count(Taxon, pattern = " "))

# duplicate genera
taxa_list %>%
  mutate(dup = duplicated(Taxon)) %>%
  filter(dup == T) %>%
  select(Taxon) %>%
  inner_join(taxa_list) %>%
  data.frame() %>%
  unique()
# all are fine to use taxon

# select species
species_list <- taxa_list %>%
  filter(words > 0) # synonyms function failed with genera

# synonyms from ITIS
# taxa_syn <- synonyms(species_list$Taxon, db = "itis")
# manually accepted duplicate names
# chose accepted ones and checked on website
# date run: 8/12/21

# make into dataframe
# taxa_syn2 <- rbindlist(lapply(taxa_syn, as.data.table), use.names = T, fill = T, idcol = "Taxon") %>%
#   select(-V1) %>%
#   as_tibble()

# manually check missing species
# taxa_syn2 %>%
#   filter(is.na(sub_tsn)) %>%
#   select(Taxon)
# one is general algae
# Pomacea insularum is the snail -- remove
# type-os: 
# Ricciocarpus natans = Ricciocarpos natans
# Ludwigia grandifolia = Ludwigia grandiflora
# Triadenum virginium = Triadenum virginicum
# synonyms:
# Symphyotrichum carolinianum = Ampelaster carolinianus
# Cyperus blepharoleptos = Oxycaryum cubense (the taxon name is accepted in FL Plant Atlas and GBIF, just not in ITIS)

# alternate names
# last one was left out of results
species_list2 <- tibble(Taxon = c("Ricciocarpus natans", "Ludwigia grandifolia", "Triadenum virginium", "Symphyotrichum carolinianum", "Cyperus blepharoleptos", "Potamogeton amplifolius"),
                        Alt_taxon = c("Ricciocarpos natans", "Ludwigia grandiflora", "Triadenum virginicum", "Ampelaster carolinianus", "Oxycaryum cubense", "Potamogeton amplifolius"))

# check for correct spelling/accepted names
species_list %>%
  filter(Taxon %in% species_list2$Alt_taxon)

# synonyms from ITIS
# taxa_syn3 <- synonyms(species_list2$Alt_taxon, db = "itis")
# manually accepted duplicate names
# chose accepted ones and checked on website
# date run: 8/12/21

# make into dataframe
# taxa_syn4 <- rbindlist(lapply(taxa_syn3, as.data.table), use.names = T, fill = T, idcol = "Taxon") %>%
#   as_tibble() %>%
#   full_join(taxa_syn2)

# no synonyms (manually double checked)
# Potamogeton amplifolius
# Schinus terebinthifolius
# Lygodium japonicum
# Lygodium microphyllum
# Cyperus alopecuroides
# Cyperus papyrus
# Scleria lacustris
# Luziola subintegra

# save
# write_csv(taxa_syn4, "intermediate-data/FWC_plant_survey_species_synonyms.csv")

# import
taxa_syn4 <- read_csv("intermediate-data/FWC_plant_survey_species_synonyms.csv")

# check for matches between synonyms and taxa
taxa_syn4 %>%
  filter(syn_name %in% species_list$Taxon & syn_name != Taxon)
# Paspalum repens is the accepted name for Paspalum fluitans

# check for other snails
species_list %>%
  filter(str_detect(Taxon, "Pomacea"))
# Pomacea paludosa, Pomacea insularum

# survey numbers
taxa_surveys <- fwc_plant %>%
  select(SpeciesName, SurveyYear, WaterbodyName) %>%
  full_join(fwc_plant_new %>%
              select(SpeciesName, SurveyYear, WaterbodyName))

# names to combine
taxa_surveys %>%
  filter(SpeciesName %in% c("Paspalum repens", "Paspalum fluitans")) %>%
  ggplot(aes(x = SurveyYear, fill = SpeciesName)) +
  geom_bar()
# changed name after first year, used once after that -- combine

taxa_surveys %>%
  filter(SpeciesName %in% c("Ludwigia grandiflora/hexapetala", "Ludwigia grandifolia")) %>%
  ggplot(aes(x = SurveyYear, fill = SpeciesName)) +
  geom_bar()
# incorrect spelling only used in first three surveys and once after
# leave as separate since it could be referring to a different species

taxa_surveys %>%
  filter(SpeciesName %in% c("Cyperus blepharoleptos", "Oxycaryum cubense")) %>%
  ggplot(aes(x = SurveyYear, fill = SpeciesName)) +
  geom_bar()
# name switched to Cyperus recently

taxa_surveys %>%
  filter(SpeciesName %in% c("Sagittaria stagnorum", "Sagittaria subulata/graminea/gracillima")) %>%
  ggplot(aes(x = SurveyYear, fill = SpeciesName)) +
  geom_bar()
# stagnorum used in same years as other, it's less common

taxa_surveys %>%
  filter(SpeciesName %in% c("Echinochloa spp.", "Echinochloa spp. (exotic)", "Echinochloa spp")) %>%
  ggplot(aes(x = SurveyYear, fill = SpeciesName)) +
  geom_bar()
# no period used in last two years, otherwise exotic and non co-occur in time

taxa_surveys %>%
  filter(SpeciesName %in% c("Bidens spp.", "Bidens spp")) %>%
  ggplot(aes(x = SurveyYear, fill = SpeciesName)) +
  geom_bar()
# no period used in last two years

taxa_surveys %>%
  filter(SpeciesName %in% c("Fuirena spp.", "Fuirena spp")) %>%
  ggplot(aes(x = SurveyYear, fill = SpeciesName)) +
  geom_bar()
# periods seem to be left out of recent surveys

taxa_surveys %>%
  filter(SpeciesName %in% c("Ipomoea sp.", "Ipomoea sp")) %>%
  ggplot(aes(x = SurveyYear, fill = SpeciesName)) +
  geom_bar()

# check for potential duplication in survey
taxa_surveys %>%
  filter(SpeciesName == "Paspalum fluitans") %>%
  select(SurveyYear, WaterbodyName) %>%
  inner_join(taxa_surveys) %>%
  filter(SpeciesName == "Paspalum repens")
# no duplication

taxa_surveys %>%
  filter(SpeciesName == "Sagittaria stagnorum") %>%
  select(SurveyYear, WaterbodyName) %>%
  inner_join(taxa_surveys) %>%
  filter(SpeciesName == "Sagittaria subulata/graminea/gracillima")
# yes, both options were on a survey

taxa_surveys %>%
  filter(SpeciesName == "Ludwigia grandifolia") %>%
  select(SurveyYear, WaterbodyName) %>%
  inner_join(taxa_surveys) %>%
  filter(SpeciesName == "Ludwigia grandiflora/hexapetala")
# yes, both options were on a survey

# combine names
# remove snails
taxa_list_mod <- fwc_plant %>%
  select(SpeciesName) %>%
  unique() %>%
  full_join(fwc_plant_new %>%
              select(SpeciesName) %>%
              unique()) %>%
  unique() %>%
  mutate(TaxonName = case_when(SpeciesName == "Oxycaryum cubense" ~ "Cyperus blepharoleptos",
                               SpeciesName == "Paspalum fluitans" ~ "Paspalum repens",
                               SpeciesName == "Ipomoea sp" ~ "Ipomoea sp.",
                               TRUE ~ SpeciesName),
         TaxonName = str_replace_all(TaxonName, "spp", "spp."),
         TaxonName = str_replace_all(TaxonName, "spp..", "spp.")) %>%
  filter(str_detect(SpeciesName, "Pomacea") == F)

# update fwc data
fwc_plant2 <- fwc_plant %>%
  left_join(taxa_list_mod)

fwc_plant_new2 <- fwc_plant_new %>%
  left_join(taxa_list_mod)

# look at missing taxa
filter(fwc_plant2, is.na(TaxonName)) %>% 
  select(SpeciesName) %>% 
  unique()
# snail

filter(fwc_plant_new2, is.na(TaxonName)) %>% 
  select(SpeciesName) %>% 
  unique()
# snails and NA


#### edit new FWC plant survey data ####

# format date
# rename columns
# add permanent ID based on AreaOfInterestID
fwc_plant_new3 <- fwc_plant_new2 %>%
  mutate(SurveyDate = as.Date(SurveyDate, "%m/%d/%y"),
         County = toupper(County)) %>%
  rename("AreaOfInterest" = "WaterbodyName") %>%
  left_join(gis %>%
              filter(CoordSource %in% c("FWC", "FWC_2")) %>%
              select(AreaOfInterestID, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique())

# duplicate rows?
fwc_plant_new3 %>%
  get_dupes()
# none

# missing ID?
fwc_plant_new3 %>%
  filter(is.na(PermanentID)) %>%
  select(AreaOfInterestID, WaterbodyType) %>%
  unique() %>%
  group_by(WaterbodyType) %>%
  count()
# 5 lakes

fwc_plant_new3 %>%
  filter(is.na(PermanentID) & WaterbodyType == "Lake") %>%
  select(AreaOfInterest) %>%
  unique()
# all are springs, creeks, run, reservoirs

# duplicate AOI ID's?
fwc_plant_new3 %>%
  distinct(AreaOfInterest, AreaOfInterestID) %>%
  get_dupes()
# no

# remove missing perm ID
fwc_plant_new4 <- fwc_plant_new3 %>%
  filter(!is.na(PermanentID))


#### FWC AreaOfInterestIDs ####

# this dataset is missing AOI ID's
# get ID's from other FWC datasets
fwc_ID2 <- ctrl_old %>%
  distinct(AreaOfInterest, County, AreaOfInterestID) %>%
  rbind(ctrl %>%
          distinct(AreaOfInterest, County, AreaOfInterestID)) %>%
  rbind(fwc_plant_new %>%
          rename("AreaOfInterest" = "WaterbodyName") %>%
          distinct(AreaOfInterest, County, AreaOfInterestID)) %>%
  rbind(fwc_ID %>%
          distinct(AreaOfInterest, County, AreaOfInterestID)) %>%
  mutate(County = toupper(County)) %>%
  unique()

# duplicate AreaOfInterest names/counties
fwc_ID2 %>%
  mutate(AOI_County = paste(AreaOfInterest, County, sep = "_")) %>%
  get_dupes("AOI_County") %>%
  select(AreaOfInterestID) %>%
  inner_join(fwc_ID2)
# 44 has two different counties, one is the same as 45
# 714 is a private lake
# 218 has two different counties, one is the same as 219
# 872 is a private lake

filter(ctrl_old, AreaOfInterestID == 44) %>% select(County) %>% unique()
filter(ctrl, AreaOfInterestID == 44) %>% select(County) %>% unique()
filter(fwc_plant_new, AreaOfInterestID == 44) %>% select(County) %>% unique()
# 44 is in Highlands, not Polk

filter(ctrl_old, AreaOfInterestID == 218) %>% select(County) %>% unique()
filter(ctrl, AreaOfInterestID == 218) %>% select(County) %>% unique()
filter(fwc_plant_new, AreaOfInterestID == 218) %>% select(County) %>% unique()
# 218 is in Alachua, not Clay

# duplicate ID's
get_dupes(fwc_ID2, AreaOfInterestID) %>%
  group_by(AreaOfInterestID) %>%
  mutate(names = n_distinct(AreaOfInterest)) %>%
  ungroup()
# 387 has multiple spellings
# 44, 218 fixed above
# 60 is in both Orange and Lake counties -- more is in Orange
# 414 has two different county spellings
# 465 two spellings
# 975 there are two Saddleback Lakes, one in each county, they're both private

filter(fwc_ID2, County %in% c("SAINT LUCIE", "ST. LUCIE"))
# 387: use saint

filter(fwc_ID2, County %in% c("SAINT JOHNS", "ST JOHNS"))
# 414: use saint

filter(fwc_plant_new4, AreaOfInterestID == 465)
# 465: use 1 "l" for 465

filter(fwc_plant_new4, AreaOfInterestID == 975)
# not in the new dataset
filter(fwc_plant2, str_detect(WaterbodyName, "Saddleback"))
# not in new dataset

# remove duplicates
fwc_ID3 <- fwc_ID2 %>%
  filter(!(AreaOfInterestID == 44 & County == "POLK") &
           !(AreaOfInterestID == 218 & County == "CLAY") &
           !(AreaOfInterestID %in% c(714, 872)) &
           !(AreaOfInterestID == 387 & County == "ST. LUCIE") &
           !(AreaOfInterestID == 387 & AreaOfInterest == "Savannas Preserve State Park") &
           !(AreaOfInterestID == 60 & County == "LAKE") &
           !(AreaOfInterestID == 414 & County == "ST JOHNS") &
           !(AreaOfInterestID == 465 & AreaOfInterest == "Watermellon Pond"))

# check duplicates
fwc_ID3 %>%
  mutate(AOI_County = paste(AreaOfInterest, County, sep = "_")) %>%
  get_dupes("AOI_County") 

get_dupes(fwc_ID3, AreaOfInterestID)

# missing IDs?
fwc_ID3 %>%
  filter(is.na(AreaOfInterestID))
# none


#### old FWC plant survey data ####

# matches in ID dataset
fwc_plant2 %>%
  rename("AreaOfInterest" = "WaterbodyName") %>%
  select(AreaOfInterest, County) %>%
  unique() %>%
  mutate(County = toupper(County)) %>%
  anti_join(fwc_ID3)
# all are matched, except Watermellon

# add ID's
# format dates
# add permanent ID based on AreaOfInterestID
# remove missing ID
# remove duplicate rows
fwc_plant3 <- fwc_plant2 %>%
  mutate(SurveyDate = as.Date(SurveyDate, "%m/%d/%Y"),
         County = toupper(County),
         WaterbodyName = str_replace(WaterbodyName, "Watermellon Pond",
                                     "Watermelon Pond")) %>%
  rename("AreaOfInterest" = "WaterbodyName") %>%
  left_join(fwc_ID3) %>%
  left_join(gis %>%
              filter(CoordSource %in% c("FWC", "FWC_2")) %>%
              select(AreaOfInterestID, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique())

# all have AOI ID's?
filter(fwc_plant3, is.na(AreaOfInterestID))

# duplicate rows?
fwc_plant3 %>%
  get_dupes()
# 415

# remove duplicate rows
fwc_plant4 <- fwc_plant3 %>%
  unique()

# missing ID?
fwc_plant4 %>%
  filter(is.na(PermanentID)) %>%
  select(AreaOfInterestID, WaterbodyType) %>%
  unique() %>%
  group_by(WaterbodyType) %>%
  count()
# 7 lakes

fwc_plant4 %>%
  filter(is.na(PermanentID) & WaterbodyType == "Lake") %>%
  select(AreaOfInterest) %>%
  unique()
# all are springs, creeks, run, reservoirs

# remove missing perm ID
fwc_plant5 <- fwc_plant4 %>%
  filter(!is.na(PermanentID))


#### combine old and new datasets ####

# overlapping surveys?
fwc_plant_new4 %>%
  select(AreaOfInterestID, SurveyDate) %>%
  inner_join(fwc_plant5 %>%
               select(AreaOfInterestID, SurveyDate))
# no

# duplicate measurements?
fwc_plant5 %>%
  select(AreaOfInterestID, SurveyDate, TaxonName) %>%
  get_dupes() %>%
  select(TaxonName) %>%
  unique() %>%
  data.frame()
# yes

# remove rows with missing species names
# sum cover for duplicate non-species
# use max cover for duplicate species
fwc_plant6 <- fwc_plant5 %>%
  full_join(fwc_plant_new4) %>%
  filter(!is.na(SpeciesName)) %>%
  group_by(AreaOfInterestID, AreaOfInterest, WaterbodyAcres, WaterbodyType, County, WMD, 
           PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource,
           SurveyYear, SurveyDate, IsUnableToSurvey, Surveyor,
           SpeciesName, TaxonName, IsDetected, IsAcreageRequire, Origin, Eppc, Habitat, HabitatShortName) %>%
  summarise(SumSpeciesAcres = sum(SpeciesAcres, na.rm = T),
            MaxSpeciesAcres = max(SpeciesAcres, na.rm = T)) %>%
  ungroup() %>%
  mutate(SpeciesAcres = if_else(str_detect(TaxonName, "Filamentous algae|other|spp.|/") == T,
                                SumSpeciesAcres,
                                MaxSpeciesAcres),
         SpeciesAcres = ifelse(SpeciesAcres == -Inf, NA_real_, SpeciesAcres)) %>%
  select(-c(SumSpeciesAcres, MaxSpeciesAcres))
 # will give warnings, but these are addressed by making -Inf into NAs

# check for duplicates
(dup_fwc_plant6 <- get_dupes(fwc_plant6))

 # number of waterbodies
n_distinct(fwc_plant6$AreaOfInterestID)


#### first records/continuous surveying ####

# first records
fwc_plant6 %>%
  filter(IsDetected == "Yes") %>%
  group_by(TaxonName) %>%
  summarise(FirstDetect = min(SurveyYear)) %>%
  ungroup() %>%
  ggplot(aes(x = FirstDetect)) +
  geom_bar()

# save for manual comparison with printed survey
write_csv(fwc_plant6 %>%
            group_by(TaxonName) %>%
            summarise(FirstDetect = min(SurveyYear)) %>%
            ungroup(),
          "intermediate-data/FWC_plant_survey_first_detection.csv")

# continued detection (in 2020 survey)
plant_detect <- read_csv("intermediate-data/FWC_plant_survey_first_detection_manual.csv")

# look at notes
plant_detect %>%
  select(Notes) %>%
  unique()

# species sampled somewhat continuously 
# remove confused origin species (not sampled continuously anyway)
plant_cont <- plant_detect %>%
  filter(FirstDetect %in% c(1982, 1983) & Survey2020 == 1 & 
           (str_detect(Notes, "meaning of this changes over time") == F | is.na(Notes)))
# 100 taxa


#### edit keys ####

# taxa with acreage
taxa_acres <- fwc_plant6 %>%
  select(TaxonName, Origin, SurveyYear, AreaOfInterest, SpeciesAcres) %>%
  filter(!is.na(SpeciesAcres))

# are any taxa not in the EPPC list?
non_eppc_taxa <- taxa_acres %>%
  filter(SurveyYear %in% c(1999:2002, 2004:2016)) %>% # years EPPC list was used
  left_join(eppc_1999 %>%
              mutate(eppc = 1)) %>%
  filter(is.na(eppc)) %>%
  select(TaxonName, Origin) %>%
  unique()
# 143 taxa, some are native

# are their synonyms on the EPPC list?
non_eppc_taxa %>%
  inner_join(taxa_syn4 %>%
               rename(TaxonName = Taxon)) %>%
  inner_join(eppc_1999 %>%
               rename(syn_name = TaxonName))
# no

# acreage for all taxa
key_all_acre2 <- key_all_acre %>%
  expand_grid(plant_cont %>%
                select(TaxonName)) %>%
  mutate(AcreageSurveyed = 1,
         PresenceSurveyed = 1)

# acreage for exotic taxa
key_exotic_acre2 <- key_exotic_acre %>%
  inner_join(plant_cont %>%
               select(TaxonName)) %>% # some of the EPPC taxa are not aquatic
  mutate(AcreageSurveyed = 1,
         PresenceSurveyed = 1)

# make sure no aquatic taxa are missing
key_exotic_acre %>%
  select(TaxonName) %>%
  unique() %>%
  anti_join(key_exotic_acre2) %>%
  data.frame()
# checked all in Google images and they all seem terrestrial

# presence for all taxa
key_all_pres2 <- key_all_pres %>%
  expand_grid(plant_cont %>%
                select(TaxonName)) %>%
  mutate(PresenceSurveyed = 1)

# identify invasive species of major concern (2017-2021)
key_major_concern <- fwc_plant2 %>%
  filter(SurveyYear %in% 2017:2021 & !is.na(SpeciesAcres) & Origin == "Exotic") %>%
  select(WaterbodyName, County, SurveyYear, TaxonName) %>%
  full_join(fwc_plant_new2 %>%
              filter(SurveyYear %in% 2017:2021 & !is.na(SpeciesAcres) & Origin == "Exotic") %>%
              select(WaterbodyName, County, SurveyYear, TaxonName)) %>%
  mutate(Waterbody = paste(WaterbodyName, County)) %>%
  group_by(TaxonName) %>%
  summarize(Waterbodies = n_distinct(Waterbody),
            Years = n_distinct(SurveyYear),
            Tot = n())

# visualize to determine cut-off
ggplot(key_major_concern, aes(x = Waterbodies, y = Tot, color = as.factor(Years))) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  geom_vline(xintercept = 100, linetype = "dashed") +
  geom_point(size = 0.75, alpha = 0.5) +
  geom_text(aes(label = TaxonName), size = 2, hjust = 0, vjust = 0) +
  coord_cartesian(xlim = c(0, 450))

# zoom in on less common taxa
ggplot(key_major_concern, aes(x = Waterbodies, y = Tot, color = as.factor(Years))) +
  geom_point(size = 0.75, alpha = 0.5) +
  geom_text(aes(label = TaxonName), size = 2, hjust = 0, vjust = 0) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100))

# key for concern
key_major_concern2 <- key_major_concern %>%
  filter(Waterbodies > 100 & Tot > 100) %>%
  select(TaxonName) %>%
  mutate(PresenceSurveyed = 1,
         AcreageSurveyed = 1) %>%
  expand_grid(SurveyYear = 2017:max(fwc_plant_new2$SurveyYear))

# save
write_csv(key_major_concern2, "intermediate-data/fwc_invasive_taxa_major_concern.csv")

# combine keys
keys <- key_all_acre2 %>%
  full_join(key_exotic_acre2) %>%
  full_join(key_major_concern2) %>%
  full_join(key_all_pres2) %>%
  mutate(AcreageSurveyed = replace_na(AcreageSurveyed, 0))


#### add in non-detects and zero abundances ####

# list of all surveys
area_yr <- fwc_plant6 %>%
  select(WaterbodyAcres, WaterbodyType, County, WMD, 
         SurveyYear, SurveyDate, Surveyor,
         AreaOfInterest, AreaOfInterestID, PermanentID, GNISID, GNISName,
         Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
  unique()

# check unable to survey
unique(fwc_plant6$IsUnableToSurvey)
# all no
# that could be the reason for some non-detects

# see if SpeciesAcres = 0 when IsDetect is NA
fwc_plant6 %>%
  filter(SpeciesAcres == 0 & is.na(IsDetected))
# no - can leave out of the case_when

# expand keys for surveys in dataset
keys2 <- keys %>% 
  left_join(fwc_plant6 %>%
              select(TaxonName, Origin, Eppc, Habitat, HabitatShortName) %>%
              unique()) %>%
  inner_join(area_yr)

# add surveyed columns
fwc_plant7 <- fwc_plant6 %>%
  filter(TaxonName %in% plant_cont$TaxonName & SurveyYear > 1982) %>% # select for continuously monitored taxa
  full_join(keys2) %>%
  mutate(IsDetected = case_when(is.na(IsDetected) & SpeciesAcres > 0 ~ "Yes",
                                is.na(IsDetected) & PresenceSurveyed == 1 ~ "No", # covers all the rest
                                TRUE ~ IsDetected),
         SpeciesAcres = case_when(is.na(SpeciesAcres) & AcreageSurveyed == 1 ~ 0,
                                  SpeciesAcres == 0 & AcreageSurveyed == 0 ~  NA_real_, # a lot of taxa have 0 species acres in years they weren't recorded as surveyed
                                  TRUE ~ SpeciesAcres),
         AcreageSurveyed = if_else(!is.na(SpeciesAcres), 1, 0),
         PresenceSurveyed = if_else(!is.na(IsDetected), 1, 0))

# check for duplication
get_dupes(fwc_plant7)

# check for missing values
filter(fwc_plant7, is.na(IsDetected))
filter(fwc_plant7, is.na(SpeciesAcres) & AcreageSurveyed == 1)


#### identify problematic surveys ####

# empty table
outliers <- tibble(AreaOfInterest = NA,
                   AreaOfInterestID = NA,
                   SurveyYear = NA,
                   DevRich = NA)

# list of AOIs
AOIs <- sort(unique(fwc_plant7$AreaOfInterestID))

# start pdf
pdf("output/plant_survey_detected_counts.pdf")

# cycle through AIOs
for(i in AOIs){
  
  # filter data
  dat_sub <- filter(fwc_plant7, AreaOfInterestID == i)
  
  # get name
  dat_name <- unique(dat_sub$AreaOfInterest)
  
  # summarize data
  dat_rich <- dat_sub %>%
    inner_join(key_all_pres) %>%
    group_by(AreaOfInterest, AreaOfInterestID, SurveyYear) %>%
    summarize(Richness = sum(IsDetected == "Yes"),
              .groups = "drop") %>%
    mutate(PrevRich = lag(Richness),
           NextRich = lead(Richness),
           MinSurvey = min(Richness),
           DevRich = case_when((PrevRich - Richness)/PrevRich >= 0.8 &
                                (NextRich - Richness)/NextRich >= 0.8 &
                                 Richness == MinSurvey ~ "outlier",
                               (PrevRich - Richness)/PrevRich >= 0.8 &
                                 (NextRich - Richness)/NextRich >= 0.8 ~ "deviation",
                              TRUE ~ "consistent"))
  
  dat_acre <- dat_sub %>%
    inner_join(key_exotic_acre2) %>%
    filter(TaxonName %in% c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"))
  
  # plot
  print(ggplot(dat_rich, aes(x = SurveyYear, y = Richness)) +
          geom_point(aes(color = DevRich)) +
          geom_line() +
          ggtitle(dat_name) +
          theme_bw())
  
  print(ggplot(dat_acre, aes(x = SurveyYear, y = SpeciesAcres, color = TaxonName)) +
          geom_point() +
          geom_line() +
          ggtitle(dat_name) +
          theme_bw())
  
  # save outliers
  dat_outliers = filter(dat_rich, DevRich != "consistent")
  
  if(nrow(dat_outliers >= 1)){
    
    outliers <- dat_outliers %>%
      select(AreaOfInterest, AreaOfInterestID, SurveyYear, DevRich) %>%
      full_join(outliers)
    
  }
  
}

dev.off()

# indicate outliers
fwc_plant8 <- fwc_plant7 %>%
  left_join(outliers) %>%
  rename(Outlier = DevRich) %>%
  mutate(Outlier = if_else(is.na(Outlier), 0, 1))


#### Permanent ID/AOI ####

# one permID per AOI?
fwc_plant8 %>%
  group_by(AreaOfInterestID) %>%
  summarize(nPerm = n_distinct(PermanentID)) %>%
  ungroup() %>%
  filter(nPerm > 1)
# yes

# are AOIs for each Perm ID consistent?
fwc_plant8 %>%
  group_by(PermanentID) %>%
  summarize(nAOI = n_distinct(AreaOfInterest),
            nAOI_ID = n_distinct(AreaOfInterestID)) %>%
  ungroup() %>%
  filter(nAOI_ID > 1) %>% # select lakes with multiple AOIs
  inner_join(fwc_plant8 %>%
               group_by(PermanentID, SurveyYear) %>%
               summarize(AOI = paste(sort(unique(AreaOfInterest)), collapse = ", ")) %>%
               ungroup() %>%
               group_by(PermanentID, AOI) %>%
               summarize(Years = paste(sort(unique(SurveyYear)), collapse = ", ")) %>%
               ungroup()) %>% # AOIs surveyed each year
  get_dupes(PermanentID) # select lakes with different AOIs over time

# dataset 1: greater temporal coverage
# remove AOI additions that create inconsistencies across time
fwc_plant8_time <- fwc_plant8 %>%
  filter(!(AreaOfInterest == "Silver Glen Springs" & PermanentID == "107881197") &
           !(AreaOfInterest == "Wauseon Bay" & PermanentID == "112029141") &
           !(AreaOfInterest == "Red Water, Lake" & PermanentID == "112047993") &
           !(AreaOfInterest == "Josephine Creek" & PermanentID == "112049879") &
           !(AreaOfInterest == "Hunt, Lake" & PermanentID == "167180956") &
           !(AreaOfInterest == "Tarpon, Lake Outfall Canal" & PermanentID == "68792760"))

# dataset2: greater spatial coverage
# remove years when both AOIs weren't sampled
fwc_plant8_space <- fwc_plant8 %>%
  filter(!(SurveyYear %in% c(1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 2017) & PermanentID == "107881197") &
           !(SurveyYear %in% c(1983, 1984, 1985) & PermanentID == "112029141") &
           !(SurveyYear %in% c(1983, 1984, 1986, 1988, 1989) & PermanentID == "112047993") &
           !(SurveyYear %in% c(1983, 1984, 1985, 1986, 1987, 1988) & PermanentID == "112049879") &
           !(SurveyYear %in% c(1983, 1984, 2000) & PermanentID == "167180956") &
           !(SurveyYear %in% c(2018, 2019) & PermanentID == "68792760"))

# reran summarizing function above with _time and _space datasets
# should return no rows


#### outputs ####
write_csv(fwc_plant8, "intermediate-data/FWC_plant_formatted.csv")
write_csv(fwc_plant8_time, "intermediate-data/FWC_plant_formatted_temporal_coverage.csv")
write_csv(fwc_plant8_space, "intermediate-data/FWC_plant_formatted_spatial_coverage.csv")
