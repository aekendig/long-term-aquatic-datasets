#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(data.table)
library(tidyverse)
library(janitor)
library(lubridate)
library(taxize)

# load data
gis <- read_csv("intermediate-data/gis_fwc_lakewatch_fwri.csv",
                col_types = list(wkt_geom = col_character(),
                                 AOI = col_character(),
                                 Lake = col_character()))
ctrl_old <- read_csv("original-data/PrePMARS_IPMData.csv")
ctrl <- read_csv("original-data/FWC_Herbicide_Treatments_mid2010-Oct2020.csv")
herb_type <- read_csv("intermediate-data/herbicide_types_all.csv") # I added a few to Candice's list, she only did ones for the three main invasive plants, her file is herbicide_types_all
plants <- read_csv("intermediate-data/FWC_plant_formatted.csv")


#### old control data ####

# years
ctrl_old %>% select(FY, Year) %>% unique()
# fiscal year starts July 1
# all treatments are 2nd year of FY range, not sure how accurate this is

# rename columns
# add permanent ID based on AreaOfInterestID
ctrl_old2 <- ctrl_old %>%
  filter(TotalAcres > 0) %>%
  rename("Longitude_FWC" = "longitude",
         "Latitude_FWC" = "latitude",
         "SpeciesOrig" = "Species_orig",
         "TotalContFWC" = "Total_cont_fwc",
         "County_FWC" = "County") %>%
  left_join(gis %>%
              filter(CoordSource %in% c("FWC", "FWC_2")) %>%
              select(AreaOfInterestID, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique()) %>%
  mutate(County_FWC = toupper(County_FWC),
         ControlMethod = "unknown",
         TreatmentID = paste("old", Year, substr(Species, 1, 1), TotalAcres, sep = "_"),
         CtrlSet = "old",
         MethodHerbicide = "unknown") %>%
  rename(TreatmentYear = Year)

# duplicate rows?
ctrl_old2 %>%
  get_dupes()
# none

# missing ID?
ctrl_old2 %>%
  filter(is.na(PermanentID)) %>%
  select(AreaOfInterestID, AreaOfInterest) %>%
  unique() %>%
  data.frame()
# 72
# not lakes

# consistent waterbody name?
ctrl_old2 %>%
  distinct(AreaOfInterestID, AreaOfInterest) %>%
  get_dupes(AreaOfInterestID)
# yes

# remove missing perm ID
# remove duplicate rows
# select relevant columns
ctrl_old3 <- ctrl_old2 %>%
  filter(!is.na(PermanentID))

# one permID per AOI?
ctrl_old3 %>%
  group_by(AreaOfInterestID) %>%
  summarize(nPerm = n_distinct(PermanentID)) %>%
  ungroup() %>%
  filter(nPerm > 1)
# yes

# AOIs per permID?
ctrl_old3 %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest) %>%
  get_dupes(PermanentID)
# these are the ones that cause temporal inconsistencies


#### new control data ####

# non-herbicide methods (from herbicide_initial_visualizations)
non_herb <- c("Mechanical Harvester", 
              "Snagging (tree removal)", 
              "Aquatic Dye (for shading)", 
              "Grass Carp", "Hand Removal", 
              "Mechanical (Other)", 
              "Mechanical Shredder", 
              "Prescribed Fire")

# format date
# rename columns
# add permanent ID based on AreaOfInterestID
# add herbicide features
ctrl2 <- ctrl %>%
  filter(TotalAcres > 0) %>%
  mutate(BeginDate = as.Date(BeginDate, "%m/%d/%y"),
         County = toupper(County),
         MethodHerbicide = case_when(ControlMethod %in% non_herb ~ "no",
                               is.na(ControlMethod) ~ "unknown",
                               TRUE ~ "yes"),
         CtrlSet = "new") %>%
  rename("TreatmentDate" = "BeginDate",
         "TreatmentYear" = "year",
         "TreatmentMonth" = "month",
         "Species" = "species",
         "Herbicide" = "herbicide",
         "TotalHerbicideUsed" = "totalherbicideused",
         "Longitude_FWC" = "longitude",
         "Latitude_FWC" = "latitude",
         "County_FWC" = "County")  %>%
  left_join(gis %>%
              filter(CoordSource %in% c("FWC", "FWC_2")) %>%
              select(AreaOfInterestID, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique()) %>%
  left_join(herb_type)

# duplicate rows?
ctrl2 %>%
  get_dupes()
# 39 rows

# missing ID?
ctrl2 %>%
  filter(is.na(PermanentID)) %>%
  select(AreaOfInterestID, AreaOfInterest) %>%
  unique() %>%
  data.frame()
# 97
# not lakes

# consistent waterbody name?
ctrl2 %>%
  distinct(AreaOfInterestID, AreaOfInterest) %>%
  get_dupes(AreaOfInterestID)
# yes

# remove missing perm ID
# remove duplicate rows
ctrl3 <- ctrl2 %>%
  filter(!is.na(PermanentID)) %>%
  unique()
# dataset includes duplication of treatment ID (single event) because of different methods used within event

# one permID per AOI?
ctrl3 %>%
  group_by(AreaOfInterestID) %>%
  summarize(nPerm = n_distinct(PermanentID)) %>%
  ungroup() %>%
  filter(nPerm > 1)
# yes

# AOIs per permID?
ctrl3 %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest) %>%
  get_dupes(PermanentID)
# these cause temporal inconsistencies


#### ctrl outputs ####

write_csv(ctrl_old3, "intermediate-data/FWC_management_old_formatted.csv")
write_csv(ctrl3, "intermediate-data/FWC_management_new_formatted.csv")


#### combine datasets ####

# duplicates within old dataset
ctrl_old3 %>%
  get_dupes(AreaOfInterestID, TreatmentYear, Species, TotalAcres) %>%
  data.frame()
# two different locations within same lake

# duplicates within new dataset
ctrl3 %>%
  get_dupes(AreaOfInterestID, TreatmentYear, Species, 
            TotalAcres, ControlMethod, TreatmentID)
# 253, look closer at a few
ctrl3 %>%
  get_dupes(AreaOfInterestID, TreatmentYear, Species, 
            TotalAcres, ControlMethod, TreatmentID) %>%
  filter(FormID == 14501) %>%
  data.frame()
# different total herbicide and amt active ingredient
ctrl3 %>%
  get_dupes(AreaOfInterestID, TreatmentYear, Species, 
            TotalAcres, ControlMethod, TreatmentID) %>%
  filter(FormID == 27276) %>%
  data.frame()
# different herbicide

# same waterbody names?
ctrl_old3 %>%
  distinct(AreaOfInterestID, AreaOfInterest) %>%
  full_join(ctrl3 %>%
              distinct(AreaOfInterestID, AreaOfInterest)) %>%
  get_dupes(AreaOfInterestID)
# change watermellon to one l like plant

# management dates in new dataset
ggplot(ctrl3, aes(x = TreatmentMonth)) +
  geom_histogram(binwidth = 1)

# combine datasets
# add fake date for old ctrl
# correct mispelled AOI
# identify focal invasive species
ctrl4 <- ctrl3 %>%
  mutate(TreatmentID = as.character(TreatmentID)) %>%
  full_join(ctrl_old3 %>%
              mutate(TreatmentDate = as.Date(paste0(TreatmentYear, "-01-01")))) %>%
  mutate(AreaOfInterest = if_else(AreaOfInterest == "Watermellon Pond" & AreaOfInterestID == 465,
                                  "Watermelon Pond",
                                  AreaOfInterest))
  
# check for duplication around 2010
# can use date from new dataset?
(ctrl4_2010_dups <- ctrl4 %>%
  filter(TreatmentYear == 2010) %>%
    group_by(AreaOfInterestID, Species, TotalAcres) %>%
    summarize(Sets = n_distinct(CtrlSet)) %>%
    ungroup() %>%
    filter(Sets > 1))
# 6 cases, two are floating plants
# floating plant acres are "1", so it could have been two different treatments


#### compare to plant surveys ####

# are AOIs for each Perm ID consistent over time?
ctrl4 %>%
  group_by(PermanentID) %>%
  summarize(nAOI = n_distinct(AreaOfInterest),
            nAOI_ID = n_distinct(AreaOfInterestID)) %>%
  ungroup() %>%
  filter(nAOI_ID > 1) %>% # select lakes with multiple AOIs
  inner_join(ctrl4 %>%
               group_by(PermanentID, TreatmentYear) %>%
               summarize(AOI = paste(sort(unique(AreaOfInterest)), collapse = ", ")) %>%
               ungroup() %>%
               group_by(PermanentID, AOI) %>%
               summarize(Years = paste(sort(unique(TreatmentYear)), collapse = ", ")) %>%
               ungroup()) %>% # AOIs treated each year
  get_dupes(PermanentID) # select lakes with different AOIs over time
# consistent with plants (fwc_plant_data_processing)
# check Little Santa Fe

# plant surveys
plant_surv <- plants %>%
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, SurveyYear, SurveyDate)

# look at Little Santa Fe
filter(plant_surv, PermanentID == "113638769") %>%
  group_by(AreaOfInterest) %>%
  summarize(Years = paste(sort(unique(SurveyYear)), collapse = ", ")) %>%
  ungroup() %>%
  data.frame()
# Little Santa Fe was only managed in one year, but it was surveyed with Santa Fe every year

# check that AOIs are consistent
ctrl4 %>%
  filter(PermanentID %in% plant_surv$PermanentID & 
           (!(AreaOfInterestID %in% plant_surv$AreaOfInterestID) |
              !(AreaOfInterest %in% plant_surv$AreaOfInterest))) %>%
  distinct(AreaOfInterest, AreaOfInterestID, PermanentID)
# DeLeon Springs St Park Wetlands

# look at DeLeon Springs
plant_surv %>%
  filter(PermanentID == "107776163") %>%
  distinct(AreaOfInterest, AreaOfInterestID, PermanentID)
ctrl4 %>%
  filter(PermanentID == "107776163") %>%
  distinct(AreaOfInterest, AreaOfInterestID, PermanentID)
# okay to remove the wetlands from the ctrl dataset
# different area than lake in both datasets

# check that waterbody names are consistent
ctrl4 %>%
  distinct(AreaOfInterestID, AreaOfInterest) %>%
  full_join(plant_surv %>%
              distinct(AreaOfInterestID, AreaOfInterest)) %>%
  get_dupes(AreaOfInterestID)
# yes

# remove wetlands that aren't in plant surveys
ctrl5 <- ctrl4 %>%
  filter(!(AreaOfInterest == "DeLeon Springs St Park Wetlands" & PermanentID == "107776163"))


#### resolve species names ####

# notes:
# AHRES = Aquatic Habitat Restoration & Enhancement Section

# list of taxa
taxa_list <- ctrl5 %>%
  distinct(Species) %>%
  mutate(Taxon = str_replace_all(Species, " spp.| spp| sp.| sp", ""), # for genus-level
         Taxon = str_replace_all(Taxon, ", sub| \\(exotic\\)| \\(other natives\\)| \\(other\\)|, natives|, emersed|, sub/floating", ""), # for origin/growth type
         Taxon = str_replace(Taxon, "\\/.*", ""), # anything after / removed
         Taxon = str_replace(Taxon, "sub\\.", "ssp."), # format sub-species
         words = str_count(Taxon, pattern = " ") + 1) %>%
  full_join(ctrl5 %>%
              distinct(Species) %>%
              filter(str_detect(Species, "\\/") & str_detect(Species, "spp") == F) %>%
              mutate(Taxon = paste(word(Species, 1, 1),
                                   str_replace(Species, ".*\\/", "")), # match genus with second species name
                     words = str_count(Taxon, pattern = " ") + 1)) %>%
  filter(words > 1 & str_detect(Taxon, "AHRES") == F & 
           !(Taxon %in% c("Mixed Grasses", "Trees (blocking navigation)",
                          "Floating Plants (Eichhornia and Pistia)"))) %>%
  arrange(Species)

# synonyms from ITIS
# taxa_syn <- synonyms(taxa_list$Taxon, db = "itis")
# manually accepted duplicate names
# chose accepted ones and checked on website
# date run: 10/26/24

# make into dataframe
# taxa_syn2 <- rbindlist(lapply(taxa_syn, as.data.table), use.names = T, fill = T, idcol = "Taxon") %>%
#   select(-V1) %>%
#   as_tibble()

# manually check missing species
# taxa_syn2 %>%
#   filter(is.na(sub_tsn)) %>%
#   select(Taxon)
# one is general algae
# Cyperus blepharoleptos not in ITIS database - Oxycarum is the accepted name; IFAS uses C. bleph
# Egeria najas not in ITIS database, neither of synonyms Elodea najas and Anacharis najas
# type-os: 
# Mormodica charantia = Momordica charantia - no synonyms


# save
# write_csv(taxa_syn2, "intermediate-data/FWC_management_species_synonyms.csv")

# import
taxa_syn2 <- read_csv("intermediate-data/FWC_management_species_synonyms.csv")

# identify synonyms to rename
(taxa_acc <- taxa_syn2 %>%
    filter(!is.na(acc_name) & Taxon != acc_name) %>%
    distinct(Taxon, acc_name) %>%
    add_row(Taxon = "Oxycaryum cubense", acc_name = "Cyperus blepharoleptos") %>%
    add_row(Taxon = "Mormodica charantia", acc_name = "Momordica charantia"))

# are accepted names used in list?
taxa_list %>%
  filter(Taxon %in% taxa_acc$acc_name)
# 2

# update synonyms
taxa_list2 <- taxa_list %>%
  left_join(taxa_acc %>%
              rename(TaxonName = acc_name)) %>%
  mutate(TaxonName = if_else(is.na(TaxonName), Species, TaxonName))

# look at taxa with multiple species
get_dupes(taxa_list2, Species)

# manually update taxa with multiple species and accepted names
taxa_list3 <- taxa_list2 %>%
  mutate(TaxonName = if_else(Species == "Schoenoplectus californicus/validus",
                             "Schoenoplectus californicus/tabernaemontani",
                             TaxonName)) %>%
  distinct(Species, TaxonName)

# add taxon name to data
ctrl6 <- ctrl5 %>%
  left_join(taxa_list3) %>%
  mutate(TaxonName = if_else(is.na(TaxonName), Species, TaxonName))
  

#### final formatting ####

# remove AOIs for temporal consistency
ctrl6_time <- ctrl6 %>%
  filter(!(AreaOfInterest == "Silver Glen Springs" & PermanentID == "107881197") &
           !(AreaOfInterest == "Wauseon Bay" & PermanentID == "112029141") &
           !(AreaOfInterest == "Red Water, Lake" & PermanentID == "112047993") &
           !(AreaOfInterest == "Josephine Creek" & PermanentID == "112049879") &
           !(AreaOfInterest == "Hunt, Lake" & PermanentID == "167180956") &
           !(AreaOfInterest == "Tarpon, Lake Outfall Canal" & PermanentID == "68792760"))

# save
write_csv(ctrl6, "intermediate-data/FWC_management_formatted.csv")
write_csv(ctrl6_time, "intermediate-data/FWC_management_formatted_temporal_coverage.csv")
