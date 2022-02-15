#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor) # used get_dupes when developing methods

# import data
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv",
                      col_types = list(JoinNotes = col_character(),
                                       PermanentID = col_character()))
plant_fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv",
                       col_types = list(Depth_ft = col_double(),
                                        PermanentID = col_character(),
                                        YearF = col_character()))
taxa_acres <- read_csv("intermediate-data/fwc_taxa_with_acres_summary.csv")

# load scripts
source("code/generic-functions/proportion_transformations.R")
source("code/data-processing/plant_abundance_formatting.R")
# this loads packages, surveyor function, Lake O materials, remove dups function
source("code/data-processing/plant_frequency_formatting.R")


#### edit data ####

# multiple AOIs per permanent ID
plant_fwc %>%
  group_by(PermanentID) %>%
  summarize(Waterbodies = n_distinct(AreaOfInterestID),
            AOI = paste(unique(AreaOfInterest), collapse = "/"),
            ID = paste(unique(AreaOfInterestID), collapse = ", ")) %>%
  ungroup() %>%
  filter(Waterbodies > 1) %>%
  left_join(plant_fwri %>%
              select(PermanentID) %>%
              unique() %>%
              mutate(FWRI = 1))
# ID's to omit if comparing total area (probably not sampled by FWRI):
# 402, 469, 42, 220, 436
# Harris (177) and Little Harris (244) are sampled as one lake in FWRI
# Red Water (365) and Little Red Water (250) are sampled as one lake in FWRI

plant_fwri %>%
  group_by(PermanentID) %>%
  summarize(Waterbodies = n_distinct(AOI),
            AOI = paste(unique(AOI), collapse = "/")) %>%
  ungroup() %>%
  filter(Waterbodies > 1)

# choose taxa with > 1000 records 
taxa_acres %>%
  filter(Origin == "Exotic" & Surveys > 1000)
# Colocasia esculenta (Wild taro, Code: ELEA) is invasive, but commonly confused with native elephant ear (grouped together in FWRI)
# Cyperus blepharoleptos has synonyms: Oxycaryum cubense (syn. Scirpus cubensis) in FWRI
# FWRI doesn't distinguish between Salvinia except for Salvinia molesta (WAFE is Salvinia sp.)
# From FWC list, did not include:
# Canna spp. (too vague)
# Schinus terebinthifolius (Brazilian peppertree, terrestrial)
# Sapium sebiferum (Chinese tallow tree, terrestrial)
# Melaleuca quinquenervia (Melaleuca, terrestrial/aquatic, not in FWRI)

# specify taxa of interest
inv_taxa <- tibble(TaxonName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes", "Panicum repens", "Colocasia esculenta", "Urochloa mutica", "Alternanthera philoxeroides", "Cyperus blepharoleptos", "Salvinia minima"),
                   CommonName = c("Hydrilla", "Water lettuce", "Water hyacinth", "Torpedograss", "Wild taro", "Para grass", "Alligator weed", "Cuban bulrush", "Water fern"),
                   Code = c("HYDR", "WALE", "WAHY", "TORP", "ELEA", "PAGR", "ALWE", "BUSE", "WAFE"))

# modify data
inv_fwc <- plant_fwc %>%
  filter(!(WaterbodyName == "Tohopekaliga, Lake" & SurveyYear == 2017)) %>% # survey seems incomplete
  plant_abun_format(inv_taxa) # warnings from min/max functions, replaced with NA
inv_fwri <- plant_freq_format(plant_fwri, inv_taxa)

# hydrilla look-alikes
# Elodea canadensis is also one, but not in datasets
hyd_looks_fwc <- plant_fwc %>%
  filter(!(WaterbodyName == "Tohopekaliga, Lake" & SurveyYear == 2017)) %>% # survey seems incomplete
  plant_abun_format(tibble(TaxonName = c("Egeria densa", "Najas guadalupensis"),
                                          CommonName = c("Edensa", "Nguad"))) %>%
  mutate(SpeciesPresent = ifelse(SpeciesAcres > 0, 1, 0)) %>%
  select(PermanentID, GSYear, CommonName, SpeciesPresent) %>%
  pivot_wider(names_from = CommonName,
              values_from = SpeciesPresent,
              names_glue = "{CommonName}_Present")

hyd_looks_fwri <- plant_freq_format(plant_fwri,
                                   tibble(Code = "SONA",
                                          CommonName = "Nguad")) %>%
  select(PermanentID, GSYear, SpeciesPresent) %>%
  rename("Nguad_Present" = "SpeciesPresent")

# remove missing years
# add hydrilla look-alikes
inv_fwc2 <- inv_fwc %>%
  filter(!is.na(PropCovered)) %>%
  left_join(hyd_looks_fwc)

inv_fwri2 <- inv_fwri %>%
  left_join(hyd_looks_fwri)


#### outputs ####
write_csv(inv_fwc2, "intermediate-data/FWC_invasive_plant_formatted.csv")
write_csv(inv_fwri2, "intermediate-data/FWRI_invasive_plant_formatted.csv")
