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

# specify taxa of interest
inv_taxa <- tibble(TaxonName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"),
                   CommonName = c("Hydrilla", "Water lettuce", "Water hyacinth"),
                   Code = c("HYDR", "WALE", "WAHY"))

# modify data
inv_fwc <- plant_abun_format(plant_fwc, inv_taxa)
inv_fwri <- plant_freq_format(plant_fwri, inv_taxa)

# hydrilla look-alikes
# Elodea canadensis is also one, but not in datasets
hyd_looks_fwc <- plant_abun_format(plant_fwc,
                                   tibble(TaxonName = c("Egeria densa", "Najas guadalupensis"),
                                          CommonName = c("Edensa", "Nguad"))) %>%
  select(PermanentID, GSYear, CommonName, SpeciesPresent) %>%
  pivot_wider(names_from = CommonName,
              values_from = SpeciesPresent,
              names_glue = "{CommonName}_Present")

hyd_looks_fwri <- plant_freq_format(plant_fwri,
                                   tibble(Code = "SONA",
                                          CommonName = "Nguad")) %>%
  select(PermanentID, GSYear, SpeciesPresent) %>%
  rename("Nguad_Present" = "SpeciesPresent")

# add hydrilla look-alikes
inv_fwc2 <- inv_fwc %>%
  left_join(hyd_looks_fwc)

inv_fwri2 <- inv_fwri %>%
  left_join(hyd_looks_fwri)


#### outputs ####
write_csv(inv_fwc2, "intermediate-data/FWC_invasive_plant_formatted.csv")
write_csv(inv_fwri2, "intermediate-data/FWRI_invasive_plant_formatted.csv")
