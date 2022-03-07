#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor) # used get_dupes when developing methods

# import data
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted_temporal_coverage.csv",
                      col_types = list(JoinNotes = col_character(),
                                       PermanentID = col_character()))
plant_fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv",
                       col_types = list(Depth_ft = col_double(),
                                        PermanentID = col_character(),
                                        YearF = col_character()))
major_concern <- read_csv("intermediate-data/fwc_invasive_taxa_major_concern.csv")

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
# ID's to omit if comparing total area (probably not sampled by FWRI): 42
# Harris (177) and Little Harris (244) are sampled as one lake in FWRI
# Red Water (365) and Little Red Water (250) are sampled as one lake in FWRI, but two in FWC

plant_fwri %>%
  group_by(PermanentID) %>%
  summarize(Waterbodies = n_distinct(AOI),
            AOI = paste(unique(AOI), collapse = "/")) %>%
  ungroup() %>%
  filter(Waterbodies > 1)

# choose high priority taxa 
major_concern %>%
  select(TaxonName) %>%
  unique()
# Colocasia esculenta (Wild taro, Code: ELEA) is invasive, but commonly confused with native elephant ear (grouped together in FWRI)
# Cyperus blepharoleptos has synonyms: Oxycaryum cubense (syn. Scirpus cubensis) in FWRI


# specify taxa of interest
inv_taxa <- tibble(TaxonName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes", "Panicum repens", "Urochloa mutica", "Cyperus blepharoleptos"),
                   CommonName = c("Hydrilla", "Water lettuce", "Water hyacinth", "Torpedograss", "Para grass", "Cuban bulrush"),
                   Code = c("HYDR", "WALE", "WAHY", "TORP", "PAGR", "BUSE"))

# visualize raw abundances
pdf("output/invasive_plant_raw_abundance_time_series.pdf", width = 18, height = 16)
plant_fwc %>%
  filter(TaxonName %in% inv_taxa$TaxonName & !is.na(SpeciesAcres)) %>%
  ggplot(aes(x = GSYear, y = (SpeciesAcres * 0.405)/(ShapeArea * 100), color = PermanentID)) +
  geom_line(size = 0.5, alpha = 0.25) +
  geom_point(size = 0.5, alpha = 0.25) +
  facet_wrap(~ TaxonName,
             scales = "free") +
  theme(legend.position = "none") + 
  labs(x = "Year", y = "Proportion area covered")

plant_fwc %>%
  filter(TaxonName %in% inv_taxa$TaxonName) %>%
  mutate(SpeciesAcres = ifelse(SpeciesAcres == 0, NA_real_, SpeciesAcres)) %>%
  filter(!is.na(SpeciesAcres)) %>%
  ggplot(aes(x = GSYear)) +
  geom_bar() +
  facet_wrap(~ TaxonName, scales = "free") +
  labs(x = "Year", y = "Number of surveys with cover > 0")
dev.off()

# modify data
inv_fwc <-  plant_fwc %>%
  filter(!(AreaOfInterestID == 476 & SurveyYear == 2017)) %>% # incomplete survey
  plant_abun_format(inv_taxa)
inv_fwri <- plant_freq_format(plant_fwri, inv_taxa)

# check patterns
pdf("output/invasive_plant_processed_abundance_time_series.pdf", width = 18, height = 16)
inv_fwc %>%
  filter(!is.na(PropCovered)) %>% # there are so many NA's that a lot of data don't show up
  ggplot(aes(x = GSYear, y = PropCovered, color = PermanentID)) +
  geom_line(size = 0.5, alpha = 0.25) +
  geom_point(size = 0.5, alpha = 0.25) + 
  facet_wrap(~ TaxonName, scales = "free") +
  theme(legend.position = "none") + 
  labs(x = "Year", y = "Proportion area covered")

inv_fwc %>%
  mutate(PropCovered = ifelse(PropCovered == 0, NA_real_, PropCovered)) %>%
  filter(!is.na(PropCovered)) %>%
  ggplot(aes(x = GSYear)) +
  geom_bar() +
  facet_wrap(~ TaxonName, scales = "free") +
  labs(x = "Year", y = "Number of surveys with cover data")
dev.off()

# hydrilla look-alikes
plant_fwc %>%
  filter(TaxonName %in% c("Egeria densa", "Najas guadalupensis", "Elodea canadensis") &
           !is.na(SpeciesAcres)) %>%
  ggplot(aes(x = SurveyYear, y = SpeciesAcres/WaterbodyAcres, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ TaxonName,
             scales = "free") +
  theme(legend.position = "none") + 
  labs(x = "Year", y = "Proportion area covered")
# no Elodea canadensis
# limited Egeria dens data
# Najas guadalupensis data ends around 2005

hyd_looks_fwri <- plant_freq_format(plant_fwri,
                                   tibble(Code = "SONA",
                                          CommonName = "Nguad")) %>%
  select(PermanentID, GSYear, SpeciesPresent) %>%
  rename("Nguad_Present" = "SpeciesPresent")

# remove missing years
# add hydrilla look-alikes
inv_fwc2 <- inv_fwc %>%
  filter(!is.na(PropCovered)) # remove rows without acreage surveyed and rows added in to calculate growth

inv_fwri2 <- inv_fwri %>%
  left_join(hyd_looks_fwri)


#### outputs ####
write_csv(inv_fwc2, "intermediate-data/FWC_invasive_plant_formatted.csv")
write_csv(inv_fwri2, "intermediate-data/FWRI_invasive_plant_formatted.csv")
