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
common_fwc <- read_csv("intermediate-data/FWC_common_native_plants_formatted.csv")

# load scripts
source("code/generic-functions/proportion_transformations.R")
source("code/data-processing/plant_abundance_formatting.R")
# this loads packages, surveyor function, Lake O materials, remove dups function


#### edit data ####

# summarize by taxon
surv_num <- plant_fwc %>%
  filter(Origin == "Native" & !is.na(SpeciesAcres) & 
           SpeciesAcres > 0) %>%
  group_by(TaxonName) %>%
  summarize(surveys = n()) %>%
  ungroup()

# visualize survey numbers
ggplot(surv_num, aes(x = surveys)) +
  geom_histogram(bins = 100)

# choose taxa with > 200 records 
surv_comm <- surv_num %>%
  filter(surveys >= 1000) %>%
  select(TaxonName) %>%
  filter(str_detect(TaxonName, "spp.") == F &# remove genus-level
           str_detect(TaxonName, "/") == F &
           TaxonName != "Filamentous algae")

# add common names
nat_taxa <- surv_comm %>%
  mutate(CommonName = case_when(TaxonName =="Cephalanthus occidentalis" ~ "common buttonbush",
                                TaxonName =="Cladium jamaicense" ~ "Jamaica swamp sawgrass",
                                TaxonName =="Nuphar advena" ~ "yellow pond-lily",
                                TaxonName =="Nymphaea odorata" ~ "American white waterlily",
                                TaxonName =="Panicum hemitomon" ~ "maidencane",
                                TaxonName =="Paspalidium geminatum" ~ "Egyptian panicgrass",
                                TaxonName =="Pontederia cordata" ~ "pickerelweed",
                                TaxonName =="Sagittaria lancifolia" ~ "bulltongue arrowhead"))

# visualize raw abundances
pdf("output/native_plant_raw_abundance_time_series.pdf", width = 20, height = 16)
plant_fwc %>%
  filter(TaxonName %in% nat_taxa$TaxonName & !is.na(SpeciesAcres)) %>%
  ggplot(aes(x = GSYear, y = (SpeciesAcres * 0.405)/(ShapeArea * 100), color = PermanentID)) +
  geom_line(size = 0.5, alpha = 0.25) +
  geom_point(size = 0.5, alpha = 0.25) + 
  facet_wrap(~ TaxonName,
             scales = "free_y") +
  theme(legend.position = "none") + 
  labs(x = "Year", y = "Proportion area covered")

plant_fwc %>%
  filter(SpeciesName %in% nat_taxa$TaxonName) %>%
  mutate(SpeciesAcres = ifelse(SpeciesAcres == 0, NA_real_, SpeciesAcres)) %>%
  filter(!is.na(SpeciesAcres)) %>%
  ggplot(aes(x = GSYear)) +
  geom_bar() +
  facet_wrap(~ TaxonName, scales = "free_y") +
  labs(x = "Year", y = "Number of surveys with cover > 0")
dev.off()

# looked at taxa with highest survey numbers after 1998
# and they are included in this list

# modify data
nat_fwc <- plant_fwc %>%
  filter(!(AreaOfInterestID == 476 & SurveyYear == 2017)) %>% # incomplete survey
  plant_abun_format(nat_taxa) 

# check patterns
pdf("output/native_plant_processed_abundance_time_series.pdf", width = 20, height = 16)
nat_fwc %>%
  filter(!is.na(PropCovered)) %>%
  ggplot(aes(x = GSYear, y = PropCovered, color = PermanentID)) +
  geom_line(size = 0.5, alpha = 0.25) +
  geom_point(size = 0.5, alpha = 0.25) + 
  facet_wrap(~ TaxonName, scales = "free") +
  theme(legend.position = "none") + 
  labs(x = "Year", y = "Proportion area covered")

nat_fwc %>%
  mutate(PropCovered = ifelse(PropCovered == 0, NA_real_, PropCovered)) %>%
  filter(!is.na(PropCovered)) %>%
  ggplot(aes(x = GSYear)) +
  geom_bar() +
  facet_wrap(~ TaxonName, scales = "free") +
  labs(x = "Year", y = "Number of surveys with cover > 0")
dev.off()

# remove missing years
# remove non-surveyed years
nat_fwc2 <- nat_fwc %>%
  filter(!is.na(PropCovered))

#### outputs ####
write_csv(nat_fwc2, "intermediate-data/FWC_native_plant_abundance_formatted.csv")
