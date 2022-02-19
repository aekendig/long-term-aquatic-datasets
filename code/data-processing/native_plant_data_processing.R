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
taxa_acres <- read_csv("intermediate-data/fwc_taxa_with_acres_summary.csv")
common_fwc <- read_csv("intermediate-data/FWC_common_native_plants_formatted.csv")

# load scripts
source("code/generic-functions/proportion_transformations.R")
source("code/data-processing/plant_abundance_formatting.R")
# this loads packages, surveyor function, Lake O materials, remove dups function


#### edit data ####

# choose taxa with > 1000 records 
nat_taxa <- taxa_acres %>%
  filter(TaxonName %in% unique(common_fwc$TaxonName) & 
           Surveys > 1000 & 
           str_detect(TaxonName, "spp.") == F &
           str_detect(TaxonName, "/") == F) %>%
  select(TaxonName) %>%
  mutate(CommonName = case_when(TaxonName == "Pontederia cordata" ~ "pickerelweed",
                                TaxonName == "Panicum hemitomon" ~ "maidencane",
                                TaxonName == "Nuphar advena" ~ "yellow pond-lily",
                                TaxonName == "Sagittaria lancifolia" ~ "bulltongue arrowhead",
                                TaxonName == "Cephalanthus occidentalis" ~ "common buttonbush",
                                TaxonName == "Nymphaea odorata" ~ "American white waterlily",
                                TaxonName == "Cladium jamaicense" ~ "Jamaica swamp sawgrass",
                                TaxonName == "Paspalidium geminatum" ~ "Egyptian panicgrass",
                                TaxonName == "Vallisneria americana" ~ "American eelgrass",
                                TaxonName == "Sacciolepis striata" ~ "American cupscale",
                                TaxonName == "Ceratophyllum demersum" ~ "coon's tail",
                                TaxonName == "Najas guadalupensis" ~ "southern waternymph",
                                TaxonName == "Spartina bakeri" ~ "sand cordgrass",
                                TaxonName == "Utricularia foliosa" ~ "leafy bladderwort",
                                TaxonName == "Nymphoides aquatica" ~ "big floatingheart",
                                TaxonName == "Luziola fluitans" ~ "southern watergrass",
                                TaxonName == "Juncus effusus" ~ "common rush",
                                TaxonName == "Bacopa caroliniana" ~ "blue waterhyssop",
                                TaxonName == "Utricularia gibba" ~ "humped bladderwort"))

# visualize raw abundances
pdf("output/native_plant_raw_abundance_time_series.pdf", width = 20, height = 16)
plant_fwc %>%
  filter(SpeciesName %in% nat_taxa$TaxonName) %>%
  ggplot(aes(x = SurveyYear, y = SpeciesAcres/WaterbodyAcres, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ TaxonName,
             scales = "free") +
  theme(legend.position = "none") + 
  labs(x = "Year", y = "Proportion area covered")

plant_fwc %>%
  filter(SpeciesName %in% nat_taxa$TaxonName) %>%
  mutate(SpeciesAcres = ifelse(SpeciesAcres == 0, NA_real_, SpeciesAcres)) %>%
  filter(!is.na(SpeciesAcres)) %>%
  ggplot(aes(x = SurveyYear)) +
  geom_bar() +
  facet_wrap(~ TaxonName, scales = "free") +
  labs(x = "Year", y = "Number of surveys with cover > 0")
dev.off()
# missing every other year 1985-1993 + 1994
# missing most years after 1995

# surveys
nat_surv <- plant_fwc %>%
  select(AreaOfInterestID) %>%
  unique() %>%
  expand_grid(SurveyYear = seq(min(plant_fwc$SurveyYear), max(plant_fwc$SurveyYear))) %>%
  left_join(plant_fwc %>%
              filter(Origin == "Native") %>%
              select(AreaOfInterestID, SurveyYear, TaxonName, SpeciesAcres, IsDetected)) %>%
  mutate(Measured = if_else(SpeciesAcres > 0, 1, 0),
         Measured = replace_na(Measured, 0)) %>%
  group_by(AreaOfInterestID, SurveyYear) %>%
  summarize(Measured = sum(Measured)) %>%
  ungroup() %>%
  mutate(Surveyed = if_else(Measured > 0, 1, 0)) %>%
  select(-Measured)

# check Tohopekaliga survey in 2017 (seems incomplete in other data exploration)
filter(nat_surv, AreaOfInterestID == 476 & SurveyYear == 2017)

# check that patterns look right
nat_surv %>%
  filter(Surveyed == 1) %>%
  ggplot(aes(x = SurveyYear)) +
  geom_bar()

# modify data
nat_fwc <- plant_fwc %>%
  plant_abun_format(nat_taxa, nat_surv) 

# check patterns
pdf("output/native_plant_processed_abundance_time_series.pdf", width = 20, height = 16)
nat_fwc %>%
  filter(!is.na(PropCovered)) %>%
  ggplot(aes(x = GSYear, y = PropCovered, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ TaxonName, scales = "free") +
  theme(legend.position = "none") + 
  labs(x = "Year", y = "Proportion area covered")

nat_fwc %>%
  filter(!is.na(PropCovered)) %>%
  ggplot(aes(x = GSYear)) +
  geom_bar() +
  facet_wrap(~ TaxonName, scales = "free") +
  labs(x = "Year", y = "Number of surveys with cover data")
dev.off()

# remove missing years
# remove non-surveyed years
nat_fwc2 <- nat_fwc %>%
  filter(!is.na(PropCovered))

#### outputs ####
write_csv(nat_fwc2, "intermediate-data/FWC_native_plant_abundance_formatted.csv")