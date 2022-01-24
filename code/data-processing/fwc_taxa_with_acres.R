#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv",
                      col_types = list(JoinNotes = col_character(),
                                       PermanentID = col_character()))

#### edit data ####

# summarize by taxa
summ_fwc <- plant_fwc %>%
  filter(!is.na(SpeciesAcres)) %>% 
  group_by(Origin, TaxonName) %>% 
  count() %>%
  ungroup() %>%
  rename(Surveys = "n") %>%
  arrange(Origin, desc(Surveys))


#### output ####
write_csv(summ_fwc, "intermediate-data/FWC_taxa_with_acres_summary.csv")
