
# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
herb_type <- read_csv("intermediate-data/herbicide_types.csv")

# non-herbicide methods
non_herb <- c("Mechanical Harvester", 
              "Snagging (tree removal)", 
              "Aquatic Dye (for shading)", 
              "Grass Carp", "Hand Removal", 
              "Mechanical (Other)", 
              "Mechanical Shredder", 
              "Prescribed Fire")

# taxa of interest
herb_taxa <- tibble(Species = c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)"))

# edit data
ai_dat <- ctrl_new %>% filter(Species %in% herb_taxa$Species &
                      TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>% 
  left_join(herb_type) %>% 
  select(HerbicideType, RateType, Unit, ActiveIngredient, ControlMethod, Herbicide) %>% 
  unique() %>% 
  arrange(HerbicideType, RateType, Unit, ActiveIngredient, ControlMethod, Herbicide) %>% 
  data.frame()

# save for manual editing
write_csv(ai_dat, "intermediate-data/active_ingredient_conversion.csv")
# edit "manual" version