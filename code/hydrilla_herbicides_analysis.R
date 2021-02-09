#### info ####

# goal: evaluate effects of herbicides on hydrilla


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)

# import data
ctrl_old <- read_csv("intermediate-data/FWC_control_old_formatted.csv")
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
plant_lw <- read_csv("intermediate-data/LW_plant_formatted.csv")


#### edit data ####

# select data associated with Hydrilla
ctrl_old_hyd <- ctrl_old %>%
  filter(Species == "Hydrilla verticillata")

ctrl_new_hyd <- ctrl_new %>%
  filter(Species == "Hydrilla verticillata")

plant_fwc_hyd <- plant_fwc %>%
  filter(SpeciesName == "Hydrilla verticillata")

plant_lw_hyd <- plant_lw %>%
  filter(GenusSpecies == "Hydrilla verticillata")