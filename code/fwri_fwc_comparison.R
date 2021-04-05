#### info ####

# goal: see how FWRI compares to FWC


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv")
fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")


