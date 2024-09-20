#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom.mixed)
library(glmmTMB)
library(DHARMa)
library(GGally)

# import data
rich_dat <- read_csv("intermediate-data/FWC_plant_management_richness_analysis_formatted.csv")


#### examine richness data ####

# correlations among explanatory variables
rich_dat %>%
  distinct(AreaOfInterestID, HydrPAC, FloatPAC, HydrTrtFreq, HydrTrtArea, FloatTrtFreq, FloatTrtArea,
           OtherTrtFreq, OtherTrtArea) %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# correlations >= 0.4 and sig
# OtherTrtFreq and FloatTrtFreq: 0.6
# OtherTrtArea and FloatTrtArea: 0.4


#### native richness model ####

# response variable
ggplot(rich_dat, aes(x = NativeRichness)) +
  geom_density()

# fit model
nat_mod1 <- glmmTMB(NativeRichness ~ Time * (),
                    data = rich_dat,
                    family = poisson)