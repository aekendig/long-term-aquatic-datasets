#### info ####

# goal: evaluate effects of herbicides on hydrilla


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)
library(zoo)
library(ggfortify)

# stan settings
source("code/stan_settings.R")

# import data
ctrl_old <- read_csv("intermediate-data/FWC_control_old_formatted.csv")
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
plant_lw <- read_csv("intermediate-data/LW_plant_formatted.csv")

# assumptions
MinHerbLag = 14 # herbicides applied within x days haven't had an effect yet
MaxHerbLag = 365 * 2 # effects of herbicide applied longer than x days ago can't be detected

# figure settings
def_theme <- theme_bw() +
  theme(axis.text = element_text(size = 8, color="black"),
        axis.title = element_text(size = 10, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.text = element_text(size = 10, color="black"),
        strip.background = element_blank())


#### edit data ####

# select data associated with Pistia
# convert all to hectares
ctrl_old_pis <- ctrl_old %>%
  filter(str_detect(Species, "Pistia") == T) %>%
  mutate(AreaTreated_ha= TotalAcres * 0.405,
         Area_ha = ShapeArea * 100)

ctrl_new_pis <- ctrl_new %>%
  filter(str_detect(Species, "Pistia") == T) %>%
  mutate(AreaTreated_ha = TotalAcres * 0.405,
         Area_ha = ShapeArea * 100)

plant_fwc_pis <- plant_fwc %>%
  filter(SpeciesName == "Pistia stratiotes") %>%
  mutate(AreaCovered_ha = SpeciesAcres * 0.405,
         Area_ha = ShapeArea * 100,
         AreaCovered_ha = case_when(AreaCovered_ha > Area_ha ~ Area_ha,
                                    TRUE ~ AreaCovered_ha),
         SpeciesFrequency_ha = AreaCovered_ha / Area_ha,
         log_AreaCovered_ha = log(AreaCovered_ha + 0.001))

plant_lw_pis <- plant_lw %>%
  filter(GenusSpecies == "Pistia stratiotes") %>%
  mutate(Area_ha = ShapeArea * 100)


#### plant_fwc_pis exploratory figures ####

pdf("output/pistia_abundance_fwc_over_time.pdf")
for(i in sort(unique(plant_fwc_pis$AreaOfInterestID))){
  subdat <- filter(plant_fwc_pis, AreaOfInterestID == i)
  lake <- unique(subdat$AreaOfInterest)
  print(ggplot(subdat, aes(SurveyDate, AreaCovered_ha)) +
    geom_point() +
    geom_line() +
    ggtitle(lake) +
    def_theme)
}
dev.off()
