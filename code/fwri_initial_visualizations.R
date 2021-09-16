#### info ####

# goal: visualize FWRI data


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)

# import data
fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv",
                 col_types = list(Depth_ft = col_double()))


#### edit data ####

# single year lakes
single_year <- fwri %>%
  group_by(Lake, AOI) %>%
  summarise(Years = n_distinct(Year),
            Year = unique(Year)) %>%
  ungroup() %>%
  filter(Years < 2)

# remove problematic samples
fwri2 <- fwri %>%
  filter(!(AOI %in% c("Orange", "Eustis2") & !(Lake == "Kings Bay" & Year == 2015 & Month == 2)))
# Orange has 2 sets of surveys: one for full lake and one for open water
# Eustis2 is probably EastToho in 2019 based on coordiantes, but that lake has a survey that year
# Kings Bay was sampled twice in 2015, but it is usually sampled in Oct, so Feb sample is removed

# remove single year lakes when we're interested in change over time
fwri3 <- fwri2 %>%
  anti_join(single_year)


#### subset data functions ####

# select a species and fill in 0 abundance
spp_fun <- function(code, dat){
  
  dat_out <- dat %>%
    select(AOI, Lake, Year, Site) %>% # start with all surveys (no species means abundance = 0)
    unique() %>%
    left_join(dat %>%  # add species information
                filter(Code == code) %>%
                select(AOI, Lake, Year, Site, Abundance)) %>%
    mutate(Abundance = replace_na(Abundance, 0)) # abundance 0 when it wasn't in a survey

  return(dat_out)
}

# proportion of sites occupied over time
prop_fun <- function(dat){
  
  dat_out <- dat %>%
    group_by(AOI, Lake, Year) %>%
    summarise(Sites = n(),
              SitesOcc1 = sum(Abundance > 0),
              SitesOcc2 = sum(Abundance > 1),
              SitesOcc3 = sum(Abundance > 2)) %>%
    ungroup() %>%
    pivot_longer(cols = starts_with("SitesOcc"),
                 names_to = "MinAbundance",
                 names_prefix = "SitesOcc",
                 values_to = "SitesOcc") %>%
    mutate(PropCovered = SitesOcc / Sites)
  
  return(dat_out)
  
}

# change in site over time
site_fun <- function(dat){
  
  dat_out <- dat %>%
    group_by(AOI, Lake, Site) %>%
    arrange(Year) %>%
    mutate(MinYear = min(Year),
           YearChange = Year - lag(Year),
           AbundanceChange = Abundance - lag(Abundance)) %>%
    ungroup() %>%
    filter(Year != MinYear)
  
  return(dat_out)
  
}


#### time series ####

# hydrilla data
fwri_hydr <- spp_fun("HYDR", fwri3)
fwri_hydr_prop <- prop_fun(fwri_hydr)
# fwri_hydr_site <- site_fun(fwri_hydr)
# realized sites are inconsistent over time for some lakes
# starting resolving this is fwri_data_formatting
# need to finish there

ggplot(fwri_hydr_prop, aes(x = Year, y = PropCovered, color = AOI)) +
  geom_point(show.legend = F) +
  geom_line(show.legend = F) +
  facet_wrap(~ MinAbundance)

