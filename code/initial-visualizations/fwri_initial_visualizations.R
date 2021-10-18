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
fwri_gis <- read_csv("gis/intermediate-data/FWRI_plants_edited.csv")
# gis needs to be checked
# see fwri_data_formatting for details
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")


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

# hydrilla control
ctrl_hydr <- ctrl_new %>%
  filter(Species == "Hydrilla verticillata" & TotalAcres > 0) %>%
  select(PermanentID, ShapeArea, Year, TreatmentID, TotalAcres) %>%
  unique() %>% # removes duplication due to multiple herbicides with one treatment
  group_by(PermanentID, ShapeArea, Year) %>%
  summarise(Control_acres = sum(TotalAcres)) %>%
  mutate(Control = 1,
         Area_acres = ShapeArea * 247.105, # convert from square km
         Control_acres = case_when(Control_acres > Area_acres ~ Area_acres, # reduce if greater than lake size
                                   TRUE ~ Control_acres)) %>%
  select(-ShapeArea)

# add control data
fwri_hydr_ctrl <- fwri_hydr_prop %>%
  left_join(fwri_gis %>%
              select(AOI, Lake, Permanent_) %>%
              rename(PermanentID = Permanent_)) %>%
  left_join(ctrl_hydr %>%
              select(PermanentID, Area_acres) %>%
              unique()) %>%
  left_join(ctrl_hydr) %>%
  mutate(Control = replace_na(Control, 0),
         Control_acres = replace_na(Control_acres, 0),
         Occ_acres = PropCovered * Area_acres,
         Control_prop = case_when(is.na(Occ_acres) & Control == 0 ~ 0, # no shape area info (not in ctrl_hydr)
                                  Occ_acres == 0 ~ 0, # dividing by 0
                                  TRUE ~ Control_acres / Occ_acres))

ggplot(fwri_hydr_ctrl, aes(x = Year, y = PropCovered, color = AOI, fill = as.factor(Control))) +
  geom_point(show.legend = F, shape = 21) +
  geom_line(show.legend = F) +
  scale_fill_manual(values = c("white", "black")) +
  facet_wrap(~ MinAbundance)

ggplot(fwri_hydr_ctrl, aes(x = as.factor(Control), y = PropCovered)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ MinAbundance)
# treated lake-year combos have higher initial hydrilla abundance

ggplot(fwri_hydr_ctrl, aes(x = PropCovered, y = Control_prop)) +
  geom_point() +
  facet_wrap(~ MinAbundance, scales = "free")
# much larger areas treated than estimated by plant survey
# survey may have occurred after treatment within the same year
# need to format control data like I did in plant_analysis.R