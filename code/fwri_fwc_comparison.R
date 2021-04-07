#### info ####

# goal: see how FWRI compares to FWC


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)

# import data
fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv",
                 col_types = list(PermanentID = col_character()))
fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")


#### subset data functions ####

fwri_fun <- function(code){
  
  dat_out <- fwri %>%
    select(AOI, Lake, PermanentID, Year, Date, Site) %>% # start with all surveys (no species means abundance = 0)
    unique() %>%  # don't need a row for each species
    left_join(fwri %>%  # add species information
                filter(Code == code) %>%
                select(AOI, Lake, PermanentID, Year, Date, Site, Abundance)) %>%
    mutate(Abundance = replace_na(Abundance, 0), # species cover 0 when it wasn't in a survey
           FWRI_Date = as.Date(Date),
           FWRI_Month = month(Date)) %>%
    group_by(AOI, Lake, PermanentID, Year, FWRI_Month, FWRI_Date) %>% # checked above for duplicates within sites - none
    summarise(Area_sites = n(),
              AreaCovered_sites = sum(Abundance > 0)) %>%
    ungroup() %>%
    mutate(AOI = case_when(AOI %in% c("NorthConway", "SouthConway") ~ "Conway", # combine these lakes (one lake in FWC data)
                           TRUE ~ AOI),
           Lake = case_when(Lake %in% c("North Conway", "South Conway") ~ "Conway",
                            TRUE ~ Lake)) %>%
    group_by(AOI, Lake, PermanentID, Year, FWRI_Month) %>%
    summarise(Area_sites = sum(Area_sites), # combine north and south conway data
              AreaCovered_sites = sum(AreaCovered_sites),
              FWRI_Date = min(FWRI_Date)) %>%
    ungroup() %>%
    mutate(FWRI_PropCovered = AreaCovered_sites / Area_sites) %>%
    select(-c(Area_sites, AreaCovered_sites)) %>%
    rename(FWRI_Year = Year) %>%
    filter(!(AOI %in% c("OrangeOpenWater", "Eustis2")))
  
  return(dat_out)
}

fwc_fun <- function(species){
  
  dat_out <- fwc %>% # start with all surveys (no species means abundance = 0)
    filter(SpeciesName != species) %>%
    select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate) %>%
    unique() %>% # don't need a row for each species
    left_join(fwc %>% # add species information
                filter(SpeciesName == species) %>%
                select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, SpeciesAcres)) %>%
    mutate(SpeciesAcres = replace_na(SpeciesAcres, 0), # species cover 0 when it wasn't in a survey
           Area_ha = ShapeArea * 100, # convert lake area from km-squared to hectares
           AreaCovered_ha = SpeciesAcres * 0.405, # convert plant cover from acres to hectares
           AreaCovered_ha = case_when(AreaCovered_ha > Area_ha ~ Area_ha,
                                      TRUE ~ AreaCovered_ha), # make plant cover the size of the lake area if it exceeds it
           FWC_PropCovered = AreaCovered_ha / Area_ha,
           FWC_Year = year(SurveyDate),
           FWC_Month = month(SurveyDate)) %>%
    group_by(AreaOfInterestID, FWC_Year) %>% # remove reports of zero when there is another report that year
    mutate(AreaCoveredAnnAvg_ha = mean(AreaCovered_ha),
           FirstSurveyPerYear = min(SurveyDate)) %>%
    ungroup() %>%
    filter(!(AreaCoveredAnnAvg_ha > 0 & AreaCovered_ha == 0) & # kept report is non-zero
             !(AreaCovered_ha == 0 & SurveyDate != FirstSurveyPerYear)) %>% # kept report is zero
    select(AreaOfInterest, AreaOfInterestID, PermanentID, FWC_Year, FWC_Month, SurveyDate, FWC_PropCovered) %>%
    rename(FWC_Date = SurveyDate) %>%
    filter(!(AreaOfInterestID %in% c(402, 469, 42, 220, 436)))
  # Harris (177) and Little Harris (244) are the same lake
  # Red Water (365) and Little Red Water (250) are the same lake
  
  return(dat_out)
}


#### compare hydrilla cover ####

# subset data
fwri_hydr <- fwri_fun("HYDR")
fwc_hydr <- fwc_fun("Hydrilla verticillata")

# check lakes that will be excluded
anti_join(fwri_hydr, fwc_hydr) %>%
  select(AOI, PermanentID) %>%
  unique()

# check date formatting
tz(fwri_hydr$FWRI_Date)
tz(fwc_hydr$FWC_Date)

# combine data
hydr <- fwri_hydr %>%
  inner_join(fwc_hydr) %>%
  mutate(TimeDiff = abs(FWC_Date - FWRI_Date)) %>%
  filter(TimeDiff <= 150)

# duplicate AOI's
hydr %>%
  group_by(PermanentID) %>%
  summarise(FWRI_AOI = length(unique(AOI)),
            FWC_AOI = length(unique(AreaOfInterestID))) %>%
  ungroup() %>%
  filter(FWRI_AOI > 1 | FWC_AOI > 1) %>%
  left_join(hydr %>%
              select(PermanentID, AOI, AreaOfInterest, AreaOfInterestID) %>%
              unique())
# two lakes that have same coordinates in FWC dataset

# duplicate FWC dates for each FWRI
hydr %>%
  group_by(PermanentID, AreaOfInterestID, FWRI_Date) %>%
  summarise(FWC_surveys = length(unique(FWC_Date))) %>%
  ungroup() %>%
  filter(FWC_surveys > 1)

# figure
ggplot(hydr, aes(FWC_PropCovered, FWRI_PropCovered, 
                 group = as.factor(AreaOfInterestID), 
                 color = as.numeric(TimeDiff))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point()


#### compare water hyacinth cover ####

# subset data
fwri_wahy <- fwri_fun("WAHY")
fwc_wahy <- fwc_fun("Eichhornia crassipes")

# combine data
wahy <- fwri_wahy %>%
  inner_join(fwc_wahy) %>%
  mutate(TimeDiff = abs(FWC_Date - FWRI_Date)) %>%
  filter(TimeDiff <= 150)

# duplicate AOI's
wahy %>%
  group_by(PermanentID) %>%
  summarise(FWRI_AOI = length(unique(AOI)),
            FWC_AOI = length(unique(AreaOfInterestID))) %>%
  ungroup() %>%
  filter(FWRI_AOI > 1 | FWC_AOI > 1) %>%
  left_join(wahy %>%
              select(PermanentID, AOI, AreaOfInterest, AreaOfInterestID) %>%
              unique())
# two lakes that have same coordinates in FWC dataset

# duplicate FWC dates for each FWRI
wahy %>%
  group_by(PermanentID, AreaOfInterestID, FWRI_Date) %>%
  summarise(FWC_surveys = length(unique(FWC_Date))) %>%
  ungroup() %>%
  filter(FWC_surveys > 1)

# figure
ggplot(wahy, aes(FWC_PropCovered, FWRI_PropCovered, 
                 group = as.factor(AreaOfInterestID), 
                 color = as.numeric(TimeDiff))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point()


#### compare water lettuce cover ####

# subset data
fwri_wale <- fwri_fun("WALE")
fwc_wale <- fwc_fun("Pistia stratiotes")

# combine data
wale <- fwri_wale %>%
  inner_join(fwc_wale) %>%
  mutate(TimeDiff = abs(FWC_Date - FWRI_Date)) %>%
  filter(TimeDiff <= 150)

# duplicate AOI's
wale %>%
  group_by(PermanentID) %>%
  summarise(FWRI_AOI = length(unique(AOI)),
            FWC_AOI = length(unique(AreaOfInterestID))) %>%
  ungroup() %>%
  filter(FWRI_AOI > 1 | FWC_AOI > 1) %>%
  left_join(wale %>%
              select(PermanentID, AOI, AreaOfInterest, AreaOfInterestID) %>%
              unique())
# two lakes that have same coordinates in FWC dataset

# duplicate FWC dates for each FWRI
wale %>%
  group_by(PermanentID, AreaOfInterestID, FWRI_Date) %>%
  summarise(FWC_surveys = length(unique(FWC_Date))) %>%
  ungroup() %>%
  filter(FWC_surveys > 1)

# figure
ggplot(wale, aes(FWC_PropCovered, FWRI_PropCovered, 
                 group = as.factor(AreaOfInterestID), 
                 color = as.numeric(TimeDiff))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point()


#### combined figure ####

# combine
compDat <- hydr %>%
  mutate(Species = "Hydrilla verticillata") %>%
  full_join(wahy %>%
              mutate(Species = "Eichhornia crassipes")) %>%
  full_join(wale %>%
              mutate(Species = "Pistia stratiotes")) %>%
  mutate(DaysDiff = as.numeric(TimeDiff))

# check duplicate lakes
filter(compDat, PermanentID == "112047993") # all 0's
filter(compDat, PermanentID == "120024301") %>% data.frame() # different values

# remove duplicates
compDat2 <- compDat %>%
  filter(AreaOfInterestID != 365)

# figure
pdf("output/fwri_fwc_comparison_invasive_cover.pdf", width = 7.5, height = 4)
ggplot(compDat, aes(FWC_PropCovered, FWRI_PropCovered, 
                 group = as.factor(AreaOfInterestID), 
                 color = DaysDiff)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  facet_wrap(~Species, scales = "free") +
  xlab("Proportion of lake (FWC)") +
  ylab("Proportion of lake (FWRI)")
dev.off()
