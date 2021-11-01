#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)
library(lubridate)

# load data
gis <- read_csv("intermediate-data/gis_fwc_lakewatch_fwri.csv",
                col_types = list(wkt_geom = col_character(),
                                 AOI = col_character(),
                                 Lake = col_character()))
ctrl_old <- read_csv("original-data/PrePMARS_IPMData.csv")
ctrl <- read_csv("original-data/FWC_Herbicide_Treatments_mid2010-Oct2020.csv")
herb_type <- read_csv("intermediate-data/herbicide_types.csv")


#### old control data ####

# years
ctrl_old %>% select(FY, Year) %>% unique()
# fiscal year starts July 1
# all treatments are before July (2nd year of FY range)

# rename columns
# add permanent ID based on AreaOfInterestID
# remove missing ID
# remove duplicate rows
ctrl_old2 <- ctrl_old %>%
  filter(TotalAcres > 0) %>%
  rename("Longitude_FWC" = "longitude",
         "Latitude_FWC" = "latitude",
         "SpeciesOrig" = "Species_orig",
         "TotalContFWC" = "Total_cont_fwc",
         "County_FWC" = "County") %>%
  left_join(gis %>%
              filter(CoordSource %in% c("FWC", "FWC_2")) %>%
              select(AreaOfInterestID, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique()) %>%
  mutate(County_FWC = toupper(County_FWC),
         AreaTreated_ha = TotalAcres * 0.405,
         Area_ha = ShapeArea * 100,
         AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha, # make full area if it exceeds it
                                    TRUE ~ AreaTreated_ha),
         PropTreated = AreaTreated_ha / Area_ha,
         TreatmentMethod = "unknown",
         TreatmentID = paste("old", Year, substr(Species, 1, 1), TotalAcres, sep = "_"),
         TreatmentYear = Year,
         CtrlSet = "old",
         GSYear = TreatmentYear, # assume treatments are after April
         MethodHerbicide = "unknown")

# duplicate rows?
ctrl_old2 %>%
  get_dupes()
# none

# missing ID?
ctrl_old2 %>%
  filter(is.na(PermanentID)) %>%
  select(AreaOfInterestID, AreaOfInterest) %>%
  unique() %>%
  data.frame()
# 72
# not lakes

# remove missing perm ID
# remove duplicate rows
# select relevant columns
ctrl_old3 <- ctrl_old2 %>%
  filter(!is.na(PermanentID)) %>%
  unique() %>%
  select(AreaOfInterestID, PermanentID, TreatmentYear, Species, Area_ha, AreaTreated_ha, PropTreated, TreatmentMethod, TreatmentID, MethodHerbicide, CtrlSet, GSYear)


#### control data ####

# non-herbicide methods (from herbicide_initial_visualizations)
non_herb <- c("Mechanical Harvester", 
              "Snagging (tree removal)", 
              "Aquatic Dye (for shading)", 
              "Grass Carp", "Hand Removal", 
              "Mechanical (Other)", 
              "Mechanical Shredder", 
              "Prescribed Fire")

# format date
# rename columns
# add permanent ID based on AreaOfInterestID
# remove missing ID
# remove duplicate rows
ctrl2 <- ctrl %>%
  filter(TotalAcres > 0) %>%
  mutate(BeginDate = as.Date(BeginDate, "%m/%d/%y"),
         County = toupper(County),
         MethodHerbicide = case_when(ControlMethod %in% non_herb ~ "no",
                               is.na(ControlMethod) ~ "unknown",
                               TRUE ~ "yes")) %>%
  rename("Year" = "year",
         "Month" = "month",
         "Species" = "species",
         "Herbicide" = "herbicide",
         "TotalHerbicideUsed" = "totalherbicideused",
         "Longitude_FWC" = "longitude",
         "Latitude_FWC" = "latitude",
         "County_FWC" = "County")  %>%
  left_join(gis %>%
              filter(CoordSource %in% c("FWC", "FWC_2")) %>%
              select(AreaOfInterestID, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource) %>%
              unique()) %>%
  left_join(herb_type %>%
              select(ControlMethod, ActiveIngredient)) # can add mechanism and whether it's contact - need to be formatted below


# duplicate rows?
ctrl2 %>%
  get_dupes()
# 39 rows

# missing ID?
ctrl2 %>%
  filter(is.na(PermanentID)) %>%
  select(AreaOfInterestID, AreaOfInterest) %>%
  unique() %>%
  data.frame()
# 97
# not lakes

# remove missing perm ID
# remove duplicate rows
# remove duplicate methods
ctrl3 <- ctrl2 %>%
  filter(!is.na(PermanentID)) %>%
  unique() %>%
  group_by(AreaOfInterestID, PermanentID, Species, BeginDate, TreatmentID, TotalAcres, ShapeArea) %>%  # captures area treated for an event without duplication due to multiple methods
  mutate(AreaTreated_ha = TotalAcres * 0.405,
         Area_ha = ShapeArea * 100,
         AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha, # make full area if it exceeds it
                                    TRUE ~ AreaTreated_ha),
         PropTreated = AreaTreated_ha / Area_ha,
         TreatmentMethod = paste(sort(unique(ControlMethod)), collapse = " + "), # combines control methods
         MethodHerbicide = case_when("yes" %in% unique(MethodHerbicide) ~ "herbicide",
                               "no" %in% unique(MethodHerbicide) & !("yes" %in% unique(MethodHerbicide)) ~ "not herbicide",
                               TRUE ~ "unknown"),
         ActiveIngredients = paste(sort(unique(ActiveIngredient)), collapse = ","),
         Endothall = if_else(str_detect("Endothall", ActiveIngredients) == T, 1, 0),
         TreatmentYear = year(BeginDate),
         TreatmentMonth = month(BeginDate),
         TreatmentID = as.character(TreatmentID),
         CtrlSet = "new",
         GSYear = case_when(TreatmentMonth >= 4 ~ TreatmentYear,
                            TreatmentMonth < 4 ~ TreatmentYear - 1)) %>%
  ungroup() %>%
  select(AreaOfInterestID, PermanentID, TreatmentYear, Species, Area_ha, AreaTreated_ha, PropTreated, TreatmentMethod, MethodHerbicide, ActiveIngredients, Endothall, TreatmentMonth, BeginDate, TreatmentID, CtrlSet, GSYear) %>%
  rename(TreatmentDate = BeginDate) %>%
  unique() %>%  # some duplication due to different TotalHerbicideUsed (but nothing else) for the same TreatmentID
  mutate(Endothall = replace_na(Endothall, 0)) # warnings returned for NA active ingredient


#### outputs ####

write_csv(ctrl_old3, "intermediate-data/FWC_control_old_formatted.csv")
write_csv(ctrl3, "intermediate-data/FWC_control_new_formatted.csv")

