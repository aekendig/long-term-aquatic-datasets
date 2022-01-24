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
inv_fwc <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")


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


#### ctrl outputs ####

write_csv(ctrl_old3, "intermediate-data/FWC_control_old_formatted.csv")
write_csv(ctrl3, "intermediate-data/FWC_control_new_formatted.csv")


#### edit old data for FWC plants ####

# species
inv_taxa <- tibble(Species = c("Hydrilla verticillata", rep("Floating Plants (Eichhornia and Pistia)", 2), # double each row that has floating plants
                               "Panicum repens", "Colocasia esculenta", "Urochloa mutica", 
                               "Alternanthera philoxeroides", "Oxycaryum cubense", "Cyperus blepharoleptos", "Salvinia minima"),
                        TaxonName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes",
                                      "Panicum repens", "Colocasia esculenta", "Urochloa mutica",
                                      "Alternanthera philoxeroides", rep("Cyperus blepharoleptos", 2), "Salvinia minima"))
# old data: "Oxycaryum cubense"
# new data: "Cyperus blepharoleptos"

# all treatments
# don't know treatment data, don't adjust GSYear
all_ctrl_old <- ctrl_old3 %>%
  group_by(PermanentID, GSYear, Area_ha) %>%
  summarize(AllAreaTreated_ha = sum(AreaTreated_ha)) %>%
  ungroup() %>%
  full_join(ctrl_old3 %>%
              select(PermanentID, Area_ha) %>%
              unique() %>%
              expand_grid(GSYear = min(ctrl_old3$GSYear):max(ctrl_old3$GSYear))) %>%
  mutate(AllAreaTreated_ha = replace_na(AllAreaTreated_ha, 0),
         AllPropTreated = AllAreaTreated_ha / Area_ha,
         CtrlSet = "old")

# focal species
foc_ctrl_old <- ctrl_old3 %>%
  filter(Species %in% inv_taxa$Species) %>%
  group_by(PermanentID, GSYear, Area_ha, Species) %>%
  summarize(AreaTreated_ha = sum(AreaTreated_ha)) %>%
  ungroup() %>%
  full_join(ctrl_old3 %>%
              select(PermanentID, Area_ha) %>%
              unique() %>%
              expand_grid(GSYear = min(ctrl_old3$GSYear):max(ctrl_old3$GSYear)) %>%
              expand_grid(Species = unique(inv_taxa$Species))) %>%
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0),
         PropTreated = AreaTreated_ha / Area_ha,
         CtrlSet = "old") %>%
  filter(Species != "Cyperus blepharoleptos") # new name of Cuban bulrush added with expand_grid


#### edit new data for FWC plants ####

# all treatments
all_ctrl <- ctrl3 %>%
  left_join(inv_fwc %>%
              select(PermanentID, SurveyDate, GSYear) %>% 
              unique()) %>% # add survey dates for each lake and year
  mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
         GSYear = case_when(is.na(SurveyTreatDays) | SurveyTreatDays >= 14 ~ GSYear, # no survey/survey after treatment -> keep year
                            SurveyTreatDays < 14 ~ GSYear + 1)) %>% # survey before or very soon after treatment -> move treatment to next year 
  full_join(ctrl3 %>%
              select(PermanentID, Area_ha) %>%
              unique() %>%
              expand_grid(GSYear = min(ctrl3$GSYear):max(ctrl3$GSYear))) %>% # could add +1 to max, but don't want to include last year anyway (likely incomplete)
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0)) %>% # add zeros for missing years
  group_by(PermanentID, GSYear, Area_ha) %>%
  summarize(AllAreaTreated_ha = sum(AreaTreated_ha)) %>%
  ungroup() %>%
  mutate(AllPropTreated = AllAreaTreated_ha / Area_ha,
         CtrlSet = "new")

# focal species
foc_ctrl <- ctrl3 %>%
  filter(Species %in% inv_taxa$Species) %>%
  left_join(inv_fwc %>%
              select(PermanentID, SurveyDate, GSYear) %>% 
              unique()) %>% # add survey dates for each lake and year
  mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
         GSYear = case_when(is.na(SurveyTreatDays) | SurveyTreatDays >= 14 ~ GSYear, # no survey/survey after treatment -> keep year
                            SurveyTreatDays < 14 ~ GSYear + 1)) %>% # survey before or very soon after treatment -> move treatment to next year 
  full_join(ctrl3 %>%
              select(PermanentID, Area_ha) %>%
              unique() %>%
              expand_grid(GSYear = min(ctrl3$GSYear):max(ctrl3$GSYear)) %>% # could add +1 to max, but don't want to include last year anyway (likely incomplete)
              expand_grid(Species = unique(inv_taxa$Species))) %>% 
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0)) %>% # add zeros for missing years
  group_by(PermanentID, GSYear, Area_ha, Species) %>%
  summarize(AreaTreated_ha = sum(AreaTreated_ha),
            TreatmentDays = case_when(AreaTreated_ha > 0 ~ as.numeric(n_distinct(TreatmentDate)),
                                      TRUE ~ 0),
            TreatmentDate = if_else(AreaTreated_ha > 0, max(TreatmentDate), NA_real_)) %>%
  ungroup() %>%
  mutate(PropTreated = AreaTreated_ha / Area_ha,
         CtrlSet = "new") %>%
  filter(Species != "Oxycaryum cubense") # old name of Cuban bulrush added with expand_grid


#### combine data ####

# combine 2010 data (assume records between old and new are distinct)
all_ctrl_2010 <- all_ctrl_old %>%
  filter(GSYear == 2010) %>%
  full_join(all_ctrl %>%
              filter(GSYear == 2010)) %>%
  group_by(PermanentID, GSYear, Area_ha) %>%
  summarize(AllAreaTreated_ha = sum(AllAreaTreated_ha)) %>%
  ungroup() %>%
  mutate(AllPropTreated = AllAreaTreated_ha / Area_ha,
         CtrlSet = "old/new")

foc_ctrl_2010 <- foc_ctrl_old %>%
  filter(GSYear == 2010) %>%
  mutate(Species = if_else(Species == "Oxycaryum cubense", # update old name to avoid duplicate records
                           "Cyperus blepharoleptos", Species)) %>%
  full_join(foc_ctrl %>%
              filter(GSYear == 2010)) %>%
  group_by(PermanentID, GSYear, Area_ha, Species) %>%
  summarize(AreaTreated_ha = sum(AreaTreated_ha)) %>%
  ungroup() %>%
  mutate(PropTreated = AreaTreated_ha / Area_ha,
         CtrlSet = "old/new") %>%
  left_join(foc_ctrl %>%
              filter(GSYear == 2010) %>%
              select(PermanentID, Species, TreatmentDays, TreatmentDate)) # add treatment info

ctrl_2010 <- all_ctrl_2010 %>%
  full_join(foc_ctrl_2010)

# combine
inv_ctrl <- all_ctrl_old %>%
  full_join(foc_ctrl_old) %>%
  full_join(all_ctrl %>%
              full_join(foc_ctrl)) %>%
  filter(GSYear != 2010) %>%
  full_join(ctrl_2010) %>%
  left_join(inv_taxa) # duplicate each floating plant row for the two species, add taxon names

# check for redundancy
inv_ctrl %>%
  select(PermanentID, GSYear, Species) %>%
  get_dupes() %>%
  select(Species, dupe_count) %>%
  unique()

inv_ctrl %>%
  select(PermanentID, GSYear, TaxonName) %>%
  get_dupes()


#### lag intervals ####

# function for cumulative treatment
ctrl_lag_fun <- function(GSYear, Lag){
  
  # change name
  year1 <- GSYear
  
  # filter dataset
  subdat <- inv_ctrl %>%
      filter(GSYear <= year1 & GSYear >= (year1 - Lag))
  
  # summarize
  outdat <- subdat %>%
    mutate(AllTreated = if_else(AllAreaTreated_ha > 0, 1, 0),
           Treated = ifelse(AreaTreated_ha > 0, 1, 0)) %>%
    group_by(PermanentID, TaxonName, Species) %>% # summarize over lag years
    summarise(AllPropTreated = mean(AllPropTreated), # average proportion treated per year 
              PropTreated = mean(PropTreated),
              AllTreated = mean(AllTreated), # proportion of years with treatment
              Treated = mean(Treated),
              NYears = n_distinct(GSYear)) %>%
    ungroup() %>%
    mutate(GSYear = year1,
           Lag = Lag)
  
  # return
  return(outdat)
}

# treatments by year and lag
inv_ctrl2 <- inv_ctrl %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 0:5)) %>% # remove repeat row for each species
  pmap(ctrl_lag_fun) %>% # summarizes for each GS, Lag, PermID, and Sp
  bind_rows() %>%
  filter(NYears == Lag + 1) %>% # all years for a lag must be available
  select(-NYears) %>%
  pivot_wider(names_from = Lag,
              values_from = c(AllPropTreated, PropTreated, AllTreated, Treated),
              names_glue = "Lag{Lag}{.value}") %>% # make treatments wide by lag
  filter(!is.na(Lag0AllPropTreated)) %>% # must have measurement in final year 
  full_join(inv_ctrl)


#### inv ctrl output ####
write_csv(inv_ctrl2, "intermediate-data/FWC_invasive_control_formatted.csv")
