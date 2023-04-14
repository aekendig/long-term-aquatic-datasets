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
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")


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

# consistent waterbody name?
ctrl_old2 %>%
  distinct(AreaOfInterestID, AreaOfInterest) %>%
  get_dupes(AreaOfInterestID)

# remove missing perm ID
# remove duplicate rows
# select relevant columns
ctrl_old3 <- ctrl_old2 %>%
  filter(!is.na(PermanentID)) %>%
  unique() %>%
  select(AreaOfInterestID, AreaOfInterest, PermanentID, TreatmentYear, 
         Species, Area_ha, AreaTreated_ha, PropTreated, TreatmentMethod, 
         TreatmentID, MethodHerbicide, CtrlSet, GSYear)

# one permID per AOI?
ctrl_old3 %>%
  group_by(AreaOfInterestID) %>%
  summarize(nPerm = n_distinct(PermanentID)) %>%
  ungroup() %>%
  filter(nPerm > 1)
# yes


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

# consistent waterbody name?
ctrl2 %>%
  distinct(AreaOfInterestID, AreaOfInterest) %>%
  get_dupes(AreaOfInterestID)

# remove missing perm ID
# remove duplicate rows
# remove duplicate methods
ctrl3 <- ctrl2 %>%
  filter(!is.na(PermanentID)) %>%
  unique() %>%
  group_by(AreaOfInterestID, AreaOfInterest, PermanentID, Species, BeginDate, TreatmentID, TotalAcres, ShapeArea) %>%  # captures area treated for an event without duplication due to multiple methods
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
         TreatmentYear = year(BeginDate),
         TreatmentMonth = month(BeginDate),
         TreatmentID = as.character(TreatmentID),
         CtrlSet = "new",
         GSYear = case_when(TreatmentMonth >= 4 ~ TreatmentYear,
                            TreatmentMonth < 4 ~ TreatmentYear - 1)) %>%
  ungroup() %>%
  select(AreaOfInterestID, AreaOfInterest, PermanentID, TreatmentYear, Species, Area_ha, AreaTreated_ha, PropTreated, TreatmentMethod, MethodHerbicide, ActiveIngredients, TreatmentMonth, BeginDate, TreatmentID, CtrlSet, GSYear) %>%
  rename(TreatmentDate = BeginDate) %>%
  unique()  # some duplication due to different TotalHerbicideUsed (but nothing else) for the same TreatmentID

# one permID per AOI?
ctrl3 %>%
  group_by(AreaOfInterestID) %>%
  summarize(nPerm = n_distinct(PermanentID)) %>%
  ungroup() %>%
  filter(nPerm > 1)
# yes


#### ctrl outputs ####

write_csv(ctrl_old3, "intermediate-data/FWC_control_old_formatted.csv")
write_csv(ctrl3, "intermediate-data/FWC_control_new_formatted.csv")


#### combine datasets ####

# duplicates within old dataset
ctrl_old3 %>%
  get_dupes(AreaOfInterestID, TreatmentYear, Species, AreaTreated_ha)
# one, but it's not an aquatic species

# duplicates within new dataset
ctrl3 %>%
  get_dupes(AreaOfInterestID, TreatmentYear, Species, 
            AreaTreated_ha, TreatmentMethod, TreatmentID)

# same waterbody names?
ctrl_old3 %>%
  distinct(AreaOfInterestID, AreaOfInterest) %>%
  full_join(ctrl3 %>%
              distinct(AreaOfInterestID, AreaOfInterest)) %>%
  get_dupes(AreaOfInterestID)

# combine datasets
ctrl4 <- ctrl3 %>%
  full_join(ctrl_old3 %>%
              mutate(TreatmentDate = as.Date(paste0(GSYear, "-01-01")))) %>%
  mutate(AreaOfInterest = if_else(AreaOfInterest == "Watermelon Pond" & AreaOfInterestID == 465,
                                  "Watermellon Pond", # correct dual-spelled name
                                  AreaOfInterest))

# waterbodies
n_distinct(ctrl4$AreaOfInterestID)

# check for duplication around 2010
ctrl4_dups <- ctrl4 %>%
  filter(TreatmentYear == 2010) %>%
  get_dupes(AreaOfInterestID, Species, AreaTreated_ha) %>%
  group_by(AreaOfInterestID, Species, AreaTreated_ha, TreatmentYear) %>%
  summarize(Sets = n_distinct(CtrlSet)) %>%
  ungroup() %>%
  filter(Sets > 1)
# only 6 examples, 3 will be included in analyses

# examine duplicates
inner_join(ctrl4, ctrl4_dups) %>%
  arrange(Species, AreaOfInterestID)
# one rep is in new dataset, one is in old, use new

# are AOIs for each Perm ID consistent?
ctrl4 %>%
  group_by(PermanentID) %>%
  summarize(nAOI = n_distinct(AreaOfInterest),
            nAOI_ID = n_distinct(AreaOfInterestID)) %>%
  ungroup() %>%
  filter(nAOI_ID > 1) %>% # select lakes with multiple AOIs
  inner_join(ctrl4 %>%
               group_by(PermanentID, GSYear) %>%
               summarize(AOI = paste(sort(unique(AreaOfInterest)), collapse = ", ")) %>%
               ungroup() %>%
               group_by(PermanentID, AOI) %>%
               summarize(Years = paste(sort(unique(GSYear)), collapse = ", ")) %>%
               ungroup()) %>% # AOIs surveyed each year
  get_dupes(PermanentID) # select lakes with different AOIs over time
# consistent with plants (fwc_plant_data_processing), but
# leave Little Santa Fe because it was surveyed for plants

# first plant survey date of growing season
plant_surveys <- plant_fwc %>%
  group_by(AreaOfInterestID, PermanentID, GSYear) %>%
  summarize(SurveyDate = min(SurveyDate)) %>%
  ungroup()

# check that AOIs are consistent
ctrl4 %>%
  filter(PermanentID %in% plant_surveys$PermanentID & 
           !(AreaOfInterestID %in% plant_surveys$AreaOfInterestID)) %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID) %>%
  unique()
# DeLeon Springs St Park Wetlands

plant_fwc %>%
  filter(PermanentID == "107776163") %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID) %>%
  unique()
# okay to remove the wetlands from the ctrl dataset
# different area than lake in both datasets

# remove additional AOIs (temporal consistency)
# add plant survey date
# adjust GSYear based on plant survey
ctrl5 <- ctrl4 %>%
  anti_join(ctrl4_dups %>%
              mutate(CtrlSet = "old")) %>% # should remove 6 rows
  filter( #!(AreaOfInterest == "Wauseon Bay" & PermanentID == "112029141") & # remove if adjusting for temporal consistency
           # !(AreaOfInterest == "Red Water, Lake" & PermanentID == "112047993") &
           # !(AreaOfInterest == "Josephine Creek" & PermanentID == "112049879") &
           # !(AreaOfInterest == "Hunt, Lake" & PermanentID == "167180956") &
           !(AreaOfInterest == "DeLeon Springs St Park Wetlands" & PermanentID == "107776163")) %>%
  left_join(plant_surveys) %>%
  mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
         GSYear = case_when(is.na(SurveyTreatDays) | SurveyTreatDays >= 14 ~ GSYear, # no survey/survey after treatment -> keep year
                            SurveyTreatDays < 14 ~ GSYear + 1), # survey before or very soon after treatment -> move treatment to next year
         Species = if_else(Species == "Oxycaryum cubense", "Cyperus blepharoleptos", Species)) # use more recent Cuban bulrush name


#### summarize data ####

# species
inv_taxa <- tibble(Species = c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)", # double each row that has floating plants
                               "Panicum repens", "Urochloa mutica", "Cyperus blepharoleptos"))

# all treatments (for plant analyses)
# fill in missing years
# summarize by permanent ID and GS year
all_ctrl <- ctrl5 %>%
  full_join(ctrl5 %>%
              select(AreaOfInterestID, AreaOfInterest, PermanentID, Area_ha) %>%
              unique() %>%
              expand_grid(GSYear = min(ctrl5$GSYear):max(ctrl5$GSYear))) %>% # could add +1 to max, but don't want to include last year anyway (likely incomplete)
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0)) %>% # add zeros for missing years (did not treat)
  group_by(AreaOfInterestID, AreaOfInterest, PermanentID, GSYear, Area_ha) %>%
  summarize(AllAreaTreated_ha = sum(AreaTreated_ha)) %>% # add all treatments within a year
  ungroup() %>%
  mutate(AllTreated = if_else(AllAreaTreated_ha > 0, 1, 0))

# focal species (for plant analyses)
# fill in missing years
# summarize by permanent ID and GSYear
foc_ctrl <- ctrl5 %>%
  filter(Species %in% inv_taxa$Species) %>%
  full_join(ctrl5 %>%
              select(AreaOfInterestID, AreaOfInterest, PermanentID, Area_ha) %>%
              unique() %>%
              expand_grid(GSYear = min(ctrl5$GSYear):max(ctrl5$GSYear)) %>% # could add +1 to max, but don't want to include last year anyway (likely incomplete)
              expand_grid(Species = unique(inv_taxa$Species))) %>% 
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0)) %>% # add zeros for missing years
  group_by(AreaOfInterestID, AreaOfInterest, PermanentID, GSYear, Area_ha, Species) %>%
  summarize(AreaTreated_ha = sum(AreaTreated_ha),  # add all treatments within a year
            TreatmentDays = case_when(AreaTreated_ha > 0 ~ as.numeric(n_distinct(TreatmentDate)),
                                      TRUE ~ 0),
            TreatmentDate = max(TreatmentDate)) %>%
  ungroup() %>%
  mutate(Treated = ifelse(AreaTreated_ha > 0, 1, 0))

# combine
inv_ctrl <- all_ctrl %>%
  full_join(foc_ctrl)

# check for redundancy
inv_ctrl %>%
  select(AreaOfInterestID, GSYear, Species) %>%
  get_dupes()

# all treatments (for quality analyses)
# fill in missing years
# summarize by permanent ID and actual year
qual_all_ctrl <- ctrl5 %>%
  full_join(ctrl5 %>%
              select(PermanentID, AreaOfInterest, Area_ha) %>%
              unique() %>%
              expand_grid(GSYear = min(ctrl5$GSYear):max(ctrl5$GSYear))) %>% # could add +1 to max, but don't want to include last year anyway (likely incomplete)
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0)) %>% # add zeros for missing years (did not treat)
  group_by(PermanentID, GSYear, Area_ha) %>%
  summarize(AreaName = paste(sort(unique(AreaOfInterest)), collapse = "/"),
            AllAreaTreated_ha = sum(AreaTreated_ha)) %>% # add all treatments within a year
  ungroup() %>%
  mutate(AllTreated = if_else(AllAreaTreated_ha > 0, 1, 0))

# focal treatments (for quality analyses)
# fill in missing years
# summarize by permanent ID and actual year
qual_foc_ctrl <- ctrl5 %>%
  filter(Species %in% inv_taxa$Species) %>%
  full_join(ctrl5 %>%
              select(PermanentID, AreaOfInterest, Area_ha) %>%
              unique() %>%
              expand_grid(GSYear = min(ctrl5$GSYear):max(ctrl5$GSYear)) %>% # could add +1 to max, but don't want to include last year anyway (likely incomplete)
              expand_grid(Species = unique(inv_taxa$Species))) %>% 
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0)) %>% # add zeros for missing years
  group_by(PermanentID, GSYear, Area_ha, Species) %>%
  summarize(AreaName = paste(sort(unique(AreaOfInterest)), collapse = "/"),
            AreaTreated_ha = sum(AreaTreated_ha),  # add all treatments within a year
            TreatmentDays = case_when(AreaTreated_ha > 0 ~ as.numeric(n_distinct(TreatmentDate)),
                                      TRUE ~ 0),
            TreatmentDate = max(TreatmentDate)) %>%
  ungroup() %>%
  mutate(Treated = ifelse(AreaTreated_ha > 0, 1, 0)) 

# combine
qual_ctrl <- qual_all_ctrl %>%
  full_join(qual_foc_ctrl)

# check for redundancy
qual_ctrl %>%
  select(PermanentID, GSYear, Species) %>%
  get_dupes() %>%
  select(Species, dupe_count) %>%
  unique()

#### lag intervals ####

# function for cumulative treatment
ctrl_lag_fun <- function(Year, Lag, dat_in){
  
  # change name
  year1 <- Year
  
  # filter dataset
  subdat <- dat_in %>%
    filter(GSYear <= year1 & GSYear > (year1 - Lag)) # average frequency over preceeding years up to lag years
  
  # summarize
  outdat <- subdat %>%
    group_by(Location, Species) %>% # summarize over lag years
    summarise(AllTreated = mean(AllTreated), # proportion of years with treatment
              Treated = mean(Treated),
              NYears = n_distinct(GSYear)) %>%
    ungroup()
  
  # return
  return(outdat)
}

# treatments by year and lag
inv_ctrl2 <- inv_ctrl %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 2:6)) %>% # remove repeat row for each species
  mutate(out = pmap(., function(GSYear, Lag)  # summarizes for each GS, Lag, PermID, and Sp
    ctrl_lag_fun(Year = GSYear, Lag = Lag, 
                 dat_in = inv_ctrl %>%
                   mutate(Location = AreaOfInterestID)))) %>%
  unnest(cols = out) %>%
  filter(NYears == Lag) %>% # all years for a lag must be available
  select(-NYears) %>%
  pivot_wider(names_from = Lag,
              values_from = c(AllTreated, Treated),
              names_glue = "Lag{Lag}{.value}") %>% # make treatments wide by lag
  rename(AreaOfInterestID = Location) %>%
  full_join(inv_ctrl)

# treatments by year and lag
qual_ctrl2 <- qual_ctrl %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 2:6)) %>% # remove repeat row for each species
  mutate(out = pmap(., function(GSYear, Lag)  # summarizes for each GS, Lag, PermID, and Sp
    ctrl_lag_fun(Year = GSYear, Lag = Lag, 
                 dat_in = qual_ctrl %>%
                   mutate(Location = PermanentID)))) %>%
  unnest(cols = out) %>%
  filter(NYears == Lag) %>% # all years for a lag must be available
  select(-NYears) %>%
  pivot_wider(names_from = Lag,
              values_from = c(AllTreated, Treated),
              names_glue = "Lag{Lag}{.value}") %>% # make treatments wide by lag
  rename(PermanentID = Location) %>%
  full_join(qual_ctrl)


#### years since treatment ####

# max lag
inv_max_lag <- max(inv_ctrl2$GSYear) - min(inv_ctrl2$GSYear)
qual_max_lag <- max(qual_ctrl2$GSYear) - min(qual_ctrl2$GSYear)

# add column for last treatment
inv_ctrl3 <- inv_ctrl2 %>%
  arrange(AreaOfInterestID, Species, GSYear) %>%
  group_by(AreaOfInterestID, Species) %>%
  mutate(LastTreatment = if_else(lag(Treated, n = inv_max_lag) > 0,
                                 inv_max_lag,
                                 NA_real_))

qual_ctrl3 <- qual_ctrl2 %>%
  arrange(PermanentID, Species, GSYear) %>%
  group_by(PermanentID, Species) %>%
  mutate(LastTreatment = if_else(lag(Treated, n = qual_max_lag) > 0,
                                 qual_max_lag,
                                 NA_real_))

# cycle through lags
for(i in (inv_max_lag - 1):0) {
  
  inv_ctrl3 <- inv_ctrl3 %>%
    mutate(LastTreatment = if_else(lag(Treated, n = i) > 0,
                                   as.numeric(i),
                                   LastTreatment))
    
}

for(i in (qual_max_lag - 1):0) {
  
  qual_ctrl3 <- qual_ctrl3 %>%
    mutate(LastTreatment = if_else(lag(Treated, n = i) > 0,
                                   as.numeric(i),
                                   LastTreatment))
  
}

# ungroup
inv_ctrl4 <- inv_ctrl3 %>%
  ungroup() %>%
  mutate(RecentTreatment = 1 / (LastTreatment + 1),
         RecentTreatment = replace_na(RecentTreatment, 0)) # no record of treatment

qual_ctrl4 <- qual_ctrl3 %>%
  ungroup() %>%
  mutate(RecentTreatment = 1 / (LastTreatment + 1),
         RecentTreatment = replace_na(RecentTreatment, 0))

# check that it worked (use hydrilla because it's managed the most)
aois <- unique(inv_ctrl4$AreaOfInterestID)

for(i in aois) {
  
  print(inv_ctrl4 %>%
          filter(Species == "Hydrilla verticillata" & AreaOfInterestID == i) %>%
          ggplot(aes(x = GSYear, color = AreaOfInterestID)) +
          geom_line(aes(y = Treated)) +
          geom_point(aes(y = LastTreatment)) +
          theme_bw() +
          theme(legend.position = "none"))
  
}

ggplot(inv_ctrl4, aes(x = LastTreatment)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ Species)


#### inv ctrl output ####
write_csv(inv_ctrl4, "intermediate-data/FWC_invasive_control_formatted.csv")
write_csv(qual_ctrl4, "intermediate-data/FWC_quality_control_formatted.csv")
