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
inv_fwc <- read_csv("intermediate-data/FWC_only_invasive_plant_formatted.csv")


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

# species
inv_taxa <- tibble(Species = c("Hydrilla verticillata", "Eichhornia crassipes", "Pistia stratiotes",
                               "Floating Plants (Eichhornia and Pistia)", # combine with above
                               "Panicum repens", "Urochloa mutica", 
                               "Cyperus blepharoleptos", "Oxycaryum cubense")) # combine
# combine datasets
ctrl4 <- ctrl3 %>%
  full_join(ctrl_old3 %>%
              mutate(TreatmentDate = as.Date(paste0(GSYear, "-01-01")))) %>%
  inner_join(inv_taxa) %>%
  mutate(AreaOfInterest = if_else(AreaOfInterest == "Watermelon Pond" & AreaOfInterestID == 465,
                                  "Watermellon Pond", # correct dual-spelled name
                                  AreaOfInterest),
         Species = if_else(Species == "Oxycaryum cubense", "Cyperus blepharoleptos", Species))

# assign floating plants to both species
ctrl5 <- ctrl4 %>%
  filter(Species == "Floating Plants (Eichhornia and Pistia)") %>%
  mutate(Species = "Eichhornia crassipes",
         Floating = 1) %>%
  full_join(ctrl4 %>%
              filter(Species == "Floating Plants (Eichhornia and Pistia)") %>%
              mutate(Species = "Pistia stratiotes",
                     Floating = 1)) %>%
  full_join(ctrl4 %>%
              filter(Species != "Floating Plants (Eichhornia and Pistia)") %>%
              mutate(Floating = 0))

# duplicates due to name differences
ctrl5 %>%
  get_dupes(AreaOfInterestID, Species, TreatmentDate, AreaTreated_ha, TreatmentID) %>%
  select(AreaOfInterestID, Species, TreatmentDate, AreaTreated_ha, CtrlSet, Floating)

# note there are multiple treatments on the same date with different TreatmentID
# decide whether to sum area or use single area if area treated is used
  
# check for duplication around 2010
# can use date from new dataset?
(ctrl5_2010_dups <- ctrl5 %>%
  filter(TreatmentYear == 2010) %>%
    group_by(AreaOfInterestID, Species, AreaTreated_ha) %>%
    summarize(Sets = n_distinct(CtrlSet)) %>%
    ungroup() %>%
    filter(Sets > 1))
# use dates from new dataset for these

# are AOIs for each Perm ID consistent?
ctrl5 %>%
  group_by(PermanentID) %>%
  summarize(nAOI = n_distinct(AreaOfInterest),
            nAOI_ID = n_distinct(AreaOfInterestID)) %>%
  ungroup() %>%
  filter(nAOI_ID > 1) %>% # select lakes with multiple AOIs
  inner_join(ctrl5 %>%
               group_by(PermanentID, GSYear) %>%
               summarize(AOI = paste(sort(unique(AreaOfInterest)), collapse = ", ")) %>%
               ungroup() %>%
               group_by(PermanentID, AOI) %>%
               summarize(Years = paste(sort(unique(GSYear)), collapse = ", ")) %>%
               ungroup()) %>% # AOIs surveyed each year
  get_dupes(PermanentID) # select lakes with different AOIs over time
# consistent with plants (fwc_plant_data_processing), but
# leave Little Santa Fe because it was surveyed for plants

# plant survey times
inv_surveys <- inv_fwc %>%
  select(AreaOfInterestID, AreaOfInterest, PermanentID, TaxonName, GSYear, SurveyDate) %>%
  rename(Species = TaxonName)

# check that AOIs are consistent
ctrl5 %>%
  filter(PermanentID %in% inv_surveys$PermanentID & 
           !(AreaOfInterestID %in% inv_surveys$AreaOfInterestID)) %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID) %>%
  unique()
# DeLeon Springs St Park Wetlands

inv_fwc %>%
  filter(PermanentID == "107776163") %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID) %>%
  unique()
# okay to remove the wetlands from the ctrl dataset
# different area than lake in both datasets

# check that waterbody names are consistent
ctrl5 %>%
  distinct(AreaOfInterestID, AreaOfInterest) %>%
  full_join(inv_surveys %>%
              distinct(AreaOfInterestID, AreaOfInterest)) %>%
  get_dupes(AreaOfInterestID)
# yes

# are all treatment areas > 0?
filter(ctrl5, AreaTreated_ha == 0 | is.na(AreaTreated_ha))
# yes

# remove additional AOIs (temporal consistency)
# select treatment occurrences (ignoring area, herbicide, etc.)
# add plant survey date
# adjust GSYear based on survey date
# summarize by permanent ID and GSYear
ctrl6 <- ctrl5 %>%
  anti_join(ctrl5_2010_dups %>%
              mutate(CtrlSet = "old")) %>% # should remove 10 rows
  filter( #!(AreaOfInterest == "Wauseon Bay" & PermanentID == "112029141") & # remove if adjusting for temporal consistency
           # !(AreaOfInterest == "Red Water, Lake" & PermanentID == "112047993") &
           # !(AreaOfInterest == "Josephine Creek" & PermanentID == "112049879") &
           # !(AreaOfInterest == "Hunt, Lake" & PermanentID == "167180956") &
           !(AreaOfInterest == "DeLeon Springs St Park Wetlands" & PermanentID == "107776163")) %>%
  distinct(AreaOfInterestID, AreaOfInterest, PermanentID, Species, TreatmentYear, TreatmentDate, GSYear) %>%
  inner_join(inv_taxa) %>% # select focal taxa
  left_join(inv_surveys) %>% # add plant survey info
  mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
         GSYear = case_when(is.na(SurveyTreatDays) | SurveyTreatDays >= 7 ~ GSYear, # no survey/survey after treatment -> keep year
                            SurveyTreatDays < 7 ~ GSYear + 1)) %>% # survey before or very soon after treatment -> move treatment to next year
  group_by(AreaOfInterestID, AreaOfInterest, PermanentID, Species, GSYear) %>%
  summarize(MaxTreatmentDate = max(TreatmentDate),
            nTreatmentDays = n_distinct(TreatmentDate),
            TreatmentDate = if_else(nTreatmentDays > 1,
                                    paste(min(TreatmentDate), MaxTreatmentDate, sep = " - "),
                                    as.character(MaxTreatmentDate)),
            MaxTreatmentYear = max(TreatmentYear)) %>%
  ungroup()

# check for duplicates
get_dupes(ctrl6, AreaOfInterestID, Species, GSYear)

# fill in missing years
# calculate days between survey and treatment
ctrl7 <- ctrl6 %>%
  mutate(Treated = 1) %>%
  full_join(ctrl6 %>%
              select(AreaOfInterestID, AreaOfInterest, PermanentID) %>%
              unique() %>%
              expand_grid(GSYear = min(ctrl6$GSYear):max(ctrl6$GSYear)) %>% # started keeping track of all activities in 2011
              expand_grid(Species = unique(ctrl6$Species))) %>%
  mutate(Treated = replace_na(Treated, 0)) %>%
  left_join(inv_surveys) %>%
  mutate(SurveyTreatDays = as.numeric(SurveyDate - MaxTreatmentDate)) %>%
  filter(!(is.na(TreatmentDate) & is.na(SurveyDate))) # remove rows without survey or treatment

# check for duplicates
get_dupes(ctrl7, AreaOfInterestID, Species, GSYear)


#### combine datasets for quality analysis ####

# this was the old method for processing data. Don't need area treated an all species treatments, so simplified above
# quality needs to be summarized by permanent ID rather than AOI ID

# # all treatments (for quality analyses)
# # fill in missing years
# # summarize by permanent ID and actual year
# qual_all_ctrl <- ctrl5 %>%
#   full_join(ctrl5 %>%
#               select(PermanentID, AreaOfInterest, Area_ha) %>%
#               unique() %>%
#               expand_grid(GSYear = min(ctrl5$GSYear):max(ctrl5$GSYear))) %>% # could add +1 to max, but don't want to include last year anyway (likely incomplete)
#   mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0)) %>% # add zeros for missing years (did not treat)
#   group_by(PermanentID, GSYear, Area_ha) %>%
#   summarize(AreaName = paste(sort(unique(AreaOfInterest)), collapse = "/"),
#             AllAreaTreated_ha = sum(AreaTreated_ha)) %>% # add all treatments within a year
#   ungroup() %>%
#   mutate(AllTreated = if_else(AllAreaTreated_ha > 0, 1, 0))

# # focal treatments (for quality analyses)
# # fill in missing years
# # summarize by permanent ID and actual year
# qual_foc_ctrl <- ctrl5 %>%
#   filter(Species %in% inv_taxa$Species) %>%
#   full_join(ctrl5 %>%
#               select(PermanentID, AreaOfInterest, Area_ha) %>%
#               unique() %>%
#               expand_grid(GSYear = min(ctrl5$GSYear):max(ctrl5$GSYear)) %>% # could add +1 to max, but don't want to include last year anyway (likely incomplete)
#               expand_grid(Species = unique(inv_taxa$Species))) %>% 
#   mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0)) %>% # add zeros for missing years
#   group_by(PermanentID, GSYear, Area_ha, Species) %>%
#   summarize(AreaName = paste(sort(unique(AreaOfInterest)), collapse = "/"),
#             AreaTreated_ha = sum(AreaTreated_ha),  # add all treatments within a year
#             TreatmentDays = case_when(AreaTreated_ha > 0 ~ as.numeric(n_distinct(TreatmentDate)),
#                                       TRUE ~ 0),
#             TreatmentDate = max(TreatmentDate)) %>%
#   ungroup() %>%
#   mutate(Treated = ifelse(AreaTreated_ha > 0, 1, 0)) 

# # combine
# qual_ctrl <- qual_all_ctrl %>%
#   full_join(qual_foc_ctrl)

# # check for redundancy
# qual_ctrl %>%
#   select(PermanentID, GSYear, Species) %>%
#   get_dupes() %>%
#   select(Species, dupe_count) %>%
#   unique()


#### lag intervals ####

# function for cumulative treatment
ctrl_lag_fun <- function(Year, Lag, dat_in){
  
  # change name
  year1 <- Year
  
  # filter dataset
  subdat <- dat_in %>%
    filter(GSYear <= year1 & GSYear > (year1 - Lag)) # average frequency over preceding years up to lag years
  
  # summarize
  outdat <- subdat %>%
    group_by(Location, Species) %>% # summarize over lag years
    summarise(Treated = mean(Treated),
              NYears = n_distinct(GSYear)) %>%
    ungroup()
  
  # return
  return(outdat)
}

# treatments by year and lag
ctrl8 <- ctrl7 %>%
  select(GSYear) %>%
  unique() %>%
  expand_grid(tibble(Lag = 2:6)) %>% # remove repeat row for each species
  mutate(out = pmap(., function(GSYear, Lag)  # summarizes for each GS, Lag, PermID, and Sp
    ctrl_lag_fun(Year = GSYear, Lag = Lag, 
                 dat_in = ctrl7 %>%
                   mutate(Location = AreaOfInterestID)))) %>%
  unnest(cols = out) %>%
  filter(NYears == Lag) %>% # all years for a lag must be available
  select(-NYears) %>%
  pivot_wider(names_from = Lag,
              values_from = Treated,
              names_glue = "Lag{Lag}Treated") %>% # make treatments wide by lag
  rename(AreaOfInterestID = Location) %>%
  full_join(ctrl7)

# # treatments by year and lag
# qual_ctrl2 <- qual_ctrl %>%
#   select(GSYear) %>%
#   unique() %>%
#   expand_grid(tibble(Lag = 2:6)) %>% # remove repeat row for each species
#   mutate(out = pmap(., function(GSYear, Lag)  # summarizes for each GS, Lag, PermID, and Sp
#     ctrl_lag_fun(Year = GSYear, Lag = Lag, 
#                  dat_in = qual_ctrl %>%
#                    mutate(Location = PermanentID)))) %>%
#   unnest(cols = out) %>%
#   filter(NYears == Lag) %>% # all years for a lag must be available
#   select(-NYears) %>%
#   pivot_wider(names_from = Lag,
#               values_from = c(AllTreated, Treated),
#               names_glue = "Lag{Lag}{.value}") %>% # make treatments wide by lag
#   rename(PermanentID = Location) %>%
#   full_join(qual_ctrl)


#### years since treatment ####

# max lag
inv_max_lag <- max(ctrl8$GSYear) - min(ctrl8$GSYear)
# qual_max_lag <- max(qual_ctrl2$GSYear) - min(qual_ctrl2$GSYear)

# add column for last treatment
ctrl9 <- ctrl8 %>%
  arrange(AreaOfInterestID, Species, GSYear) %>%
  group_by(AreaOfInterestID, Species) %>%
  mutate(LastTreatment = if_else(lag(Treated, n = inv_max_lag) > 0,
                                 inv_max_lag,
                                 NA_real_))

# qual_ctrl3 <- qual_ctrl2 %>%
#   arrange(PermanentID, Species, GSYear) %>%
#   group_by(PermanentID, Species) %>%
#   mutate(LastTreatment = if_else(lag(Treated, n = qual_max_lag) > 0,
#                                  qual_max_lag,
#                                  NA_real_))

# cycle through lags
for(i in (inv_max_lag - 1):0) {
  
  ctrl9 <- ctrl9 %>%
    mutate(LastTreatment = if_else(lag(Treated, n = i) > 0,
                                   as.numeric(i),
                                   LastTreatment))
    
}

# for(i in (qual_max_lag - 1):0) {
#   
#   qual_ctrl3 <- qual_ctrl3 %>%
#     mutate(LastTreatment = if_else(lag(Treated, n = i) > 0,
#                                    as.numeric(i),
#                                    LastTreatment))
#   
# }

# ungroup
ctrl10 <- ctrl9 %>%
  ungroup() %>%
  mutate(RecentTreatment = 1 / (LastTreatment + 1),
         RecentTreatment = replace_na(RecentTreatment, 0)) # no record of treatment

# qual_ctrl4 <- qual_ctrl3 %>%
#   ungroup() %>%
#   mutate(RecentTreatment = 1 / (LastTreatment + 1),
#          RecentTreatment = replace_na(RecentTreatment, 0))

# check that it worked (use hydrilla because it's managed the most)
aois <- unique(ctrl10$AreaOfInterestID)

for(i in aois) {
  
  print(ctrl10 %>%
          filter(Species == "Hydrilla verticillata" & AreaOfInterestID == i) %>%
          ggplot(aes(x = GSYear, color = AreaOfInterestID)) +
          geom_line(aes(y = Treated)) +
          geom_point(aes(y = LastTreatment)) +
          theme_bw() +
          theme(legend.position = "none"))
  
}

ggplot(ctrl10, aes(x = LastTreatment)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ Species)


#### inv ctrl output ####
write_csv(ctrl10, "intermediate-data/FWC_invasive_control_formatted.csv")
# write_csv(qual_ctrl4, "intermediate-data/FWC_quality_control_formatted.csv")
