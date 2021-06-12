
#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)

# import data
ctrl_old <- read_csv("intermediate-data/FWC_control_old_formatted.csv")
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
dayDat <- read_csv("intermediate-data/LakeO_day_data_for_model.csv")

# load models
load("output/lakeO_intraannual_veg_change_mod.rda")

# figure settings
def_theme <- theme_bw() +
  theme(axis.text = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.text = element_text(size = 14, color="black"),
        strip.background = element_blank())


#### lake O parameters ####

# extract coefficients
lakeO_beta1 <- coef(summary(lake0_mod))[2, "Estimate"] # days
lakeO_beta2 <- coef(summary(lake0_mod))[3, "Estimate"] # days^2

# date of max abundance (-b/2a)
lakeO_days <- -lakeO_beta1 / (2 * lakeO_beta2)


#### edit FWC plant data ####

# Perm IDs per AOI
plant_fwc %>%
  group_by(AreaOfInterestID) %>%
  summarise(IDs = length(unique(PermanentID))) %>%
  filter(IDs > 1)
# one Perm ID per AOI

# function to remove duplicates
rem_dups_fun <- function(dat) {
  
  if (nrow(dat) == 1){
    dat2 <- dat
  } else if (var(dat$AreaCovered_ha) > 0){
    dat2 <- dat %>%
      filter(AreaCovered_ha == max(AreaCovered_ha)) # choose survey with maximum area
  } else if (var(dat$AreaCovered_ha) == 0){
    dat2 <- dat %>%
      mutate(DaysDiff = abs(lakeO_days - Days)) %>%
      filter(DaysDiff == min(DaysDiff)) %>% # choose survey closes to max. abundance time
      select(-DaysDiff)
  }
  
  return(dat2)
}

plant_fwc2 <- plant_fwc %>% # start with all surveys
  select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate) %>%
  unique() %>% # one row per survey
  expand_grid(tibble(SpeciesName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"))) %>% # one row per species per survey
  full_join(plant_fwc %>% # add invasive plant information
              filter(SpeciesName %in% c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes")) %>%
              select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, SpeciesName, SpeciesAcres)) %>%
  mutate(SpeciesAcres = replace_na(SpeciesAcres, 0), # cover 0 when it wasn't in a survey
         Area_ha = ShapeArea * 100, # convert lake area from km-squared to hectares
         AreaCovered_ha = SpeciesAcres * 0.405, # convert plant cover from acres to hectares
         SurveyMonth = month(SurveyDate),
         SurveyDay = day(SurveyDate),
         SurveyYear = case_when(SurveyMonth >= 4 ~ year(SurveyDate),
                                SurveyMonth < 4 ~ year(SurveyDate) - 1), # assume growing season starts in April
         MonthDay = case_when(SurveyMonth >= 4 ~ as.Date(paste("2020", SurveyMonth, SurveyDay, sep = "-")), # start "year" in April (2020/2021 are arbitrary)
                              SurveyMonth < 4 ~ as.Date(paste("2021", SurveyMonth, SurveyDay, sep = "-")))) %>%
  left_join(dayDat) %>% # add standardized days (proportion between April 1 and March 31)
  nest(data = c(SurveyDate, SpeciesAcres, AreaCovered_ha, SurveyMonth, SurveyDay, MonthDay, Days)) %>% # find duplicates within lake, year, and species
  mutate(newdata = map(data, ~rem_dups_fun(.))) %>% # remove duplicates
  select(-data) %>% # removes 132 rows of data
  unnest(newdata) %>%
  mutate(AreaChangeSD = lakeO_beta1 * (lakeO_days-Days) + lakeO_beta2 * (lakeO_days^2 - Days^2),
         CommonName = case_when(SpeciesName == "Eichhornia crassipes" ~ "Water hyacinth", 
                                SpeciesName == "Hydrilla verticillata" ~ "Hydrilla", 
                                SpeciesName == "Pistia stratiotes" ~ "Water lettuce")) %>% # calculate the number of sd's to change to get est. max abundance
  group_by(AreaOfInterestID, SpeciesName) %>%
  mutate(EstAreaCoveredRaw_ha = AreaCovered_ha + AreaChangeSD * sd(AreaCovered_ha), # calculate est. max abundance
         EstAreaCovered_ha = case_when(EstAreaCoveredRaw_ha > Area_ha ~ Area_ha, # reduce areas covered to total area
                                       TRUE ~ EstAreaCoveredRaw_ha),
         log_EstAreaCovered_ha = log(EstAreaCovered_ha + 1e-9),  # log-transform (min EstAreaCovered_ha (> 0) = 2x10^-8)
         PropCovered = EstAreaCovered_ha / Area_ha,
         log_PropCovered = log((EstAreaCovered_ha + 1e-9) / Area_ha),
         Detected = as.numeric(sum(SpeciesAcres) > 0)) %>%
  ungroup()

# duplicates in same year
plant_fwc2 %>%
  group_by(AreaOfInterestID, SurveyYear, SpeciesName) %>%
  count() %>%
  filter(n > 1)

# save data
write_csv(plant_fwc2, "intermediate-data/FWC_hydrilla_pistia_eichhornia_survey_formatted.csv")


#### edit control data ####

# non-herbicide methods (from herbicide_initial_visualizations)
non_herb <- c("Mechanical Harvester", 
              "Snagging (tree removal)", 
              "Aquatic Dye (for shading)", 
              "Grass Carp", "Hand Removal", 
              "Mechanical (Other)", 
              "Mechanical Shredder", 
              "Prescribed Fire")

# duplicates per year/date?
ctrl_old %>%
  filter(Species %in% c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)") & TotalAcres > 0) %>%
  group_by(AreaOfInterestID, Year, Species) %>%
  summarise(apps = n(),
            acres = paste(TotalAcres, collapse = ", ")) %>%
  filter(apps > 1)
# yes, 106
# different acres in same year, probably different events

ctrl_new %>%
  filter(Species %in% c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)") & 
           TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
  group_by(AreaOfInterestID, BeginDate, Species) %>%
  summarise(apps = n(),
            acres = paste(TotalAcres, collapse = ", "),
            method = paste(ControlMethod, collapse = ", "),
            id = paste(TreatmentID, collapse = ", ")) %>%
  filter(apps > 1)
# yes, 5866
# for treatments involving more than one herbicide, separate lines for each herbicide type, but they share a treatment ID and acres
ctrl_new %>%
  group_by(TreatmentID) %>%
  summarise(acres = length(unique(TotalAcres)),
            species = length(unique(Species)),
            dates = length(unique(BeginDate))) %>%
  filter(acres > 1 | species > 1 | dates > 1)
# each treatment ID is a unique application for each species, area, and date
sum(is.na(ctrl_new$TreatmentID))

# old herbicide data
ctrl_old2 <- ctrl_old %>%
  filter(Species %in% c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)") & TotalAcres > 0) %>%
  mutate(AreaTreated_ha= TotalAcres * 0.405,
         Area_ha = ShapeArea * 100,
         TreatmentMethod = "unknown",
         TreatmentMonth = 12,
         TreatmentDate = as.Date(paste(as.character(Year), "-12-01", sep = "")), # no date given for these surveys
         TreatmentID = paste("old", Year, substr(Species, 1, 1), TotalAcres, sep = "_")) %>%
  select(AreaOfInterestID, PermanentID, Year, Species, Area_ha, AreaTreated_ha, TreatmentMethod, TreatmentMonth, TreatmentDate, TreatmentID) %>%
  rename(TreatmentYear = Year)

# new herbicide data
ctrl_new2 <- ctrl_new %>%
  filter(Species %in% c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)") & 
           TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>% # herbicide control only
  group_by(AreaOfInterestID, PermanentID, Species, BeginDate, TreatmentID, TotalAcres, ShapeArea) %>%  # captures area treated for an event without duplication due to multiple herbicides
  mutate(AreaTreated_ha = TotalAcres * 0.405,
         Area_ha = ShapeArea * 100,
         TreatmentMethod = paste(ControlMethod, collapse = ", "),
         TreatmentYear = year(BeginDate),
         TreatmentMonth = month(BeginDate),
         TreatmentID = as.character(TreatmentID)) %>%
  ungroup() %>%
  select(AreaOfInterestID, PermanentID, TreatmentYear, Species, Area_ha, AreaTreated_ha, TreatmentMethod, TreatmentMonth, BeginDate, TreatmentID) %>%
  rename(TreatmentDate = BeginDate)

# combine herbicide data
ctrl <- ctrl_old2 %>%
  full_join(ctrl_new2) %>%
  full_join(tibble(Species = c("Hydrilla verticillata", rep("Floating Plants (Eichhornia and Pistia)", 2)), # double each row that has floating plants - one row for each species
                   SpeciesName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"))) %>%
  mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha, # make treatment area size of lake if it exceeds it
                                    TRUE ~ AreaTreated_ha)) 

# add zero for missing years (assume no application)
ctrl2 <- ctrl_old %>% # one row for every lake
  mutate(Area_ha = ShapeArea * 100) %>%
  select(AreaOfInterestID, PermanentID, Area_ha) %>%
  unique() %>%
  full_join(ctrl_new %>%
              mutate(Area_ha = ShapeArea * 100) %>%
              select(AreaOfInterestID, PermanentID, Area_ha) %>%
              unique()) %>%
  expand_grid(tibble(TreatmentYear = min(ctrl_old2$TreatmentYear):max(ctrl_new2$TreatmentYear)) %>% # one row per year and species
                expand_grid(tibble(Species = c("Hydrilla verticillata", rep("Floating Plants (Eichhornia and Pistia)", 2)),
                                   SpeciesName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes")))) %>% 
  full_join(ctrl) %>%
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0), # make area treated zero if the year wasn't included
         TreatmentEvent = TreatmentID, # NA when no treatment applied
         TreatmentID = case_when(is.na(TreatmentID) ~ paste("none", TreatmentYear, sep = "_"),
                                 TRUE ~ TreatmentID))

# save data
write_csv(ctrl2, "intermediate-data/FWC_hydrilla_pistia_eichhornia_herbicide_formatted.csv")


#### combine FWC plant and ctrl data ####

plant_ctrl <- ctrl2 %>%
  mutate(SurveyYear = TreatmentYear) %>% # same year (for difference calcs of abundance, which use lead)
  group_by(AreaOfInterestID, PermanentID, SurveyYear, Area_ha, Species, SpeciesName) %>%
  summarise(TotalAreaTreated_ha = sum(AreaTreated_ha)) %>%
  ungroup() %>%
  full_join(ctrl2 %>%
              mutate(SurveyYear = TreatmentYear + 1) %>% # lag 1 year
              group_by(AreaOfInterestID, PermanentID, SurveyYear, Species, SpeciesName) %>%
              summarise(TotalAreaTreatedLag1_ha = sum(AreaTreated_ha)) %>%
              ungroup()) %>%
  full_join(ctrl2 %>%
              mutate(SurveyYear = TreatmentYear + 2) %>% # lag 2 years
              group_by(AreaOfInterestID, PermanentID, SurveyYear, Species, SpeciesName) %>%
              summarise(TotalAreaTreatedLag2_ha = sum(AreaTreated_ha)) %>%
              ungroup()) %>%
  full_join(ctrl2 %>%
              mutate(SurveyYear = TreatmentYear + 3) %>% # lag 3 years
              group_by(AreaOfInterestID, PermanentID, SurveyYear, Species, SpeciesName) %>%
              summarise(TotalAreaTreatedLag3_ha = sum(AreaTreated_ha)) %>%
              ungroup()) %>%
  mutate(TotalPropTreated = TotalAreaTreated_ha / Area_ha,
         TotalPropTreatedLag1 = TotalAreaTreatedLag1_ha / Area_ha,
         TotalPropTreatedLag2 = TotalAreaTreatedLag2_ha / Area_ha,
         TotalPropTreatedLag3 = TotalAreaTreatedLag3_ha / Area_ha,
         Treated = ifelse(TotalAreaTreated_ha > 0, 1, 0),
         TreatedLag1 = ifelse(TotalAreaTreatedLag1_ha > 0, 1, 0),
         TreatedLag2 = ifelse(TotalAreaTreatedLag2_ha > 0, 1, 0),
         TreatedLag3 = ifelse(TotalAreaTreatedLag3_ha > 0, 1, 0)) %>%
  left_join(plant_fwc2) %>% # only include years with treatment info
  group_by(AreaOfInterestID, PermanentID, SpeciesName) %>%
  arrange(SurveyYear) %>%
  mutate(log_NextPropCovered = lead(log_PropCovered), # includes correction for zero abundance
         log_DiffPropCovered = log_NextPropCovered - log_PropCovered,
         NextPropCovered = lead(PropCovered),
         RatioPropCovered = NextPropCovered/PropCovered,
         Detected = as.numeric(sum(Detected, na.rm = T) > 0)) %>% # fills in NA's for detected
  ungroup()

# missing data
plant_ctrl %>%
  filter(is.na(PropCovered) | is.na(NextPropCovered))
# > 10,000 (~ 1/3 data)

# treatment without detection?
treat_no_detect <- plant_ctrl %>%
  group_by(Species, AreaOfInterestID) %>%
  mutate(GroupDetected = as.numeric(sum(Detected) > 0)) %>% # were either of the floating plants detected?
  ungroup() %>%
  filter(GroupDetected == 0 & TotalPropTreated > 0)
# 229 examples

treat_no_detect %>%
  ggplot(aes(x = TotalPropTreated)) +
  geom_histogram(binwidth = 0.1) +
  facet_wrap(~ Species)
# 292 instances of herbicide applied to a treat a species without it being detected by surveys

# look at some of these
treat_no_detect %>%
  full_join(tibble(Species = c("Hydrilla verticillata", rep("Floating Plants (Eichhornia and Pistia)", 2)), # double each row that has floating plants - one row for each species
                   SpeciesName = c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes"))) %>%
  rename(TreatmentYear = SurveyYear) %>%
  select(SpeciesName, AreaOfInterestID, TreatmentYear) %>%
  unique() %>%
  inner_join(plant_fwc2 %>%
               select(AreaOfInterestID, SurveyYear, SpeciesName, SpeciesAcres, PropCovered, Detected)) %>%
  ggplot(aes(x = SurveyYear, y = SpeciesAcres)) +
  geom_line(aes(color = as.factor(AreaOfInterestID)), show.legend = F) +
  geom_vline(aes(xintercept = TreatmentYear, color = as.factor(AreaOfInterestID)), show.legend = F) +
  facet_wrap(~ SpeciesName)

# remove missing data for difference analysis
plant_ctrl2 <- plant_ctrl %>%
  filter(!is.na(PropCovered) & !is.na(NextPropCovered) & Detected == 1)
# Detected == 1: plant must have been detected in the lake at some point


#### summary stats ####

# overall dataset
plant_ctrl2 %>%
  summarise(Lakes = length(unique(AreaOfInterestID)),
            Years = length(unique(SurveyYear)))

# detected lakes
lake_sum <- plant_ctrl2 %>%
  group_by(CommonName) %>%
  summarise(Lakes = length(unique(AreaOfInterestID)),
            YearLakes = n()) %>%
  mutate(NameLakes = paste(CommonName, "\n(", Lakes, " lakes, ", YearLakes, " data points)", sep = ""))

# year ranges
pdf("output/year_ranges_herbicide_analysis.pdf", width = 11.5, height = 5.75)
plant_ctrl2 %>%
  group_by(CommonName, AreaOfInterestID) %>%
  summarise(Years = max(SurveyYear) - min(SurveyYear) + 1) %>%
  ungroup() %>%
  left_join(lake_sum) %>%
  ggplot(aes(x = Years)) +
  geom_bar() +
  facet_wrap(~ NameLakes) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  labs(x = "Years of data", y = "Number of lakes") +
  def_theme
dev.off()

# population sizes
prop_sum <- plant_ctrl2 %>%
  group_by(CommonName) %>%
  summarise(mean = mean(NextPropCovered))

pdf("output/pop_sizes_herbicide_analysis.pdf", width = 11.5, height = 5.75)
plant_ctrl2 %>%
  ggplot(aes(NextPropCovered)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(data = prop_sum, aes(xintercept = mean), color = "blue", linetype = "dashed") +
  geom_text(data = prop_sum, aes(x =  mean, label = paste("mean = ", as.character(round(mean, 3)), sep = "")), 
            y = 5800, hjust = -0.05, size = 4, color = "blue") +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  labs(x = "Proportion of lake inhabited", y = "Data points") +
  def_theme
dev.off()

pdf("output/pop_changes_herbicide_analysis.pdf", width = 11.5, height = 5.75)
plant_ctrl2 %>%
  ggplot(aes(log_DiffPropCovered)) +
  geom_histogram() +
  facet_wrap(~ CommonName) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf)  +
  labs(x = expression(paste("ln(", N[t+1], "/", N[t], ")", sep = "")),
       y = "Data points") +
  def_theme
dev.off()