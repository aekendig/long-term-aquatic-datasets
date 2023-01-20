#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
# qual_ctrl <- read_csv("intermediate-data/FWC_quality_control_formatted.csv")
ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv") # this one is summarized by GS year instead of treatment year
# using the GS year one because summarizing quality values by year led to a max of 9 months, I think because GS year of the invasive species data were mismatched with treatment year of the control data
lw <- read_csv("intermediate-data/LW_quality_formatted.csv")
atlas <- read_csv("intermediate-data/water_atlas_quality_formatted.csv")

# import functions
source("code/generic-functions/continuous_time_interval.R")


#### invasive plant data ####

# combine water hyacinth and lettuce percent covered
floating_cover <- inv_plant %>%
  filter(CommonName %in% c("Water hyacinth", "Water lettuce")) %>%
  group_by(PermanentID, GSYear) %>%
  summarize(across(.cols = ends_with ("PropCovered"), sum),
            InitPercCovered = sum(InitPercCovered),
            SpeciesAcres = sum(SpeciesAcres),
            EstAreaCoveredRaw_ha = sum(EstAreaCoveredRaw_ha)) %>%
  ungroup() %>%
  mutate(across(.cols = ends_with ("PropCovered"), ~ if_else(.x > 1, 1, .x)),
         InitPercCovered = if_else(InitPercCovered > 1, 1, InitPercCovered),
         CommonName = "floating plants")

# add floating cover
inv_plant2 <- inv_plant %>%
  full_join(floating_cover) %>%
  mutate(across(ends_with("PropCovered"), ~ .x * 100), # change proportion to percent
         Lag2APCsq = Lag2AvgPropCovered^2) %>% # square perc covered
  rename_with(str_replace, pattern = "PropCovered", replacement = "PercCovered") %>%
  select(CommonName, PermanentID, GSYear, ends_with("PercCovered"), Lag2APCsq, SpeciesAcres, EstAreaCoveredRaw_ha)

# Species names
inv_taxa <- tibble(Species = c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)", "Panicum repens", "Urochloa mutica", "Cyperus blepharoleptos"),
                   CommonName = c("Hydrilla", "floating plants", "Torpedograss", "Para grass", "Cuban bulrush"))


#### management data ####

# use "Species" and "unique" to get floating plants combined
ctrl2 <- ctrl %>%
  select(Species, PermanentID, GSYear, LastTreatment, RecentTreatment, ends_with("Treated")) %>%
  unique() %>%
  left_join(inv_taxa)

#### identify waterbodies to use ####

# has the plant ever been established?
# by using EstAreaCoveredRaw_ha, there had to be more than one year per permanentID
perm_plant <- inv_plant2 %>%
  group_by(PermanentID, CommonName) %>%
  summarize(Established = as.numeric(sum(EstAreaCoveredRaw_ha) > 0)) %>%
  ungroup() %>%
  filter(Established > 0)

# has the plant ever been treated?
perm_ctrl <- ctrl2 %>%
  group_by(PermanentID, Species) %>%
  summarize(Treatments = as.numeric(sum(Treated, na.rm = T) > 0)) %>%
  ungroup() %>%
  filter(Treatments > 0) %>%
  left_join(inv_taxa)

# waterbodies that have the species present and been managed at least once
perm_plant_ctrl <- inner_join(perm_plant, perm_ctrl) %>%
  select(PermanentID, CommonName)


#### quality data ####

# choose/rename columns
lw2 <- lw %>%
  select(PermanentID, GSYear, Month, QualityMetric, QualityValue) %>%
  filter(QualityMetric %in% c("CHL_ug_L", "TN_ug_L", "TP_ug_L", "Secchi_ft"))

atlas2 <- atlas %>%
  mutate(GSYear = case_when(Month >= 4 ~ Year, # match to plant survey data from the same growing season and treatment from previous year
                            Month < 4 ~ Year - 1)) %>%
  select(PermanentID, GSYear, Month, Parameter, Result_Value) %>%
  mutate(Parameter = fct_recode(Parameter,
                                "TN_ug_L" = "tn_ugl",
                                "Secchi_ft" = "secchi_ft",
                                "CHL_ug_L" = "chla_ugl",
                                "TP_ug_L" = "tp_ugl") %>%
           as.character()) %>%
  rename(QualityMetric = Parameter,
         QualityValue = Result_Value)

# combine data from lake watch and water atlas
# select lakes with invasive plant
qual <- lw2 %>%
  full_join(atlas2) %>%
  group_by(PermanentID, GSYear, Month, QualityMetric) %>%
  summarize(QualityValue = mean(QualityValue)) %>%
  ungroup()

# month distribution
qual %>%
  group_by(QualityMetric, GSYear, Month) %>%
  summarize(WB = n_distinct(PermanentID)) %>%
  ungroup() %>%
  ggplot(aes(x = Month, y = WB, color = as.factor(GSYear))) +
  geom_point() +
  geom_line() +
  facet_wrap(~ QualityMetric) +
  theme(legend.position = "none")

# quarters
quarters <- tibble(Month = 1:12,
                   Quarter = rep(1:4, each = 3))

# quarter distribution
qual %>%
  left_join(quarters) %>%
  distinct(QualityMetric, GSYear, Quarter, PermanentID) %>%
  group_by(QualityMetric, GSYear, Quarter) %>%
  summarize(WB = n_distinct(PermanentID)) %>%
  ungroup() %>%
  ggplot(aes(x = Quarter, y = WB, color = as.factor(GSYear))) +
  geom_point() +
  geom_line() +
  facet_wrap(~ QualityMetric) +
  theme(legend.position = "none")

# annual metrics
qual_ann <- qual %>%
  left_join(quarters) %>%
  group_by(PermanentID, GSYear, QualityMetric) %>%
  mutate(Quarters = n_distinct(Quarter)) %>%
  ungroup() %>%
  filter(Quarters == 4) %>% # require one month per quarter
  group_by(PermanentID, GSYear, QualityMetric) %>%
  summarize(Months = n_distinct(Month),
            MonthsOmitted = setdiff(1:12, c(Month)) %>%
              paste(collapse = ", "),
            QualityMean = mean(QualityValue),
            QualityMax = max(QualityValue),
            QualityMin = min(QualityValue),
            QualityCV = sd(QualityValue)/mean(QualityValue)) %>%
  ungroup() %>%
  filter(Months > 8) %>% # require 8 months
  group_by(PermanentID, GSYear) %>%
  mutate(Metrics = n_distinct(QualityMetric)) %>%
  ungroup() %>%
  filter(Metrics == 4) # require all four metrics


#### combine data and select waterbodies/years ####

# check recent treatment values (to modify missing values below)
ctrl2 %>%
  filter(is.na(RecentTreatment)) %>%
  distinct(LastTreatment)
# no last treatment values 
# (because dataset started with no treatment and continued without treatment until one was applied)

range(ctrl2$RecentTreatment, na.rm = T)

# add invasive plant and control data
qual_ann2 <- qual_ann %>%
  inner_join(inv_plant2) %>% # select waterbodies and years in both datasets
  inner_join(ctrl2) %>%
  left_join(perm_plant) %>%
  mutate(LogQualityMean = log(QualityMean),
         LogQualityMax = log(QualityMax),
         LogQualityMin = log(QualityMin),
         LogQualityCV = log(QualityCV),
         RecentTreatment = replace_na(RecentTreatment, 0),
         Established = replace_na(Established, 0))

# select waterbodies that have had species and at least one year of management
qual_inv <- qual_ann2 %>%
  inner_join(perm_plant_ctrl)

# check for issues with log-transformation
filter(qual_ann2, is.na(LogQualityMean) | LogQualityMean == -Inf)
filter(qual_ann2, is.na(LogQualityMax) | LogQualityMax == -Inf)
filter(qual_ann2, is.na(LogQualityMin) | LogQualityMin == -Inf) %>% distinct(QualityMetric) # all chlorophyll
filter(qual_ann2, is.na(LogQualityCV) | LogQualityCV == -Inf)

# look at minimum chlorophyll
min_chl <- qual_ann2 %>%
  filter(QualityMetric == "CHL_ug_L" & QualityMin > 0) %>%
  pull(QualityMin) %>%
  min()

# update QualityMin
qual_ann3 <- qual_ann2 %>%
  mutate(LogQualityMin = if_else(QualityMetric == "CHL_ug_L", log(QualityMin + min_chl),
                                 log(QualityMin)))

qual_inv2 <- qual_inv %>%
  mutate(LogQualityMin = if_else(QualityMetric == "CHL_ug_L", log(QualityMin + min_chl),
                                 log(QualityMin)))

# list of lakes, years, and SpeciesAcres (to detect NAs)
qual_perm_year <- qual_ann3 %>%
  distinct(PermanentID, GSYear, CommonName, SpeciesAcres)

# lake/year coverage
qual_perm_year2 <- qual_perm_year %>%
  distinct(PermanentID, GSYear, CommonName) %>% # remove acres, shouldn't change row number
  group_by(CommonName) %>%
  mutate(Perm_tot = n_distinct(PermanentID), # waterbodies per species
         Year_tot = n_distinct(GSYear)) %>% # years per species
  ungroup() %>%
  group_by(CommonName, PermanentID) %>%
  mutate(Years_per_perm = n_distinct(GSYear)) %>% # years per waterbody per species
  ungroup() %>%
  group_by(CommonName, GSYear) %>%
  mutate(Perms_per_year = n_distinct(PermanentID)) %>% # waterbodies per year per species
  ungroup() %>%
  mutate(Prop_years_per_perm = Years_per_perm / Year_tot,
         Prop_perms_per_year = Perms_per_year / Perm_tot)

# look at distiribution
qual_perm_year2 %>%
  distinct(CommonName, PermanentID, Prop_years_per_perm) %>%
  ggplot(aes(x = Prop_years_per_perm)) +
  geom_histogram(binwidth = 0.1) +
  facet_wrap(~ CommonName, scales = "free")

qual_perm_year2 %>%
  distinct(CommonName, GSYear, Prop_perms_per_year) %>%
  ggplot(aes(x = Prop_perms_per_year)) +
  geom_histogram(binwidth = 0.1) +
  facet_wrap(~ CommonName, scales = "free")

# try different cutoffs
qual_perm_year2 %>%
  expand_grid(MinProp = seq(0, 1, by = 0.1)) %>%
  filter(Prop_years_per_perm >= MinProp &
           Prop_perms_per_year >= MinProp) %>%
  group_by(CommonName, MinProp) %>%
  summarize(Years = n_distinct(GSYear),
            WBs = n_distinct(PermanentID)) %>%
  ungroup() %>%
  ggplot(aes(x = Years, y = WBs, color = as.factor(MinProp))) +
  geom_point() +
  facet_wrap(~ CommonName, scales = "free")

# visualize data availability
ggplot(qual_perm_year2, aes(x = GSYear, y = PermanentID)) +
  geom_line() +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none",
        axis.text.y = element_blank())

# look at months omitted
sort(unique(qual_ann3$MonthsOmitted))


#### check values ####

# check for duplicates
qual_ann3 %>%
  get_dupes(GSYear, PermanentID, QualityMetric, CommonName)

# check for outliers
qual_ann3 %>%
  distinct(GSYear, PermanentID, QualityMetric, QualityMean) %>% # remove duplication by inv taxa
  ggplot(aes(x = GSYear, y = QualityMean, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ QualityMetric, scales = "free_y") +
  theme(legend.position = "none")

# select only quality data and make wide
qual_ann_wide <- qual_ann3 %>%
  distinct(GSYear, PermanentID, QualityMetric, QualityMean) %>%
  mutate(LogQualityMean = log(QualityMean)) %>%
  select(-QualityMean) %>%
  pivot_wider(names_from = QualityMetric,
              values_from = LogQualityMean)

# look at correlations
qual_ann_wide %>%
  select(-c(PermanentID, GSYear)) %>%
  GGally::ggpairs()


#### save data ####

# save data
write_csv(qual_ann3, "intermediate-data/water_quality_data_formatted.csv")
write_csv(qual_inv2, "intermediate-data/water_quality_invaded_data_formatted.csv")


#### waterbody counts ####
qual_ann3 %>%
  group_by(CommonName, Established) %>%
  summarize(waterbodies = n_distinct(PermanentID)) %>%
  ungroup() %>%
  mutate(Established = fct_recode(as.factor(Established),
                                  "Invaded" = "1",
                                  "Uninvaded" = "0"),
         CommonName = tolower(CommonName)) %>%
  pivot_wider(names_from = Established,
              values_from = waterbodies) %>%
  write_csv("intermediate-data/water_quality_invaded_uninvaded_counts.csv")

# previously had uninvaded waterbodies, but these included waterbodies
# that weren't in the plant survey dataset, so it's not actually
# certain that they were uninvaded