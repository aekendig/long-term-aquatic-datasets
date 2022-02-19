#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(plotly)
library(tidyverse)
library(fixest) # FE models
library(modelsummary)
library(inspectdf) # for inspect_cor

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")
nat_plant <- read_csv("intermediate-data/FWC_native_plant_abundance_formatted.csv")


#### edit native plant data ####

# Avg prop covered columns
inv_avg_cols <- tibble(cols = colnames(inv_plant)) %>%
  filter(str_detect(cols, "AvgProp") == T) %>%
  pull(cols)

# make long by lag
# make wide by inv. plant species
inv_plant2 <- inv_plant %>%
  filter(CommonName != "Wild taro") %>%
  select(PermanentID, GSYear, CommonName, all_of(inv_avg_cols)) %>%
  pivot_longer(cols = all_of(inv_avg_cols),
               names_to = "Lag",
               values_to = "AvgPropCovered") %>%
  filter(!is.na(AvgPropCovered)) %>% # remove missing data
  mutate(CommonName = fct_recode(CommonName,
                                 "WaterHyacinth" = "Water hyacinth",
                                 "WaterLettuce" = "Water lettuce",
                                 "AlligatorWeed" = "Alligator weed",
                                 "CubanBulrush" = "Cuban bulrush",
                                 "WaterFern" = "Water fern",
                                 "ParaGrass" = "Para grass"),
         Lag = as.numeric(str_sub(Lag, 4, 4)),
         AvgPercCovered = AvgPropCovered * 100) %>%
  select(-AvgPropCovered) %>%
  pivot_wider(names_from = CommonName,
              values_from = AvgPercCovered,
              names_glue = "{CommonName}_AvgPercCovered") %>%
  mutate(Floating_AvgPercCovered = WaterHyacinth_AvgPercCovered + WaterLettuce_AvgPercCovered)

# Avg prop covered columns
inv_ctrl_cols <- tibble(cols = colnames(inv_ctrl)) %>%
  filter(str_detect(cols, "Treated") == T & 
           str_detect(cols, "Lag") == T & 
           str_detect(cols, "PropTreated") == F & 
           str_detect(cols, "All") == F) %>%
  pull(cols)

# make long by lag
# make wide by control target
inv_ctrl2 <- inv_ctrl %>%
  filter(Species != "Colocasia esculenta") %>%
  select(PermanentID, Species, GSYear, all_of(inv_ctrl_cols)) %>% 
  unique() %>% # remove duplication of floating plant treatment
  pivot_longer(cols = all_of(inv_ctrl_cols),
               names_to = "Lag",
               values_to = "TreatFreq") %>%
  filter(!is.na(TreatFreq)) %>% # remove missing data
  mutate(Species = fct_recode(Species,
                              "Floating" = "Floating Plants (Eichhornia and Pistia)",
                              "Hydrilla" = "Hydrilla verticillata",
                              "AlligatorWeed" = "Alternanthera philoxeroides",
                              "CubanBulrush" = "Cyperus blepharoleptos",
                              "Torpedograss" = "Panicum repens",
                              "WaterFern" = "Salvinia minima",
                              "ParaGrass" = "Urochloa mutica"),
         Lag = as.numeric(str_sub(Lag, 4, 4))) %>%
  pivot_wider(names_from = Species,
              values_from = TreatFreq,
              names_glue = "{Species}_TreatFreq")


#### combine data ####

# filter invasion for all lags
inv_plant3 <- inv_plant2 %>%
  filter(Lag == 6) %>%
  select(PermanentID, GSYear) %>%
  inner_join(inv_plant2)

# filter control for all lags
inv_ctrl3 <- inv_ctrl2 %>%
  filter(Lag == 6) %>%
  select(PermanentID, GSYear) %>%
  inner_join(inv_ctrl2)

# native and invasive data
nat_inv <- nat_plant %>%
  filter(!is.na(Lag1LogRatioCovered)) %>%
  select(TaxonName, CommonName, PermanentID, GSYear, PropCovered,
         Lag1InitPercCovered, Lag1LogRatioCovered, Lag1MinSurveyorExperience) %>%
  inner_join(inv_plant3)

# combine native, invasive, control
# remove missing data
nat_dat <- nat_inv %>%
  inner_join(inv_ctrl3)


#### initial visualizations ####

# surveys
nat_dat %>%
  select(PermanentID, GSYear) %>%
  unique() %>%
  ggplot(aes(x = GSYear)) +
  geom_bar()

# native abundance over space and time
nat_dat %>%
  filter(Lag == 1) %>%
  ggplot(aes(x = GSYear, y = PropCovered, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")
# not enough data to analyze

nat_inv %>%
  filter(Lag == 1) %>%
  ggplot(aes(x = GSYear, y = PropCovered, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")
# nat_plant is missing every other year between 1986 and 1994
# nat_inv then omits 1986-1994 (no previous year to build growth rate)
# inv_plant omits pre-1988 because it needs 6 years of data
# constrains nat_inv to 1995+