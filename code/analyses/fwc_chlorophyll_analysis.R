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
lw_chl <- read_csv("intermediate-data/LW_chlorophyll_formatted.csv")
lwwa_chl <- read_csv("intermediate-data/LW_water_atlas_chlorophyll_formatted.csv")


#### edit data ####

# Avg prop covered columns
inv_avg_cols <- tibble(cols = colnames(inv_plant)) %>%
  filter(str_detect(cols, "AvgProp") == T) %>%
  pull(cols)

# check invasive plant data availability
inv_plant %>%
  inner_join(lwwa_chl %>%
               select(PermanentID, GSYear) %>%
               unique()) %>%
  inner_join(inv_ctrl %>%
               select(PermanentID, GSYear) %>%
               unique()) %>%
  filter(Lag6AvgPropCovered > 0) %>%
  ggplot(aes(x = GSYear, y = Lag6AvgPropCovered, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ CommonName, scales = "free_y") +
  theme(legend.position = "none")
# not enough data for: Alligator weed, water fern
# Cuban bulrush is missing a lot of data before 2013
# wild taro is difficult to distinguish from elephant ear, surveys may be inaccurate

inv_taxa <- inv_plant %>%
  filter(CommonName %in% c("Hydrilla", "Water lettuce", "Water hyacinth",
                           "Torpedograss", "Para grass", "Cuban bulrush")) %>%
  select(CommonName, TaxonName) %>%
  unique()

# make long by lag
# make wide by inv. plant species
inv_plant2 <- inv_plant %>%
  filter(CommonName %in% inv_taxa$CommonName) %>%
  select(PermanentID, GSYear, CommonName, all_of(inv_avg_cols)) %>%
  pivot_longer(cols = all_of(inv_avg_cols),
               names_to = "Lag",
               values_to = "AvgPropCovered") %>%
  filter(!is.na(AvgPropCovered)) %>% # remove missing data
  mutate(CommonName = fct_recode(CommonName,
                                 "WaterHyacinth" = "Water hyacinth",
                                 "WaterLettuce" = "Water lettuce",
                                 "ParaGrass" = "Para grass",
                                 "CubanBulrush" = "Cuban bulrush"),
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
  filter(TaxonName %in% inv_taxa$TaxonName) %>%
  select(PermanentID, Species, GSYear, all_of(inv_ctrl_cols)) %>% 
  unique() %>% # remove duplication of floating plant treatment
  pivot_longer(cols = all_of(inv_ctrl_cols),
               names_to = "Lag",
               values_to = "TreatFreq") %>%
  filter(!is.na(TreatFreq)) %>% # remove missing data
  mutate(Species = fct_recode(Species,
                              "Floating" = "Floating Plants (Eichhornia and Pistia)",
                              "Hydrilla" = "Hydrilla verticillata",
                              "Torpedograss" = "Panicum repens",
                              "ParaGrass" = "Urochloa mutica",
                              "CubanBulrush" = "Cyperus blepharoleptos"),
         Lag = as.numeric(str_sub(Lag, 4, 4))) %>%
  pivot_wider(names_from = Species,
              values_from = TreatFreq,
              names_glue = "{Species}_TreatFreq")

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

# combine chlorophyll, invasive, control
chl_dat <- lwwa_chl %>%
  inner_join(inv_plant3) %>%
  inner_join(inv_ctrl3)

lw_chl_dat <- lw_chl %>%
  inner_join(inv_plant3) %>%
  inner_join(inv_ctrl3)


#### initial visualizations ####

# covariate correlations
chl_dat %>%
  select(ends_with("TreatFreq"), 
         ends_with("AvgPercCovered"),
         PrevValue,
         Lag) %>%
  group_by(Lag) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & corr >= 0.4 &
           !(col_1 %in% c("Floating_AvgPercCovered", "WaterLettuce_AvgPercCovered", "WaterHyacinth_AvgPercCovered") &
               col_2 %in% c("Floating_AvgPercCovered", "WaterLettuce_AvgPercCovered", "WaterHyacinth_AvgPercCovered")))
# only floating plants

lw_chl_dat %>%
  select(ends_with("TreatFreq"), 
         ends_with("AvgPercCovered"),
         PrevValue,
         Lag) %>%
  group_by(Lag) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & corr >= 0.4 &
           !(col_1 %in% c("Floating_AvgPercCovered", "WaterLettuce_AvgPercCovered", "WaterHyacinth_AvgPercCovered") &
               col_2 %in% c("Floating_AvgPercCovered", "WaterLettuce_AvgPercCovered", "WaterHyacinth_AvgPercCovered")))
# para grass with water hyacinth
# para grass with its own treatment (0.4)

# prev and change in richness
chl_dat %>%
  filter(Lag == 1) %>%
  ggplot(aes(x = PrevValue, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm")

lw_chl_dat %>%
  filter(Lag == 1) %>%
  ggplot(aes(x = PrevValue, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm")

# invasive plant and change in richness
ggplot(chl_dat, aes(x = Hydrilla_AvgPercCovered, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(chl_dat, aes(x = Floating_AvgPercCovered, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(chl_dat, aes(x = ParaGrass_AvgPercCovered, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(chl_dat, aes(x = Torpedograss_AvgPercCovered, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(chl_dat, aes(x = CubanBulrush_AvgPercCovered, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(lw_chl_dat, aes(x = Hydrilla_AvgPercCovered, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(lw_chl_dat, aes(x = Floating_AvgPercCovered, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(lw_chl_dat, aes(x = ParaGrass_AvgPercCovered, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(lw_chl_dat, aes(x = Torpedograss_AvgPercCovered, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(lw_chl_dat, aes(x = CubanBulrush_AvgPercCovered, y = LogRatioQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)


#### fit LW + water atlas models ####

# subset data
chl_dat_mod1 <- filter(chl_dat, Lag == 1) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))
chl_dat_mod2 <- filter(chl_dat, Lag == 2) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))
chl_dat_mod3 <- filter(chl_dat, Lag == 3) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))
chl_dat_mod4 <- filter(chl_dat, Lag == 4) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))
chl_dat_mod5 <- filter(chl_dat, Lag == 5) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))
chl_dat_mod6 <- filter(chl_dat, Lag == 6) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))  

# subset data for Cuban bulrush
chl_dat_rec1 <- filter(chl_dat, Lag == 1 & GSYear > 2013) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))

# fit models for focal invasive plants
chl_foc_mod1 <- feols(LogRatioQual ~ PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod1)
chl_foc_mod2 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod2)
chl_foc_mod3 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod3)
chl_foc_mod4 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod4)
chl_foc_mod5 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod5)
chl_foc_mod6 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod6)  

# add non-focal invasive plants to model
chl_all_mod1 <- feols(LogRatioQual ~ PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod1)
chl_all_mod2 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod2)
chl_all_mod3 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod3)
chl_all_mod4 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod4)
chl_all_mod5 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod5)
chl_all_mod6 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = chl_dat_mod6) 

# add Cuban bulrush
chl_rec_mod1 <- feols(LogRatioQual ~ PrevChl_s + CubanBulrush_TreatFreq + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + CubanBulrush_AvgPercCovered + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = chl_dat_rec1)


#### LW + water model coefficient figures ####

# combine models
chl_foc_mods <- list(chl_foc_mod1, chl_foc_mod2, chl_foc_mod3, chl_foc_mod4, chl_foc_mod5, chl_foc_mod6)
chl_all_mods <- list(chl_all_mod1, chl_all_mod2, chl_all_mod3, chl_all_mod4, chl_all_mod5, chl_all_mod6)

# name models
names(chl_foc_mods) <- names(chl_all_mods) <- c("1", "2", "3", "4", "5", "6")

# rename coefficients
mgmt_coef_names <- c("CubanBulrush_TreatFreq" = "Cuban bulrush",
                     "Torpedograss_TreatFreq" = "Torpedograss",
                     "ParaGrass_TreatFreq" = "Para grass",
                     "Hydrilla_TreatFreq" = "Hydrilla",
                     "Floating_TreatFreq" = "Floating plant")

inv_coef_names <- c("CubanBulrush_AvgPercCovered" = "Cuban bulrush",
                    "Torpedograss_AvgPercCovered" = "Torpedograss",
                    "ParaGrass_AvgPercCovered" = "Para grass",
                    "Hydrilla_AvgPercCovered" = "Hydrilla",
                    "Floating_AvgPercCovered" = "Floating plant")

# figure
foc_mgmt_fig <- modelplot(chl_foc_mods,
                          coef_map = mgmt_coef_names,
                          background = list(geom_vline(xintercept = 0, color = "black",
                                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Management\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Treatment frequency",
       title = "Change in chlorophyll") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

foc_inv_fig <- modelplot(chl_foc_mods,
                         coef_map = inv_coef_names,
                         background = list(geom_vline(xintercept = 0, color = "black",
                                                      size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Invasion\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered",
       title = "Change in chlorophyll") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

all_mgmt_fig <- modelplot(chl_all_mods,
                          coef_map = mgmt_coef_names,
                          background = list(geom_vline(xintercept = 0, color = "black",
                                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Management\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Treatment frequency",
       title = "Change in chlorophyll") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

all_inv_fig <- modelplot(chl_all_mods,
                         coef_map = inv_coef_names,
                         background = list(geom_vline(xintercept = 0, color = "black",
                                                      size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Invasion\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered",
       title = "Change in chlorophyll") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

rec_mgmt_fig <- modelplot(chl_rec_mod1,
                          coef_map = mgmt_coef_names,
                          background = list(geom_vline(xintercept = 0, color = "black",
                                                       size = 0.5, linetype = "dashed"))) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Treatment frequency",
       title = "Change in chlorophyll") +
  def_theme_paper

rec_inv_fig <- modelplot(chl_rec_mod1,
                         coef_map = inv_coef_names,
                         background = list(geom_vline(xintercept = 0, color = "black",
                                                      size = 0.5, linetype = "dashed"))) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered",
       title = "Change in chlorophyll") +
  def_theme_paper


#### export LW + water figures and models ####

# figures
ggsave("output/LW_water_atlas_chlorophyll_focal_invasive_model.eps", foc_inv_fig,
       device = "eps", width = 3.5, height = 3, units = "in")
ggsave("output/LW_water_atlas_chlorophyll_focal_control_model.eps", foc_mgmt_fig,
       device = "eps", width = 3.5, height = 3, units = "in")
ggsave("output/LW_water_atlas_chlorophyll_all_invasive_model.eps", all_inv_fig,
       device = "eps", width = 3.5, height = 5, units = "in")
ggsave("output/LW_water_atlas_chlorophyll_all_control_model.eps", all_mgmt_fig,
       device = "eps", width = 3.5, height = 5, units = "in")
ggsave("output/LW_water_atlas_chlorophyll_recent_invasive_model.eps", rec_inv_fig,
       device = "eps", width = 3.5, height = 4, units = "in")
ggsave("output/LW_water_atlas_chlorophyll_recent_control_model.eps", rec_mgmt_fig,
       device = "eps", width = 3.5, height = 4, units = "in")

# model objects
save(chl_foc_mods, file = "output/LW_water_atlas_chlorophyll_focal_invasive_models.rda")
save(chl_all_mods, file = "output/LW_water_atlas_chlorophyll_all_invasive_models.rda")
save(chl_rec_mod1, file = "output/LW_water_atlas_chlorophyll_recent_invasive_model.rda")


#### fit LW models ####

# subset data
lw_chl_dat_mod1 <- filter(lw_chl_dat, Lag == 1) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))
lw_chl_dat_mod2 <- filter(lw_chl_dat, Lag == 2) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))
lw_chl_dat_mod3 <- filter(lw_chl_dat, Lag == 3) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))
lw_chl_dat_mod4 <- filter(lw_chl_dat, Lag == 4) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))
lw_chl_dat_mod5 <- filter(lw_chl_dat, Lag == 5) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))
lw_chl_dat_mod6 <- filter(lw_chl_dat, Lag == 6) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))  

# subset data for Cuban bulrush
lw_chl_dat_rec1 <- filter(lw_chl_dat, Lag == 1 & GSYear > 2013) %>%
  mutate(PrevChl_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))

# fit models for focal invasive plants
lw_chl_foc_mod1 <- feols(LogRatioQual ~ PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod1)
lw_chl_foc_mod2 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod2)
lw_chl_foc_mod3 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod3)
lw_chl_foc_mod4 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod4)
lw_chl_foc_mod5 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod5)
lw_chl_foc_mod6 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod6)  

# add non-focal invasive plants to model
lw_chl_all_mod1 <- feols(LogRatioQual ~ PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod1)
lw_chl_all_mod2 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod2)
lw_chl_all_mod3 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod3)
lw_chl_all_mod4 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod4)
lw_chl_all_mod5 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod5)
lw_chl_all_mod6 <- feols(LogRatioQual ~  PrevChl_s + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_mod6) 

# add Cuban bulrush
lw_chl_rec_mod1 <- feols(LogRatioQual ~ PrevChl_s + CubanBulrush_TreatFreq + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + CubanBulrush_AvgPercCovered + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = lw_chl_dat_rec1)


#### LW coefficient figures ####

# combine models
lw_chl_foc_mods <- list(lw_chl_foc_mod1, lw_chl_foc_mod2, lw_chl_foc_mod3, lw_chl_foc_mod4, lw_chl_foc_mod5, lw_chl_foc_mod6)
lw_chl_all_mods <- list(lw_chl_all_mod1, lw_chl_all_mod2, lw_chl_all_mod3, lw_chl_all_mod4, lw_chl_all_mod5, lw_chl_all_mod6)

# name models
names(lw_chl_foc_mods) <- names(lw_chl_all_mods) <- c("1", "2", "3", "4", "5", "6")

# rename coefficients
mgmt_coef_names <- c("CubanBulrush_TreatFreq" = "Cuban bulrush",
                     "Torpedograss_TreatFreq" = "Torpedograss",
                     "ParaGrass_TreatFreq" = "Para grass",
                     "Hydrilla_TreatFreq" = "Hydrilla",
                     "Floating_TreatFreq" = "Floating plant")

inv_coef_names <- c("CubanBulrush_AvgPercCovered" = "Cuban bulrush",
                    "Torpedograss_AvgPercCovered" = "Torpedograss",
                    "ParaGrass_AvgPercCovered" = "Para grass",
                    "Hydrilla_AvgPercCovered" = "Hydrilla",
                    "Floating_AvgPercCovered" = "Floating plant")

# figure
lw_foc_mgmt_fig <- modelplot(lw_chl_foc_mods,
                          coef_map = mgmt_coef_names,
                          background = list(geom_vline(xintercept = 0, color = "black",
                                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Management\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Treatment frequency",
       title = "Change in chlorophyll") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

lw_foc_inv_fig <- modelplot(lw_chl_foc_mods,
                         coef_map = inv_coef_names,
                         background = list(geom_vline(xintercept = 0, color = "black",
                                                      size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Invasion\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered",
       title = "Change in chlorophyll") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

lw_all_mgmt_fig <- modelplot(lw_chl_all_mods,
                          coef_map = mgmt_coef_names,
                          background = list(geom_vline(xintercept = 0, color = "black",
                                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Management\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Treatment frequency",
       title = "Change in chlorophyll") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

lw_all_inv_fig <- modelplot(lw_chl_all_mods,
                         coef_map = inv_coef_names,
                         background = list(geom_vline(xintercept = 0, color = "black",
                                                      size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Invasion\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered",
       title = "Change in chlorophyll") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

lw_rec_mgmt_fig <- modelplot(lw_chl_rec_mod1,
                          coef_map = mgmt_coef_names,
                          background = list(geom_vline(xintercept = 0, color = "black",
                                                       size = 0.5, linetype = "dashed"))) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Treatment frequency",
       title = "Change in chlorophyll") +
  def_theme_paper

lw_rec_inv_fig <- modelplot(lw_chl_rec_mod1,
                         coef_map = inv_coef_names,
                         background = list(geom_vline(xintercept = 0, color = "black",
                                                      size = 0.5, linetype = "dashed"))) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered",
       title = "Change in chlorophyll") +
  def_theme_paper


#### export LW figures and models ####

# figures
ggsave("output/LW_chlorophyll_focal_invasive_model.eps", lw_foc_inv_fig,
       device = "eps", width = 3.5, height = 3, units = "in")
ggsave("output/LW_chlorophyll_focal_control_model.eps", lw_foc_mgmt_fig,
       device = "eps", width = 3.5, height = 3, units = "in")
ggsave("output/LW_chlorophyll_all_invasive_model.eps", lw_all_inv_fig,
       device = "eps", width = 3.5, height = 5, units = "in")
ggsave("output/LW_chlorophyll_all_control_model.eps", lw_all_mgmt_fig,
       device = "eps", width = 3.5, height = 5, units = "in")
ggsave("output/LW_chlorophyll_recent_invasive_model.eps", lw_rec_inv_fig,
       device = "eps", width = 3.5, height = 4, units = "in")
ggsave("output/LW_chlorophyll_recent_control_model.eps", lw_rec_mgmt_fig,
       device = "eps", width = 3.5, height = 4, units = "in")

# model objects
save(lw_chl_foc_mods, file = "output/LW_chlorophyll_focal_invasive_models.rda")
save(lw_chl_all_mods, file = "output/LW_chlorophyll_all_invasive_models.rda")
save(lw_chl_rec_mod1, file = "output/LW_chlorophyll_recent_invasive_model.rda")
