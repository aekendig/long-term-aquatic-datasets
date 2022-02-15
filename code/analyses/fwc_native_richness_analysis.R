#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(plotly)
library(tidyverse)
library(GGally)
library(fixest) # FE models
library(modelsummary)
library(inspectdf) # for inspect_cor

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")
nat_plant <- read_csv("intermediate-data/FWC_common_native_plants_formatted.csv")


#### richness-area relationship ####

# richness per waterbody for all years
nat_rich <- nat_plant %>%
  filter(Detected == 1) %>%
  group_by(PermanentID, AreaName, Area_ha) %>%
  summarize(Richness = n_distinct(TaxonName)) %>%
  ungroup() %>%
  mutate(LogRich = log(Richness),
         LogArea = log(Area_ha))

# are richness and area linearly related?
ggplot(nat_rich, aes(x = Area_ha, y = Richness)) +
  geom_point()

ggplot(nat_rich, aes(x = LogArea, y = LogRich)) +
  geom_point()

nat_rich %>%
  filter(Area_ha < 90000) %>%
  ggplot(aes(x = Area_ha, y = Richness)) +
  geom_point()
# no, richness saturates with area

# fit species richness-area relationship
rich_area_mod <- lm(LogRich ~ LogArea, data = nat_rich)
summary(rich_area_mod)
c <- exp(as.numeric(coef(rich_area_mod)[1]))
z <- as.numeric(coef(rich_area_mod)[2])

# simulate relationship
nat_rich_sim <- tibble(Area_ha = seq(min(nat_rich$Area_ha), 
                                     max(nat_rich$Area_ha), 
                                     length.out = 100)) %>%
  mutate(Richness = c*Area_ha^z)

ggplot(nat_rich, aes(x = Area_ha, y = Richness)) +
  geom_point() +
  geom_line(data = nat_rich_sim) +
  coord_cartesian(xlim = c(0, 20000))
# suggests that richness isn't saturated
# not a ton of data to evaluate it


#### edit native plant data ####

# average richness per waterbody
wat_rich <- nat_plant %>%
  group_by(PermanentID, GSYear, Area_ha) %>%
  summarize(Richness = sum(Detected)) %>% # richness by year
  ungroup() %>%
  group_by(PermanentID) %>%
  summarize(AvgRich = mean(Richness, na.rm = T)) %>% # average waterbody richness
  ungroup()

# select data with previous year's detection data
# summarize richness by waterbody and year
nat_plant2 <-  nat_plant %>%
  filter(!is.na(PrevDetected)) %>%
  mutate(Loss = ifelse(Detected == 0 & PrevDetected == 1, 1, 0),
         Gain = ifelse(Detected == 1 & PrevDetected == 0, 1, 0)) %>%
  group_by(PermanentID, GSYear, Area_ha) %>%
  summarize(Richness = sum(Detected),
            PrevRich = sum(PrevDetected),
            Losses = sum(Loss),
            Gains = sum(Gain)) %>%
  ungroup() %>%
  mutate(RatioRich = Richness/PrevRich,
         LogRatioRich = log(RatioRich),
         LogArea = log(Area_ha),
         LogRich = log(Richness)) %>%
  left_join(wat_rich)

# initial visualizations
plot_ly(nat_plant2, x = ~GSYear, y = ~Richness, color = ~PermanentID) %>%
  add_lines() %>% 
  layout(showlegend = FALSE)

plot_ly(nat_plant2, x = ~GSYear, y = ~Gains, color = ~PermanentID) %>%
  add_lines() %>% 
  layout(showlegend = FALSE)

plot_ly(nat_plant2, x = ~GSYear, y = ~Losses, color = ~PermanentID) %>%
  add_lines() %>% 
  layout(showlegend = FALSE)

# one huge outlier for gains/losses in one year
# nat_plant2 %>%
#   filter(Losses > 40 | Gains > 40) %>%
#   inner_join(nat_plant %>%
#                filter(Detected == 1))
# removed this survey (Tohopekaliga, Lake in 2017) from data processing scripts
# survey only had 6 taxa detected and average was around 40
# 5/6 were the firsts in alphabetical order -- seemed incomplete


#### edit other data ####

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

# combine native, invasive, control
# remove missing data
# center and scale variables
nat_dat <- nat_plant2 %>%
  inner_join(inv_plant3) %>%
  inner_join(inv_ctrl3)


#### initial visualizations ####

# covariate correlations
nat_dat %>%
  select(ends_with("TreatFreq"), 
         ends_with("AvgPercCovered"),
         PrevRich,
         Lag) %>%
  group_by(Lag) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & corr >= 0.4)
# only floating plants

# prev and change in richness
nat_dat %>%
  filter(Lag == 1) %>%
  ggplot(aes(x = PrevRich, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm")

# invasive plant and change in richness
ggplot(nat_dat, aes(x = Hydrilla_AvgPercCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(nat_dat, aes(x = Floating_AvgPercCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(nat_dat, aes(x = AlligatorWeed_AvgPercCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(nat_dat, aes(x = CubanBulrush_AvgPercCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(nat_dat, aes(x = ParaGrass_AvgPercCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(nat_dat, aes(x = Torpedograss_AvgPercCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(nat_dat, aes(x = WaterFern_AvgPercCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)


#### fit models ####

# subset data
nat_dat_mod1 <- filter(nat_dat, Lag == 1) %>%
  mutate(PrevRich_s = (PrevRich - mean(PrevRich)) / sd(PrevRich))
nat_dat_mod2 <- filter(nat_dat, Lag == 2) %>%
  mutate(PrevRich_s = (PrevRich - mean(PrevRich)) / sd(PrevRich))
nat_dat_mod3 <- filter(nat_dat, Lag == 3) %>%
  mutate(PrevRich_s = (PrevRich - mean(PrevRich)) / sd(PrevRich))
nat_dat_mod4 <- filter(nat_dat, Lag == 4) %>%
  mutate(PrevRich_s = (PrevRich - mean(PrevRich)) / sd(PrevRich))
nat_dat_mod5 <- filter(nat_dat, Lag == 5) %>%
  mutate(PrevRich_s = (PrevRich - mean(PrevRich)) / sd(PrevRich))
nat_dat_mod6 <- filter(nat_dat, Lag == 6) %>%
  mutate(PrevRich_s = (PrevRich - mean(PrevRich)) / sd(PrevRich))  
  
# fit models for focal invasive plants
nat_foc_mod1 <- feols(LogRatioRich ~ PrevRich_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod1)
nat_foc_mod2 <- feols(LogRatioRich ~  PrevRich_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod2)
nat_foc_mod3 <- feols(LogRatioRich ~  PrevRich_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod3)
nat_foc_mod4 <- feols(LogRatioRich ~  PrevRich_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod4)
nat_foc_mod5 <- feols(LogRatioRich ~  PrevRich_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod5)
nat_foc_mod6 <- feols(LogRatioRich ~  PrevRich_s + Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod6)  
  
# fit models for all invasive plants
nat_all_mod1 <- feols(LogRatioRich ~ PrevRich_s + AlligatorWeed_TreatFreq + CubanBulrush_TreatFreq + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + WaterFern_TreatFreq + AlligatorWeed_AvgPercCovered + CubanBulrush_AvgPercCovered + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered + WaterFern_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod1)
nat_all_mod2 <- feols(LogRatioRich ~  PrevRich_s + AlligatorWeed_TreatFreq + CubanBulrush_TreatFreq + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + WaterFern_TreatFreq + AlligatorWeed_AvgPercCovered + CubanBulrush_AvgPercCovered + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered + WaterFern_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod2)
nat_all_mod3 <- feols(LogRatioRich ~  PrevRich_s + AlligatorWeed_TreatFreq + CubanBulrush_TreatFreq + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + WaterFern_TreatFreq + AlligatorWeed_AvgPercCovered + CubanBulrush_AvgPercCovered + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered + WaterFern_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod3)
nat_all_mod4 <- feols(LogRatioRich ~  PrevRich_s + AlligatorWeed_TreatFreq + CubanBulrush_TreatFreq + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + WaterFern_TreatFreq + AlligatorWeed_AvgPercCovered + CubanBulrush_AvgPercCovered + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered + WaterFern_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod4)
nat_all_mod5 <- feols(LogRatioRich ~  PrevRich_s + AlligatorWeed_TreatFreq + CubanBulrush_TreatFreq + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + WaterFern_TreatFreq + AlligatorWeed_AvgPercCovered + CubanBulrush_AvgPercCovered + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered + WaterFern_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod5)
nat_all_mod6 <- feols(LogRatioRich ~  PrevRich_s + AlligatorWeed_TreatFreq + CubanBulrush_TreatFreq + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + WaterFern_TreatFreq + AlligatorWeed_AvgPercCovered + CubanBulrush_AvgPercCovered + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered + WaterFern_AvgPercCovered | PermanentID + GSYear, data = nat_dat_mod6)  


#### model coefficient figures ####

# combine models
nat_foc_mods <- list(nat_foc_mod1, nat_foc_mod2, nat_foc_mod3, nat_foc_mod4, nat_foc_mod5, nat_foc_mod6)
nat_all_mods <- list(nat_all_mod1, nat_all_mod2, nat_all_mod3, nat_all_mod4, nat_all_mod5, nat_all_mod6)

# name models
names(nat_foc_mods) <- names(nat_all_mods) <- c("1", "2", "3", "4", "5", "6")

# rename coefficients
mgmt_coef_names <- c("WaterFern_TreatFreq" = "Water fern",
                     "Torpedograss_TreatFreq" = "Torpedograss",
                     "ParaGrass_TreatFreq" = "Para grass",
                     "Hydrilla_TreatFreq" = "Hydrilla",
                     "Floating_TreatFreq" = "Floating plant",
                     "CubanBulrush_TreatFreq" = "Cuban bulrush",
                     "AlligatorWeed_TreatFreq" = "Alligator weed")

inv_coef_names <- c("WaterFern_AvgPercCovered" = "Water fern",
                    "Torpedograss_AvgPercCovered" = "Torpedograss",
                    "ParaGrass_AvgPercCovered" = "Para grass",
                    "Hydrilla_AvgPercCovered" = "Hydrilla",
                    "Floating_AvgPercCovered" = "Floating plant",
                    "CubanBulrush_AvgPercCovered" = "Cuban bulrush",
                    "AlligatorWeed_AvgPercCovered" = "Alligator weed")

# figure
foc_mgmt_fig <- modelplot(nat_foc_mods,
                          coef_map = mgmt_coef_names,
                          background = list(geom_vline(xintercept = 0, color = "black",
                                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Management\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Treatment frequency",
       title = "Change in native richness") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

foc_inv_fig <- modelplot(nat_foc_mods,
                         coef_map = inv_coef_names,
                         background = list(geom_vline(xintercept = 0, color = "black",
                                                      size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Invasion\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered",
       title = "Change in native richness") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

all_mgmt_fig <- modelplot(nat_all_mods,
                          coef_map = mgmt_coef_names,
                          background = list(geom_vline(xintercept = 0, color = "black",
                                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Management\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Treatment frequency",
       title = "Change in native richness") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

all_inv_fig <- modelplot(nat_all_mods,
                         coef_map = inv_coef_names,
                         background = list(geom_vline(xintercept = 0, color = "black",
                                                      size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Invasion\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered",
       title = "Change in native richness") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))


#### export figures and models ####

# figures
ggsave("output/fwc_native_richness_focal_invasive_model.eps", foc_inv_fig,
       device = "eps", width = 3.5, height = 3, units = "in")
ggsave("output/fwc_native_richness_focal_control_model.eps", foc_mgmt_fig,
       device = "eps", width = 3.5, height = 3, units = "in")
ggsave("output/fwc_native_richness_all_invasive_model.eps", all_inv_fig,
       device = "eps", width = 3.5, height = 5, units = "in")
ggsave("output/fwc_native_richness_all_control_model.eps", all_mgmt_fig,
       device = "eps", width = 3.5, height = 5, units = "in")

# model objects
save(nat_foc_mods, file = "output/fwc_native_richness_focal_invasive_models.rda")
save(nat_all_mods, file = "output/fwc_native_richness_all_invasive_models.rda")
