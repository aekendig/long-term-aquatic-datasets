#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(inspectdf) # inspect_cor
library(modelsummary) # modelplot
library(patchwork) # combining figures
library(plm) # panel data models
library(pglm) # panel data models
library(sandwich) # vcovHC
library(lmtest) # coeftest
library(glmmTMB) # random effects
library(pals) # color palettes               

# figure settings
source("code/settings/figure_settings.R")

# functions
source("code/generic-functions/model_structure_comparison.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_analysis_formatted.csv") # plant and control data, continuous data
lw_chl <- read_csv("intermediate-data/LW_chlorophyll_formatted.csv")
lwwa_chl <- read_csv("intermediate-data/LW_water_atlas_chlorophyll_formatted.csv")


#### edit data ####

# add max years to invasion data
inv_plant2 <- inv_plant %>%
  filter(!is.na(Lag3AvgPropCovered) & !is.na(Lag3Treated)) %>%
  group_by(CommonName) %>%
  mutate(maxYears = n_distinct(GSYear)) %>%
  ungroup()

# check
inv_plant2 %>%
  select(CommonName, maxYears) %>%
  unique()


#### start here ####
# this method omits torpedograss
# don't just want max years within chlorophyll data because they may not be the same years

# combine chlorophyll, invasive, control
# select waterbodies sampled throughout
chl_dat <- lwwa_chl %>%
  mutate(ValueDiff = QualityValue - PrevValue) %>% # change over time
  filter(!is.na(PrevValue)) %>%
  inner_join(inv_plant2) %>%
  group_by(CommonName, PermanentID, Quarter) %>%
  mutate(nYears = n_distinct(GSYear)) %>% # years per waterbody
  ungroup() %>%
  filter(nYears == maxYears)

# sample sizes
chl_dat %>%
  group_by(CommonName, Quarter) %>%
  summarize(Years = n_distinct(GSYear),
            Waterbodies = n_distinct(PermanentID))
  
  # %>% # remove incomplete time-series
  # group_by(CommonName, Quarter) %>%
  # mutate(nWBs = n_distinct(PermanentID)) %>% # waterbodies per quarter
  # ungroup() %>%
  # group_by(CommonName) %>%
  # mutate(maxWBs = max(nWBs)) %>% # max waterbodies per taxon
  # ungroup() %>%
  # filter(nWBs == maxWBs) # select quarter with max

# taxa
inv_taxa <- sort(unique(chl_dat$CommonName))

# loop through taxa
pdf("output/chlorophyll_continuous_time_series_by_taxon.pdf")

for(i in 1:length(inv_taxa)){
  
  # subset data
  subdat <- chl_dat %>% filter(CommonName == inv_taxa[i])
  subdat_ctrl <- subdat %>% filter(Lag1Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = QualityValue, color = PermanentID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          facet_wrap(~ Quarter) + 
          labs(x = "Year", y = "Chlorophyll a (ug/L)", title = inv_taxa[i]) +
          def_theme_paper +
          theme(plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()


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
