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
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv",
                      col_types = list(PrevPropCovered = col_double(),
                                       PrevAreaCoveredRaw_ha = col_double(),
                                       SurveyDays = col_double(),
                                       RatioCovered = col_double(),
                                       LogRatioCovered = col_double(),
                                       LogitPrevPropCovered = col_double(),
                                       LogRatioCovered = col_double(),
                                       LogitPrevPropCovered = col_double(),
                                       LogRatioCovered = col_double()))
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")
nat_plant <- read_csv("intermediate-data/FWC_common_native_plants_formatted.csv")
lw_qual <- read_csv("intermediate-data/LW_quality_formatted.csv")
wa_qual <- read_csv("intermediate-data/water_atlas_quality_formatted.csv")


#### edit native data ####

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

# select data with previous year's detection data
# summarize richness by waterbody and year
nat_plant2 <- nat_plant %>%
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
         LogRich = log(Richness))

# initial visualizations
plot_ly(nat_plant2, x = ~GSYear, y = ~Richness, color = ~PermanentID) %>%
  add_lines() %>% 
  layout(showlegend = FALSE)

# several cases of zero richness
nat_plant2 %>%
  filter(GSYear %in% c(1985, 1987, 1989, 1991, 1993) & Richness == 0) %>%
  group_by(PermanentID) %>%
  count() %>%
  ungroup()
  
# richness-area with year variation
ggplot(nat_plant2, aes(x = Area_ha, y = Richness)) +
  geom_point()

ggplot(nat_plant2, aes(x = LogArea, y = LogRich)) +
  geom_point()

ggplot(nat_plant2, aes(x = LogArea, y = Richness)) +
  geom_point()
# need to account for area in the model

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


#### edit other data ####

# make wide by inv. plant species
inv_plant2 <- inv_plant %>%
  filter(!is.na(PrevPropCovered)) %>% # remove lakes/years without surveys
  mutate(CommonName = fct_recode(CommonName,
                                 "WaterHyacinth" = "Water hyacinth",
                                 "WaterLettuce" = "Water lettuce",
                                 "AlligatorWeed" = "Alligator weed",
                                 "WildTaro" = "Wild taro",
                                 "CubanBulrush" = "Cuban bulrush",
                                 "WaterFern" = "Water fern",
                                 "ParaGrass" = "Para grass")) %>%
  select(PermanentID, GSYear, CommonName, LogRatioCovered, PrevPropCovered, MinSurveyorExperience) %>%
  pivot_wider(names_from = CommonName,
              values_from = c(PrevPropCovered, LogRatioCovered),
              names_glue = "{CommonName}_{.value}") %>%
  mutate(across(.cols = ends_with("PrevPropCovered"), ~ .x * 100),
         Floating_PrevPropCovered = WaterHyacinth_PrevPropCovered + WaterLettuce_PrevPropCovered)

# make wide by control target
inv_ctrl2 <- inv_ctrl %>%
  select(PermanentID, Species, GSYear, starts_with("Lag")) %>% # remove duplication of floating plant treatment
  select(-contains("All")) %>% # remove to make wide
  unique() %>%
  mutate(Species = fct_recode(Species,
                              "Floating" = "Floating Plants (Eichhornia and Pistia)",
                              "Hydrilla" = "Hydrilla verticillata",
                              "AlligatorWeed" = "Alternanthera philoxeroides",
                              "WildTaro" = "Colocasia esculenta",
                              "CubanBulrush" = "Oxycaryum cubense",
                              "CubanBulrush" = "Cyperus blepharoleptos",
                              "Torpedograss" = "Panicum repens",
                              "WaterFern" = "Salvinia minima",
                              "ParaGrass" = "Urochloa mutica")) %>%
  pivot_wider(names_from = Species,
              values_from = starts_with("Lag"),
              names_glue = "{Species}_{.value}") %>%
  full_join(inv_ctrl %>% # add "all treatment" back in
              select(PermanentID, GSYear, contains("All")) %>%
              unique())

# check for missing data
filter(inv_ctrl2, across(ends_with("Lag0Treated"), ~ is.na(.x)))
# none

# water quality
qual <- lw_qual %>%
  full_join(wa_qual) %>%
  filter(QualityMetric == "Secchi_ft") %>%
  group_by(PermanentID, GSYear) %>%
  summarize(QualityValue = mean(QualityValue),
            MonthsSampled = mean(MonthsSampled))


#### combine data ####

# combine native, invasive, control
# remove missing data
# center and scale variables
nat_dat <- nat_plant2 %>%
  inner_join(inv_plant2) %>% # this should maintain nat_plant2's row number
  inner_join(inv_ctrl2) %>%
  filter(!is.na(MinSurveyorExperience) & !(LogRatioRich %in% c(-Inf, Inf))) %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
         PrevRich_s = (PrevRich - mean(PrevRich)) / sd(PrevRich),
         across(ends_with("PrevPropCovered"), ~ .x - mean(.x))) # center percent changes

# missing surveyor data (checked before removing)
sum(is.na(nat_dat$MinSurveyorExperience)) # 5 rows

# zero richness (checked before removing)
filter(nat_dat, LogRatioRich %in% c(-Inf, Inf)) # 1 row

# add quality data
# re-center and scale variables (smaller dataset)
nat_qual_dat <- nat_dat %>%
  inner_join(qual %>%
               select(PermanentID, GSYear, QualityValue, MonthsSampled)) %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
         PrevRich_s = (PrevRich - mean(PrevRich)) / sd(PrevRich),
         across(ends_with("PrevPropCovered"), ~ .x - mean(.x)),
         Clarity_s = (QualityValue - mean(QualityValue)) / sd(QualityValue))


#### initial visualizations ####

# covariate correlations
nat_dat %>%
  select(ends_with("Lag0Treated"), "Lag0AllTreated") %>%
  inspect_cor() %>% 
  filter(p_value < 0.05 & corr >= 0.4)
# all and floating: 0.6
# all and hydrilla: 0.4

nat_dat %>%
  select(ends_with("Lag0Treated"), 
         ends_with("LogRatioCovered"),
         PrevRich, SurveyorExperience_s) %>%
  inspect_cor() %>% 
  filter(p_value < 0.05 & corr >= 0.4)
# water fern and alligator weed cover (0.8)

nat_dat %>%
  select(ends_with("Lag0Treated"),
         ends_with("PrevPercCovered"),
         PrevRich, SurveyorExperience_s) %>%
  inspect_cor() %>% 
  filter(p_value < 0.05 & corr >= 0.4)
# none

nat_qual_dat %>%
  select(ends_with("Lag0Treated"), 
         ends_with("LogRatioCovered"),
         PrevRich, SurveyorExperience_s, Clarity_s) %>%
  inspect_cor() %>% 
  filter(p_value < 0.05 & corr >= 0.4)
# water fern and alligator weed cover (0.8)

nat_qual_dat %>%
  select(ends_with("Lag0Treated"),
         ends_with("PrevPercCovered"),
         PrevRich, SurveyorExperience_s, Clarity_s) %>%
  inspect_cor() %>% 
  filter(p_value < 0.05 & corr >= 0.4)
# none

# prev and change in richness
ggplot(nat_dat, aes(x = PrevRich, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm")

# hydrilla and change in richness
ggplot(nat_dat, aes(x = Hydrilla_LogRatioCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(nat_dat, aes(x = Hydrilla_PrevPropCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm")

# water hyacinth and change in richness
ggplot(nat_dat, aes(x = WaterHyacinth_LogRatioCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(nat_dat, aes(x = WaterHyacinth_PrevPropCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm")

# water lettuce and change in richness
ggplot(nat_dat, aes(x = WaterLettuce_LogRatioCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm")
# all 3 log ratio are positive
# conditions promoting invasive plant growth promote native richness

ggplot(nat_dat, aes(x = WaterLettuce_PrevPropCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm")

# floating and change in richness
ggplot(nat_dat, aes(x = Floating_PrevPropCovered, y = LogRatioRich)) +
  geom_point() +
  geom_smooth(method = "lm")


#### fit models ####

# log(richness_t/richness_t-1) ~ richness_t-1 + treatment + log(inv_t/inv_t-1)
# does richness increase or decrease?
# initial richness constrains change: can lose more if you have more, can gain more if you have fewer
# control performed between t-1 and t may reduce richness
# control performed before t-1 may increase richness by reducing the invasive plant pop
# an increase in invasive cover may exclude plants/make them rarer and a decrease in cover may allow them to recover
# initial invasive cover may lead to losses when higher

# models with just hydrilla, water hyacinth, and water lettuce

# rich_mod <- feols(LogRatioRich ~ PrevRich_s + Floating_Lag0Treated + Hydrilla_Lag0Treated + Hydrilla_LogRatioCovered + WaterHyacinth_LogRatioCovered + WaterLettuce_LogRatioCovered + SurveyorExperience_s | PermanentID + GSYear, data = nat_dat)
# summary(rich_mod)
# 
# rich_modb <- feols(LogRatioRich ~ PrevRich_s + Floating_Lag0Treated + Hydrilla_Lag0Treated + Hydrilla_PrevPercCovered + Floating_PrevPercCovered + SurveyorExperience_s | PermanentID + GSYear, data = nat_dat)
# summary(rich_modb)
# 
# rich_qual_mod <- feols(LogRatioRich ~ PrevRich_s + Floating_Lag0Treated + Hydrilla_Lag0Treated + Hydrilla_LogRatioCovered + WaterHyacinth_LogRatioCovered + WaterLettuce_LogRatioCovered + SurveyorExperience_s + Clarity_s | PermanentID + GSYear, data = nat_qual_dat)
# summary(rich_qual_mod)
# 
# rich_qual_modb <- feols(LogRatioRich ~ PrevRich_s + Floating_Lag0Treated + Hydrilla_Lag0Treated + Hydrilla_PrevPercCovered + Floating_PrevPercCovered + SurveyorExperience_s + Clarity_s | PermanentID + GSYear, data = nat_qual_dat)
# summary(rich_qual_modb)
# 
# # hydrilla by itself
# rich_hydr_mod <- feols(LogRatioRich ~ PrevRich_s + Hydrilla_Lag0Treated + Hydrilla_LogRatioCovered + SurveyorExperience_s | PermanentID + GSYear, data = nat_dat)
# summary(rich_hydr_mod)
# 
# rich_hydr_modb <- feols(LogRatioRich ~ PrevRich_s + Hydrilla_Lag0Treated + Hydrilla_PrevPercCovered + SurveyorExperience_s | PermanentID + GSYear, data = nat_dat)
# summary(rich_hydr_modb)
# # estimates are all very similar to full model
# 
# # floating by itself
# rich_float_mod <- feols(LogRatioRich ~ PrevRich_s + Floating_Lag0Treated + WaterHyacinth_LogRatioCovered + WaterLettuce_LogRatioCovered + SurveyorExperience_s | PermanentID + GSYear, data = nat_dat)
# summary(rich_float_mod)
# 
# rich_float_modb <- feols(LogRatioRich ~ PrevRich_s + Floating_Lag0Treated + Floating_PrevPercCovered + SurveyorExperience_s | PermanentID + GSYear, data = nat_dat)
# summary(rich_float_modb)
# # estimates are all very similar to full model

# full dataset
rich_all_mod <- feols(LogRatioRich ~ PrevRich_s + AlligatorWeed_PrevPropCovered + CubanBulrush_PrevPropCovered + Floating_PrevPropCovered + Hydrilla_PrevPropCovered + ParaGrass_PrevPropCovered + Torpedograss_PrevPropCovered + WaterFern_PrevPropCovered + WildTaro_PrevPropCovered + AlligatorWeed_Lag0Treated + CubanBulrush_Lag0Treated + Floating_Lag0Treated + Hydrilla_Lag0Treated + ParaGrass_Lag0Treated + Torpedograss_Lag0Treated + WaterFern_Lag0Treated + WildTaro_Lag0Treated + SurveyorExperience_s | PermanentID + GSYear, data = nat_dat)
summary(rich_all_mod)

# water quality
rich_all_qual_mod <- feols(LogRatioRich ~ PrevRich_s + AlligatorWeed_PrevPropCovered + CubanBulrush_PrevPropCovered + Floating_PrevPropCovered + Hydrilla_PrevPropCovered + ParaGrass_PrevPropCovered + Torpedograss_PrevPropCovered + WaterFern_PrevPropCovered + WildTaro_PrevPropCovered + AlligatorWeed_Lag0Treated + CubanBulrush_Lag0Treated + Floating_Lag0Treated + Hydrilla_Lag0Treated + ParaGrass_Lag0Treated + Torpedograss_Lag0Treated + WaterFern_Lag0Treated + WildTaro_Lag0Treated + SurveyorExperience_s + Clarity_s | PermanentID + GSYear, data = nat_qual_dat)
summary(rich_all_qual_mod)


#### model coefficient figure ####

# combine models
rich_mods <- list(rich_all_mod, rich_all_qual_mod)

# name models
names(rich_mods) <- c("no", "yes")

# rename coefficients
coef_names <- c("Clarity_s" = "Water clarity",
                "SurveyorExperience_s" = "Surveyor experience",
                "WildTaro_PrevPropCovered" = "Wild taro abundance",
                "WaterFern_PrevPropCovered" = "Water fern abundance",
                "Torpedograss_PrevPropCovered" = "Torpedograss abundance",
                "ParaGrass_PrevPropCovered" = "Para grass abundance",
                "CubanBulrush_PrevPropCovered" = "Cuban bulrush abundance",
                "AlligatorWeed_PrevPropCovered" = "Alligator weed abundance",
                "Floating_PrevPropCovered" = "Floating plant abundance",
                "Hydrilla_PrevPropCovered" = "Hydrilla abundance",
                "WildTaro_Lag0Treated" = "Wild taro treatment",
                "WaterFern_Lag0Treated" = "Water fern treatment",
                "Torpedograss_Lag0Treated" = "Torpedograss treatment",
                "ParaGrass_Lag0Treated" = "Para grass treatment",
                "CubanBulrush_Lag0Treated" = "Cuban bulrush treatment",
                "AlligatorWeed_Lag0Treated" = "Alligator weed treatment",
                "Floating_Lag0Treated" = "Floating plant treatment",
                "Hydrilla_Lag0Treated" = "Hydrilla treatment",
                "PrevRich_s" = "Initial richness")

# figure
modelplot(rich_mods,
          coef_map = coef_names,
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(option = "plasma", end = 0.7, direction = -1,
                        name = "Model\nincludes\nclarity") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = ""))) +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10))
# models estimates are similar

pdf("output/fwc_native_richness_treatment_clarity_all_inv_model.pdf", width = 4, height = 6)
modelplot(rich_all_qual_mod,
          coef_map = coef_names,
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed")),
          size = 0.3, fatten = 2, color = "blue") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       title = "Native plant taxonomic richness") +
  def_theme_paper +
  theme(plot.title = element_text(size = 11, hjust = 0.5))
dev.off()
