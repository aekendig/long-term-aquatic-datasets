#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)
library(GGally)
library(broom)
library(janitor)

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
plant_com <- read_csv("intermediate-data/FWC_common_plant_community_formatted.csv")


#### edit data ####

# add SurveyYear to invasive plant
inv_plant$SurveyYear = year(inv_plant$SurveyDate)

# check that invasive plant and plant community data match
inv_plant %>%
  filter(!is.na(SpeciesAcres) & !(SurveyYear %in% c(2000, 2001))) %>%
  select(PermanentID, SurveyDate) %>%
  unique() %>%
  anti_join(plant_com %>%
              select(PermanentID, SurveyDate) %>%
              unique())
# water hyacinth was surveyd 9 days before other plants (or type-o?)

plant_com %>%
  select(PermanentID, SurveyDate) %>%
  unique() %>%
  anti_join(inv_plant %>%
              select(PermanentID, SurveyDate))

# check that management data doesn't go later than plant surveys
inv_ctrl %>%
  left_join(plant_com %>%
              group_by(PermanentID) %>%
              summarize(MaxYearPlant = max(GSYear)) %>%
              ungroup()) %>%
  filter(GSYear > MaxYearPlant) %>%
  select(PermanentID, GSYear, MaxYearPlant) %>%
  unique()
# it does for 400 cases

# summarize data
plant_com2 <- plant_com %>%
  filter(PreCtrl == "post ctrl data") %>%
  group_by(PermanentID, TaxonName, Origin, Habitat) %>%
  summarize(YearsDetected = sum(Detected),
            YearsSurveyed = n()) %>%
  ungroup() %>%
  mutate(YearsUndetected = YearsSurveyed - YearsDetected,
         PropDetected = YearsDetected / YearsSurveyed,
         Origin = fct_recode(Origin, "non-native" = "Exotic",
                             "native" = "Native") %>%
           fct_relevel("native"))

inv_plant2 <- inv_plant %>%
  filter(!is.na(SpeciesAcres) & GSYear >= min(inv_ctrl$GSYear)) %>%
  mutate(CommonName = fct_recode(CommonName,
                                 "WaterHyacinth" = "Water hyacinth",
                                 "WaterLettuce" = "Water lettuce")) %>%
  select(PermanentID, GSYear, CommonName, PropCovered) %>%
  pivot_wider(names_from = CommonName,
              values_from = PropCovered) %>%
  mutate(Floating = WaterHyacinth + WaterLettuce) %>%
  group_by(PermanentID) %>%
  summarize(Hydrilla = mean(Hydrilla) * 100,
            WaterHyacinth = mean(WaterHyacinth) * 100,
            WaterLettuce = mean(WaterLettuce) * 100,
            Floating = mean(Floating) * 100) %>%
  ungroup()

inv_ctrl2 <- inv_ctrl %>%
  left_join(plant_com %>%
              group_by(PermanentID) %>%
              summarize(MaxYearPlant = max(GSYear)) %>%
              ungroup()) %>%
  filter(GSYear <= MaxYearPlant) %>%
  group_by(PermanentID, TaxonName) %>%
  summarize(TreatmentFreq = mean(Lag0Treated),
            YearsTreatment = n()) %>%
  ungroup() %>%
  mutate(TaxonName = fct_recode(TaxonName, 
                                "FloatingTrt" = "Eichhornia crassipes",
                                "HydrillaTrt" = "Hydrilla verticillata",
                                "FloatingTrt2" = "Pistia stratiotes")) %>%
  filter(TaxonName != "FloatingTrt2") %>%
  pivot_wider(names_from = TaxonName,
              values_from = TreatmentFreq)

# combine data
plant_com3 <- plant_com2 %>%
  inner_join(inv_plant2) %>%
  inner_join(inv_ctrl2)


#### initial visualizations ####

ggplot(plant_com3, aes(Hydrilla, PropDetected, color = Origin)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

ggplot(plant_com3, aes(WaterHyacinth, PropDetected, color = Origin)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

ggplot(plant_com3, aes(WaterLettuce, PropDetected, color = Origin)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

ggplot(plant_com3, aes(Floating, PropDetected, color = Origin)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

ggplot(plant_com3, aes(FloatingTrt, PropDetected, color = Origin)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

ggplot(plant_com3, aes(HydrillaTrt, PropDetected)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

plant_com3 %>%
  select(Hydrilla, HydrillaTrt, WaterHyacinth, WaterLettuce, Floating, FloatingTrt) %>%
  unique() %>%
  ggpairs()
# all are less than 0.4 except water lettuce and water hyacinth (0.7) - use floating

# update data, remove high floating point
# prevented models from converging (see below)
plant_com3b <- plant_com3 %>%
  filter(Floating < 4)

plant_com3b %>%
  select(Hydrilla, HydrillaTrt, WaterHyacinth, WaterLettuce, Floating, FloatingTrt) %>%
  unique() %>%
  ggpairs()
# water hyacinth and lettuce are no less correlated


#### models ####

# species with too many 1's or 0's
plant_com3 %>%
  group_by(TaxonName) %>%
  summarize(YearsDetected = sum(YearsDetected),
            YearsSurveyed = sum(YearsSurveyed)) %>%
  ungroup() %>%
  filter(YearsDetected >= (YearsSurveyed - 20) | YearsDetected <= 20)
# nothing obvious

# apply model to each taxon
plant_mods <- plant_com3 %>%
  select(TaxonName, YearsDetected, YearsUndetected, 
         Hydrilla, Floating, HydrillaTrt, FloatingTrt) %>%
  nest(data = c(YearsDetected, YearsUndetected, Hydrilla, Floating, HydrillaTrt, 
                FloatingTrt)) %>%
  mutate(fit = map(data, ~glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + Floating + FloatingTrt, 
                              data = ., family = binomial)))
# multiple models with errors

# model summaries
plant_mods %>%
  mutate(glanced = map(fit, glance)) %>%
  select(TaxonName, glanced) %>%
  unnest(glanced) %>%
  data.frame()
# Nuphar advena has a much lower loglik

# examine species
nu_ad_dat <- plant_com3 %>%
  filter(TaxonName == "Nuphar advena")

nu_ad_mod <- glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + Floating + FloatingTrt, 
                 data = nu_ad_dat, family = binomial)
summary(nu_ad_mod)
# very large estimates

# remove high floating point
nu_ad_dat2 <- plant_com3b %>%
  filter(TaxonName == "Nuphar advena")

# refit model
nu_ad_mod2 <- glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + Floating + FloatingTrt, 
                  data = nu_ad_dat2, family = binomial)
summary(nu_ad_mod2)
# much more reasonable estimates without warnings

# update models
plant_mods2 <- plant_com3b %>%
  select(TaxonName, YearsDetected, YearsUndetected, 
         Hydrilla, WaterHyacinth, WaterLettuce, Floating,
         HydrillaTrt, FloatingTrt) %>%
  nest(data = c(YearsDetected, YearsUndetected, 
                Hydrilla, WaterHyacinth, WaterLettuce, Floating,
                HydrillaTrt, FloatingTrt)) %>%
  mutate(fit = map(data, ~glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + WaterHyacinth + WaterLettuce + FloatingTrt, 
                              data = ., family = binomial)),
         fit2 = map(data, ~glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + Floating + FloatingTrt, 
                              data = ., family = binomial)))
# no errors

# model summaries
plant_aic <- plant_mods2 %>%
  mutate(glanced = map(fit, glance)) %>%
  select(TaxonName, glanced) %>%
  unnest(glanced) %>%
  select(TaxonName, AIC) %>%
  rename(AICSep = AIC) %>%
  full_join(plant_mods2 %>%
              mutate(glanced = map(fit2, glance)) %>%
              select(TaxonName, glanced) %>%
              unnest(glanced) %>%
              select(TaxonName, AIC) %>%
              rename(AICAdd = AIC)) %>%
  mutate(AICDiff = AICSep - AICAdd)

ggplot(plant_aic, aes(x = AICDiff)) +
  geom_vline(xintercept = 4, color = "blue") +
  geom_vline(xintercept = -4, color = "blue") +
  geom_histogram(binwidth = 1)
# when models differ, AICSep is a better fit

# coefficients
plant_coef <- plant_mods2 %>%
  mutate(tidied = map(fit, tidy)) %>%
  select(TaxonName, tidied) %>%
  unnest(tidied) %>%
  mutate(sig = case_when(p.value < 0.1 & p.value >= 0.05 ~ "< 0.1",
                         p.value < 0.05 & p.value >= 0.01 ~ "< 0.05",
                         p.value < 0.01 & p.value >= 0.001 ~ "< 0.01",
                         p.value < 0.001 ~ "< 0.001",
                         TRUE ~ "\u2265 0.1"),
         term = fct_recode(term, "intercept" = "(Intercept)",
                           "floating\ntreat." = "FloatingTrt",
                           "hydrilla" = "Hydrilla",
                           "hydrilla\ntreat." = "HydrillaTrt",
                           "water\nhyacinth" = "WaterHyacinth",
                           "water\nlettuce" = "WaterLettuce") %>%
           fct_relevel("intercept", "hydrilla", "water\nhyacinth",
                       "water\nlettuce", "hydrilla\ntreat.", "floating\ntreat.")) %>%
  left_join(plant_com3 %>%
              select(TaxonName, Origin, Habitat) %>%
              unique()) %>%
  mutate(HabitatLabel = fct_recode(Habitat,
                                   "(B) emersed taxa" = "Emersed",
                                   "(C) floating taxa" = "Floating",
                                   "(D) submersed taxa" = "Submersed"))

# divide coefficients by group
plant_coef2 <- plant_coef %>%
  mutate(HabitatLabel = "(A) all taxa") %>%
  full_join(plant_coef)


#### figure ####
cairo_pdf("output/fwc_plant_community_coefficients.pdf", width = 6.5, height = 6.5)
ggplot(plant_coef2, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = NA, outlier.shape = 4, size = 0.25) +
  geom_point(aes(fill = sig, shape = Origin), 
             position = position_jitter(width = 0.1), size = 1, stroke = 0.3) +
  scale_fill_viridis_d(name = expression(paste(italic(P), " value", sep = ""))) +
  scale_shape_manual(values = c(21, 24)) +
  labs(y = "Estimate (log-odds)") +
  facet_wrap(~ HabitatLabel, scales = "free_y") +
  def_theme_paper +
  theme(strip.text = element_text(size = 11, color="black", hjust = 0),
        axis.title.x = element_blank(),
        legend.position = c(0.12, 0.93),
        legend.box = "horizontal",
        legend.spacing.x = unit(-0.1, "cm"),
        legend.box.background = element_rect(fill = NA, color = NA),
        legend.key.height = unit(2, "mm")) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))
dev.off()


#### table ####

plant_table <- plant_coef %>%
  relocate(TaxonName, Habitat, Origin) %>%
  mutate(term = str_replace(term, "\n", " ") %>%
           fct_relevel("intercept", "hydrilla", "water hyacinth",
                       "water lettuce", "hydrilla treat.", "floating treat."),
         Habitat = tolower(Habitat),
         across(where(is.double) & !p.value, ~ round_half_up(., digits = 3))) %>%
  rename("Taxon" = "TaxonName",
         "Parameter" = "term",
         "Estimate" = "estimate",
         "Std.error" = "std.error",
         "Z" = "statistic",
         "P" = "p.value",
         "Significance" = "sig") %>%
  select(-c(HabitatLabel)) %>%
  arrange(Parameter, Taxon)

write_csv(plant_table, "output/fwc_plant_community_coefficients.csv")
