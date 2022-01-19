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
nat_plant <- read_csv("intermediate-data/FWC_common_native_plants_formatted.csv")


#### edit data ####

# add SurveyYear to invasive plant
inv_plant$SurveyYear = year(inv_plant$SurveyDate)

# check that invasive plant and plant community data match
inv_plant %>%
  filter(!is.na(SpeciesAcres) & !(SurveyYear %in% c(2000, 2001))) %>%
  select(PermanentID, SurveyDate) %>%
  unique() %>%
  anti_join(nat_plant %>%
              select(PermanentID, SurveyDate) %>%
              unique())
# water hyacinth was surveyed 9 days before other plants (or type-o?)

nat_plant %>%
  select(PermanentID, SurveyDate) %>%
  unique() %>%
  anti_join(inv_plant %>%
              select(PermanentID, SurveyDate))

# check that management data doesn't go later than plant surveys
inv_ctrl %>%
  left_join(nat_plant %>%
              group_by(PermanentID) %>%
              summarize(MaxYearPlant = max(GSYear)) %>%
              ungroup()) %>%
  filter(GSYear > MaxYearPlant) %>%
  select(PermanentID, GSYear, MaxYearPlant) %>%
  unique()
# it does for 400 cases

# summarize data
nat_plant2 <- nat_plant %>%
  filter(PreCtrl == "post ctrl data") %>%
  group_by(PermanentID, TaxonName, Habitat) %>%
  summarize(YearsDetected = sum(Detected),
            YearsSurveyed = n()) %>%
  ungroup() %>%
  mutate(YearsUndetected = YearsSurveyed - YearsDetected,
         PropDetected = YearsDetected / YearsSurveyed)

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
  left_join(nat_plant %>%
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
nat_plant3 <- nat_plant2 %>%
  inner_join(inv_plant2) %>%
  inner_join(inv_ctrl2)


#### initial visualizations ####

ggplot(nat_plant3, aes(Hydrilla, PropDetected, color = Habitat)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

ggplot(nat_plant3, aes(WaterHyacinth, PropDetected, color = Habitat)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

ggplot(nat_plant3, aes(WaterLettuce, PropDetected, color = Habitat)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

ggplot(nat_plant3, aes(Floating, PropDetected, color = Habitat)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

ggplot(nat_plant3, aes(FloatingTrt, PropDetected, color = Habitat)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

ggplot(nat_plant3, aes(HydrillaTrt, PropDetected)) +
  geom_point(size = 0.75, alpha = 0.5) +
  facet_wrap(~ TaxonName)

nat_plant3 %>%
  select(Hydrilla, HydrillaTrt, WaterHyacinth, WaterLettuce, Floating, FloatingTrt) %>%
  unique() %>%
  ggpairs()
# all are less than 0.4 except water lettuce and water hyacinth (0.7) - use floating

# update data, remove high floating point
# prevented models from converging (see below)
nat_plant3b <- nat_plant3 %>%
  filter(Floating < 40)

nat_plant3b %>%
  select(Hydrilla, HydrillaTrt, WaterHyacinth, WaterLettuce, Floating, FloatingTrt) %>%
  unique() %>%
  ggpairs()
# water hyacinth and lettuce less correlated


#### models ####

# species with too many 1's or 0's
nat_plant3 %>%
  group_by(TaxonName) %>%
  summarize(YearsDetected = sum(YearsDetected),
            YearsSurveyed = sum(YearsSurveyed)) %>%
  ungroup() %>%
  filter(YearsDetected >= (YearsSurveyed - 20) | YearsDetected <= 20)
# nothing obvious

# apply model to each taxon
plant_mods <- nat_plant3 %>%
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
nu_ad_dat <- nat_plant3 %>%
  filter(TaxonName == "Nuphar advena")

nu_ad_mod <- glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + Floating + FloatingTrt, 
                 data = nu_ad_dat, family = binomial)
summary(nu_ad_mod)
# very large estimates

# remove high floating point
nu_ad_dat2 <- nat_plant3b %>%
  filter(TaxonName == "Nuphar advena")

# refit model
nu_ad_mod2 <- glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + Floating + FloatingTrt, 
                  data = nu_ad_dat2, family = binomial)
summary(nu_ad_mod2)
# much more reasonable estimates without warnings

# update models
plant_mods2 <- nat_plant3b %>%
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
# two errors

# taxa list
taxa <- sort(unique(nat_plant3b$TaxonName))
taxa2 <- taxa[taxa != "Mayaca fluviatilis"] 
# use to run loop after identifying M. fluviatilis as the first warning

# convert warnings to errors to break loop
options(warn = 2)

# find issue model
for(i in 1:length(taxa)){
  
  print(taxa[i])
  
  print(glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + WaterHyacinth + WaterLettuce + FloatingTrt,
      data = filter(nat_plant3b, TaxonName == taxa[i]), family = binomial))
  
}

# refit model
ma_fl_dat <- nat_plant3b %>%
  filter(TaxonName == "Mayaca fluviatilis")
ma_fl_mod <- glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + WaterHyacinth + FloatingTrt, 
                  data = ma_fl_dat, family = binomial)
summary(ma_fl_mod)
# water lettuce triggers warning

my_la_dat <- nat_plant3b %>%
  filter(TaxonName == "Myriophyllum laxum/pinnatum")
my_la_mod <- glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + WaterHyacinth + FloatingTrt, 
                 data = my_la_dat, family = binomial)
summary(my_la_mod)
# water lettuce triggers warning

# reset warning
options(warn = 1)

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

# update models to remove errors
plant_mods3 <- nat_plant3b %>%
  select(TaxonName, YearsDetected, YearsUndetected, 
         Hydrilla, WaterHyacinth, WaterLettuce, Floating,
         HydrillaTrt, FloatingTrt) %>%
  nest(data = c(YearsDetected, YearsUndetected, 
                Hydrilla, WaterHyacinth, WaterLettuce, Floating,
                HydrillaTrt, FloatingTrt)) %>%
  mutate(fit = map(data, ~glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + WaterHyacinth + WaterLettuce + FloatingTrt, 
                              data = ., family = binomial)),
         fit2 = map(data, ~glm(cbind(YearsDetected, YearsUndetected) ~ Hydrilla + HydrillaTrt + WaterHyacinth + FloatingTrt, 
                               data = ., family = binomial)))

#### coefficients ####
plant_coef1 <- plant_mods3 %>%
  mutate(tidied = map(fit, tidy)) %>%
  select(TaxonName, tidied) %>%
  unnest(tidied)

plant_coef2 <- plant_mods3 %>%
  mutate(tidied = map(fit2, tidy)) %>%
  select(TaxonName, tidied) %>%
  unnest(tidied)

plant_coef <- plant_coef1 %>%
  filter(!(TaxonName %in% c("Mayaca fluviatilis", "Myriophyllum laxum/pinnatum"))) %>%
  full_join(plant_coef2 %>%
              filter(TaxonName %in% c("Mayaca fluviatilis", "Myriophyllum laxum/pinnatum"))) %>%
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
  left_join(nat_plant3 %>%
              select(TaxonName, Habitat) %>%
              unique()) %>%
  mutate(HabitatLabel = fct_recode(Habitat,
                                   "(B) emersed taxa" = "Emersed",
                                   "(C) floating taxa" = "Floating",
                                   "(D) submersed taxa" = "Submersed"))

# divide coefficients by group
plant_coef_fig <- plant_coef %>%
  mutate(HabitatLabel = "(A) all taxa") %>%
  full_join(plant_coef)


#### one-sample t-test ####

# use all coefficients
plant_coef_test <- plant_coef %>%
  mutate(HabitatLabel = "(A) all taxa")

# t-tests
plant_coef_test2 <- plant_coef_test %>%
  select(term, estimate) %>%
  nest(data = estimate) %>%
  mutate(fit = map(data, ~ t.test(.)))

# summary
plant_coef_test3 <- plant_coef_test2 %>%
  mutate(tidied = map(fit, tidy)) %>%
  select(term, tidied) %>%
  unnest(tidied) %>%
  mutate(sig = case_when(p.value < 0.1 & p.value >= 0.05 ~ ".",
                           p.value < 0.05 & p.value >= 0.01 ~ "*",
                           p.value < 0.01 & p.value >= 0.001 ~ "**",
                           p.value < 0.001 ~ "***",
                           TRUE ~ NA_character_),
         HabitatLabel = "(A) all taxa")

# save
write_csv(plant_coef_test3, "output/fwc_native_plant_ttests.csv")


#### figure ####
cairo_pdf("output/fwc_native_plant_coefficients.pdf", width = 6.5, height = 6.5)
ggplot(plant_coef_fig , aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = NA, outlier.shape = 4, size = 0.25) +
  geom_point(aes(fill = sig), shape = 21, 
             position = position_jitter(width = 0.1), size = 1, stroke = 0.3) +
  geom_text(data = plant_coef_test3, aes(y = 5.1, label = sig)) +
  scale_fill_viridis_d(name = expression(paste(italic(P), " value", sep = ""))) +
  labs(y = "Estimate (log-odds)") +
  facet_wrap(~ HabitatLabel, scales = "free_y") +
  def_theme_paper +
  theme(strip.text = element_text(size = 11, color="black", hjust = 0),
        axis.title.x = element_blank(),
        legend.position = c(0.05, 0.4),
        legend.box = "horizontal",
        legend.spacing.x = unit(-0.1, "cm"),
        legend.box.background = element_rect(fill = NA, color = NA),
        legend.key.height = unit(2, "mm"))
dev.off()


#### table ####

plant_table <- plant_coef %>%
  relocate(TaxonName, Habitat) %>%
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

write_csv(plant_table, "output/fwc_native_plant_coefficients.csv")
