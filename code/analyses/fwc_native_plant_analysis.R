#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(inspectdf) # inspect_cor
library(broom) # glance
library(pals) # color palettes

# figure settings
source("code/settings/figure_settings.R")

# functions
source("code/generic-functions/proportion_transformations.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_analysis_formatted.csv") # plant and control data, continuous data
nat_plant <- read_csv("intermediate-data/FWC_common_native_plants_formatted.csv")


#### edit data ####

# select relevant invasive plant columns
# join with native plant data
# summarize
comb_dat <- inv_plant %>%
  select(PermanentID, GSYear, CommonName, PropCovered, Lag1Treated) %>%
  inner_join(nat_plant) %>%
  group_by(CommonName, PermanentID, TaxonName, Habitat) %>%
  summarize(AvgPAC = mean(PropCovered * 100),
            TreatFreq = mean(Lag1Treated),
            YearsDetected = sum(Detected),
            YearsSurveyed = n_distinct(GSYear)) %>%
  ungroup() %>%
  group_by(CommonName, TaxonName) %>%
  mutate(AvgPAC_c = AvgPAC - mean(AvgPAC)) %>%
  ungroup() %>%
  mutate(YearsUndetected = YearsSurveyed - YearsDetected)

# check that invasive plant values are repeated by permanentID
comb_dat %>%
  group_by(CommonName, PermanentID) %>%
  summarize(nPAC = n_distinct(AvgPAC),
            nTreat = n_distinct(TreatFreq)) %>%
  ungroup() %>%
  filter(nPAC > 1 | nTreat > 1)
# yes

# check that all lakes were surveyed the same number of times
comb_dat %>%
  group_by(CommonName) %>%
  summarize(nYears = n_distinct(YearsSurveyed)) %>%
  ungroup() %>%
  filter(nYears > 1)
# yes
  

#### initial visualizations ####

# covariate correlations
comb_dat %>%
  select(CommonName, PermanentID, TreatFreq, AvgPAC_c) %>%
  unique() %>%
  select(-PermanentID) %>%
  group_by(CommonName) %>%
  inspect_cor() %>% 
  ungroup()
# treatment not significantly correlated with PAC 
# except for hydrilla, but corr = 0.2

ggplot(comb_dat, aes(AvgPAC, YearsDetected, color = TaxonName)) +
  geom_point(size = 0.75, alpha = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")

ggplot(comb_dat, aes(TreatFreq, YearsDetected, color = TaxonName)) +
  geom_point(size = 0.75, alpha = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")


#### fit models ####

# species with too many 1's or 0's
comb_dat %>%
  group_by(TaxonName, CommonName) %>%
  summarize(YearsDetected = sum(YearsDetected),
            YearsSurveyed = sum(YearsSurveyed)) %>%
  ungroup() %>%
  filter(YearsDetected >= (YearsSurveyed - 20) | YearsDetected <= 20)
# para grass models might not work, some have a few or a lot of detections

# apply model to each invasive plant and taxon
plant_mods <- comb_dat %>%
  select(CommonName, TaxonName, Habitat, YearsDetected, YearsUndetected, 
         AvgPAC_c, TreatFreq) %>%
  nest(data = c(YearsDetected, YearsUndetected, AvgPAC_c, TreatFreq)) %>%
  mutate(fit = map(data, ~glm(cbind(YearsDetected, YearsUndetected) ~ AvgPAC_c + TreatFreq, 
                              data = ., family = binomial)))
# model fit error

# convert warnings to errors to break loop
options(warn = 2)

# common/taxa list
# update as pairs break loop
comm_tax <- comb_dat %>%
  select(CommonName, TaxonName) %>%
  unique() %>%
  mutate(SppPair = paste(CommonName, TaxonName, sep = ", ")) %>%
  filter(!(SppPair %in% c("Para grass, Cicuta maculata",
                          "Para grass, Crinum americanum",
                          "Para grass, Salix spp.",
                          "Water hyacinth, Mayaca fluviatilis",
                          "Water hyacinth, Myriophyllum laxum/pinnatum",
                          "Water lettuce, Fuirena spp.",
                          "Water lettuce, Mayaca fluviatilis",
                          "Water lettuce, Micranthemum umbrosum",
                          "Water lettuce, Myriophyllum laxum/pinnatum",
                          "Water lettuce, Nitella spp.",
                          "Water lettuce, Nymphoides aquatica")))

# find issue model
for(i in 1:nrow(comm_tax)) {
    
    print(comm_tax$SppPair[i])
    
    print(glm(cbind(YearsDetected, YearsUndetected) ~ AvgPAC_c + TreatFreq,
              data = filter(comb_dat, 
                            CommonName == comm_tax$CommonName[i] & 
                              TaxonName == comm_tax$TaxonName[i]), 
              family = binomial))
  
}

# remove problematic species pairs
comb_dat2 <- comb_dat %>%
  inner_join(comm_tax)

# refit models
plant_mods2 <- comb_dat2 %>%
  select(CommonName, TaxonName, Habitat, YearsDetected, YearsUndetected, 
         AvgPAC_c, TreatFreq) %>%
  nest(data = c(YearsDetected, YearsUndetected, AvgPAC_c, TreatFreq)) %>%
  mutate(fit = map(data, ~glm(cbind(YearsDetected, YearsUndetected) ~ AvgPAC_c + TreatFreq, 
                              data = ., family = binomial)))


#### coefficients ####

# all coefficients
plant_coef <- plant_mods2 %>%
  mutate(tidied = map(fit, tidy)) %>%
  select(CommonName, TaxonName, Habitat, tidied) %>%
  unnest(tidied)

# models with significant PAC effects
plant_coef_PAC <- plant_coef %>%
  mutate(Sig = if_else(term == "AvgPAC_c" & p.value < 0.05, "yes", "no")) %>%
  filter(term == "AvgPAC_c") %>%
  rename(Coef = estimate) %>%
  select(CommonName, TaxonName, Habitat, Coef, Sig) %>%
  full_join(plant_coef %>%
              filter(term == "(Intercept)") %>%
              select(CommonName, TaxonName, estimate) %>%
              rename(Intercept = estimate))

# models with significant treatment effects
plant_coef_treat <- plant_coef %>%
  mutate(Sig = if_else(term == "TreatFreq" & p.value < 0.05, "yes", "no")) %>%
  filter(term == "TreatFreq") %>%
  rename(Coef = estimate) %>%
  select(CommonName, TaxonName, Habitat, Coef, Sig) %>%
  full_join(plant_coef %>%
              filter(term == "(Intercept)") %>%
              select(CommonName, TaxonName, estimate) %>%
              rename(Intercept = estimate))


#### predicted values ####

# add coefficients to full data
comb_dat_PAC <- comb_dat2 %>%
  left_join(plant_coef_PAC) %>%
  mutate(Pred = logit2prob(Intercept + Coef * AvgPAC_c) * YearsSurveyed,
         Sig = fct_relevel(Sig, "yes"),
         PanelName = case_when(CommonName == "Hydrilla" ~ "(A) hydrilla",
                               CommonName == "Water hyacinth" ~ "(B) water hyacinth",
                               CommonName == "Water lettuce" ~ "(C) water lettuce",
                               CommonName == "Cuban bulrush" ~ "(A) Cuban bulrush",
                               CommonName == "Para grass" ~ "(B) para grass",
                               CommonName == "Torpedograss" ~ "(C) torpedograss")) 

comb_dat_treat <- comb_dat2 %>%
  left_join(plant_coef_treat) %>%
  mutate(Pred = logit2prob(Intercept + Coef * TreatFreq) * YearsSurveyed,
         Sig = fct_relevel(Sig, "yes"),
         Treat = TreatFreq * 3,
         PanelName = case_when(CommonName == "Hydrilla" ~ "(A) hydrilla management",
                               CommonName == "Water hyacinth" ~ "(B) water hyacinth management",
                               CommonName == "Water lettuce" ~ "(C) water lettuce management",
                               CommonName == "Cuban bulrush" ~ "(A) Cuban bulrush management",
                               CommonName == "Para grass" ~ "(B) para grass management",
                               CommonName == "Torpedograss" ~ "(C) torpedograss management"))

# split by focal/non-focal
foc_dat_PAC <- comb_dat_PAC %>%
  filter(CommonName %in% c("Hydrilla", "Water hyacinth", "Water lettuce"))

foc_dat_treat <- comb_dat_treat %>%
  filter(CommonName %in% c("Hydrilla", "Water hyacinth", "Water lettuce"))

non_foc_dat_PAC <- comb_dat_PAC %>%
  filter(CommonName %in% c("Cuban bulrush", "Para grass", "Torpedograss"))

non_foc_dat_treat <- comb_dat_treat %>%
  filter(CommonName %in% c("Cuban bulrush", "Para grass", "Torpedograss"))



#### figures ####

foc_fig_PAC <- ggplot(foc_dat_PAC, aes(x = AvgPAC, y = Pred, group = TaxonName, color = Habitat)) +
  geom_line(aes(alpha = Sig, size = Sig)) +
  facet_wrap(~ PanelName, scales = "free") +
  scale_alpha_manual(values = c(1, 0.5), name = "P < 0.05") +
  scale_size_manual(values = c(0.75, 0.5), name = "P < 0.05") +
  labs(x = "Invasive plant PAC", y = "Years detected") +
  scale_color_manual(values = brewer.set2(n = 3)) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0)) +
  guides(color = guide_legend(override.aes = list(size = 0.75, alpha = 1)))

foc_fig_treat <- ggplot(foc_dat_treat, aes(x = Treat, y = Pred, group = TaxonName, color = Habitat)) +
  geom_line(aes(alpha = Sig, size = Sig)) +
  facet_wrap(~ PanelName, scales = "free") +
  scale_alpha_manual(values = c(1, 0.5), name = "P < 0.05") +
  scale_size_manual(values = c(0.75, 0.5), name = "P < 0.05") +
  labs(x = "Years managed (out of 3)", y = "Years detected") +
  scale_color_manual(values = brewer.set2(n = 3)) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0)) +
  guides(color = guide_legend(override.aes = list(size = 0.75, alpha = 1)))

non_foc_fig_PAC <- ggplot(non_foc_dat_PAC, aes(x = AvgPAC, y = Pred, group = TaxonName, color = Habitat)) +
  geom_line(aes(alpha = Sig, size = Sig)) +
  facet_wrap(~ PanelName, scales = "free") +
  scale_alpha_manual(values = c(1, 0.5), name = "P < 0.05") +
  scale_size_manual(values = c(0.75, 0.5), name = "P < 0.05") +
  labs(x = "Invasive plant PAC", y = "Years detected") +
  scale_color_manual(values = brewer.set2(n = 3)) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0)) +
  guides(color = guide_legend(override.aes = list(size = 0.75, alpha = 1)))

non_foc_fig_treat <- ggplot(non_foc_dat_treat, aes(x = Treat, y = Pred, group = TaxonName, color = Habitat)) +
  geom_line(aes(alpha = Sig, size = Sig)) +
  facet_wrap(~ PanelName, scales = "free") +
  scale_alpha_manual(values = c(1, 0.5), name = "P < 0.05") +
  scale_size_manual(values = c(0.75, 0.5), name = "P < 0.05") +
  labs(x = "Years managed (out of 3)", y = "Years detected") +
  scale_color_manual(values = brewer.set2(n = 3)) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0)) +
  guides(color = guide_legend(override.aes = list(size = 0.75, alpha = 1)))



#### older code ####

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
