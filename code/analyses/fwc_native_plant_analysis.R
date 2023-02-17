#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(pals) # color palettes
library(patchwork)
library(Hmsc)
library(corrplot)
library(tidybayes)

# figure settings
source("code/settings/figure_settings.R")

# import data
dat <- read_csv("intermediate-data/FWC_common_native_plants_invaded_data_formatted.csv")
dat_full <- read_csv("intermediate-data/FWC_common_native_plants_invasive_species_data_formatted.csv")


#### edit data ####

# split data by invasive species
hydr_dat <- filter(dat, CommonName == "Hydrilla") %>%
  arrange(PermanentID, GSYear, TaxonName)
flpl_dat <- filter(dat, CommonName == "floating plants") %>%
  arrange(PermanentID, GSYear, TaxonName)
cubu_dat <- filter(dat, CommonName == "Cuban bulrush") %>%
  arrange(PermanentID, GSYear, TaxonName)
pagr_dat <- filter(dat, CommonName == "Para grass") %>%
  arrange(PermanentID, GSYear, TaxonName)
torp_dat <- filter(dat, CommonName == "Torpedograss") %>%
  arrange(PermanentID, GSYear, TaxonName)

# site-by-species matrices
# remove taxa with no presences
hydr_mat <- hydr_dat %>%
  select(PermanentID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(PermanentID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
flpl_mat <- flpl_dat %>%
  select(PermanentID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(PermanentID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
cubu_mat <- cubu_dat %>%
  select(PermanentID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(PermanentID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
pagr_mat <- pagr_dat %>%
  select(PermanentID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(PermanentID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
torp_mat <- torp_dat %>%
  select(PermanentID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(PermanentID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()

# covariates
hydr_cov <- hydr_dat %>%
  distinct(PermanentID, GSYear, PercCovered, RecentTreatment) %>%
  arrange(PermanentID, GSYear) %>%
  select(-c(PermanentID, GSYear)) %>%
  data.frame()
flpl_cov <- flpl_dat %>%
  distinct(PermanentID, GSYear, PercCovered, RecentTreatment) %>%
  arrange(PermanentID, GSYear) %>%
  select(-c(PermanentID, GSYear)) %>%
  data.frame()
cubu_cov <- cubu_dat %>%
  distinct(PermanentID, GSYear, PercCovered, RecentTreatment) %>%
  arrange(PermanentID, GSYear) %>%
  select(-c(PermanentID, GSYear)) %>%
  data.frame()
pagr_cov <- pagr_dat %>%
  distinct(PermanentID, GSYear, PercCovered, RecentTreatment) %>%
  arrange(PermanentID, GSYear) %>%
  select(-c(PermanentID, GSYear)) %>%
  data.frame()
torp_cov <- torp_dat %>%
  distinct(PermanentID, GSYear, PercCovered, RecentTreatment) %>%
  arrange(PermanentID, GSYear) %>%
  select(-c(PermanentID, GSYear)) %>%
  data.frame()

# study design
hydr_stud <- hydr_dat %>%
  distinct(PermanentID, GSYear) %>%
  arrange(PermanentID, GSYear) %>%
  mutate(PermanentID = as.factor(PermanentID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
flpl_stud <- flpl_dat %>%
  distinct(PermanentID, GSYear) %>%
  arrange(PermanentID, GSYear) %>%
  mutate(PermanentID = as.factor(PermanentID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
cubu_stud <- cubu_dat %>%
  distinct(PermanentID, GSYear) %>%
  arrange(PermanentID, GSYear) %>%
  mutate(PermanentID = as.factor(PermanentID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
pagr_stud <- pagr_dat %>%
  distinct(PermanentID, GSYear) %>%
  arrange(PermanentID, GSYear) %>%
  mutate(PermanentID = as.factor(PermanentID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
torp_stud <- torp_dat %>%
  distinct(PermanentID, GSYear) %>%
  arrange(PermanentID, GSYear) %>%
  mutate(PermanentID = as.factor(PermanentID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()

# traits
hydr_hab <- hydr_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(hydr_hab) <- hydr_hab$TaxonName
hydr_trait <- select(hydr_hab, -TaxonName)
flpl_hab <- flpl_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(flpl_hab) <- flpl_hab$TaxonName
flpl_trait <- select(flpl_hab, -TaxonName)
cubu_hab <- cubu_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(cubu_hab) <- cubu_hab$TaxonName
cubu_trait <- select(cubu_hab, -TaxonName)
pagr_hab <- pagr_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(pagr_hab) <- pagr_hab$TaxonName
pagr_trait <- select(pagr_hab, -TaxonName)
torp_hab <- torp_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(torp_hab) <- torp_hab$TaxonName
torp_trait <- select(torp_hab, -TaxonName)


#### define models ####

hydr_mod = Hmsc(Y = hydr_mat,
                XData = hydr_cov, 
                XFormula = ~PercCovered * RecentTreatment,
                TrData = hydr_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = hydr_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = hydr_stud$GSYear), 
                                 PermanentID = HmscRandomLevel(units = hydr_stud$PermanentID)))

flpl_mod = Hmsc(Y = flpl_mat,
                XData = flpl_cov, 
                XFormula = ~PercCovered * RecentTreatment,
                TrData = flpl_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = flpl_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = flpl_stud$GSYear), 
                                 PermanentID = HmscRandomLevel(units = flpl_stud$PermanentID)))

cubu_mod = Hmsc(Y = cubu_mat,
                XData = cubu_cov, 
                XFormula = ~PercCovered * RecentTreatment,
                TrData = cubu_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = cubu_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = cubu_stud$GSYear), 
                                 PermanentID = HmscRandomLevel(units = cubu_stud$PermanentID)))

pagr_mod = Hmsc(Y = pagr_mat,
                XData = pagr_cov, 
                XFormula = ~PercCovered * RecentTreatment,
                TrData = pagr_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = pagr_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = pagr_stud$GSYear), 
                                 PermanentID = HmscRandomLevel(units = pagr_stud$PermanentID)))

torp_mod = Hmsc(Y = torp_mat,
                XData = torp_cov, 
                XFormula = ~PercCovered * RecentTreatment,
                TrData = torp_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = torp_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = torp_stud$GSYear), 
                                 PermanentID = HmscRandomLevel(units = torp_stud$PermanentID)))

#### fit models ####

# MCMC settings
nChains = 3
thin = 10
samples = 500
transient = 500
# compared 3 chains/500 samples/200 transient to
# 2 chains/1000 samples/500 transient and the distribution
# of effective sample sizes and PSRF didn't improve

# fit model
hydr_fit = sampleMcmc(hydr_mod,
                      thin = thin,
                      samples = samples,
                      transient = transient,
                      nChains = nChains,
                      nParallel = nChains)
save(hydr_fit, file = "output/fwc_native_plant_hydrilla_hmsc.rda")

flpl_fit = sampleMcmc(flpl_mod,
                      thin = thin,
                      samples = samples,
                      transient = transient,
                      nChains = nChains,
                      nParallel = nChains)
save(flpl_fit, file = "output/fwc_native_plant_floating_plants_hmsc.rda")

cubu_fit = sampleMcmc(cubu_mod,
                      thin = thin,
                      samples = samples,
                      transient = transient,
                      nChains = nChains,
                      nParallel = nChains)
save(cubu_fit, file = "output/fwc_native_plant_cuban_bulrush_hmsc.rda")

pagr_fit = sampleMcmc(pagr_mod,
                      thin = thin,
                      samples = samples,
                      transient = transient,
                      nChains = nChains,
                      nParallel = nChains)
save(pagr_fit, file = "output/fwc_native_plant_paragrass_hmsc.rda")

torp_fit = sampleMcmc(torp_mod,
                      thin = thin,
                      samples = samples,
                      transient = transient,
                      nChains = nChains,
                      nParallel = nChains)
save(torp_fit, file = "output/fwc_native_plant_torpedograss_hmsc.rda")


#### reload models ####

load("output/fwc_native_plant_hydrilla_hmsc.rda")
load("output/fwc_native_plant_floating_plants_hmsc.rda")
load("output/fwc_native_plant_Cuban_bulrush_hmsc.rda")
load("output/fwc_native_plant_paragrass_hmsc.rda")
load("output/fwc_native_plant_torpedograss_hmsc.rda")


#### evaluate fit ####

# posterior samples
hydr_post = convertToCodaObject(hydr_fit)
flpl_post = convertToCodaObject(flpl_fit)
cubu_post = convertToCodaObject(cubu_fit)
pagr_post = convertToCodaObject(pagr_fit)
torp_post = convertToCodaObject(torp_fit)

# effective sample sizes
# total samples = samples * chains
hist(effectiveSize(hydr_post$Beta))
hist(effectiveSize(flpl_post$Beta))
hist(effectiveSize(cubu_post$Beta))
hist(effectiveSize(pagr_post$Beta)) # almost uniform
hist(effectiveSize(torp_post$Beta))

# Gelman diagnostics (potential scale reduction factors)
# one indicates better convergence among chains
hist(gelman.diag(hydr_post$Beta, multivariate=FALSE)$psrf)
hist(gelman.diag(flpl_post$Beta, multivariate=FALSE)$psrf)
hist(gelman.diag(cubu_post$Beta, multivariate=FALSE)$psrf)
hist(gelman.diag(pagr_post$Beta, multivariate=FALSE)$psrf)
hist(gelman.diag(torp_post$Beta, multivariate=FALSE)$psrf)

# explanatory power, compares posterior distribution to observed values
hydr_preds = computePredictedValues(hydr_fit)
hydr_eval = evaluateModelFit(hM = hydr_fit, predY = hydr_preds)
hist(hydr_eval$RMSE)
hist(hydr_eval$TjurR2)
hist(hydr_eval$AUC)

flpl_preds = computePredictedValues(flpl_fit,
                                    thin = 10,
                                    nParallel = nChains)
# had to use settings because vector memory was exhausted
flpl_eval = evaluateModelFit(hM = flpl_fit, predY = flpl_preds)
hist(flpl_eval$RMSE)
hist(flpl_eval$TjurR2)
hist(flpl_eval$AUC)

cubu_preds = computePredictedValues(cubu_fit)
cubu_eval = evaluateModelFit(hM = cubu_fit, predY = cubu_preds)
hist(cubu_eval$RMSE)
hist(cubu_eval$TjurR2)
hist(cubu_eval$AUC)

pagr_preds = computePredictedValues(pagr_fit)
pagr_eval = evaluateModelFit(hM = pagr_fit, predY = pagr_preds)
hist(pagr_eval$RMSE)
hist(pagr_eval$TjurR2)
hist(pagr_eval$AUC)

torp_preds = computePredictedValues(torp_fit)
torp_eval = evaluateModelFit(hM = torp_fit, predY = torp_preds)
hist(torp_eval$RMSE)
hist(torp_eval$TjurR2)
hist(torp_eval$AUC)


#### coefficient figures ####

# summaries
hydr_sum <- as_tibble(summary(hydr_post$Beta)$statistics) %>%
  mutate(species = rep(colnames(hydr_mat), each = 4),
         covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                         ncol(hydr_mat))) %>%
  full_join(as_tibble(summary(hydr_post$Beta)$quantiles) %>%
              mutate(species = rep(colnames(hydr_mat), each = 4),
                     covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                                     ncol(hydr_mat)))) %>%
  mutate(invader = "hydrilla")

flpl_sum <- as_tibble(summary(flpl_post$Beta)$statistics) %>%
  mutate(species = rep(colnames(flpl_mat), each = 4),
         covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                         ncol(flpl_mat))) %>%
  full_join(as_tibble(summary(flpl_post$Beta)$quantiles) %>%
              mutate(species = rep(colnames(flpl_mat), each = 4),
                     covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                                     ncol(flpl_mat)))) %>%
  mutate(invader = "floating plants")

cubu_sum <- as_tibble(summary(cubu_post$Beta)$statistics) %>%
  mutate(species = rep(colnames(cubu_mat), each = 4),
         covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                         ncol(cubu_mat))) %>%
  full_join(as_tibble(summary(cubu_post$Beta)$quantiles) %>%
              mutate(species = rep(colnames(cubu_mat), each = 4),
                     covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                                     ncol(cubu_mat)))) %>%
  mutate(invader = "Cuban bulrush")

pagr_sum <- as_tibble(summary(pagr_post$Beta)$statistics) %>%
  mutate(species = rep(colnames(pagr_mat), each = 4),
         covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                         ncol(pagr_mat))) %>%
  full_join(as_tibble(summary(pagr_post$Beta)$quantiles) %>%
              mutate(species = rep(colnames(pagr_mat), each = 4),
                     covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                                     ncol(pagr_mat)))) %>%
  mutate(invader = "paragrass")

torp_sum <- as_tibble(summary(torp_post$Beta)$statistics) %>%
  mutate(species = rep(colnames(torp_mat), each = 4),
         covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                         ncol(torp_mat))) %>%
  full_join(as_tibble(summary(torp_post$Beta)$quantiles) %>%
              mutate(species = rep(colnames(torp_mat), each = 4),
                     covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                                     ncol(torp_mat)))) %>%
  mutate(invader = "torpedograss")

# combine
foc_coef_sum <- full_join(hydr_sum, flpl_sum) %>%
  relocate(invader, species, covariate) %>%
  rename(naive_se = "Naive SE",
         time_series_se = "Time-series SE",
         quant_2_5 = "2.5%",
         quant_25 = "25%",
         quant_50 = "50%",
         quant_75 = "75%",
         quant_97_5 = "97.5%") %>%
  left_join(dat %>%
              distinct(TaxonName, Habitat) %>%
              rename(species = TaxonName)) %>%
  mutate(covariate = fct_relevel(covariate,
                                 "intercept", "invasive plant PAC", "management"),
         invader = fct_relevel(invader, "hydrilla"),
         habitat_num = as.numeric(as.factor(Habitat)),
         sp_num = habitat_num * 100 + as.numeric(as.factor(species)),
         species = fct_reorder(species, sp_num, .desc = T),
         sig = case_when(quant_2_5 > 0 & quant_97_5 > 0 ~ "sig",
                         quant_2_5 < 0 & quant_97_5 < 0 ~ "sig",
                         TRUE ~ "not sig"))

non_foc_coef_sum <- full_join(cubu_sum, pagr_sum) %>%
  full_join(torp_sum) %>%
  relocate(invader, species, covariate) %>%
  rename(naive_se = "Naive SE",
         time_series_se = "Time-series SE",
         quant_2_5 = "2.5%",
         quant_25 = "25%",
         quant_50 = "50%",
         quant_75 = "75%",
         quant_97_5 = "97.5%") %>%
  left_join(dat %>%
              distinct(TaxonName, Habitat) %>%
              rename(species = TaxonName)) %>%
  mutate(covariate = fct_relevel(covariate,
                                 "intercept", "invasive plant PAC", "management"),
         cov2 = if_else(covariate == "invasive plant PAC",
                        paste(invader, "PAC"),
                        as.character(covariate)) %>%
           fct_relevel("intercept", "Cuban bulrush PAC", "paragrass PAC",
                       "torpedograss PAC", "management"),
         habitat_num = as.numeric(as.factor(Habitat)),
         sp_num = habitat_num * 100 + as.numeric(as.factor(species)),
         species = fct_reorder(species, sp_num, .desc = T),
         sig = case_when(quant_2_5 > 0 & quant_97_5 > 0 ~ "sig",
                         quant_2_5 < 0 & quant_97_5 < 0 ~ "sig",
                         TRUE ~ "not sig"))

# coefficient plots
foc_coef_fig <- foc_coef_sum %>%
  filter(covariate != "intercept") %>%
  ggplot(aes(x = Mean, y = species, color = Habitat, alpha = sig)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = quant_2_5,
                     xmax = quant_97_5),
                 height = 0) +
  geom_point() +
  facet_grid(invader ~ covariate, scales = "free") +
  scale_color_manual(values = brewer.set2(n = 3)) +
  scale_alpha_manual(values = c(0.2, 1), guide = "none") +
  labs(x = "Model estimate ± 95% CI") +
  def_theme_paper +
  theme(legend.position = c(0.8, 0.95),
        axis.text.y = element_text(size = 6.5, face = "italic"),
        axis.title.y = element_blank())

hydr_coef_fig <- foc_coef_sum %>%
  filter(covariate != "intercept" & invader == "hydrilla") %>%
  ggplot(aes(x = Mean, y = species, color = Habitat, alpha = sig)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = quant_2_5,
                     xmax = quant_97_5),
                 height = 0) +
  geom_point() +
  facet_grid(~ covariate, scales = "free") +
  scale_color_manual(values = brewer.set2(n = 3)) +
  scale_alpha_manual(values = c(0.2, 1), guide = "none") +
  labs(x = "Hydrilla model estimate ± 95% CI") +
  def_theme_paper +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.y = element_text(size = 6.5, face = "italic"),
        axis.title.y = element_blank())

flpl_coef_fig <- foc_coef_sum %>%
  filter(covariate != "intercept" & invader == "floating plants") %>%
  ggplot(aes(x = Mean, y = species, color = Habitat, alpha = sig)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = quant_2_5,
                     xmax = quant_97_5),
                 height = 0) +
  geom_point() +
  facet_grid(~ covariate, scales = "free") +
  scale_color_manual(values = brewer.set2(n = 3)) +
  scale_alpha_manual(values = c(0.2, 1), guide = "none") +
  labs(x = "Floating plant model estimate ± 95% CI") +
  def_theme_paper +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.y = element_text(size = 6.5, face = "italic"),
        axis.title.y = element_blank())

cubu_coef_fig <- non_foc_coef_sum %>%
  filter(covariate != "intercept" & invader == "Cuban bulrush") %>%
  ggplot(aes(x = Mean, y = species, color = Habitat, alpha = sig)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = quant_2_5,
                     xmax = quant_97_5),
                 height = 0) +
  geom_point() +
  facet_grid(~ cov2, scales = "free") +
  scale_color_manual(values = brewer.set2(n = 3)) +
  scale_alpha_manual(values = c(0.2, 1), guide = "none") +
  labs(x = "Cuban bulrush model estimate ± 95% CI") +
  def_theme_paper +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.y = element_text(size = 6.5, face = "italic"),
        axis.title.y = element_blank())

pagr_coef_fig <- non_foc_coef_sum %>%
  filter(covariate != "intercept" & invader == "paragrass") %>%
  ggplot(aes(x = Mean, y = species, color = Habitat, alpha = sig)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = quant_2_5,
                     xmax = quant_97_5),
                 height = 0) +
  geom_point() +
  facet_grid(~ cov2, scales = "free") +
  scale_color_manual(values = brewer.set2(n = 3)) +
  scale_alpha_manual(values = c(0.2, 1), guide = "none") +
  labs(x = "Paragrass model estimate ± 95% CI") +
  def_theme_paper +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.y = element_text(size = 6.5, face = "italic"),
        axis.title.y = element_blank())

torp_coef_fig <- non_foc_coef_sum %>%
  filter(covariate != "intercept" & invader == "torpedograss") %>%
  ggplot(aes(x = Mean, y = species, color = Habitat, alpha = sig)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = quant_2_5,
                     xmax = quant_97_5),
                 height = 0) +
  geom_point() +
  facet_grid(~ cov2, scales = "free") +
  scale_color_manual(values = brewer.set2(n = 3)) +
  scale_alpha_manual(values = c(0.2, 1), guide = "none") +
  labs(x = "Torpedograss model estimate ± 95% CI") +
  def_theme_paper +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.y = element_text(size = 6.5, face = "italic"),
        axis.title.y = element_blank())

ggsave("output/fwc_native_plant_coefficients_focal_invaders.png", foc_coef_fig,
       device = "png", width = 6.5, height = 10, units = "in")
ggsave("output/fwc_native_plant_coefficients_hydrilla.png", hydr_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_floating_plants.png", flpl_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_cuban_bulrush.png", cubu_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_paragrass.png", pagr_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_torpedograss.png", torp_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")

write_csv(foc_coef_sum, "output/fwc_native_plant_coefficients_focal_invaders.csv")
write_csv(non_foc_coef_sum, "output/fwc_native_plant_coefficients_non_focal_invaders.csv")


#### species richness figures ####

# code from: https://github.com/hmsc-r/HMSC/issues/48

# Percentiles used in calculation
p <- c(.025,.5,.975)
p_names <- paste0(p*100)
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>%
  set_names(nm = p_names)

# ggplot of gradient
plot_grad_fun <- function(mod, raw_dat, plot_title){
  
  # variables
  v_expV <- mod$covNames %>% .[-1] # remove intercept
  df_CIs <-NULL
  df_raw <-NULL
  schange <- NULL
  npoints <- 20
  
  # loop over all exp vars
  for (i in 1:2){
    
    # gradients
    grad <- constructGradient(mod, focalVariable = v_expV[i],
                              ngrid = npoints,
                              non.focalVariables = 2) # (2 = net effect)
    
    # use gradient to make predictions
    pred <- predict(mod, Gradient = grad, expected = T) # (probit: expected = T)
    
    # dataframe with 1500 pred * npoints units of exp vars
    predGrad <- pred %>%
      lapply(rowSums) %>% # richness of all pred values (see line 53 in plotGradient function for CWM)
      rlist::list.rbind() %>% as_tibble() # merge list elements to rows
    
    # dataframe with quantiles as cols
    CIs <- predGrad %>%
      pivot_longer(everything(), names_to = "grad", values_to = "pred") %>%
      mutate(across(everything(), type.convert, as.is = T)) %>%
      group_by(grad) %>%
      # mean_hdi() %>% # tried this too, similar results
      summarise(across(everything(), p_funs)) %>%
      select(-1) %>%
      rename_with(~c("Qua_low","Qua_mid","Qua_high"))
    
    # change in richness
    schange <- predGrad %>%
      rename(unit0 = "1",
             unit1 = as.character(npoints)) %>%
      mutate(change = unit1 - unit0) %>%
      median_qi(change) %>%
      add_column(covNames = v_expV[i], .before = 1) %>%
      bind_rows(schange)
    
    # dataframe with all exp vars
    df_CIs <- CIs %>% 
      add_column(gradient = grad$XDataNew[, v_expV[i]], .before = 1) %>%
      add_column(covNames = v_expV[i], .before = 1) %>%
      bind_rows(df_CIs)
    
    # summarize raw data
    rich_dat <- raw_dat %>%
      group_by(across(all_of(c("PermanentID", "GSYear", v_expV[i])))) %>%
      summarize(Qua_mid = sum(Detected)) %>%
      ungroup() %>%
      rename(gradient = !!v_expV[i])
    
    # combine raw data 
    df_raw <- rich_dat %>%
      add_column(covNames = v_expV[i], .before = 1) %>%
      bind_rows(df_raw)
    
  }
  
  # change variable names
  df_CIs2 <- df_CIs %>%
    mutate(covNames = fct_recode(covNames,
                                 "invasive plant PAC" = "PercCovered",
                                 "management" = "RecentTreatment"))
  
  df_raw2 <- df_raw %>%
    mutate(covNames = fct_recode(covNames,
                                 "invasive plant PAC" = "PercCovered",
                                 "management" = "RecentTreatment"))
  
  fig_out <- ggplot(df_raw2, aes(x = gradient, y = Qua_mid)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_ribbon(data = df_CIs2,
                aes(ymin = Qua_low, ymax = Qua_high), alpha = 0.5) +
    geom_line(data = df_CIs2) +
    facet_wrap(~ covNames, scales = "free_x",
               strip.position = "bottom") +
    labs(y = "Native richness",
         title = plot_title) +
    def_theme_paper +
    theme(axis.title.x = element_blank(),
          strip.placement = "outside")
  
  # return
  return(list(fig_out, schange))
  
}

# run function
hydr_rich_out <- plot_grad_fun(hydr_fit, hydr_dat, "(A) hydrilla")
flpl_rich_out <- plot_grad_fun(flpl_fit, flpl_dat, "(B) floating plants")
cubu_rich_out <- plot_grad_fun(cubu_fit, cubu_dat, "(A) Cuban bulrush")
pagr_rich_out <- plot_grad_fun(pagr_fit, pagr_dat, "(B) paragrass")
torp_rich_out <- plot_grad_fun(torp_fit, torp_dat, "(C) torpedograss")

# focal figure
hydr_rich_fig <- hydr_rich_out[[1]]
flpl_rich_fig <- flpl_rich_out[[1]]
foc_rich_fig <- hydr_rich_fig / flpl_rich_fig

ggsave("output/fwc_native_plant_richness_focal_invaders.png", foc_rich_fig,
       device = "png", width = 5, height = 5, units = "in")

# non-focal figure
cubu_rich_fig <- cubu_rich_out[[1]]
pagr_rich_fig <- pagr_rich_out[[1]]
torp_rich_fig <- torp_rich_out[[1]]
non_foc_rich_fig <- cubu_rich_fig / pagr_rich_fig / torp_rich_fig

ggsave("output/fwc_native_plant_richness_non_focal_invaders.png", non_foc_rich_fig,
       device = "png", width = 5, height = 7.5, units = "in")

# focal table of richness effects
foc_rich_tab <- hydr_rich_out[[2]] %>%
  mutate(Invasive = "hydrilla") %>%
  full_join(flpl_rich_out[[2]] %>%
              mutate(Invasive = "floating plants"))

write_csv(foc_rich_tab, "output/fwc_native_plant_richness_focal_invaders.csv")

non_foc_rich_tab <- cubu_rich_out[[2]] %>%
  mutate(Invasive = "Cuban bulrush") %>%
  full_join(pagr_rich_out[[2]] %>%
              mutate(Invasive = "paragrass")) %>%
  full_join(torp_rich_out[[2]] %>%
              mutate(Invasive = "torpedograss"))

write_csv(non_foc_rich_tab, "output/fwc_native_plant_richness_non_focal_invaders.csv")
