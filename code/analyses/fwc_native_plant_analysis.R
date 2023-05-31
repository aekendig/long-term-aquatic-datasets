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


#### edit data ####

# split data by invasive species
hydr_dat <- filter(dat, CommonName == "Hydrilla") %>%
  arrange(AreaOfInterestID, GSYear, TaxonName)
wale_dat <- filter(dat, CommonName == "Water lettuce") %>%
  arrange(AreaOfInterestID, GSYear, TaxonName)
wahy_dat <- filter(dat, CommonName == "Water hyacinth") %>%
  arrange(AreaOfInterestID, GSYear, TaxonName)
cubu_dat <- filter(dat, CommonName == "Cuban bulrush") %>%
  arrange(AreaOfInterestID, GSYear, TaxonName)
pagr_dat <- filter(dat, CommonName == "Para grass") %>%
  arrange(AreaOfInterestID, GSYear, TaxonName)
torp_dat <- filter(dat, CommonName == "Torpedograss") %>%
  arrange(AreaOfInterestID, GSYear, TaxonName)

# site-by-species matrices
# remove taxa with no presences
hydr_mat <- hydr_dat %>%
  select(AreaOfInterestID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
wale_mat <- wale_dat %>%
  select(AreaOfInterestID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
wahy_mat <- wahy_dat %>%
  select(AreaOfInterestID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
cubu_mat <- cubu_dat %>%
  select(AreaOfInterestID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
pagr_mat <- pagr_dat %>%
  select(AreaOfInterestID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
torp_mat <- torp_dat %>%
  select(AreaOfInterestID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()

# covariates
hydr_cov <- hydr_dat %>%
  distinct(AreaOfInterestID, GSYear, PrevPercCovered, Treated) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  data.frame()
wale_cov <- wale_dat %>%
  distinct(AreaOfInterestID, GSYear, PrevPercCovered, Treated) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  data.frame()
wahy_cov <- wahy_dat %>%
  distinct(AreaOfInterestID, GSYear, PrevPercCovered, Treated) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  data.frame()
cubu_cov <- cubu_dat %>%
  distinct(AreaOfInterestID, GSYear, PrevPercCovered, Treated) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  data.frame()
pagr_cov <- pagr_dat %>%
  distinct(AreaOfInterestID, GSYear, PrevPercCovered, Treated) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  data.frame()
torp_cov <- torp_dat %>%
  distinct(AreaOfInterestID, GSYear, PrevPercCovered, Treated) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  select(-c(AreaOfInterestID, GSYear)) %>%
  data.frame()

# study design
hydr_stud <- hydr_dat %>%
  distinct(AreaOfInterestID, GSYear) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  mutate(AreaOfInterestID = as.factor(AreaOfInterestID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
wale_stud <- wale_dat %>%
  distinct(AreaOfInterestID, GSYear) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  mutate(AreaOfInterestID = as.factor(AreaOfInterestID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
wahy_stud <- wahy_dat %>%
  distinct(AreaOfInterestID, GSYear) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  mutate(AreaOfInterestID = as.factor(AreaOfInterestID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
cubu_stud <- cubu_dat %>%
  distinct(AreaOfInterestID, GSYear) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  mutate(AreaOfInterestID = as.factor(AreaOfInterestID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
pagr_stud <- pagr_dat %>%
  distinct(AreaOfInterestID, GSYear) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  mutate(AreaOfInterestID = as.factor(AreaOfInterestID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
torp_stud <- torp_dat %>%
  distinct(AreaOfInterestID, GSYear) %>%
  arrange(AreaOfInterestID, GSYear) %>%
  mutate(AreaOfInterestID = as.factor(AreaOfInterestID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()

# traits
hydr_hab <- hydr_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(hydr_hab) <- hydr_hab$TaxonName
hydr_trait <- select(hydr_hab, -TaxonName)
wale_hab <- wale_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(wale_hab) <- wale_hab$TaxonName
wale_trait <- select(wale_hab, -TaxonName)
wahy_hab <- wahy_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(wahy_hab) <- wahy_hab$TaxonName
wahy_trait <- select(wahy_hab, -TaxonName)
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
                XFormula = ~PrevPercCovered * Treated,
                TrData = hydr_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = hydr_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = hydr_stud$GSYear), 
                                 AreaOfInterestID = HmscRandomLevel(units = hydr_stud$AreaOfInterestID)))

wale_mod = Hmsc(Y = wale_mat,
                XData = wale_cov, 
                XFormula = ~PrevPercCovered * Treated,
                TrData = wale_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = wale_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = wale_stud$GSYear), 
                                 AreaOfInterestID = HmscRandomLevel(units = wale_stud$AreaOfInterestID)))

wahy_mod = Hmsc(Y = wahy_mat,
                XData = wahy_cov, 
                XFormula = ~PrevPercCovered * Treated,
                TrData = wahy_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = wahy_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = wahy_stud$GSYear), 
                                 AreaOfInterestID = HmscRandomLevel(units = wahy_stud$AreaOfInterestID)))

cubu_mod = Hmsc(Y = cubu_mat,
                XData = cubu_cov, 
                XFormula = ~PrevPercCovered * Treated,
                TrData = cubu_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = cubu_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = cubu_stud$GSYear), 
                                 AreaOfInterestID = HmscRandomLevel(units = cubu_stud$AreaOfInterestID)))

pagr_mod = Hmsc(Y = pagr_mat,
                XData = pagr_cov, 
                XFormula = ~PrevPercCovered * Treated,
                TrData = pagr_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = pagr_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = pagr_stud$GSYear), 
                                 AreaOfInterestID = HmscRandomLevel(units = pagr_stud$AreaOfInterestID)))

torp_mod = Hmsc(Y = torp_mat,
                XData = torp_cov, 
                XFormula = ~PrevPercCovered * Treated,
                TrData = torp_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = torp_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = torp_stud$GSYear), 
                                 AreaOfInterestID = HmscRandomLevel(units = torp_stud$AreaOfInterestID)))

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

wale_fit = sampleMcmc(wale_mod,
                      thin = thin,
                      samples = samples,
                      transient = transient,
                      nChains = nChains,
                      nParallel = nChains)
save(wale_fit, file = "output/fwc_native_plant_water_lettuce_hmsc.rda")

wahy_fit = sampleMcmc(wahy_mod,
                      thin = thin,
                      samples = samples,
                      transient = transient,
                      nChains = nChains,
                      nParallel = nChains)
save(wahy_fit, file = "output/fwc_native_plant_water_hyacinth_hmsc.rda")

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
load("output/fwc_native_plant_water_lettuce_hmsc.rda")
load("output/fwc_native_plant_water_hyacinth_hmsc.rda")
load("output/fwc_native_plant_Cuban_bulrush_hmsc.rda")
load("output/fwc_native_plant_paragrass_hmsc.rda")
load("output/fwc_native_plant_torpedograss_hmsc.rda")


#### evaluate fit ####

# posterior samples
hydr_post = convertToCodaObject(hydr_fit)
wale_post = convertToCodaObject(wale_fit)
wahy_post = convertToCodaObject(wahy_fit)
cubu_post = convertToCodaObject(cubu_fit)
pagr_post = convertToCodaObject(pagr_fit)
torp_post = convertToCodaObject(torp_fit)

# effective sample sizes
# total samples = samples * chains
hist(effectiveSize(hydr_post$Beta))
hist(effectiveSize(wale_post$Beta))
hist(effectiveSize(wahy_post$Beta))
hist(effectiveSize(cubu_post$Beta))
hist(effectiveSize(pagr_post$Beta)) # almost uniform
hist(effectiveSize(torp_post$Beta))

# Gelman diagnostics (potential scale reduction factors)
# one indicates better convergence among chains
hist(gelman.diag(hydr_post$Beta, multivariate=FALSE)$psrf)
hist(gelman.diag(wale_post$Beta, multivariate=FALSE)$psrf)
hist(gelman.diag(wahy_post$Beta, multivariate=FALSE)$psrf)
hist(gelman.diag(cubu_post$Beta, multivariate=FALSE)$psrf)
hist(gelman.diag(pagr_post$Beta, multivariate=FALSE)$psrf)
hist(gelman.diag(torp_post$Beta, multivariate=FALSE)$psrf)

# explanatory power, compares posterior distribution to observed values
hydr_preds = computePredictedValues(hydr_fit)
hydr_eval = evaluateModelFit(hM = hydr_fit, predY = hydr_preds)
hist(hydr_eval$RMSE)
hist(hydr_eval$TjurR2)
hist(hydr_eval$AUC)

wale_preds = computePredictedValues(wale_fit)
wale_eval = evaluateModelFit(hM = wale_fit, predY = wale_preds)
hist(wale_eval$RMSE)
hist(wale_eval$TjurR2)
hist(wale_eval$AUC)

wahy_preds = computePredictedValues(wahy_fit)
wahy_eval = evaluateModelFit(hM = wahy_fit, predY = wahy_preds)
hist(wahy_eval$RMSE)
hist(wahy_eval$TjurR2)
hist(wahy_eval$AUC)

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

wale_sum <- as_tibble(summary(wale_post$Beta)$statistics) %>%
  mutate(species = rep(colnames(wale_mat), each = 4),
         covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                         ncol(wale_mat))) %>%
  full_join(as_tibble(summary(wale_post$Beta)$quantiles) %>%
              mutate(species = rep(colnames(wale_mat), each = 4),
                     covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                                     ncol(wale_mat)))) %>%
  mutate(invader = "water lettuce")

wahy_sum <- as_tibble(summary(wahy_post$Beta)$statistics) %>%
  mutate(species = rep(colnames(wahy_mat), each = 4),
         covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                         ncol(wahy_mat))) %>%
  full_join(as_tibble(summary(wahy_post$Beta)$quantiles) %>%
              mutate(species = rep(colnames(wahy_mat), each = 4),
                     covariate = rep(c("intercept", "invasive plant PAC", "management", "interaction"), 
                                     ncol(wahy_mat)))) %>%
  mutate(invader = "water hyacinth")

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
foc_coef_sum <- full_join(hydr_sum, wale_sum) %>%
  full_join(wahy_sum) %>%
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
         cov2 = case_when(covariate == "invasive plant PAC" ~ paste(invader, "PAC"),
                          covariate == "management" ~ paste(invader, "mgmt."),
                          TRUE ~ as.character(covariate)) %>%
           fct_relevel("intercept", "hydrilla PAC", "water hyacinth PAC", "water lettuce PAC",
                       "hydrilla mgmt.", "water hyacinth mgmt.", "water lettuce mgmt."),
         invader = fct_relevel(invader, "hydrilla"),
         habitat_num = as.numeric(as.factor(Habitat)),
         sp_num = habitat_num * 100 + as.numeric(as.factor(species)),
         species = fct_reorder(species, sp_num, .desc = T),
         sig = case_when(quant_2_5 > 0 & quant_97_5 > 0 ~ "yes",
                         quant_2_5 < 0 & quant_97_5 < 0 ~ "yes",
                         TRUE ~ "no"),
         Habitat = tolower(Habitat))

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
         cov2 = case_when(covariate == "invasive plant PAC" ~ paste(invader, "PAC"),
                          covariate == "management" ~ paste(invader, "mgmt."),
                        TRUE ~ as.character(covariate)) %>%
           fct_relevel("intercept", "Cuban bulrush PAC", 
                       "paragrass PAC",
                       "torpedograss PAC", 
                       "Cuban bulrush mgmt.", 
                       "paragrass mgmt.",
                       "torpedograss mgmt."),
         habitat_num = as.numeric(as.factor(Habitat)),
         sp_num = habitat_num * 100 + as.numeric(as.factor(species)),
         species = fct_reorder(species, sp_num, .desc = T),
         sig = case_when(quant_2_5 > 0 & quant_97_5 > 0 ~ "yes",
                         quant_2_5 < 0 & quant_97_5 < 0 ~ "yes",
                         TRUE ~ "no"),
         Habitat = tolower(Habitat))

# coefficient plots
hydr_coef_fig <- foc_coef_sum %>%
  filter(covariate != "intercept" & invader == "hydrilla") %>%
  ggplot(aes(x = Mean, y = species, color = Habitat, alpha = sig)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = quant_2_5,
                     xmax = quant_97_5),
                 height = 0) +
  geom_point() +
  facet_grid(~ cov2, scales = "free") +
  scale_color_manual(values = brewer.set2(n = 3), name = "Growth form") +
  scale_alpha_manual(values = c(0.2, 1), name = "CI omits zero") +
  labs(x = "Change in detection z-score ± 95% CI") +
  def_theme_paper +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.spacing.y = unit(-0.4, "cm"),
        axis.text.y = element_text(size = 6.5, face = "italic"),
        axis.title.y = element_blank())

wahy_coef_fig <- hydr_coef_fig %+%
  filter(foc_coef_sum,
         covariate != "intercept" & invader == "water hyacinth")

wale_coef_fig <- hydr_coef_fig %+%
  filter(foc_coef_sum,
         covariate != "intercept" & invader == "water lettuce")

cubu_coef_fig <- hydr_coef_fig %+%
  filter(non_foc_coef_sum,
         covariate != "intercept" & invader == "Cuban bulrush") +
  theme(strip.text = element_text(size = 8))

pagr_coef_fig <- hydr_coef_fig %+%
  filter(non_foc_coef_sum,
         covariate != "intercept" & invader == "paragrass")

torp_coef_fig <- hydr_coef_fig %+%
  filter(non_foc_coef_sum,
         covariate != "intercept" & invader == "torpedograss")

ggsave("output/fwc_native_plant_coefficients_hydrilla_PAC.png", hydr_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_water_hyacinth_PAC.png", wahy_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_water_lettuce_PAC.png", wale_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_cuban_bulrush_PAC.png", cubu_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_paragrass_PAC.png", pagr_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_torpedograss_PAC.png", torp_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")

write_csv(foc_coef_sum, "output/fwc_native_plant_coefficients_focal_invaders_PAC.csv")
write_csv(non_foc_coef_sum, "output/fwc_native_plant_coefficients_non_focal_invaders_PAC.csv")
