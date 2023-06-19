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

# coefficients and net effect of management and cover
coef_fun <- function(post_in, inv_name) {
  
  post_long <- do.call("rbind", post_in$Beta)  %>% # collapse three chain dataframes into one
    as_tibble() %>%
    mutate(sample = 1:n()) %>%
    pivot_longer(cols = -sample, # make covariates and taxa into rows
                 names_to = "var",
                 values_to = "est") %>%
    mutate(covariate = case_when(str_detect(var, "(Intercept)") == T ~ "intercept", # column for covariate
                                 str_detect(var, ":Treated") == T ~ "interaction",
                                 str_detect(var, "PrevPercCovered") == T ~ "cover",
                                 str_detect(var, "Treated") == T ~ "management")) %>%
    mutate(TaxonName = str_split_i(var, "\\), ", 2), # column for taxa
           TaxonName = str_split_i(TaxonName, " \\(S", 1)) %>%
    select(-var) %>%
    pivot_wider(names_from = "covariate", # make covariates into columns
                values_from = "est") %>%
    mutate(cover = 10 * cover, # assume 10% increase in cover
           interaction = 10 * interaction,
           net = cover + management + interaction) %>% # net effect
    pivot_longer(cols = c(intercept, cover, management, interaction, net), # put covariates back into rows
                 names_to = "covariate",
                 values_to = "est") %>%
    group_by(TaxonName, covariate) %>%
    summarize(Mean = mean(est), # summarize by covariate and taxon
              quant_2_5 = quantile(est, 0.025),
              quant_97_5 = quantile(est, 0.975),
              .groups = "drop") %>%
    mutate(invader = inv_name)
  
  return(post_long)
  
}

hydr_sum <- coef_fun(hydr_post, "hydrilla")
wale_sum <- coef_fun(wale_post, "water lettuce")
wahy_sum <- coef_fun(wahy_post, "water hyacinth")
cubu_sum <- coef_fun(cubu_post, "Cuban bulrush")
pagr_sum <- coef_fun(pagr_post, "para grass")
torp_sum <- coef_fun(torp_post, "torpedograss")

# combine
foc_coef_sum <- full_join(hydr_sum, wale_sum) %>%
  full_join(wahy_sum) %>%
  left_join(dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(cov2 = case_when(covariate == "cover" ~ paste(invader, "cover"), # switch between PAC and cover
                          covariate == "management" ~ paste(invader, "mgmt."),
                          covariate == "net" ~ "mgmt. + cover",
                          TRUE ~ as.character(covariate)) %>%
           fct_relevel("intercept",
                       "hydrilla mgmt.", 
                       "water hyacinth mgmt.", 
                       "water lettuce mgmt.", 
                       "hydrilla cover", 
                       "water hyacinth cover", 
                       "water lettuce cover"),
         invader = fct_relevel(invader, "hydrilla"),
         habitat_num = as.numeric(as.factor(Habitat)),
         tax_num = habitat_num * 100 + as.numeric(as.factor(TaxonName)),
         TaxonName = fct_reorder(TaxonName, tax_num, .desc = T),
         sig = case_when(quant_2_5 > 0 & quant_97_5 > 0 ~ "yes",
                         quant_2_5 < 0 & quant_97_5 < 0 ~ "yes",
                         TRUE ~ "no"),
         Habitat = tolower(Habitat)) %>%
  filter(covariate %in% c("management", "cover", "net"))

non_foc_coef_sum <- full_join(cubu_sum, pagr_sum) %>%
  full_join(torp_sum) %>%
  left_join(dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(cov2 = case_when(covariate == "cover" ~ paste(invader, "cover"),
                          covariate == "management" ~ paste(invader, "mgmt."),
                          covariate == "net" ~ "mgmt. + cover",
                        TRUE ~ as.character(covariate)) %>%
           fct_relevel("intercept", 
                       "Cuban bulrush mgmt.", 
                       "para grass mgmt.",
                       "torpedograss mgmt.", 
                       "Cuban bulrush cover", 
                       "para grass cover",
                       "torpedograss cover"),
         habitat_num = as.numeric(as.factor(Habitat)),
         tax_num = habitat_num * 100 + as.numeric(as.factor(TaxonName)),
         TaxonName = fct_reorder(TaxonName, tax_num, .desc = T),
         sig = case_when(quant_2_5 > 0 & quant_97_5 > 0 ~ "yes",
                         quant_2_5 < 0 & quant_97_5 < 0 ~ "yes",
                         TRUE ~ "no"),
         Habitat = tolower(Habitat)) %>%
  filter(covariate %in% c("management", "cover", "net"))

# coefficient plots
hydr_coef_fig <- foc_coef_sum %>%
  filter(invader == "hydrilla") %>%
  ggplot(aes(x = Mean, y = TaxonName, color = Habitat, alpha = sig)) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = quant_2_5,
                     xmax = quant_97_5),
                 height = 0) +
  geom_point() +
  facet_grid(~ cov2, scales = "free") +
  scale_color_manual(values = brewer.set2(n = 3), name = "Growth form") +
  scale_alpha_manual(values = c(0.2, 1), name = "CI omits zero") +
  labs(x = "Change in z-score ± 95% CI") +
  def_theme_paper +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.spacing.y = unit(-0.4, "cm"),
        panel.spacing.x = unit(4, "mm"),
        axis.text.y = element_text(size = 6.5, face = "italic"),
        axis.title.y = element_blank())

wahy_coef_fig <- hydr_coef_fig %+%
  filter(foc_coef_sum, invader == "water hyacinth")

wale_coef_fig <- hydr_coef_fig %+%
  filter(foc_coef_sum, invader == "water lettuce" &
           TaxonName != "Mayaca fluviatilis")

cubu_coef_fig <- hydr_coef_fig %+%
  filter(non_foc_coef_sum, invader == "Cuban bulrush") +
  theme(strip.text = element_text(size = 8))

pagr_coef_fig <- hydr_coef_fig %+%
  filter(non_foc_coef_sum, invader == "para grass")

torp_coef_fig <- hydr_coef_fig %+%
  filter(non_foc_coef_sum, invader == "torpedograss")

ggsave("output/fwc_native_plant_coefficients_hydrilla_cover.png", hydr_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_water_hyacinth_cover.png", wahy_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_water_lettuce_cover.png", wale_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_cuban_bulrush_cover.png", cubu_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_paragrass_cover.png", pagr_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")
ggsave("output/fwc_native_plant_coefficients_torpedograss_cover.png", torp_coef_fig,
       device = "png", width = 6.5, height = 6.5, units = "in")

write_csv(foc_coef_sum, "output/fwc_native_plant_coefficients_focal_invaders_cover.csv")
write_csv(non_foc_coef_sum, "output/fwc_native_plant_coefficients_non_focal_invaders_cover.csv")


#### summary figure ####

# color palette
pal <- c("#0072B2", "#D55E00")

# summarize data
foc_sum_dat <- foc_coef_sum %>%
  mutate(direction = if_else(Mean > 0, "+", "-")) %>%
  group_by(invader, covariate, cov2, sig, direction) %>%
  summarize(taxa = n_distinct(TaxonName)) %>%
  ungroup() %>%
  mutate(taxa_plot = if_else(direction == "-", -1 * taxa, taxa),
         panel_name = case_when(invader == "hydrilla" ~ "(A) hydrilla",
                                invader == "water hyacinth" ~ "(B) water hyacinth",
                                invader == "water lettuce" ~ "(C) water lettuce"),
         # covariate = str_replace(covariate, "cover", "PAC") %>%
         #   fct_relevel("management", "PAC"),
         covariate = fct_relevel(covariate, "management", "cover"),
         sig = fct_relevel(sig, "yes")) %>%
  group_by(invader, covariate, direction) %>%
  mutate(taxa_col = sum(taxa_plot)) %>%
  ungroup() %>%
  mutate(lab_y = if_else(sig == "no", taxa_plot/2,
                         taxa_col - taxa_plot/2))

non_foc_sum_dat <- non_foc_coef_sum %>%
  mutate(direction = if_else(Mean > 0, "+", "-")) %>%
  group_by(invader, covariate, cov2, sig, direction) %>%
  summarize(taxa = n_distinct(TaxonName)) %>%
  ungroup() %>%
  filter(covariate != "intercept") %>%
  mutate(taxa_plot = if_else(direction == "-", -1 * taxa, taxa),
         panel_name = case_when(invader == "Cuban bulrush" ~ "(A) Cuban bulrush",
                                invader == "para grass" ~ "(B) para grass",
                                invader == "torpedograss" ~ "(C) torpedograss"),
         # covariate = str_replace(covariate, "cover", "PAC") %>%
         #   fct_relevel("management", "PAC"),
         covariate = fct_relevel(covariate, "management", "cover"),
         sig = fct_relevel(sig, "yes")) %>%
  group_by(invader, covariate, direction) %>%
  mutate(taxa_col = sum(taxa_plot)) %>%
  ungroup() %>%
  mutate(lab_y = if_else(sig == "no", taxa_plot/2,
                         taxa_col - taxa_plot/2))

# figure function
foc_sum_fig <- ggplot(foc_sum_dat, aes(x = covariate, y = taxa_plot, 
                          fill = direction)) +
    geom_col(aes(alpha = sig)) +
    geom_hline(yintercept = 0) +
  geom_text(aes(label = taxa, y = lab_y), size = 3) +
    facet_wrap(~ panel_name) +
    scale_fill_manual(values = pal, name = "Response direction") +
  scale_alpha_manual(values = c(1, 0.5), name = "CI omits zero") +
  scale_y_continuous(breaks = c(-50, -25, 0, 25, 50),
                    labels = c(50, 25, 0, 25, 50)) +
    labs(y = "Number of native taxa") +
    def_theme_paper +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.spacing.y = unit(-0.4, "cm"),
          axis.title.x = element_blank(),
          strip.text = element_text(hjust = 0))

ggsave("output/fwc_native_plant_focal_summary_figure_cover.png", foc_sum_fig,
       device = "png", width = 6.5, height = 4, units = "in")

non_foc_sum_fig <- ggplot(non_foc_sum_dat, aes(x = covariate, y = taxa_plot, 
                                               fill = direction)) +
  geom_col(aes(alpha = sig)) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = taxa, y = lab_y), size = 3) +
  facet_wrap(~ panel_name) +
  scale_fill_manual(values = pal, name = "Response direction") +
  scale_alpha_manual(values = c(1, 0.5), name = "CI omits zero") +
  scale_y_continuous(breaks = c(-50, -25, 0, 25, 50),
                     labels = c(50, 25, 0, 25, 50)) +
  labs(y = "Number of native taxa") +
  def_theme_paper +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.spacing.y = unit(-0.4, "cm"),
        axis.title.x = element_blank(),
        strip.text = element_text(hjust = 0))

ggsave("output/fwc_native_plant_non_focal_summary_figure_cover.png", non_foc_sum_fig,
       device = "png", width = 6.5, height = 4, units = "in")


#### time between management and survey ####

# calculate time difference
time_diff <- dat %>%
  filter(Treated == 1) %>%
  distinct(AreaOfInterestID, GSYear, CommonName, SurveyDate, MaxTreatmentDate) %>%
  mutate(SurveyTreatDays = as.numeric(SurveyDate - MaxTreatmentDate))

# check that it makes sense
janitor::get_dupes(time_diff, AreaOfInterestID, GSYear, CommonName)
# should return no duplicates

count(time_diff, CommonName)
# should be similar, but higher than richness N
# higher because richness required previous year

# focal taxa
foc_time_diff <- time_diff %>%
  filter(CommonName %in% c("Hydrilla", "Water hyacinth", "Water lettuce")) %>%
  mutate(pan_lab = case_when(CommonName == "Hydrilla" ~ "(A) hydrilla",
                             CommonName == "Water hyacinth" ~ "(B) water hyacinth",
                             CommonName == "Water lettuce" ~ "(C) water lettuce"))



# figure
foc_time_diff_fig <- ggplot(foc_time_diff, aes(x = SurveyTreatDays)) +
  geom_histogram(binwidth = 7) +
  facet_wrap(~ pan_lab, scales = "free") +
  labs(x = "Days between management and plant survey",
       y = "Number of waterbody-year combinations") +
  def_theme_paper +
  theme(strip.text = element_text(hjust = 0))

ggsave("output/fwc_focal_survey_management_time_difference.png", foc_time_diff_fig,
       device = "png", width = 6.5, height = 3, units = "in")
