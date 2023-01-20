#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(inspectdf) # inspect_cor
library(plm) # panel data models
library(sandwich) # vcovHC
library(lmtest) # coeftest
library(car) # logit
library(glmmTMB) # glmmTMB
library(broom) # tidy

# figure settings
source("code/settings/figure_settings.R")

# functions
source("code/generic-functions/model_structure_comparison.R")

# import data
dat <- read_csv("intermediate-data/FWC_common_native_richness_invaded_data_formatted.csv")
dat_full <- read_csv("intermediate-data/FWC_common_native_richness_invasive_species_data_formatted.csv")


#### edit data ####

# taxa
inv_taxa <- sort(unique(dat$CommonName))

# loop through taxa
pdf("output/native_richness_time_series_by_taxon.pdf")

for(i in inv_taxa){
  
  # subset data
  subdat <- dat %>% filter(CommonName == i)
  subdat_ctrl <- subdat %>% filter(Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = Richness, color = PermanentID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          labs(x = "Year", y = "Native richness", title = i) +
          def_theme_paper +
          theme(plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()

# split by species (para grass doesn't have enough data)
hydr_dat <- filter(dat, CommonName == "Hydrilla") %>%
  mutate(PercCovered_c = PercCovered - mean(PercCovered))
torp_dat <- filter(dat, CommonName == "Torpedograss") %>%
  mutate(PercCovered_c = PercCovered - mean(PercCovered))
cubu_dat <- filter(dat, CommonName == "Cuban bulrush") %>%
  mutate(PercCovered_c = PercCovered - mean(PercCovered))
flpl_dat <- filter(dat, CommonName == "floating plants") %>%
  mutate(PercCovered_c = PercCovered - mean(PercCovered))

hydr_dat_full <- filter(dat_full, CommonName == "Hydrilla") %>%
  mutate(PercCovered_c = PercCovered - mean(PercCovered))
torp_dat_full <- filter(dat_full, CommonName == "Torpedograss") %>%
  mutate(PercCovered_c = PercCovered - mean(PercCovered))
cubu_dat_full <- filter(dat_full, CommonName == "Cuban bulrush") %>%
  mutate(PercCovered_c = PercCovered - mean(PercCovered))
flpl_dat_full <- filter(dat_full, CommonName == "floating plants") %>%
  mutate(PercCovered_c = PercCovered - mean(PercCovered))


#### initial visualizations ####

# covariate correlations
dat %>%
  select(CommonName, PercCovered, Treated,
         PropTreated, RecentTreatment) %>%
  group_by(CommonName) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & abs(corr) >= 0.4) %>%
  data.frame()

# response distributions
ggplot(dat, aes(x = Richness)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# skewed

ggplot(dat, aes(x = log(Richness + 1))) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = logit(PercCovered), y = LogRichness, color = PermanentID)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")

dat %>%
  filter(logit(PercCovered) < -3) %>%
  ggplot(aes(x = logit(PercCovered), y = LogRichness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = Treated, y = LogRichness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = Lag3Treated, y = LogRichness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# all missing values are from beginning of dataset

ggplot(dat, aes(x = RecentTreatment, y = LogRichness)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  # geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = RecentTreatment, y = LogRichness, color = PermanentID)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")

# established
ggplot(dat_full, aes(x = as.factor(Established), y = LogRichness)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ CommonName, scales = "free")


#### evaluate model structure ####

# function to fit models for each species
mod_structure_fits <- function(dat_in){
  
  # create fixed effects data frame
  # individual-time rows
  # each waterbody is an individual
  dat_fix <- dat_in %>%
    pdata.frame(index = c("PermanentID", "GSYear"))
  
  # simple lm
  mod_lm <- lm(LogRichness ~ PercCovered_c + Treated, data = dat_fix)
  
  # random effects
  mod_ran_loc <- glmmTMB(LogRichness ~ PercCovered_c + Treated + (1|PermanentID), data = dat_fix)
  mod_ran_yr <- glmmTMB(LogRichness ~ PercCovered_c + Treated + (1|GSYear), data = dat_fix)
  mod_ran_loc_yr <- glmmTMB(LogRichness ~ PercCovered_c + Treated + (1|PermanentID) + (1|GSYear), data = dat_fix)
  
  # fixed effects
  mod_fix_loc <- plm(LogRichness ~ PercCovered_c + Treated, data = dat_fix,
                     model = "within")
  mod_fix_loc_yr <- plm(LogRichness ~ PercCovered_c + Treated, data = dat_fix,
                        model = "within", effect = "twoways")
  
  # use last treatment
  mod_fix_rec <- plm(LogRichness ~ PercCovered_c + RecentTreatment, data = dat_fix,
                     model = "within", effect = "twoways")
  
  # treatment-plant interactions
  mod_fix_int <- plm(LogRichness ~ PercCovered_c * Treated, data = dat_fix,
                     model = "within", effect = "twoways")
  mod_fix_rec_int <- plm(LogRichness ~ PercCovered_c * RecentTreatment, data = dat_fix,
                         model = "within", effect = "twoways")
  
  # return list of models
  return(list(lm = mod_lm,
              ran_loc = mod_ran_loc,
              ran_yr = mod_ran_yr,
              ran_loc_yr = mod_ran_loc_yr,
              fix_loc = mod_fix_loc,
              fix_loc_yr = mod_fix_loc_yr,
              rec = mod_fix_rec,
              int = mod_fix_int,
              rec_int = mod_fix_rec_int))
  
}

# fit models for each species
hydr_mod_struc <- mod_structure_fits(hydr_dat)
flpl_mod_struc <- mod_structure_fits(flpl_dat)
cubu_mod_struc <- mod_structure_fits(cubu_dat)
torp_mod_struc <- mod_structure_fits(torp_dat)

# compare model estimates
hydr_mod_comp <- mod_structure_comp(simp_mods = hydr_mod_struc[1], 
                                    ran_mods = hydr_mod_struc[2:4],
                                    fix_mods = hydr_mod_struc[5:9])
flpl_mod_comp <- mod_structure_comp(simp_mods = flpl_mod_struc[1], 
                                    ran_mods = flpl_mod_struc[2:4],
                                    fix_mods = flpl_mod_struc[5:9])
cubu_mod_comp <- mod_structure_comp(simp_mods = cubu_mod_struc[1], 
                                    ran_mods = cubu_mod_struc[2:4],
                                    fix_mods = cubu_mod_struc[5:9])
torp_mod_comp <- mod_structure_comp(simp_mods = torp_mod_struc[1], 
                                    ran_mods = torp_mod_struc[2:4],
                                    fix_mods = torp_mod_struc[5:9])

# combine species
mod_comp <- hydr_mod_comp %>%
  mutate(Species = "hydrilla") %>%
  full_join(flpl_mod_comp %>%
              mutate(Species = "floating plants")) %>%
  full_join(cubu_mod_comp %>%
              mutate(Species = "Cuban bulrush")) %>%
  full_join(torp_mod_comp %>%
              mutate(Species = "torpedograss")) %>%
  mutate(across(!c(coefficients, Species), ~ round(.x, digits = 3))) %>%
  relocate(Species)

write_csv(mod_comp, "output/fwc_native_richness_model_structure_comparison.csv")

# model comparison notes:
# hydrilla treatment increases richness until location is accounted for, then decrease
# hydrilla slightly increases richness
# recent treatment has a weaker effect than treated
# interaction is zero

# accounting for location reduces the effects of floating plants and treatment
# floating decreases then increases richness
# treatment increases it
# very small negative interactions

# accounting for location reduces effects of Cuban bulrush and treatment
# PAC decreases, treatment decreases when location and year are accounted for
# positive interaction

# accounting for location reduces effects of torpedograss
# both increase except for treatment in interaction model
# interaction is very small

# test fixed effects
# have to refit because data need to be accessible (not "dat_fix")
hydr_mod_diff_fix_loc_yr <- plm(LogRichness ~ PercCovered_c * RecentTreatment, 
                                data = hydr_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(hydr_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(hydr_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

flpl_mod_diff_fix_loc_yr <- plm(LogRichness ~ PercCovered_c * RecentTreatment, 
                                data = flpl_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(flpl_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
  plmtest(flpl_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

cubu_mod_diff_fix_loc_yr <- plm(LogRichness ~ PercCovered_c * RecentTreatment, 
                                data = cubu_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(cubu_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(cubu_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

torp_mod_diff_fix_loc_yr <- plm(LogRichness ~ PercCovered_c * RecentTreatment, 
                                data = torp_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(torp_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(torp_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

# use time and waterbody fixed effects


#### finalize models ####

# format data
hydr_dat_fix <- hydr_dat %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
flpl_dat_fix <- flpl_dat %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
cubu_dat_fix <- cubu_dat %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
torp_dat_fix <- torp_dat %>%
  pdata.frame(index = c("PermanentID", "GSYear"))

# fit models
hydr_mod <- plm(LogRichness ~ PercCovered_c * RecentTreatment, data = hydr_dat_fix,
                     index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
flpl_mod <- update(hydr_mod, data = flpl_dat_fix)
cubu_mod <- update(hydr_mod, data = cubu_dat_fix)
torp_mod <- update(hydr_mod, data = torp_dat_fix)

# add fitted values to pdata.frame (important to match rows)
# convert to regular dataframe
hydr_fit <- mutate(hydr_dat_fix, Fitted = as.numeric(hydr_dat_fix$LogRichness - hydr_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_fit <- mutate(flpl_dat_fix, Fitted = as.numeric(flpl_dat_fix$LogRichness - flpl_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit <- mutate(cubu_dat_fix, Fitted = as.numeric(cubu_dat_fix$LogRichness - cubu_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit <- mutate(torp_dat_fix, Fitted = as.numeric(torp_dat_fix$LogRichness - torp_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)

# fitted vs. observed
ggplot(hydr_fit, aes(x = Fitted, y = LogRichness)) + geom_point()
ggplot(flpl_fit, aes(x = Fitted, y = LogRichness)) + geom_point()
ggplot(cubu_fit, aes(x = Fitted, y = LogRichness)) + geom_point()
ggplot(torp_fit, aes(x = Fitted, y = LogRichness)) + geom_point()

# export models
save(hydr_mod, file = "output/fwc_hydrilla_native_richness_model.rda")
save(flpl_mod, file = "output/fwc_floating_plant_native_richness_model.rda")
save(cubu_mod, file = "output/fwc_cuban_bulrush_native_richness_model.rda")
save(torp_mod, file = "output/fwc_torpedograss_native_richness_model.rda")

# process model SE
mod_se_fun <- function(model, dat, spp, var_type = "HC3"){
  
  dat_out <- tidy(coeftest(model, vcov = vcovHC, type = var_type)) %>%
    mutate(R2 = r.squared(model),
           Invasive = spp,
           Waterbodies = n_distinct(dat$PermanentID),
           Years = n_distinct(dat$GSYear),
           N = Waterbodies * Years)
  
  return(dat_out)
  
}

# combine SE tables
foc_mod_se <- mod_se_fun(hydr_mod, hydr_dat, "hydrilla") %>%
  full_join(mod_se_fun(flpl_mod, flpl_dat, "floating plants")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "PercCovered_c",
                           "management" = "RecentTreatment",
                           "invasive PAC:management" = "PercCovered_c:RecentTreatment")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive)

non_foc_mod_se <- mod_se_fun(cubu_mod, cubu_dat, "Cuban bulrush") %>%
  full_join(mod_se_fun(torp_mod, torp_dat, "torpedograss")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "PercCovered_c",
                           "management" = "RecentTreatment",
                           "invasive PAC:management" = "PercCovered_c:RecentTreatment")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive)


#### binary models ####

# format data
hydr_dat_full_fix <- hydr_dat_full %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
flpl_dat_full_fix <- flpl_dat_full %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
cubu_dat_full_fix <- cubu_dat_full %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
torp_dat_full_fix <- torp_dat_full %>%
  pdata.frame(index = c("PermanentID", "GSYear"))

# fit models
hydr_mod_full <- plm(LogRichness ~ Established, data = hydr_dat_full_fix, 
                          model = "random", random.method="walhus")
flpl_mod_full <- update(hydr_mod_full, data = flpl_dat_full)
cubu_mod_full <- update(hydr_mod_full, data = cubu_dat_full)
torp_mod_full <- update(hydr_mod_full, data = torp_dat_full)

# add fitted values to pdata.frame (important to match rows)
# convert to regular dataframe
hydr_fit_full <- mutate(hydr_dat_full_fix, Fitted = as.numeric(hydr_dat_full_fix$LogRichness - hydr_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_fit_full <- mutate(flpl_dat_full_fix, Fitted = as.numeric(flpl_dat_full_fix$LogRichness - flpl_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_full <- mutate(cubu_dat_full_fix, Fitted = as.numeric(cubu_dat_full_fix$LogRichness - cubu_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_full <- mutate(torp_dat_full_fix, Fitted = as.numeric(torp_dat_full_fix$LogRichness - torp_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)

# fitted vs. observed
ggplot(hydr_fit_full, aes(x = Fitted, y = LogRichness)) + geom_point()
ggplot(flpl_fit_full, aes(x = Fitted, y = LogRichness)) + geom_point()
ggplot(cubu_fit_full, aes(x = Fitted, y = LogRichness)) + geom_point()
ggplot(torp_fit_full, aes(x = Fitted, y = LogRichness)) + geom_point()

# export models
save(hydr_mod_full, file = "output/fwc_hydrilla_native_richness_binary_model.rda")
save(flpl_mod_full, file = "output/fwc_floating_plant_native_richness_binary_model.rda")
save(cubu_mod_full, file = "output/fwc_cuban_bulrush_native_richness_binary_model.rda")
save(torp_mod_full, file = "output/fwc_torpedograss_native_richness_binary_model.rda")

# combine SE tables
foc_mod_full_se <- mod_se_fun(hydr_mod_full, hydr_dat_full, "hydrilla") %>%
  full_join(mod_se_fun(flpl_mod_full, flpl_dat_full, "floating plants")) %>%
  mutate(term = fct_recode(term,
                           "established" = "Established",
                           "intercept" = "(Intercept)")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive)

non_foc_mod_full_se <- mod_se_fun(cubu_mod_full, cubu_dat_full, "Cuban bulrush") %>%
  full_join(mod_se_fun(torp_mod_full, torp_dat_full, "torpedograss")) %>%
  mutate(term = fct_recode(term,
                           "established" = "Established",
                           "intercept" = "(Intercept)")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive)

# save
write_csv(foc_mod_full_se, "output/fwc_focal_native_richness_binary_model_summary.csv")
write_csv(non_foc_mod_full_se, "output/fwc_non_focal_native_richness_binary_model_summary.csv")


#### values for text ####

# extract intercept from fixed effects models
intercept_fun <- function(models, spp){
  
  dat_out <- tibble(Invasive = spp,
                    Intercept = sapply(models, function(x) mean(fixef(x)))) %>%
    mutate(Intercept = exp(Intercept))
  
  return(dat_out)
  
}

# combine model summary info
foc_mod_sum <- foc_mod_se %>%
  left_join(tibble(Invasive = "hydrilla",
                   Intercept = exp(mean(fixef(hydr_mod)))) %>%
              full_join(tibble(Invasive = "floating plants",
                               Intercept = exp(mean(fixef(flpl_mod)))))) %>%
  relocate(Intercept, .before = "R2")

write_csv(foc_mod_sum, "output/fwc_focal_native_richness_model_summary.csv")

non_foc_mod_sum <- non_foc_mod_se %>%
  left_join(tibble(Invasive = "Cuban bulrush",
                   Intercept = exp(mean(fixef(cubu_mod)))) %>%
              full_join(tibble(Invasive = "torpedograss",
                               Intercept = exp(mean(fixef(torp_mod)))))) %>%
  relocate(Intercept, .before = "R2")

write_csv(non_foc_mod_sum, "output/fwc_non_focal_native_richness_model_summary.csv")

# combine binary and continuous models
foc_mean_summary <- foc_mod_full_se %>%
  select(Invasive, Term, Estimate) %>%
  pivot_wider(names_from = Term,
              values_from = Estimate) %>%
  rename(Intercept = intercept,
         Estimate = established) %>%
  mutate(Term = "established",
         Model = "binary",
         Intercept = exp(Intercept)) %>%
  left_join(foc_mod_full_se) %>%
  full_join(foc_mod_sum %>%
              mutate(Model = "continuous")) %>%
  relocate(c(Model, Term), .after = "Invasive") %>%
  relocate(Intercept, .before = "R2") %>%
  arrange(Invasive, Model, Term)

non_foc_mean_summary <- non_foc_mod_full_se %>%
  select(Invasive, Term, Estimate) %>%
  pivot_wider(names_from = Term,
              values_from = Estimate) %>%
  rename(Intercept = intercept,
         Estimate = established) %>%
  mutate(Term = "established",
         Model = "binary",
         Intercept = exp(Intercept)) %>%
  left_join(non_foc_mod_full_se) %>%
  full_join(non_foc_mod_sum %>%
              mutate(Model = "continuous")) %>%
  relocate(c(Model, Term), .after = "Invasive") %>%
  relocate(Intercept, .before = "R2") %>%
  arrange(Invasive, Model, Term)

# save
write_csv(foc_mean_summary, "output/fwc_focal_native_richness_models.csv")
write_csv(non_foc_mean_summary, "output/fwc_non_focal_native_richness_models.csv")


#### model prediction figures ####

# combine data
# extract fixed effects and coefficients from models
# separate predictions for invasive-only and treatment-only
foc_fit_dat <- hydr_dat %>%
  full_join(flpl_dat) %>%
  select(CommonName, PermanentID, GSYear, RecentTreatment, PercCovered, PercCovered_c, Richness, LogRichness) %>%
  full_join(tibble(PermanentID = names(fixef(hydr_mod)),
                   fixef = as.numeric(fixef(hydr_mod)),
                   coefPAC = as.numeric(coef(hydr_mod)[1]),
                   coefTreat = as.numeric(coef(hydr_mod)[2]),
                   CommonName = "Hydrilla",
                   PanelNamePAC = "(A) hydrilla",
                   PanelNameTreat = "(A) hydrilla management") %>%
              full_join(tibble(PermanentID = names(fixef(flpl_mod)),
                               fixef = as.numeric(fixef(flpl_mod)),
                               coefPAC = as.numeric(coef(flpl_mod)[1]),
                               coefTreat = as.numeric(coef(flpl_mod)[2]),
                               CommonName = "floating plants",
                               PanelNamePAC = "(B) floating plants",
                               PanelNameTreat = "(B) floating plant management"))) %>%
  mutate(FittedPAC = fixef + coefPAC * PercCovered_c, # PAC-only effect
         RichnessPAC = exp(FittedPAC) - 1,
         FittedTreat = fixef + coefTreat * RecentTreatment, # treatment-only effect
         RichnessTreat = exp(FittedTreat) - 1) # %>%
  # group_by(CommonName) %>%
  # mutate(BinPAC = cut_number(log(PercCovered + 1), n = 3)) %>%
  # group_by(CommonName, BinPAC) %>%
  # mutate(BinPACMean = cut_mean(BinPAC)) %>%
  # ungroup()

  non_foc_fit_dat <- cubu_dat %>%
  full_join(torp_dat) %>%
  select(CommonName, PermanentID, GSYear, RecentTreatment, PercCovered, PercCovered_c, Richness, LogRichness) %>%
  full_join(tibble(PermanentID = names(fixef(cubu_mod)),
                   fixef = as.numeric(fixef(cubu_mod)),
                   coefPAC = as.numeric(coef(cubu_mod)[1]),
                   coefTreat = as.numeric(coef(cubu_mod)[2]),
                   CommonName = "Cuban bulrush",
                   PanelNamePAC = "(A) Cuban bulrush",
                   PanelNameTreat = "(A) Cuban bulrush management") %>%
              full_join(tibble(PermanentID = names(fixef(torp_mod)),
                               fixef = as.numeric(fixef(torp_mod)),
                               coefPAC = as.numeric(coef(torp_mod)[1]),
                               coefTreat = as.numeric(coef(torp_mod)[2]),
                               CommonName = "Torpedograss",
                               PanelNamePAC = "(B) torpedograss",
                               PanelNameTreat = "(B) torpedograss management"))) %>%
  mutate(FittedPAC = fixef + coefPAC * PercCovered_c, # PAC-only effect
         RichnessPAC = exp(FittedPAC) - 1,
         FittedTreat = fixef + coefTreat * RecentTreatment, # treatment-only effect
         RichnessTreat = exp(FittedTreat) - 1) # %>%
  # group_by(CommonName) %>%
  # mutate(BinPAC = cut_number(log(PercCovered + 1), n = 3)) %>%
  # group_by(CommonName, BinPAC) %>%
  # mutate(BinPACMean = cut_mean(BinPAC)) %>%
  # ungroup()

# raw data (each waterbody represented)
ggplot(foc_fit_dat, aes(x = PercCovered, color = PermanentID)) +
  geom_point(aes(y = Richness), alpha = 0.3) +
  geom_line(aes(y = RichnessPAC)) +
  facet_wrap(~ PanelNamePAC, scales = "free") +
  labs(x = "PAC", y = "Native richness") +
  def_theme_paper +
  theme(legend.position = "none")

ggplot(foc_fit_dat, aes(x = log(PercCovered + 1), color = PermanentID)) +
  geom_point(aes(y = LogRichness), alpha = 0.3) +
  geom_line(aes(y = FittedPAC)) +
  facet_wrap(~ PanelNamePAC, scales = "free") +
  labs(x = "Log PAC", y = "Log Native richness") +
  def_theme_paper +
  theme(legend.position = "none")

ggplot(foc_fit_dat, aes(x = RecentTreatment, color = PermanentID)) +
  geom_point(aes(y = Richness), alpha = 0.3) +
  geom_line(aes(y = RichnessTreat)) +
  facet_wrap(~ PanelNameTreat, scales = "free") +
  labs(x = "Recent treatment", y = "Native richness") +
  def_theme_paper +
  theme(legend.position = "none")

#### START HERE ####
# binning didn't work above -- too small of a range or something

# summarize raw data
foc_pred_PAC_fig <- ggplot(foc_fit_dat, aes(x = log(Lag3AvgPercCovered + 1))) +
  geom_point(aes(y = FittedPAC, color = PermanentID), size = 0.5, alpha = 0.1) +
  geom_line(aes(y = FittedPAC, color = PermanentID), alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0,
               aes(x = BinPACMean, y = RichnessDiff)) +
  stat_summary(geom = "point", fun = "mean", size = 2,
               aes(x = BinPACMean, y = RichnessDiff)) +
  facet_wrap(~ PanelNamePAC, scales = "free") +
  labs(x = "3-year average PAC (log[x + 1])", y = "Annual difference in native richness") +
  scale_color_manual(values = rep(kelly(), 9)) +
  def_theme_paper +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, color = "black", hjust = 0))

foc_pred_treat_fig <- ggplot(foc_fit_dat, aes(x = Treated)) +
  geom_point(aes(y = FittedTreat, color = PermanentID), size = 0.5, alpha = 0.1) +
  geom_line(aes(y = FittedTreat, color = PermanentID), alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0,
               aes(y = RichnessDiff)) +
  stat_summary(geom = "point", fun = "mean", size = 2,
               aes(y = RichnessDiff)) +
  facet_wrap(~ PanelNameTreat, scales = "free") +
  labs(x = "Years managed (out of 3)", y = "Annual difference in native richness") +
  scale_color_manual(values = rep(kelly(), 9)) +
  def_theme_paper +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, color = "black", hjust = 0))

non_foc_pred_PAC_fig <- ggplot(non_foc_fit_dat, aes(x = log(Lag3AvgPercCovered + 1))) +
  geom_point(aes(y = FittedPAC, color = PermanentID), size = 0.5, alpha = 0.1) +
  geom_line(aes(y = FittedPAC, color = PermanentID), alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0,
               aes(x = BinPACMean, y = RichnessDiff)) +
  stat_summary(geom = "point", fun = "mean", size = 2,
               aes(x = BinPACMean, y = RichnessDiff)) +
  facet_wrap(~ PanelNamePAC, scales = "free") +
  labs(x = "3-year average PAC (log[x + 1])", y = "Annual difference in native richness") +
  scale_color_manual(values = rep(kelly(), 9)) +
  def_theme_paper +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, color = "black", hjust = 0))

non_foc_pred_treat_fig <- ggplot(non_foc_fit_dat, aes(x = Treated)) +
  geom_point(aes(y = FittedTreat, color = PermanentID), size = 0.5, alpha = 0.1) +
  geom_line(aes(y = FittedTreat, color = PermanentID), alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0,
               aes(y = RichnessDiff)) +
  stat_summary(geom = "point", fun = "mean", size = 2,
               aes(y = RichnessDiff)) +
  facet_wrap(~ PanelNameTreat, scales = "free") +
  labs(x = "Years managed (out of 3)", y = "Annual difference in native richness") +
  scale_color_manual(values = rep(kelly(), 9)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3)) +
  def_theme_paper +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, color = "black", hjust = 0))

# save
ggsave("output/fwc_focal_invasive_native_richness_PAC_prediction.png", foc_pred_PAC_fig,
       device = "png", width = 6.5, height = 2.5, units = "in")

ggsave("output/fwc_focal_invasive_native_richness_treatment_prediction.png", foc_pred_treat_fig,
       device = "png", width = 6.5, height = 2.5, units = "in")

ggsave("output/fwc_non_focal_invasive_native_richness_PAC_prediction.png", non_foc_pred_PAC_fig,
       device = "png", width = 6.5, height = 2.5, units = "in")

ggsave("output/fwc_non_focal_invasive_native_richness_treatment_prediction.png", non_foc_pred_treat_fig,
       device = "png", width = 6.5, height = 2.5, units = "in")

