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
qual <- read_csv("intermediate-data/water_quality_invaded_data_formatted.csv")
qual_full <- read_csv("intermediate-data/water_quality_data_formatted.csv")


#### edit data ####

# focal water quality
dat <- filter(qual, QualityMetric == "TP_ug_L")
dat_full <- filter(qual_full, QualityMetric == "TP_ug_L")

# taxa
inv_taxa <- sort(unique(dat$CommonName))

# loop through taxa
pdf("output/phosphorus_time_series_by_taxon.pdf")

for(i in inv_taxa){
  
  # subset data
  subdat <- dat %>% filter(CommonName == i)
  subdat_ctrl <- subdat %>% filter(Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = QualityMean, color = PermanentID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          labs(x = "Year", y = "total phosphorus (ug/L)", title = i) +
          def_theme_paper +
          theme(plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()

# save data
write_csv(dat, "intermediate-data/FWC_phosphorus_analysis_formatted.csv")

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
ggplot(dat, aes(x = QualityMean)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# skewed

ggplot(dat, aes(x = log(QualityMean))) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = QualityMax)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# skewed

ggplot(dat, aes(x = log(QualityMax))) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = QualityMin)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# skewed

ggplot(dat, aes(x = log(QualityMin))) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = QualityCV)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# skewed

ggplot(dat, aes(x = log(QualityCV))) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

# coefficients and response values
ggplot(dat, aes(x = logit(PercCovered), y = LogQualityMean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = logit(PercCovered), y = LogQualityMean, color = PermanentID)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")

dat %>%
  filter(logit(PercCovered) < -3) %>%
  ggplot(aes(x = logit(PercCovered), y = LogQualityMean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = Treated, y = LogQualityMean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = Lag3Treated, y = LogQualityMean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# all missing values are from beginning of dataset

ggplot(dat, aes(x = RecentTreatment, y = LogQualityMean)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  # geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# P increases the more recently lake has been treated for floating plants

ggplot(dat, aes(x = RecentTreatment, y = LogQualityMean, color = PermanentID)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")
# floating plant trend may be driven by lake/treatment confounding

# established
ggplot(dat_full, aes(x = as.factor(Established), y = LogQualityMean)) +
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
  mod_lm <- lm(LogQualityMean ~ PercCovered_c + Treated, data = dat_fix)
  
  # random effects
  mod_ran_loc <- glmmTMB(LogQualityMean ~ PercCovered_c + Treated + (1|PermanentID), data = dat_fix)
  mod_ran_yr <- glmmTMB(LogQualityMean ~ PercCovered_c + Treated + (1|GSYear), data = dat_fix)
  mod_ran_loc_yr <- glmmTMB(LogQualityMean ~ PercCovered_c + Treated + (1|PermanentID) + (1|GSYear), data = dat_fix)
  
  # fixed effects
  mod_fix_loc <- plm(LogQualityMean ~ PercCovered_c + Treated, data = dat_fix,
                      model = "within")
  mod_fix_loc_yr <- plm(LogQualityMean ~ PercCovered_c + Treated, data = dat_fix,
                         model = "within", effect = "twoways")
  
  # add quadratic term to account for multiple processes
  # mod_fix_quad <- plm(LogQualityMean ~ PercCovered_c + Lag3APCsq + Lag3Treated, data = dat_fix,
  #                     model = "within", effect = "twoways")
  # tried this, need more years to do a squared term
  
  # use last treatment
  mod_fix_rec <- plm(LogQualityMean ~ PercCovered_c + RecentTreatment, data = dat_fix,
                      model = "within", effect = "twoways")
  
  # treatment-plant interactions
  mod_fix_int <- plm(LogQualityMean ~ PercCovered_c * Treated, data = dat_fix,
                      model = "within", effect = "twoways")
  mod_fix_rec_int <- plm(LogQualityMean ~ PercCovered_c * RecentTreatment, data = dat_fix,
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

write_csv(mod_comp, "output/fwc_phosphorus_model_structure_comparison.csv")

# model comparison notes:
# hydrilla and treatment increase P until location is accounted for, then they decrease
# recent treatment has a weaker effect than treated
# interaction is minimal

# accounting for location reduces the effects of floating plants and treatment on P
# both have positive effect
# recent treatment has a stronger effect
# treatment eliminates the positive effect of the plants

# accounting for location and time don't seem to affect estimates much
# recent treatment is also similar to treated
# treatment and plants both reduce P
# synergistic effect of the two (even more reduction)

# accounting for location reduces the effects
# both reduce P
# recent similar to treated
# interaction accounts for most of the treated effects

# test fixed effects (seems like year isn't necessary)
# have to refit because data need to be accessible (not "dat_fix")
hydr_mod_diff_fix_loc_yr <- plm(LogQualityMean ~ PercCovered_c * RecentTreatment, 
                                data = hydr_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(hydr_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(hydr_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

flpl_mod_diff_fix_loc_yr <- plm(LogQualityMean ~ PercCovered_c * RecentTreatment, 
                                data = flpl_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(flpl_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(flpl_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

cubu_mod_diff_fix_loc_yr <- plm(LogQualityMean ~ PercCovered_c * RecentTreatment, 
                                data = cubu_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(cubu_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(cubu_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

torp_mod_diff_fix_loc_yr <- plm(LogQualityMean ~ PercCovered_c * RecentTreatment, 
                                data = torp_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(torp_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(torp_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

# use waterbody fixed effects
# decided to leave in year effects because there should in theory be a year effect


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
hydr_mean_mod <- plm(LogQualityMean ~ PercCovered_c * RecentTreatment, data = hydr_dat_fix,
                       index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
flpl_mean_mod <- update(hydr_mean_mod, data = flpl_dat_fix)
cubu_mean_mod <- update(hydr_mean_mod, data = cubu_dat_fix)
torp_mean_mod <- update(hydr_mean_mod, data = torp_dat_fix)

hydr_max_mod <- plm(LogQualityMax ~ PercCovered_c * RecentTreatment, data = hydr_dat_fix,
                     index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
flpl_max_mod <- update(hydr_max_mod, data = flpl_dat_fix)
cubu_max_mod <- update(hydr_max_mod, data = cubu_dat_fix)
torp_max_mod <- update(hydr_max_mod, data = torp_dat_fix)

hydr_min_mod <- plm(LogQualityMin ~ PercCovered_c * RecentTreatment, data = hydr_dat_fix,
                     index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
flpl_min_mod <- update(hydr_min_mod, data = flpl_dat_fix)
cubu_min_mod <- update(hydr_min_mod, data = cubu_dat_fix)
torp_min_mod <- update(hydr_min_mod, data = torp_dat_fix)

hydr_cv_mod <- plm(LogQualityCV ~ PercCovered_c * RecentTreatment, data = hydr_dat_fix,
                     index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
flpl_cv_mod <- update(hydr_cv_mod, data = flpl_dat_fix)
cubu_cv_mod <- update(hydr_cv_mod, data = cubu_dat_fix)
torp_cv_mod <- update(hydr_cv_mod, data = torp_dat_fix)

# add fitted values to pdata.frame (important to match rows)
# convert to regular dataframe
hydr_mean_fit <- mutate(hydr_dat_fix, Fitted = as.numeric(hydr_dat_fix$LogQualityMean - hydr_mean_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_mean_fit <- mutate(flpl_dat_fix, Fitted = as.numeric(flpl_dat_fix$LogQualityMean - flpl_mean_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_mean_fit <- mutate(cubu_dat_fix, Fitted = as.numeric(cubu_dat_fix$LogQualityMean - cubu_mean_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_mean_fit <- mutate(torp_dat_fix, Fitted = as.numeric(torp_dat_fix$LogQualityMean - torp_mean_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_max_fit <- mutate(hydr_dat_fix, Fitted = as.numeric(hydr_dat_fix$LogQualityMax - hydr_max_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_max_fit <- mutate(flpl_dat_fix, Fitted = as.numeric(flpl_dat_fix$LogQualityMax - flpl_max_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_max_fit <- mutate(cubu_dat_fix, Fitted = as.numeric(cubu_dat_fix$LogQualityMax - cubu_max_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_max_fit <- mutate(torp_dat_fix, Fitted = as.numeric(torp_dat_fix$LogQualityMax - torp_max_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_min_fit <- mutate(hydr_dat_fix, Fitted = as.numeric(hydr_dat_fix$LogQualityMin - hydr_min_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_min_fit <- mutate(flpl_dat_fix, Fitted = as.numeric(flpl_dat_fix$LogQualityMin - flpl_min_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_min_fit <- mutate(cubu_dat_fix, Fitted = as.numeric(cubu_dat_fix$LogQualityMin - cubu_min_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_min_fit <- mutate(torp_dat_fix, Fitted = as.numeric(torp_dat_fix$LogQualityMin - torp_min_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_cv_fit <- mutate(hydr_dat_fix, Fitted = as.numeric(hydr_dat_fix$LogQualityCV - hydr_cv_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_cv_fit <- mutate(flpl_dat_fix, Fitted = as.numeric(flpl_dat_fix$LogQualityCV - flpl_cv_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_cv_fit <- mutate(cubu_dat_fix, Fitted = as.numeric(cubu_dat_fix$LogQualityCV - cubu_cv_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_cv_fit <- mutate(torp_dat_fix, Fitted = as.numeric(torp_dat_fix$LogQualityCV - torp_cv_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)

# fitted vs. observed
ggplot(hydr_mean_fit, aes(x = Fitted, y = LogQualityMean)) + geom_point()
ggplot(flpl_mean_fit, aes(x = Fitted, y = LogQualityMean)) + geom_point()
ggplot(cubu_mean_fit, aes(x = Fitted, y = LogQualityMean)) + geom_point()
ggplot(torp_mean_fit, aes(x = Fitted, y = LogQualityMean)) + geom_point()

ggplot(hydr_max_fit, aes(x = Fitted, y = LogQualityMax)) + geom_point()
ggplot(flpl_max_fit, aes(x = Fitted, y = LogQualityMax)) + geom_point()
ggplot(cubu_max_fit, aes(x = Fitted, y = LogQualityMax)) + geom_point()
ggplot(torp_max_fit, aes(x = Fitted, y = LogQualityMax)) + geom_point()

ggplot(hydr_min_fit, aes(x = Fitted, y = LogQualityMin)) + geom_point()
ggplot(flpl_min_fit, aes(x = Fitted, y = LogQualityMin)) + geom_point()
ggplot(cubu_min_fit, aes(x = Fitted, y = LogQualityMin)) + geom_point()
ggplot(torp_min_fit, aes(x = Fitted, y = LogQualityMin)) + geom_point()

ggplot(hydr_cv_fit, aes(x = Fitted, y = LogQualityCV)) + geom_point()
ggplot(flpl_cv_fit, aes(x = Fitted, y = LogQualityCV)) + geom_point()
ggplot(cubu_cv_fit, aes(x = Fitted, y = LogQualityCV)) + geom_point()
ggplot(torp_cv_fit, aes(x = Fitted, y = LogQualityCV)) + geom_point()
# much weaker relationships with CV

# export models
save(hydr_mean_mod, file = "output/fwc_hydrilla_mean_phosphorus_model.rda")
save(flpl_mean_mod, file = "output/fwc_floating_plant_mean_phosphorus_model.rda")
save(cubu_mean_mod, file = "output/fwc_cuban_bulrush_mean_phosphorus_model.rda")
save(torp_mean_mod, file = "output/fwc_torpedograss_mean_phosphorus_model.rda")

save(hydr_max_mod, file = "output/fwc_hydrilla_max_phosphorus_model.rda")
save(flpl_max_mod, file = "output/fwc_floating_plant_max_phosphorus_model.rda")
save(cubu_max_mod, file = "output/fwc_cuban_bulrush_max_phosphorus_model.rda")
save(torp_max_mod, file = "output/fwc_torpedograss_max_phosphorus_model.rda")

save(hydr_min_mod, file = "output/fwc_hydrilla_min_phosphorus_model.rda")
save(flpl_min_mod, file = "output/fwc_floating_plant_min_phosphorus_model.rda")
save(cubu_min_mod, file = "output/fwc_cuban_bulrush_min_phosphorus_model.rda")
save(torp_min_mod, file = "output/fwc_torpedograss_min_phosphorus_model.rda")

save(hydr_cv_mod, file = "output/fwc_hydrilla_cv_phosphorus_model.rda")
save(flpl_cv_mod, file = "output/fwc_floating_plant_cv_phosphorus_model.rda")
save(cubu_cv_mod, file = "output/fwc_cuban_bulrush_cv_phosphorus_model.rda")
save(torp_cv_mod, file = "output/fwc_torpedograss_cv_phosphorus_model.rda")

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
foc_mod_se <- mod_se_fun(hydr_mean_mod, hydr_dat, "hydrilla") %>%
  mutate(Response = "mean") %>%
  full_join(mod_se_fun(hydr_max_mod, hydr_dat, "hydrilla") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(hydr_min_mod, hydr_dat, "hydrilla") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(hydr_cv_mod, hydr_dat, "hydrilla") %>%
              mutate(Response = "cv")) %>%
  full_join(mod_se_fun(flpl_mean_mod, flpl_dat, "floating plants") %>%
              mutate(Response = "mean")) %>%
  full_join(mod_se_fun(flpl_max_mod, flpl_dat, "floating plants") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(flpl_min_mod, flpl_dat, "floating plants") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(flpl_cv_mod, flpl_dat, "floating plants") %>%
              mutate(Response = "cv")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "PercCovered_c",
                           "management" = "RecentTreatment",
                           "invasive PAC:management" = "PercCovered_c:RecentTreatment")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Response)

non_foc_mod_se <- mod_se_fun(cubu_mean_mod, cubu_dat, "Cuban bulrush") %>%
  mutate(Response = "mean") %>%
  full_join(mod_se_fun(cubu_max_mod, cubu_dat, "Cuban bulrush") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(cubu_min_mod, cubu_dat, "Cuban bulrush") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(cubu_cv_mod, cubu_dat, "Cuban bulrush") %>%
              mutate(Response = "cv")) %>%
  full_join(mod_se_fun(torp_mean_mod, torp_dat, "torpedograss") %>%
              mutate(Response = "mean")) %>%
  full_join(mod_se_fun(torp_max_mod, torp_dat, "torpedograss") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(torp_min_mod, torp_dat, "torpedograss") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(torp_cv_mod, torp_dat, "torpedograss") %>%
              mutate(Response = "cv")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "PercCovered_c",
                           "management" = "RecentTreatment",
                           "invasive PAC:management" = "PercCovered_c:RecentTreatment")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Response)


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
hydr_mean_mod_full <- plm(LogQualityMean ~ Established, data = hydr_dat_full_fix, 
                          model = "random", random.method="walhus")
flpl_mean_mod_full <- update(hydr_mean_mod_full, data = flpl_dat_full)
cubu_mean_mod_full <- update(hydr_mean_mod_full, data = cubu_dat_full)
torp_mean_mod_full <- update(hydr_mean_mod_full, data = torp_dat_full)

hydr_max_mod_full <- plm(LogQualityMax ~ Established, data = hydr_dat_full_fix, 
                         model = "random", random.method="walhus")
flpl_max_mod_full <- update(hydr_max_mod_full, data = flpl_dat_full)
cubu_max_mod_full <- update(hydr_max_mod_full, data = cubu_dat_full)
torp_max_mod_full <- update(hydr_max_mod_full, data = torp_dat_full)

hydr_min_mod_full <- plm(LogQualityMin ~ Established, data = hydr_dat_full_fix, 
                         model = "random", random.method="walhus")
flpl_min_mod_full <- update(hydr_min_mod_full, data = flpl_dat_full)
cubu_min_mod_full <- update(hydr_min_mod_full, data = cubu_dat_full)
torp_min_mod_full <- update(hydr_min_mod_full, data = torp_dat_full)

hydr_cv_mod_full <- plm(LogQualityCV ~ Established, data = hydr_dat_full_fix, 
                        model = "random", random.method="walhus")
flpl_cv_mod_full <- update(hydr_cv_mod_full, data = flpl_dat_full)
cubu_cv_mod_full <- update(hydr_cv_mod_full, data = cubu_dat_full)
torp_cv_mod_full <- update(hydr_cv_mod_full, data = torp_dat_full)

# add fitted values to pdata.frame (important to match rows)
# convert to regular dataframe
hydr_mean_fit_full <- mutate(hydr_dat_full_fix, Fitted = as.numeric(hydr_dat_full_fix$LogQualityMean - hydr_mean_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_mean_fit_full <- mutate(flpl_dat_full_fix, Fitted = as.numeric(flpl_dat_full_fix$LogQualityMean - flpl_mean_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_mean_fit_full <- mutate(cubu_dat_full_fix, Fitted = as.numeric(cubu_dat_full_fix$LogQualityMean - cubu_mean_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_mean_fit_full <- mutate(torp_dat_full_fix, Fitted = as.numeric(torp_dat_full_fix$LogQualityMean - torp_mean_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_max_fit_full <- mutate(hydr_dat_full_fix, Fitted = as.numeric(hydr_dat_full_fix$LogQualityMax - hydr_max_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_max_fit_full <- mutate(flpl_dat_full_fix, Fitted = as.numeric(flpl_dat_full_fix$LogQualityMax - flpl_max_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_max_fit_full <- mutate(cubu_dat_full_fix, Fitted = as.numeric(cubu_dat_full_fix$LogQualityMax - cubu_max_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_max_fit_full <- mutate(torp_dat_full_fix, Fitted = as.numeric(torp_dat_full_fix$LogQualityMax - torp_max_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_min_fit_full <- mutate(hydr_dat_full_fix, Fitted = as.numeric(hydr_dat_full_fix$LogQualityMin - hydr_min_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_min_fit_full <- mutate(flpl_dat_full_fix, Fitted = as.numeric(flpl_dat_full_fix$LogQualityMin - flpl_min_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_min_fit_full <- mutate(cubu_dat_full_fix, Fitted = as.numeric(cubu_dat_full_fix$LogQualityMin - cubu_min_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_min_fit_full <- mutate(torp_dat_full_fix, Fitted = as.numeric(torp_dat_full_fix$LogQualityMin - torp_min_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_cv_fit_full <- mutate(hydr_dat_full_fix, Fitted = as.numeric(hydr_dat_full_fix$LogQualityCV - hydr_cv_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_cv_fit_full <- mutate(flpl_dat_full_fix, Fitted = as.numeric(flpl_dat_full_fix$LogQualityCV - flpl_cv_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_cv_fit_full <- mutate(cubu_dat_full_fix, Fitted = as.numeric(cubu_dat_full_fix$LogQualityCV - cubu_cv_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_cv_fit_full <- mutate(torp_dat_full_fix, Fitted = as.numeric(torp_dat_full_fix$LogQualityCV - torp_cv_mod_full$residuals)) %>%
  as.data.frame(keep.attributes = F)

# fitted vs. observed
ggplot(hydr_mean_fit_full, aes(x = Fitted, y = LogQualityMean)) + geom_point()
ggplot(flpl_mean_fit_full, aes(x = Fitted, y = LogQualityMean)) + geom_point()
ggplot(cubu_mean_fit_full, aes(x = Fitted, y = LogQualityMean)) + geom_point()
ggplot(torp_mean_fit_full, aes(x = Fitted, y = LogQualityMean)) + geom_point()

ggplot(hydr_max_fit_full, aes(x = Fitted, y = LogQualityMax)) + geom_point()
ggplot(flpl_max_fit_full, aes(x = Fitted, y = LogQualityMax)) + geom_point()
ggplot(cubu_max_fit_full, aes(x = Fitted, y = LogQualityMax)) + geom_point()
ggplot(torp_max_fit_full, aes(x = Fitted, y = LogQualityMax)) + geom_point()

ggplot(hydr_min_fit_full, aes(x = Fitted, y = LogQualityMin)) + geom_point()
ggplot(flpl_min_fit_full, aes(x = Fitted, y = LogQualityMin)) + geom_point()
ggplot(cubu_min_fit_full, aes(x = Fitted, y = LogQualityMin)) + geom_point()
ggplot(torp_min_fit_full, aes(x = Fitted, y = LogQualityMin)) + geom_point()

ggplot(hydr_cv_fit_full, aes(x = Fitted, y = LogQualityCV)) + geom_point()
ggplot(flpl_cv_fit_full, aes(x = Fitted, y = LogQualityCV)) + geom_point()
ggplot(cubu_cv_fit_full, aes(x = Fitted, y = LogQualityCV)) + geom_point()
ggplot(torp_cv_fit_full, aes(x = Fitted, y = LogQualityCV)) + geom_point()
# much weaker relationships with CV

# export models
save(hydr_mean_mod_full, file = "output/fwc_hydrilla_mean_phosphorus_binary_model.rda")
save(flpl_mean_mod_full, file = "output/fwc_floating_plant_mean_phosphorus_binary_model.rda")
save(cubu_mean_mod_full, file = "output/fwc_cuban_bulrush_mean_phosphorus_binary_model.rda")
save(torp_mean_mod_full, file = "output/fwc_torpedograss_mean_phosphorus_binary_model.rda")

save(hydr_max_mod_full, file = "output/fwc_hydrilla_max_phosphorus_binary_model.rda")
save(flpl_max_mod_full, file = "output/fwc_floating_plant_max_phosphorus_binary_model.rda")
save(cubu_max_mod_full, file = "output/fwc_cuban_bulrush_max_phosphorus_binary_model.rda")
save(torp_max_mod_full, file = "output/fwc_torpedograss_max_phosphorus_binary_model.rda")

save(hydr_min_mod_full, file = "output/fwc_hydrilla_min_phosphorus_binary_model.rda")
save(flpl_min_mod_full, file = "output/fwc_floating_plant_min_phosphorus_binary_model.rda")
save(cubu_min_mod_full, file = "output/fwc_cuban_bulrush_min_phosphorus_binary_model.rda")
save(torp_min_mod_full, file = "output/fwc_torpedograss_min_phosphorus_binary_model.rda")

save(hydr_cv_mod_full, file = "output/fwc_hydrilla_cv_phosphorus_binary_model.rda")
save(flpl_cv_mod_full, file = "output/fwc_floating_plant_cv_phosphorus_binary_model.rda")
save(cubu_cv_mod_full, file = "output/fwc_cuban_bulrush_cv_phosphorus_binary_model.rda")
save(torp_cv_mod_full, file = "output/fwc_torpedograss_cv_phosphorus_binary_model.rda")

# combine SE tables
foc_mod_full_se <- mod_se_fun(hydr_mean_mod_full, hydr_dat_full, "hydrilla") %>%
  mutate(Response = "mean") %>%
  full_join(mod_se_fun(hydr_max_mod_full, hydr_dat_full, "hydrilla") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(hydr_min_mod_full, hydr_dat_full, "hydrilla") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(hydr_cv_mod_full, hydr_dat_full, "hydrilla") %>%
              mutate(Response = "cv")) %>%
  full_join(mod_se_fun(flpl_mean_mod_full, flpl_dat_full, "floating plants") %>%
              mutate(Response = "mean")) %>%
  full_join(mod_se_fun(flpl_max_mod_full, flpl_dat_full, "floating plants") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(flpl_min_mod_full, flpl_dat_full, "floating plants") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(flpl_cv_mod_full, flpl_dat_full, "floating plants") %>%
              mutate(Response = "cv")) %>%
  mutate(term = fct_recode(term,
                           "established" = "Established",
                           "intercept" = "(Intercept)")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Response)

non_foc_mod_full_se <- mod_se_fun(cubu_mean_mod_full, cubu_dat_full, "Cuban bulrush") %>%
  mutate(Response = "mean") %>%
  full_join(mod_se_fun(cubu_max_mod_full, cubu_dat_full, "Cuban bulrush") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(cubu_min_mod_full, cubu_dat_full, "Cuban bulrush") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(cubu_cv_mod_full, cubu_dat_full, "Cuban bulrush") %>%
              mutate(Response = "cv")) %>%
  full_join(mod_se_fun(torp_mean_mod_full, torp_dat_full, "torpedograss", var_type = "HC2") %>%
              mutate(Response = "mean")) %>%
  full_join(mod_se_fun(torp_max_mod_full, torp_dat_full, "torpedograss", var_type = "HC2") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(torp_min_mod_full, torp_dat_full, "torpedograss", var_type = "HC0") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(torp_cv_mod_full, torp_dat_full, "torpedograss", var_type = "HC2") %>%
              mutate(Response = "cv")) %>%
  mutate(term = fct_recode(term,
                           "established" = "Established",
                           "intercept" = "(Intercept)")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Response)

# save
write_csv(foc_mod_full_se, "output/fwc_focal_phosphorus_binary_model_summary.csv")
write_csv(non_foc_mod_full_se, "output/fwc_non_focal_phosphorus_binary_model_summary.csv")


#### values for text ####

# extract intercept from fixed effects models
intercept_fun <- function(models, spp){
  
  dat_out <- tibble(Invasive = spp,
                    Response = c("mean", "max", "min", "cv"),
                    Intercept = sapply(models, function(x) mean(fixef(x)))) %>%
    mutate(Intercept = exp(Intercept))
  
  return(dat_out)
  
}

# combine model summary info
foc_mod_sum <- foc_mod_se %>%
  left_join(intercept_fun(list(hydr_mean_mod, hydr_max_mod, hydr_min_mod, hydr_cv_mod), "hydrilla") %>%
              full_join(intercept_fun(list(flpl_mean_mod, flpl_max_mod, flpl_min_mod, flpl_cv_mod), "floating plants"))) %>%
  mutate(Metric = "phosphorus") %>%
  relocate(Metric) %>%
  relocate(Intercept, .before = "R2")

write_csv(foc_mod_sum, "output/fwc_focal_phosphorus_model_summary.csv")

non_foc_mod_sum <- non_foc_mod_se %>%
  left_join(intercept_fun(list(cubu_mean_mod, cubu_max_mod, cubu_min_mod, cubu_cv_mod), "Cuban bulrush") %>%
              full_join(intercept_fun(list(torp_mean_mod, torp_max_mod, torp_min_mod, torp_cv_mod), "torpedograss"))) %>%
  mutate(Metric = "phosphorus") %>%
  relocate(Metric) %>%
  relocate(Intercept, .before = "R2")

write_csv(non_foc_mod_sum, "output/fwc_non_focal_phosphorus_model_summary.csv")

# combine binary and continuous models
foc_mean_summary <- foc_mod_full_se %>%
  select(Invasive, Response, Term, Estimate) %>%
  pivot_wider(names_from = Term,
              values_from = Estimate) %>%
  rename(Intercept = intercept,
         Estimate = established) %>%
  mutate(Term = "established",
         Model = "binary",
         Metric = "phosphorus",
         Intercept = exp(Intercept)) %>%
  left_join(foc_mod_full_se) %>%
  full_join(foc_mod_sum %>%
              mutate(Model = "continuous")) %>%
  filter(Response == "mean") %>%
  select(-Response) %>%
  relocate(Metric) %>%
  relocate(c(Model, Term), .after = "Invasive") %>%
  relocate(Intercept, .before = "R2") %>%
  arrange(Invasive, Model, Term)

non_foc_mean_summary <- non_foc_mod_full_se %>%
  select(Invasive, Response, Term, Estimate) %>%
  pivot_wider(names_from = Term,
              values_from = Estimate) %>%
  rename(Intercept = intercept,
         Estimate = established) %>%
  mutate(Term = "established",
         Model = "binary",
         Metric = "phosphorus",
         Intercept = exp(Intercept)) %>%
  left_join(non_foc_mod_full_se) %>%
  full_join(non_foc_mod_sum %>%
              mutate(Model = "continuous")) %>%
  filter(Response == "mean") %>%
  select(-Response) %>%
  relocate(Metric) %>%
  relocate(c(Model, Term), .after = "Invasive") %>%
  relocate(Intercept, .before = "R2") %>%
  arrange(Invasive, Model, Term)

# save
write_csv(foc_mean_summary, "output/fwc_focal_phosphorus_mean_models.csv")
write_csv(non_foc_mean_summary, "output/fwc_non_focal_phosphorus_mean_models.csv")
