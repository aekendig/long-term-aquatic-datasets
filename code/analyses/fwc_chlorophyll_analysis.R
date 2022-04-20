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
uninv <- read_csv("output/fwc_uninvaded_permID.csv") # lakes with no recorded invasion


#### edit data ####

# require longest reasonable lag
# longer lags mean fewer years are needed from chlorophyll dataset
# add max years column
inv_plant2 <- inv_plant %>%
  filter((CommonName != "Cuban bulrush" & !is.na(Lag6AvgPropCovered) & !is.na(Lag6Treated)) |
           (CommonName == "Cuban bulrush" & !is.na(Lag3AvgPropCovered) & !is.na(Lag3Treated))) %>%
  group_by(CommonName) %>%
  mutate(maxYears = n_distinct(GSYear)) %>%
  ungroup()

# check
inv_plant2 %>%
  group_by(CommonName) %>%
  summarize(maxYears = unique(maxYears),
            lastYear = max(GSYear),
            firstYear = min(GSYear),
            locations = n_distinct(PermanentID))

# year distribution
ggplot(lwwa_chl, aes(x = GSYear)) +
  geom_bar() +
  facet_wrap(~ Quarter)
# some are low for 2019

# combine chlorophyll, invasive, control
# select waterbodies sampled throughout
chl_dat <- lwwa_chl %>%
  # filter(!is.na(PrevValue)) %>%
  inner_join(inv_plant2 %>%
               filter(GSYear >= 2005 & GSYear < 2019) %>% # year cut-offs from data exploration below
               group_by(CommonName) %>%
               mutate(maxYears = n_distinct(GSYear)) %>% # recalculate max years
               ungroup()) %>%
  group_by(CommonName, PermanentID, Quarter) %>%
  mutate(nYears = n_distinct(GSYear)) %>% # years per waterbody
  ungroup() %>%
  filter(nYears == maxYears) %>%
  mutate(ValueDiff = QualityValue - PrevValue,  # change over time
         across(ends_with("AvgPropCovered"), ~ .x * 100),
         logQual = log(QualityValue)) %>%
  rename_with(str_replace, pattern = "AvgPropCovered", replacement = "AvgPercCovered")

# sample sizes
(chl_samp_sum <- chl_dat %>%
  group_by(CommonName, Quarter) %>%
  summarize(Years = n_distinct(GSYear),
            Waterbodies = n_distinct(PermanentID),
            minYear = min(GSYear),
            maxYear = max(GSYear)) %>%
  data.frame())
# torpedograss was completely missing before date cut-offs
# its distribution is highly concentrated in center of state
# https://nas.er.usgs.gov/viewer/omap.aspx?SpeciesID=1124
# para grass only has one or two waterbodies, but it only has 10 total

# why are there no data matching with torpedograss?
inv_plant2 %>%
  filter(CommonName == "Torpedograss") %>%
  select(PermanentID, GSYear) %>%
  unique() %>%
  inner_join(lwwa_chl %>%
               select(PermanentID, GSYear) %>%
               unique()) %>%
  group_by(GSYear) %>%
  count() # a lot of waterbodies each year

inv_plant2 %>%
  filter(CommonName == "Torpedograss") %>%
  select(PermanentID, GSYear) %>%
  unique() %>%
  inner_join(lwwa_chl %>%
               select(PermanentID, GSYear) %>%
               unique()) %>%
  group_by(PermanentID) %>%
  summarize(Years = n_distinct(GSYear),
            firstYear = min(GSYear),
            lastYear = max(GSYear),
            yearDiff = lastYear - firstYear + 1) %>%
  data.frame() 
# most are missing sampling in 2004

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

# remove para grass (only 1-2 waterbodies)
# use same waterbodies in all 4 quarters
chl_dat2 <- chl_dat %>%
  filter(CommonName != "Para grass") %>%
  group_by(CommonName, PermanentID, GSYear) %>%
  mutate(nQuart = n_distinct(Quarter)) %>%
  ungroup() %>%
  filter(nQuart == 4)

# sample sizes
chl_dat2 %>%
  group_by(CommonName, Quarter) %>%
  summarize(Years = n_distinct(GSYear),
            Waterbodies = n_distinct(PermanentID),
            N = Years * Waterbodies) %>%
  data.frame()

# split by species
hydr_dat <- filter(chl_dat2, CommonName == "Hydrilla")
wale_dat <- filter(chl_dat2, CommonName == "Water lettuce")
wahy_dat <- filter(chl_dat2, CommonName == "Water hyacinth")
torp_dat <- filter(chl_dat2, CommonName == "Torpedograss")
cubu_dat <- filter(chl_dat2, CommonName == "Cuban bulrush")

# add water quality to uninvaded dataset
# select years to match invasion dataset
uninv2 <- lwwa_chl %>%
  # filter(!is.na(PrevValue)) %>%
  inner_join(uninv) %>%
  left_join(chl_samp_sum %>%
              select(CommonName, minYear, maxYear) %>%
              unique()) %>%
  filter(GSYear >= minYear & GSYear <= maxYear) %>%
  mutate(logQual = log(QualityValue))


#### initial visualizations ####

# covariate correlations
chl_dat2 %>%
  select(Quarter, CommonName, 
         Lag3Treated, Lag3AvgPercCovered, PrevValue) %>%
  group_by(Quarter, CommonName) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & abs(corr) >= 0.4) %>%
  data.frame()
# water lettuce: prev value is correlated with PAC

# response distributions
ggplot(chl_dat2, aes(x = ValueDiff)) +
  geom_histogram() +
  facet_grid(CommonName ~ Quarter, scales = "free")
# normal and relatively wide around zero

ggplot(chl_dat2, aes(x = QualityValue)) +
  geom_histogram() +
  facet_grid(CommonName ~ Quarter, scales = "free")
# skewed

ggplot(chl_dat2, aes(x = logQual)) +
  geom_histogram() +
  facet_grid(CommonName ~ Quarter, scales = "free")

# coefficients and QualityValue
ggplot(chl_dat2, aes(x = Lag3AvgPercCovered, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# generally close to zero

ggplot(chl_dat2, aes(x = Lag3AvgPercCovered, y = logQual)) +
  geom_point(aes(color = PermanentID)) +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free") +
  theme(legend.position = "none")
# strong positive for water hyacinth and lettuce (driven by outliers)
# negative for hydrilla and torpedograss
# switches over time for Cuban bulrush

ggplot(chl_dat2, aes(x = Lag3AvgPercCovered, y = logQual, color = PermanentID)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_grid(Quarter ~ CommonName, scales = "free") +
  theme(legend.position = "none")
# all the high floating plant values are teh same lake, which has high chlorophyll
# hydrilla values are more spread out over lakes

ggplot(chl_dat2, aes(x = Lag3AvgPercCovered, y = logQual, color = as.factor(round(Lag3Treated, 1)))) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_grid(Quarter ~ CommonName, scales = "free") +
  scale_color_viridis_d(name = "Mgmt", direction = -1)

ggplot(chl_dat2, aes(x = Lag3Treated, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# negligible

ggplot(chl_dat2, aes(x = Lag3Treated, y = logQual)) +
  geom_point(aes(color = PermanentID)) +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free") +
  theme(legend.position = "none")
# shallow slopes
# negative for newer invasive spp
# slightly positive for hydrilla

ggplot(chl_dat2, aes(x = Lag3Treated, y = logQual, color = PermanentID)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_grid(Quarter ~ CommonName, scales = "free") +
  theme(legend.position = "none")

ggplot(chl_dat2, aes(x = PrevValue, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(CommonName ~ Quarter, scales = "free")
# consistently negative

ggplot(chl_dat2, aes(x = PrevValue, y = logQual)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(CommonName ~ Quarter, scales = "free")
# consistently positive


#### evaluate model structure ####

# function to fit models for each species
mod_structure_fits <- function(dat_in, quarter){
  
  # create fixed effects data frame
  # choose a quarter
  # cannot include quarter as an interaction in fixed-effect models because then there are duplicate
  # individual-time rows
  dat_fix <- dat_in %>%
    filter(Quarter == quarter) %>%
    ungroup() %>%
    pdata.frame(index = c("PermanentID", "GSYear"))
  # each waterbody is an individual
  
  # simple lm
  mod_lm <- lm(logQual ~ Lag3AvgPercCovered + Lag3Treated, data = dat_fix)
  
  # random effects
  mod_ran_loc <- glmmTMB(logQual ~ Lag3AvgPercCovered + Lag3Treated + (1|PermanentID), data = dat_fix)
  mod_ran_yr <- glmmTMB(logQual ~ Lag3AvgPercCovered + Lag3Treated + (1|GSYear), data = dat_fix)
  mod_ran_loc_yr <- glmmTMB(logQual ~ Lag3AvgPercCovered + Lag3Treated + (1|PermanentID) + (1|GSYear), data = dat_fix)
  
  # fixed effects
  mod_fix_loc <- plm(logQual ~ Lag3AvgPercCovered + Lag3Treated, data = dat_fix,
                      model = "within")
  mod_fix_loc_yr <- plm(logQual ~ Lag3AvgPercCovered + Lag3Treated, data = dat_fix,
                         model = "within", effect = "twoways")
  
  # use difference to account for reverse causality
  mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag3AvgPercCovered + Lag3Treated, data = dat_fix,
                             index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
  
  # return list of models
  return(list(lm = mod_lm,
              ran_loc = mod_ran_loc,
              ran_yr = mod_ran_yr,
              ran_loc_yr = mod_ran_loc_yr,
              fix_loc = mod_fix_loc,
              fix_loc_yr = mod_fix_loc_yr,
              diff_fix_loc_yr = mod_diff_fix_loc_yr))
  
}

# fit models for each species
hydr_mod_struc_1 <- mod_structure_fits(hydr_dat, 1)
wahy_mod_struc_1 <- mod_structure_fits(wahy_dat, 1)
wale_mod_struc_1 <- mod_structure_fits(wale_dat, 1)

hydr_mod_struc_2 <- mod_structure_fits(hydr_dat, 2)
wahy_mod_struc_2 <- mod_structure_fits(wahy_dat, 2)
wale_mod_struc_2 <- mod_structure_fits(wale_dat, 2)

hydr_mod_struc_3 <- mod_structure_fits(hydr_dat, 3)
wahy_mod_struc_3 <- mod_structure_fits(wahy_dat, 3)
wale_mod_struc_3 <- mod_structure_fits(wale_dat, 3)

hydr_mod_struc_4 <- mod_structure_fits(hydr_dat, 4)
wahy_mod_struc_4 <- mod_structure_fits(wahy_dat, 4)
wale_mod_struc_4 <- mod_structure_fits(wale_dat, 4)

# compare model estimates
hydr_mod_comp_1 <- mod_structure_comp(simp_mods = hydr_mod_struc_1[1], 
                                      ran_mods = hydr_mod_struc_1[2:4],
                                      fix_mods = hydr_mod_struc_1[5:7])
wahy_mod_comp_1 <- mod_structure_comp(simp_mods = wahy_mod_struc_1[1], 
                                      ran_mods = wahy_mod_struc_1[2:4],
                                      fix_mods = wahy_mod_struc_1[5:7]) 
wale_mod_comp_1 <- mod_structure_comp(simp_mods = wale_mod_struc_1[1], 
                                      ran_mods = wale_mod_struc_1[2:4],
                                      fix_mods = wale_mod_struc_1[5:7]) 

hydr_mod_comp_2 <- mod_structure_comp(simp_mods = hydr_mod_struc_2[1], 
                                      ran_mods = hydr_mod_struc_2[2:4],
                                      fix_mods = hydr_mod_struc_2[5:7])
wahy_mod_comp_2 <- mod_structure_comp(simp_mods = wahy_mod_struc_2[1], 
                                      ran_mods = wahy_mod_struc_2[2:4],
                                      fix_mods = wahy_mod_struc_2[5:7]) 
wale_mod_comp_2 <- mod_structure_comp(simp_mods = wale_mod_struc_2[1], 
                                      ran_mods = wale_mod_struc_2[2:4],
                                      fix_mods = wale_mod_struc_2[5:7]) 

hydr_mod_comp_3 <- mod_structure_comp(simp_mods = hydr_mod_struc_3[1], 
                                      ran_mods = hydr_mod_struc_3[2:4],
                                      fix_mods = hydr_mod_struc_3[5:7])
wahy_mod_comp_3 <- mod_structure_comp(simp_mods = wahy_mod_struc_3[1], 
                                      ran_mods = wahy_mod_struc_3[2:4],
                                      fix_mods = wahy_mod_struc_3[5:7]) 
wale_mod_comp_3 <- mod_structure_comp(simp_mods = wale_mod_struc_3[1], 
                                      ran_mods = wale_mod_struc_3[2:4],
                                      fix_mods = wale_mod_struc_3[5:7]) 

hydr_mod_comp_4 <- mod_structure_comp(simp_mods = hydr_mod_struc_4[1], 
                                      ran_mods = hydr_mod_struc_4[2:4],
                                      fix_mods = hydr_mod_struc_4[5:7])
wahy_mod_comp_4 <- mod_structure_comp(simp_mods = wahy_mod_struc_4[1], 
                                      ran_mods = wahy_mod_struc_4[2:4],
                                      fix_mods = wahy_mod_struc_4[5:7]) 
wale_mod_comp_4 <- mod_structure_comp(simp_mods = wale_mod_struc_4[1], 
                                      ran_mods = wale_mod_struc_4[2:4],
                                      fix_mods = wale_mod_struc_4[5:7]) 

# combine species
mod_comp <- hydr_mod_comp_1 %>%
  mutate(Quarter = 1) %>%
  full_join(hydr_mod_comp_2 %>%
              mutate(Quarter = 2)) %>%
  full_join(hydr_mod_comp_3 %>%
              mutate(Quarter = 3)) %>%
  full_join(hydr_mod_comp_4 %>%
              mutate(Quarter = 4)) %>%
  mutate(Species = "hydrilla") %>%
  full_join(wahy_mod_comp_1 %>%
              mutate(Quarter = 1) %>%
              full_join(wahy_mod_comp_2 %>%
                          mutate(Quarter = 2)) %>%
              full_join(wahy_mod_comp_3 %>%
                          mutate(Quarter = 3)) %>%
              full_join(wahy_mod_comp_4 %>%
                          mutate(Quarter = 4)) %>%
              mutate(Species = "water hyacinth")) %>%
  full_join(wale_mod_comp_1 %>%
              mutate(Quarter = 1) %>%
              full_join(wale_mod_comp_2 %>%
                          mutate(Quarter = 2)) %>%
              full_join(wale_mod_comp_3 %>%
                          mutate(Quarter = 3)) %>%
              full_join(wale_mod_comp_4 %>%
                          mutate(Quarter = 4)) %>%
              mutate(Species = "water lettuce")) %>%
  mutate(coefficients = str_replace(coefficients, "Lag3Treated", "management"),
         coefficients = str_replace(coefficients, "Lag3AvgPercCovered", "PAC"),
         across(!c(coefficients, Species), ~ round(.x, digits = 3))) %>%
  relocate(Species, Quarter)

write_csv(mod_comp, "output/fwc_chlorophyll_model_structure_comparison.csv")

# model comparison notes:
# Hydrilla PAC negatively affects chlorophyll until location-specific intercepts
# are included, and then the estimate becomes ~ 0, indicating that lakes with
# more hydrilla have lower chlorophyll, but the same isn't true over time
# within a lake. The difference model suggests hydrilla PAC accelerates chlorophyll
# change. Hydrilla management has an inverse pattern: more management increases
# chlorophyll unless there are location-specific intercepts, and then it decreases
# chlorophyll, indicating lakes with more management have higher chlorophyll, but
# more management within a lake decreases chlorophyll (opposite of expected).
# Difference model estimates large positive hydrilla management effect.

# Water hyacinth and water lettuce PAC increase chlorophyll. The magnitude of the
# effect decreases with location-specific intercepts, indicating lakes with more 
# floating plants have higher chlorophyll, but temporal variation in floating
# plant coverage is less influential (and also less variable). Water hyacinth
# management increases chlorophyll with minimal effects of location/year intercepts.
# water lettuce management decreases chlorophyll, which becomes a positive effect with
# location-specific intercepts, indicating lakes with more management have lower
# chlorophyll, but management within a lake increases chlorophyll

# test fixed effects (seems like year isn't necessary)
# have to refit because data need to be accessible (not "dat_fix")
hydr_mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag3AvgPercCovered + Lag3Treated, 
                                data = filter(hydr_dat, Quarter == 1), 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(hydr_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(hydr_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

wahy_mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag3AvgPercCovered + Lag3Treated, 
                                data = filter(wahy_dat, Quarter == 1), 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(wahy_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(wahy_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

wale_mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag3AvgPercCovered + Lag3Treated, 
                                data = filter(wale_dat, Quarter == 1), 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(wale_mod_diff_fix_loc_yr, effect = "time", type = "bp") # marginal
plmtest(wale_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

# use log quality without initial value and with waterbody and year fixed effects


#### model-fitting functions ####

# data filter function
dat_mod_filt <- function(treat_col, inv_col, dat_in){
  
  dat_mod <- dat_in %>%
    mutate(Treated = !!sym(treat_col),
           AvgPercCovered = !!sym(inv_col))
  
  return(dat_mod)
  
}

# function to fit models
mod_fit <- function(dat_in){
  
  # focal species
  foc_sp <- unique(dat_in$CommonName)
  
  # subset data
  dat_mod1 <- dat_mod_filt("Lag1Treated", "Lag1AvgPercCovered", dat_in)
  dat_mod2 <- dat_mod_filt("Lag2Treated", "Lag2AvgPercCovered", dat_in)
  dat_mod3 <- dat_mod_filt("Lag3Treated", "Lag3AvgPercCovered", dat_in)
  dat_mod4 <- dat_mod_filt("Lag4Treated", "Lag4AvgPercCovered", dat_in)
  dat_mod5 <- dat_mod_filt("Lag5Treated", "Lag5AvgPercCovered", dat_in)
  dat_mod6 <- dat_mod_filt("Lag6Treated", "Lag6AvgPercCovered", dat_in)
  
  # fit models
  if(foc_sp == "Cuban bulrush") {
    
    mod1 <- plm(logQual ~ AvgPercCovered + Treated, data = dat_mod1,
                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
    mod2 <- update(mod1, data = dat_mod2)
    mod3 <- update(mod1, data = dat_mod3)
    
    mods_out <- list(mod1, mod2, mod3)
    
  } else {
    
    mod1 <- plm(logQual ~ AvgPercCovered + Treated, data = dat_mod1,
                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
    mod2 <- update(mod1, data = dat_mod2)
    mod3 <- update(mod1, data = dat_mod3)
    mod4 <- update(mod1, data = dat_mod4)
    mod5 <- update(mod1, data = dat_mod5)
    mod6 <- update(mod1, data = dat_mod6)
    
    mods_out <- list(mod1, mod2, mod3, mod4, mod5, mod6)
    
  }
  
  # output
  return(mods_out)
  
}


#### fit models ####

# fit models with all lags
hydr_mods_q1 <- mod_fit(filter(hydr_dat, Quarter == 1))
wahy_mods_q1 <- mod_fit(filter(wahy_dat, Quarter == 1))
wale_mods_q1 <- mod_fit(filter(wale_dat, Quarter == 1))
torp_mods_q1 <- mod_fit(filter(torp_dat, Quarter == 1))
cubu_mods_q1 <- mod_fit(filter(cubu_dat, Quarter == 1))

hydr_mods_q2 <- mod_fit(filter(hydr_dat, Quarter == 2))
wahy_mods_q2 <- mod_fit(filter(wahy_dat, Quarter == 2))
wale_mods_q2 <- mod_fit(filter(wale_dat, Quarter == 2))
torp_mods_q2 <- mod_fit(filter(torp_dat, Quarter == 2))
cubu_mods_q2 <- mod_fit(filter(cubu_dat, Quarter == 2))

hydr_mods_q3 <- mod_fit(filter(hydr_dat, Quarter == 3))
wahy_mods_q3 <- mod_fit(filter(wahy_dat, Quarter == 3))
wale_mods_q3 <- mod_fit(filter(wale_dat, Quarter == 3))
torp_mods_q3 <- mod_fit(filter(torp_dat, Quarter == 3))
cubu_mods_q3 <- mod_fit(filter(cubu_dat, Quarter == 3))

hydr_mods_q4 <- mod_fit(filter(hydr_dat, Quarter == 4))
wahy_mods_q4 <- mod_fit(filter(wahy_dat, Quarter == 4))
wale_mods_q4 <- mod_fit(filter(wale_dat, Quarter == 4))
torp_mods_q4 <- mod_fit(filter(torp_dat, Quarter == 4))
cubu_mods_q4 <- mod_fit(filter(cubu_dat, Quarter == 4))

# name models
names(hydr_mods_q1) <- names(wahy_mods_q1) <- names(wale_mods_q1) <- names(torp_mods_q1) <- names(hydr_mods_q2) <- names(wahy_mods_q2) <- names(wale_mods_q2) <- names(torp_mods_q2) <- names(hydr_mods_q3) <- names(wahy_mods_q3) <- names(wale_mods_q3) <- names(torp_mods_q3) <- names(hydr_mods_q4) <- names(wahy_mods_q4) <- names(wale_mods_q4) <- names(torp_mods_q4) <- c("1", "2", "3", "4", "5", "6")
  
names(cubu_mods_q1) <- names(cubu_mods_q2) <- names(cubu_mods_q3) <- names(cubu_mods_q4) <- c("1", "2", "3")


#### coefficient figures and tables ####

# rename coefficients
coef_names <- c("Treated" = "Management", 
                "AvgPercCovered" = "Invasive PAC")

# ggplot function
plot_fun <- function(models){
  
  plot_out <- modelplot(models,
                        coef_map = coef_names,
                        background = list(geom_vline(xintercept = 0, color = "black",
                                                     size = 0.5, linetype = "dashed"))) +
    scale_color_viridis_d(direction = -1) +
    scale_x_continuous(labels = scale_fun_1) +
    def_theme_paper +
    theme(plot.title = element_text(size = 9))
  
  return(plot_out)
  
}

# panel plot function
panel_plot_fun <- function(mods1, mods2, mods3, 
                           spp1, spp2, spp3,
                           filename){
  
  # focal panels
  fig1 <- plot_fun(mods1) +
    labs(x = "",
         title = paste("(A)", spp1)) +
    def_theme_paper +
    theme(legend.position = "none")
  
  fig2 <- plot_fun(mods2) +
    labs(x = expression(paste("Estimate"%+-%" 95% CI", sep = "")),
         title = paste("(B)", spp2)) +
    theme(legend.position = "none",
          axis.text.y = element_blank())
  
  fig3 <- plot_fun(mods3) +
    labs(x = "",
         title = paste("(C)", spp3)) +
    theme(axis.text.y = element_blank(),
          legend.box.margin = margin(-10, 0, -10, -10)) +
    scale_color_viridis_d(direction = -1, name = "Lag\n(years)") +
    guides(color = guide_legend(reverse = TRUE))
  
  comb_fig <- fig1 + fig2 + fig3 + plot_annotation(
    theme = theme(plot.margin = margin(5, -5, 0, -10),
                  plot.title = element_text(size = 10, hjust = 0.5)),
    title = "Effects on chlorophyll a")
  
  ggsave(filename, comb_fig,
         device = "eps", width = 6.5, height = 3, units = "in")
  
}

panel_plot_non_foc_fun <- function(mods1, mods2,
                                   spp1, spp2,
                                   filename){
  
  # focal panels
  fig1 <- plot_fun(mods1) +
    labs(x = "",
         title = paste("(A)", spp1)) +
    theme(legend.position = "none")
  
  fig2 <- plot_fun(mods2) +
    labs(x = "",
         title = paste("(B)", spp2)) +
    theme(axis.text.y = element_blank(),
          legend.box.margin = margin(-10, 0, -10, -10)) +
    scale_color_viridis_d(direction = -1, name = "Lag\n(years)") +
    guides(color = guide_legend(reverse = TRUE))
  
  comb_fig <- (fig1 + fig2) + 
    # labs(x = expression(paste("Estimate"%+-%" 95% CI", sep = ""))) +
    # xlab(label = "Estimate2") +
    plot_annotation(
      caption = expression(paste("Estimate"%+-%" 95% CI", sep = "")),
      theme = theme(plot.caption = element_text(size = 9, color="black", hjust = 0.6, vjust = 10),
                    plot.margin = margin(5, -5, -5, -10),
                    plot.title = element_text(size = 10, hjust = 0.5)),
      title = "Effects on chlorophyll a")
  
  ggsave(filename, comb_fig,
         device = "eps", width = 4.7, height = 3, units = "in")
  
}

# figures
panel_plot_fun(hydr_mods_q1, wahy_mods_q1, wale_mods_q1,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_chlorophyll_quarter1_diff_model.eps")
panel_plot_fun(hydr_mods_q2, wahy_mods_q2, wale_mods_q2,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_chlorophyll_quarter2_diff_model.eps")
panel_plot_fun(hydr_mods_q3, wahy_mods_q3, wale_mods_q3,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_chlorophyll_quarter3_diff_model.eps")
panel_plot_fun(hydr_mods_q4, wahy_mods_q4, wale_mods_q4,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_chlorophyll_quarter4_diff_model.eps")

panel_plot_non_foc_fun(cubu_mods_q1, torp_mods_q1,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_chlorophyll_quarter1_diff_model.eps")

panel_plot_non_foc_fun(cubu_mods_q2, torp_mods_q2,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_chlorophyll_quarter2_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q3, torp_mods_q3,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_chlorophyll_quarter3_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q4, torp_mods_q4,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_chlorophyll_quarter4_diff_model.eps")

# lag/quarter notes
# floating plants and management have stronger positive effects on chlorophyll with greater lag times
# hydrilla management has weaker negative effect on chlorophyll with greater lag times
# torpedograss has stronger positive effect on chlorophyll with greater lag times


#### finalize models ####

# data format function
dat_mod_fin <- function(dat_in, quarter){
  
  # center invasive PAC
  # make p data frame
  dat_out <- dat_in %>%
    filter(Quarter == quarter) %>%
    mutate(AvgPercCovered = Lag3AvgPercCovered,
           Treated = Lag3Treated) %>%
    pdata.frame(index = c("PermanentID", "GSYear"))
  
  return(dat_out)
    
}

# format data
hydr_dat3_q1 <- dat_mod_fin(hydr_dat, 1)
wahy_dat3_q1 <- dat_mod_fin(wahy_dat, 1)
wale_dat3_q1 <- dat_mod_fin(wale_dat, 1)
cubu_dat3_q1 <- dat_mod_fin(cubu_dat, 1)
torp_dat3_q1 <- dat_mod_fin(torp_dat, 1)

hydr_dat3_q2 <- dat_mod_fin(hydr_dat, 2)
wahy_dat3_q2 <- dat_mod_fin(wahy_dat, 2)
wale_dat3_q2 <- dat_mod_fin(wale_dat, 2)
cubu_dat3_q2 <- dat_mod_fin(cubu_dat, 2)
torp_dat3_q2 <- dat_mod_fin(torp_dat, 2)

hydr_dat3_q3 <- dat_mod_fin(hydr_dat, 3)
wahy_dat3_q3 <- dat_mod_fin(wahy_dat, 3)
wale_dat3_q3 <- dat_mod_fin(wale_dat, 3)
cubu_dat3_q3 <- dat_mod_fin(cubu_dat, 3)
torp_dat3_q3 <- dat_mod_fin(torp_dat, 3)

hydr_dat3_q4 <- dat_mod_fin(hydr_dat, 4)
wahy_dat3_q4 <- dat_mod_fin(wahy_dat, 4)
wale_dat3_q4 <- dat_mod_fin(wale_dat, 4)
cubu_dat3_q4 <- dat_mod_fin(cubu_dat, 4)
torp_dat3_q4 <- dat_mod_fin(torp_dat, 4)

# fit models
hydr_chl_mod_q1 <- plm(logQual ~ AvgPercCovered + Treated, data = hydr_dat3_q1,
                       index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
wahy_chl_mod_q1 <- update(hydr_chl_mod_q1, data = wahy_dat3_q1)
wale_chl_mod_q1 <- update(hydr_chl_mod_q1, data = wale_dat3_q1)
cubu_chl_mod_q1 <- update(hydr_chl_mod_q1, data = cubu_dat3_q1)
torp_chl_mod_q1 <- update(hydr_chl_mod_q1, data = torp_dat3_q1)

hydr_chl_mod_q2 <- update(hydr_chl_mod_q1, data = hydr_dat3_q2)
wahy_chl_mod_q2 <- update(hydr_chl_mod_q2, data = wahy_dat3_q2)
wale_chl_mod_q2 <- update(hydr_chl_mod_q2, data = wale_dat3_q2)
cubu_chl_mod_q2 <- update(hydr_chl_mod_q2, data = cubu_dat3_q2)
torp_chl_mod_q2 <- update(hydr_chl_mod_q2, data = torp_dat3_q2)

hydr_chl_mod_q3 <- update(hydr_chl_mod_q1, data = hydr_dat3_q3)
wahy_chl_mod_q3 <- update(hydr_chl_mod_q3, data = wahy_dat3_q3)
wale_chl_mod_q3 <- update(hydr_chl_mod_q3, data = wale_dat3_q3)
cubu_chl_mod_q3 <- update(hydr_chl_mod_q3, data = cubu_dat3_q3)
torp_chl_mod_q3 <- update(hydr_chl_mod_q3, data = torp_dat3_q3)

hydr_chl_mod_q4 <- update(hydr_chl_mod_q1, data = hydr_dat3_q4)
wahy_chl_mod_q4 <- update(hydr_chl_mod_q4, data = wahy_dat3_q4)
wale_chl_mod_q4 <- update(hydr_chl_mod_q4, data = wale_dat3_q4)
cubu_chl_mod_q4 <- update(hydr_chl_mod_q4, data = cubu_dat3_q4)
torp_chl_mod_q4 <- update(hydr_chl_mod_q4, data = torp_dat3_q4)

# SE with heteroscedasticity and autocorrelation
coeftest(hydr_chl_mod_q1, vcov = vcovHC, type = "HC3") # -treat
coeftest(wahy_chl_mod_q1, vcov = vcovHC, type = "HC3")
coeftest(wale_chl_mod_q1, vcov = vcovHC, type = "HC3") # +treat
coeftest(cubu_chl_mod_q1, vcov = vcovHC, type = "HC3") 
coeftest(torp_chl_mod_q1, vcov = vcovHC, type = "HC3") 

coeftest(hydr_chl_mod_q2, vcov = vcovHC, type = "HC3")
coeftest(wahy_chl_mod_q2, vcov = vcovHC, type = "HC3") # +PAC
coeftest(wale_chl_mod_q2, vcov = vcovHC, type = "HC3") # +PAC
coeftest(cubu_chl_mod_q2, vcov = vcovHC, type = "HC3")
coeftest(torp_chl_mod_q2, vcov = vcovHC, type = "HC3") 

coeftest(hydr_chl_mod_q3, vcov = vcovHC, type = "HC3") # -treat
coeftest(wahy_chl_mod_q3, vcov = vcovHC, type = "HC3") 
coeftest(wale_chl_mod_q3, vcov = vcovHC, type = "HC3") 
coeftest(cubu_chl_mod_q3, vcov = vcovHC, type = "HC3") 
coeftest(torp_chl_mod_q3, vcov = vcovHC, type = "HC3") # +PAC

coeftest(hydr_chl_mod_q4, vcov = vcovHC, type = "HC3")
coeftest(wahy_chl_mod_q4, vcov = vcovHC, type = "HC3") # +PAC 
coeftest(wale_chl_mod_q4, vcov = vcovHC, type = "HC3") # +PAC
coeftest(cubu_chl_mod_q4, vcov = vcovHC, type = "HC3") 
coeftest(torp_chl_mod_q4, vcov = vcovHC, type = "HC3") 

# add fitted values to pdata.frame (important to match rows)
# convert to regular dataframe
hydr_fit_q1 <- mutate(hydr_dat3_q1, Fitted = as.numeric(hydr_dat3_q1$logQual - hydr_chl_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q1 <- mutate(wahy_dat3_q1, Fitted = as.numeric(wahy_dat3_q1$logQual - wahy_chl_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q1 <- mutate(wale_dat3_q1, Fitted = as.numeric(wale_dat3_q1$logQual - wale_chl_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q1 <- mutate(cubu_dat3_q1, Fitted = as.numeric(cubu_dat3_q1$logQual - cubu_chl_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q1 <- mutate(torp_dat3_q1, Fitted = as.numeric(torp_dat3_q1$logQual - torp_chl_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_fit_q2 <- mutate(hydr_dat3_q2, Fitted = as.numeric(hydr_dat3_q2$logQual - hydr_chl_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q2 <- mutate(wahy_dat3_q2, Fitted = as.numeric(wahy_dat3_q2$logQual - wahy_chl_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q2 <- mutate(wale_dat3_q2, Fitted = as.numeric(wale_dat3_q2$logQual - wale_chl_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q2 <- mutate(cubu_dat3_q2, Fitted = as.numeric(cubu_dat3_q2$logQual - cubu_chl_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q2 <- mutate(torp_dat3_q2, Fitted = as.numeric(torp_dat3_q2$logQual - torp_chl_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_fit_q3 <- mutate(hydr_dat3_q3, Fitted = as.numeric(hydr_dat3_q3$logQual - hydr_chl_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q3 <- mutate(wahy_dat3_q3, Fitted = as.numeric(wahy_dat3_q3$logQual - wahy_chl_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q3 <- mutate(wale_dat3_q3, Fitted = as.numeric(wale_dat3_q3$logQual - wale_chl_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q3 <- mutate(cubu_dat3_q3, Fitted = as.numeric(cubu_dat3_q3$logQual - cubu_chl_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q3 <- mutate(torp_dat3_q3, Fitted = as.numeric(torp_dat3_q3$logQual - torp_chl_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_fit_q4 <- mutate(hydr_dat3_q4, Fitted = as.numeric(hydr_dat3_q4$logQual - hydr_chl_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q4 <- mutate(wahy_dat3_q4, Fitted = as.numeric(wahy_dat3_q4$logQual - wahy_chl_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q4 <- mutate(wale_dat3_q4, Fitted = as.numeric(wale_dat3_q4$logQual - wale_chl_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q4 <- mutate(cubu_dat3_q4, Fitted = as.numeric(cubu_dat3_q4$logQual - cubu_chl_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q4 <- mutate(torp_dat3_q4, Fitted = as.numeric(torp_dat3_q4$logQual - torp_chl_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)

# fitted vs. observed
ggplot(hydr_fit_q1, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(wahy_fit_q1, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(wale_fit_q1, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(cubu_fit_q1, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(torp_fit_q1, aes(x = Fitted, y = logQual)) + geom_point()

ggplot(hydr_fit_q2, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(wahy_fit_q2, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(wale_fit_q2, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(cubu_fit_q2, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(torp_fit_q2, aes(x = Fitted, y = logQual)) + geom_point()

ggplot(hydr_fit_q3, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(wahy_fit_q3, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(wale_fit_q3, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(cubu_fit_q3, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(torp_fit_q3, aes(x = Fitted, y = logQual)) + geom_point()

ggplot(hydr_fit_q4, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(wahy_fit_q4, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(wale_fit_q4, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(cubu_fit_q4, aes(x = Fitted, y = logQual)) + geom_point()
ggplot(torp_fit_q4, aes(x = Fitted, y = logQual)) + geom_point()
# pretty good model fits

# combine models
hydr_chl_mods <- list(hydr_chl_mod_q1, hydr_chl_mod_q2, hydr_chl_mod_q3, hydr_chl_mod_q4)
wahy_chl_mods <- list(wahy_chl_mod_q1, wahy_chl_mod_q2, wahy_chl_mod_q3, wahy_chl_mod_q4)
wale_chl_mods <- list(wale_chl_mod_q1, wale_chl_mod_q2, wale_chl_mod_q3, wale_chl_mod_q4)
cubu_chl_mods <- list(cubu_chl_mod_q1, cubu_chl_mod_q2, cubu_chl_mod_q3, cubu_chl_mod_q4)
torp_chl_mods <- list(torp_chl_mod_q1, torp_chl_mod_q2, torp_chl_mod_q3, torp_chl_mod_q4)

# export models
save(hydr_chl_mods, file = "output/fwc_hydrilla_chlorophyll_models.rda")
save(wahy_chl_mods, file = "output/fwc_water_hyacinth_chlorophyll_models.rda")
save(wale_chl_mods, file = "output/fwc_water_lettuce_chlorophyll_models.rda")
save(cubu_chl_mods, file = "output/fwc_cuban_bulrush_chlorophyll_models.rda")
save(torp_chl_mods, file = "output/fwc_torpedograss_chlorophyll_models.rda")

# load models
load("output/fwc_hydrilla_chlorophyll_models.rda")
load("output/fwc_water_hyacinth_chlorophyll_models.rda")
load("output/fwc_water_lettuce_chlorophyll_models.rda")
load("output/fwc_cuban_bulrush_chlorophyll_models.rda")
load("output/fwc_torpedograss_chlorophyll_models.rda")

# process model SE
mod_se_fun <- function(models, dat, spp){
  
  dat_out <- tidy(coeftest(models[[1]], vcov = vcovHC, type = "HC3")) %>%
    mutate(Quarter = "Apr-Jun",
           R2 = r.squared(models[[1]])) %>%
    full_join(tidy(coeftest(models[[2]], vcov = vcovHC, type = "HC3")) %>%
                mutate(Quarter = "Jul-Sep",
                       R2 = r.squared(models[[2]]))) %>%
    full_join(tidy(coeftest(models[[3]], vcov = vcovHC, type = "HC3")) %>%
                mutate(Quarter = "Oct-Dec",
                       R2 = r.squared(models[[3]]))) %>%
    full_join(tidy(coeftest(models[[4]], vcov = vcovHC, type = "HC3")) %>%
                mutate(Quarter = "Jan-Mar",
                       R2 = r.squared(models[[4]]))) %>%
    mutate(Invasive = spp,
           Waterbodies = n_distinct(dat$PermanentID),
           Years = n_distinct(dat$GSYear),
           N = Waterbodies * Years)
  
  return(dat_out)
  
}

# combine SE tables
foc_mod_se <- mod_se_fun(hydr_chl_mods, hydr_dat, "hydrilla") %>%
  full_join(mod_se_fun(wahy_chl_mods, wahy_dat, "water hyacinth")) %>%
  full_join(mod_se_fun(wale_chl_mods, wale_dat, "water lettuce")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "AvgPercCovered",
                           "management" = "Treated")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Quarter)

non_foc_mod_se <- mod_se_fun(cubu_chl_mods, cubu_dat, "Cuban bulrush") %>%
  full_join(mod_se_fun(torp_chl_mods, torp_dat, "torpedograss")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "AvgPercCovered",
                           "management" = "Treated")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Quarter)

# export
write_csv(foc_mod_se, "output/fwc_focal_chlorophyll_model_summary.csv")
write_csv(non_foc_mod_se, "output/fwc_non_focal_chlorophyll_model_summary.csv")


#### values for text ####

# summarize uninvaded
uninv_sum <- uninv2 %>%
  group_by(CommonName, Quarter) %>%
  summarize(UninvAvg = mean(QualityValue),
            UninvN = n()) %>%
  ungroup()

# translate model coefficients
mod_coef_fun <- function(models, spp){
  
  dat_out <- tibble(Invasive = spp,
                    Quarter = c("Apr-Jun", "Jul-Sep", "Oct-Dec", "Jan-Mar"),
                    Intercept = sapply(models, function(x) mean(fixef(x))),
                    BetaPAC = sapply(models, function(x) coef(x)[1]),
                    BetaTreat = sapply(models, function(x) coef(x)[2]))
  
  return(dat_out)
  
}

# identify significant effects
foc_sig <- foc_mod_se %>%
  filter(P < 0.1) %>%
  select(Invasive, Quarter, Term) %>%
  left_join(mod_coef_fun(hydr_chl_mods, "hydrilla") %>%
              full_join(mod_coef_fun(wahy_chl_mods, "water hyacinth")) %>%
              full_join(mod_coef_fun(wale_chl_mods, "water lettuce"))) %>%
  mutate(BetaPAC = if_else(Term == "management", NA_real_, BetaPAC),
         BetaTreat = if_else(Term == "invasive PAC", NA_real_, BetaTreat),
         Metric = "chlorophyll a") %>%
  left_join(uninv_sum %>%
              mutate(Quarter = case_when(Quarter == 1 ~ "Apr-Jun", 
                                         Quarter == 2 ~ "Jul-Sep", 
                                         Quarter == 3 ~ "Oct-Dec", 
                                         Quarter == 4 ~ "Jan-Mar"),
                     CommonName = tolower(CommonName)) %>%
              rename(Invasive = CommonName))

write_csv(foc_sig, "output/fwc_focal_invasive_chlorophyll_significant.csv")

non_foc_sig <- non_foc_mod_se %>%
  filter(P < 0.1) %>%
  select(Invasive, Quarter, Term) %>%
  left_join(mod_coef_fun(cubu_chl_mods, "Cuban bulrush") %>%
              full_join(mod_coef_fun(torp_chl_mods, "torpedograss"))) %>%
  mutate(BetaPAC = if_else(Term == "management", NA_real_, BetaPAC),
         BetaTreat = if_else(Term == "invasive PAC", NA_real_, BetaTreat),
         Metric = "chlorophyll a") %>%
  left_join(uninv_sum %>%
              mutate(Quarter = case_when(Quarter == 1 ~ "Apr-Jun", 
                                         Quarter == 2 ~ "Jul-Sep", 
                                         Quarter == 3 ~ "Oct-Dec", 
                                         Quarter == 4 ~ "Jan-Mar"),
                     CommonName = fct_recode(CommonName, 
                                             "torpedograss" = "Torpedograss")) %>%
              rename(Invasive = CommonName))

write_csv(non_foc_sig, "output/fwc_non_focal_invasive_chlorophyll_significant.csv")
