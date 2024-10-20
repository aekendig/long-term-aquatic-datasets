#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom.mixed)
library(glmmTMB)
library(DHARMa)
library(GGally)
library(betareg)
library(janitor)

# figure settings
source("code/settings/figure_settings.R")

# import data
target_dat <- read_csv("intermediate-data/FWC_plant_management_target_analysis_formatted.csv")
target_taxa_dat <- read_csv("intermediate-data/FWC_plant_management_target_taxa_analysis_formatted.csv")
methods_dat <- read_csv("intermediate-data/FWC_plant_management_methods_analysis_formatted.csv")
methods_taxa_dat <- read_csv("intermediate-data/FWC_plant_management_methods_taxa_analysis_formatted.csv")


#### examine full dataset ####
target_dat_var <- target_dat %>%
  select(AreaOfInterestID, HydrPAC, FloatPAC, HydrTrtFreq, HydrTrtArea, FloatTrtFreq, FloatTrtArea,
           OtherTrtFreq, OtherTrtArea, ends_with("F")) %>%
  distinct() %>%
  mutate(HydrTrt = if_else(HydrTrtFreq == 0, 0, 1),
         FloatTrt = if_else(FloatTrtFreq == 0, 0, 1),
         across(.cols = ends_with("F"), 
                .fns = ~ fct_relevel(.x, "none", "low")))

# correlations among explanatory variables
target_dat_var %>%
  select(-c(AreaOfInterestID, HydrTrt, FloatTrt, ends_with("F"))) %>%
  ggpairs()
# correlations >= 0.4 and sig
# OtherTrtFreq and FloatTrtFreq: 0.6
# OtherTrtArea and FloatTrtArea: 0.4
# neither looks like a particularly strong relationship

# order factors
target_dat2 <- target_dat %>%
  mutate(across(.cols = ends_with("F"), 
                .fns = ~ fct_relevel(.x, "none", "low")))


#### examine newer data ####

# correlations among explanatory variables
methods_dat %>%
  select(c(AreaOfInterestID, ends_with(c("Con", "Sys", "Non", "Month")))) %>%
  distinct() %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# correlations >= 0.4 and sig
# TrtAreaSys & TrtAreaCon = 0.4
# TrtFreqSys & TrtFreqCon: 0.6
# visual patterns aren't super strong

methods_dat %>%
  select(c(AreaOfInterestID, ends_with(c("Q1", "Q2", "Q3", "Q4")))) %>%
  distinct() %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# correlations >= 0.4 and sig
# TrtFreqQ1 and Q2: 0.7 (visible)
# TrtFreqQ1 and Q3: 0.7 (visible)
# TrtFreqQ1 and Q4: 0.7 (visible)
# lots more --> just use average month

# quarter to set as intercept in categorical model
methods_dat %>%
  distinct(AreaOfInterestID, TrtMonthF) %>%
  count(TrtMonthF)

methods_dat2 <- methods_dat %>%
  mutate(across(.cols = (!starts_with("TrtMonth") & ends_with("F")), 
                .fns = ~ fct_relevel(.x, "none", "low")),
         TrtMonthF = fct_relevel(TrtMonthF, "Q2", "Q3", "Q4", "Q1"))


#### management efficacy models ####

# format PAC for beta distribution
target_dat_var2 <- target_dat_var %>%
  mutate(HydrProp = HydrPAC / 100,
         FloatProp = FloatPAC / 100)

# hydrilla binary model
hydr_mod1 <- betareg(HydrProp ~ HydrTrt, data = target_dat_var2)
plot(hydr_mod1)
summary(hydr_mod1)

# hydrilla categorical models
hydr_mod2 <- betareg(HydrProp ~ HydrTrtFreqF, data = target_dat_var2)
summary(hydr_mod2)
# linear
hydr_mod3 <- betareg(HydrProp ~ HydrTrtAreaF, data = target_dat_var2)
summary(hydr_mod3)
# linear

# hydrilla continuous model
hydr_mod4 <- betareg(HydrProp ~ HydrTrtFreq + HydrTrtArea, data = target_dat_var2)
plot(hydr_mod4)
summary(hydr_mod4)

# floating plant binary model
float_mod1 <- betareg(FloatProp ~ FloatTrt, data = target_dat_var2)
plot(float_mod1)
summary(float_mod1)

# floating plant categorical models
float_mod2 <- betareg(FloatProp ~ FloatTrtFreqF, data = target_dat_var2)
summary(float_mod2)
# linear
float_mod3 <- betareg(FloatProp ~ FloatTrtAreaF, data = target_dat_var2)
summary(float_mod3)
# linear

# floating plant continuous model
float_mod4 <- betareg(FloatProp ~ FloatTrtFreq + FloatTrtArea, data = target_dat_var2)
plot(float_mod4)
summary(float_mod4)


#### native richness methods model ####

#### start here: remove management intercepts ####

# response variable
ggplot(methods_dat, aes(x = NativeRichness)) +
  geom_density()

# fit frequency model
method_mod1 <- glmmTMB(NativeRichness ~ Time + Time:(TrtFreqConF + TrtFreqSysF + TrtFreqNonF +
                                                       TrtMonthF) + (1|AreaOfInterestID),
                    data = methods_dat2,
                    family = poisson)
method_res1 <- simulateResiduals(method_mod1, n = 1000)
plot(method_res1)
summary(method_mod1)
# TrtFreqCon same for low and high
# TrtFreqSys is negative then positive, but negative effect is very small
# TrtFreqNon higher with low
# month is maybe quadratic

# fit area model
method_mod2 <- glmmTMB(NativeRichness ~ Time + Time:(TrtAreaConF + TrtAreaSysF + TrtAreaNonF +
                                                       TrtMonthF) + (1|AreaOfInterestID),
                       data = methods_dat2,
                       family = poisson)
method_res2 <- simulateResiduals(method_mod2, n = 1000)
plot(method_res2)
summary(method_mod2)
# TrtAreaCon lower with low than high, but somewhat similar
# TrtAreaSys higher with low than high

# get average treatment month to standardize
trt_month_avg <- methods_dat2 %>%
  distinct(AreaOfInterestID, TrtMonth) %>%
  pull(TrtMonth) %>%
  mean()

# transform variables
methods_dat3 <- methods_dat2 %>%
  mutate(TrtFreqConLog = log(TrtFreqCon + 0.01),
         TrtFreqNonSq = TrtFreqNon^2,
         TrtAreaConLog = log(TrtAreaCon + 0.01),
         TrtAreaSysSq = TrtAreaSys^2,
         TrtMonthStd = TrtMonth - trt_month_avg,
         TrtMonthStdSq = TrtMonthStd^2)

# check correlations
methods_dat3 %>%
  distinct(AreaOfInterestID, TrtAreaConLog, TrtAreaSys, TrtAreaSysSq, TrtAreaNon,
           TrtFreqConLog, TrtFreqSys, TrtFreqNon, TrtFreqNonSq,
           TrtMonthStd, TrtMonthStdSq) %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# TreatAreaConLog + TrtFreqConLog: 0.5, driven by zeros

# fit continuous model
method_mod3 <- glmmTMB(NativeRichness ~ Time + Time:(TrtAreaConLog + TrtAreaSys + TrtAreaSysSq + TrtAreaNon +
                                                TrtFreqConLog + TrtFreqSys + TrtFreqNon + TrtFreqNonSq +
                                                TrtMonthStd + TrtMonthStdSq) + (1|AreaOfInterestID),
                       data = methods_dat3,
                       family = poisson,
                       control=glmmTMBControl(optimizer=optim,
                                              optArgs=list(method="BFGS")))
method_res3 <- simulateResiduals(method_mod3, n = 1000)
plot(method_res3)
summary(method_mod3)

# better fit with negative binomial?
method_mod3a <- update(method_mod3, family = "nbinom2")

# save model (slow to fit)
save(method_mod3a, file = "output/native_richness_methods_model.rda")
load("output/native_richness_methods_model.rda")

# model diagnostics
method_res3a <- simulateResiduals(method_mod3a, n = 1000)
plot(method_res3a)
summary(method_mod3a)


#### taxon-specific methods models ####

# select native taxa
# standardize treatment month using average from above
methods_taxa_dat2 <- methods_taxa_dat %>%
  filter(Origin == "Native") %>%
  mutate(TrtFreqConLog = log(TrtFreqCon + 0.01),
         TrtFreqNonSq = TrtFreqNon^2,
         TrtAreaConLog = log(TrtAreaCon + 0.01),
         TrtAreaSysSq = TrtAreaSys^2,
         TrtMonthStd = TrtMonth - trt_month_avg,
         TrtMonthStdSq = TrtMonthStd^2)

# check taxa for all 1's or 0's
taxa_methods_sum <- methods_taxa_dat2 %>%
  group_by(TaxonName) %>%
  summarize(Detections = sum(Detected) / n(),
            .groups = "drop")

arrange(taxa_methods_sum, Detections)
arrange(taxa_methods_sum, desc(Detections))
# check the two extremes

# list of taxa
methods_taxa <- sort(unique(methods_taxa_dat2$TaxonName))

# set i
i <- 1

# subset data
methods_taxa_sub <- filter(methods_taxa_dat2, TaxonName == methods_taxa[i])

# fit model
methods_taxa_mod <- glmmTMB(Detected ~ Time + Time:(TrtAreaConLog + TrtAreaSys + TrtAreaSysSq + TrtAreaNon +
                                               TrtFreqConLog + TrtFreqSys + TrtFreqNon + TrtFreqNonSq +
                                               TrtMonthStd + TrtMonthStdSq) + (1|AreaOfInterestID),
                                  data = methods_taxa_sub, family = "binomial")

# use this to manually cycle through taxa (some models had warnings below)
# i <- i + 1
# taxa with model convergence issues:
# 22: Fontinalis spp.
# 29: Juncus roemerianus
# 43: Najas marina
# 50: Orontium aquaticum
# 60: Proserpinaca spp.
# 62: Ruppia maritima
# 64: Sagittaria kurziana
# 70: Sparganium americanum
# 72: Stuckenia pectinata

# # look at models
summary(methods_taxa_mod)
methods_taxa_res <- simulateResiduals(methods_taxa_mod, n = 1000)
plot(methods_taxa_res, title = methods_taxa[i])

# save
save(methods_taxa_mod, file = paste0("output/methods_model_",
                                           str_to_lower(methods_taxa[i]) %>%
                                             str_replace_all(" ", "_"),
                                           ".rda"))

# save results when i <- 1
methods_taxa_coefs <- tidy(methods_taxa_mod) %>%
  filter(str_detect(term, "Time:")) %>%
  mutate(TaxonName = methods_taxa[1]) %>%
  relocate(TaxonName)

# loop through taxa
pdf("output/taxa_specific_methods_models.pdf")

for(i in methods_taxa[-c(22, 29, 43, 50, 60, 62, 64, 70, 72)]) {
  
  # subset data
  methods_taxa_sub <- filter(methods_taxa_dat2, TaxonName == i)
  
  # fit model
  methods_taxa_mod <- glmmTMB(Detected ~ Time + Time:(TrtAreaConLog + TrtAreaSys + TrtAreaSysSq + TrtAreaNon +
                                                 TrtFreqConLog + TrtFreqSys + TrtFreqNon + TrtFreqNonSq +
                                                 TrtMonthStd + TrtMonthStdSq) + (1|AreaOfInterestID),
                              data = methods_taxa_sub, family = "binomial")
  
  # extract and plot residuals
  methods_taxa_res <- simulateResiduals(methods_taxa_mod, n = 1000)
  plot(methods_taxa_res, title = i)
  
  # save model
  save(methods_taxa_mod, file = paste0("output/methods_model_",
                                       str_to_lower(i) %>%
                                         str_remove_all("\\.|\\(|\\)")%>%
                                         str_replace_all("/|, | ", "_"),
                                       ".rda"))
  
  # save model coefficients
  methods_taxa_coefs <- methods_taxa_coefs %>%
    full_join(tidy(methods_taxa_mod) %>%
                filter(str_detect(term, "Time:")) %>%
                mutate(TaxonName = i) %>%
                relocate(TaxonName))
}

dev.off()

# correct p-values
# remove taxon that the estimates didn't work
methods_taxa_coefs2 <- methods_taxa_coefs %>%
  mutate(q.value = p.adjust(p.value, method = "fdr"),
         lower = estimate - 1.96 * std.error,
         upper = estimate + 1.96 * std.error,
         odds_perc = 100 * (exp(estimate) - 1),
         lower_odds_perc = 100 * (exp(lower) - 1),
         upper_odds_perc = 100 * (exp(upper) - 1))

# save
write_csv(methods_taxa_coefs2, "output/methods_model_taxa_interaction_coefficients.csv")
# import if needed
# methods_taxa_coefs2 <- read_csv("output/methods_model_taxa_interaction_coefficients.csv")

# contact herbicide frequency
methods_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:TrtFreqConLog") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(methods_taxa_dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(Habitat = tolower(Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
write_csv("output/methods_model_taxa_contact_frequency_interaction_table.csv")

# systemic herbicide intensity
methods_taxa_coefs2 %>%
  filter(term == "Time:TrtAreaSys") %>%
  select(TaxonName, estimate, std.error, q.value) %>%
  rename_with(~ paste0("intensity_", .), .cols = -TaxonName) %>%
  right_join(methods_taxa_coefs2 %>%
               filter(q.value < 0.05 & term == "Time:TrtAreaSysSq") %>% 
               select(TaxonName, estimate, std.error, q.value) %>%
               rename_with(~ paste0("intensity2_", .), .cols = -TaxonName)) %>%
  full_join(methods_taxa_coefs2 %>%
              filter(q.value < 0.05 & term == "Time:TrtAreaSys") %>%
              select(TaxonName, estimate, std.error, q.value) %>%
              rename_with(~ paste0("intensity_", .), .cols = -TaxonName) %>%
              left_join(methods_taxa_coefs2 %>%
                           filter(term == "Time:TrtAreaSysSq") %>% 
                           select(TaxonName, estimate, std.error, q.value) %>%
                           rename_with(~ paste0("intensity2_", .), .cols = -TaxonName))) %>%
  left_join(methods_taxa_dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(Habitat = tolower(Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
  write_csv("output/methods_model_taxa_systemic_intensity_interaction_table.csv")

# non-herbicide frequency
methods_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:TrtFreqNon") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  rename_with(~ paste0("frequency_", .), .cols = -TaxonName) %>%
  left_join(methods_taxa_coefs2 %>%
              filter(term == "Time:TrtFreqNonSq") %>% 
              select(TaxonName, estimate, std.error, q.value) %>%
              rename_with(~ paste0("frequency2_", .), .cols = -TaxonName)) %>%
  left_join(methods_taxa_dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(Habitat = tolower(Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
  write_csv("output/methods_model_taxa_nonherbicide_frequency_interaction_table.csv")


#### native richness target model ####

# response variable
ggplot(target_dat2, aes(x = NativeRichness)) +
  geom_density()

# fit frequency model
target_mod1 <- glmmTMB(NativeRichness ~ Time + Time:(HydrPACF + HydrTrtFreqF + FloatPACF + FloatTrtFreqF + 
                                             OtherTrtFreqF) + (1|AreaOfInterestID),
                    data = target_dat2,
                    family = poisson)
target_res1 <- simulateResiduals(target_mod1, n = 1000)
plot(target_res1)
summary(target_mod1)
# interactions with time:
# floating PAC changes directions
# Hydr PAC similar low and high
# others seem linear

# fit area model
target_mod2 <- glmmTMB(NativeRichness ~ Time + Time:(HydrPACF + HydrTrtAreaF + FloatPACF + FloatTrtAreaF + 
                                                    OtherTrtAreaF) + (1|AreaOfInterestID),
                    data = target_dat2,
                    family = poisson)
target_res2 <- simulateResiduals(target_mod2, n = 1000)
plot(target_res2)
summary(target_mod2)
# interactions with time:
# floating area similar with low and high
# other area positive then very low

# transform variables with nonlinearities
target_dat3 <- target_dat2 %>%
  mutate(FloatPACSq = FloatPAC^2,
         HydrPACLog = log(HydrPAC + 0.01),
         FloatTrtAreaLog = log(FloatTrtArea + 0.01),
         OtherTrtAreaSq = OtherTrtArea^2)

# check correlations
target_dat3 %>%
  distinct(AreaOfInterestID, HydrPACLog, FloatPAC, FloatPACSq,
             HydrTrtFreq,FloatTrtFreq, OtherTrtFreq, 
             HydrTrtArea, FloatTrtAreaLog, 
             OtherTrtArea, OtherTrtAreaSq) %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# HydrTrtFreq and HydrPACLog: 0.6
# HydrTrtArea and HydrPACLog: 0.6 (visible)
# OtherTrtFreq and FloatTrtFreq: 0.6
# FloatTrtAreaLog and FloatTrtFreq: 0.5

# take out transformation of hydrPAC
target_dat3 %>%
  distinct(AreaOfInterestID, HydrPAC, FloatPAC, FloatPACSq,
           HydrTrtFreq,FloatTrtFreq, OtherTrtFreq, 
           HydrTrtArea, FloatTrtAreaLog, 
           OtherTrtArea, OtherTrtAreaSq) %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# no longer correlated

# fit model with continuous variables
target_mod3 <- glmmTMB(NativeRichness ~ Time + Time:(HydrPAC + FloatPAC + FloatPACSq +
                                                       HydrTrtFreq +FloatTrtFreq + OtherTrtFreq + 
                                                       HydrTrtArea + FloatTrtAreaLog + 
                                                       OtherTrtArea + OtherTrtAreaSq) +
                      (1|AreaOfInterestID),
                    data = target_dat3,
                    family = "poisson",
                    control=glmmTMBControl(optimizer=optim,
                                           optArgs=list(method="BFGS")))

target_res3 <- simulateResiduals(target_mod3, n = 1000)
plot(target_res3)
summary(target_mod3)

# better fit with negative binomial?
target_mod3a <- update(target_mod3, family = "nbinom2")

# save model (slow to fit)
save(target_mod3a, file = "output/native_richness_target_model.rda")
load("output/native_richness_target_model.rda")

# model diagnostics
target_res3a <- simulateResiduals(target_mod3a, n = 1000)
plot(target_res3a)
summary(target_mod3a)
# very large dispersion parameter
# results are very similar


#### taxon-specific target models ####

# select native taxa
# add columns
target_taxa_dat2 <- target_taxa_dat %>%
  filter(Origin == "Native") %>%
  mutate(FloatPACSq = FloatPAC^2,
         HydrPACLog = log(HydrPAC + 0.01),
         FloatTrtAreaLog = log(FloatTrtArea + 0.01),
         OtherTrtAreaSq = OtherTrtArea^2)

# check taxa for all 1's or 0's
taxa_target_sum <- target_taxa_dat2 %>%
  group_by(TaxonName) %>%
  summarize(Detections = sum(Detected) / n(),
            .groups = "drop")

arrange(taxa_target_sum, Detections)
arrange(taxa_target_sum, desc(Detections))
# check the two extremes

# list of taxa
target_taxa <- sort(unique(target_taxa_dat2$TaxonName))

# set i
i <- 1

# subset data
target_taxa_sub <- filter(target_taxa_dat2, TaxonName == target_taxa[i])

# fit model
target_taxa_mod <- glmmTMB(Detected ~ Time + Time:(HydrPAC + FloatPAC + FloatPACSq +
                                              HydrTrtFreq +FloatTrtFreq + OtherTrtFreq + 
                                              HydrTrtArea + FloatTrtAreaLog + 
                                              OtherTrtArea + OtherTrtAreaSq) + (1|AreaOfInterestID),
                            data = target_taxa_sub, family = "binomial")

# use this to manually cycle through taxa (some models had warnings below)
i <- i + 1
# taxa with model convergence issues:


# # look at models
summary(target_taxa_mod)
target_taxa_res <- simulateResiduals(target_taxa_mod, n = 1000)
plot(target_taxa_res, title = target_taxa[i])

# save
save(target_taxa_mod, file = paste0("output/target_model_",
                                     str_to_lower(target_taxa[i]) %>%
                                       str_replace_all(" ", "_"),
                                     ".rda"))

# save results when i <- 1
target_taxa_coefs <- tidy(target_taxa_mod) %>%
  filter(str_detect(term, "Time:")) %>%
  mutate(TaxonName = target_taxa[1]) %>%
  relocate(TaxonName)

# loop through taxa
pdf("output/taxa_specific_target_models.pdf")

for(i in target_taxa) {
  
  # subset data
  target_taxa_sub <- filter(target_taxa_dat2, TaxonName == i)
  
  # fit model
  target_taxa_mod <- glmmTMB(Detected ~ Time + Time:(HydrPAC + FloatPAC + FloatPACSq +
                                                HydrTrtFreq +FloatTrtFreq + OtherTrtFreq + 
                                                HydrTrtArea + FloatTrtAreaLog + 
                                                OtherTrtArea + OtherTrtAreaSq) + (1|AreaOfInterestID),
                              data = target_taxa_sub, family = "binomial")
  
  # extract and plot residuals
  target_taxa_res <- simulateResiduals(target_taxa_mod, n = 1000)
  plot(target_taxa_res, title = i)
  
  # save model
  save(target_taxa_mod, file = paste0("output/target_model_",
                                       str_to_lower(i) %>%
                                         str_remove_all("\\.|\\(|\\)")%>%
                                         str_replace_all("/|, | ", "_"),
                                       ".rda"))
  
  # save model coefficients
  target_taxa_coefs <- target_taxa_coefs %>%
    full_join(tidy(target_taxa_mod) %>%
                filter(str_detect(term, "Time:")) %>%
                mutate(TaxonName = i) %>%
                relocate(TaxonName))
}

dev.off()

# correct p-values
# remove taxon that the estimates didn't work
target_taxa_coefs2 <- target_taxa_coefs %>%
  mutate(q.value = p.adjust(p.value, method = "fdr"),
         lower = estimate - 1.96 * std.error,
         upper = estimate + 1.96 * std.error,
         odds_perc = 100 * (exp(estimate) - 1),
         lower_odds_perc = 100 * (exp(lower) - 1),
         upper_odds_perc = 100 * (exp(upper) - 1))

# save
write_csv(target_taxa_coefs2, "output/target_model_taxa_interaction_coefficients.csv")
# import if needed
# target_taxa_coefs2 <- read_csv("output/target_model_taxa_interaction_coefficients.csv")

# hydrilla treatment frequency
target_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:HydrTrtFreq") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(target_taxa_dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(Habitat = tolower(Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
  write_csv("output/target_model_taxa_hydrilla_frequency_interaction_table.csv")

# floating PAC
target_taxa_coefs2 %>%
  filter(term == "Time:FloatPAC") %>%
  select(TaxonName, estimate, std.error, q.value) %>%
  rename_with(~ paste0("pac_", .), .cols = -TaxonName) %>%
  right_join(target_taxa_coefs2 %>%
               filter(q.value < 0.05 & term == "Time:FloatPACSq") %>% 
               select(TaxonName, estimate, std.error, q.value) %>%
               rename_with(~ paste0("pac2_", .), .cols = -TaxonName)) %>%
  full_join(target_taxa_coefs2 %>%
              filter(q.value < 0.05 & term == "Time:FloatPAC") %>%
              select(TaxonName, estimate, std.error, q.value) %>%
              rename_with(~ paste0("pac_", .), .cols = -TaxonName) %>%
              left_join(target_taxa_coefs2 %>%
                          filter(term == "Time:FloatPACSq") %>% 
                          select(TaxonName, estimate, std.error, q.value) %>%
                          rename_with(~ paste0("pac2_", .), .cols = -TaxonName))) %>%
  left_join(target_taxa_dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(Habitat = tolower(Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
  write_csv("output/target_model_taxa_floating_PAC_interaction_table.csv")

# floating treatment frequency
target_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:FloatTrtFreq") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(target_taxa_dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(Habitat = tolower(Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
  write_csv("output/target_model_taxa_floating_frequency_interaction_table.csv")

# floating treatment intensity
target_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:FloatTrtAreaLog") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(target_taxa_dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(Habitat = tolower(Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
  write_csv("output/target_model_taxa_floating_intensity_interaction_table.csv")

# other treatment frequency
target_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:OtherTrtFreq") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(target_taxa_dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(Habitat = tolower(Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
  write_csv("output/target_model_taxa_other_frequency_interaction_table.csv")

# floating PAC
target_taxa_coefs2 %>%
  filter(term == "Time:OtherTrtArea") %>%
  select(TaxonName, estimate, std.error, q.value) %>%
  rename_with(~ paste0("intensity_", .), .cols = -TaxonName) %>%
  right_join(target_taxa_coefs2 %>%
               filter(q.value < 0.05 & term == "Time:OtherTrtAreaSq") %>% 
               select(TaxonName, estimate, std.error, q.value) %>%
               rename_with(~ paste0("intensity2_", .), .cols = -TaxonName)) %>%
  full_join(target_taxa_coefs2 %>%
              filter(q.value < 0.05 & term == "Time:OtherTrtArea") %>%
              select(TaxonName, estimate, std.error, q.value) %>%
              rename_with(~ paste0("intensity_", .), .cols = -TaxonName) %>%
              left_join(target_taxa_coefs2 %>%
                          filter(term == "Time:OtherTrtAreaSq") %>% 
                          select(TaxonName, estimate, std.error, q.value) %>%
                          rename_with(~ paste0("intensity2_", .), .cols = -TaxonName))) %>%
  left_join(target_taxa_dat %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(Habitat = tolower(Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
  write_csv("output/target_model_taxa_other_intensity_interaction_table.csv")


#### native richness predicted values ####

# remove time component
mod_dat_var <- methods_dat3 %>%
  distinct(AreaOfInterestID, TrtAreaCon, TrtAreaSys, TrtAreaNon,
           TrtFreqCon, TrtFreqSys, TrtFreqNon, 
           TrtMonth) %>%
  full_join(target_dat3 %>%
  distinct(AreaOfInterestID, HydrPAC, FloatPAC,
           HydrTrtFreq, FloatTrtFreq, OtherTrtFreq,
           HydrTrtArea, FloatTrtArea, OtherTrtArea))

# function for predictions
pred_fun <- function(variable, model, dat){
  
  pred_dat <- tibble(AreaOfInterestID = "A",
                     Time = min(dat$Time):max(dat$Time),
                     TrtAreaCon = 0, 
                     TrtAreaSys = 0, 
                     TrtAreaNon = 0,
                     TrtFreqCon = 0, 
                     TrtFreqSys = 0, 
                     TrtFreqNon = 0, 
                     TrtMonthStd = 0,
                     HydrPAC = mean(mod_dat_var$HydrPAC), 
                     FloatPAC = mean(mod_dat_var$FloatPAC),
                     HydrTrtFreq = 0, 
                     FloatTrtFreq = 0, 
                     OtherTrtFreq = 0,
                     HydrTrtArea = 0, 
                     FloatTrtArea = 0, 
                     OtherTrtArea = 0) %>%
    select(-{{variable}}) %>%
    expand_grid(mod_dat_var%>%
                  filter(!is.na({{variable}}))  %>%
                  select({{variable}})) %>%
    mutate(across(.cols = -c(AreaOfInterestID, Time), 
                  .fns = list(Sq = ~.x^2, Log = ~log(.x + 0.01)),
                  .names = "{.col}{.fn}")) %>%
    mutate(Pred = predict(model, newdata = ., type = "response",
                          allow.new.levels = T),
           PredSE = predict(model, newdata = ., type = "response",
                          allow.new.levels = T, se.fit = T)$se.fit)
  
  return(pred_dat)
  
}

# contact herb freq
pred_cont_freq <- pred_fun(TrtFreqCon, method_mod3a, methods_dat3)
method_fig_1 <- ggplot(pred_cont_freq, aes(x = Time, y = Pred, 
                              color = TrtFreqCon,
                              group = TrtFreqCon)) +
  geom_line() +
  scale_color_viridis_c(option = "A", 
                        name = "Contact\nherbicide\nfrequency") +
  labs(x = "Years of monitoring", y = "Native taxonomic richness") +
  def_theme_paper

# systemic herb intensity
pred_sys_ints <- pred_fun(TrtAreaSys, method_mod3a, methods_dat3)
method_fig_2 <- ggplot(pred_sys_ints, aes(x = Time, y = Pred, 
                                           color = TrtAreaSys,
                                           group = TrtAreaSys)) +
  geom_line() +
  scale_color_viridis_c(option = "H", 
                        name = "Systemic\nherbicide\nintensity") +
  labs(x = "Years of monitoring", y = "Native taxonomic richness") +
  def_theme_paper

# hydr PAC
nat_pred_hydr_pac <- pred_fun(HydrPAC, target_mod3a, target_dat3)
ggplot(nat_pred_hydr_pac, aes(x = Time, y = Pred, 
                              color = HydrPAC,
                              group = HydrPAC)) +
  geom_line()

# floating PAC
nat_pred_float_pac <- pred_fun(FloatPAC, target_mod3a, target_dat3)
ggplot(nat_pred_float_pac, aes(x = Time, y = Pred, 
                              color = FloatPAC,
                              group = FloatPAC)) +
  geom_line()

# floating PAC trt freq
nat_pred_float_freq <- pred_fun(FloatTrtFreq, target_mod3a, target_dat3)
ggplot(nat_pred_float_freq, aes(x = Time, y = Pred, 
                               color = FloatTrtFreq,
                               group = FloatTrtFreq)) +
  geom_line()

# other PAC trt area
nat_pred_other_area <- pred_fun(OtherTrtArea, target_mod3a, target_dat3)
ggplot(nat_pred_other_area, aes(x = Time, y = Pred, 
                                color = OtherTrtArea,
                                group = OtherTrtArea)) +
  geom_line()

# area treated with systemic herbicide
nat_pred_sys_area <- pred_fun(TrtAreaSys, method_mod3a, methods_dat3)
ggplot(nat_pred_sys_area, aes(x = Time, y = Pred, 
                                color = TrtAreaSys,
                                group = TrtAreaSys)) +
  geom_line()

# frequency treated with non-herbicide
nat_pred_non_freq <- pred_fun(TrtFreqNon, method_mod3a, methods_dat3)
ggplot(nat_pred_non_freq, aes(x = Time, y = Pred, 
                              color = TrtFreqNon,
                              group = TrtFreqNon)) +
  geom_line()




#### native richness coefficient plots ####

change_sys_ints <- confint(method_mod3a, parm = c("Time:TrtAreaSys", "Time:TrtAreaSysSq")) %>%
   cbind(tibble(Term = c("b", "a"))) %>%
   select(Term, Estimate) %>%
   pivot_wider(names_from = Term, values_from = Estimate) %>%
   expand_grid(tibble(Intensity = sort(unique(methods_dat3$TrtAreaSys)))) %>%
   mutate(ChangeOverTime = a * Intensity^2 + b * Intensity)

change_sys_ints_fig <- ggplot(change_sys_ints, aes(x = Intensity, y = ChangeOverTime)) +
  geom_point() +
  labs(x = "Systemic herbicide intensity", y = "Relative change in native richness over time") + 
  def_theme_paper


