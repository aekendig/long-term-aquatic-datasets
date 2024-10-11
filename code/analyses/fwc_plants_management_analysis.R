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
methods_dat <- read_csv("intermediate-data/FWC_plant_management_methods_analysis_formatted.csv")
rich_dat <- read_csv("intermediate-data/FWC_plant_management_richness_analysis_formatted.csv")
taxa_dat <- read_csv("intermediate-data/FWC_plant_management_taxa_analysis_formatted.csv")
methods_taxa_dat <- read_csv("intermediate-data/FWC_plant_management_methods_taxa_analysis_formatted.csv")


#### examine full dataset ####
rich_dat_var <- rich_dat %>%
  select(AreaOfInterestID, HydrPAC, FloatPAC, HydrTrtFreq, HydrTrtArea, FloatTrtFreq, FloatTrtArea,
           OtherTrtFreq, OtherTrtArea, ends_with("F")) %>%
  distinct() %>%
  mutate(HydrTrt = if_else(HydrTrtFreq == 0, 0, 1),
         FloatTrt = if_else(FloatTrtFreq == 0, 0, 1),
         across(.cols = ends_with("F"), 
                .fns = ~ fct_relevel(.x, "none", "low")))

# correlations among explanatory variables
rich_dat_var %>%
  select(-c(AreaOfInterestID, HydrTrt, FloatTrt, ends_with("F"))) %>%
  ggpairs()
# correlations >= 0.4 and sig
# OtherTrtFreq and FloatTrtFreq: 0.6
# OtherTrtArea and FloatTrtArea: 0.4
# neither looks like a particularly strong relationship

# order factors
rich_dat2 <- rich_dat %>%
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
rich_dat_var2 <- rich_dat_var %>%
  mutate(HydrProp = HydrPAC / 100,
         FloatProp = FloatPAC / 100)

# hydrilla binary model
hydr_mod1 <- betareg(HydrProp ~ HydrTrt, data = rich_dat_var2)
plot(hydr_mod1)
summary(hydr_mod1)

# hydrilla categorical models
hydr_mod2 <- betareg(HydrProp ~ HydrTrtFreqF, data = rich_dat_var2)
summary(hydr_mod2)
# linear
hydr_mod3 <- betareg(HydrProp ~ HydrTrtAreaF, data = rich_dat_var2)
summary(hydr_mod3)
# linear

# hydrilla continuous model
hydr_mod4 <- betareg(HydrProp ~ HydrTrtFreq + HydrTrtArea, data = rich_dat_var2)
plot(hydr_mod4)
summary(hydr_mod4)

# floating plant binary model
float_mod1 <- betareg(FloatProp ~ FloatTrt, data = rich_dat_var2)
plot(float_mod1)
summary(float_mod1)

# floating plant categorical models
float_mod2 <- betareg(FloatProp ~ FloatTrtFreqF, data = rich_dat_var2)
summary(float_mod2)
# linear
float_mod3 <- betareg(FloatProp ~ FloatTrtAreaF, data = rich_dat_var2)
summary(float_mod3)
# linear

# floating plant continuous model
float_mod4 <- betareg(FloatProp ~ FloatTrtFreq + FloatTrtArea, data = rich_dat_var2)
plot(float_mod4)
summary(float_mod4)


#### native richness methods model ####

# response variable
ggplot(methods_dat, aes(x = NativeRichness)) +
  geom_density()

# fit frequency model
method_mod1 <- glmmTMB(NativeRichness ~ Time*(TrtFreqConF + TrtFreqSysF + TrtFreqNonF +
                                                       TrtMonthF) + (1|AreaOfInterestID),
                    data = methods_dat2,
                    family = poisson)
method_res1 <- simulateResiduals(method_mod1, n = 1000)
plot(method_res1)
summary(method_mod1)
# TrtFreqCon same for low and high
# TrtFreqSys is negative then positive, but negative effect is very small
# TrtFreqNon only positive with low
# month is maybe quadratic

# fit area model
method_mod2 <- glmmTMB(NativeRichness ~ Time*(TrtAreaConF + TrtAreaSysF + TrtAreaNonF +
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
method_mod3 <- glmmTMB(NativeRichness ~ Time*(TrtAreaConLog + TrtAreaSys + TrtAreaSysSq + TrtAreaNon +
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
methods_taxa_mod <- glmmTMB(Detected ~ Time*(TrtAreaConLog + TrtAreaSys + TrtAreaSysSq + TrtAreaNon +
                                               TrtFreqConLog + TrtFreqSys + TrtFreqNon + TrtFreqNonSq +
                                               TrtMonthStd + TrtMonthStdSq) + (1|AreaOfInterestID),
                                  data = methods_taxa_sub, family = "binomial")

# use this to manually cycle through taxa (some models had warnings below)
# i <- i + 1
# taxa with model convergence issues:
# 22: Fontinalis spp.
# 43: Najas marina
# 50: Orontium aquaticum
# 60: Proserpinaca spp.
# 62: Ruppia maritima
# 64: Sagittaria kurziana
# 70: Sparganium americanum
# 72: Stuckenia pectinata

# other notes: 
# covariance matrix estimated for Juncus roemerianus uses potentiall less
# accurate method, so standard errors may be off

# # look at models
summary(methods_taxa_mod)
methods_taxa_res <- simulateResiduals(methods_taxa_mod, n = 1000)
plot(methods_taxa_res, title = methods_taxa[i])

# # save
# save(methods_taxa_mod, file = paste0("output/methods_model_",
#                                            str_to_lower(methods_taxa[i]) %>%
#                                              str_replace_all(" ", "_"),
#                                            ".rda"))

# save results when i <- 1
methods_taxa_coefs <- tidy(methods_taxa_mod) %>%
  filter(str_detect(term, "Time:")) %>%
  mutate(TaxonName = methods_taxa[1]) %>%
  relocate(TaxonName)

# loop through taxa
pdf("output/taxa_specific_methods_models.pdf")

for(i in methods_taxa[-c(22, 43, 50, 60, 62, 64, 70, 72)]) {
  
  # subset data
  methods_taxa_sub <- filter(methods_taxa_dat2, TaxonName == i)
  
  # fit model
  methods_taxa_mod <- glmmTMB(Detected ~ Time*(TrtAreaConLog + TrtAreaSys + TrtAreaSysSq + TrtAreaNon +
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
  filter(TaxonName != "Juncus roemerianus") %>%
  mutate(q.value = p.adjust(p.value, method = "fdr"),
         lower = estimate - 1.96 * std.error,
         upper = estimate + 1.96 * std.error,
         odds_perc = 100 * (exp(estimate) - 1),
         lower_odds_perc = 100 * (exp(lower) - 1),
         upper_odds_perc = 100 * (exp(upper) - 1))

# save
write_csv(methods_taxa_coefs2, "output/methods_model_taxa_interaction_coefficients.csv")

# contact herbicide frequency
methods_taxa_coefs2 %>%
  filter(q.value < 0.05 & estimate < 0 & term == "Time:TrtFreqConLog") %>% 
  select(TaxonName, estimate, std.error, q.value)

# systemic herbicide intensity
methods_taxa_coefs2 %>%
  filter(q.value < 0.05 & 
           (estimate < 0 & term == "Time:TrtAreaSys") |
           (estimate > 0 & term == "Time:TrtAreaSysSq"))%>% 
  select(TaxonName, term, estimate, std.error, q.value) %>%
  group_by(TaxonName) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 2) %>%
  select(-n) %>%
  left_join(distinct(methods_taxa_dat2, TaxonName, Habitat)) %>%
  mutate(term = str_replace_all(term, "Time:TrtAreaSysSq", "intensity2") %>%
           str_replace_all("Time:TrtAreaSys", "intensity"),
         estimate = round_half_up(estimate, 4),
         std.error = round_half_up(std.error, 4),
         q.value = round_half_up(q.value, 3),
         Habitat = str_to_lower(Habitat)) %>%
  pivot_wider(names_from = term, 
              values_from = c(estimate, std.error, q.value),
              names_glue = "{term}_{.value}") %>%
  write_csv("output/methods_model_taxa_systemic_intensity_interaction_table.csv")

# non-herbicide frequency
methods_taxa_coefs2 %>%
  filter(q.value < 0.05 & 
           (estimate > 0 & term == "Time:TrtFreqNon") |
           (estimate < 0 & term == "Time:TrtFreqNonSq"))%>% 
  select(TaxonName, term, estimate, std.error, q.value) %>%
  group_by(TaxonName) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 2) %>%
  select(-n) %>%
  left_join(distinct(methods_taxa_dat2, TaxonName, Habitat)) %>%
  mutate(term = str_replace_all(term, "Time:TrtFreqNonSq", "frequency2") %>%
           str_replace_all("Time:TrtFreqNon", "frequency"),
         estimate = round_half_up(estimate, 4),
         std.error = round_half_up(std.error, 4),
         q.value = round_half_up(q.value, 3),
         Habitat = str_to_lower(Habitat)) %>%
  pivot_wider(names_from = term, 
              values_from = c(estimate, std.error, q.value),
              names_glue = "{term}_{.value}") %>%
  write_csv("output/methods_model_taxa_nonherbicide_frequency_interaction_table.csv")


#### native richness target model ####

# response variable
ggplot(rich_dat2, aes(x = NativeRichness)) +
  geom_density()

# fit frequency model
target_mod1 <- glmmTMB(NativeRichness ~ Time + Time*(HydrPACF + HydrTrtFreqF + FloatPACF + FloatTrtFreqF + 
                                             OtherTrtFreqF) + (1|AreaOfInterestID),
                    data = rich_dat2,
                    family = poisson)
target_res1 <- simulateResiduals(target_mod1, n = 1000)
plot(target_res1)
summary(target_mod1)
# interactions with time:
# floating PAC changes directions
# Hydr PAC similar low and high
# others seem linear

# fit area model
target_mod2 <- glmmTMB(NativeRichness ~ Time + Time*(HydrPACF + HydrTrtAreaF + FloatPACF + FloatTrtAreaF + 
                                                    OtherTrtAreaF) + (1|AreaOfInterestID),
                    data = rich_dat2,
                    family = poisson)
target_res2 <- simulateResiduals(target_mod2, n = 1000)
plot(target_res2)
summary(target_mod2)
# interactions with time:
# floating area similar with low and high
# other area changes directions

# transform variables with nonlinearities
rich_dat3 <- rich_dat2 %>%
  mutate(FloatPACSq = FloatPAC^2,
         HydrPACLog = log(HydrPAC + 0.01),
         FloatTrtAreaLog = log(FloatTrtArea + 0.01),
         OtherTrtAreaSq = OtherTrtArea^2)

# check correlations
rich_dat3 %>%
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
rich_dat3 %>%
  distinct(AreaOfInterestID, HydrPAC, FloatPAC, FloatPACSq,
           HydrTrtFreq,FloatTrtFreq, OtherTrtFreq, 
           HydrTrtArea, FloatTrtAreaLog, 
           OtherTrtArea, OtherTrtAreaSq) %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# no longer correlated

# fit model with continuous variables
target_mod3 <- glmmTMB(NativeRichness ~ Time + Time*(HydrPAC + FloatPAC + FloatPACSq +
                                                       HydrTrtFreq +FloatTrtFreq + OtherTrtFreq + 
                                                       HydrTrtArea + FloatTrtAreaLog + 
                                                       OtherTrtArea + OtherTrtAreaSq) +
                      (1|AreaOfInterestID),
                    data = rich_dat3,
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


#### native richness predicted values ####

# remove time component
mod_dat_var <- methods_dat3 %>%
  distinct(AreaOfInterestID, TrtAreaCon, TrtAreaSys, TrtAreaSysSq, TrtAreaNon,
           TrtFreqCon, TrtFreqConLog, TrtFreqSys, TrtFreqNon, TrtFreqNonSq,
           TrtMonth, TrtMonthSq) %>%
  full_join(rich_dat3 %>%
  distinct(AreaOfInterestID, HydrPAC, FloatPAC, FloatPACSq,
           HydrTrtFreq, FloatTrtFreq, OtherTrtFreq,
           HydrTrtArea, FloatTrtArea, FloatTrtAreaLog, OtherTrtArea, OtherTrtAreaSq))

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

# hydr PAC
nat_pred_hydr_pac <- pred_fun(HydrPAC, target_mod3a, rich_dat3)
ggplot(nat_pred_hydr_pac, aes(x = Time, y = Pred, 
                              color = HydrPAC,
                              group = HydrPAC)) +
  geom_line()

# floating PAC
nat_pred_float_pac <- pred_fun(FloatPAC, target_mod3a, rich_dat3)
ggplot(nat_pred_float_pac, aes(x = Time, y = Pred, 
                              color = FloatPAC,
                              group = FloatPAC)) +
  geom_line()

# floating PAC trt freq
nat_pred_float_freq <- pred_fun(FloatTrtFreq, target_mod3a, rich_dat3)
ggplot(nat_pred_float_freq, aes(x = Time, y = Pred, 
                               color = FloatTrtFreq,
                               group = FloatTrtFreq)) +
  geom_line()

# other PAC trt area
nat_pred_other_area <- pred_fun(OtherTrtArea, target_mod3a, rich_dat3)
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




#### native richness coefficients ####

# update: just do tidy table for supp and confint for coefficient figure (no transformations)

# extract coefficients and 95% CI
# transform to response scale by exponentiating (incidence rate ratio)
nat_est <- confint(nat_mod1a) %>%
  as_tibble() %>%
  mutate(Coef = rownames(confint(nat_mod1a)),
         Mod = "Freq") %>%
  left_join(tidy(nat_mod1a) %>%
              rename(Coef = term,
                     PValue = p.value) %>%
              select(Coef, PValue)) %>%
  full_join(confint(nat_mod2a) %>%
              as_tibble() %>%
              mutate(Coef = rownames(confint(nat_mod2a)),
                     Mod = "Area") %>%
              left_join(tidy(nat_mod2a) %>%
                          rename(Coef = term,
                                 PValue = p.value) %>%
                          select(Coef, PValue))) %>%
  filter(str_detect(Coef, "Time")) %>%
  rename(Lower = `2.5 %`,
         Upper = `97.5 %`) %>%
  mutate(CoefType = case_when(str_detect(Coef, "PAC") ~ "PAC",
                              str_detect(Coef, "Trt") ~ "Trt"),
         Level = case_when(str_detect(Coef, "low") ~ "Low",
                           str_detect(Coef, "high") ~ "High",
                           TRUE ~ "None") %>%
           fct_relevel("None", "Low", "High"),
         Target = case_when(str_detect(Coef, "Hydr") ~ "hydrilla",
                            str_detect(Coef, "Float") ~ "floating plants",
                            str_detect(Coef, "Other") ~ "other",
                            TRUE ~ "none") %>%
           fct_relevel("none", "hydrilla", "floating plants"),
         IRR = exp(Estimate), # assumes time = 1, relative to estimate with this variable = "none"
         LowerIRR = exp(Lower),
         UpperIRR = exp(Upper))


#### individual species models ####

# for each species, only choose waterbodies where it has been detected at least once within this dataset
# model is same as richness: change over time in presence/absence, variation in change over time due to PAC, trt freq, and trt area
# use same variable transformations as native richness model to keep things simpler
# bonferoni correction all p-values
# require X lakes to analyze species

