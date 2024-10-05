#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom.mixed)
library(glmmTMB)
library(DHARMa)
library(GGally)

# figure settings
source("code/settings/figure_settings.R")

# import data
methods_dat <- read_csv("intermediate-data/FWC_plant_management_methods_analysis_formatted.csv")
rich_dat <- read_csv("intermediate-data/FWC_plant_management_richness_analysis_formatted.csv")


#### examine methods data ####

# correlations among explanatory variables
methods_dat %>%
  select(c(AreaOfInterestID, ends_with(c("Con", "Sys", "Non", "Month")))) %>%
  distinct() %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# correlations >= 0.4 and sig
# TrtFreqSys & TrtFreqCon: 0.6
# higher values are correlated - driven by high values

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


#### native richness methods model ####

# response variable
ggplot(methods_dat, aes(x = NativeRichness)) +
  geom_density()

# fit frequency model
method_mod1 <- glmmTMB(NativeRichness ~ Time + Time*(TrtFreqConF + TrtFreqSysF + TrtFreqNonF +
                                                       TrtMonthF) + (1|AreaOfInterestID),
                    data = methods_dat2,
                    family = poisson)
method_res1 <- simulateResiduals(method_mod1, n = 1000)
plot(method_res1)
summary(method_mod1)
# TrtFreqCon same for low and high
# TrtFreqNon only positive with low
# month is maybe quadratic

# fit area model
method_mod2 <- glmmTMB(NativeRichness ~ Time + Time*(TrtAreaConF + TrtAreaSysF + TrtAreaNonF +
                                                       TrtMonthF) + (1|AreaOfInterestID),
                       data = methods_dat2,
                       family = poisson)
method_res2 <- simulateResiduals(method_mod2, n = 1000)
plot(method_res2)
summary(method_mod2)
# TrtAreaSys higher with low than high

# transform variables
methods_dat3 <- methods_dat2 %>%
  mutate(TrtFreqConLog = log(TrtFreqCon + 0.01),
         TrtFreqNonSq = TrtFreqNon^2,
         TrtAreaSysSq = TrtAreaSys^2,
         TrtMonthSq = TrtMonth^2)

# check correlations
methods_dat3 %>%
  distinct(AreaOfInterestID, TrtAreaCon, TrtAreaSys, TrtAreaSysSq, TrtAreaNon,
           TrtFreqConLog, TrtFreqSys, TrtFreqNon, TrtFreqNonSq,
           TrtMonth, TrtMonthSq) %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# only squared values

# fit continuous model
method_mod3 <- glmmTMB(NativeRichness ~ Time + Time*(TrtAreaCon + TrtAreaSys + TrtAreaSysSq + TrtAreaNon +
                                                       TrtFreqConLog + TrtFreqSys + TrtFreqNon + TrtFreqNonSq +
                                                       TrtMonth + TrtMonthSq) + (1|AreaOfInterestID),
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


#### examine target data ####

# correlations among explanatory variables
rich_dat %>%
  distinct(AreaOfInterestID, HydrPAC, FloatPAC, HydrTrtFreq, HydrTrtArea, FloatTrtFreq, FloatTrtArea,
           OtherTrtFreq, OtherTrtArea) %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# correlations >= 0.4 and sig
# OtherTrtFreq and FloatTrtFreq: 0.6
# OtherTrtArea and FloatTrtArea: 0.4
# neither looks like a particularly strong relationship

# order factors
rich_dat2 <- rich_dat %>%
  mutate(across(.cols = ends_with("F"), 
                .fns = ~ fct_relevel(.x, "none", "low")))


#### native richness target model ####

# response variable
ggplot(rich_dat, aes(x = NativeRichness)) +
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
rich_dat_var <- rich_dat3 %>%
  distinct(AreaOfInterestID, HydrPAC, HydrTrtFreq, HydrTrtArea, HydrTrtAreaLog,
           FloatPAC, FloatPACLog, FloatTrtFreq, FloatTrtArea,
           OtherTrtFreq, OtherTrtArea, OtherTrtAreaLog)

# function for predictions
pred_fun <- function(variable, model){
  
  pred_dat <- tibble(AreaOfInterestID = "A",
                     Time = min(rich_dat3$Time):max(rich_dat3$Time),
                     HydrPAC = mean(rich_dat_var$HydrPAC),
                     HydrTrtFreq = 0,
                     HydrTrtArea = 0,
                     HydrTrtAreaLog = log(0.01),
                     FloatPAC = mean(rich_dat_var$FloatPAC),
                     FloatPACLog = mean(rich_dat_var$FloatPACLog),
                     FloatTrtFreq = 0,
                     FloatTrtArea = 0,
                     OtherTrtFreq = 0,
                     OtherTrtArea = 0,
                     OtherTrtAreaLog = log(0.01)) %>%
    select(-{{variable}}) %>%
    expand_grid(rich_dat_var %>%
                  select({{variable}}))%>%
    mutate(Pred = predict(model, newdata = ., type = "response",
                          allow.new.levels = T),
           PredSE = predict(model, newdata = ., type = "response",
                          allow.new.levels = T, se.fit = T)$se.fit)
  
  return(pred_dat)
  
}

# hydr PAC
nat_pred_hydr_pac <- pred_fun(HydrPAC, target_mod3a)
ggplot(nat_pred_hydr_pac, aes(x = Time, y = Pred, 
                              color = HydrPAC,
                              group = HydrPAC)) +
  geom_line()

# floating PAC
nat_pred_float_pac <- pred_fun(FloatPACLog, target_mod3a) %>%
  mutate(FloatPAC = exp(FloatPACLog) - 0.01)
ggplot(nat_pred_float_pac, aes(x = Time, y = Pred, 
                              color = FloatPAC,
                              group = FloatPAC)) +
  geom_line()

# floating PAC trt freq
nat_pred_float_freq <- pred_fun(FloatTrtFreq, target_mod3a)
ggplot(nat_pred_float_freq, aes(x = Time, y = Pred, 
                               color = FloatTrtFreq,
                               group = FloatTrtFreq)) +
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

