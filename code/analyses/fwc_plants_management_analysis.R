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
rich_dat <- read_csv("intermediate-data/FWC_plant_management_richness_analysis_formatted.csv")


#### examine richness data ####

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
  mutate(HydrPACF = fct_relevel(HydrPACF, "none", "low"),
         FloatPACF = fct_relevel(FloatPACF, "none", "low"),
         HydrTrtFreqF = fct_relevel(HydrTrtFreqF, "none", "low"),
         HydrTrtAreaF = fct_relevel(HydrTrtAreaF, "none", "low"),
         FloatTrtFreqF = fct_relevel(FloatTrtFreqF, "none", "low"),
         FloatTrtAreaF = fct_relevel(FloatTrtAreaF, "none", "low"),
         OtherTrtFreqF = fct_relevel(OtherTrtFreqF, "none", "low"),
         OtherTrtAreaF = fct_relevel(OtherTrtAreaF, "none", "low"),)


#### native richness model ####

# response variable
ggplot(rich_dat, aes(x = NativeRichness)) +
  geom_density()

# fit frequency model
nat_mod1 <- glmmTMB(NativeRichness ~ Time + Time*(HydrPACF + HydrTrtFreqF + FloatPACF + FloatTrtFreqF + 
                                             OtherTrtFreqF) + (1|AreaOfInterestID),
                    data = rich_dat2,
                    family = poisson)
nat_res1 <- simulateResiduals(nat_mod1, n = 1000)
plot(nat_res1)
summary(nat_mod1)
# floating PAC changes directions

# fit area model
nat_mod2 <- glmmTMB(NativeRichness ~ Time + Time*(HydrPACF + HydrTrtAreaF + FloatPACF + FloatTrtAreaF + 
                                                    OtherTrtAreaF) + (1|AreaOfInterestID),
                    data = rich_dat2,
                    family = poisson)
nat_res2 <- simulateResiduals(nat_mod2, n = 1000)
plot(nat_res2)
summary(nat_mod2)
# floating PAC changes directions
# hydr area not better with high than low
# other area changes directions

# transform variables with nonlinearities
rich_dat3 <- rich_dat2 %>%
  mutate(FloatPACLog = log(FloatPAC + 0.01),
         HydrTrtAreaLog = log(HydrTrtArea + 0.01),
         OtherTrtAreaLog = log(OtherTrtArea + 0.01))

# fit model with continuous variables
nat_mod3 <- glmmTMB(NativeRichness ~ Time + Time*(HydrPAC + FloatPACLog +
                                                    FloatTrtFreq + OtherTrtFreq + HydrTrtFreq + 
                                                    HydrTrtAreaLog + FloatTrtArea + OtherTrtAreaLog) +
                      (1|AreaOfInterestID),
                    data = rich_dat3,
                    family = "poisson",
                    control=glmmTMBControl(optimizer=optim,
                                           optArgs=list(method="BFGS")))
# okay to ignore function evaluation warnings if there are no convergence warnings
nat_res3 <- simulateResiduals(nat_mod3, n = 1000)
plot(nat_res3)
summary(nat_mod3)

# better fit with negative binomial?
nat_mod3a <- update(nat_mod3, family = "nbinom2")

# save model (slow to fit)
save(nat_mod3a, file = "output/native_richness_continuous_model.rda")
load("output/native_richness_continuous_model.rda")

# model diagnostics
nat_res3a <- simulateResiduals(nat_mod3a, n = 1000)
plot(nat_res3a)
summary(nat_mod3a)
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
nat_pred_hydr_pac <- pred_fun(HydrPAC, nat_mod3a)
ggplot(nat_pred_hydr_pac, aes(x = Time, y = Pred, 
                              color = HydrPAC,
                              group = HydrPAC)) +
  geom_line()

# floating PAC
nat_pred_float_pac <- pred_fun(FloatPACLog, nat_mod3a) %>%
  mutate(FloatPAC = exp(FloatPACLog) - 0.01)
ggplot(nat_pred_float_pac, aes(x = Time, y = Pred, 
                              color = FloatPAC,
                              group = FloatPAC)) +
  geom_line()

# floating PAC trt freq
nat_pred_float_freq <- pred_fun(FloatTrtFreq, nat_mod3a)
ggplot(nat_pred_float_freq, aes(x = Time, y = Pred, 
                               color = FloatTrtFreq,
                               group = FloatTrtFreq)) +
  geom_line()


#### native richness coefficients ####

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



