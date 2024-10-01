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
nat_mod1 <- glmmTMB(NativeRichness ~ Time + Time:(HydrPACF + HydrTrtFreqF + FloatPACF + FloatTrtFreqF + 
                                             OtherTrtFreqF) + (1|AreaOfInterestID),
                    data = rich_dat2,
                    family = poisson)
nat_res1 <- simulateResiduals(nat_mod1, n = 1000)
plot(nat_res1)
summary(nat_mod1)
# floating PAC changes directions

# better fit with negative binomial?
nat_mod1a <- update(nat_mod1, family = "nbinom2",
                    control=glmmTMBControl(optimizer=optim,
                                           optArgs=list(method="BFGS")))

# save model (slow to fit)
save(nat_mod1a, file = "output/native_richness_frequency_model.rda")
load("output/native_richness_frequency_model.rda")

# model diagnostics
nat_res1a <- simulateResiduals(nat_mod1a, n = 1000)
plot(nat_res1a)
summary(nat_mod1a)
# very large dispersion parameter

# fit area model
nat_mod2 <- glmmTMB(NativeRichness ~ Time + Time:(HydrPACF + HydrTrtAreaF + FloatPACF + FloatTrtAreaF + 
                                                    OtherTrtAreaF) + (1|AreaOfInterestID),
                    data = rich_dat2,
                    family = poisson)
nat_res2 <- simulateResiduals(nat_mod2, n = 1000)
plot(nat_res2)
summary(nat_mod2)
# hydr area not better with high than low
# other area changes directions

# better fit with negative binomial?
nat_mod2a <- update(nat_mod2, family = "nbinom2",
                    control=glmmTMBControl(optimizer=optim,
                                           optArgs=list(method="BFGS")))

# save model (slow to fit)
save(nat_mod2a, file = "output/native_richness_area_model.rda")
load("output/native_richness_area_model.rda")

# model diagnostics
nat_res2a <- simulateResiduals(nat_mod2a, n = 1000)
plot(nat_res2a)
summary(nat_mod2a)

# fit model with non-linear relationships?
# nat_mod3 <- glmmTMB(NativeRichness ~ Time + Time:(HydrPAC + HydrTrtFreq + FloatPAC + I(FloatPAC^2) + 
#                                                     FloatTrtFreq + OtherTrtFreq + HydrTrtArea + I(HydrTrtArea^2) +
#                                                     FloatTrtArea + OtherTrtArea + I(OtherTrtArea^2)) + 
#                       (1|AreaOfInterestID),
#                     data = rich_dat2,
#                     family = "nbinom2",
#                     control=glmmTMBControl(optimizer=optim,
#                                            optArgs=list(method="BFGS")))
# can't converge


#### native richness predicted values ####

# START HERE: try out function below

# function for predictions
pred_fun <- function(variable, model){
  
  pred_dat <- tibble(AreaOfInterestID = "A",
                     Time = min(rich_dat2$Time):max(rich_dat2$Time),
                     HydrPACF = "none",
                     HydrTrtFreqF = "none",
                     HydrTrtAreaF = "none",
                     FloatPACF = "none",
                     FloatTrtFreqF = "none",
                     FloatTrtAreaF = "none",
                     OtherTrtFreqF = "none",
                     OtherTrtAreaF = "none") %>%
    select(-{{variabe}}) %>%
    expand_grid({{variable}} = c("none", "low", "high"))%>%
    mutate(Pred = predict(model, newdata = ., type = "response",
                          allow.new.levels = T),
           PredSE = predict(model, newdata = ., type = "response",
                          allow.new.levels = T, se.fit = T)$se.fit)
  
  return(pred_dat)
  
  
}

# predicted richness for different levels of floating plant PAC
nat_pred_float_pac <- tibble(AreaOfInterestID = "A",
                             Time = 0:500,
                             HydrTrtFreqF = "none",
                             HydrPACF = "none",
                             FloatTrtFreqF = "none",
                             OtherTrtFreqF = "none") %>%
  expand_grid(FloatPACF = c("none", "low", "high")) %>%
  mutate(Pred = predict(nat_mod1a, newdata = ., type = "response",
                        allow.new.levels = T),
         PredLog = predict(nat_mod1a, newdata = ., allow.new.levels = T)) %>%
  left_join(filter(., FloatPACF == "none") %>% 
              select(Time, Pred) %>%
              rename(PredNone = Pred)) %>%
  mutate(PercChange = 100 * (Pred - PredNone) / PredNone, # % change relative to none
         Level = str_to_sentence(FloatPACF) %>%
           fct_relevel("None", "Low"))

# richness over time on log scale
ggplot(nat_pred_float_pac, aes(x = Time, y = PredLog, color = Level)) +
  geom_line()

# richness over time
ggplot(nat_pred_float_pac, aes(x = Time, y = Pred, color = Level)) +
  geom_line()
# not linear (exponential increase)

# change in richness relative to "none" over time
ggplot(nat_pred_float_pac, aes(x = Time, y = PercChange, color = Level)) +
  geom_line()
# percent change isn't stable over time


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



