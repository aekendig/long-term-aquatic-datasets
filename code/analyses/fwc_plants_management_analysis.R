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
library(patchwork)
library(khroma)

# figure settings
source("code/settings/figure_settings.R")

# import data
target_dat <- read_csv("intermediate-data/FWC_plant_management_target_analysis_formatted.csv")
target_taxa_dat <- read_csv("intermediate-data/FWC_plant_management_target_taxa_analysis_formatted.csv")
methods_dat <- read_csv("intermediate-data/FWC_plant_management_methods_analysis_formatted.csv")
methods_taxa_dat <- read_csv("intermediate-data/FWC_plant_management_methods_taxa_analysis_formatted.csv")
target_dat_new <- read_csv("intermediate-data/FWC_plant_management_new_target_analysis_formatted.csv")


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
                .fns = ~ fct_relevel(.x, "none", "low")),
         SurveyDay = yday(SurveyDate))


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

# quarter-specific frequency/intensity
methods_dat %>%
  select(c(AreaOfInterestID, ends_with(c("Q1", "Q2", "Q3", "Q4")))) %>%
  distinct() %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# correlations >= 0.4 and sig
# TrtFreqQ1 and Q2: 0.7 (visible)
# TrtFreqQ1 and Q3: 0.7 (visible)
# TrtFreqQ1 and Q4: 0.7 (visible)
# lots more correlated --> just use average month

# decide quarter to set as intercept in categorical model
methods_dat %>%
  distinct(AreaOfInterestID, TrtMonthF) %>%
  count(TrtMonthF)

# rearrange factor levels
methods_dat2 <- methods_dat %>%
  mutate(across(.cols = (!starts_with("TrtMonth") & ends_with("F")), 
                .fns = ~ fct_relevel(.x, "none", "low")),
         TrtMonthF = fct_relevel(TrtMonthF, "Q2", "Q3", "Q4", "Q1"))


#### survey date ####

# change in survey dates over time
day_mod1 <- glmmTMB(SurveyDay ~ Time + (1|AreaOfInterestID), data = target_dat2,
                    family = poisson)
day_res1 <- simulateResiduals(day_mod1, n = 1000)
plot(day_res1) # need to refit

day_mod2 <- update(day_mod1, family = "nbinom2")
day_res2 <- simulateResiduals(day_mod2, n = 1000)
plot(day_res2) # better

summary(day_mod2)

write_csv(tidy(day_mod2), "output/survey_date_model_summary.csv")


#### management efficacy models ####

# format PAC for beta distribution
target_dat_var2 <- target_dat_var %>%
  mutate(HydrProp = HydrPAC / 100,
         FloatProp = FloatPAC / 100)

# hydrilla binary model
hydr_mod1 <- betareg(HydrProp ~ HydrTrt, data = target_dat_var2)
plot(hydr_mod1)
summary(hydr_mod1)

# save
write_csv(tidy(hydr_mod1), "output/hydrilla_mgmt_efficacy_binary_model_summary.csv")

# # hydrilla categorical models
# hydr_mod2 <- betareg(HydrProp ~ HydrTrtFreqF, data = target_dat_var2)
# summary(hydr_mod2)
# # linear
# hydr_mod3 <- betareg(HydrProp ~ HydrTrtAreaF, data = target_dat_var2)
# summary(hydr_mod3)
# # linear

# hydrilla continuous model
hydr_mod2 <- betareg(HydrProp ~ HydrTrtFreq + HydrTrtArea, data = target_dat_var2)
plot(hydr_mod2)
summary(hydr_mod2)

# save
write_csv(tidy(hydr_mod2), "output/hydrilla_mgmt_efficacy_continuous_model_summary.csv")

# floating plant binary model
float_mod1 <- betareg(FloatProp ~ FloatTrt, data = target_dat_var2)
plot(float_mod1)
summary(float_mod1)

# save
write_csv(tidy(float_mod1), "output/floating_mgmt_efficacy_binary_model_summary.csv")

# # floating plant categorical models
# float_mod2 <- betareg(FloatProp ~ FloatTrtFreqF, data = target_dat_var2)
# summary(float_mod2)
# # linear
# float_mod3 <- betareg(FloatProp ~ FloatTrtAreaF, data = target_dat_var2)
# summary(float_mod3)
# # linear

# floating plant continuous model
float_mod2 <- betareg(FloatProp ~ FloatTrtFreq + FloatTrtArea, data = target_dat_var2)
plot(float_mod2)
summary(float_mod2)

# save
write_csv(tidy(float_mod2), "output/floating_mgmt_efficacy_continuous_model_summary.csv")


#### native richness model function ####

# combine model variables without time component
mod_dat_var <- methods_dat2 %>%
  distinct(AreaOfInterestID, TrtAreaCon, TrtAreaSys, TrtAreaNon,
           TrtFreqCon, TrtFreqSys, TrtFreqNon, 
           TrtMonth) %>%
  full_join(target_dat2 %>%
              distinct(AreaOfInterestID, HydrPAC, FloatPAC,
                       HydrTrtFreq, FloatTrtFreq, OtherTrtFreq,
                       HydrTrtArea, FloatTrtArea, OtherTrtArea))

# function for native plant model predictions
pred_fun <- function(variable, model, dat){
  
  pred_dat <- tibble(AreaOfInterestID = "A",
                     Time = max(dat$Time),
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
                  distinct({{variable}})) %>%
    mutate(across(.cols = -c(AreaOfInterestID, Time), 
                  .fns = list(Sq = ~.x^2, Log = ~log(.x + 0.01)),
                  .names = "{.col}{.fn}")) %>%
    mutate(Pred = predict(model, newdata = ., type = "response",
                          allow.new.levels = T),
           PredSE = predict(model, newdata = ., type = "response",
                            allow.new.levels = T, se.fit = T)$se.fit)
  
  return(pred_dat)
  
}


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
# Hydr PAC similar low and high
# floating PAC changes directions

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
# hydr area similar with low and high
# floating area similar with low and high
# other area positive then very low
# hydr pac now linear

# transform variables with nonlinearities
target_dat3 <- target_dat2 %>%
  mutate(HydrPACLog = log(HydrPAC + 0.01),
         FloatPACSq = FloatPAC^2,
         HydrTrtAreaLog = log(HydrTrtArea + 0.01),
         FloatTrtAreaLog = log(FloatTrtArea + 0.01),
         OtherTrtAreaSq = OtherTrtArea^2)

# check correlations
target_dat3 %>%
  distinct(AreaOfInterestID, HydrPACLog, FloatPAC, FloatPACSq,
             HydrTrtFreq, FloatTrtFreq, OtherTrtFreq, 
             HydrTrtAreaLog, FloatTrtAreaLog, 
             OtherTrtArea, OtherTrtAreaSq) %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# HydrTrtFreq and HydrPACLog: 0.6
# HydrTrtAreaLog and HydrPACLog: 0.7 (visible)
# HydrTrtAreaLog and HydrTrtFreq: 0.7
# OtherTrtFreq and FloatTrtFreq: 0.6
# FloatTrtAreaLog and FloatTrtFreq: 0.5

# take out transformation of hydrPAC
target_dat3 %>%
  distinct(AreaOfInterestID, HydrPAC, FloatPAC, FloatPACSq,
           HydrTrtFreq, FloatTrtFreq, OtherTrtFreq, 
           HydrTrtAreaLog, FloatTrtAreaLog, 
           OtherTrtArea, OtherTrtAreaSq) %>%
  select(-AreaOfInterestID) %>%
  ggpairs()
# no longer correlated

# fit model with continuous variables
target_mod3 <- glmmTMB(NativeRichness ~ Time + Time:(HydrPAC + FloatPAC + FloatPACSq +
                                                       HydrTrtFreq +FloatTrtFreq + OtherTrtFreq + 
                                                       HydrTrtAreaLog + FloatTrtAreaLog + 
                                                       OtherTrtArea + OtherTrtAreaSq) +
                      (1|AreaOfInterestID),
                    data = target_dat3,
                    family = "poisson",
                    control=glmmTMBControl(optimizer=optim,
                                           optArgs=list(method="BFGS")))

target_res3 <- simulateResiduals(target_mod3, n = 1000)
plot(target_res3)
summary(target_mod3)

# remove non-significant squared terms
target_mod3a <- update(target_mod3, .~. - Time:OtherTrtAreaSq)
summary(target_mod3a)

# evaluate non-linear terms
target_pred_FloatPAC <- pred_fun(FloatPAC, target_mod3a, target_dat3)
target_pred_FloatPAC %>%
  ggplot(aes(FloatPAC, Pred)) +
  geom_point()
target_dat %>%
  filter(FloatPAC > 7 & Time == max(Time)) %>%
  select(AreaOfInterestID, FloatPAC, NativeRichness)
# 2 waterbodies

target_pred_HydrTrtArea <- pred_fun(HydrTrtArea, target_mod3a, target_dat3)
target_pred_HydrTrtArea %>%
  ggplot(aes(HydrTrtArea, Pred)) +
  geom_point()
target_dat %>%
  filter(HydrTrtArea < 10 & Time == max(Time)) %>%
  select(AreaOfInterestID, HydrTrtArea, NativeRichness)
# a lot of points

target_pred_FloatTrtArea <- pred_fun(FloatTrtArea, target_mod3a, target_dat3)
target_pred_FloatTrtArea %>%
  ggplot(aes(FloatTrtArea, Pred)) +
  geom_point()
target_dat %>%
  filter(FloatTrtArea < 5 & Time == max(Time)) %>%
  select(AreaOfInterestID, FloatTrtArea, NativeRichness)
# a lot of points

# remove squared term
target_mod3b <- update(target_mod3a, .~. - Time:FloatPACSq)
summary(target_mod3b)

# change marginally sig log to linear?
target_mod3c <- update(target_mod3b, .~. - Time:FloatTrtAreaLog + Time:FloatTrtArea)
summary(target_mod3c)
# keep as log

# better fit with negative binomial?
target_mod3d <- update(target_mod3b, family = "nbinom2")
target_res3d <- simulateResiduals(target_mod3d, n = 1000)
plot(target_res3d)
summary(target_mod3d)
# almost identical to poisson

# save model (slow to fit)
save(target_mod3d, file = "output/native_richness_target_model.rda")
load("output/native_richness_target_model.rda")

# data for figures
flt_pac_fig_dat <- pred_fun(FloatPAC, target_mod3d, target_dat3)
hyd_frq_fig_dat <- pred_fun(HydrTrtFreq , target_mod3d, target_dat3)
flt_frq_fig_dat <- pred_fun(FloatTrtFreq, target_mod3d, target_dat3)
oth_frq_fig_dat <- pred_fun(OtherTrtFreq, target_mod3d, target_dat3)
hyd_ext_fig_dat <- pred_fun(HydrTrtArea, target_mod3d, target_dat3)
oth_ext_fig_dat <- pred_fun(OtherTrtArea, target_mod3d, target_dat3)
flt_ext_fig_dat <- pred_fun(FloatTrtArea, target_mod3d, target_dat3)

# figures
hyd_ext_fig <- ggplot(hyd_ext_fig_dat, aes(x = HydrTrtArea, y = Pred)) +
  geom_ribbon(aes(ymin = Pred - PredSE, ymax = Pred + PredSE), alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Hydrilla mgmt. extent",
       y = "Native plant richness") +
  def_theme_paper

flt_ext_fig <- ggplot(flt_ext_fig_dat, aes(x = FloatTrtArea, y = Pred)) +
  geom_ribbon(aes(ymin = Pred - PredSE, ymax = Pred + PredSE), alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Floating plant mgmt. extent",
       y = "Native plant richness") +
  def_theme_paper

# the log-transformed values have very little change in richness
target_mod3e <- glmmTMB(NativeRichness ~ Time + Time:(HydrPAC + FloatPAC +
                                                        HydrTrtFreq + FloatTrtFreq + OtherTrtFreq + 
                                                        HydrTrtArea + FloatTrtArea + OtherTrtArea) +
                          (1|AreaOfInterestID),
                        data = target_dat3,
                        family = "nbinom2",
                        control=glmmTMBControl(optimizer=optim,
                                               optArgs=list(method="BFGS")))
target_res3e <- simulateResiduals(target_mod3e, n = 1000)
plot(target_res3e)
summary(target_mod3e)

# save model (slow to fit)
save(target_mod3e, file = "output/native_richness_target_model.rda")
write_csv(tidy(target_mod3e), "output/native_richness_target_model_summary.csv")
load("output/native_richness_target_model.rda")

# data for figures
flt_pac_fig_dat <- pred_fun(FloatPAC, target_mod3e, target_dat3)
hyd_frq_fig_dat <- pred_fun(HydrTrtFreq , target_mod3e, target_dat3)
flt_frq_fig_dat <- pred_fun(FloatTrtFreq, target_mod3e, target_dat3)
oth_frq_fig_dat <- pred_fun(OtherTrtFreq, target_mod3e, target_dat3)
hyd_ext_fig_dat <- pred_fun(HydrTrtArea, target_mod3e, target_dat3)

# figures
flt_pac_fig <- ggplot(flt_pac_fig_dat, aes(x = FloatPAC, y = Pred)) +
  geom_ribbon(aes(ymin = Pred - PredSE, ymax = Pred + PredSE), alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Floating plant cover",
       y = "Native plant richness") +
  def_theme_paper

hyd_frq_fig <- ggplot(hyd_frq_fig_dat, aes(x = HydrTrtFreq, y = Pred)) +
  geom_ribbon(aes(ymin = Pred - PredSE, ymax = Pred + PredSE), alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Hydrilla mgmt. frequency",
       y = "Native plant richness") +
  def_theme_paper

flt_frq_fig <- ggplot(flt_frq_fig_dat, aes(x = FloatTrtFreq, y = Pred)) +
  geom_ribbon(aes(ymin = Pred - PredSE, ymax = Pred + PredSE), alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Floating plant mgmt. frequency",
       y = "Native plant richness") +
  def_theme_paper

oth_frq_fig <- ggplot(oth_frq_fig_dat, aes(x = OtherTrtFreq, y = Pred)) +
  geom_ribbon(aes(ymin = Pred - PredSE, ymax = Pred + PredSE), alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Other plant mgmt. frequency",
       y = "Native plant richness") +
  def_theme_paper

hyd_ext_fig <- ggplot(hyd_ext_fig_dat, aes(x = HydrTrtArea, y = Pred)) +
  geom_ribbon(aes(ymin = Pred - PredSE, ymax = Pred + PredSE), alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Hydrilla mgmt. extent",
       y = "Native plant richness") +
  def_theme_paper


#### taxon-specific target models ####

# select native taxa
# add columns
target_taxa_dat2 <- target_taxa_dat %>%
  filter(Origin == "Native") %>%
  mutate(HydrTrtAreaLog = log(HydrTrtArea + 0.01),
         FloatTrtAreaLog = log(FloatTrtArea + 0.01),
         Habitat = tolower(Habitat) %>%
           fct_recode("emergent" = "emersed"))

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
target_taxa_mod <- glmmTMB(Detected ~ Time + Time:(HydrPAC + FloatPAC +
                                               HydrTrtFreq +FloatTrtFreq + OtherTrtFreq + 
                                               HydrTrtArea + FloatTrtArea + OtherTrtArea) + 
                             (1|AreaOfInterestID),
                            data = target_taxa_sub, family = "binomial")

# use this to manually cycle through taxa (some models had warnings below)
# i <- i + 1
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
  target_taxa_mod <- glmmTMB(Detected ~ Time + Time:(HydrPAC + FloatPAC +
                                                       HydrTrtFreq +FloatTrtFreq + OtherTrtFreq + 
                                                       HydrTrtArea + FloatTrtArea + OtherTrtArea) + 
                               (1|AreaOfInterestID),
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
  
  target_taxa_coefs <- target_taxa_coefs %>%
    full_join(tidy(target_taxa_mod) %>%
                mutate(TaxonName = i) %>%
                relocate(TaxonName))
}

dev.off()

# save
write_csv(target_taxa_coefs, "output/target_model_taxa_coefficients.csv")
# import if needed
# target_taxa_coefs <- read_csv("output/target_model_taxa_coefficients.csv")

# correct p-values
# remove taxon that the estimates didn't work
target_taxa_coefs2 <- target_taxa_coefs %>%
  filter(str_detect(term, "Time:")) %>%
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

# significant interactions
target_taxa_coefs3 <- target_taxa_coefs2 %>%
  filter(q.value < 0.05) %>% 
  select(TaxonName, estimate, std.error, q.value, term) %>%
  left_join(target_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat) %>%
  mutate(Habitat = fct_relevel(Habitat, "floating"))

# significant coefficients table
target_taxa_tab <- target_taxa_coefs3 %>%
  mutate(across(.cols = c(estimate, std.error), 
                .fns = ~round_half_up(.x, 4) %>%
                  format(scientific = F)),
         q.value = round_half_up(q.value, 3),
         q.value = if_else(q.value == 0, "< 0.001", as.character(q.value)),
         response = paste0(estimate, " (", std.error, ")\n", q.value),
         term = str_remove(term, "Time:")) %>%
  select(Habitat, TaxonName, response, term) %>%
  pivot_wider(names_from = "term",
              values_from = "response") %>%
  arrange(Habitat, TaxonName) %>%
  relocate(c(HydrPAC, HydrTrtFreq, HydrTrtArea, 
             FloatPAC, FloatTrtFreq, FloatTrtArea, 
             OtherTrtFreq, OtherTrtArea), 
           .after = TaxonName) %>%
  mutate(across(.cols = c(HydrPAC, HydrTrtFreq, HydrTrtArea, 
                          FloatPAC, FloatTrtFreq, FloatTrtArea, 
                          OtherTrtFreq, OtherTrtArea),
                .fns = ~replace_na(.x, "")),
         TaxonName = str_replace_all(TaxonName, " ", "\n") %>%
           str_replace_all("\\/", "\\/\n"))

write_csv(target_taxa_tab, "output/target_model_taxa_interactions_table.csv")

# format estimates for figure
target_taxa_est <- target_taxa_coefs3 %>%
  mutate(term = fct_recode(term,
                           "hydrilla\nPAC" = "Time:HydrPAC",
                           "hydrilla\nmgmt.\nfrequency" = "Time:HydrTrtFreq",
                           "hydrilla\nmgmt.\nintensity" = "Time:HydrTrtArea",
                           "floating\nplant\nPAC" = "Time:FloatPAC",
                           "floating\nplant\nmgmt.\nfrequency" = "Time:FloatTrtFreq",
                           "floating\nplant\nmgmt.\nintensity" = "Time:FloatTrtArea",
                           "other\nplant\nmgmt.\nfrequency" = "Time:OtherTrtFreq",
                           "other\nplant\nmgmt.\nintensity" = "Time:OtherTrtArea"),
         sign = if_else(estimate > 0, "pos", "neg")) %>%
  count(term, sign, Habitat) %>%
  complete(term, sign, Habitat) %>%
  mutate(n = replace_na(n, 0),
         n = if_else(sign == "neg", -1 * n, n),
         term = fct_relevel(term, "hydrilla\nPAC",
                            "hydrilla\nmgmt.\nfrequency",
                            "hydrilla\nmgmt.\nintensity",
                            "floating\nplant\nPAC",
                            "floating\nplant\nmgmt.\nfrequency",
                            "floating\nplant\nmgmt.\nintensity",
                            "other\nplant\nmgmt.\nfrequency",
                            "other\nplant\nmgmt.\nintensity"))

# figure
taxa_target_fig <- ggplot(target_taxa_est,
                          aes(x = term, y = n, fill = Habitat)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  scale_fill_manual(values = colour("muted")(7)[c(7, 5, 3)]) +
  labs(x = "Covariate interaction with time",
       y = "Number of native taxa\n(in direction of response)") +
  def_theme_paper +
  theme(legend.position = "inside",
        legend.position.inside = c(0.5, 0.78))

# combine figures
target_fig1 <- hyd_frq_fig + hyd_ext_fig + flt_pac_fig
target_fig2 <- plot_spacer() + flt_frq_fig + oth_frq_fig + plot_spacer() +
  plot_layout(widths = c(1/4, 1, 1, 1/4))
target_fig <- target_fig1 / target_fig2 / taxa_target_fig +
  plot_layout(nrow = 3, heights = c(0.6, 0.6, 1)) &
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

ggsave("output/target_figure.png", target_fig,
       width = 6.5, height = 8)


#### native richness methods model ####

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
         TrtAreaSysSq = TrtAreaSys^2,
         TrtMonthStd = TrtMonth - trt_month_avg,
         TrtMonthStdSq = TrtMonthStd^2)

# check correlations
methods_dat3 %>%
  distinct(AreaOfInterestID, TrtAreaCon, TrtAreaSys, TrtAreaSysSq, TrtAreaNon,
           TrtFreqConLog, TrtFreqSys, TrtFreqNon, TrtFreqNonSq,
           TrtMonthStd, TrtMonthStdSq) %>%
  select(-AreaOfInterestID) %>%
  ggpairs()

# fit continuous model
method_mod3 <- glmmTMB(NativeRichness ~ Time + Time:(TrtAreaCon + TrtAreaSys + TrtAreaSysSq + TrtAreaNon +
                                                       TrtFreqConLog + TrtFreqSys + TrtFreqNon + TrtFreqNonSq +
                                                       TrtMonthStd + TrtMonthStdSq) + (1|AreaOfInterestID),
                       data = methods_dat3,
                       family = poisson,
                       control=glmmTMBControl(optimizer=optim,
                                              optArgs=list(method="BFGS")))
method_res3 <- simulateResiduals(method_mod3, n = 1000)
plot(method_res3)
summary(method_mod3)

# remove non-significant squares
method_mod3a <- update(method_mod3, .~. -Time:TrtFreqNonSq - Time:TrtMonthStdSq)
method_res3a <- simulateResiduals(method_mod3a, n = 1000)
plot(method_res3a)
summary(method_mod3a)

# evaluate non-linear effect of systemic area
method_pred_TrtAreaSys <- pred_fun(TrtAreaSys, method_mod3a, methods_dat3)
method_pred_TrtAreaSys %>%
  ggplot(aes(TrtAreaSys, Pred)) +
  geom_point()
methods_dat %>%
  filter(TrtAreaSys > 50 & Time == max(Time)) %>%
  select(AreaOfInterestID, TrtAreaSys, NativeRichness)
# 3 points

# evaluate non-linear effect of contact frequency
method_pred_TrtFreqCon <- pred_fun(TrtFreqCon, method_mod3a, methods_dat3)
method_pred_TrtFreqCon %>%
  ggplot(aes(TrtFreqCon, Pred)) +
  geom_point()
methods_dat %>%
  filter(TrtFreqCon < 25) %>%
  group_by(TrtFreqCon & Time == max(Time)) %>%
  summarize(n = n(),
            NativeRichness = mean(NativeRichness))
# there are a lot of points, but it's non significant

# remove squared term
method_mod3b <- update(method_mod3a, .~. -Time:TrtAreaSysSq)
method_res3b <- simulateResiduals(method_mod3b, n = 1000)
plot(method_res3b)
summary(method_mod3b)

# remove log-transformation
method_mod3c <- update(method_mod3b, .~. -Time:TrtFreqConLog + Time:TrtFreqCon)
method_res3c <- simulateResiduals(method_mod3c, n = 1000)
plot(method_res3c)
summary(method_mod3c)

# negative binomial response dist.
method_mod3d <- update(method_mod3c, family = "nbinom2") 
method_res3d <- simulateResiduals(method_mod3d, n = 1000)
plot(method_res3d)
summary(method_mod3d)
# similar diagnostics and estimates

# save model (slow to fit)
save(method_mod3d, file = "output/native_richness_method_model.rda")
write_csv(tidy(method_mod3d), "output/native_richness_method_model_summary.csv")
load("output/native_richness_method_model.rda")

# data for figures
cont_ext_fig_dat <- pred_fun(TrtAreaCon, method_mod3d, methods_dat3)
sys_ext_fig_dat <- pred_fun(TrtAreaSys, method_mod3d, methods_dat3)
non_freq_fig_dat <- pred_fun(TrtFreqNon, method_mod3d, methods_dat3)

# figures
cont_ext_fig <- ggplot(cont_ext_fig_dat, aes(x = TrtAreaCon, y = Pred)) +
  geom_ribbon(aes(ymin = Pred - PredSE, ymax = Pred + PredSE), alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Contact herbicide extent",
       y = "Native plant richness") +
  def_theme_paper

sys_ext_fig <- ggplot(sys_ext_fig_dat, aes(x = TrtAreaSys, y = Pred)) +
  geom_ribbon(aes(ymin = Pred - PredSE, ymax = Pred + PredSE), alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Systemic herbicide extent",
       y = "Native plant richness") +
  def_theme_paper

non_freq_fig <- ggplot(non_freq_fig_dat, aes(x = TrtFreqNon, y = Pred)) +
  geom_ribbon(aes(ymin = Pred - PredSE, ymax = Pred + PredSE), alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Non-herbicide frequency",
       y = "Native plant richness") +
  def_theme_paper


#### taxon-specific methods models ####

# select native taxa
# standardize treatment month using average from above
methods_taxa_dat2 <- methods_taxa_dat %>%
  filter(Origin == "Native") %>%
  mutate(TrtFreqConLog = log(TrtFreqCon + 0.01),
         TrtMonthStd = TrtMonth - trt_month_avg,
         Habitat = tolower(Habitat) %>%
           fct_recode("emergent" = "emersed"))

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
methods_taxa_mod <- glmmTMB(Detected ~ Time + Time:(TrtAreaCon + TrtAreaSys + TrtAreaNon +
                                                      TrtFreqCon + TrtFreqSys + TrtFreqNon +
                                                      TrtMonthStd) + (1|AreaOfInterestID),
                            data = methods_taxa_sub, family = "binomial")

# use this to manually cycle through taxa (some models had warnings below)
# i <- i + 1
# taxa with model convergence issues:
# 43: Najas marina

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
  mutate(TaxonName = methods_taxa[1]) %>%
  relocate(TaxonName)

# loop through taxa
pdf("output/taxa_specific_methods_models.pdf")

for(i in methods_taxa[-43]) {
  
  # subset data
  methods_taxa_sub <- filter(methods_taxa_dat2, TaxonName == i)
  
  # fit model
  methods_taxa_mod <- glmmTMB(Detected ~ Time + Time:(TrtAreaCon + TrtAreaSys + TrtAreaNon +
                                                        TrtFreqCon + TrtFreqSys + TrtFreqNon +
                                                        TrtMonthStd) + (1|AreaOfInterestID),
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
                mutate(TaxonName = i) %>%
                relocate(TaxonName))
}

dev.off()

# save
write_csv(methods_taxa_coefs, "output/methods_model_taxa_coefficients.csv")
# import if needed
# methods_taxa_coefs <- read_csv("output/methods_model_taxa_coefficients.csv")

# correct p-values for interaction terms
methods_taxa_coefs2 <- methods_taxa_coefs %>%
  filter(str_detect(term, "Time:")) %>%
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

# significant interactions
methods_taxa_coefs3 <- methods_taxa_coefs2 %>%
  filter(q.value < 0.05) %>% 
  select(TaxonName, estimate, std.error, q.value, term) %>%
  left_join(methods_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat) %>%
  mutate(Habitat = fct_relevel(Habitat, "floating"))

# significant coefficients table
methods_taxa_tab <- methods_taxa_coefs3 %>%
  mutate(across(.cols = c(estimate, std.error, q.value), .fns = ~round_half_up(.x, 3)),
         q.value = if_else(q.value == 0, "< 0.001", as.character(q.value)),
         response = paste0(estimate, " (", std.error, ")\n", q.value),
         term = str_remove(term, "Time:Trt")) %>%
  select(Habitat, TaxonName, response, term) %>%
  pivot_wider(names_from = "term",
              values_from = "response") %>%
  arrange(Habitat, TaxonName) %>%
  relocate(c(FreqCon, AreaCon, FreqSys, AreaSys, FreqNon, AreaNon, MonthStd), .after = TaxonName) %>%
  mutate(across(.cols = c(FreqCon, AreaCon, FreqSys, AreaSys, FreqNon, AreaNon, MonthStd),
                .fns = ~replace_na(.x, "")),
         TaxonName = str_replace_all(TaxonName, " ", "\n") %>%
           str_replace_all("\\/", "\\/\n"))

write_csv(methods_taxa_tab, "output/methods_model_taxa_interactions_table.csv")

# format estimates for figure
method_taxa_est <- methods_taxa_coefs3 %>%
  mutate(term = fct_recode(term,
                           "contact\nherbicide\nextent" = "Time:TrtAreaCon",
                           "systemic\nherbicide\nextent" = "Time:TrtAreaSys",
                           "non-\nherbicide\nextent" = "Time:TrtAreaNon",
                           "contact\nherbicide\nfrequency" = "Time:TrtFreqCon",
                           "systemic\nherbicide\nfrequency" = "Time:TrtFreqSys",
                           "non-\nherbicide\nfrequency" = "Time:TrtFreqNon",
                           "average\nmanagement\nmonth" = "Time:TrtMonthStd"),
         sign = if_else(estimate > 0, "pos", "neg")) %>%
  count(term, sign, Habitat) %>%
  complete(term, sign, Habitat) %>%
  mutate(n = replace_na(n, 0),
         n = if_else(sign == "neg", -1 * n, n),
         term = fct_relevel(term, 
                            "contact\nherbicide\nfrequency",
                            "contact\nherbicide\nextent",
                            "systemic\nherbicide\nfrequency",
                            "systemic\nherbicide\nextent",
                            "non-\nherbicide\nfrequency",
                            "non-\nherbicide\nextent",
                            "average\nmanagement\nmonth"))

# figure
taxa_method_fig <- ggplot(method_taxa_est,
                          aes(x = term, y = n, fill = Habitat)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  scale_fill_manual(values = colour("muted")(7)[c(7, 5, 3)]) +
  # scale_fill_brewer(type = "qual", palette = "Dark2",
  #                   name = "Growth form") +
  labs(x = "Covariate interaction with time",
       y = "Number of native taxa\n(in direction of response)") +
  def_theme_paper +
  theme(legend.position = "inside",
        legend.position.inside = c(0.24, 0.78))

# combine figures
method_fig <- (cont_ext_fig + sys_ext_fig + non_freq_fig) /
  taxa_method_fig +
  plot_layout(heights = c(0.6, 1)) &
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

ggsave("output/method_figure.png", method_fig,
       width = 18, height = 16, units = "cm", dpi = 600)


# #### native richness predicted values ####
# 
# # contact herb intensity
# pred_cont_ints <- pred_fun(TrtAreaCon, method_mod3a, methods_dat3)
# method_fig1 <- ggplot(pred_cont_ints, aes(x = Time, y = Pred, 
#                               color = TrtAreaCon,
#                               group = TrtAreaCon)) +
#   geom_line() +
#   scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
#                         midpoint = mean(c(min(mod_dat_var$TrtAreaCon, na.rm = T),
#                                           max(mod_dat_var$TrtAreaCon, na.rm = T))),
#                         name = "Contact\nherbicide\nintensity") +
#   labs(x = "Years of monitoring", y = "Native taxonomic richness") +
#   def_theme_paper
# 
# # systemic herb intensity
# pred_sys_ints <- pred_fun(TrtAreaSys, method_mod3a, methods_dat3)
# method_fig2 <- ggplot(pred_sys_ints, aes(x = Time, y = Pred, 
#                                            color = TrtAreaSys,
#                                            group = TrtAreaSys)) +
#   geom_line() +
#   scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
#                         midpoint = mean(c(min(mod_dat_var$TrtAreaSys, na.rm = T),
#                                           max(mod_dat_var$TrtAreaSys, na.rm = T))),
#                         name = "Systemic\nherbicide\nintensity") +
#   labs(x = "Years of monitoring", y = "Native taxonomic richness") +
#   def_theme_paper
# 
# # non-herb frequency
# pred_non_freq <- pred_fun(TrtFreqNon, method_mod3a, methods_dat3)
# method_fig3 <- ggplot(pred_non_freq, aes(x = Time, y = Pred, 
#                                           color = TrtFreqNon,
#                                           group = TrtFreqNon)) +
#   geom_line() +
#   scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
#                         midpoint = mean(c(min(mod_dat_var$TrtFreqNon, na.rm = T),
#                                           max(mod_dat_var$TrtFreqNon, na.rm = T))),
#                         name = "Non-\nherbicide\nfrequency") +
#   labs(x = "Years of monitoring", y = "Native taxonomic richness") +
#   def_theme_paper
# 
# # hydr trt frequency
# pred_hydr_freq <- pred_fun(HydrTrtFreq, target_mod3b, target_dat3)
# target_fig1 <- ggplot(pred_hydr_freq, aes(x = Time, y = Pred, 
#                                           color = HydrTrtFreq,
#                                           group = HydrTrtFreq)) +
#   geom_line() +
#   scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
#                         midpoint = mean(c(min(mod_dat_var$HydrTrtFreq, na.rm = T),
#                                           max(mod_dat_var$HydrTrtFreq, na.rm = T))),
#                         name = "Hydrilla\nmanagement\nfrequency") +
#   labs(x = "Years of monitoring", y = "Native taxonomic richness") +
#   def_theme_paper
# 
# # hydr trt intensity
# pred_hydr_ints <- pred_fun(HydrTrtArea, target_mod3b, target_dat3)
# target_fig2 <- ggplot(pred_hydr_ints, aes(x = Time, y = Pred, 
#                                           color = HydrTrtArea,
#                                           group = HydrTrtArea)) +
#   geom_line() +
#   scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
#                         midpoint = mean(c(min(mod_dat_var$HydrTrtArea, na.rm = T),
#                                           max(mod_dat_var$HydrTrtArea, na.rm = T))),
#                         name = "Hydrilla\nmanagement\nintensity") +
#   labs(x = "Years of monitoring", y = "Native taxonomic richness") +
#   def_theme_paper
# 
# # floating PAC
# pred_float_pac <- pred_fun(FloatPAC, target_mod3b, target_dat3)
# target_fig3 <- ggplot(pred_float_pac, aes(x = Time, y = Pred, 
#                                          color = FloatPAC,
#                                          group = FloatPAC)) +
#   geom_line() +
#   scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
#                         midpoint = mean(c(min(mod_dat_var$FloatPAC, na.rm = T),
#                                           max(mod_dat_var$FloatPAC, na.rm = T))),
#                         name = "Floating\nplant\ncover") +
#   labs(x = "Years of monitoring", y = "Native taxonomic richness") +
#   def_theme_paper
# 
# # floating trt freq
# pred_float_freq <- pred_fun(FloatTrtFreq, target_mod3b, target_dat3)
# target_fig4 <- ggplot(pred_float_freq, aes(x = Time, y = Pred, 
#                                color = FloatTrtFreq,
#                                group = FloatTrtFreq)) +
#   geom_line() +
#   scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
#                         midpoint = mean(c(min(mod_dat_var$FloatTrtFreq, na.rm = T),
#                                           max(mod_dat_var$FloatTrtFreq, na.rm = T))),
#                         name = "Floating\nmanagement\nfrequency") +
#   labs(x = "Years of monitoring", y = "Native taxonomic richness") +
#   def_theme_paper
# 
# # other trt freq
# pred_other_freq <- pred_fun(OtherTrtFreq, target_mod3b, target_dat3)
# target_fig5 <- ggplot(pred_other_freq, aes(x = Time, y = Pred, 
#                                 color = OtherTrtFreq,
#                                 group = OtherTrtFreq)) +
#   geom_line() +
#   geom_line() +
#   scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
#                         midpoint = mean(c(min(mod_dat_var$OtherTrtFreq, na.rm = T),
#                                           max(mod_dat_var$OtherTrtFreq, na.rm = T))),
#                         name = "Other\nmanagement\nfrequency") +
#   labs(x = "Years of monitoring", y = "Native taxonomic richness") +
#   def_theme_paper
# 
# 
# #### native richness coefficient plots ####
# 
# # systemic intensity
# change_sys_ints <- confint(method_mod3a, parm = c("Time", "Time:TrtAreaSys", "Time:TrtAreaSysSq")) %>%
#    cbind(tibble(Term = c("c", "b", "a"))) %>%
#    select(Term, Estimate) %>%
#    pivot_wider(names_from = Term, values_from = Estimate) %>%
#    expand_grid(tibble(Intensity = sort(unique(methods_dat3$TrtAreaSys)))) %>%
#    mutate(Change = a * Intensity^2 + b * Intensity + c,
#           ExpChange = exp(Change))
# 
# change_sys_ints_fig <- ggplot(change_sys_ints, aes(x = Intensity, y = ExpChange)) +
#   geom_point() +
#   scale_y_continuous(breaks = 1) +
#   labs(x = "Intensity", y = "Rate of change") + 
#   def_theme_paper +
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
#         axis.title = element_text(size = 8))
# 
# # check understanding of meaning
# pred_sys_ints %>%
#   select(Time, TrtAreaSys, Pred) %>%
#   full_join(pred_sys_ints %>%
#               filter(Time == 0) %>%
#               select(TrtAreaSys, Pred) %>%
#               rename(Pred0 = Pred)) %>%
#   mutate(RelChange = Pred/Pred0) %>%
#   ggplot(aes(x = TrtAreaSys, y = RelChange)) +
#   geom_point() +
#   facet_wrap(~ Time, scales = "free")
# 
# # floating plant PAC
# change_float_pac <- confint(target_mod3b, parm = c("Time", "Time:FloatPAC", "Time:FloatPACSq")) %>%
#   cbind(tibble(Term = c("c", "b", "a"))) %>%
#   select(Term, Estimate) %>%
#   pivot_wider(names_from = Term, values_from = Estimate) %>%
#   expand_grid(tibble(PAC = sort(unique(target_dat3$FloatPAC)))) %>%
#   mutate(Change = a * PAC^2 + b * PAC + c,
#          ExpChange = exp(Change))
# 
# change_float_pac_fig <- ggplot(change_float_pac, aes(x = PAC, y = ExpChange)) +
#   geom_point() +
#   scale_y_continuous(breaks = 1) +
#   labs(x = "Cover", y = "Rate of change") + 
#   def_theme_paper +
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
#         axis.title = element_text(size = 8))
# 
# # check understanding of meaning
# pred_float_pac %>%
#   select(Time, FloatPAC, Pred) %>%
#   full_join(pred_float_pac %>%
#               filter(Time == 0) %>%
#               select(FloatPAC, Pred) %>%
#               rename(Pred0 = Pred)) %>%
#   mutate(RelChange = Pred/Pred0) %>%
#   ggplot(aes(x = FloatPAC, y = RelChange)) +
#   geom_point() +
#   facet_wrap(~ Time, scales = "free")
# 
# 
# #### combine figures ####
# 
# method_fig <- method_fig1 + 
#   method_fig2 + inset_element(change_sys_ints_fig, left = 0.01, bottom = 0.01, right = 0.4, top = 0.45) +
#   method_fig3 +
#   plot_layout(ncol = 1) +
#   plot_annotation(tag_levels = list(c("A", "B", "", "C")))
# 
# ggsave("output/methods_figure.png", method_fig,
#        width = 4, height = 9)
# 
# target_fig <- target_fig1 + target_fig2 + target_fig3 + inset_element(change_float_pac_fig, left = 0.01, bottom = 0.01, right = 0.38, top = 0.42) + target_fig4 + target_fig5 +
#   plot_layout(ncol = 2) +
#   plot_annotation(tag_levels = list(c("A", "B", "C", "", "D", "E")))
# 
# ggsave("output/target_figure.png", target_fig,
#        width = 7.5, height = 9)
# 
# 
#### variations of target model ####

# remove waterbodies with no management
target_dat3_mgmt <- target_dat3 %>%
  filter(!(HydrTrtFreq == 0 & FloatTrtFreq == 0 & OtherTrtFreq == 0))

# compare number of waterbodies
n_distinct(target_dat3_mgmt$AreaOfInterestID)
n_distinct(target_dat3$AreaOfInterestID)

# fit model
target_mod3_mgmt <- update(target_mod3e, data = target_dat3_mgmt)
summary(target_mod3_mgmt)

# save model (slow to fit)
save(target_mod3_mgmt, file = "output/native_richness_target_model_mgmt_only.rda")
write_csv(tidy(target_mod3_mgmt), "output/target_model_all_mgmt_summary.csv")

# target dataset from "new" management period (used for methods model)
# remove waterbodies with no management
target_dat_new2 <- target_dat_new %>%
  filter(!(HydrTrtFreq == 0 & FloatTrtFreq == 0 & OtherTrtFreq == 0))

# compare number of waterbodies
n_distinct(target_dat_new2$AreaOfInterestID)
n_distinct(methods_dat3$AreaOfInterestID)
# same

# compare sample size
nrow(target_dat_new2)
nrow(methods_dat3)
# same

# fit model
target_mod3_new <- update(target_mod3e, data = target_dat_new2) # convergence error
summary(target_mod3_new) # very large dispersion parameter -- > small overdispersion
target_mod3_new2 <- update(target_mod3e, data = target_dat_new2, 
                           family = "poisson")
summary(target_mod3_new2) # this seems to change the estimates substantially
target_mod3_new3 <- update(target_mod3_new, # try default optimizer
                           control=glmmTMBControl(optimizer=nlminb,
                                                  optArgs=list())) # convergence error

# poisson version of full dataset model
target_mod3f <- update(target_mod3e, family = "poisson")

# model summaries
summary(target_mod3_new2)
summary(target_mod3f)

# save model
save(target_mod3_new2, file = "output/native_richness_target_model_newer_data.rda")
write_csv(tidy(target_mod3_new2) %>%
            mutate(dataset = "detailed methods") %>%
            full_join(tidy(target_mod3f) %>%
                        mutate(dataset = "full")), "output/target_model_new_mgmt_summary.csv")
