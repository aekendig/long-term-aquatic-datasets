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


#### survey date ####

# change in survey dates over time
day_mod1 <- glmmTMB(SurveyDay ~ Time + (1|AreaOfInterestID), data = target_dat2,
                    family = poisson)
day_res1 <- simulateResiduals(day_mod1, n = 1000)
plot(day_res1) # need to refit

day_mod2 <- update(day_mod1, family = "nbinom1")
day_res2 <- simulateResiduals(day_mod2, n = 1000)
plot(day_res2) # better

day_mod3 <- update(day_mod1, family = "nbinom2")
day_res3 <- simulateResiduals(day_mod3, n = 1000)
plot(day_res3) # same as above

summary(day_mod2)

# visualize
day_pred <- tibble(Time = min(target_dat2$Time):max(target_dat2$Time),
                   AreaOfInterestID = "A") %>%
  mutate(PredSurveyDay = predict(day_mod2, newdata = ., type = "response",
                                 allow.new.levels = T),
         PredSurveyDaySE = predict(day_mod2, newdata = ., type = "response",
                                   allow.new.levels = T, se.fit = T)$se.fit)

ggplot(day_pred, aes(x = Time, y = PredSurveyDay)) +
  geom_ribbon(aes(ymin = PredSurveyDay - PredSurveyDaySE,
                  ymax = PredSurveyDay + PredSurveyDaySE),
              alpha = 0.5) +
  geom_line() +
  def_theme_paper


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

# better fit with negative binomial?
method_mod3b <- update(method_mod3a, family = "nbinom1") 
summary(method_mod3b)
# convergence problem with both nbinom1 and 2
# note dispersion parameter explanation in family_glmmTMB
# couldn't fit with method_mod3 either
# dispersion parameter is very small (1*10^-7), which could cause it


#### taxon-specific methods models ####

# select native taxa
# standardize treatment month using average from above
methods_taxa_dat2 <- methods_taxa_dat %>%
  filter(Origin == "Native") %>%
  mutate(TrtFreqConLog = log(TrtFreqCon + 0.01),
         TrtAreaSysSq = TrtAreaSys^2,
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
methods_taxa_mod <- glmmTMB(Detected ~ Time + Time:(TrtAreaCon + TrtAreaSys + TrtAreaSysSq + TrtAreaNon +
                                                      TrtFreqConLog + TrtFreqSys + TrtFreqNon +
                                                      TrtMonthStd) + (1|AreaOfInterestID),
                                  data = methods_taxa_sub, family = "binomial")

# use this to manually cycle through taxa (some models had warnings below)
i <- i + 1
# taxa with model convergence issues:
# 22: Fontinalis spp.
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

for(i in methods_taxa[-c(22, 43)]) {
  
  # subset data
  methods_taxa_sub <- filter(methods_taxa_dat2, TaxonName == i)
  
  # fit model
  methods_taxa_mod <- glmmTMB(Detected ~ Time + Time:(TrtAreaCon + TrtAreaSys + TrtAreaSysSq + TrtAreaNon +
                                                        TrtFreqConLog + TrtFreqSys + TrtFreqNon +
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

# correct p-values
# remove taxon that the estimates didn't work
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

# contact herbicide intensity
taxa_cont_ints <- methods_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:TrtAreaCon") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(methods_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat, .after = 1)

write_csv(taxa_cont_ints, "output/methods_model_taxa_contact_intensity_interaction_table.csv")

# systemic herbicide intensity
taxa_sys_ints <- methods_taxa_coefs2 %>%
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
  left_join(methods_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
  arrange(TaxonName)

write_csv(taxa_sys_ints , "output/methods_model_taxa_systemic_intensity_interaction_table.csv")

# non-herbicide frequency
taxa_non_freq <- methods_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:TrtFreqNon") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(methods_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat, .after = 1)

write_csv(taxa_non_freq, "output/methods_model_taxa_nonherbicide_frequency_interaction_table.csv")

# format estimates
# squared: choose squared estimate if significant, linear if not; make long
# add estimates that don't have tables
method_taxa_est <- taxa_cont_ints %>%
  mutate(term = "contact\nherbicide\nintensity") %>%
  full_join(taxa_sys_ints %>%
              mutate(intensity_estimate = if_else(intensity2_q.value < 0.05, NA_real_,
                                                  intensity_estimate),
                     intensity2_estimate = if_else(intensity2_q.value < 0.05, intensity2_estimate,
                                                   NA_real_)) %>%
              pivot_longer(cols = starts_with("intensity"),
                          names_to = c("term", ".value"),
                          names_sep = "_") %>%
              mutate(term = fct_recode(term,
                                       "systemic\nherbicide\nintensity" = "intensity",
                                       "squared\nsystemic\nherbicide\nintensity" = "intensity2")) %>%
              filter(!is.na(estimate))) %>%
  full_join(taxa_non_freq %>%
              mutate(term = "non-\nherbicide\nfrequency")) %>%
  full_join(methods_taxa_coefs2 %>%
              filter(q.value < 0.05 & term %in% c("Time:TrtAreaNon",
                                                  "Time:TrtFreqConLog",
                                                  "Time:TrtFreqSys",
                                                  "Time:TrtMonthStd")) %>% 
              select(TaxonName, estimate, std.error, q.value, term) %>%
              left_join(methods_taxa_dat2 %>%
                          distinct(TaxonName, Habitat)) %>%
              relocate(Habitat, .after = 1) %>%
              mutate(term = fct_recode(term,
                                       "non-\nherbicide\nintensity" = "Time:TrtAreaNon",
                                       "contact\nherbicide\nfrequency" = "Time:TrtFreqConLog",
                                       "systemic\nherbicide\nfrequency" = "Time:TrtFreqSys",
                                       "average\nmanagement\nmonth" = "Time:TrtMonthStd"))) %>%
  mutate(sign = if_else(estimate > 0, "pos", "neg")) %>%
  count(term, sign, Habitat) %>%
  mutate(Habitat = fct_relevel(Habitat, "floating"),
         n = if_else(sign == "neg", -1 * n, n),
         term = fct_relevel(term, 
                            "contact\nherbicide\nfrequency",
                            "contact\nherbicide\nintensity",
                            "systemic\nherbicide\nfrequency",
                            "systemic\nherbicide\nintensity",
                            "squared\nsystemic\nherbicide\nintensity",
                            "non-\nherbicide\nfrequency",
                            "non-\nherbicide\nintensity",
                            "average\nmanagement\nmonth"))

taxa_method_fig <- ggplot(method_taxa_est,
       aes(x = term, y = n, fill = Habitat)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  scale_fill_brewer(type = "qual", palette = "Dark2",
                    name = "Growth form") +
  labs(x = "Covariate interaction with time",
       y = "Number of taxa (in direction of response)") +
  def_theme_paper +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.2))

ggsave("output/taxa_method_figure.png", taxa_method_fig,
       width = 6, height = 4)


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

# better fit with negative binomial?
target_mod3a <- update(target_mod3, family = "nbinom2")

# model diagnostics
target_res3a <- simulateResiduals(target_mod3a, n = 1000)
plot(target_res3a)
summary(target_mod3a)
# very small dispersion parameter
# results are very similar

# remove squared terms that aren't needed
target_mod3b <- update(target_mod3a, .~. -Time:OtherTrtAreaSq)

# save model (slow to fit)
save(target_mod3b, file = "output/native_richness_target_model.rda")
load("output/native_richness_target_model.rda")

# model diagnostics
target_res3b <- simulateResiduals(target_mod3b, n = 1000)
plot(target_res3b)
summary(target_mod3b)


#### taxon-specific target models ####

# select native taxa
# add columns
target_taxa_dat2 <- target_taxa_dat %>%
  filter(Origin == "Native") %>%
  mutate(FloatPACSq = FloatPAC^2,
         HydrTrtAreaLog = log(HydrTrtArea + 0.01),
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
target_taxa_mod <- glmmTMB(Detected ~ Time + Time:(HydrPAC + FloatPAC + FloatPACSq +
                                               HydrTrtFreq +FloatTrtFreq + OtherTrtFreq + 
                                               HydrTrtAreaLog + FloatTrtAreaLog + OtherTrtArea) + 
                             (1|AreaOfInterestID),
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
                                                       HydrTrtAreaLog + FloatTrtAreaLog + OtherTrtArea) + 
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

# hydrilla treatment frequency
taxa_hydr_freq <- target_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:HydrTrtFreq") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(target_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat, .after = 1)

write_csv(taxa_hydr_freq, "output/target_model_taxa_hydrilla_frequency_interaction_table.csv")

# hydrilla treatment intensity
taxa_hydr_ints <- target_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:HydrTrtAreaLog") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(target_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat, .after = 1)

write_csv(taxa_hydr_ints, "output/target_model_taxa_hydrilla_intensity_interaction_table.csv")

# floating PAC
taxa_float_pac <- target_taxa_coefs2 %>%
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
  left_join(target_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat, .after = 1) %>%
  arrange(TaxonName)

write_csv(taxa_float_pac, "output/target_model_taxa_floating_PAC_interaction_table.csv")

# floating treatment frequency
taxa_float_freq <- target_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:FloatTrtFreq") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(target_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat, .after = 1)

write_csv(taxa_float_freq, "output/target_model_taxa_floating_frequency_interaction_table.csv")

# other treatment frequency
taxa_other_freq <- target_taxa_coefs2 %>%
  filter(q.value < 0.05 & term == "Time:OtherTrtFreq") %>% 
  select(TaxonName, estimate, std.error, q.value) %>%
  left_join(target_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat, .after = 1)

write_csv(taxa_other_freq, "output/target_model_taxa_other_frequency_interaction_table.csv")

# format estimates
# squared: choose squared estimate if significant, linear if not; make long
# add estimates that don't have tables
target_taxa_est <- taxa_hydr_freq %>%
  mutate(term = "hydrilla\nmgmt.\nfrequency") %>%
  full_join(taxa_hydr_ints %>%
              mutate(term = "hydrilla\nmgmt.\nintensity")) %>%
  full_join(taxa_float_pac %>%
              mutate(pac_estimate = if_else(pac2_q.value < 0.05, NA_real_,
                                                  pac_estimate),
                     pac2_estimate = if_else(pac2_q.value < 0.05, pac2_estimate,
                                                   NA_real_)) %>%
              pivot_longer(cols = starts_with("pac"),
                           names_to = c("term", ".value"),
                           names_sep = "_") %>%
              mutate(term = fct_recode(term,
                                       "floating\nplant\nPAC" = "pac",
                                       "squared\nfloating\nplant\nPAC" = "pac2")) %>%
              filter(!is.na(estimate))) %>%
  full_join(taxa_float_freq %>%
              mutate(term = "floating\nplant\nmgmt.\nfrequency")) %>%
    full_join(taxa_other_freq %>%
                mutate(term = "other\nplant\nmgmt.\nfrequency")) %>%
  full_join(target_taxa_coefs2 %>%
              filter(q.value < 0.05 & term %in% c("Time:HydrPAC",
                                                  "Time:FloatTrtAreaLog",
                                                  "Time:OtherTrtArea")) %>% 
              select(TaxonName, estimate, std.error, q.value, term) %>%
              left_join(target_taxa_dat2 %>%
                          distinct(TaxonName, Habitat)) %>%
              relocate(Habitat, .after = 1) %>%
              mutate(term = fct_recode(term,
                                       "hydrilla\nPAC" = "Time:HydrPAC",
                                       "floating\nplant\nmgmt.\nintensity" = "Time:FloatTrtAreaLog",
                                       "other\nplant\nmgmt.\nintensity" = "Time:OtherTrtArea"))) %>%
  mutate(sign = if_else(estimate > 0, "pos", "neg")) %>%
  count(term, sign, Habitat) %>%
  mutate(Habitat = fct_relevel(Habitat, "floating"),
         n = if_else(sign == "neg", -1 * n, n),
         term = fct_relevel(term, "hydrilla\nPAC",
                            "hydrilla\nmgmt.\nfrequency",
                            "hydrilla\nmgmt.\nintensity",
                            "floating\nplant\nPAC",
                            "squared\nfloating\nplant\nPAC",
                            "floating\nplant\nmgmt.\nfrequency",
                            "floating\nplant\nmgmt.\nintensity",
                            "other\nplant\nmgmt.\nfrequency",
                            "other\nplant\nmgmt.\nintensity"))

taxa_target_fig <- ggplot(target_taxa_est,
                          aes(x = term, y = n, fill = Habitat)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  scale_fill_brewer(type = "qual", palette = "Dark2", name = "Growth form") +
  labs(x = "Covariate interaction with time",
       y = "Number of taxa (in direction of response)") +
  def_theme_paper +
  theme(legend.position = "inside",
        legend.position.inside = c(0.12, 0.85))

ggsave("output/taxa_target_figure.png", taxa_target_fig,
       width = 6, height = 4)


#### native richness predicted values ####

# remove time component
mod_dat_var <- methods_dat2 %>%
  distinct(AreaOfInterestID, TrtAreaCon, TrtAreaSys, TrtAreaNon,
           TrtFreqCon, TrtFreqSys, TrtFreqNon, 
           TrtMonth) %>%
  full_join(target_dat2 %>%
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

# contact herb intensity
pred_cont_ints <- pred_fun(TrtAreaCon, method_mod3a, methods_dat3)
method_fig1 <- ggplot(pred_cont_ints, aes(x = Time, y = Pred, 
                              color = TrtAreaCon,
                              group = TrtAreaCon)) +
  geom_line() +
  scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
                        midpoint = mean(c(min(mod_dat_var$TrtAreaCon, na.rm = T),
                                          max(mod_dat_var$TrtAreaCon, na.rm = T))),
                        name = "Contact\nherbicide\nintensity") +
  labs(x = "Years of monitoring", y = "Native taxonomic richness") +
  def_theme_paper

# systemic herb intensity
pred_sys_ints <- pred_fun(TrtAreaSys, method_mod3a, methods_dat3)
method_fig2 <- ggplot(pred_sys_ints, aes(x = Time, y = Pred, 
                                           color = TrtAreaSys,
                                           group = TrtAreaSys)) +
  geom_line() +
  scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
                        midpoint = mean(c(min(mod_dat_var$TrtAreaSys, na.rm = T),
                                          max(mod_dat_var$TrtAreaSys, na.rm = T))),
                        name = "Systemic\nherbicide\nintensity") +
  labs(x = "Years of monitoring", y = "Native taxonomic richness") +
  def_theme_paper

# non-herb frequency
pred_non_freq <- pred_fun(TrtFreqNon, method_mod3a, methods_dat3)
method_fig3 <- ggplot(pred_non_freq, aes(x = Time, y = Pred, 
                                          color = TrtFreqNon,
                                          group = TrtFreqNon)) +
  geom_line() +
  scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
                        midpoint = mean(c(min(mod_dat_var$TrtFreqNon, na.rm = T),
                                          max(mod_dat_var$TrtFreqNon, na.rm = T))),
                        name = "Non-\nherbicide\nfrequency") +
  labs(x = "Years of monitoring", y = "Native taxonomic richness") +
  def_theme_paper

# hydr trt frequency
pred_hydr_freq <- pred_fun(HydrTrtFreq, target_mod3b, target_dat3)
target_fig1 <- ggplot(pred_hydr_freq, aes(x = Time, y = Pred, 
                                          color = HydrTrtFreq,
                                          group = HydrTrtFreq)) +
  geom_line() +
  scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
                        midpoint = mean(c(min(mod_dat_var$HydrTrtFreq, na.rm = T),
                                          max(mod_dat_var$HydrTrtFreq, na.rm = T))),
                        name = "Hydrilla\nmanagement\nfrequency") +
  labs(x = "Years of monitoring", y = "Native taxonomic richness") +
  def_theme_paper

# hydr trt intensity
pred_hydr_ints <- pred_fun(HydrTrtArea, target_mod3b, target_dat3)
target_fig2 <- ggplot(pred_hydr_ints, aes(x = Time, y = Pred, 
                                          color = HydrTrtArea,
                                          group = HydrTrtArea)) +
  geom_line() +
  scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
                        midpoint = mean(c(min(mod_dat_var$HydrTrtArea, na.rm = T),
                                          max(mod_dat_var$HydrTrtArea, na.rm = T))),
                        name = "Hydrilla\nmanagement\nintensity") +
  labs(x = "Years of monitoring", y = "Native taxonomic richness") +
  def_theme_paper

# floating PAC
pred_float_pac <- pred_fun(FloatPAC, target_mod3b, target_dat3)
target_fig3 <- ggplot(pred_float_pac, aes(x = Time, y = Pred, 
                                         color = FloatPAC,
                                         group = FloatPAC)) +
  geom_line() +
  scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
                        midpoint = mean(c(min(mod_dat_var$FloatPAC, na.rm = T),
                                          max(mod_dat_var$FloatPAC, na.rm = T))),
                        name = "Floating\nplant\ncover") +
  labs(x = "Years of monitoring", y = "Native taxonomic richness") +
  def_theme_paper

# floating trt freq
pred_float_freq <- pred_fun(FloatTrtFreq, target_mod3b, target_dat3)
target_fig4 <- ggplot(pred_float_freq, aes(x = Time, y = Pred, 
                               color = FloatTrtFreq,
                               group = FloatTrtFreq)) +
  geom_line() +
  scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
                        midpoint = mean(c(min(mod_dat_var$FloatTrtFreq, na.rm = T),
                                          max(mod_dat_var$FloatTrtFreq, na.rm = T))),
                        name = "Floating\nmanagement\nfrequency") +
  labs(x = "Years of monitoring", y = "Native taxonomic richness") +
  def_theme_paper

# other trt freq
pred_other_freq <- pred_fun(OtherTrtFreq, target_mod3b, target_dat3)
target_fig5 <- ggplot(pred_other_freq, aes(x = Time, y = Pred, 
                                color = OtherTrtFreq,
                                group = OtherTrtFreq)) +
  geom_line() +
  geom_line() +
  scale_color_gradient2(low = "#1b9e77", mid = "#d95f02", high = "#7570b3",
                        midpoint = mean(c(min(mod_dat_var$OtherTrtFreq, na.rm = T),
                                          max(mod_dat_var$OtherTrtFreq, na.rm = T))),
                        name = "Other\nmanagement\nfrequency") +
  labs(x = "Years of monitoring", y = "Native taxonomic richness") +
  def_theme_paper


#### native richness coefficient plots ####

# systemic intensity
change_sys_ints <- confint(method_mod3a, parm = c("Time", "Time:TrtAreaSys", "Time:TrtAreaSysSq")) %>%
   cbind(tibble(Term = c("c", "b", "a"))) %>%
   select(Term, Estimate) %>%
   pivot_wider(names_from = Term, values_from = Estimate) %>%
   expand_grid(tibble(Intensity = sort(unique(methods_dat3$TrtAreaSys)))) %>%
   mutate(Change = a * Intensity^2 + b * Intensity + c,
          ExpChange = exp(Change))

change_sys_ints_fig <- ggplot(change_sys_ints, aes(x = Intensity, y = ExpChange)) +
  geom_point() +
  scale_y_continuous(breaks = 1) +
  labs(x = "Intensity", y = "Rate of change") + 
  def_theme_paper +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title = element_text(size = 8))

# check understanding of meaning
pred_sys_ints %>%
  select(Time, TrtAreaSys, Pred) %>%
  full_join(pred_sys_ints %>%
              filter(Time == 0) %>%
              select(TrtAreaSys, Pred) %>%
              rename(Pred0 = Pred)) %>%
  mutate(RelChange = Pred/Pred0) %>%
  ggplot(aes(x = TrtAreaSys, y = RelChange)) +
  geom_point() +
  facet_wrap(~ Time, scales = "free")

# floating plant PAC
change_float_pac <- confint(target_mod3b, parm = c("Time", "Time:FloatPAC", "Time:FloatPACSq")) %>%
  cbind(tibble(Term = c("c", "b", "a"))) %>%
  select(Term, Estimate) %>%
  pivot_wider(names_from = Term, values_from = Estimate) %>%
  expand_grid(tibble(PAC = sort(unique(target_dat3$FloatPAC)))) %>%
  mutate(Change = a * PAC^2 + b * PAC + c,
         ExpChange = exp(Change))

change_float_pac_fig <- ggplot(change_float_pac, aes(x = PAC, y = ExpChange)) +
  geom_point() +
  scale_y_continuous(breaks = 1) +
  labs(x = "Cover", y = "Rate of change") + 
  def_theme_paper +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title = element_text(size = 8))

# check understanding of meaning
pred_float_pac %>%
  select(Time, FloatPAC, Pred) %>%
  full_join(pred_float_pac %>%
              filter(Time == 0) %>%
              select(FloatPAC, Pred) %>%
              rename(Pred0 = Pred)) %>%
  mutate(RelChange = Pred/Pred0) %>%
  ggplot(aes(x = FloatPAC, y = RelChange)) +
  geom_point() +
  facet_wrap(~ Time, scales = "free")


#### combine figures ####

method_fig <- method_fig1 + 
  method_fig2 + inset_element(change_sys_ints_fig, left = 0.01, bottom = 0.01, right = 0.4, top = 0.45) +
  method_fig3 +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = list(c("A", "B", "", "C")))

ggsave("output/methods_figure.png", method_fig,
       width = 4, height = 9)

target_fig <- target_fig1 + target_fig2 + target_fig3 + inset_element(change_float_pac_fig, left = 0.01, bottom = 0.01, right = 0.38, top = 0.42) + target_fig4 + target_fig5 +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "", "D", "E")))

ggsave("output/target_figure.png", target_fig,
       width = 7.5, height = 9)


#### variations of target model ####

# remove waterbodies with no management
target_dat3_mgmt <- target_dat3 %>%
  filter(!(HydrTrtFreq == 0 & FloatTrtFreq == 0 & OtherTrtFreq == 0))

# compare number of waterbodies
n_distinct(target_dat3_mgmt$AreaOfInterestID)
n_distinct(target_dat3$AreaOfInterestID)

# fit model
target_mod3_mgmt <- glmmTMB(NativeRichness ~ Time + Time:(HydrPAC + FloatPAC + FloatPACSq +
                                                       HydrTrtFreq +FloatTrtFreq + OtherTrtFreq + 
                                                       HydrTrtAreaLog + FloatTrtAreaLog + 
                                                       OtherTrtArea + OtherTrtAreaSq) +
                         (1|AreaOfInterestID),
                       data = target_dat3_mgmt,
                       family = "poisson",
                       control=glmmTMBControl(optimizer=optim,
                                              optArgs=list(method="BFGS")))
summary(target_mod3_mgmt)

# refit with negative binomial
# target_mod3a_mgmt <- update(target_mod3_mgmt, family = "nbinom2")
# couldn't converge

# remove squared terms that aren't needed
target_mod3b_mgmt <- update(target_mod3_mgmt, .~. -Time:OtherTrtAreaSq)

# compare estimates
summary(target_mod3b_mgmt)
summary(target_mod3b)
# same direction, similar magnitude

# target dataset from "new" management period (used for methods model)
# remove waterbodies with no management
target_dat_new2 <- target_dat_new %>%
  filter(!(HydrTrtFreq == 0 & FloatTrtFreq == 0 & OtherTrtFreq == 0)) %>%
  mutate(HydrPACLog = log(HydrPAC + 0.01),
         FloatPACSq = FloatPAC^2,
         HydrTrtAreaLog = log(HydrTrtArea + 0.01),
         FloatTrtAreaLog = log(FloatTrtArea + 0.01),
         OtherTrtAreaSq = OtherTrtArea^2)

# compare number of waterbodies
n_distinct(target_dat_new2$AreaOfInterestID)
n_distinct(methods_dat3$AreaOfInterestID)
# same

# compare sample size
nrow(target_dat_new2)
nrow(methods_dat3)
# same

# fit model
target_mod3_new <- glmmTMB(NativeRichness ~ Time + Time:(HydrPAC + FloatPAC + FloatPACSq +
                                                            HydrTrtFreq +FloatTrtFreq + OtherTrtFreq + 
                                                            HydrTrtAreaLog + FloatTrtAreaLog + 
                                                            OtherTrtArea + OtherTrtAreaSq) +
                              (1|AreaOfInterestID),
                            data = target_dat_new2,
                            family = "poisson",
                            control=glmmTMBControl(optimizer=optim,
                                                   optArgs=list(method="BFGS")))

# refit with negative binomial
# target_mod3a_new <- update(target_mod3_new, family = "nbinom2")
# couldn't converge

# squared term significance are different - keep all in for comparison
summary(target_mod3_new)
summary(target_mod3b)

# comparison with same terms
target_mod3b_new <- update(target_mod3_new, .~. -Time:OtherTrtAreaSq)
summary(target_mod3b_new)
summary(target_mod3b)

# export tables
tidy(target_mod3b_mgmt) %>%
  write_csv("output/target_model_all_mgmt_summary.csv")

tidy(target_mod3b_new) %>%
  write_csv("output/target_model_new_mgmt_summary.csv")


#### list of native taxa ####

# both datasets should be 83 species
methods_taxa_dat2 %>%
  distinct(TaxonName, Habitat) %>%
  full_join(target_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  mutate(Habitat = fct_relevel(Habitat, "submersed")) %>%
  arrange(Habitat, TaxonName) %>%
  write_csv("output/native_taxa_table.csv")
