#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom.mixed)
library(glmmTMB)
library(DHARMa)
library(GGally)
library(msm)

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
nat_mod1a <- update(nat_mod1, family = "nbinom2")
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
nat_res2a <- simulateResiduals(nat_mod2a, n = 1000)
plot(nat_res2a)
summary(nat_mod2a)

# fit model with non-linear relationships?
nat_mod3 <- glmmTMB(NativeRichness ~ Time + Time:(HydrPAC + HydrTrtFreq + FloatPAC + I(FloatPAC^2) + 
                                                    FloatTrtFreq + OtherTrtFreq + HydrTrtArea + I(HydrTrtArea^2) +
                                                    FloatTrtArea + OtherTrtArea + I(OtherTrtArea^2)) + 
                      (1|AreaOfInterestID),
                    data = rich_dat2,
                    family = poisson)
# can't converge


#### native richness coefficients ####

# extract coefficients and 95% CI
# transform to response scale by exponentiating and subtracting 1 (proportional change in richness when everything else is 0 and time is 1)
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
         PercChangeY1 = 100 * (exp(Estimate) - 1),
         PercChangeY20 = 100 * (exp(Estimate * 20) - 1),
         LowerPercChangeY1 = 100 * (exp(Lower) - 1),
         UpperPercChangeY1 = 100 * (exp(Upper) - 1))

# # extract estimates for delta methods
# nat_coef1 <- nat_est %>%
#   filter(Mod == "Freq") %>%
#   pull(Estimate)
# nat_coef2 <- nat_est %>%
#   filter(Mod == "Area") %>%
#   pull(Estimate)
# 
# # extract variance covariance matrix
# nat_vcov1 <- vcov(nat_mod1a)[[1]][2:12,2:12]
# nat_vcov2 <- vcov(nat_mod2a)[[1]][2:12,2:12]
# 
# # delta method using percent change transformation
# nat_sd1 <- deltamethod(list(~100 * (exp(x1) - 1), ~100 * (exp(x2) - 1), ~100 * (exp(x3) - 1), ~100 * (exp(x4) - 1),
#                             ~100 * (exp(x5) - 1), ~100 * (exp(x6) - 1), ~100 * (exp(x7) - 1), ~100 * (exp(x8) - 1),
#                             ~100 * (exp(x9) - 1), ~100 * (exp(x10) - 1), ~100 * (exp(x11) - 1)),
#                        nat_coef1, nat_vcov1)
# nat_sd2 <- deltamethod(list(~100 * (exp(x1) - 1), ~100 * (exp(x2) - 1), ~100 * (exp(x3) - 1), ~100 * (exp(x4) - 1),
#                             ~100 * (exp(x5) - 1), ~100 * (exp(x6) - 1), ~100 * (exp(x7) - 1), ~100 * (exp(x8) - 1),
#                             ~100 * (exp(x9) - 1), ~100 * (exp(x10) - 1), ~100 * (exp(x11) - 1)),
#                        nat_coef2, nat_vcov2)
# 
# # update etimates
# nat_est2 <- nat_est %>%
#   mutate(SDPercChange = c(nat_sd1, nat_sd2),
#          EstLowerPercChange = PercChangeY1 - SDPercChange * 1.26,
#          EstUpperPercChange = PercChangeY1 + SDPercChange * 1.26)
# # confidence intervals less conservative than transformed confint and not always consistent with p-value


#### visualize native richness model ####

# figure
nat_est %>%
  filter(Mod == "Freq" & (CoefType == "PAC" | is.na(CoefType))) %>%
  ggplot(aes(x = Level, y = PercChangeY1, color = Target)) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  geom_errorbar(aes(ymin = LowerPercChangeY1, ymax = UpperPercChangeY1), 
                width = 0.1, position = position_dodge(0.3)) +
  geom_point(size = 2, position = position_dodge(0.3)) +
  labs(x = "Abundance (PAC)", y = "Annual percent change in native richness") +
  def_theme_paper

nat_est %>%
  filter(Mod == "Area" & (CoefType == "PAC" | is.na(CoefType))) %>%
  ggplot(aes(x = Level, y = PercChangeY1, color = Target)) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  geom_errorbar(aes(ymin = LowerPercChangeY1, ymax = UpperPercChangeY1), 
                width = 0.1, position = position_dodge(0.3)) +
  geom_point(size = 2, position = position_dodge(0.3)) +
  labs(x = "Percentage of waterbody covered", y = "Annual percent change in native richness") +
  def_theme_paper
# pretty similar

nat_est %>%
  filter(Mod == "Freq" & (CoefType == "Trt" | is.na(CoefType))) %>%
  ggplot(aes(x = Level, y = PercChangeY1, color = Target)) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  geom_errorbar(aes(ymin = LowerPercChangeY1, ymax = UpperPercChangeY1), 
                width = 0.1, position = position_dodge(0.3)) +
  geom_point(size = 2, position = position_dodge(0.3)) +
  labs(x = "Percentage of years treated", 
       y = "Annual percent change in native richness") +
  def_theme_paper

nat_est %>%
  filter(Mod == "Area" & (CoefType == "Trt" | is.na(CoefType))) %>%
  ggplot(aes(x = Level, y = PercChangeY1, color = Target)) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  geom_errorbar(aes(ymin = LowerPercChangeY1, ymax = UpperPercChangeY1), 
                width = 0.1, position = position_dodge(0.3)) +
  geom_point(size = 2, position = position_dodge(0.3)) +
  labs(x = "Percentage of waterbody treated", 
       y = "Annual percent change in native richness") +
  def_theme_paper


#### to do - check estimates against predicted values ####