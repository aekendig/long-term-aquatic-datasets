#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom.mixed)
library(glmmTMB)
library(DHARMa)
library(GGally)

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


#### native richness model ####

# response variable
ggplot(rich_dat, aes(x = NativeRichness)) +
  geom_density()

# fit model
nat_mod1 <- glmmTMB(NativeRichness ~ Time + Time:(HydrPAC + HydrTrtFreq + HydrTrtArea + FloatPAC + FloatTrtFreq + 
                                             FloatTrtArea + OtherTrtFreq + OtherTrtArea) + (1|AreaOfInterestID),
                    data = rich_dat,
                    family = poisson)
# warnings - these are okay when there's no convergence warning
nat_res1 <- simulateResiduals(nat_mod1, n = 1000)
plot(nat_res1) # aggregation of points around one prediction value (intercept?)
summary(nat_mod1)

# try different distribution for residual patterns
nat_mod2 <- update(nat_mod1, family = nbinom2)
# didn't converge

# try different distribution for residual patterns
nat_mod3 <- update(nat_mod1, family = nbinom1)
# didn't converge

# average values for prediction
rich_cov <- rich_dat %>%
  distinct(AreaOfInterestID, HydrPAC, FloatPAC, HydrTrtFreq, HydrTrtArea, FloatTrtFreq, FloatTrtArea,
           OtherTrtFreq, OtherTrtArea) %>%
  select(-AreaOfInterestID) %>%
  summarize(across(.cols = everything(), .fn = mean))

# frequency predictions
nat_hydr_freq_pred <- rich_cov %>%
  select(-HydrTrtFreq) %>%
  expand_grid(tibble(HydrTrtFreq = c(0, 25, 50, 100))) %>%
  expand_grid(rich_dat %>%
                distinct(Time)) %>%
  mutate(AreaOfInterestID = "A") %>%
  mutate(NativeRichness = predict(nat_mod1, newdata = ., type = "response", allow.new.levels = T)) %>%
  group_by() %>%
  mutate(BaselineRichness = min(NativeRichness)) %>%
  ungroup() %>%
  mutate(Frequency = as.factor(HydrTrtFreq),
         RichnessChange = NativeRichness - BaselineRichness)

# visualize
ggplot(nat_hydr_freq_pred, aes(x = Time, y = NativeRichness, color = Frequency)) +
  geom_line()

ggplot(nat_hydr_freq_pred, aes(x = Time, y = RichnessChange, color = Frequency)) +
  geom_line()

rich_dat %>%
  filter(HydrTrtFreq == 0) %>%
  ggplot(aes(x = Time, y = NativeRichness, color = as.character(AreaOfInterestID))) +
  geom_smooth(formula = 'y ~ x', method = "glm", se = F) +
  theme(legend.position = "none")

rich_dat %>%
  filter(HydrTrtFreq == 100) %>%
  full_join(rich_dat %>%
              filter(HydrTrtFreq == 100 & Time == 0) %>%
              select(AreaOfInterestID, NativeRichness) %>%
              rename(BaselineRichness = NativeRichness)) %>%
  mutate(RichnessChange = NativeRichness - BaselineRichness) %>%
  ggplot(aes(x = Time, y = RichnessChange, color = as.character(AreaOfInterestID))) +
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = "glm", se = F) +
  theme(legend.position = "none")
