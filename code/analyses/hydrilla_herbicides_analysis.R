#### info ####

# goal: evaluate effects of herbicides on hydrilla
# model structure: https://nwfsc-timeseries.github.io/atsa-labs/sec-marss-fitting-with-stan.html


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)
library(reshape2) # for melt
library(cowplot)
library(nlme)
library(piecewiseSEM)
library(lme4)
library(car)
# library(MARSS) # for marss model
# library(zoo) # for uneven time series
# library(ggfortify)

# stan settings
source("code/stan_settings.R")

# import data
ctrl_old <- read_csv("intermediate-data/FWC_control_old_formatted.csv")
ctrl_new <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
plant_fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")
dayDat <- read_csv("intermediate-data/LakeO_day_data_for_model.csv")

# load models
load("output/lakeO_intraannual_veg_change_mod.rda")

# assumptions
MinTime = 5 # minimum number of time points needed per waterbody

# figure settings
def_theme <- theme_bw() +
  theme(axis.text = element_text(size = 8, color="black"),
        axis.title = element_text(size = 10, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.text = element_text(size = 10, color="black"),
        strip.background = element_blank())


#### lake O parameters ####

# extract coefficients
lakeO_beta1 <- coef(summary(lake0_mod))[2, "Estimate"] # days
lakeO_beta2 <- coef(summary(lake0_mod))[3, "Estimate"] # days^2

# date of max abundance (-b/2a)
lakeO_days <- -lakeO_beta1 / (2 * lakeO_beta2)


#### edit FWC plant data ####

# Perm IDs per AOI
plant_fwc %>%
  group_by(AreaOfInterestID) %>%
  summarise(IDs = length(unique(PermanentID))) %>%
  filter(IDs > 1)
# one Perm ID per AOI

# function to remove duplicates
rem_dups_fun <- function(dat) {
  
  if (nrow(dat) == 1){
    dat2 <- dat
  } else if (var(dat$AreaCovered_ha) > 0){
    dat2 <- dat %>%
      filter(AreaCovered_ha == max(AreaCovered_ha)) # choose survey with maximum area
  } else if (var(dat$AreaCovered_ha) == 0){
    dat2 <- dat %>%
      mutate(DaysDiff = abs(lakeO_days - Days)) %>%
      filter(DaysDiff == min(DaysDiff)) %>% # choose survey closes to max. abundance time
      select(-DaysDiff)
  }
  
  return(dat2)
}

plant_fwc_hyd <- plant_fwc %>% # start with all surveys (no hydrilla means abundance = 0)
  filter(SpeciesName != "Hydrilla verticillata") %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate) %>%
  unique() %>% # don't need a row for each species
  left_join(plant_fwc %>% # add hydrilla information
              filter(SpeciesName == "Hydrilla verticillata") %>%
              select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, SpeciesAcres)) %>%
  mutate(SpeciesAcres = replace_na(SpeciesAcres, 0), # hydrilla cover 0 when it wasn't in a survey
         Area_ha = ShapeArea * 100, # convert lake area from km-squared to hectares
         AreaCovered_ha = SpeciesAcres * 0.405, # convert plant cover from acres to hectares
         AreaCovered_ha = case_when(AreaCovered_ha > Area_ha ~ Area_ha,
                                    TRUE ~ AreaCovered_ha), # make plant cover the size of the lake area if it exceeds it
         SurveyMonth = month(SurveyDate),
         SurveyDay = day(SurveyDate),
         SurveyYear = case_when(SurveyMonth >= 4 ~ year(SurveyDate),
                                SurveyMonth < 4 ~ year(SurveyDate) - 1), # assume growing season starts in April
         MonthDay = case_when(SurveyMonth >= 4 ~ as.Date(paste("2020", SurveyMonth, SurveyDay, sep = "-")), # start "year" in April (2020/2021 are arbitrary)
                              SurveyMonth < 4 ~ as.Date(paste("2021", SurveyMonth, SurveyDay, sep = "-")))) %>%
  left_join(dayDat) %>% # add standardized days (proportion between April 1 and March 31)
  mutate(AreaChangeSD = lakeO_beta1 * (lakeO_days-Days) + lakeO_beta2 * (lakeO_days^2 - Days^2)) %>% # calculate the number of sd's to change to get est. max abundance
  group_by(AreaOfInterestID) %>%
  mutate(EstAreaCoveredRaw_ha = AreaCovered_ha + AreaChangeSD * sd(AreaCovered_ha), # calculate est. max abundance
         EstAreaCovered_ha = case_when(EstAreaCoveredRaw_ha > Area_ha ~ Area_ha, # reduce areas covered to total area
                                       TRUE ~ EstAreaCoveredRaw_ha),
         log_EstAreaCovered_ha = log(EstAreaCovered_ha + 1e-10),  # log-transform (min EstAreaCovered_ha (> 0) = 1x10^-9)
         PropCovered = EstAreaCovered_ha / Area_ha,
         log_PropCovered = log((EstAreaCovered_ha + 1e-10) / Area_ha),
         Present = as.numeric(AreaCovered_ha > 0)) %>%
  ungroup() %>%
  nest(-c(AreaOfInterestID, SurveyYear)) %>% # find duplicates within year and waterbody
  mutate(newdata = map(data, ~rem_dups_fun(.))) %>% # remove duplicates
  select(-data) %>% # removes 41 rows of data
  unnest(newdata)

# duplicates in same year
plant_fwc_hyd %>%
  group_by(AreaOfInterestID, SurveyYear) %>%
  count() %>%
  filter(n > 1)

# save data
write_csv(plant_fwc_hyd, "intermediate-data/FWC_hydrilla_survey_formatted.csv")


#### edit control data ####

# non-herbicide methods (from herbicide_initial_visualizations)
non_herb <- c("Mechanical Harvester", 
              "Snagging (tree removal)", 
              "Aquatic Dye (for shading)", 
              "Grass Carp", "Hand Removal", 
              "Mechanical (Other)", 
              "Mechanical Shredder", 
              "Prescribed Fire")

# duplicates per year/date?
ctrl_old %>%
  filter(Species == "Hydrilla verticillata" & TotalAcres > 0) %>%
  group_by(AreaOfInterestID, Year) %>%
  summarise(apps = n()) %>%
  filter(apps > 1)
# yes, 34
# different acres in same year, probably different events

ctrl_new %>%
  filter(Species == "Hydrilla verticillata" & TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
  group_by(AreaOfInterestID, BeginDate) %>%
  summarise(apps = n()) %>%
  filter(apps > 1)
# yes, 1059
ctrl_new %>%
  filter(AreaOfInterestID == 8 & BeginDate == "2019-05-20")
# for treatments involving more than one herbicide, separate lines for each herbicide type, but they share a treatment ID and acres

# old herbicide data
ctrl_old_hyd <- ctrl_old %>%
  filter(Species == "Hydrilla verticillata" & TotalAcres > 0) %>%
  mutate(AreaTreated_ha= TotalAcres * 0.405,
         Area_ha = ShapeArea * 100,
         TreatmentMonth = 12,
         TreatmentDate = as.Date(paste(as.character(Year), "-12-01", sep = "")), # no date given for these surveys
         TreatmentID = paste("old", seq(1:length(Area_ha)), sep = "_")) %>%
  select(AreaOfInterestID, PermanentID, Year, Area_ha, AreaTreated_ha, TreatmentMonth, TreatmentDate, TreatmentID) %>%
  rename(TreatmentYear = Year)

# new herbicide data
ctrl_new_hyd <- ctrl_new %>%
  filter(Species == "Hydrilla verticillata" & TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>% # herbicide control only
  select(AreaOfInterestID, PermanentID, BeginDate, TreatmentID, TotalAcres, ShapeArea) %>%
  unique() %>% # captures area treated for an event without duplication due to multiple herbicides
  mutate(AreaTreated_ha = TotalAcres * 0.405,
         Area_ha = ShapeArea * 100,
         TreatmentYear = year(BeginDate),
         TreatmentMonth = month(BeginDate),
         TreatmentID = as.character(TreatmentID)) %>%
  select(AreaOfInterestID, PermanentID, TreatmentYear, Area_ha, AreaTreated_ha, TreatmentMonth, BeginDate, TreatmentID) %>%
  rename(TreatmentDate = BeginDate)

# combine herbicide data
ctrl_hyd <- ctrl_old_hyd %>%
  full_join(ctrl_new_hyd) %>%
  mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha, # make treatment area size of lake if it exceeds it
                                    TRUE ~ AreaTreated_ha))

# add zero for missing years (assume no application)
ctrl_hyd2 <- ctrl_old %>%
  mutate(Area_ha = ShapeArea * 100) %>%
  select(AreaOfInterestID, PermanentID, Area_ha) %>%
  unique() %>%
  full_join(ctrl_new %>%
              mutate(Area_ha = ShapeArea * 100) %>%
              select(AreaOfInterestID, PermanentID, Area_ha) %>%
              unique()) %>%
  expand_grid(tibble(TreatmentYear = min(ctrl_old_hyd$TreatmentYear):max(ctrl_new_hyd$TreatmentYear))) %>% # one row per year
  full_join(ctrl_hyd) %>%
  mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0), # make area treated zero if the year wasn't included
         TreatmentEvent = TreatmentID, # NA when no treatment applied
         TreatmentID = case_when(is.na(TreatmentID) ~ paste("none", TreatmentYear, sep = "_"),
                                 TRUE ~ TreatmentID))

# save data
write_csv(ctrl_hyd2, "intermediate-data/FWC_hydrilla_herbicide_formatted.csv")


#### long-term linear trends by lake ####

# find trend and herbicide usage for each lake
trend_hyd <- plant_fwc_hyd %>%
  group_by(AreaOfInterestID, PermanentID) %>%
  summarise(Trend = as.numeric(coef(lm(log_EstAreaCovered_ha ~ SurveyYear))[2])) %>%
  ungroup() %>%
  left_join(ctrl_hyd2 %>%
              group_by(AreaOfInterestID, PermanentID) %>%
              summarise(TreatEvents = length(unique(na.omit(TreatmentEvent))),
                        TreatFreq = TreatEvents/length(unique(TreatmentYear)),
                        TreatInt = mean(AreaTreated_ha/Area_ha)))

# correlation
cor.test(~ TreatFreq + TreatInt, data = trend_hyd)
# marginal, low value

# standardize?
ggplot(trend_hyd, aes(x = TreatFreq)) +
  geom_histogram()

ggplot(trend_hyd, aes(x = TreatInt)) +
  geom_histogram()

# frequency model type?
ggplot(trend_hyd, aes(x = TreatFreq, y = Trend)) +
  geom_point()

freqFig <- trend_hyd %>%
  filter(Trend < 4) %>%
  ggplot(aes(x = TreatFreq, y = Trend)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  xlab("Herbicide frequency \n(treatments per year)") +
  ylab("Linear change in log abundance over time") +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10, hjust = 1),
        axis.text = element_text(size = 8))

trend_hyd %>%
  filter(Trend < 4 & TreatFreq < 5) %>%
  ggplot(aes(x = TreatFreq, y = Trend)) +
  geom_point() +
  geom_smooth(method = "lm")

# intensity model type?
ggplot(trend_hyd, aes(x = TreatInt, y = Trend)) +
  geom_point()

intFig <- trend_hyd %>%
  filter(Trend < 4) %>%
  ggplot(aes(x = TreatInt, y = Trend)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  xlab("Herbicide intensity \n(proportion lake treated per year)") +
  ylab("") +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

# combine and export figures
jpeg("output/linear_abundance_herbicide_trends.jpeg", 
    width = 6, height = 3,
    res = 450,
    units = "in")
plot_grid(freqFig, intFig,
          nrow = 1)
dev.off()

# look at a specific case
examp <- trend_hyd %>%
  filter(TreatInt > 0.2 & Trend > 0.18) %>%
  pull(AreaOfInterestID)

jpeg("output/Maude_lake_time_series.jpeg",
     width = 3, height = 3,
     res = 450,
     units = "in")
filter(plant_fwc_hyd, AreaOfInterestID == examp) %>%
  ggplot(aes(x = SurveyYear, y = EstAreaCovered_ha/Area_ha)) +
  geom_line() +
  geom_line(data = filter(ctrl_hyd2, AreaOfInterestID == examp), 
            aes(x = TreatmentYear, y = AreaTreated_ha/Area_ha), color = "red", linetype = "dashed") +
  xlab("Year") +
  ylab("Proportion of lake") +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))
dev.off()
# treatment did seem to manage it: high application in years around outbreak, one more after it crashes
# big treatment before outbreak seems off: large amount of lake treated in 2004 and 2005 when pop is just starting, but no treatment in 2006, when the pop explodes


#### combine FWC plant and ctrl data ####

plant_ctrl_hyd <- ctrl_hyd2 %>%
  mutate(SurveyYear = TreatmentYear) %>% # same year (for difference calcs of abundance, which use lead)
  group_by(AreaOfInterestID, PermanentID, SurveyYear, Area_ha) %>%
  summarise(TotalAreaTreated_ha = sum(AreaTreated_ha)) %>%
  ungroup() %>%
  full_join(ctrl_hyd2 %>%
              mutate(SurveyYear = TreatmentYear + 1) %>% # lag 1 year
              group_by(AreaOfInterestID, PermanentID, SurveyYear) %>%
              summarise(TotalAreaTreatedLag1_ha = sum(AreaTreated_ha)) %>%
              ungroup()) %>%
  full_join(ctrl_hyd2 %>%
              mutate(SurveyYear = TreatmentYear + 2) %>% # lag 2 years
              group_by(AreaOfInterestID, PermanentID, SurveyYear) %>%
              summarise(TotalAreaTreatedLag2_ha = sum(AreaTreated_ha)) %>%
              ungroup()) %>%
  full_join(ctrl_hyd2 %>%
              mutate(SurveyYear = TreatmentYear + 3) %>% # lag 3 years
              group_by(AreaOfInterestID, PermanentID, SurveyYear) %>%
              summarise(TotalAreaTreatedLag3_ha = sum(AreaTreated_ha)) %>%
              ungroup()) %>%
  mutate(TotalPropTreated = TotalAreaTreated_ha / Area_ha,
         TotalPropTreatedLag1 = TotalAreaTreatedLag1_ha / Area_ha,
         TotalPropTreatedLag2 = TotalAreaTreatedLag2_ha / Area_ha,
         TotalPropTreatedLag3 = TotalAreaTreatedLag3_ha / Area_ha) %>%
  left_join(plant_fwc_hyd) %>% # only use data when ctrl occurred (every year included)
  mutate(log_Area_ha = log(Area_ha), # log-transform for offset
         Time = SurveyYear - min(SurveyYear)) %>% # relative time variable
  group_by(AreaOfInterestID, PermanentID) %>%
  arrange(SurveyYear) %>%
  mutate(log_NextPropCovered = lead(log_PropCovered), # includes correction for zero abundance
         log_DiffPropCovered = log_NextPropCovered - log_PropCovered,
         NextPropCovered = lead(PropCovered),
         RatioPropCovered = NextPropCovered/PropCovered) %>%
  ungroup()

# remove missing data for difference analysis
plant_ctrl_hyd2 <- plant_ctrl_hyd %>%
  filter(!is.na(log_NextPropCovered) & !is.na(log_DiffPropCovered))


#### mixed-effects model ####

# visualizations
ggplot(plant_ctrl_hyd2, aes(x = log_DiffPropCovered)) +
  geom_histogram()

ggplot(plant_ctrl_hyd2, aes(x = TotalPropTreated)) +
  geom_histogram()
# zero heavy

ggplot(plant_ctrl_hyd2, aes(x = log_PropCovered)) +
  geom_histogram()
# binary (because of zeros)

ggplot(plant_ctrl_hyd, aes(x = PropCovered)) +
  geom_histogram()
# zero heavy

# treatment and cover
ggplot(plant_ctrl_hyd2, aes(TotalPropTreated, log_DiffPropCovered)) +
  geom_point(alpha = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm")

ggplot(plant_ctrl_hyd2, aes(TotalPropTreatedLag1, log_DiffPropCovered)) +
  geom_point(alpha = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm")

ggplot(plant_ctrl_hyd2, aes(TotalPropTreatedLag2, log_DiffPropCovered)) +
  geom_point(alpha = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm")

ggplot(plant_ctrl_hyd, aes(TotalPropTreatedLag1, PropCovered)) +
  geom_point(alpha = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm")

ggplot(plant_ctrl_hyd, aes(TotalPropTreatedLag2, PropCovered)) +
  geom_point(alpha = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm")

ggplot(plant_ctrl_hyd, aes(TotalPropTreatedLag3, PropCovered)) +
  geom_point(alpha = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm")

# initial cover and change in cover
plant_ctrl_hyd2 %>%
  mutate(ZeroAbund = case_when(EstAreaCovered_ha == 0 ~ "zero",
                               TRUE ~ "> zero")) %>%
  ggplot(aes(x = log_PropCovered, y = log_DiffPropCovered, color = ZeroAbund)) +
  geom_point(alpha = 0.5)
# zeros throw off trend
# constrains change in prop (prop = 1, can only decrease)

plant_ctrl_hyd2 %>%
  filter(EstAreaCovered_ha == 0) %>%
  ggplot(aes(x = log_PropCovered, y = log_NextPropCovered)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1)
# in every case where abundance goes to zero, it stays at zero

plant_ctrl_hyd2 %>%
  filter(EstAreaCovered_ha > 0) %>%
  ggplot(aes(x = log_PropCovered, y = TotalPropTreated)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", formula = y ~ x)
# positively correlated (non-linear)

plant_ctrl_hyd2 %>%
  filter(EstAreaCovered_ha > 0) %>%
  ggplot(aes(x = log_PropCovered, y = log(TotalPropTreated + 1e-8))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x)

# remove zeros
# go back to formatting code to remove zero correction (very small, shouldn't affect outcome)
# add treated variable
plant_ctrl_hyd3 <- plant_ctrl_hyd2 %>%
  filter(SpeciesAcres > 0) %>% # use raw data as cut-off
  mutate(Treated = ifelse(TotalPropTreated == 0, 0, 1),
         log_PropCoveredC = log_PropCovered - mean(log_PropCovered))

# distribution of centered log_PropCovered
ggplot(plant_ctrl_hyd3, aes(x = log_PropCoveredC)) +
  geom_histogram()
  
# random intercept model
mod1 <- lme(log_DiffPropCovered ~ log_PropCoveredC + TotalPropTreated, random =  ~1|AreaOfInterestID, data = plant_ctrl_hyd3)
summary(mod1)
# treatment effect is negative, but not sig

# random slope model
mod2 <- lme(log_DiffPropCovered ~ log_PropCoveredC + TotalPropTreated, random =  ~1 + log_PropCoveredC|AreaOfInterestID, data = plant_ctrl_hyd3)
summary(mod2)
rsquared(mod2) # https://jonlefcheck.net/2013/03/13/r2-for-linear-mixed-effects-models/

# compare models
AIC(mod1, mod2)
# mod2 better

# add autocorrelation to this?

# predictions
pred_lme_hyd <- plant_ctrl_hyd3 %>%
  mutate(log_PropCoveredC = 0) %>%
  mutate(PredictedDiff = predict(mod2, newdata = ., level = 0))

# figure
ggplot(pred_lme_hyd, aes(x = TotalPropTreated, y = log_DiffPropCovered)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = PredictedDiff))


#### piecewise sem ####

# using TotalPropTreated would be nice, but it's zero-inflated
# glmmTMB has zero-inflated models for count data
# could use lake area as an offset, but area treated is continuous
# can't do a zero-inflated gaussian model with log-transformed area treated
# can't do a zero-inflated beta model because no 0/1's allowed - can do censored beta model

# effect of initial population on herbicide
mod3 <- glmer(Treated ~ log_PropCoveredC + (1|AreaOfInterestID), family = binomial, data = plant_ctrl_hyd3)
summary(mod3)

# update mod2 with treated
mod4 <- lme(log_DiffPropCovered ~ log_PropCoveredC + Treated, random =  ~1+log_PropCoveredC|AreaOfInterestID, data = plant_ctrl_hyd3)
summary(mod4)

# fit model
mod5 <- psem(lme(log_DiffPropCovered ~ log_PropCoveredC + Treated, random =  ~1+log_PropCoveredC|AreaOfInterestID, data = plant_ctrl_hyd3),
             glmer(Treated ~ log_PropCoveredC + (1|AreaOfInterestID), family = binomial, data = plant_ctrl_hyd3))

# summary
summary(mod5)
# same estimates as component models


#### whole state model ####

# remove extra years from the lags
plant_ctrl_hyd4 <- plant_ctrl_hyd %>%
  filter(SurveyYear <= max(plant_fwc_hyd$SurveyYear)) # remove extra years from the lags

# function to find longest time interval with the most lakes
time_int_fun <- function(year1){
  
  dat <- plant_ctrl_hyd4
  
  dat2 <- dat %>%
    filter(SurveyYear >= year1 & is.na(EstAreaCovered_ha)) %>% # select missing years
    group_by(AreaOfInterestID, PermanentID) %>%
    summarise(year2 = min(SurveyYear)) %>% # identify first year missing data
    ungroup() %>%
    full_join(dat %>%
                select(AreaOfInterestID, PermanentID) %>%
                unique()) %>% # add all lakes (in case some had no missing data)
    mutate(year2 = replace_na(year2, max(dat$SurveyYear) + 1),   # assign last year (add one to count that year)
           years = year2 - year1) %>%
    group_by(years) %>%
    count() %>% # summarise number of lakes per timespan
    ungroup() %>%
    expand_grid(tibble(years_out = 0:(max(dat$SurveyYear) + 1 - year1))) %>% # expand for every years_out value
    filter(years >= years_out) %>%
    group_by(years_out) %>%
    summarise(lakes = sum(n))
  
  return(dat2)
}

# complete time intervals
plant_ctrl_times <- plant_ctrl_hyd4 %>%
  select(SurveyYear) %>%
  unique() %>%
  mutate(out = map(SurveyYear, time_int_fun)) %>%
  unnest(cols = out) %>%
  mutate(data_points = years_out * lakes) 

ggplot(plant_ctrl_times, aes(x = data_points)) +
  geom_histogram()

# select largest number of datapoints
plant_ctrl_times %>%
  filter(data_points == max(data_points))

# filter 
plant_ctrl_hyd5 <- plant_ctrl_hyd4 %>%
  filter(SurveyYear >= 1998 & SurveyYear <= (1998 + 21)) %>%
  group_by(AreaOfInterestID, PermanentID) %>%
  mutate(NAVals = sum(is.na(EstAreaCovered_ha))) %>%
  ungroup() %>%
  filter(NAVals == 0)

# make sure it worked
length(unique(plant_ctrl_hyd5$AreaOfInterestID))

# visualize
ggplot(plant_ctrl_hyd5, aes(x = SurveyYear, y = PropCovered)) +
  geom_line(aes(color = PermanentID), show.legend = F)

# summarize by year
plant_ctrl_year <- plant_ctrl_hyd5 %>%
  group_by(SurveyYear) %>%
  summarise(across(.cols = c(Area_ha, SpeciesAcres, EstAreaCovered_ha, TotalAreaTreated_ha, TotalAreaTreatedLag1_ha, TotalAreaTreatedLag2_ha, TotalAreaTreatedLag3_ha), sum)) %>%
  ungroup() %>%
  arrange(SurveyYear) %>%
  mutate(PropCovered = EstAreaCovered_ha / Area_ha,
         PropCoveredChange = lead(PropCovered) - PropCovered,
         TotalPropTreated = TotalAreaTreated_ha / Area_ha,
         TotalPropTreatedLag1 = TotalAreaTreatedLag1_ha / Area_ha,
         TotalPropTreatedLag2 = TotalAreaTreatedLag2_ha / Area_ha,
         TotalPropTreatedLag3 = TotalAreaTreatedLag3_ha / Area_ha)

# visualize
ggplot(plant_ctrl_year, aes(x = SurveyYear, y = PropCovered)) +
  geom_line()

ggplot(plant_ctrl_year, aes(x = SurveyYear, y = TotalPropTreated)) +
  geom_line()

ggplot(plant_ctrl_year, aes(x = TotalPropTreated, y = PropCoveredChange)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm")

ggplot(plant_ctrl_year, aes(x = TotalPropTreatedLag1, y = PropCoveredChange)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm")

ggplot(plant_ctrl_year, aes(x = TotalPropTreatedLag1, y = PropCoveredChange)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm")

ggplot(plant_ctrl_year, aes(x = TotalAreaTreatedLag2_ha, y = EstAreaCovered_ha)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm")
# these don't look like Hiatt et al. because the years prior to 1998 influence the trend a lot


#### binary herbicide model ####

# add column to data
plant_ctrl_hyd6 <- plant_ctrl_hyd %>%
  mutate(Treated = case_when(TotalPropTreated > 0 ~ 1,
                             TRUE ~ 0),
         TreatedF = ifelse(Treated == 1, "yes", "no"),
         TreatedLag1 = case_when(TotalPropTreatedLag1 > 0 ~ 1,
                             TRUE ~ 0),
         TreatedLag2 = case_when(TotalPropTreatedLag2 > 0 ~ 1,
                             TRUE ~ 0),
         TreatedLag3 = case_when(TotalPropTreatedLag3 > 0 ~ 1,
                             TRUE ~ 0))

# visualize
ggplot(plant_ctrl_hyd6, aes(x = log_PropCovered, color = TreatedF)) +
  geom_density()
# lakes with no invasion are typically not treated
# treated lakes tend to have larger invasions
# may be good to use data with zero cover removed
# these are more common in the untreated than treated lakes

ggplot(plant_ctrl_hyd6, aes(x = log_DiffPropCovered, color = TreatedF)) +
  geom_density()
# untreated lakes do not change much over time
# this may be due to no invasion
# treated lakes symmetrically increase and decrease

# remove zero cover
# transform for t-test
plant_ctrl_hyd7 <- plant_ctrl_hyd6 %>%
  filter(SpeciesAcres > 0) %>% # use SpeciesAcres so that if it was absent, it's removed
  mutate(logit_PropCovered = logit(PropCovered, adjust = 0.001))

# visualize
ggplot(plant_ctrl_hyd7, aes(x = log_DiffPropCovered, color = TreatedF)) +
  geom_density()
# the change in treated lakes is slightly negative

ggplot(plant_ctrl_hyd7, aes(x = RatioPropCovered, color = TreatedF)) +
  geom_density()

ggplot(plant_ctrl_hyd7, aes(x = logit_PropCovered, color = TreatedF)) +
  geom_density()

ggplot(plant_ctrl_hyd7, aes(x = log10(SpeciesAcres), color = TreatedF)) +
  geom_density() +
  geom_vline(xintercept = log10(0.1))

# sample sizes
plant_ctrl_hyd7 %>%
  group_by(TreatedF) %>%
  count()

# t-test
t.test(formula = logit_PropCovered ~ TreatedF, 
       data = plant_ctrl_hyd7)
# significantly lower prop cover in untreated lakes

# remove low values?
plant_ctrl_hyd7 %>%
  filter(SpeciesAcres > 0.1) %>% # minimum "presence" value
  ggplot(aes(x = log10(SpeciesAcres), color = TreatedF)) +
  geom_density()

plant_ctrl_hyd8 <- plant_ctrl_hyd7 %>%
  filter(SpeciesAcres > 0.1) %>%
  mutate(log_PropCovered = log(PropCovered),  # recalculate to remove zero correction
         log_PropCoveredC = log_PropCovered - mean(log_PropCovered),
         log_DiffPropCovered = log(NextPropCovered) - log_PropCovered) %>%
  filter(!is.na(log_DiffPropCovered))

# visualize
ggplot(plant_ctrl_hyd8, aes(x = log_DiffPropCovered, color = TreatedF)) +
  geom_density()
# the change in treated lakes is slightly negative

ggplot(plant_ctrl_hyd8, aes(x = RatioPropCovered, color = TreatedF)) +
  geom_density() +
  xlim(c(0, 100))

# t-test
t.test(formula = logit_PropCovered ~ TreatedF, 
       data = plant_ctrl_hyd8)
# significantly initial lower prop cover in untreated lakes

# model
mod6 <- lme(log_DiffPropCovered ~ Treated * log_PropCoveredC, random = ~1|AreaOfInterestID, data = plant_ctrl_hyd8)
summary(mod6)

mod6b <- lme(log_DiffPropCovered ~ TreatedLag1 * log_PropCoveredC, random = ~1|AreaOfInterestID, data = plant_ctrl_hyd8)
summary(mod6b)

mod6c <- lme(log_DiffPropCovered ~ TreatedLag2 * log_PropCoveredC, random = ~1|AreaOfInterestID, data = plant_ctrl_hyd8)
summary(mod6c)

mod6d <- lme(log_DiffPropCovered ~ TreatedLag3 * log_PropCoveredC, random = ~1|AreaOfInterestID, data = plant_ctrl_hyd8)
summary(mod6d)

# remove high values?
plant_ctrl_hyd8 %>%
  filter(EstAreaCoveredRaw_ha <= Area_ha) %>%
  ggplot(aes(x = logit_PropCovered, color = TreatedF)) +
  geom_density()
# tried AreaCovered_ha <= Area_ha, but it made less of a difference

# remove unrealistic abundances
# calculate change over time
plant_ctrl_hyd9 <- plant_ctrl_hyd8 %>%
  filter(EstAreaCoveredRaw_ha <= Area_ha) %>%
  mutate(log_PropCovered = log(PropCovered),  # recalculate to remove zero correction
         log_PropCoveredC = log_PropCovered - mean(log_PropCovered),
         log_DiffPropCovered = log(NextPropCovered) - log_PropCovered,
         ChgPropCovered = (NextPropCovered - PropCovered) / PropCovered,
         log_PropCoveredBin = cut_number(log_PropCoveredC, 3,
                                         labels = c("low", "medium", "high"))) %>%
  filter(!is.na(log_DiffPropCovered))

# t-test
t.test(formula = logit_PropCovered ~ TreatedF, 
       data = plant_ctrl_hyd9)
# still significantly different, but  much closer

# variation within year
ggplot(plant_ctrl_hyd9, aes(x = logit_PropCovered, color = TreatedF)) +
  geom_density() +
  facet_wrap(~SurveyYear)

# sample size
plant_ctrl_hyd9 %>%
  group_by(TreatedF) %>%
  count()

# response variable
ggplot(plant_ctrl_hyd9, aes(x = RatioPropCovered)) +
  geom_histogram()

ggplot(plant_ctrl_hyd9, aes(x = log_DiffPropCovered)) +
  geom_histogram()

# treatment effect on change?
ggplot(plant_ctrl_hyd9, aes(x = TreatedF, y = log_DiffPropCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")
# yes

# initial prop?
ggplot(plant_ctrl_hyd9, aes(x = log_PropCoveredC, y = log_DiffPropCovered)) +
  geom_point() +
  geom_smooth(method = "lm")
# less severe than with full dataset

# model
mod7 <- lme(log_DiffPropCovered ~ Treated * log_PropCoveredC, random = ~1|AreaOfInterestID, data = plant_ctrl_hyd9)
summary(mod7)

mod7b <- lme(log_DiffPropCovered ~ TreatedLag1 * log_PropCoveredC, random = ~1|AreaOfInterestID, data = plant_ctrl_hyd9)
summary(mod7b)

mod7c <- lme(log_DiffPropCovered ~ TreatedLag2 * log_PropCoveredC, random = ~1|AreaOfInterestID, data = plant_ctrl_hyd9)
summary(mod7c)

mod7d <- lme(log_DiffPropCovered ~ TreatedLag3 * log_PropCoveredC, random = ~1|AreaOfInterestID, data = plant_ctrl_hyd9)
summary(mod7d)

# change in prop covered
ggplot(plant_ctrl_hyd9, aes(x = ChgPropCovered)) +
  geom_histogram()
# way too skewed
# log-transformation loses 1235 rows of no change

# initial prop by category
ggplot(plant_ctrl_hyd9, aes(x = log_PropCoveredC, y = log_DiffPropCovered)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~log_PropCoveredBin, scales = "free_x")
# still a large effect of initial abundance at high abundances

# sample size
plant_ctrl_hyd9 %>%
  group_by(log_PropCoveredBin, TreatedF) %>%
  count()

# treatment effect by category
ggplot(plant_ctrl_hyd9, aes(x = TreatedF, y = log_DiffPropCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~log_PropCoveredBin, scales = "free_x")

# model
mod8 <- lme(log_DiffPropCovered ~ log_PropCoveredBin * (TreatedLag1 + log_PropCoveredC), random =  ~1|AreaOfInterestID, data = plant_ctrl_hyd9)
summary(mod8)

# predicted values
mod8_pred <- plant_ctrl_hyd9 %>%
  group_by(log_PropCoveredBin) %>%
  summarise(log_PropCoveredC = mean(log_PropCoveredC)) %>% # use overall mean for each group
  ungroup() %>%
  expand_grid(plant_ctrl_hyd9 %>%
                select(Treated, TreatedF) %>%
                unique()) %>%
  mutate(log_DiffPropCovered = predict(mod7, newdata = ., level = 0))

# below from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q1/003336.html
Designmat <- model.matrix(eval(eval(mod7$call$fixed)[-2]), mod8_pred[-4])
predvar <- diag(Designmat %*% mod8$varFix %*% t(Designmat))
mod8_pred$SE <- sqrt(predvar)
mod8_predt$SE2 <- sqrt(predvar+mod8$sigma^2)

# predicted effect by category
ggplot(mod8_pred, aes(x = TreatedF, y = log_DiffPropCovered)) +
  geom_point() +
  facet_wrap(~log_PropCoveredBin, scales = "free_x")


#### binary hydrilla model ####

# lake-year combos for presence/absence
ggplot(plant_fwc_hyd, aes(x = Present)) +
  geom_bar()

ggplot(plant_ctrl_hyd, aes(x = Present)) +
  geom_bar()

# lake counts for proportion of years present
plant_fwc_hyd %>%
  group_by(AreaOfInterestID) %>%
  summarise(Presence = round(sum(Present)/n(), 1)) %>%
  ggplot(aes(x = Presence)) +
  geom_bar()

plant_ctrl_hyd %>%
  filter(!is.na(Present)) %>%
  group_by(AreaOfInterestID) %>%
  summarise(Presence = round(sum(Present)/n(), 1)) %>%
  ggplot(aes(x = Presence)) +
  geom_bar()

# control in previous year and presence
ggplot(plant_ctrl_hyd, aes(TotalPropTreatedLag3, Present)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "glm", method.args = list(family = "binomial"))
# lakes with larger treatment areas tend to have persistent populations
  

#### pre-herbicide data FWC surveys ####

# y matrix: each lake is a row and each year is a column, missing values

# select years
plant_fwc_hyd_pre <- plant_fwc_hyd %>%
  filter(SurveyYear <= (min(ctrl_old$Year))) %>% # goes up to year of first application
  group_by(AreaOfInterestID) %>%
  mutate(NumSurveys = n()) %>%
  ungroup() %>%
  filter(NumSurveys >= MinTime) # remove lakes with too few surveys

# make wide
plant_fwc_hyd_pre2 <- plant_fwc_hyd_pre %>%
  select(AreaOfInterestID, SurveyYear, log_EstAreaCovered_ha) %>%
  arrange(SurveyYear) %>%
  pivot_wider(names_from = SurveyYear, # make years the column names
              values_from = log_EstAreaCovered_ha) %>%
  arrange(AreaOfInterestID) %>%
  mutate(`1982` = case_when(is.na(`1982`) ~ `1983`, # need a first-year value to estimate initial pop size
                            TRUE ~ `1982`))

# extract AOI ID's
plant_fwc_hyd_pre_AOI <- plant_fwc_hyd_pre2 %>%
  select(AreaOfInterestID)

# make matrix
plant_fwc_hyd_pre_mat <- plant_fwc_hyd_pre2 %>%
  select(-AreaOfInterestID) %>%
  as.matrix()

# non-NA values (vector)
plant_fwc_hyd_pre_pos <- plant_fwc_hyd_pre_mat[!is.na(plant_fwc_hyd_pre_mat)]
plant_fwc_hyd_pre_n <- length(plant_fwc_hyd_pre_pos) # number of non-NAs
plant_fwc_hyd_pre_indx <- which(!is.na(plant_fwc_hyd_pre_mat), arr.ind = TRUE)  # index of the non-NAs
plant_fwc_hyd_pre_col <- as.vector(plant_fwc_hyd_pre_indx[, "col"])
plant_fwc_hyd_pre_row <- as.vector(plant_fwc_hyd_pre_indx[, "row"])


#### pre-herbicide model ####

# model file
pre_herb_mod <- 'models/pre_herbicide_model.stan'
cat(paste(readLines(pre_herb_mod)), sep = '\n')

# time model fit
time1 <- Sys.time()

# fit model
pre_herb_fit <- rstan::stan(file = pre_herb_mod, 
                            data = list(y = plant_fwc_hyd_pre_pos, 
                                        TT = ncol(plant_fwc_hyd_pre_mat), 
                                        N = nrow(plant_fwc_hyd_pre_mat), 
                                        n_pos = plant_fwc_hyd_pre_n, 
                                        col_indx_pos = plant_fwc_hyd_pre_col, 
                                        row_indx_pos = plant_fwc_hyd_pre_row), 
                            pars = c("sd_q", "x", "sd_r","u", "x0"), 
                            iter = 1000, chains = 3, thin = 1)

time2 <- Sys.time()
# a little over 20 min
# 516 divergent transitions
# 984 transitions exceeded maximum treedepth
# two chains had low Bayesian Fraction of Missing Information
# large R-hat
# low bulk ESS
# low tail ESS

# update model
time3 <- Sys.time()

pre_herb_fit2 <- rstan::stan(file = pre_herb_mod, 
                             data = list(y = plant_fwc_hyd_pre_pos, 
                                         TT = ncol(plant_fwc_hyd_pre_mat), 
                                         N = nrow(plant_fwc_hyd_pre_mat), 
                                         n_pos = plant_fwc_hyd_pre_n, 
                                         col_indx_pos = plant_fwc_hyd_pre_col, 
                                         row_indx_pos = plant_fwc_hyd_pre_row), 
                             pars = c("sd_q", "x", "sd_r","u", "x0"),
                             iter = 5000, chains = 3, warmup = 2000,
                             control = list(adapt_delta = 0.99, max_treedepth = 15))

time4 <- Sys.time()
# this model can't progress

# fit model with MARSS

# define matrices
marss.list <- list(B = "identity", U = "unequal", Q = "diagonal and equal", 
                   Z = "identity", A = "scaling", R = "diagonal and unequal", 
                   x0 = "unequal", tinitx = 0)

time5 <- Sys.time()
pre_herb_fit3 <- MARSS(plant_fwc_hyd_pre_mat, model = marss.list, silent = 2)
time6 <- Sys.time()

# Error: Stopped in MARSS() before fitting because MARSSkem returned errors.  Try control$trace=1 for more information as the reported error may not be helpful. You can also try method='BFGS' if you are seeing a 'chol' error.
# Error : vector memory exhausted (limit reached?)
# memore: 4.45 GB


 #### post-herbicide data FWC suveys ####

# select data
plant_fwc_hyd_post <-  plant_fwc_hyd %>%
  filter(SurveyYearAdj > min(ctrl_old$Year)) %>% # starts with year after first application
  group_by(AreaOfInterestID) %>%
  mutate(NumSurveys = n()) %>%
  ungroup() %>%
  filter(NumSurveys >= MinTime) # remove lakes with too few surveys

# missing data?
plant_fwc_hyd_post %>%
  group_by(AreaOfInterestID) %>%
  summarise(YearRange = max(SurveyYearAdj) - min(SurveyYearAdj) + 1,
            SurveyYears = n()) %>%
  ungroup() %>%
  filter(YearRange != SurveyYears)
# 109 lakes

# add -99 for missing years
plant_fwc_hyd_post2 <- plant_fwc_hyd_post %>%
  group_by(AreaOfInterestID, PermanentID, AreaOfInterest, Area_ha) %>%
  summarise(MinYear = min(SurveyYearAdj),
            MaxYear = max(SurveyYearAdj)) %>% # year range for a waterbody
  ungroup() %>%
  mutate(SurveyYearAdj = map2(MinYear, MaxYear, seq, by = 1)) %>%
  unnest(cols = SurveyYearAdj) %>% # populate all years for a waterbody
  select(-c(MinYear, MaxYear)) %>%
  full_join(plant_fwc_hyd_post) %>% # missing years will have NA for area metrics
  group_by(AreaOfInterestID) %>%
  arrange(SurveyYearAdj) %>%
  mutate(log_AreaCoveredChange_ha = log_AreaCovered_ha - lag(log_AreaCovered_ha), # change in cover since last survey (NA if missing data)
         log_AreaCoveredMis_ha = replace_na(log_AreaCovered_ha, -99),
         MonthDisplacementMis = replace_na(abs(MonthDisplacementAdj), 0), # make month displacement positive
         SigmaJ = case_when(MonthDisplacementMis != 0 ~ log(1 + MonthDisplacementMis * 0.1), # 10% error by month
                            MonthDisplacementMis == 0 ~ 0.001)) %>% # small error if in survey month
  ungroup() %>%
  arrange(AreaOfInterestID, SurveyYearAdj)

#### post-herbicide surveys ####

# post herbicide surveys with control data
plant_fwc_ctrl_hyd_post <-  plant_fwc_hyd_post2 %>% # plant data has NA for missing values
  rename(YearAdj = SurveyYearAdj) %>%
  left_join(ctrl_hyd_adj %>%
              rename(YearAdj = TreatmentYearAdj) %>%
              mutate(YearAdj = YearAdj + 1)) %>% # add year to control data so that it's in the same row as the survey that followed it
  mutate(TreatmentFrequency = replace_na(TreatmentFrequency, 0), # make values 0 if no control data exists for a waterbody (revisit)
         TreatmentIntensity = replace_na(TreatmentIntensity, 0),
         Treatment = replace_na(Treatment, 0))


#### featured lakes ####

# select six lakes with varying levels of missing data
featured_ID <- plant_fwc_hyd_pre2 %>%
  group_by(AreaOfInterestID) %>%
  summarise(OldYears = n(),
            OldMissingYears = sum(log_AreaCoveredMis_ha == -99),
            OldVarArea = sd(log_AreaCovered_ha, na.rm = T)) %>%
  inner_join(plant_fwc_hyd_post2 %>%
               group_by(AreaOfInterestID) %>%
               summarise(NewYears = n(),
                         NewMissingYears = sum(log_AreaCoveredMis_ha == -99),
                         NewVarArea = sd(log_AreaCovered_ha, na.rm = T))) %>%
  mutate(TotalYears = OldYears + NewYears,
         TotalMissingYears = OldMissingYears + NewMissingYears,
         TotalVarArea = OldVarArea + NewVarArea) %>%
  filter(OldYears >= 10 & NewYears >= 10) %>% # at least 10 years in each dataset
  group_by(TotalMissingYears) %>%
  mutate(MaxVar = max(TotalVarArea)) %>% # choose highest variation
  ungroup() %>%
  filter(TotalMissingYears %in% c(0, 2, 3, 4, 6, 10) & TotalVarArea == MaxVar)

# select lakes from pre dataset
lakes_fwc_hyd_pre <- plant_fwc_hyd_pre2 %>%
  inner_join(featured_ID %>%
              select(AreaOfInterestID, TotalMissingYears)) %>%
  mutate(LakeSize = paste(str_replace(AreaOfInterest, ", Lake", ""), " (", round(Area_ha), " ha)", sep = "") %>%
           fct_reorder(TotalMissingYears))

lakes_fwc_hyd_post <- plant_fwc_ctrl_hyd_post %>%
  inner_join(featured_ID %>%
               select(AreaOfInterestID, TotalMissingYears)) %>%
  mutate(LakeSize = paste(str_replace(AreaOfInterest, ", Lake", ""), " (", round(Area_ha), " ha)", sep = "") %>%
           fct_reorder(TotalMissingYears))


#### edit LW plant data ####

# select hydrilla data
# create date range to merge with FWC
plant_lw_hyd <- plant_lw %>%
  filter(GenusSpecies == "Hydrilla verticillata") %>% # add hydrilla information
  select(County_LW, Lake, PermanentID, Date, SpeciesFrequency) %>%
  mutate(SurveyMin = Date - SurveyDist,
         SurveyMax = Date + SurveyDist) %>%
  rename(SurveyDate_LW = Date,
         PropCover_LW = SpeciesFrequency)


#### combine FWC and LW plant data ####

# merge based on permanentID
# merge if fwc survey is within SurveyDist days of lw survey
plant_lw_fwc_hyd <- plant_fwc_hyd %>%
  rename(SurveyDate_FWC = SurveyDate) %>%
  inner_join(plant_lw_hyd) %>%
  filter(SurveyDate_FWC >= SurveyMin & SurveyDate_FWC <= SurveyMax) %>%
  mutate(PropCover_FWC = AreaCovered_ha / Area_ha,
         SurveyDiff = SurveyDate_LW - SurveyDate_FWC)


#### visualizations ####

# original month displacement
ggplot(plant_fwc_hyd, aes(x = MonthDisplacement)) +
  geom_histogram(binwidth = 1)

# adjusted month displacement
pdf("output/adjusted_month_displacement_hydrilla_fwc.pdf", width = 4, height = 4)
ggplot(plant_fwc_hyd, aes(x = MonthDisplacementAdj)) +
  geom_histogram(binwidth = 1) +
  xlab("Month displacement (J)") +
  ylab("Surveys") +
  def_theme
dev.off()

# adjusted month displacement
pdf("output/month_sampled_hydrilla_fwc.pdf", width = 4, height = 4)
ggplot(plant_fwc_hyd, aes(x = factor(month.abb[MostFreqMonth], levels = month.abb))) +
  geom_bar(stat = "count") +
  xlab("Month surveyed") +
  ylab("Surveys") +
  def_theme
dev.off()

# treatment
pdf("output/herbicide_distribution_hydrilla.pdf", width = 4, height = 4)
ggplot(plant_fwc_ctrl_hyd_post, aes(x = Treatment)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Herbicide (H, proportion lake x applications)") +
  ylab("Surveys") +
  def_theme
dev.off()

# initial treatment visualization
ggplot(plant_fwc_ctrl_hyd_post, aes(x = Treatment, y = log_AreaCoveredChange_ha)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (H, proportion lake x applications)") +
  ylab("Change in hydrilla cover (log(final/initial))") +
  def_theme
# missing data when a year was not surveyed

# featured lakes
pdf("output/featured_lakes_hydrilla_herbicides.pdf", width = 7, height = 5)
lakes_fwc_hyd_pre%>%
  ggplot(aes(SurveyYearAdj, log_AreaCovered_ha)) +
  geom_vline(xintercept = min(ctrl_old$Year), linetype = "dashed") +
  geom_line() +
  geom_point() +
  geom_line(data = lakes_fwc_hyd_post, aes(x = YearAdj)) +
  geom_point(data = lakes_fwc_hyd_post, aes(x = YearAdj, color = Treatment)) +
  geom_text(x = 2018, y = 8, aes(label = paste("missing = ", TotalMissingYears, sep = "")), check_overlap = T, size = 3, hjust = 1) +
  facet_wrap(~ LakeSize) +
  scale_color_viridis_c(name = "Herbicide") +
  xlab("Year") +
  ylab("Hydrilla cover (log hectares)") +
  def_theme
dev.off()

# LW and FWC surveys
pdf("output/lw_fwc_survey_comparison_hydrilla.pdf", width = 4, height = 4)
ggplot(plant_lw_fwc_hyd, aes(PropCover_LW, PropCover_FWC, color = abs(as.numeric(SurveyDiff)))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c(name = "Days\nbetween\nsurveys") +
  xlab("Proportion of lake covered (LW)") +
  ylab("Proportion of lake covered (FWC)") +
  def_theme

ggplot(plant_lw_fwc_hyd, aes(PropCover_LW, PropCover_FWC, color = log(Area_ha))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c(name = "Lake\narea\n(log ha)") +
  xlab("Proportion of lake covered (LW)") +
  ylab("Proportion of lake covered (FWC)") +
  def_theme
dev.off()
  

#### pre-herbicide model ####

# stan code based on example downloaded from: https://laptrinhx.com/smooth-poll-aggregation-using-state-space-modeling-in-stan-from-jim-savage-3654694539/
# observational error is not identifiable when growth rate is a random effect across lakes (strong corelation with lp__ in the pairs plots and error about BFMI low). it may be because of this issue, but I don't totally understand it: https://discourse.mc-stan.org/t/dynamic-panel-data-models-with-stan/5136/16
# when I made growth rate year-specific, process error was more strongly correlated with lp__. BFMI went away, but estimates aren't super close. Maybe they can get closer with stronger priors on the observational error made with other datasets?

# model file
pre_herb_mod <- 'models/hydrilla_pre_herbicide.stan'
cat(paste(readLines(pre_herb_mod)), sep = '\n')

# generate fake data
set.seed(85)
pre_sig_proc = abs(rnorm(1)) # sd of process error
pre_sig_obs = abs(rnorm(1)) # sd of observation error
pre_sig_r = abs(rnorm(1)) # sd of growth rate

pre_test_dat <- tibble(lake = rep(1:10, each = 20), # 10 lakes
                       time = rep(1:20, 10), # 20 years
                       month = sample(plant_fwc_hyd_pre2$MonthDisplacementMis, 200, replace = T)) %>% # sampling time
  rowwise() %>%
  mutate(x = case_when(time == 1 ~ rnorm(1), # random value for log abundance at t = 1
                       TRUE ~ NA_real_),
         sigmaJ = case_when(month != 0 ~ log(1 + month * 0.1), # 10% error by month
                            month == 0 ~ 0.001),
         sig_month = rnorm(1, 0, sigmaJ),
         sig_proc = rnorm(1, 0, pre_sig_proc),
         sig_obs = rnorm(1, 0, pre_sig_obs)) %>%
  ungroup() %>%
  group_by(time) %>%
  mutate(br = rnorm(1, 0.25, pre_sig_r)) %>% # year-specific growth rate
  ungroup()

for(i in 2:20){
  pre_test_dat <- pre_test_dat %>%
    mutate(x = case_when(time == i ~ lag(x) + br + sig_proc, # fill in remaining x
                         TRUE ~ x))
}

pre_test_dat <- pre_test_dat %>%
  mutate(ysurv = x + sig_obs, # add observational error to log-abundance
         y = ysurv + sig_month, # add error for month to observed abundance
         missing = sample(c(1, rep(0, 10)), 200, replace = T), # add indicator variable for missing
         y = case_when(missing == 1 ~ -99,
                       TRUE ~ y)) %>% 
  group_by(lake) %>%
  mutate(ymis = mean(y),
         time_non_mis = case_when(y == -99 ~ NA_integer_,
                                  TRUE ~ time),
         x0 = case_when(time == min(time_non_mis, na.rm = T) ~ y,
                        TRUE ~ NA_real_)) %>%
  ungroup()

pre_test_dat %>%
  filter(y != -99) %>%
  ggplot(aes(x = time, y = x, color = as.factor(lake))) +
  geom_line() +
  geom_point(aes(y = y)) +
  def_theme +
  theme(legend.position = "none")

# create stan data
pre_test_stan_dat <- within(list(), {
  n <- length(unique(pre_test_dat$time))
  l <- length(unique(pre_test_dat$lake))
  y <- matrix(pre_test_dat$y, nrow = length(unique(pre_test_dat$time)), ncol = length(unique(pre_test_dat$lake)))
  sigmaJ <- matrix(pre_test_dat$sigmaJ, nrow = length(unique(pre_test_dat$time)), ncol = length(unique(pre_test_dat$lake)))
  x0 <- filter(pre_test_dat, !is.na(x0))$x0
  ymis <- unique(pre_test_dat$ymis)
})

# fit model
pre_test_fit <- stan(file = pre_herb_mod, data = pre_test_stan_dat, 
                     iter = 5000, chains = 3, warmup = 2000)

# model output
pairs(pre_test_fit, pars = c("sigma_obs", "sigma_proc", "sigma_r", "br_hat", "lp__"))
summary(pre_test_fit, pars = c("sigma_obs", "sigma_proc", "sigma_r", "br_hat", "eta"))$summary

# extract model fit
pre_test_x <- extract(pre_test_fit, pars = "x")[[1]] %>% 
  as_tibble() %>% 
  melt() %>% 
  mutate(variable = as.character(variable)) %>%
  left_join(pre_test_dat %>%
              select(time, lake) %>%
              unique() %>%
              mutate(variable = paste(time, lake, sep = "."))) %>%
  group_by(lake, time) %>% 
  summarise(mean = mean(value),
            median = median(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>%
  ungroup()

# look at results
pdf("output/pre_herbicide_test_model.pdf", width = 5, height = 4)
pre_test_dat %>%
  filter(y != -99) %>%
  full_join(pre_test_x) %>%
  ggplot(aes(time, median, color = as.factor(lake), fill = as.factor(lake))) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, color = NA) +
  geom_line() +
  geom_point(aes(y = y)) +
  xlab("Time") +
  ylab("Area covered (log ha)") +
  def_theme + 
  theme(legend.position = "none")
dev.off()

pdf("output/pre_herbicide_test_model_true_values.pdf", width = 4, height = 4)
pre_test_dat %>%
  filter(y != -99) %>%
  full_join(pre_test_x) %>%
  ggplot(aes(x, median, color = as.factor(lake))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point() +
  xlab("True values") +
  ylab("Estimated true values") +
  def_theme + 
  theme(legend.position = "none")
dev.off()

# create real stan data
pre_herb_stan_dat <- within(list(), {
  n <- length(unique(plant_fwc_hyd_pre2$SurveyYearAdj))
  l <- length(unique(plant_fwc_hyd_pre2$AreaOfInterestID))
  y <- matrix(plant_fwc_hyd_pre2$log_AreaCoveredMis_ha, 
              nrow = length(unique(plant_fwc_hyd_pre2$SurveyYearAdj)), 
              ncol = length(unique(plant_fwc_hyd_pre2$AreaOfInterestID)))
  sigmaJ <- matrix(plant_fwc_hyd_pre2$SigmaJ, 
                   nrow = length(unique(plant_fwc_hyd_pre2$SurveyYearAdj)), 
                   ncol = length(unique(plant_fwc_hyd_pre2$AreaOfInterestID)))
  x0 <- filter(plant_fwc_hyd_pre2, !is.na(x0))$x0
  ymis <- filter(plant_fwc_hyd_pre2, !is.na(ymis))$ymis
})

#### start here: model below is craaazy slow ####

# fit model
pre_herb_fit <- stan(file = pre_herb_mod, data = pre_herb_stan_dat, 
                     iter = 5000, chains = 3, warmup = 2000)

# model output
pairs(pre_herb_fit, pars = c("sigma_obs", "sigma_proc", "sigma_r", "br_hat", "lp__"))
summary(pre_herb_fit, pars = c("sigma_obs", "sigma_proc", "sigma_r", "br_hat", "eta"))$summary


# initiate lists
pre_test_x <- list()
pre_test_y_surv <- list()
pre_test_br <- list()
pre_test_sigma_obs <- list()
pre_test_sigma_proc <- list()
pre_test_x_last <- list()

# loop through each featured lake
for(i in 1:nrow(featured_ID)) {
  
  # ID
  ID <- featured_ID$AreaOfInterestID[i]
  
  # subset data
  pre_test_dat <- lakes_fwc_hyd_pre %>%
    filter(AreaOfInterestID == ID)
  
  # create stan data
  pre_test_stan_data <- within(list(), {
    n <- nrow(pre_test_dat)
    y <- pre_test_dat$log_AreaCoveredMis_ha
    sigmaJ <- pre_test_dat$SigmaJ
    x0 <- filter(pre_test_dat, !is.na(log_AreaCovered_ha))$log_AreaCoveredMis_ha[1] # first non-NA
    ymis <- mean(pre_test_dat$log_AreaCovered_ha, na.rm = T)
  })
  
  # fit model
  pre_test_fit <- stan(file = pre_herb_mod, data = pre_test_stan_data, 
              iter = 5000, chains = 3, warmup = 2000)
  
  # extract process values
  pre_test_x[[i]] <- extract(pre_test_fit, pars = "x")[[1]] %>% 
    as_tibble() %>% 
    melt() %>% 
    mutate(SurveyYearAdj = min(pre_test_dat$SurveyYearAdj) + parse_number(as.character(variable)) - 1) %>% 
    group_by(SurveyYearAdj) %>% 
    summarise(median = median(value),
              lower = quantile(value, 0.025),
              upper = quantile(value, 0.975)) %>%
    mutate(AreaOfInterestID = ID)
  
  pre_test_y_surv[[i]] <- extract(pre_test_fit, pars = "y_surv")[[1]] %>% 
    as_tibble() %>% 
    melt() %>% 
    mutate(SurveyYearAdj = min(pre_test_dat$SurveyYearAdj) + parse_number(as.character(variable)) - 1) %>% 
    group_by(SurveyYearAdj) %>% 
    summarise(median = median(value),
              lower = quantile(value, 0.025),
              upper = quantile(value, 0.975)) %>%
    mutate(AreaOfInterestID = ID)
  
  pre_test_br[[i]] <- summary(pre_test_fit, pars = "br")$summary %>%
    as_tibble() %>%
    mutate(AreaOfInterestID = ID)
  
  pre_test_sigma_obs[[i]] <- summary(pre_test_fit, pars = "sigma_obs")$summary %>%
    as_tibble() %>%
    mutate(AreaOfInterestID = ID)
  
  pre_test_sigma_proc[[i]] <- summary(pre_test_fit, pars = "sigma_proc")$summary %>%
    as_tibble() %>%
    mutate(AreaOfInterestID = ID)
  
  pre_test_x_last[[i]] <- summary(pre_test_fit, pars = "x")$summary %>%
    as_tibble() %>%
    tail(n = 1) %>%
    mutate(AreaOfInterestID = ID)
  
  # save model
  assign(paste("pre_test_fit", ID, sep = "_"), pre_test_fit)
}

# examine models
print(pre_test_fit_290)
print(pre_test_fit_301)
print(pre_test_fit_303)
print(pre_test_fit_311)
print(pre_test_fit_390)
print(pre_test_fit_475)

# combine estimates
pre_test_est <- pre_test_x %>%
  bind_rows() %>%
  rename(median_x = median, lower_x = lower, upper_x = upper) %>%
  full_join(pre_test_y_surv %>%
              bind_rows() %>%
              rename(median_y_surv = median, lower_y_surv = lower, upper_y_surv = upper))

# look at results
pdf("output/featured_lakes_pre_herbicide_fit.pdf", width = 7, height = 5)
lakes_fwc_hyd_pre %>%
  full_join(pre_test_est) %>%
  ggplot(aes(SurveyYearAdj, median_x)) +
  geom_ribbon(aes(ymin = lower_x, ymax = upper_x), alpha = 0.5) +
  geom_line() +
  # geom_point(aes(y = median_y_surv), color = "red") +
  geom_point(aes(y = log_AreaCovered_ha)) +
  facet_wrap(~ LakeSize) +
  xlab("Year") +
  ylab("Hydrilla cover (log hectares)") +
  def_theme
dev.off()

# priors for next model
pre_test_prior <- pre_test_br %>%
  bind_rows() %>%
  select(AreaOfInterestID, mean, se_mean) %>%
  rename(mean_br = mean, se_br = se_mean) %>%
  full_join(pre_test_x_last %>%
              bind_rows() %>%
              select(AreaOfInterestID, mean, se_mean) %>%
              rename(mean_x0 = mean, se_x0 = se_mean)) %>%
  full_join(pre_test_sigma_obs %>%
              bind_rows() %>%
              select(AreaOfInterestID, mean, se_mean) %>%
              rename(mean_obs = mean, se_obs = se_mean)) %>%
  full_join(pre_test_sigma_proc %>%
              bind_rows() %>%
              select(AreaOfInterestID, mean, se_mean) %>%
              rename(mean_proc = mean, se_proc = se_mean)) %>%
  left_join(lakes_fwc_hyd_pre %>%
              select(AreaOfInterest, AreaOfInterestID) %>%
              unique())

write_csv(pre_test_prior, "output/featured_lakes_pre_herbicide_params.csv")


#### post-herbicide model ###

# model file
post_herb_mod <- 'models/hydrilla_post_herbicide.stan'
cat(paste(readLines(post_herb_mod)), sep = '\n')

# initiate lists
post_test_x <- list()
post_test_y_surv <- list()
post_test_br <- list()
post_test_sigma_obs <- list()
post_test_sigma_proc <- list()
post_test_bH <- list()

# loop through each featured lake
for(i in 1:nrow(featured_ID)) {
  
  # ID
  ID <- featured_ID$AreaOfInterestID[i]
  
  # subset data
  post_test_dat <- lakes_fwc_hyd_post %>%
    filter(AreaOfInterestID == ID)
  
  # subset priors
  post_test_prior <- pre_test_prior %>%
    filter(AreaOfInterestID == ID)
  
  # create stan data
  post_test_stan_data <- within(list(), {
    n <- nrow(post_test_dat)
    y <- post_test_dat$log_AreaCoveredMis_ha
    sigmaJ <- post_test_dat$SigmaJ
    H <- post_test_dat$Treatment
    x0 <- post_test_prior$mean_x0
    ymis <- mean(post_test_dat$log_AreaCovered_ha, na.rm = T)
    br_mean <- post_test_prior$mean_br
    #br_sd <- post_test_prior$se_br * 10
    br_sd <- 1
    obs_mean <- post_test_prior$mean_obs
    #obs_sd <- post_test_prior$se_obs * 10
    obs_sd <- 1
    proc_mean <- post_test_prior$mean_proc
    # proc_sd <- post_test_prior$se_proc * 10
    proc_sd <- 1
  })
  
  # fit model
  post_test_fit <- stan(file = post_herb_mod, data = post_test_stan_data, 
                       iter = 5000, chains = 3, warmup = 2000)
  
  # extract process values
  post_test_x[[i]] <- extract(post_test_fit, pars = "x")[[1]] %>% 
    as_tibble() %>% 
    melt() %>% 
    mutate(YearAdj = min(post_test_dat$YearAdj) + parse_number(as.character(variable)) - 1) %>% 
    group_by(YearAdj) %>% 
    summarise(median = median(value),
              lower = quantile(value, 0.025),
              upper = quantile(value, 0.975)) %>%
    mutate(AreaOfInterestID = ID)
  
  post_test_y_surv[[i]] <- extract(post_test_fit, pars = "y_surv")[[1]] %>% 
    as_tibble() %>% 
    melt() %>% 
    mutate(YearAdj = min(post_test_dat$YearAdj) + parse_number(as.character(variable)) - 1) %>% 
    group_by(YearAdj) %>% 
    summarise(median = median(value),
              lower = quantile(value, 0.025),
              upper = quantile(value, 0.975)) %>%
    mutate(AreaOfInterestID = ID)
  
  post_test_br[[i]] <- summary(post_test_fit, pars = "br")$summary %>%
    as_tibble() %>%
    mutate(AreaOfInterestID = ID)
  
  post_test_bH[[i]] <- summary(post_test_fit, pars = "bH")$summary %>%
    as_tibble() %>%
    mutate(AreaOfInterestID = ID)
  
  post_test_sigma_obs[[i]] <- summary(post_test_fit, pars = "sigma_obs")$summary %>%
    as_tibble() %>%
    mutate(AreaOfInterestID = ID)
  
  post_test_sigma_proc[[i]] <- summary(post_test_fit, pars = "sigma_proc")$summary %>%
    as_tibble() %>%
    mutate(AreaOfInterestID = ID)
  
  # save model
  assign(paste("post_test_fit", ID, sep = "_"), post_test_fit)
}

# combine estimates
post_test_est <- post_test_x %>%
  bind_rows() %>%
  rename(median_x = median, lower_x = lower, upper_x = upper) %>%
  full_join(post_test_y_surv %>%
              bind_rows() %>%
              rename(median_y_surv = median, lower_y_surv = lower, upper_y_surv = upper))

# look at results
pdf("output/featured_lakes_post_herbicide_fit.pdf", width = 7, height = 5)
lakes_fwc_hyd_post %>%
  full_join(post_test_est) %>%
  ggplot(aes(YearAdj, median_x)) +
  geom_ribbon(aes(ymin = lower_x, ymax = upper_x), alpha = 0.5) +
  geom_line() +
  # geom_point(aes(y = median_y_surv), color = "red") +
  geom_point(aes(y = log_AreaCovered_ha, color = Treatment)) +
  facet_wrap(~ LakeSize) +
  scale_color_viridis_c(name = "Herbicide") +
  xlab("Year") +
  ylab("Hydrilla cover (log hectares)") +
  def_theme
dev.off()

# parameter estimates
post_test_params <- post_test_br %>%
  bind_rows() %>%
  select(AreaOfInterestID, mean, se_mean) %>%
  rename(mean_br = mean, se_br = se_mean) %>%
  full_join(post_test_bH %>%
              bind_rows() %>%
              select(AreaOfInterestID, mean, se_mean) %>%
              rename(mean_bH = mean, se_bH = se_mean)) %>%
  full_join(post_test_sigma_obs %>%
              bind_rows() %>%
              select(AreaOfInterestID, mean, se_mean) %>%
              rename(mean_obs = mean, se_obs = se_mean)) %>%
  full_join(post_test_sigma_proc %>%
              bind_rows() %>%
              select(AreaOfInterestID, mean, se_mean) %>%
              rename(mean_proc = mean, se_proc = se_mean)) %>%
  left_join(lakes_fwc_hyd_post %>%
              select(AreaOfInterest, AreaOfInterestID, Area_ha) %>%
              unique())

# figure
pdf("output/featured_lakes_herbicide_size.pdf", width = 4, height = 4)
ggplot(post_test_params, aes(x = log(Area_ha), y = mean_bH)) +
  geom_errorbar(aes(ymin = mean_bH - se_bH, ymax = mean_bH + se_bH), width = 0) +
  geom_point(size = 2) +
  xlab("Waterbody area (log hectares)") +
  ylab("Herbicide effect") +
  def_theme
dev.off()



#### old code: edit data ####

# select data associated with Hydrilla
# convert all to hectares




plant_lw_hyd <- plant_lw %>%
  filter(GenusSpecies == "Hydrilla verticillata") %>%
  mutate(Area_ha = ShapeArea * 100)

# see if treated area exceeds total area
ctrl_old_hyd %>%
  filter(AreaTreated_ha > Area_ha) %>%
  select(AreaOfInterest, County_FWC, AreaOfInterestID, PermanentID, AreaTreated_ha, Area_ha) %>%
  unique() %>%
  data.frame()
# checked a few of the large differences and the areas exceed reported areas on websites - round down to 1

# how many lw plant surveys are relevant?
plant_lw_hyd %>%
  select(PermanentID, Date) %>%
  unique() %>%
  inner_join(ctrl_old_hyd %>%
               select(PermanentID) %>%
               unique() %>%
               full_join(ctrl_new_hyd %>%
                           select(PermanentID) %>%
                           unique())) %>% # 184 surveys
  select(PermanentID) %>%
  unique() %>% # 51 lakes
  anti_join(plant_fwc_hyd %>%
               select(PermanentID) %>%
               unique) # 1 lake was not sampled by FWC, the other 50 were

# how consistent are LW and FWC plant surveys?
plant_lw_hyd %>%
  select(Year, PermanentID, SpeciesFrequency) %>%
  inner_join(plant_fwc_hyd %>%
               select(SurveyYear, PermanentID, SpeciesFrequency_ha) %>%
               rename(Year = SurveyYear)) %>%
  ggplot(aes(x = SpeciesFrequency, y = SpeciesFrequency_ha)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  xlab("LW frequency") +
  ylab("FWC frequency")
# generally much lower estimates from FWC than LW

# FWC average plant change
plant_change_fwc_hyd <- plant_fwc_hyd %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha) %>%
  summarise(AreaInitial_ha = AreaCovered_ha[SurveyDate == min(SurveyDate)],
            AreaChange_ha = coef(lm(AreaCovered_ha ~ SurveyDate))[2] * 365,
            FreqChange = coef(lm(SpeciesFrequency_ha ~ SurveyDate))[2] * 365,
            FirstSurvey = min(SurveyDate),
            LastSurvey = max(SurveyDate),
            DaysChange = max(SurveyDate) - min(SurveyDate)) %>%
  ungroup()

# FWC survey-to-survey plant change
survey_change_fwc_hyd <- plant_fwc_hyd %>%
  select(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, SurveyDate, AreaCovered_ha, log_AreaCovered_ha) %>%
  arrange(AreaOfInterestID, SurveyDate) %>%
  group_by(AreaOfInterestID, PermanentID) %>%
  mutate(AreaNotCovered_ha = Area_ha - AreaCovered_ha,
         DaysChange = SurveyDate - lag(SurveyDate),
         OrigAreaCovered_ha = lag(AreaCovered_ha),
         Lag2AreaCovered_ha = lag(AreaCovered_ha, n = 2),
         AreaChange_haPerYear = (AreaCovered_ha - lag(AreaCovered_ha))/as.numeric(DaysChange) * 365,
         RelAreaChange = (AreaCovered_ha - lag(AreaCovered_ha))/lag(AreaCovered_ha),
         SpeciesFrequency_ha = AreaCovered_ha/Area_ha,
         OrigFreq = lag(SpeciesFrequency_ha),
         FreqChange = (SpeciesFrequency_ha - lag(SpeciesFrequency_ha)),
         FreqChangePerYear = (SpeciesFrequency_ha - lag(SpeciesFrequency_ha))/as.numeric(DaysChange) * 365,
         AbsFreqChangePerYear = abs(FreqChangePerYear),
         FreqChangePreFreqPerYear = (SpeciesFrequency_ha - lag(SpeciesFrequency_ha))/(lag(SpeciesFrequency_ha) * as.numeric(DaysChange)) * 365) %>%
  ungroup()

# check ctrl dataset for unique ID-name combos
ctrl_old_hyd %>%
  select(AreaOfInterest, AreaOfInterestID) %>%
  unique() %>%
  mutate(dup = duplicated(AreaOfInterestID)) %>%
  filter(dup == T)

ctrl_new_hyd %>%
  select(AreaOfInterest, AreaOfInterestID) %>%
  unique() %>%
  mutate(dup = duplicated(AreaOfInterestID)) %>%
  filter(dup == T)


# FWC average herbicide and plants
# assign 12/31 to the year herbicide was applied for old dataset (no date info)
# the earliest last survey is in April
# only include herbicide treatments before the last survey (4 month+ lag)
# for known dates, required MinHerbLag days before survey
# herbicide treatments more than MaxHerbLag days before the first survey are not included
plant_herb_fwc_hyd <- plant_change_fwc_hyd %>%
  left_join(ctrl_old_hyd %>%
              filter(AreaTreated_ha > 0) %>%
              group_by(AreaOfInterestID, PermanentID, Year, Area_ha) %>%
              summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
              ungroup() %>%
              mutate(TreatmentDate = paste0(Year, "-12-31") %>% as.Date("%Y-%m-%d")) %>%
              full_join(ctrl_new_hyd %>%
                          filter(AreaTreated_ha > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
                          group_by(AreaOfInterestID, PermanentID, BeginDate, Year, Area_ha) %>%
                          summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
                          ungroup() %>%
                          rename(TreatmentDate = BeginDate)) %>%
              mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha,
                                                TRUE ~ AreaTreated_ha),
                     MinYear = min(Year)) %>%
              select(-c(Area_ha, Year))) %>%
  mutate(AreaTreated_ha = case_when(FirstSurvey - TreatmentDate > MaxHerbLag ~ 0,
                                    LastSurvey - TreatmentDate <= MinHerbLag ~ 0,
                                    is.na(TreatmentDate) ~ 0,
                                    TRUE ~ AreaTreated_ha),
         TreatmentDate = case_when(AreaTreated_ha == 0 ~ NA_Date_,
                                 TRUE ~ TreatmentDate),
         MinYear = case_when(AreaTreated_ha == 0 ~ NA_real_,
                             TRUE ~ MinYear),
         YearsTotal = case_when(AreaTreated_ha > 0 ~ year(LastSurvey) - MinYear,
                                TRUE ~ 1)) %>%
  unique() %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID) %>%
  mutate(rows = n()) %>%
  ungroup() %>%
  filter(!(AreaTreated_ha == 0 & rows > 1)) %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, AreaInitial_ha, AreaChange_ha, FreqChange, DaysChange) %>%
  summarise(TreatmentYears = case_when(sum(AreaTreated_ha) > 0 ~ unique(YearsTotal),
                                       TRUE ~ 0),
            TreatmentFrequency = sum(AreaTreated_ha > 0)/unique(YearsTotal),
            TreatmentIntensity = mean(AreaTreated_ha/AreaInitial_ha),
            TreatmentIntensityLake = mean(AreaTreated_ha/Area_ha)) %>%
  ungroup()

# FWC survey-to-survey herbicide and plants
# herbicide treatments that occurred between the two surveys are included
# unless it is within MinHerbLag days of the last survey
survey_herb_fwc_hyd <- survey_change_fwc_hyd %>%
  left_join(ctrl_old_hyd %>%
              filter(AreaTreated_ha > 0) %>%
              group_by(AreaOfInterestID, PermanentID, Year, Area_ha) %>%
              summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
              ungroup() %>%
              mutate(TreatmentDate = paste0(Year, "-12-31") %>% as.Date("%Y-%m-%d")) %>%
              full_join(ctrl_new_hyd %>%
                          filter(AreaTreated_ha > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
                          group_by(AreaOfInterestID, PermanentID, BeginDate, Area_ha) %>%
                          summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
                          ungroup() %>%
                          rename(TreatmentDate = BeginDate)) %>%
              mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha,
                                                TRUE ~ AreaTreated_ha)) %>%
              select(-c(Area_ha, Year))) %>%
  mutate(SurvTreatDays = SurveyDate - TreatmentDate,
         AreaTreated_ha = case_when(SurvTreatDays > MinHerbLag & SurvTreatDays <= (DaysChange + MinHerbLag) ~ AreaTreated_ha,
                                    TRUE ~ 0),
         SurvTreatDays = case_when(AreaTreated_ha == 0 ~ NA_real_,
                                   TRUE ~ as.numeric(SurvTreatDays)),
         TreatmentDate = case_when(AreaTreated_ha == 0 ~ NA_Date_,
                                 TRUE ~ TreatmentDate)) %>%
  unique() %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, SurveyDate) %>%
  mutate(rows = n()) %>%
  ungroup() %>%
  filter(!(AreaTreated_ha == 0 & rows > 1)) %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, SurveyDate, AreaCovered_ha, OrigAreaCovered_ha, log_AreaCovered_ha, DaysChange, AreaChange_haPerYear, SpeciesFrequency_ha, RelAreaChange, FreqChange, FreqChangePerYear, AbsFreqChangePerYear) %>%
  summarise(TreatmentYears = case_when(sum(AreaTreated_ha) > 0 ~ unique((as.numeric(DaysChange) - MinHerbLag)/365),
                                       TRUE ~ 0),
            TreatmentFrequency = sum(AreaTreated_ha > 0)/unique((as.numeric(DaysChange) - MinHerbLag) * 365),
            TreatmentIntensity = mean(AreaTreated_ha/OrigAreaCovered_ha),
            TreatmentIntensityLake = mean(AreaTreated_ha/Area_ha)) %>%
  ungroup() %>%
  mutate(Herbicide = TreatmentFrequency * TreatmentIntensityLake)

# subset with minimum time points
survey_ts_fwc_hyd <- survey_herb_fwc_hyd %>%
  group_by(AreaOfInterestID, PermanentID) %>%
  mutate(TimePoints = n()) %>%
  ungroup() %>%
  filter(TimePoints >= 5)

# FWC survey-to-survey cumulative herbicide and plants
cltv_herb_fwc_hyd <- survey_change_fwc_hyd %>%
  left_join(ctrl_old_hyd %>%
              filter(AreaTreated_ha > 0) %>%
              group_by(AreaOfInterestID, PermanentID, Year, Area_ha) %>%
              summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
              ungroup() %>%
              mutate(TreatmentDate = paste0(Year, "-12-31") %>% as.Date("%Y-%m-%d")) %>%
              full_join(ctrl_new_hyd %>%
                          filter(AreaTreated_ha > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
                          group_by(AreaOfInterestID, PermanentID, BeginDate, Year, Area_ha) %>%
                          summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
                          ungroup() %>%
                          rename(TreatmentDate = BeginDate)) %>%
              mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha,
                                                TRUE ~ AreaTreated_ha),
                     MinYear = min(Year)) %>%
              select(-c(Area_ha, Year))) %>%
  mutate(SurvTreatDays = SurveyDate - TreatmentDate,
         AreaTreated_ha = case_when(SurvTreatDays > MinHerbLag & SurvTreatDays < MaxHerbLag & !is.na(DaysChange) ~ AreaTreated_ha,
                                    TRUE ~ 0),
         SurvTreatDays = case_when(AreaTreated_ha == 0 ~ NA_real_,
                                   TRUE ~ as.numeric(SurvTreatDays)),
         TreatmentDate = case_when(AreaTreated_ha == 0 ~ NA_Date_,
                                 TRUE ~ TreatmentDate),
         MinYear = case_when(AreaTreated_ha == 0 ~ NA_real_,
                             TRUE ~ MinYear)) %>%
  unique() %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, SurveyDate) %>%
  mutate(rows = n()) %>%
  ungroup() %>%
  filter(!(AreaTreated_ha == 0 & rows > 1)) %>%
  group_by(AreaOfInterest, AreaOfInterestID, PermanentID, County_FWC, Area_ha, SurveyDate, AreaCovered_ha, OrigAreaCovered_ha, DaysChange, AreaChange_haPerYear, SpeciesFrequency_ha, FreqChangePerYear) %>%
  summarise(TreatmentFrequency = sum(AreaTreated_ha > 0)/MaxHerbLag,
            TreatmentIntensity = mean(AreaTreated_ha/OrigAreaCovered_ha),
            TreatmentIntensityLake = mean(AreaTreated_ha/Area_ha)) %>%
  ungroup()

# check LW data for freq > 1
plant_lw_hyd %>%
  filter(SpeciesFrequency > 1)

# LW survey-to-survey plant change
survey_change_lw_hyd <- plant_lw_hyd %>%
  select(Lake, PermanentID, County_LW, Area_ha, Date, SpeciesFrequency) %>%
  arrange(PermanentID, Lake, Date) %>%
  group_by(PermanentID, Lake, County_LW, Area_ha) %>%
  mutate(DaysChange = Date - lag(Date),
         OrigFreq = lag(SpeciesFrequency),
         FreqChange = (SpeciesFrequency - lag(SpeciesFrequency)),
         FreqChangePerYear = (SpeciesFrequency - lag(SpeciesFrequency))/as.numeric(DaysChange) * 365,
         FreqChangePreFreqPerYear = (SpeciesFrequency - lag(SpeciesFrequency))/(lag(SpeciesFrequency) * as.numeric(DaysChange)) * 365) %>%
  ungroup() %>%
  mutate(uniqueID = paste(Lake, County_LW, PermanentID, sep = "_"))

# LW survey-to-survey herbicide and plants
# herbicide treatments that occurred between the two surveys are included
# unless it is within MinHerbLag days of the last survey
survey_herb_lw_hyd <- survey_change_lw_hyd %>%
  left_join(ctrl_old_hyd %>%
              filter(AreaTreated_ha > 0) %>%
              group_by(AreaOfInterestID, PermanentID, Year, Area_ha) %>%
              summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
              ungroup() %>%
              mutate(TreatmentDate = paste0(Year, "-12-31") %>% as.Date("%Y-%m-%d")) %>%
              full_join(ctrl_new_hyd %>%
                          filter(AreaTreated_ha > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>%
                          group_by(AreaOfInterestID, PermanentID, BeginDate, Area_ha) %>%
                          summarise(AreaTreated_ha = sum(AreaTreated_ha)) %>%
                          ungroup() %>%
                          rename(TreatmentDate = BeginDate)) %>%
              mutate(AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha,
                                                TRUE ~ AreaTreated_ha)) %>%
              select(-c(Area_ha, Year))) %>%
  mutate(SurvTreatDays = Date - TreatmentDate,
         AreaTreated_ha = case_when(SurvTreatDays > MinHerbLag & SurvTreatDays <= DaysChange ~ AreaTreated_ha,
                                    TRUE ~ 0),
         SurvTreatDays = case_when(AreaTreated_ha == 0 ~ NA_real_,
                                   TRUE ~ as.numeric(SurvTreatDays)),
         TreatmentDate = case_when(AreaTreated_ha == 0 ~ NA_Date_,
                                   TRUE ~ TreatmentDate)) %>%
  unique() %>%
  group_by(uniqueID, Date) %>%
  mutate(rows = n()) %>%
  ungroup() %>%
  filter(!(AreaTreated_ha == 0 & rows > 1)) %>%
  group_by(uniqueID, PermanentID, Lake, County_LW, Area_ha, Date, DaysChange, SpeciesFrequency, FreqChangePerYear) %>%
  summarise(TreatmentYears = case_when(sum(AreaTreated_ha) > 0 ~ unique((as.numeric(DaysChange) - MinHerbLag)/365),
                                       TRUE ~ 0),
            TreatmentFrequency = sum(AreaTreated_ha > 0)/unique((as.numeric(DaysChange) - MinHerbLag) * 365),
            TreatmentIntensityLake = mean(AreaTreated_ha/Area_ha)) %>%
  ungroup()

  

#### plant_fwc_hyd exploratory figures ####

# small invasions over time
plant_fwc_hyd %>%
  filter(SpeciesAcres < 10) %>%
  ggplot(aes(x = SurveyDate, y = SpeciesFrequency_ha, color = as.factor(AreaOfInterestID))) +
  geom_line(alpha = 0.5) +
  theme(legend.position = "none")
# lots of low values
# non-linear changes over time

# acres by surveyor
plant_fwc_hyd %>%
  ggplot(aes(x = SpeciesAcres)) +
  geom_histogram() +
  facet_wrap(~ Surveyor, scales = "free")

# acres by surveyor when Hydrilla present
plant_fwc_hyd %>%
  filter(SpeciesAcres > 0) %>%
  ggplot(aes(x = SpeciesAcres)) +
  geom_histogram() +
  facet_wrap(~ Surveyor, scales = "free")

# distribution of low values
plant_fwc_hyd %>%
  filter(SpeciesAcres < 20) %>%
  ggplot(aes(x = SpeciesAcres)) +
  geom_histogram(binwidth = 0.1)

# distribution of medium values
plant_fwc_hyd %>%
  filter(SpeciesAcres > 20 & SpeciesAcres < 500) %>%
  ggplot(aes(x = SpeciesAcres)) +
  geom_histogram(binwidth = 10)

# distribution of high values
plant_fwc_hyd %>%
  filter(SpeciesAcres > 500) %>%
  ggplot(aes(x = SpeciesAcres)) +
  geom_histogram(binwidth = 100)

# see if species frequency depends on area source
ggplot(plant_fwc_hyd, aes(x = SpeciesFrequency, y = SpeciesFrequency_ha)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point()
# a lot more are greater than 1
# cluster around zero: probably total shape area includes larger waterbody connected to smaller one


#### plant_change_fwc_hyd exploratory figures ####

ggplot(plant_change_fwc_hyd, aes(DaysChange, AreaChange_ha)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(AreaInitial_ha/Area_ha, AreaChange_ha)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(Area_ha, AreaChange_ha)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(AreaInitial_ha/Area_ha, FreqChange)) +
  geom_point()

ggplot(plant_change_fwc_hyd, aes(DaysChange, FreqChange)) +
  geom_point()


#### survey_change_fwc_hyd exploratory figures ####

# lagged population regulation
pdf("output/hydrilla_survey_lag_by_lake.pdf")
for(i in sort(unique(survey_change_fwc_hyd$AreaOfInterestID))){
  subDat <- filter(survey_change_fwc_hyd, AreaOfInterestID == i)
  lakeName <- unique(subDat$AreaOfInterest)
  print(ggplot(subDat, aes(OrigAreaCovered_ha, AreaCovered_ha)) +
          geom_point() +
          geom_smooth(method = lm) +
          xlab("Area covered by hydrilla at t-1 (log-ha)") +
          ylab("Area covered by hydrilla at t (log-ha)") +
          ggtitle(lakeName) +
          theme_bw())
}
dev.off()

ggplot(survey_change_fwc_hyd, aes(DaysChange, AreaChange_ha)) +
  geom_vline(xintercept = 365, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.5)
# smaller changes with longer times between surveys

survey_change_fwc_hyd %>%
  filter(!is.na(AreaChange_ha)) %>%
  ggplot(aes(DaysChange, log(AreaCovered_ha))) +
  geom_hline(yintercept = log(0.0405), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 365, color = "red", linetype = "dashed") +
  geom_point(alpha = 0.5)
# not all small changes are 0.01 acres (common entry)
# many surveys taken 1 year apart

ggplot(survey_change_fwc_hyd, aes(x = FreqChangePerYear)) +
  geom_histogram()

ggplot(survey_change_fwc_hyd, aes(x = AreaChange_haPerYear)) +
  geom_histogram()

ggplot(survey_change_fwc_hyd, aes(log(Area_ha), FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm")
# no strong relationship between lake size and freq change

ggplot(survey_change_fwc_hyd, aes(log(Area_ha), abs(AreaChange_haPerYear))) +
  geom_point() +
  geom_smooth(method = "lm")
# the absolute value in area change increases with area (expected, constrained variable)

ggplot(survey_change_fwc_hyd, aes(OrigFreq, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm")
# frequency change is constrained by initial frequency

ggplot(survey_change_fwc_hyd, aes(OrigFreq, SpeciesFrequency_ha)) +
  geom_point() +
  geom_smooth(method = "lm")
# trend of temporal autocorrelation, but tons of variation

survey_change_fwc_hyd %>%
  mutate(OrigFreqBin = cut_interval(OrigFreq, n = 10)) %>%
  filter(!is.na(OrigFreqBin))  %>%
  ggplot(aes(x = FreqChangePerYear)) +
  geom_histogram() +
  facet_wrap(~ OrigFreqBin, scales = "free_y")
# the distribution of frequency changes shifts from positive to negative with increasing original frequency
# this is partially because of the constraint noted above

ggplot(survey_change_fwc_hyd, aes(DaysChange, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm")


#### plant_herb_fwc_hyd  exploratory figures ####

ggplot(plant_herb_fwc_hyd, aes(TreatmentFrequency, AreaChange_ha)) +
  geom_point()

ggplot(plant_herb_fwc_hyd, aes(TreatmentIntensity, AreaChange_ha)) +
  geom_point()

ggplot(plant_herb_fwc_hyd, aes(TreatmentFrequency, FreqChange)) +
  geom_point()

ggplot(plant_herb_fwc_hyd, aes(TreatmentIntensity, FreqChange)) +
  geom_point()

plant_herb_fwc_hyd %>%
  filter(!is.na(FreqChange)) %>%
  ggplot(aes(TreatmentFrequency, TreatmentIntensity)) +
  geom_point()

plant_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentFrequency * TreatmentIntensity, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm")


#### survey_herb_fwc_hyd exploratory figures ####

# hydrilla change over time
pdf("output/hydrilla_survey_herb_over_time_by_lake.pdf")
for(i in sort(unique(survey_herb_fwc_hyd$AreaOfInterestID))){
  subDat <- filter(survey_herb_fwc_hyd, AreaOfInterestID == i)
  lakeName <- unique(subDat$AreaOfInterest)
  print(ggplot(subDat, aes(SurveyDate, SpeciesFrequency_ha)) +
          geom_line() +
          geom_point(aes(color = Herbicide), size = 2) +
          xlab("Date") +
          ylab(expression(paste(italic("Hydrilla"), " abundance (prop. lake)", sep = ""))) +
          ggtitle(lakeName) +
          scale_colour_viridis_c(name = "Herbicide\n(prop. lake\ntreated/year)") +
          coord_cartesian(ylim = c(0, 1)) +
          theme_bw())
}
dev.off()

ggplot(survey_herb_fwc_hyd, aes(x = Herbicide)) +
  geom_histogram()

ggplot(survey_herb_fwc_hyd, aes(x = TreatmentFrequency)) +
  geom_histogram()

ggplot(survey_herb_fwc_hyd, aes(x = TreatmentIntensityLake)) +
  geom_histogram()

ggplot(survey_herb_fwc_hyd, aes(Herbicide, AbsFreqChangePerYear, color = as.factor(AreaOfInterestID))) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme(legend.position = "none")

ggplot(survey_herb_fwc_hyd, aes(TreatmentFrequency, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(survey_herb_fwc_hyd, aes(TreatmentIntensity, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(survey_herb_fwc_hyd, aes(TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

survey_herb_fwc_hyd %>%
  filter(TreatmentFrequency > 0) %>%
  ggplot(aes(TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(survey_herb_fwc_hyd, aes(TreatmentIntensity, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

survey_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

pdf("output/hydrilla_change_survey_herbicide.pdf")
survey_herb_fwc_hyd %>%
  #filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()

survey_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()
dev.off()

survey_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentFrequency, TreatmentIntensity)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)


#### survey_ts_fwc_hyd exploratory figures ####

# longest time series?
max_ts <- survey_ts_fwc_hyd %>%
  filter(TimePoints == max(TimePoints))
# max = 38

pdf("output/longest_time_series_lakes.pdf")
survey_ts_fwc_hyd %>%
  filter(AreaOfInterestID %in% max_ts$AreaOfInterestID) %>%
  mutate(Lake_size = paste(str_replace(AreaOfInterest, ", Lake", ""), " (", round(Area_ha), " ha)", sep = "")) %>%
  ggplot(aes(SurveyDate, log_AreaCovered_ha)) +
  geom_vline(xintercept = as.Date("1998-12-31", "%Y-%m-%d"), linetype = "dashed") +
  geom_line() +
  geom_point(aes(color = Herbicide)) +
  facet_wrap(~ Lake_size, scales = "free_x") +
  scale_color_viridis_c(name = "Herbicide\n(prop. lake\ntreated/year)") +
  xlab("Date") +
  ylab("Area covered by hydrilla (log-ha)") +
  def_theme
dev.off()

# minimum time interval
min(survey_ts_fwc_hyd$DaysChange, na.rm = T)
# 162 days
162/30 # ~5.5 months

# time interval distribution
ggplot(survey_ts_fwc_hyd, aes(x = DaysChange)) +
  geom_histogram()

# month samples
survey_ts_fwc_hyd %>%
  mutate(Month = month(SurveyDate)) %>%
  ggplot(aes(x = Month)) +
  geom_histogram(bins = 12)
# clustered between May and October

# samples per year
survey_ts_fwc_hyd %>%
  mutate(Year = year(SurveyDate)) %>%
  group_by(AreaOfInterestID, PermanentID, Year) %>%
  count() %>%
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 1)
# most are once per year

survey_ts_fwc_hyd %>%
  mutate(Year = year(SurveyDate)) %>%
  group_by(AreaOfInterestID, PermanentID, Year) %>%
  count() %>%
  filter(n > 1)
# one was sampled twice in one year - maybe use average


#### cltv_herb_fwc_hyd exploratory figures ####

# hydrilla change over time
pdf("output/hydrilla_cumulative_herb_over_time_by_lake.pdf")
for(i in sort(unique(cltv_herb_fwc_hyd$AreaOfInterestID))){
  subDat <- filter(cltv_herb_fwc_hyd, AreaOfInterestID == i)
  lakeName <- unique(subDat$AreaOfInterest)
  print(ggplot(subDat, aes(SurveyDate, SpeciesFrequency_ha)) +
          geom_line() +
          geom_point(aes(color = TreatmentIntensityLake*TreatmentFrequency), size = 2) +
          xlab("Date") +
          ylab(expression(paste(italic("Hydrilla"), " abundance (prop. lake)", sep = ""))) +
          ggtitle(lakeName) +
          scale_colour_viridis_c(name = "Herbicide\n(prop. lake\ntreated/year)") +
          coord_cartesian(ylim = c(0, 1)) +
          theme_bw())
}
dev.off()

ggplot(cltv_herb_fwc_hyd, aes(TreatmentFrequency, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(cltv_herb_fwc_hyd, aes(TreatmentIntensity, AreaChange_ha)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(cltv_herb_fwc_hyd, aes(TreatmentFrequency, FreqChange)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

cltv_herb_fwc_hyd %>%
  filter(TreatmentFrequency > 0) %>%
  ggplot(aes(TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(cltv_herb_fwc_hyd, aes(TreatmentIntensity, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

cltv_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

pdf("output/hydrilla_change_cumulative_herbicide.pdf")
cltv_herb_fwc_hyd %>%
  #filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()

cltv_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()
dev.off()

cltv_herb_fwc_hyd %>%
  filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentFrequency, TreatmentIntensity)) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  geom_point() +
  geom_smooth(method = "lm", se = F)

cltv_herb_fwc_hyd %>%
  filter(TreatmentFrequency <= 1 & TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensity, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

#### survey_herb_lw_hyd exploratory figures ####

# hydrilla change over time
pdf("output/hydrilla_lw_survey_herb_over_time_by_lake.pdf")
for(i in sort(unique(survey_herb_lw_hyd$uniqueID))){
  subDat <- filter(survey_herb_lw_hyd, uniqueID == i)
  lakeName <- unique(subDat$Lake)
  print(ggplot(subDat, aes(Date, SpeciesFrequency)) +
          geom_line() +
          geom_point(aes(color = TreatmentIntensityLake*TreatmentFrequency), size = 2) +
          xlab("Date") +
          ylab(expression(paste(italic("Hydrilla"), " abundance (prop. lake)", sep = ""))) +
          ggtitle(lakeName) +
          scale_colour_viridis_c(name = "Herbicide\n(prop. lake\ntreated/year)") +
          coord_cartesian(ylim = c(0, 1)) +
          theme_bw())
}
dev.off()

pdf("output/hydrilla_change_lw_survey_herbicide.pdf")
survey_herb_lw_hyd %>%
  #filter(TreatmentIntensity > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()

survey_herb_lw_hyd %>%
  filter(TreatmentIntensityLake > 0) %>%
  ggplot(aes(TreatmentIntensityLake * TreatmentFrequency, FreqChangePerYear)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Herbicide (prop. lake treated/year)") +
  ylab(expression(paste("Change in ", italic("Hydrilla"), " abundance (prop. lake/year)", sep = ""))) +
  theme_bw()
dev.off()


#### test model ####

# longest time series?
max_ts <- plant_fwc_hyd %>%
  group_by(AreaOfInterestID) %>%
  count() %>%
  ungroup() %>%
  filter(n == max(n))

# visualize
plant_fwc_hyd %>%
  filter(AreaOfInterestID %in% max_ts$AreaOfInterestID) %>%
  ggplot(aes(SurveyDate, SpeciesFrequency_ha)) + 
  geom_line() +
  facet_wrap(~AreaOfInterestID)

# test time series
test_dat <- filter(plant_fwc_hyd, AreaOfInterestID == 206) %>%
  mutate(MinDate = min(SurveyDate),
         DateNum = as.numeric(SurveyDate - MinDate))
ggplot(test_dat, aes(DateNum, SpeciesFrequency_ha)) + geom_line()
ggplot(test_dat, aes(DateNum, log_AreaCovered_ha)) + geom_line()
test_ts <- zoo(test_dat$log_AreaCovered_ha, test_dat$DateNum)
plot(test_ts)

# stan data
stan_data <- within(list(), {
  y <- as.vector(test_ts)
  n <- length(test_ts)
})

# herbicide time series
test_dat2 <- filter(survey_ts_fwc_hyd, AreaOfInterestID == 206) %>%
  mutate(MinDate = min(SurveyDate),
         DateNum = as.numeric(SurveyDate - MinDate)) %>%
  filter(!is.na(FreqChangePerYear))
ggplot(test_dat2, aes(DateNum, FreqChangePerYear)) + geom_line()
ggplot(test_dat2, aes(DateNum, Herbicide)) + geom_line()
ggplot(test_dat2, aes(DateNum, log_AreaCovered_ha)) + 
  geom_line() +
  geom_point(aes(color = Herbicide)) +
  scale_color_viridis_c()
ggplot(test_dat2, aes(DateNum, abs(FreqChangePerYear))) + geom_line()
ggplot(test_dat2, aes(Herbicide, abs(FreqChangePerYear))) + geom_point()
ggplot(test_dat2, aes(Herbicide, log_AreaCovered_ha)) + geom_point()
test_ts2 <- zoo(test_dat2$log_AreaCovered_ha, test_dat2$DateNum)
test_x2 <- zoo(test_dat2$Herbicide, test_dat2$DateNum)

# stan data
stan_data2 <- within(list(), {
  y <- as.vector(test_ts2)
  x <- as.vector(test_x2)
  n <- length(test_ts2)
})

# simulated data
test_dat3 <- tibble(days = 0:30,
                    yhat = 10,
                    y = yhat + rnorm(1, 0, 1))

for(i in 2:nrow(test_dat3)){
  test_dat3$yhat[i] <- test_dat3$yhat[i-1] * exp(1.1 - 0.055 * test_dat3$yhat[i-1])
  test_dat3$y[i] = test_dat3$yhat[i] + rnorm(1, 0, 1)
}

ggplot(test_dat3, aes(x = days, y = yhat)) +
  geom_line() +
  geom_point() + 
  geom_line(aes(y = y))

test_ts3 <- zoo(test_dat3$y, test_dat3$days)

# stan data
stan_data3 <- within(list(), {
  y <- as.vector(test_ts3)
  n <- length(test_ts3)
  mu0 <- log(test_dat3$y[1])
})


# model 1: level only, observation error only

# model file
model_file1 <- 'models/test_level_only_obs_error_only.stan'
cat(paste(readLines(model_file1)), sep = '\n')

# fit model
test_fit1 <- stan(file = model_file1, data = stan_data,
            iter = 2000, chains = 4)

mu1 <- get_posterior_mean(test_fit1, par = 'mu')[, 'mean-all chains']

# figures
autoplot(test_ts) + 
  geom_hline(yintercept = mu1, linetype = "dashed") +
  ggtitle("level only, observation error only")

autoplot(test_ts - mu1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggtitle("observation error")

# model 2: level only, observation and process error

# model file
model_file2 <- 'models/test_level_only_obs_error_proc_error.stan'
cat(paste(readLines(model_file2)), sep = '\n')

# fit model
test_fit2 <- stan(file = model_file2, data = stan_data,
                 iter = 2000, chains = 4)

mu2 <- get_posterior_mean(test_fit2, par = 'mu')[, 'mean-all chains']
plot_dat2 <- tibble(model = mu2, 
               Days = test_dat$DateNum,
               observation = test_dat$log_AreaCovered_ha) %>%
  mutate(ObsError = observation - model,
         ProcError = model - lag(model)) %>%
  pivot_longer(-Days, names_to = "yType", values_to = "y")

# figures
plot_dat2 %>%
  filter(!(yType %in% c("ObsError", "ProcError"))) %>%
  ggplot(aes(Days, y, color = yType)) +
  geom_line() +
  ylab("log(Area Covered (ha))") +
  theme_bw() +
  theme(legend.title = element_blank())

plot_dat2 %>%
  filter(yType == "ObsError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Observation error") +
  theme_bw()

plot_dat2 %>%
  filter(yType == "ProcError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Process error") +
  theme_bw()

# model 3: level and slope, observation and process error

# model file
model_file3 <- 'models/test_level_slope_obs_error_proc_error.stan'
cat(paste(readLines(model_file3)), sep = '\n')

# fit model
test_fit3 <- stan(file = model_file3, data = stan_data,
                  iter = 2000, chains = 4)

# extract values
mu3 <- get_posterior_mean(test_fit3, par = 'mu')[, 'mean-all chains']
v3 <- get_posterior_mean(test_fit3, par = 'v')[, 'mean-all chains']

plot_dat3 <- tibble(model = mu3, 
                    slope = c(v3, NA),
                    Days = test_dat$DateNum,
                    observation = test_dat$log_AreaCovered_ha) %>%
  mutate(ObsError = observation - model,
         SlopeProcError = slope - lag(slope),
         LevelProcError = model - lag(model) - lag(slope)) %>%
  pivot_longer(-Days, names_to = "yType", values_to = "y")

# figures
plot_dat3 %>%
  filter(yType %in% c("model", "observation")) %>%
  ggplot(aes(Days, y, color = yType)) +
  geom_line() +
  ylab("log(Area Covered (ha))") +
  theme_bw() +
  theme(legend.title = element_blank())

plot_dat3 %>%
  filter(yType == "ObsError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Observation error") +
  theme_bw()

plot_dat3 %>%
  filter(yType %in% c("SlopeProcError", "LevelProcError")) %>%
  ggplot(aes(Days, y, color = yType)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Process error") +
  theme_bw()

# model 4: deterministic model explanatory variable (regression)

# model file
model_file4 <- 'models/test_level_explanatory_obs_error_only.stan'
cat(paste(readLines(model_file4)), sep = '\n')

# fit model
test_fit4 <- stan(file = model_file4, data = stan_data2,
                  iter = 2000, chains = 4)

# extract values
yhat4 <- get_posterior_mean(test_fit4, par = 'yhat')[, 'mean-all chains']

plot_dat4 <- tibble(model = yhat4, 
                    Days = test_dat2$DateNum,
                    Herbicide = test_dat2$Herbicide,
                    observation = test_dat2$FreqChangePerYear) %>%
  mutate(ObsError = observation - model) %>%
  pivot_longer(-c(Days, Herbicide), names_to = "yType", values_to = "y")

# figures
plot_dat4 %>%
  filter(yType %in% c("model", "observation")) %>%
  ggplot(aes(Days, y, color = yType)) +
  geom_line() +
  ylab("Frequency change per year") +
  theme_bw() +
  theme(legend.title = element_blank())

plot_dat4 %>%
  filter(yType == "ObsError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Observation error") +
  theme_bw()


# model 5: stochastic model with explanatory variable

# model file
model_file5 <- 'models/test_level_explanatory_obs_error_proc_error.stan'
cat(paste(readLines(model_file5)), sep = '\n')

# fit model
test_fit5 <- stan(file = model_file5, data = stan_data2,
                  iter = 2000, chains = 4)

# extract values
yhat5 <- get_posterior_mean(test_fit5, par = 'yhat')[, 'mean-all chains']
mu5 <- get_posterior_mean(test_fit5, par = 'mu')[, 'mean-all chains']

plot_dat5 <- tibble(mu = mu5,
                    model = yhat5, 
                    Days = test_dat2$DateNum,
                    Herbicide = test_dat2$Herbicide,
                    observation = test_dat2$FreqChangePerYear) %>%
  mutate(ObsError = observation - model,
         ProcError = mu - lag(mu)) %>%
  pivot_longer(-c(Days, Herbicide), names_to = "yType", values_to = "y")

# figures
plot_dat5 %>%
  filter(yType %in% c("model", "observation")) %>%
  ggplot(aes(Days, y, color = yType)) +
  geom_line() +
  ylab("Frequency change per year") +
  theme_bw() +
  theme(legend.title = element_blank())

plot_dat5 %>%
  filter(yType == "ObsError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Observation error") +
  theme_bw()

plot_dat5 %>%
  filter(yType == "ProcError") %>%
  ggplot(aes(Days, y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Process error") +
  theme_bw()

# model 6: Ricker model, observation error only

# model file
model_file6 <- 'models/test_ricker_obs_error_only.stan'
cat(paste(readLines(model_file6)), sep = '\n')

# fit model
test_fit6 <- stan(file = model_file6, data = stan_data3,
                  iter = 2000, chains = 4)

print(test_fit6, max = 50)
# generally bad at estimating parameters, could give stricter priors...

# try simple model
test_fit6b <- stan(file = model_file1, data = stan_data3,
                  iter = 2000, chains = 4)

mu6b <- get_posterior_mean(test_fit6b, par = 'mu')[, 'mean-all chains']

# figures
autoplot(test_ts3) + 
  geom_hline(yintercept = mu6b, linetype = "dashed") +
  ggtitle("level only, observation error only")

autoplot(test_ts3 - mu6b) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggtitle("observation error")
