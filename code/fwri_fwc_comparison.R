#### info ####

# goal: see how FWRI compares to FWC


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)

# import data
fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv",
                 col_types = list(PermanentID = col_character()))
fwc <- read_csv("intermediate-data/FWC_plant_formatted.csv")


#### subset data functions ####

fwri_fun <- function(code){
  
  dat_out <- fwri %>%
    select(AOI, Lake, PermanentID, Year, Date, Site) %>% # start with all surveys (no species means abundance = 0)
    unique() %>%  # don't need a row for each species
    left_join(fwri %>%  # add species information
                filter(Code == code) %>%
                select(AOI, Lake, PermanentID, Year, Date, Site, Abundance)) %>%
    mutate(Abundance = replace_na(Abundance, 0), # species cover 0 when it wasn't in a survey
           FWRI_Date = as.Date(Date),
           FWRI_Month = month(Date)) %>%
    group_by(AOI, Lake, PermanentID, Year, FWRI_Month, FWRI_Date) %>% # checked above for duplicates within sites - none
    summarise(Area_sites = n(),
              AreaCovered_sites1 = sum(Abundance > 0),
              AreaCovered_sites2 = sum(Abundance > 1),
              AreaCovered_sites3 = sum(Abundance > 2)) %>%
    ungroup() %>%
    mutate(AOI = case_when(AOI %in% c("NorthConway", "SouthConway") ~ "Conway", # combine these lakes (one lake in FWC data)
                           TRUE ~ AOI),
           Lake = case_when(Lake %in% c("North Conway", "South Conway") ~ "Conway",
                            TRUE ~ Lake)) %>%
    group_by(AOI, Lake, PermanentID, Year, FWRI_Month) %>%
    summarise(Area_sites = sum(Area_sites), # combine north and south conway data
              AreaCovered_sites1 = sum(AreaCovered_sites1),
              AreaCovered_sites2 = sum(AreaCovered_sites2),
              AreaCovered_sites3 = sum(AreaCovered_sites3),
              FWRI_Date = min(FWRI_Date)) %>%
    ungroup() %>%
    mutate(FWRI_PropCovered1 = AreaCovered_sites1 / Area_sites,
           FWRI_PropCovered2 = AreaCovered_sites2 / Area_sites,
           FWRI_PropCovered3 = AreaCovered_sites3 / Area_sites,
           FWRI_SitesUnOcc1 = Area_sites - AreaCovered_sites1,
           FWRI_SitesUnOcc2 = Area_sites - AreaCovered_sites2,
           FWRI_SitesUnOcc3 = Area_sites - AreaCovered_sites3) %>%
    rename(FWRI_Year = Year, FWRI_SitesOcc1 = AreaCovered_sites1, FWRI_SitesOcc2 = AreaCovered_sites2, FWRI_SitesOcc3 = AreaCovered_sites3) %>%
    select(-Area_sites) %>%
    filter(!(AOI %in% c("Orange", "Eustis2")))
  
  return(dat_out)
}

fwc_fun <- function(species){
  
  dat_out <- fwc %>% # start with all surveys (no species means abundance = 0)
    filter(SpeciesName != species) %>%
    select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate) %>%
    unique() %>% # don't need a row for each species
    left_join(fwc %>% # add species information
                filter(SpeciesName == species) %>%
                select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, SpeciesAcres)) %>%
    mutate(SpeciesAcres = replace_na(SpeciesAcres, 0), # species cover 0 when it wasn't in a survey
           Area_ha = ShapeArea * 100, # convert lake area from km-squared to hectares
           AreaCovered_ha = SpeciesAcres * 0.405, # convert plant cover from acres to hectares
           AreaCovered_ha = case_when(AreaCovered_ha > Area_ha ~ Area_ha,
                                      TRUE ~ AreaCovered_ha), # make plant cover the size of the lake area if it exceeds it
           FWC_PropCovered = AreaCovered_ha / Area_ha,
           FWC_AreaCovered_ha = round(AreaCovered_ha),
           FWC_AreaUnCovered_ha = round(Area_ha) - FWC_AreaCovered_ha,
           FWC_Year = year(SurveyDate),
           FWC_Month = month(SurveyDate)) %>%
    group_by(AreaOfInterestID, FWC_Year) %>% # remove reports of zero when there is another report that year
    mutate(AreaCoveredAnnAvg_ha = mean(AreaCovered_ha),
           FirstSurveyPerYear = min(SurveyDate)) %>%
    ungroup() %>%
    filter(!(AreaCoveredAnnAvg_ha > 0 & AreaCovered_ha == 0) & # kept report is non-zero
             !(AreaCovered_ha == 0 & SurveyDate != FirstSurveyPerYear)) %>% # kept report is zero
    select(AreaOfInterest, AreaOfInterestID, PermanentID, FWC_Year, FWC_Month, SurveyDate, FWC_AreaCovered_ha, FWC_AreaUnCovered_ha, FWC_PropCovered) %>%
    rename(FWC_Date = SurveyDate) %>%
    filter(!(AreaOfInterestID %in% c(402, 469, 42, 220, 436)))
  # Harris (177) and Little Harris (244) are the same lake
  # Red Water (365) and Little Red Water (250) are the same lake
  
  return(dat_out)
}


#### process data functions ####

# change over time with glm
chg_fun <- function(dat){
  
  dat_out <- dat %>%
    group_by(AreaOfInterestID) %>%
    mutate(FWC_MinDate = min(FWC_Date),
           FWC_Days = as.numeric(FWC_Date - FWC_MinDate),
           FWRI_MinDate = min(FWRI_Date),
           FWRI_Days = as.numeric(FWRI_Date - FWRI_MinDate),
           Dates = n()) %>%
    ungroup() %>%
    filter(Dates > 1) %>%
    nest(-AreaOfInterestID) %>%
    mutate(FWC_Change = map(data, ~ coef(glm(cbind(FWC_AreaCovered_ha, FWC_AreaUnCovered_ha) ~ FWC_Days, data = ., family = "binomial"))["FWC_Days"]),
           FWRI_Change1 = map(data, ~ coef(glm(cbind(FWRI_SitesOcc1, FWRI_SitesUnOcc1) ~ FWRI_Days, data = ., family = "binomial"))["FWRI_Days"]),
           FWRI_Change2 = map(data, ~ coef(glm(cbind(FWRI_SitesOcc2, FWRI_SitesUnOcc2) ~ FWRI_Days, data = ., family = "binomial"))["FWRI_Days"]),
           FWRI_Change3 = map(data, ~ coef(glm(cbind(FWRI_SitesOcc3, FWRI_SitesUnOcc3) ~ FWRI_Days, data = ., family = "binomial"))["FWRI_Days"])) %>%
    unnest(FWC_Change, FWRI_Change1, FWRI_Change2, FWRI_Change3) %>%
    select(-data)
  
  return(dat_out)
}


#### compare hydrilla cover ####

# subset data
fwri_hydr <- fwri_fun("HYDR")
fwc_hydr <- fwc_fun("Hydrilla verticillata")

# check lakes that will be excluded
anti_join(fwri_hydr, fwc_hydr) %>%
  select(AOI, PermanentID) %>%
  unique()

# check date formatting
tz(fwri_hydr$FWRI_Date)
tz(fwc_hydr$FWC_Date)

# combine data
hydr <- fwri_hydr %>%
  inner_join(fwc_hydr) %>%
  mutate(TimeDiff = abs(FWC_Date - FWRI_Date)) %>%
  filter(TimeDiff <= 150)

# duplicate AOI's
hydr %>%
  group_by(PermanentID) %>%
  summarise(FWRI_AOI = length(unique(AOI)),
            FWC_AOI = length(unique(AreaOfInterestID))) %>%
  ungroup() %>%
  filter(FWRI_AOI > 1 | FWC_AOI > 1) %>%
  left_join(hydr %>%
              select(PermanentID, AOI, AreaOfInterest, AreaOfInterestID) %>%
              unique())
# two lakes that have same coordinates in FWC dataset

# duplicate FWC dates for each FWRI
hydr %>%
  group_by(PermanentID, AreaOfInterestID, FWRI_Date) %>%
  summarise(FWC_surveys = length(unique(FWC_Date))) %>%
  ungroup() %>%
  filter(FWC_surveys > 1)

# figure
ggplot(hydr, aes(FWC_PropCovered, FWRI_PropCovered1, 
                 group = as.factor(AreaOfInterestID), 
                 color = as.numeric(TimeDiff))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point()

# check duplicate lakes
filter(hydr, PermanentID == "112047993") %>% data.frame()  # all 0's
filter(hydr, PermanentID == "120024301") %>% data.frame() # different values

# remove duplicates
hydr2 <- hydr %>%
  filter(AreaOfInterestID != 365)

# remove zeros
hydrP <- hydr2 %>%
  filter(!(FWC_PropCovered == 0 & FWRI_PropCovered1 == 0))

# change over time
hydrC <- chg_fun(hydr2)


#### compare water hyacinth cover ####

# subset data
fwri_wahy <- fwri_fun("WAHY")
fwc_wahy <- fwc_fun("Eichhornia crassipes")

# combine data
wahy <- fwri_wahy %>%
  inner_join(fwc_wahy) %>%
  mutate(TimeDiff = abs(FWC_Date - FWRI_Date)) %>%
  filter(TimeDiff <= 150)

# duplicate AOI's
wahy %>%
  group_by(PermanentID) %>%
  summarise(FWRI_AOI = length(unique(AOI)),
            FWC_AOI = length(unique(AreaOfInterestID))) %>%
  ungroup() %>%
  filter(FWRI_AOI > 1 | FWC_AOI > 1) %>%
  left_join(wahy %>%
              select(PermanentID, AOI, AreaOfInterest, AreaOfInterestID) %>%
              unique())
# two lakes that have same coordinates in FWC dataset

# duplicate FWC dates for each FWRI
wahy %>%
  group_by(PermanentID, AreaOfInterestID, FWRI_Date) %>%
  summarise(FWC_surveys = length(unique(FWC_Date))) %>%
  ungroup() %>%
  filter(FWC_surveys > 1)

# figure
ggplot(wahy, aes(FWC_PropCovered, FWRI_PropCovered1, 
                 group = as.factor(AreaOfInterestID), 
                 color = as.numeric(TimeDiff))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point()

# check duplicate lakes
filter(wahy, PermanentID == "112047993") %>% data.frame()  # all 0's
filter(wahy, PermanentID == "120024301") %>% data.frame() # all 0's

# remove duplicates
wahy2 <- wahy %>%
  filter(!(AreaOfInterestID %in% c(365, 244)))

# remove zeros
wahyP <- wahy2 %>%
  filter(!(FWC_PropCovered == 0 & FWRI_PropCovered1 == 0))

# change over time
wahyC <- chg_fun(wahy2)


#### compare water lettuce cover ####

# subset data
fwri_wale <- fwri_fun("WALE")
fwc_wale <- fwc_fun("Pistia stratiotes")

# combine data
wale <- fwri_wale %>%
  inner_join(fwc_wale) %>%
  mutate(TimeDiff = abs(FWC_Date - FWRI_Date)) %>%
  filter(TimeDiff <= 150)

# duplicate AOI's
wale %>%
  group_by(PermanentID) %>%
  summarise(FWRI_AOI = length(unique(AOI)),
            FWC_AOI = length(unique(AreaOfInterestID))) %>%
  ungroup() %>%
  filter(FWRI_AOI > 1 | FWC_AOI > 1) %>%
  left_join(wale %>%
              select(PermanentID, AOI, AreaOfInterest, AreaOfInterestID) %>%
              unique())
# two lakes that have same coordinates in FWC dataset

# duplicate FWC dates for each FWRI
wale %>%
  group_by(PermanentID, AreaOfInterestID, FWRI_Date) %>%
  summarise(FWC_surveys = length(unique(FWC_Date))) %>%
  ungroup() %>%
  filter(FWC_surveys > 1)

# figure
ggplot(wale, aes(FWC_PropCovered, FWRI_PropCovered1, 
                 group = as.factor(AreaOfInterestID), 
                 color = as.numeric(TimeDiff))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point()

# check duplicate lakes
filter(wale, PermanentID == "112047993") %>% data.frame()  # all 0's
filter(wale, PermanentID == "120024301") %>% data.frame() # all 0's

# remove duplicates
wale2 <- wale %>%
  filter(!(AreaOfInterestID %in% c(365, 244))) 

# remove zeros
waleP <- wale2 %>%
  filter(!(FWC_PropCovered == 0 & FWRI_PropCovered1 == 0))

# change over time
waleC <- chg_fun(wale2)


#### proportion covered figures ####

# make long
hydrP2 <- hydrP %>%
  pivot_longer(cols = c(FWRI_PropCovered1, FWRI_PropCovered2, FWRI_PropCovered3),
               names_to = "FWRI_MinAbundance",
               values_to = "FWRI_PropCovered",
               names_prefix = "FWRI_PropCovered") %>%
  mutate(DaysDiff = as.numeric(TimeDiff))

waleP2 <- waleP %>%
  pivot_longer(cols = c(FWRI_PropCovered1, FWRI_PropCovered2, FWRI_PropCovered3),
               names_to = "FWRI_MinAbundance",
               values_to = "FWRI_PropCovered",
               names_prefix = "FWRI_PropCovered") %>%
  mutate(DaysDiff = as.numeric(TimeDiff))

wahyP2 <- wahyP %>%
  pivot_longer(cols = c(FWRI_PropCovered1, FWRI_PropCovered2, FWRI_PropCovered3),
               names_to = "FWRI_MinAbundance",
               values_to = "FWRI_PropCovered",
               names_prefix = "FWRI_PropCovered") %>%
  mutate(DaysDiff = as.numeric(TimeDiff))

# deviations from 1:1 line
hydrP_diff <- tibble(FWRI_MinAbundance = c(1, 2, 3)) %>%
  mutate(MeanDiff = NA,
         AreaOfInterestID = NA_real_,
         DaysDiff = NA_real_)
waleP_diff <- tibble(FWRI_MinAbundance = c(1, 2, 3)) %>%
  mutate(MeanDiff = NA,
         AreaOfInterestID = NA_real_,
         DaysDiff = NA_real_)
wahyP_diff <- tibble(FWRI_MinAbundance = c(1, 2, 3)) %>%
  mutate(MeanDiff = NA,
         AreaOfInterestID = NA_real_,
         DaysDiff = NA_real_)

# calculate estimates
for(i in 1:3){
  
  # subset data
  hydr_temp <- filter(hydrP2, FWRI_MinAbundance == i) %>%
    mutate(Diff = abs(FWRI_PropCovered - FWC_PropCovered))
  wale_temp <- filter(waleP2, FWRI_MinAbundance == i) %>%
    mutate(Diff = abs(FWRI_PropCovered - FWC_PropCovered))
  wahy_temp <- filter(wahyP2, FWRI_MinAbundance == i) %>%
    mutate(Diff = abs(FWRI_PropCovered - FWC_PropCovered))
  
  # enter value
  hydrP_diff$MeanDiff[i] <- round(mean(hydr_temp$Diff), 3)
  waleP_diff$MeanDiff[i] <- round(mean(wale_temp$Diff), 3)
  wahyP_diff$MeanDiff[i] <- round(mean(wahy_temp$Diff), 3)
  
}

# figure
pdf("output/fwri_fwc_comparison_invasive_cover.pdf", width = 7.5, height = 4)
ggplot(hydrP2, aes(FWC_PropCovered, FWRI_PropCovered, 
                 group = as.factor(AreaOfInterestID), 
                 color = DaysDiff)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_text(data = hydrP_diff, x = (max(hydrP2$FWC_PropCovered) - min(hydrP2$FWC_PropCovered)) / 2, y = max(hydrP2$FWRI_PropCovered), hjust = 0.5, color = "black", aes(label = paste("diff = ", MeanDiff, sep = ""))) +
  facet_wrap(~FWRI_MinAbundance) +
  xlab("Proportion of lake (FWC)") +
  ylab("Proportion of lake (FWRI)") +
  ggtitle("Hydrilla")

ggplot(waleP2, aes(FWC_PropCovered, FWRI_PropCovered, 
                     group = as.factor(AreaOfInterestID), 
                     color = DaysDiff)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_text(data = waleP_diff, x = (max(waleP2$FWC_PropCovered) - min(waleP2$FWC_PropCovered)) / 2, y = max(waleP2$FWRI_PropCovered), hjust = 0.5, color = "black", aes(label = paste("diff = ", MeanDiff, sep = ""))) +
  facet_wrap(~FWRI_MinAbundance) +
  xlab("Proportion of lake (FWC)") +
  ylab("Proportion of lake (FWRI)") +
  ggtitle("Water lettuce")

ggplot(wahyP2, aes(FWC_PropCovered, FWRI_PropCovered, 
                     group = as.factor(AreaOfInterestID), 
                     color = DaysDiff)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_text(data = wahyP_diff, x = (max(wahyP2$FWC_PropCovered) - min(wahyP2$FWC_PropCovered))/2, y = max(wahyP2$FWRI_PropCovered), hjust = 0.5, color = "black", aes(label = paste("diff = ", MeanDiff, sep = ""))) +
  facet_wrap(~FWRI_MinAbundance) +
  xlab("Proportion of lake (FWC)") +
  ylab("Proportion of lake (FWRI)") +
  ggtitle("Water hyacinth")
dev.off()


#### specific points ####

waleP2 %>%
  filter(FWC_PropCovered > 0.2) %>%
  data.frame()

wahyP2 %>%
  filter(FWC_PropCovered > 0.2) %>%
  data.frame()


#### change over time figures ####

# make long
hydrC2 <- hydrC %>%
  pivot_longer(cols = c(FWRI_Change1, FWRI_Change2, FWRI_Change3),
               names_to = "FWRI_MinAbundance",
               values_to = "FWRI_Change",
               names_prefix = "FWRI_Change")

waleC2 <- waleC %>%
  pivot_longer(cols = c(FWRI_Change1, FWRI_Change2, FWRI_Change3),
               names_to = "FWRI_MinAbundance",
               values_to = "FWRI_Change",
               names_prefix = "FWRI_Change")

wahyC2 <- wahyC %>%
  pivot_longer(cols = c(FWRI_Change1, FWRI_Change2, FWRI_Change3),
               names_to = "FWRI_MinAbundance",
               values_to = "FWRI_Change",
               names_prefix = "FWRI_Change")

# deviations from 1:1 line
hydrC_diff <- tibble(FWRI_MinAbundance = c(1, 2, 3)) %>%
  mutate(MeanDiff = NA,
         AreaOfInterestID = NA_real_,
         DaysDiff = NA_real_)
waleC_diff <- tibble(FWRI_MinAbundance = c(1, 2, 3)) %>%
  mutate(MeanDiff = NA,
         AreaOfInterestID = NA_real_,
         DaysDiff = NA_real_)
wahyC_diff <- tibble(FWRI_MinAbundance = c(1, 2, 3)) %>%
  mutate(MeanDiff = NA,
         AreaOfInterestID = NA_real_,
         DaysDiff = NA_real_)

# calculate estimates
for(i in 1:3){
  
  # subset data
  hydr_temp <- filter(hydrC2, FWRI_MinAbundance == i) %>%
    mutate(Diff = abs(FWRI_Change - FWC_Change))
  wale_temp <- filter(waleC2, FWRI_MinAbundance == i) %>%
    mutate(Diff = abs(FWRI_Change - FWC_Change))
  wahy_temp <- filter(wahyC2, FWRI_MinAbundance == i) %>%
    mutate(Diff = abs(FWRI_Change - FWC_Change))
  
  # enter value
  hydrC_diff$MeanDiff[i] <- round(mean(hydr_temp$Diff), 3)
  waleC_diff$MeanDiff[i] <- round(mean(wale_temp$Diff), 3)
  wahyC_diff$MeanDiff[i] <- round(mean(wahy_temp$Diff), 3)
  
}

# figure
pdf("output/fwri_fwc_comparison_invasive_cover_change.pdf", width = 7.5, height = 4)
ggplot(hydrC2, aes(FWC_Change, FWRI_Change, 
                   group = as.factor(AreaOfInterestID))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_text(data = hydrC_diff, x = 0, y = max(hydrC2$FWRI_Change), hjust = 0.5, color = "black", aes(label = paste("diff = ", MeanDiff, sep = ""))) +
  facet_wrap(~FWRI_MinAbundance) +
  xlab("Change in cover over time (FWC)") +
  ylab("Change in cover over time (FWRI)") +
  ggtitle("Hydrilla")

ggplot(waleC2, aes(FWC_Change, FWRI_Change, 
                   group = as.factor(AreaOfInterestID))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_text(data = waleC_diff, x = -0.015, y = max(waleC2$FWRI_Change), hjust = 0.5, color = "black", aes(label = paste("diff = ", MeanDiff, sep = ""))) +
  facet_wrap(~FWRI_MinAbundance) +
  xlab("Change in cover over time (FWC)") +
  ylab("Change in cover over time (FWRI)") +
  ggtitle("Water lettuce")

ggplot(wahyC2, aes(FWC_Change, FWRI_Change, 
                   group = as.factor(AreaOfInterestID))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_text(data = wahyC_diff, x = 0, y = max(wahyC2$FWRI_Change), hjust = 0.5, color = "black", aes(label = paste("diff = ", MeanDiff, sep = ""))) +
  facet_wrap(~FWRI_MinAbundance) +
  xlab("Change in cover over time (FWC)") +
  ylab("Change in cover over time (FWRI)") +
  ggtitle("Water hyacinth")
dev.off()