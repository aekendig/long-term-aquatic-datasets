#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)
library(fuzzyjoin)

# figure settings
source("code/settings/figure_settings.R")

# import data
fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv")
fwc <- read_csv("intermediate-data/FWC_plant_formatted_temporal_coverage.csv")
fwc_syn <- read_csv("intermediate-data/FWC_plant_survey_species_synonyms.csv")


#### find overlap between surveys ####

# rename columns
# format columns to match better
fwri2 <- fwri %>%
  rename_with(.cols = -PermanentID, .fn = ~ paste0("fwri_", .x)) %>%
  mutate(fwri_Date = as.Date(fwri_Date),
         fwri_ScientificName = str_replace(fwri_ScientificName, "sp.", "spp."))
fwc2 <- fwc %>%
  rename(Year = SurveyYear,
         Date = SurveyDate,
         AOI = AreaOfInterest) %>%
  rename_with(.cols = -PermanentID, .fn = ~ paste0("fwc_", .x)) %>%
  mutate(fwc_AOI = str_remove(fwc_AOI, ", Lake") %>%
           str_remove(" Lake"))

# isolate surveys 
fwri_surveys <- fwri2 %>%
  distinct(PermanentID, fwri_Year, fwri_AOI, fwri_Date)
fwc_surveys <- fwc2 %>%
  distinct(PermanentID, fwc_Year, fwc_AOI, fwc_Date)

# overlapping waterbodies
# calculate date difference
# make fwc name same structure as fwri
# require surveys to be within one-half year
survey_overlap <- inner_join(fwri_surveys, fwc_surveys, 
                             relationship = "many-to-many") %>%
  mutate(date_diff = abs(fwri_Date - fwc_Date)) %>%
  filter(date_diff <= 183)


#### troubleshoot spatial differences ####

# are any referring to different areas?
filter(survey_overlap, fwri_AOI != fwc_AOI) %>%
  distinct(PermanentID, fwri_AOI, fwc_AOI) %>%
  data.frame()

# North and South Conway split in FWRI
filter(fwc2, fwc_AOI == "Conway") %>% 
  distinct(fwc_AOI, fwc_WaterbodyAcres, fwc_ShapeArea*247.105)

# Boggy Creek is included as part of East Toho in FWC
filter(fwri2, fwri_AOI == "EastToho") %>% 
  distinct(fwri_X, fwri_Y) %>%
  ggplot(aes(x = fwri_X, y = fwri_Y)) +
  geom_point()
# omitted from FWRI

# Eustis 2 was a survey of East Toho saved under Eustis, but in the same year as another East Toho survey - omit

# Harris = Little Harris?
filter(fwc2, PermanentID == "120024301") %>% distinct(fwc_AOI, fwc_WaterbodyAcres)
# no
filter(fwri2, fwri_AOI == "Harris") %>% 
  distinct(fwri_X, fwri_Y) %>%
  ggplot(aes(x = fwri_X, y = fwri_Y)) +
  geom_point()
# Little Harris is included in the FWRI survey

# two Orange surveys?
filter(fwri2, str_detect(fwri_AOI, "Orange") == T) %>% distinct(fwri_Year, fwri_AOI)
# yes
filter(fwri2, str_detect(fwri_AOI, "Orange") == T) %>% 
  distinct(fwri_AOI, fwri_Year, fwri_X, fwri_Y) %>%
  ggplot(aes(x = fwri_X, y = fwri_Y)) +
  geom_point(size = 0.2) +
  facet_grid(fwri_Year ~ fwri_AOI)
# they're the same shape except in 2017 (Orange covers more than Open Water)
filter(fwri2, str_detect(fwri_AOI, "Orange") == T) %>%
  group_by(fwri_Year, fwri_AOI) %>%
  summarize(n_taxa = n_distinct(fwri_Code),
            .groups = "drop") %>%
  pivot_wider(names_from = fwri_AOI,
              values_from = n_taxa)
# likely same surveys except for 2017 - use Open water for more consistency over time

# remove locations that do not overlap
survey_overlap2 <- survey_overlap %>%
  filter(fwc_AOI != "Boggy Creek" & !(fwri_AOI %in% c("Eustis2", "Orange")))

# combine AOI's within PermID
fwri_survey_overlap <- survey_overlap2 %>%
  distinct(PermanentID, fwri_Year, fwri_AOI, fwri_Date) %>%
  group_by(PermanentID) %>%
  mutate(fwri_AOIs = paste(sort(unique(fwri_AOI)), collapse = ", ")) %>%
  ungroup()

fwc_survey_overlap <- survey_overlap2 %>%
  distinct(PermanentID, fwc_Year, fwc_AOI, fwc_Date) %>%
  group_by(PermanentID) %>%
  mutate(fwc_AOIs = paste(sort(unique(fwc_AOI)), collapse = ", ")) %>%
  ungroup()

# update larger datasets based on overlap
fwri3 <- inner_join(fwri2, fwri_survey_overlap)
fwc3 <- inner_join(fwc2, fwc_survey_overlap)


#### troubleshoot temporal issues ####

# number of days between surveys within one year
fwri_survey_dates <- fwri3 %>%
  distinct(PermanentID, fwri_Year, fwri_AOI, fwri_Date, fwri_X, fwri_Y) %>%
  count(PermanentID, fwri_AOI, fwri_Year, fwri_Date) %>%
  arrange(fwri_AOI, PermanentID, fwri_Year) %>%
  group_by(fwri_AOI, PermanentID, fwri_Year) %>%
  mutate(day_diff = fwri_Date - min(fwri_Date),
         tot_n = sum(n)) %>%
  ungroup()

fwc_survey_dates <- fwc3 %>%
  distinct(PermanentID, fwc_Year, fwc_AOI, fwc_Date) %>%
  arrange(fwc_AOI, PermanentID, fwc_Year) %>%
  group_by(fwc_AOI, PermanentID, fwc_Year) %>%
  mutate(day_diff = fwc_Date - min(fwc_Date)) %>%
  ungroup()

# do any later FWRI surveys add a ton of points?
filter(fwri_survey_dates, day_diff > 7 & (n/tot_n) > 0.1)
# manually checked all surveys within year and compared totals to other years
# all seem okay
# surveys that spanned December-January are all given the December year

# multiple surveys within one year for FWC?
filter(fwc_survey_dates, day_diff > 0)
# no


#### match taxa names ####

# taxa lists
fwri_tax <- distinct(fwri3, fwri_ScientificName) %>%
  rename(TaxonName = fwri_ScientificName)
fwc_tax <- fwc3 %>%
  filter(fwc_IsDetected == "Yes") %>%
  distinct(fwc_TaxonName) %>%
  rename(TaxonName = fwc_TaxonName)

# overlapping taxa
taxon_overlap <- inner_join(fwri_tax, fwc_tax)

# add synonyms to fwc taxa without matches
fwc_syn2 <- anti_join(fwc_tax, taxon_overlap) %>%
  left_join(fwc_syn %>%
              distinct(Taxon, syn_name) %>%
              rename(TaxonName = Taxon,
                     SynName = syn_name))

# are any taxa names in synonym list?
filter(fwc_syn2, TaxonName %in% fwc_syn$syn_name)
# yes, but it's also in the taxa list

# overlap with synonyms
taxon_overlap2 <- fwc_syn2 %>%
  filter(!is.na(SynName)) %>%
  inner_join(fwri_tax %>%
               rename(SynName = TaxonName))
# none

#### start here: look at names manually ####


#### summarize presence/absence ####

# all FWRI records are "present"
# summarize within 
fwri_occ <- fwri3 %>%
  distinct(PermanentID, fwri_AOIs, fwri_Year, fwri_CommonName, fwri_ScientificName) %>%
  rename(Year = fwri_Year,
         TaxonName = fwri_ScientificName)

# FWC records are comprehensive
fwc_occ <- fwc3 %>%
  filter(fwc_IsDetected == "Yes") %>%
  distinct(PermanentID, fwc_AOIs, fwc_Year, fwc_TaxonName) %>%
  rename(Year = fwc_Year,
         TaxonName = fwc_TaxonName)

# combine overlapping taxa
fwri_fwc_occ <- fwri_occ %>%
  inner_join(taxon_overlap) %>%
  mutate(fwri_occ = 1) %>%
  inner_join(fwc_occ %>%
               inner_join(taxon_overlap) %>%
               mutate(fwc_occ = 1)) %>%
  mutate(fwri_occ = replace_na(fwri_occ, 0),
         fwc_occ = replace_na(fwc_occ, 0))

# taxa observed by each survey
fwri_fwc_occ_sum <- fwri_fwc_occ %>%
  group_by(PermanentID, fwri_AOIs, fwc_AOIs, Year) %>%
  summarize(n_taxa = n_distinct(TaxonName),
            present_both = sum(fwri_occ == 1 & fwc_occ == 1),
            present_fwri = sum(fwri_occ == 1 & fwc_occ == 0),
            present_fwc = sum(fwri_occ == 0 & fwc_occ == 1)) %>%
  ungroup()

# differences between surveys
filter(fwri_fwc_occ_sum, present_both != n_taxa)
# occurrences are the same


#### older code ####

# combine fwc and fwri
comb <- fwc %>%
  select(PermanentID, AreaName, Area_ha, GSYear, SurveyDate, CommonName, PropCovered, Surveyor, SurveyorExperience, SurveyorExperienceB, WaterbodyList_ha, WaterbodySum_ha) %>%
  rename_with(.cols = c(AreaName, Area_ha, SurveyDate, PropCovered, Surveyor, SurveyorExperience, SurveyorExperienceB, WaterbodyList_ha, WaterbodySum_ha), 
              ~ paste0("FWC_", .x)) %>%
  inner_join(fwri %>%
               select(PermanentID, AreaName, Area_ha, GSYear, SurveyDate, CommonName, 
                      PropCovered1, PropCovered2, PropCovered3, Nguad_Present) %>%
               filter(!is.na(PropCovered1)) %>%
               rename_with(.cols = c(AreaName, Area_ha, SurveyDate, Nguad_Present, 
                                     PropCovered1, PropCovered2, PropCovered3), 
                           ~ paste0("FWRI_", .x))) %>%
  left_join(qual %>%
              select(PermanentID, GSYear, QualityValue, MonthsSampled)) %>% # current year turbidity
  mutate(DateDiff = as.numeric(difftime(FWC_SurveyDate,
                                        FWRI_SurveyDate, 
                                        units = "days"))) %>%
  pivot_longer(cols = starts_with("FWRI_PropCovered"),
               names_to = "FWRI_PropType",
               names_prefix = "FWRI_PropCovered",
               values_to = "FWRI_PropCovered") %>%
  mutate(FWRI_PropType = fct_recode(as.character(FWRI_PropType),
                                    "present" = "1",
                                    "moderate-to-dense" = "2",
                                    "dense" = "3"),
         FWC_SurveyMonth = 11 - month(FWC_SurveyDate), # max month = 11
         CommonName = tolower(CommonName),
         CoverDiff = FWC_PropCovered - FWRI_PropCovered,
         AbsCoverDiff = abs(CoverDiff))

# data points
comb %>% 
  mutate(lake_year = paste0(PermanentID, GSYear)) %>% 
  summarise(points = n_distinct(lake_year),
            lakes = n_distinct(PermanentID))

# check that areas match
comb %>%
  select(PermanentID, FWRI_Area_ha, FWC_Area_ha) %>%
  unique() %>%
  filter(FWRI_Area_ha != FWC_Area_ha)
# yes

# quality samples
comb %>%
  select(PermanentID, GSYear, MonthsSampled) %>%
  unique() %>%
  mutate(MonthsSampled = replace_na(MonthsSampled, 0)) %>%
  ggplot(aes(x = MonthsSampled)) +
  geom_histogram(binwidth = 1)
# most are at least four

# missing current Secchi data
comb %>%
  filter(is.na(QualityValue)) %>%
  select(PermanentID, FWC_AreaName, FWRI_AreaName) %>%
  unique() # 18 lakes

# initial visualizations
comb %>%
  filter(FWRI_PropType == "present") %>%
  ggplot(aes(x = GSYear, y = FWC_PropCovered, color = PermanentID))+
  geom_line() +
  geom_line(aes(y = FWRI_PropCovered), linetype = "dashed") +
  facet_wrap(~ CommonName, scales = "free") +
  theme_bw() +
  theme(legend.position = "none")

# zero abundance
comb %>%
  filter(FWRI_PropType == "present") %>%
  mutate(FWC_zero = if_else(FWC_PropCovered == 0, 1, 0),
         FWRI_zero = if_else(FWRI_PropCovered == 0, 1, 0)) %>%
  group_by(CommonName) %>%
  summarize(FWC_prop_zero = sum(FWC_zero)/length(FWC_zero),
            FWRI_prop_zero = sum(FWRI_zero)/length(FWRI_zero))
# missing from > 70% FWC: alligator weed, water fern
# missing from > 70% FWRI: cuban bulrush, para grass, wild taro
# keep: hydrilla, torpedograss, water hyacinth and lettuce

# select species with enough data
comb2 <- comb %>%
  filter(CommonName %in% c("hydrilla", "torpedograss", "water hyacinth", "water lettuce"))


#### correlations ####

comb_cor <- comb2 %>%
  group_by(CommonName, FWRI_PropType) %>%
  summarise(Cor = cor.test(FWC_PropCovered, FWRI_PropCovered)$estimate,
            P_value = cor.test(FWC_PropCovered, FWRI_PropCovered)$p.value) %>%
  ungroup() %>%
  mutate(Cor_text = as.character(round_half_up(Cor, digits = 2)),
         Star = if_else(P_value < 0.05, "*", ""),
         Cor_P = paste0(Cor_text, Star))


#### correlation figure ####

pdf("output/fwri_fwc_correlations.pdf", width = 4.5, height = 4)
comb_cor %>%
  filter(CommonName != "water fern") %>%
  ggplot(aes(x = CommonName, y = Cor, 
                     fill = as.factor(FWRI_PropType))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Cor_P), 
            position = position_dodge(width = 0.9),
            vjust = 0,
            size = paper_text_size) +
  scale_fill_viridis_d(name = "FWRI population type") +
  labs(x = "Plant species", y = "Pearson correlation coefficient") +
  def_theme_paper +
  theme(legend.position = c(0.77, 0.83))
dev.off()


#### cover difference figure ####

pdf("output/fwri_fwc_proportion_difference_histogram.pdf", width = 4.5, height = 4)
comb2 %>%
  ggplot(aes(x = CoverDiff)) +
  geom_histogram(binwidth = 0.01) +
  facet_grid(FWRI_PropType ~ CommonName, scales = "free") +
  labs(x = "Annual survey - FWRI survey",
       y = "Lake-year combinations") +
  def_theme_paper
dev.off()


#### process data for regressions ####

# function to process data
mod_dat_filt <- function(Species, PropType){
  
  dat_out <- comb2 %>%
    filter(CommonName == Species & 
             FWRI_PropType == PropType) %>%
    mutate(LakeArea_s = (FWRI_Area_ha - min(FWRI_Area_ha)) / sd(FWRI_Area_ha),
           DateDiff_s = (DateDiff - min(DateDiff)) / sd(DateDiff),
           SurveyorExperience_s = (FWC_SurveyorExperience - mean(FWC_SurveyorExperience)) / sd(FWC_SurveyorExperience),
           SurveyorExperienceB = fct_relevel(FWC_SurveyorExperienceB, "medium", "low"))
  
  dat_out_qual <- dat_out %>%
    filter(!is.na(QualityValue)) %>% # limits dataset to use yearly Secchi
    mutate(Turbidity_s = (max(QualityValue) - QualityValue) / sd(QualityValue))
  
  return(list(dat_out, dat_out_qual))
  
}

# function to fit models
mod_fit <- function(dat_in, qual = F){
  
  # fit model depending on species and water quality data
  if(unique(dat_in$CommonName) == "hydrilla" & qual == F){
    mod <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + FWRI_Nguad_Present + (1|GSYear) + (1|PermanentID), data = dat_in)
  } else if(unique(dat_in$CommonName) == "hydrilla" & qual == T){
    mod <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + FWRI_Nguad_Present + Turbidity_s + (1|GSYear) + (1|PermanentID), data = dat_in)
  } else if(qual == F){
    mod <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID), data = dat_in)
  } else{
    mod <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + Turbidity_s + (1|GSYear) + (1|PermanentID), data = dat_in)
  }
  
  return(mod)
  
}

# filter for hydrilla and each pop type
hydr1 <- mod_dat_filt("hydrilla", "present")[[1]]
hydr2 <- mod_dat_filt("hydrilla", "moderate-to-dense")[[1]]
hydr3 <- mod_dat_filt("hydrilla", "dense")[[1]]

hydr1_qual <- mod_dat_filt("hydrilla", "present")[[2]]
hydr2_qual <- mod_dat_filt("hydrilla", "moderate-to-dense")[[2]]
hydr3_qual <- mod_dat_filt("hydrilla", "dense")[[2]]

# visualize scaled variables
ggplot(hydr1, aes(x = LakeArea_s)) +
  geom_histogram()

ggplot(hydr1, aes(x = DateDiff_s)) +
  geom_histogram()

ggplot(hydr1, aes(x = SurveyorExperience_s)) +
  geom_histogram()

ggplot(hydr1, aes(x = FWC_SurveyMonth)) +
  geom_histogram()

ggplot(hydr1, aes(x = FWRI_Nguad_Present)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = LakeArea_s)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = DateDiff_s)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = SurveyorExperience_s)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = Turbidity_s)) +
  geom_histogram()

# check correlation between variables
cor.test(~ LakeArea_s + DateDiff_s, data = hydr1) # not sig
cor.test(~ LakeArea_s + SurveyorExperience_s, data = hydr1) # 0.2
cor.test(~ SurveyorExperience_s + DateDiff_s, data = hydr1) # 0.2
cor.test(~ LakeArea_s + FWC_SurveyMonth, data = hydr1) # -0.3
cor.test(~ FWC_SurveyMonth + SurveyorExperience_s, data = hydr1) # -0.2
cor.test(~ FWC_SurveyMonth + DateDiff_s, data = hydr1) # -0.7

cor.test(~ LakeArea_s + DateDiff_s, data = hydr1_qual) # not sig
cor.test(~ LakeArea_s + SurveyorExperience_s, data = hydr1_qual) # 0.2
cor.test(~ SurveyorExperience_s + DateDiff_s, data = hydr1_qual) # 0.2
cor.test(~ LakeArea_s + FWC_SurveyMonth, data = hydr1_qual) # -0.3
cor.test(~ FWC_SurveyMonth + SurveyorExperience_s, data = hydr1_qual) # -0.2
cor.test(~ FWC_SurveyMonth + DateDiff_s, data = hydr1_qual) # -0.7
cor.test(~ LakeArea_s + Turbidity_s, data = hydr1_qual) # 0.3
cor.test(~ Turbidity_s + SurveyorExperience_s, data = hydr1_qual) # not sig
cor.test(~ Turbidity_s + DateDiff_s, data = hydr1_qual) # -0.1

# filter other species
wahy1 <- mod_dat_filt("water hyacinth", "present")[[1]]
wahy2 <- mod_dat_filt("water hyacinth", "moderate-to-dense")[[1]]
wahy3 <- mod_dat_filt("water hyacinth", "dense")[[1]]

wale1 <- mod_dat_filt("water lettuce", "present")[[1]]
wale2 <- mod_dat_filt("water lettuce", "moderate-to-dense")[[1]]
wale3 <- mod_dat_filt("water lettuce", "dense")[[1]]

torp1 <- mod_dat_filt("torpedograss", "present")[[1]]
torp2 <- mod_dat_filt("torpedograss", "moderate-to-dense")[[1]]
torp3 <- mod_dat_filt("torpedograss", "dense")[[1]]

# filter for water quality data
wahy1_qual <- mod_dat_filt("water hyacinth", "present")[[2]]
wahy2_qual <- mod_dat_filt("water hyacinth", "moderate-to-dense")[[2]]
wahy3_qual <- mod_dat_filt("water hyacinth", "dense")[[2]]

wale1_qual <- mod_dat_filt("water lettuce", "present")[[2]]
wale2_qual <- mod_dat_filt("water lettuce", "moderate-to-dense")[[2]]
wale3_qual <- mod_dat_filt("water lettuce", "dense")[[2]]

torp1_qual <- mod_dat_filt("torpedograss", "present")[[2]]
torp2_qual <- mod_dat_filt("torpedograss", "moderate-to-dense")[[2]]
torp3_qual <- mod_dat_filt("torpedograss", "dense")[[2]]


#### cover difference glmmTMB regressions ####

# models
# opted for glmmTMB because we can't include time-invariant variables like LakeArea_s
# with a fixed effect model (collinear with lake fixed effect)
# glmmTMB also has flexible families if needed

# fit models
hydr_mod1 <- mod_fit(hydr1)
hydr_mod2 <- mod_fit(hydr2)
hydr_mod3 <- mod_fit(hydr3)

wahy_mod1 <- mod_fit(wahy1)
wahy_mod2 <- mod_fit(wahy2)
wahy_mod3 <- mod_fit(wahy3)

wale_mod1 <- mod_fit(wale1)
wale_mod2 <- mod_fit(wale2)
wale_mod3 <- mod_fit(wale3)

torp_mod1 <- mod_fit(torp1)
torp_mod2 <- mod_fit(torp2)
torp_mod3 <- mod_fit(torp3)

# fit quality models
hydr_mod1_qual <- mod_fit(hydr1_qual, qual = T)
hydr_mod2_qual <- mod_fit(hydr2_qual, qual = T)
hydr_mod3_qual <- mod_fit(hydr3_qual, qual = T)

wahy_mod1_qual <- mod_fit(wahy1_qual, qual = T)
wahy_mod2_qual <- mod_fit(wahy2_qual, qual = T)
wahy_mod3_qual <- mod_fit(wahy3_qual, qual = T)

wale_mod1_qual <- mod_fit(wale1_qual, qual = T)
wale_mod2_qual <- mod_fit(wale2_qual, qual = T)
wale_mod3_qual <- mod_fit(wale3_qual, qual = T)

torp_mod1_qual <- mod_fit(torp1_qual, qual = T)
torp_mod2_qual <- mod_fit(torp2_qual, qual = T)
torp_mod3_qual <- mod_fit(torp3_qual, qual = T)

# review models
summary(hydr_mod1) # surveyor
summary(hydr_mod2) # surveyor
summary(hydr_mod3) # surveyor (marg)

summary(wahy_mod1) # surveyor (marg)
summary(wahy_mod2) # date
summary(wahy_mod3) # date

summary(wale_mod1) # date
summary(wale_mod2)
summary(wale_mod3)

summary(torp_mod1)
summary(torp_mod2) # surveyor
summary(torp_mod3) # surveyor

summary(hydr_mod1_qual) # surveyor, turbidity (marg)
summary(hydr_mod2_qual) # surveyor
summary(hydr_mod3_qual) # surveyor

summary(wahy_mod1_qual) # turbidity and surveyor (marg)
summary(wahy_mod2_qual)
summary(wahy_mod3_qual) # date

summary(wale_mod1_qual) # date
summary(wale_mod2_qual)
summary(wale_mod3_qual)

summary(torp_mod1_qual) # surveyor (marg)
summary(torp_mod2_qual) # surveyor, turbidity (marg)
summary(torp_mod3_qual) # surveyor, turbidity

# check model fits
# all have issues

plot(simulateResiduals(hydr_mod1))
plot(simulateResiduals(hydr_mod2))
plot(simulateResiduals(hydr_mod3))

plot(simulateResiduals(wahy_mod1))
plot(simulateResiduals(wahy_mod2))
plot(simulateResiduals(wahy_mod3))

plot(simulateResiduals(wale_mod1))
plot(simulateResiduals(wale_mod2))
plot(simulateResiduals(wale_mod3))

plot(simulateResiduals(torp_mod1))
plot(simulateResiduals(torp_mod2))
plot(simulateResiduals(torp_mod3))

plot(simulateResiduals(hydr_mod1_qual))
plot(simulateResiduals(hydr_mod2_qual))
plot(simulateResiduals(hydr_mod3_qual))

plot(simulateResiduals(wahy_mod1_qual))
plot(simulateResiduals(wahy_mod2_qual))
plot(simulateResiduals(wahy_mod3_qual))

plot(simulateResiduals(wale_mod1_qual))
plot(simulateResiduals(wale_mod2_qual))
plot(simulateResiduals(wale_mod3_qual))

plot(simulateResiduals(torp_mod1_qual))
plot(simulateResiduals(torp_mod2_qual))
plot(simulateResiduals(torp_mod3_qual))


#### glmmTMB regression coefficient plot ####

# terms
mod_terms <- tibble(term = c("(Intercept)", "DateDiff_s", "LakeArea_s",
                             "SurveyorExperience_s", "FWRI_Nguad_Present",
                             "Turbidity_s"),
                    variable = c("intercept", "days after", "area", "surveyor",
                                 "so. naiad", "turbidity"))
# extract model summary
hydr_qual_tid2 <- tidy(hydr_mod2_qual) %>%
  mutate(model = "hydrilla")

wahy_qual_tid2 <- tidy(wahy_mod2_qual) %>%
  mutate(model = "water hyacinth")

wale_qual_tid2 <- tidy(wale_mod2_qual) %>%
  mutate(model = "water lettuce")

torp_qual_tid2 <- tidy(torp_mod2_qual) %>%
  mutate(model = "torpedograss")

mod_qual_tid <- hydr_qual_tid2 %>%
  full_join(wahy_qual_tid2) %>%
  full_join(wale_qual_tid2) %>%
  full_join(torp_qual_tid2) %>%
  left_join(mod_terms)%>%
  filter(effect == "fixed") %>%
  mutate(term = fct_relevel(variable, "intercept", "area",
                            "days after", "surveyor",
                            "turbidity", "so. naiad"),
         conf.int = std.error * 1.96)

# figure
pdf("output/fwri_fwc_comparison_coefficients.pdf", width = 4, height = 5)
dwplot(mod_qual_tid,
       vline = geom_vline(
         xintercept = 0,
         colour = "grey60",
         linetype = 2)) +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  def_theme_paper +
  theme(legend.title = element_blank(),
        legend.position = c(0.82, 0.82)) +
  xlab(expression(paste("Est. (Annual - FWRI) "%+-%" 95% CI", sep = "")))
dev.off()


#### values for text ####

# survyeor experience
sd(hydr2_qual$DateDiff)
mean(hydr2_qual$FWC_SurveyorExperience)
sd(hydr2_qual$FWC_SurveyorExperience)
sd(hydr2_qual$FWC_Area_ha)
sd(wahy2_qual$FWC_Area_ha)
sd(hydr2_qual$QualityValue)

# model summaries
summary(hydr_mod2_qual)
summary(wahy_mod2_qual)
summary(wale_mod2_qual)
summary(torp_mod2_qual)


#### surveyor experience ####

pdf("output/fwri_fwc_comparison_surveyor_experience.pdf", width = 4.4, height = 4)
ggplot(hydr2_qual, 
       aes(x = FWC_SurveyorExperience, y = CoverDiff, 
           color = PermanentID, shape = as.factor(GSYear))) +
  geom_hline(yintercept = 0) +
  geom_point() +
  scale_color_manual(values = rainbow(n_distinct(hydr2_qual$PermanentID)), guide = "none") +
  scale_shape(name = "Growing\nseason") +
  labs(x = "Surveyor experience (surveys)", y = "Annual survey - FWRI survey") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, -5, -10, -10))
dev.off()

hydr2_qual %>%
  filter(CoverDiff > 0.3) %>%
  select(FWC_AreaName, FWRI_AreaName, PermanentID, GSYear, 
         FWRI_PropCovered, FWC_PropCovered, FWC_Surveyor, FWC_SurveyorExperience)
# two lakes all years
# the lake shapes are very consistent between GIS (used to get FWC area) and FWRI sampling maps

filter(fwri, str_detect(AreaName, "Toho")) %>% 
  select(AreaName, PermanentID) %>% unique()
# two Toho's are separate, which is what we want

ggplot(hydr2, 
       aes(x = FWC_SurveyorExperience, y = CoverDiff, 
           color = PermanentID, shape = as.factor(GSYear))) +
  geom_hline(yintercept = 0) +
  geom_point() +
  scale_color_manual(values = rainbow(n_distinct(hydr1$PermanentID)), guide = "none") +
  scale_shape(name = "Growing\nseason") +
  labs(x = "Surveyor experience (surveys)", y = "Annual survey - FWRI survey") +
  def_theme_paper
# similar pattern when all data are included

hydr2 %>%
  filter(PermanentID %in% c("112039563", "112040787")) %>%
  select(FWC_AreaName, GSYear, FWC_PropCovered, FWRI_PropCovered, FWC_Area_ha, FWC_WaterbodyList_ha, FWC_WaterbodySum_ha)

# area estimates off?
hydr2_cor <- hydr2 %>%
  mutate(FWC_PropCovered_cor = (FWC_PropCovered * FWC_Area_ha) / FWC_WaterbodySum_ha,
         FWC_PropCovered = if_else(!is.na(FWC_PropCovered_cor), FWC_PropCovered_cor, FWC_PropCovered),
         CoverDiff = FWC_PropCovered - FWRI_PropCovered)

# correct lake sizes
hydr2_cor %>%
  ggplot(aes(x = FWC_SurveyorExperience, y = CoverDiff, 
             color = PermanentID, shape = as.factor(GSYear))) +
  geom_hline(yintercept = 0) +
  geom_point() +
  scale_color_manual(values = rainbow(n_distinct(hydr2_cor$PermanentID)), guide = "none") +
  scale_shape(name = "Growing\nseason") +
  labs(x = "Surveyor experience (surveys)", y = "Annual survey - FWRI survey") +
  def_theme_paper +
  theme(legend.position = c(0.1, 0.85))
# same pattern

# model without these high lakes
hydr2b <- hydr2 %>%
  filter(!(PermanentID %in% c("112039563", "112040787")))

hydr_mod2b <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + FWRI_Nguad_Present + (1|GSYear) + (1|PermanentID), data = hydr2b)
summary(hydr_mod2b)
# surveyor experience is marginal

# surveyors
# change to hydr1, hydr3 to see differences
hydr2 %>%
  group_by(FWC_Surveyor) %>%
  mutate(SampleSize = n(),
         FWC_SurveyorExperience = max(FWC_SurveyorExperience)) %>%
  ggplot(aes(x = FWC_Surveyor, y = CoverDiff, color = FWC_SurveyorExperience)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = mean(hydr2$CoverDiff), linetype = "dashed") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  scale_color_viridis_c() +
  def_theme_paper +
  theme(axis.text.x = element_text(size = 8, color="black", angle = 45, hjust = 1))
# most surveyors seem to be capturing hydr2
# the overestimators may be tring to capture hydr1, but they're still off

hydr2 %>%
  group_by(FWC_Surveyor) %>%
  summarize(SampleSize = n(),
            FWC_SurveyorExperience = max(FWC_SurveyorExperience)) %>%
  arrange(FWC_SurveyorExperience)

fwc %>%
  group_by(Surveyor) %>%
  summarize(SurveyorExperience = max(SurveyorExperience)) %>%
  arrange(SurveyorExperience) %>%
  data.frame()

# see if units might be wrong
# they should be in acres
# what if they are in hectares?
hydr2 %>%
  mutate(FWC_PropCovered = ((FWC_PropCovered * FWC_Area_ha)/0.405)/FWC_Area_ha, # convert back to original value
         CoverDiff = FWC_PropCovered - FWRI_PropCovered) %>%
  group_by(FWC_Surveyor) %>%
  mutate(SampleSize = n()) %>%
  ggplot(aes(x = FWC_Surveyor, y = CoverDiff, color = SampleSize)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = mean(hydr2$CoverDiff), linetype = "dashed") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  scale_color_viridis_c() +
  def_theme_paper +
  theme(axis.text.x = element_text(size = 8, color="black", angle = 45, hjust = 1))
# no, makes it worse

# what if they are in square m?
hydr2 %>%
  mutate(FWC_PropCovered = (((FWC_PropCovered * FWC_Area_ha)/0.405) * 1e-4)/FWC_Area_ha, # convert back to original value, then conver to ha
         CoverDiff = FWC_PropCovered - FWRI_PropCovered) %>%
  group_by(FWC_Surveyor) %>%
  mutate(SampleSize = n()) %>%
  ggplot(aes(x = FWC_Surveyor, y = CoverDiff, color = SampleSize)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = mean(hydr2$CoverDiff), linetype = "dashed") +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  scale_color_viridis_c() +
  def_theme_paper +
  theme(axis.text.x = element_text(size = 8, color="black", angle = 45, hjust = 1))
# underestimate

# surveyors across all FWC dataset
fwc %>%
  filter(CommonName == "Hydrilla") %>%
  ggplot(aes(x = Surveyor)) +
  geom_bar() +
  def_theme_paper +
  theme(axis.text.x = element_text(size = 8, color="black", angle = 45, hjust = 1))

fwc %>%
  filter(CommonName == "Hydrilla" & str_detect(Surveyor, "Daniela Alviz") == T)
# only four lakes and all in 2019

# model without surveyor
hydr2c <- hydr2 %>%
  filter(!(FWC_Surveyor == "Ed Harris"))

hydr_mod2c <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + FWRI_Nguad_Present + (1|GSYear) + (1|PermanentID), data = hydr2c)
summary(hydr_mod2c)
# no effect


#### older code ####

#### cover difference feols regressions ####

# PermanentID is a dummy variable for: lake size, bathymetry, terrestrial plant encroachment, biologist
# GSYear is a dummy variable for: funding cycle, management priorities, weather
# decided not to go with this method: we want to explain some of these variables and have the information

# hyrilla models full dataset
hydr_fe_mod1 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + FWRI_Nguad_Present | PermanentID + GSYear, data = hydr1)
summary(hydr_fe_mod1)

hydr_fe_mod2 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + FWRI_Nguad_Present | PermanentID + GSYear, data = hydr2)
summary(hydr_fe_mod2)

hydr_fe_mod3 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + FWRI_Nguad_Present | PermanentID + GSYear, data = hydr3)
summary(hydr_fe_mod3)

# hydrilla models water Turbidity dataset
hydr_qual_fe_mod1 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + FWRI_Nguad_Present + Turbidity_s | PermanentID + GSYear, data = hydr1_qual)
summary(hydr_qual_fe_mod1)

hydr_qual_fe_mod2 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + FWRI_Nguad_Present + Turbidity_s | PermanentID + GSYear, data = hydr2_qual)
summary(hydr_qual_fe_mod2)

hydr_qual_fe_mod3 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + FWRI_Nguad_Present + Turbidity_s | PermanentID + GSYear, data = hydr3_qual)
summary(hydr_qual_fe_mod3)

# water hyacinth models full dataset
wahy_fe_mod1 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s | PermanentID + GSYear , data = wahy1)
summary(wahy_fe_mod1)

wahy_fe_mod2 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s | PermanentID + GSYear , data = wahy2)
summary(wahy_fe_mod2)

wahy_fe_mod3 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s | PermanentID + GSYear , data = wahy3)
summary(wahy_fe_mod3)

# water hyacinth models water Turbidity dataset
wahy_qual_fe_mod1 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + Turbidity_s | PermanentID + GSYear, data = wahy1_qual)
summary(wahy_qual_fe_mod1)

wahy_qual_fe_mod2 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + Turbidity_s | PermanentID + GSYear, data = wahy2_qual)
summary(wahy_qual_fe_mod2)

wahy_qual_fe_mod3 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + Turbidity_s | PermanentID + GSYear, data = wahy3_qual)
summary(wahy_qual_fe_mod3)
# date diff sig (negative)

# water lettuce models full dataset
wale_fe_mod1 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s | PermanentID + GSYear , data = wale1)
summary(wale_fe_mod1)
# date diff sig (positive)

wale_fe_mod2 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s | PermanentID + GSYear , data = wale2)
summary(wale_fe_mod2)

wale_fe_mod3 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s | PermanentID + GSYear , data = wale3)
summary(wale_fe_mod3)

# water lettuce models water Turbidity dataset
wale_qual_fe_mod1 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + Turbidity_s | PermanentID + GSYear, data = wale1_qual)
summary(wale_qual_fe_mod1)
# date diff sig (positive)

wale_qual_fe_mod2 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + Turbidity_s | PermanentID + GSYear, data = wale2_qual)
summary(wale_qual_fe_mod2)

wale_qual_fe_mod3 <- feols(CoverDiff ~ DateDiff_s + SurveyorExperience_s + Turbidity_s | PermanentID + GSYear, data = wale3_qual)
summary(wale_qual_fe_mod3)


#### older code, not currently using ####

#### feols regression coefficient plot ####

# rename coefficients
coef_qual_names <- c("(Intercept)" = "Intercept",
                     "Turbidity_s" = "Turbidity", 
                     "SurveyorExperience_s" = "Surveyor experience",
                     "DateDiff_s" = "Days after")

# panels
hydr_qual_fig <- modelplot(hydr_qual_fe_mod2,
                           coef_map = coef_qual_names,
                           background = list(geom_vline(xintercept = 0, color = "black",
                                                        size = 0.5, linetype = "dashed"))) +
  labs(x = "",
       title = "(A) hydrilla") +
  def_theme_paper +
  theme(legend.position = c(0.3, 0.25))

wahy_qual_fig <- modelplot(wahy_qual_fe_mod2,
                           coef_map = coef_qual_names,
                           background = list(geom_vline(xintercept = 0, color = "black",
                                                        size = 0.5, linetype = "dashed"))) +
  labs(x = expression(paste("Est. (FWC - FWRI) "%+-%" 95% CI", sep = "")),
       title = "(B) water hyacinth") +
  def_theme_paper +
  theme(legend.position = "none",
        axis.text.y = element_blank())

wale_qual_fig <- modelplot(wale_qual_fe_mod2,
                           coef_map = coef_qual_names,
                           background = list(geom_vline(xintercept = 0, color = "black",
                                                        size = 0.5, linetype = "dashed"))) +
  labs(x = "",
       title = "(C) water lettuce") +
  def_theme_paper +
  theme(legend.position = "none",
        axis.text.y = element_blank())

plot_grid(hydr_qual_fig, wahy_qual_fig, wale_qual_fig,
          nrow = 1,
          rel_widths = c(1, 0.6, 0.6))



#### feols regression fixed effects plot ####

fe_hydr_qual_mod <- fixef(hydr_qual_fe_mod2)
summary(fe_hydr_qual_mod)
plot(fe_hydr_qual_mod)



#### correlation with lake area figure ####

# # FWRI present
# area_fig1 <- ggplot(filter(comb, FWRI_PropType == "present"), 
#        aes(FWRI_PropCovered, FWC_PropCovered)) +
#   geom_abline(intercept = 0, slope = 1) +
#   geom_point(alpha = 0.6, aes(color = log(FWRI_Area_ha))) +
#   facet_wrap(~ CommonName, scales = "free") +
#   scale_color_distiller(name = "Lake area\n(log-ha)", 
#                         direction = 1,
#                         palette = "Blues") +
#   labs(x = "FWRI survey: proportion occupied (present)") +
#   def_theme_paper +
#   theme(axis.title.y = element_blank())
# 
# # FWRI moderate
# area_fig2 <- ggplot(filter(comb, FWRI_PropType == "moderate-to-dense"), 
#        aes(FWRI_PropCovered, FWC_PropCovered)) +
#   geom_abline(intercept = 0, slope = 1) +
#   geom_point(alpha = 0.6, aes(color = log(FWRI_Area_ha))) +
#   facet_wrap(~ CommonName, scales = "free") +
#   scale_color_distiller(name = "Lake area\n(log-ha)", 
#                         direction = 1,
#                         palette = "Blues") +
#   labs(x = "FWRI survey: proportion occupied (moderate-to-dense)") +
#   def_theme_paper +
#   theme(legend.position = "none",
#         axis.title.y = element_blank())
# 
# # FWRI dense
# area_fig3 <- ggplot(filter(comb, FWRI_PropType == "dense"), 
#        aes(FWRI_PropCovered, FWC_PropCovered)) +
#   geom_abline(intercept = 0, slope = 1) +
#   geom_point(alpha = 0.6, aes(color = log(FWRI_Area_ha))) +
#   facet_wrap(~ CommonName, scales = "free") +
#   scale_color_distiller(name = "Lake area\n(log-ha)", 
#                         direction = 1,
#                         palette = "Blues") +
#   labs(x = "FWRI survey: proportion occupied (dense)") +
#   def_theme_paper +
#   theme(legend.position = "none",
#         axis.title.y = element_blank())
# 
# # axes titles
# area_figy <- textGrob("Annual survey: proportion occupied",
#                       gp = gpar(fontsize = 10), rot = 90)
# 
# # combine panels
# area_pan <- plot_grid(area_fig1 + theme(legend.position = "none"), 
#                       area_fig2, area_fig3,
#                       nrow = 3)
# 
# # add legend
# area_leg <- get_legend(area_fig1)
# area_fig <- plot_grid(area_pan, area_leg,
#                       nrow = 1,
#                       rel_widths = c(1, 0.13))
# 
# # add axes
# pdf("output/fwri_fwc_correlations_lake_area.pdf", 
#     width = 6.5, height = 6)
# grid.arrange(arrangeGrob(area_fig, 
#                          left = area_figy))
# dev.off()


