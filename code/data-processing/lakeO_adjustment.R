#### set-up ####

# load packages
library(tidyverse)
library(readxl)
library(lubridate)

# Lake O data and model values
# source("code/data-processing/okeechobee_growth.R") # can run this if you have dependencies in this script stored
dayDat <- read_csv("intermediate-data/LakeO_day_data_for_model.csv")
lakeO_area_beta1 <- 2.23485
lakeO_area_beta2 <- -1.827721
lakeO_area_days <- -lakeO_area_beta1 / (2 * lakeO_area_beta2)

# import data
dat <- read_excel("intermediate-data/Plant-coverage-Community-analysis.xlsx")


#### edit data ####

dat2 <- dat %>%
  mutate(SurveyMonth = month(SurveyDate), # extract month
         SurveyDay = day(SurveyDate), # extract day
         SurveyYear = year(SurveyDate), # extract year
         MonthDay = case_when(SurveyMonth >= 4 ~ as.Date(paste("2020", SurveyMonth, SurveyDay, sep = "-")), # start "year" in April (2020/2021 are arbitrary)
                              SurveyMonth < 4 ~ as.Date(paste("2021", SurveyMonth, SurveyDay, sep = "-")))) %>% # this is for joining dayDat
  left_join(dayDat) %>%
  mutate(AreaChangeSD = lakeO_area_beta1 * (lakeO_area_days-Days) + lakeO_area_beta2 * (lakeO_area_days^2 - Days^2)) %>% # calculate the number of sd's to change to get est. max abundance
  group_by(County, Lake, Species) %>% # take standard deviation by survey area and species
  mutate(nYears = n_distinct(SurveyYear), # count years to know sample size for sd
         SpeciesAdj_ha = Species_ha + AreaChangeSD * sd(Species_ha, na.rm = T)) %>% # adjust by sd
  ungroup() %>%
  mutate(PAC_raw = 100 * Species_ha / Area_ha,
         PAC_adj = 100 * SpeciesAdj_ha / Area_ha)

#### evaluate conversion ####

# look at sample sizes
dat2 %>%
  distinct(County, Lake, Species, nYears) %>%
  pull(nYears) %>%
  hist()
# may want to consider using raw PAC for lower sd species/lakes

# look at adjustments
ggplot(dat2, aes(x = PAC_raw, y = PAC_adj)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point()
