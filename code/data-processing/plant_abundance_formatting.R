plant_abun_format <- function(dat, taxa){
  
  # load packages
  library(magrittr)
  library(lubridate)
  library(car)
  
  # surveyor experience script
  source("code/data-processing/surveyor_experience.R")
  
  # Lake O data and model values
  source("code/data-processing/okeechobee_growth.R")
  
  # remove duplicates
  source("code/data-processing/remove_abundance_duplicates.R")
  
  # surveyor experience
  surv_out <- surv_exp(dat)
  surveyor <- surv_out[[1]]
  surveyor_bins <- surv_out[[2]]

  # plant abundance dataset
  abun_out <- dat %>% # start with all surveys
    select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, Surveyor) %>%
    unique() %>% # one row per survey
    expand_grid(tibble(TaxonName = taxa$TaxonName)) %>% # one row per species per survey
    full_join(dat %>% # add plant information
                filter(TaxonName %in% taxa$TaxonName) %>%
                select(AreaOfInterest, AreaOfInterestID, PermanentID, ShapeArea, SurveyDate, Surveyor, TaxonName, SpeciesAcres)) %>%
    mutate(SpeciesAcres = replace_na(SpeciesAcres, 0), # cover 0 when it wasn't in a survey
           Area_ha = ShapeArea * 100, # convert lake area from km-squared to hectares
           AreaCovered_ha = SpeciesAcres * 0.405, # convert plant cover from acres to hectares
           SurveyMonth = month(SurveyDate),
           SurveyDay = day(SurveyDate),
           SurveyYear = year(SurveyDate),
           GSYear = case_when(SurveyMonth >= 4 ~ SurveyYear,
                              SurveyMonth < 4 ~ SurveyYear - 1), # assume growing season starts in April
           MonthDay = case_when(SurveyMonth >= 4 ~ as.Date(paste("2020", SurveyMonth, SurveyDay, sep = "-")), # start "year" in April (2020/2021 are arbitrary)
                                SurveyMonth < 4 ~ as.Date(paste("2021", SurveyMonth, SurveyDay, sep = "-")))) %>% # this is for joining dayDat
    left_join(taxa) %>%
    left_join(dayDat) %>% # add standardized days (proportion between April 1 and March 31)
    mutate(AreaChangeSD = lakeO_area_beta1 * (lakeO_area_days-Days) + lakeO_area_beta2 * (lakeO_area_days^2 - Days^2)) %>% # calculate the number of sd's to change to get est. max abundance
    group_by(AreaOfInterestID, TaxonName) %>% # take standard deviation by survey area and species
    mutate(EstAreaCoveredRaw_ha = AreaCovered_ha + AreaChangeSD * sd(AreaCovered_ha)) %>% # calculate est. max abundance, NA if only one value is available
    ungroup() %>%
    left_join(surveyor) # surveyor experience
  
  # remove duplicates
  # summarize by waterbody
  abun_out2 <- abun_out %>%
    nest(data = c(SurveyDate, Surveyor, SurveyorExperience, SpeciesAcres, AreaCovered_ha, SurveyMonth, SurveyDay, SurveyYear, MonthDay, Days, AreaChangeSD, EstAreaCoveredRaw_ha)) %>% # find multiple surveys within area of interest, growing season year, and species
    mutate(newdata = map(data, ~rem_abun_dups(.))) %>% # remove duplicates
    select(-data) %>% # removes 123 rows of data
    unnest(newdata) %>%
    group_by(PermanentID, Area_ha, GSYear, TaxonName, CommonName) %>% # summarize for multiple AOIs in one PermanentID (i.e., waterbody)
    summarise(AreaName = paste(AreaOfInterest, collapse = "/"),
              SurveyDate = max(SurveyDate),
              SurveyorExperience = mean(SurveyorExperience, na.rm = T),
              SpeciesAcres = sum(SpeciesAcres),
              AreaCovered_ha = sum(AreaCovered_ha),
              EstAreaCoveredRaw_ha = sum(EstAreaCoveredRaw_ha)) %>%
    ungroup() %>%
    mutate(EstAreaCovered_ha = case_when(EstAreaCoveredRaw_ha > Area_ha ~ Area_ha, # reduce areas covered to total area
                                         TRUE ~ EstAreaCoveredRaw_ha),
           PropCovered = EstAreaCovered_ha / Area_ha,
           PropCoveredAdj = case_when(PropCovered < 1e-3 ~ 1e-3, # avoid super small values that skew ratios
                                      TRUE ~ PropCovered),
           SpeciesPresent = case_when(SpeciesAcres > 0 ~ 1,
                                      SpeciesAcres == 0 ~ 0),
           SurveyorExperienceB = case_when(SurveyorExperience <= surveyor_bins$low_Max ~ "low",
                                           SurveyorExperience > surveyor_bins$low_Max & SurveyorExperience <= surveyor_bins$medium_Max ~ "medium",
                                           SurveyorExperience > surveyor_bins$medium_Max ~ "high"),
           SurveyorExperienceB =  fct_relevel(SurveyorExperienceB, "high", "medium", "low")) %>%
    full_join(abun_out %>% # add row for every year for each site/species combo (NA's for missing surveys)
                select(PermanentID, TaxonName) %>%
                unique() %>%
                expand_grid(GSYear = min(abun_out$GSYear):max(abun_out$GSYear))) %>%
    group_by(PermanentID, TaxonName) %>%
    arrange(GSYear) %>% 
    mutate(PrevPropCovered = lag(PropCovered), # previous year's abundance
           PrevPropCoveredAdj = lag(PropCoveredAdj),
           PrevAreaCoveredRaw_ha = lag(EstAreaCoveredRaw_ha),
           PrevSpeciesPresent = as.numeric(lag(SpeciesAcres) > 0),
           SurveyDays = as.numeric(SurveyDate - lag(SurveyDate))) %>%
    ungroup() %>%
    mutate(RatioCovered = EstAreaCoveredRaw_ha/PrevAreaCoveredRaw_ha,
           LogRatioCovered = log(RatioCovered),
           LogitPropCovered = logit(PropCovered, adjust = prop_adjust),
           LogitPrevPropCovered = logit(PrevPropCovered, adjust = prop_adjust))
  
  return(abun_out2)
  
}