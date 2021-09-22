abundance_dataset <- function(dat, taxa){
  
  # surveyor experience
  surveyor <- surv_exp(dat)
  
  # surveyor bins
  surveyor_bins <- surveyor %>%
    select(SurveyorExperienceB) %>%
    unique() %>%
    rowwise() %>%
    mutate(Min = str_remove_all(SurveyorExperienceB, "\\[|\\]|\\(|\\)") %>% 
             str_split(",") %>%
             unlist() %>% first() %>%
             as.numeric(),
           Max = str_remove_all(SurveyorExperienceB, "\\[|\\]|\\(|\\)") %>% 
             str_split(",") %>%
             unlist() %>% last() %>%
             as.numeric())
  
  # categorize surveyor bins
  levels(surveyor_bins$SurveyorExperienceB) <- c("low", "medium", "high")
  
  # make wide
  surveyor_bins %<>% 
    pivot_wider(names_from = SurveyorExperienceB,
                values_from = c(Min, Max),
                names_glue = "{SurveyorExperienceB}_{.value}")
  
  # plant abundance dataset
  abu_out <- dat %>% # start with all surveys
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
    mutate(AreaChangeSD = lakeO_beta1 * (lakeO_days-Days) + lakeO_beta2 * (lakeO_days^2 - Days^2)) %>% # calculate the number of sd's to change to get est. max abundance
    group_by(AreaOfInterestID, TaxonName) %>% # take standard deviation by survey area and species
    mutate(EstAreaCoveredRaw_ha = AreaCovered_ha + AreaChangeSD * sd(AreaCovered_ha)) %>% # calculate est. max abundance, NA if only one value is available
    ungroup() %>%
    left_join(surveyor) # surveyor experience
  
  # remove duplicates
  # summarize by waterbody
  abu_out2 <- abu_out %>%
    nest(data = c(SurveyDate, Surveyor, SurveyorExperience, SpeciesAcres, AreaCovered_ha, SurveyMonth, SurveyDay, SurveyYear, MonthDay, Days, AreaChangeSD, EstAreaCoveredRaw_ha)) %>% # find multiple surveys within area of interest, growing season year, and species
    mutate(newdata = map(data, ~rem_dups_fun(.))) %>% # remove duplicates
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
    full_join(abu_out %>% # add row for every year for each site/species combo (NA's for missing surveys)
                select(PermanentID, TaxonName) %>%
                unique() %>%
                expand_grid(GSYear = min(abu_out$GSYear):max(abu_out$GSYear))) %>%
    group_by(PermanentID, TaxonName) %>%
    arrange(GSYear) %>% 
    mutate(PrevPropCovered = lag(PropCovered), # previous year's abundance
           PrevPropCoveredAdj = lag(PropCoveredAdj),
           PrevAreaCoveredRaw_ha = lag(EstAreaCoveredRaw_ha),
           PrevSpeciesPresent = as.numeric(lag(SpeciesAcres) > 0)) %>%
    ungroup() %>%
    mutate(RatioCovered = EstAreaCoveredRaw_ha/PrevAreaCoveredRaw_ha,
           LogRatioCovered = log(RatioCovered))
  
  return(abu_out2)
  
}