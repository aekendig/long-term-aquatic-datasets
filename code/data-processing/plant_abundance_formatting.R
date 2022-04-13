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
  abun_out <- dat %>% 
    filter(TaxonName %in% taxa$TaxonName) %>% # select desired taxa
    select(AreaOfInterest, AreaOfInterestID, County, PermanentID, ShapeArea, WaterbodyAcres, 
           SurveyDate, Surveyor, SurveyMonth, SurveyDay, SurveyYear, GSYear,
           TaxonName, SpeciesAcres) %>%
    mutate(Area_ha = ShapeArea * 100, # convert lake area from km-squared to hectares
           Waterbody_ha = WaterbodyAcres * 0.405, # convert FWC waterbody size to hectares
           AreaCovered_ha = SpeciesAcres * 0.405, # convert plant cover from acres to hectares
           MonthDay = case_when(SurveyMonth >= 4 ~ as.Date(paste("2020", SurveyMonth, SurveyDay, sep = "-")), # start "year" in April (2020/2021 are arbitrary)
                                SurveyMonth < 4 ~ as.Date(paste("2021", SurveyMonth, SurveyDay, sep = "-")))) %>% # this is for joining dayDat
    left_join(taxa) %>% # add common name and taxa code (FWRI)
    left_join(dayDat) %>% # add standardized days (proportion between April 1 and March 31)
    mutate(AreaChangeSD = lakeO_area_beta1 * (lakeO_area_days-Days) + lakeO_area_beta2 * (lakeO_area_days^2 - Days^2)) %>% # calculate the number of sd's to change to get est. max abundance
    group_by(AreaOfInterestID, TaxonName) %>% # take standard deviation by survey area and species
    mutate(EstAreaCoveredRaw_ha = case_when(SpeciesAcres > 0.01 ~ AreaCovered_ha + AreaChangeSD * sd(AreaCovered_ha, na.rm = T), # adjust by sd
                                            TRUE ~ AreaCovered_ha)) %>%
    ungroup() %>%
    left_join(surveyor) # surveyor experience
  
  # minimum area covered
  min_area <- abun_out %>%
    filter(EstAreaCoveredRaw_ha > 0 & !is.na(EstAreaCoveredRaw_ha)) %>%
    pull(EstAreaCoveredRaw_ha) %>%
    min()

  # remove duplicates
  # summarize by waterbody
  # calculate change over one year (add in missing data)
  abun_out2 <- abun_out %>%
    nest(data = c(SurveyDate, Surveyor, SurveyorExperience,
                  SpeciesAcres, AreaCovered_ha,
                  SurveyMonth, SurveyDay, SurveyYear, MonthDay, Days,
                  AreaChangeSD, EstAreaCoveredRaw_ha)) %>% # find multiple surveys within area of interest, growing season year, and species
    mutate(newdata = map(data, ~rem_abun_dups(.))) %>% # chose maximum abundance
    select(-data) %>%
    unnest(newdata) %>%
    group_by(PermanentID, Area_ha, GSYear, TaxonName, CommonName) %>% # summarize for multiple AOIs in one PermanentID (i.e., waterbody)
    summarise(AreaName = paste(sort(unique(AreaOfInterest)), collapse = "/"),
              County = paste(sort(unique(County)), collapse = "/"),
              SurveyDate = max(SurveyDate),
              Surveyor = paste(sort(unique(Surveyor)), collapse = ", "),
              SurveyorExperience = max(SurveyorExperience, na.rm = T), # assign the more experienced surveyor's exp.
              WaterbodyList_ha = paste(sort(unique(Waterbody_ha)), collapse = ", "),
              WaterbodySum_ha = sum(Waterbody_ha, na.rm = T),
              SpeciesAcres = if_else(sum(!is.na(SpeciesAcres)) == 0, NA_real_, # assign NA if all values are NA
                                     sum(SpeciesAcres, na.rm = T)), # sum otherwise
              EstAreaCoveredRaw_ha = if_else(sum(!is.na(EstAreaCoveredRaw_ha)) == 0, NA_real_, # assign NA if all values are NA
                                             sum(EstAreaCoveredRaw_ha, na.rm = T))) %>% # sum otherwise
    ungroup() %>%
    mutate(EstAreaCovered_ha = case_when(EstAreaCoveredRaw_ha > Area_ha ~ Area_ha, # reduce areas covered to total area
                                         TRUE ~ EstAreaCoveredRaw_ha),
           PropCovered = EstAreaCovered_ha / Area_ha,
           SurveyorExperience = ifelse(SurveyorExperience == -Inf, # max returns -Inf if no value is available and na.rm = T
                                       NA_real_,
                                       SurveyorExperience),
           SurveyorExperienceB = case_when(SurveyorExperience <= surveyor_bins$low_Max ~ "low",
                                           SurveyorExperience > surveyor_bins$low_Max & SurveyorExperience <= surveyor_bins$medium_Max ~ "medium",
                                           SurveyorExperience > surveyor_bins$medium_Max ~ "high"),
           SurveyorExperienceB =  fct_relevel(SurveyorExperienceB, "high", "medium", "low")) %>%
    full_join(abun_out %>% # add row for every year for each site/species combo (NA's for missing surveys)
                select(PermanentID, CommonName, TaxonName) %>%
                unique() %>%
                expand_grid(GSYear = min(abun_out$GSYear):max(abun_out$GSYear))) %>%
    group_by(PermanentID, TaxonName) %>%
    arrange(GSYear) %>%
    mutate(InitPercCovered = lag(PropCovered) * 100, # previous year's abundance
           PrevAreaCovered_ha = lag(EstAreaCovered_ha),
           PrevSurveyorExperience = lag(SurveyorExperience),
           SurveyDays = as.numeric(SurveyDate - lag(SurveyDate))) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(MinSurveyorExperience = min(c(SurveyorExperience, PrevSurveyorExperience), na.rm = T)) %>%
    ungroup() %>%
    mutate(RatioCovered = case_when(EstAreaCovered_ha == 0 & PrevAreaCovered_ha == 0 ~ 1,
                                    EstAreaCovered_ha != 0 & PrevAreaCovered_ha == 0 ~ EstAreaCovered_ha / min_area, # lower limit of SpeciesAcres
                                    EstAreaCovered_ha == 0 & PrevAreaCovered_ha != 0 ~ min_area / PrevAreaCovered_ha,
                                    TRUE ~ EstAreaCovered_ha / PrevAreaCovered_ha),
           InitPercCoveredAdj = if_else(PrevAreaCovered_ha == 0, 100 * min_area/Area_ha, InitPercCovered),
           PercCoveredAdj = if_else(EstAreaCovered_ha == 0, 100 * min_area/Area_ha, 100 * PropCovered),
           LogRatioCovered = log(RatioCovered),
           MinSurveyorExperience = ifelse(MinSurveyorExperience == Inf,  # min returns Inf if no value is available and na.rm = T
                                          NA_real_,
                                          MinSurveyorExperience))





  #### lag intervals ####

  # commented out code is to calculate growth over lag interval
  # I found this wasn't useful because the initial size isn't that closely related to the final size

  # function for cumulative treatment
  inv_lag_fun <- function(GSYear, Lag){

    # change name
    year1 <- GSYear

    # filter dataset
    subdat <- abun_out2 %>%
      filter(GSYear <= year1 & GSYear >= (year1 - Lag))

    # summarize
    outdat <- subdat %>%
      group_by(PermanentID, TaxonName) %>%
      # arrange(GSYear) %>%
      # mutate(InitPercCovered = lag(PropCovered, n = Lag) * 100, # previous year's abundance
      #        PrevAreaCovered_ha = lag(EstAreaCovered_ha, n = Lag),
      #        PrevSurveyorExperience = lag(SurveyorExperience, n = Lag),
      #        NYears = n_distinct(GSYear),
      #        AvgPropCovered = mean(PropCovered)) %>%
      summarize(AvgPropCovered = mean(PropCovered),
                NYears = n_distinct(GSYear)) %>%
      ungroup() %>%
      # filter(GSYear == year1) %>%
      # rowwise() %>%
      # mutate(MinSurveyorExperience = min(SurveyorExperience, PrevSurveyorExperience)) %>%
      # ungroup() %>%
      # mutate(RatioCovered = case_when(EstAreaCovered_ha == 0 & PrevAreaCovered_ha == 0 ~ 1,
      #                                 EstAreaCovered_ha != 0 & PrevAreaCovered_ha == 0 ~ EstAreaCovered_ha / min_area, # lower limit of SpeciesAcres
      #                                 EstAreaCovered_ha == 0 & PrevAreaCovered_ha != 0 ~ min_area / PrevAreaCovered_ha,
      #                                 TRUE ~ EstAreaCovered_ha / PrevAreaCovered_ha),
      #        LogRatioCovered = log(RatioCovered)/Lag,
      #        MinSurveyorExperience = ifelse(MinSurveyorExperience == Inf,  # min returns Inf if no value is available and na.rm = T
      #                                       NA_real_,
      #                                       MinSurveyorExperience),
      #        Lag = Lag) %>%
      mutate(Lag = Lag,
             GSYear = year1) # %>%
      # select(PermanentID, TaxonName, GSYear, Lag, NYears, AvgPropCovered)
      # select(PermanentID, TaxonName, GSYear, Lag, NYears, InitPercCovered, MinSurveyorExperience, LogRatioCovered, AvgPropCovered)

    # return
    return(outdat)
  }

  # apply lag
  abun_out3 <- abun_out2 %>%
    select(GSYear) %>%
    unique() %>%
    expand_grid(tibble(Lag = 1:6)) %>% # remove repeat row for each species
    pmap(inv_lag_fun) %>% # summarizes for each GS, Lag, PermID, and Sp
    bind_rows() %>%
    filter(NYears == Lag+1) %>% # all years for a lag must be available
    select(-NYears) %>%
    # pivot_wider(names_from = Lag,
    #             values_from = c(InitPercCovered, MinSurveyorExperience, LogRatioCovered, AvgPropCovered),
    #             names_glue = "Lag{Lag}{.value}") %>% # make treatments wide by lag
    pivot_wider(names_from = Lag,
                values_from = AvgPropCovered,
                names_glue = "Lag{Lag}{.value}") %>% # make treatments wide by lag
    full_join(abun_out2)

  return(abun_out3)
  
}
