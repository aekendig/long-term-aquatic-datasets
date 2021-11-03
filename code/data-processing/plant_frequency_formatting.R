plant_freq_format <- function(dat, taxa){
  
  # load packages
  library(lubridate)
  library(car)
  
  # Lake O data and model values
  source("code/data-processing/okeechobee_growth.R")
  
  # plant frequency dataset
  # there's no problem with grouping by PermID before making the max cover adjustment because grouped lakes were surveyed within a week of each other
  freq_out <- dat %>% # start with all surveys
    select(AOI, Lake, PermanentID, ShapeArea, Year, YearF, Date, Site) %>%
    unique() %>% # one row per site and survey
    expand_grid(tibble(Code = taxa$Code)) %>% # one row per species per survey and site
    full_join(dat %>%  # add species information
                filter(Code %in% taxa$Code) %>%
                select(AOI, Lake, PermanentID, ShapeArea, Year, YearF, Date, Site, Code, Abundance)) %>%
    filter(!(AOI %in% c("Eustis2", "OrangeOpenWater"))) %>% 
    # Eustis2 is probably EastToho in 2019 based on coordinates, but that lake has a survey that year
    # Have surveys for Orange (includes open water)
    mutate(Area_ha = ShapeArea * 100, # convert lake area from km-squared to hectares
           Abundance = replace_na(Abundance, 0)) %>%  # species cover 0 when it wasn't in a survey 
    group_by(PermanentID, Area_ha, Year, YearF, Code) %>% # summarize across sites, checked above for duplicates within sites - none; Year accounts for surveys that straddle two years (Jan sampling in year before)
    summarise(AreaName = paste(unique(AOI), collapse = "/"), # North and South Conway are the only ones combined
              Sites = n(),
              SitesOcc1 = sum(Abundance > 0),
              SitesOcc2 = sum(Abundance > 1),
              SitesOcc3 = sum(Abundance > 2),
              DateMin = min(Date),
              DateMax = max(Date)) %>%
    ungroup() %>%
    rename(SurveyDate = DateMax,
           SurveyYear = Year) %>%
    mutate(PropCovered1 = SitesOcc1 / Sites,
           PropCovered2 = SitesOcc2 / Sites,
           PropCovered3 = SitesOcc3 / Sites,
           LogitPropCovered1 = logit(PropCovered1, adjust = prop_adjust),
           LogitPropCovered2 = logit(PropCovered2, adjust = prop_adjust),
           LogitPropCovered3 = logit(PropCovered3, adjust = prop_adjust),
           SurveyDays = difftime(SurveyDate, DateMin, units = "days"),
           SurveyDay = day(SurveyDate),
           SurveyMonth = month(SurveyDate),
           GSYear = case_when(SurveyMonth >= 4 ~ year(SurveyDate),
                              SurveyMonth < 4 ~ year(SurveyDate) - 1), # assume growing season starts in April
           MonthDay = case_when(SurveyMonth >= 4 ~ as.Date(paste("2020", SurveyMonth, SurveyDay, sep = "-")), # start "year" in April (2020/2021 are arbitrary)
                                SurveyMonth < 4 ~ as.Date(paste("2021", SurveyMonth, SurveyDay, sep = "-")))) %>% # this is for joining dayDat
    left_join(dayDat) %>% # add standardized days (proportion between April 1 and March 31)
    mutate(PropCoveredChange = lakeO_prop_beta1 * (lakeO_prop_days-Days) + lakeO_prop_beta2 * (lakeO_prop_days^2 - Days^2),
           EstPropCovered1 = PropCovered1 + PropCoveredChange, # calculate the proportion change to get est. max abundance
           EstPropCovered2 = PropCovered2 + PropCoveredChange,
           EstPropCovered3 = PropCovered3 + PropCoveredChange,
           SpeciesPresent = case_when(PropCovered1 > 0 ~ 1,
                                      PropCovered1 == 0 ~ 0)) %>%
    left_join(taxa) %>% # add common names
    select(-DateMin)

  freq_out2 <- freq_out %>%
    full_join(freq_out %>% # add row for every year for each site/species combo (NA's for missing surveys)
                select(PermanentID, Code, CommonName) %>%
                unique() %>%
                expand_grid(GSYear = min(freq_out$GSYear):max(freq_out$GSYear))) %>%
    group_by(PermanentID, Code, CommonName) %>%
    arrange(GSYear) %>%
    mutate(PrevPropCovered1 = lag(PropCovered1), # previous year's abundance
           PrevPropCovered2 = lag(PropCovered2),
           PrevPropCovered3 = lag(PropCovered3),
           PrevSpeciesPresent = lag(SpeciesPresent),
           SurveyDaysLag = as.numeric(SurveyDate - lag(SurveyDate))) %>%
    ungroup() %>%
    mutate(RatioCovered1 = case_when(PropCovered1 == 0 & PrevPropCovered1 == 0 ~ 1,
                                     PropCovered1 != 0 & PrevPropCovered1 == 0 ~ PropCovered1 / (1/Sites),
                                     PropCovered1 == 0 & PrevPropCovered1 != 0 ~ (1/Sites) / PrevPropCovered1,
                                     TRUE ~ PropCovered1/PrevPropCovered1),
           RatioCovered2 = case_when(PropCovered2 == 0 & PrevPropCovered2 == 0 ~ 1,
                                     PropCovered2 != 0 & PrevPropCovered2 == 0 ~ PropCovered2 / (1/Sites),
                                     PropCovered2 == 0 & PrevPropCovered2 != 0 ~ (1/Sites) / PrevPropCovered2,
                                     TRUE ~ PropCovered2/PrevPropCovered2),
           RatioCovered3 = case_when(PropCovered3 == 0 & PrevPropCovered3 == 0 ~ 1,
                                     PropCovered3 != 0 & PrevPropCovered3 == 0 ~ PropCovered3 / (1/Sites),
                                     PropCovered3 == 0 & PrevPropCovered3 != 0 ~ (1/Sites) / PrevPropCovered3,
                                     TRUE ~ PropCovered3/PrevPropCovered3),
           LogRatioCovered1 = log(RatioCovered1),
           LogRatioCovered2 = log(RatioCovered2),
           LogRatioCovered3 = log(RatioCovered3),
           LogitPrevPropCovered1 = logit(PrevPropCovered1, adjust = prop_adjust),
           LogitPrevPropCovered2 = logit(PrevPropCovered2, adjust = prop_adjust),
           LogitPrevPropCovered3 = logit(PrevPropCovered3, adjust = prop_adjust))
  
  return(freq_out2)
  
}