herbicide_abundance_dataset <- function(ctrl_dat, abu_dat, taxa){
  
  # modify control data to join invasion data
  # if plant survey occurred before treatment, move treatment to following year
  ctrl_abu <- ctrl_dat %>%
    full_join(taxa) %>%
    left_join(abu_dat %>%
                select(PermanentID, TaxonName, SurveyDate, GSYear) %>% 
                unique()) %>% # add survey dates for each lake and year
    mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
           GSYear = case_when(is.na(SurveyTreatDays) | SurveyTreatDays >= 14 ~ GSYear, # no survey/survey after treatment -> keep year
                              SurveyTreatDays < 14 ~ GSYear + 1)) # survey before or very soon after treatment -> move treatment to next year
  
  # year-lag combos
  year_lag <- ctrl_abu %>%
    select(GSYear) %>%
    unique() %>%  # remove repeat row for each species
    expand_grid(tibble(Lag = 0:5))
  
  # treatments by year and lag
  ctrl_abu2 <- map2_dfr(.x = year_lag$GSYear,
                        .y = year_lag$Lag,
                        .f = ~(ctrl_abu %>%
                                 filter(GSYear <= .x & GSYear >= (.x - .y)) %>%
                                 ctrl_lag_fun() %>%
                                 mutate(GSYear = .x,
                                        Lag = .y))) %>%
    pivot_wider(names_from = Lag,
                values_from = c(PropTreated, Treated, AreaTreated_ha, TreatmentDays),
                names_glue = "Lag{Lag}{.value}") # make treatments wide by lag
  
  # combine treatment and invasion datasets
  abu_ctrl <- abu_dat %>%
    inner_join(ctrl_abu2) %>% # only include data from both datasets
    group_by(TaxonName) %>%
    mutate(PropCoveredBeta = transform01(PropCovered),  # uses sample size within species, leaving out NA's
           PrevPropCoveredBeta = transform01(PrevPropCovered)) %>%
    ungroup() %>%
    filter(!is.na(PrevPropCovered)) #  missing initial pop abundance
  
  return(abu_ctrl)
  
}

