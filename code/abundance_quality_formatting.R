# similar to abundance_herbicide_formatting -- could combine these
quality_abundance_dataset <- function(qual_dat, abu_dat){
  
  # modify quality data to join invasion data
  # if plant survey occurred before quality measurement, move measurement to following year
  qual_abu <- qual_dat %>%
    left_join(abu_dat %>%
                select(PermanentID, SurveyDate, GSYear) %>% 
                unique()) %>% # add survey dates for each lake and year
    mutate(SurveyTreatDays = as.numeric(SurveyDate - Date),
           GSYear = case_when(is.na(SurveyTreatDays) | SurveyTreatDays >= 14 ~ GSYear, # no survey/survey after measurement -> keep year
                              SurveyTreatDays < 14 ~ GSYear + 1)) # survey before or very soon after measurement -> move measurement to next year
  
  # year-lag combos
  year_lag <- qual_abu %>%
    select(GSYear) %>%
    unique() %>%  # remove repeat rows
    expand_grid(tibble(Lag = 0:5))
  
  # treatments by year and lag
  qual_abu2 <- map2_dfr(.x = year_lag$GSYear,
                        .y = year_lag$Lag,
                        .f = ~(qual_abu %>%
                                 filter(GSYear <= .x & GSYear >= (.x - .y)) %>%
                                 qual_lag_fun() %>%
                                 mutate(GSYear = .x,
                                        Lag = .y))) %>%
    pivot_wider(names_from = Lag,
                values_from = c(P, N, chlor, secchi),
                names_glue = "Lag{Lag}{.value}") # make treatments wide by lag
  
  # combine treatment and invasion datasets
  abu_qual <- abu_dat %>%
    inner_join(qual_abu2) %>% # only include data from both datasets
    filter(!is.na(PrevPropCovered)) #  missing initial pop abundance
  
  return(abu_qual)
  
}
