old_ctrl_abundance_dataset <- function(ctrl_dat, abu_dat, taxa){
  
  # modify control data to join invasion data
  # if plant survey occurred before treatment, move treatment to following year
  ctrl_abu <- ctrl_dat %>%
    full_join(taxa) %>%
    left_join(abu_dat %>%
                select(PermanentID, TaxonName, SurveyDate, GSYear) %>% 
                unique()) %>% # add survey dates for each lake and year
    mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
           GSYear = case_when(is.na(SurveyTreatDays) | SurveyTreatDays >= 14 ~ GSYear, # no survey/survey after treatment -> keep year
                              SurveyTreatDays < 14 ~ GSYear + 1)) %>% # survey before or very soon after treatment -> move treatment to next year
    filter(AreaTreated_ha > 0) %>% # treatments applied
    select(PermanentID, TaxonName, Species, GSYear) %>% # annual summary
    unique()
  
  return(ctrl_abu)
}

