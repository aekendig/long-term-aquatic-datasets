herbicide_abundance_dataset <- function(orig_ctrl, unknown_ctrl, herb_dat, abu_dat, taxa){
  
  # modify control data to join invasion data
  # if plant survey occurred before treatment, move treatment to following year
  herb_abu <- herb_dat %>%
    full_join(taxa) %>%
    left_join(abu_dat %>%
                select(PermanentID, TaxonName, SurveyDate, GSYear) %>% 
                unique()) %>% # add survey dates for each lake and year
    mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
           GSYear = case_when(is.na(SurveyTreatDays) | SurveyTreatDays >= 14 ~ GSYear, # no survey/survey after treatment -> keep year
                              SurveyTreatDays < 14 ~ GSYear + 1)) # survey before or very soon after treatment -> move treatment to next year
  
  # include zeros for no herbicide applied
  herb_abu2 <- herb_abu %>%
    full_join(orig_ctrl %>%
                select(PermanentID) %>%
                unique() %>%
                expand_grid(GSYear = min(orig_ctrl$Year):max(herb_abu$GSYear)) %>%
                expand_grid(taxa)) %>% # one row for each species
    mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0),
           PropTreated = replace_na(PropTreated, 0)) # no herbicide applied that year
  
  # treatments by year
  herb_abu3 <- herb_abu2 %>%
    group_by(PermanentID, TaxonName, Species, GSYear) %>% # annual summary
    summarise(PropHerbTreated = sum(PropTreated), # add proportion lake treated for all treatments in a year (can exceed 1)
              AreaHerbTreated_ha = sum(AreaTreated_ha),
              HerbTreatmentDays = case_when(PropHerbTreated > 0 ~ as.numeric(n_distinct(TreatmentDate)),
                                            TRUE ~ 0),
              HerbTreated = as.numeric(PropHerbTreated > 0),
              TreatmentDate = if_else(HerbTreated == 1, max(TreatmentDate), NA_real_)) %>% 
    ungroup() %>%
    anti_join(unknown_ctrl) # remove year/lake combos in older dataset, where control type is unknown
  
  # combine treatment and invasion datasets
  abu_herb <- abu_dat %>%
    left_join(herb_abu3) %>% # only include data for lakes and years when surveys occurred
    mutate(SurveyTreatDays = if_else(HerbTreated == 1, as.numeric(SurveyDate - TreatmentDate), SurveyDays))
  
  return(abu_herb)
  
}

