non_herb_abundance_dataset <- function(orig_ctrl, unknown_ctrl, non_herb_dat, abu_herb_dat, taxa){
  
  # modify control data to join invasion data
  # if plant survey occurred before treatment, move treatment to following year
  non_herb_abu <- non_herb_dat %>%
    full_join(taxa) %>%
    left_join(abu_herb_dat %>%
                select(PermanentID, TaxonName, SurveyDate, GSYear) %>% 
                unique()) %>% # add survey dates for each lake and year
    mutate(SurveyTreatDays = as.numeric(SurveyDate - TreatmentDate),
           GSYear = case_when(is.na(SurveyTreatDays) | SurveyTreatDays >= 14 ~ GSYear, # no survey/survey after treatment -> keep year
                              SurveyTreatDays < 14 ~ GSYear + 1)) # survey before or very soon after treatment -> move treatment to next year
  
  # fill in zeros
  non_herb_abu2 <- non_herb_abu %>%
    full_join(orig_ctrl %>%
                select(PermanentID) %>%
                unique() %>%
                expand_grid(GSYear = min(orig_ctrl$Year):max(non_herb_abu$GSYear)) %>%
                expand_grid(taxa)) %>% # one row for each species
    mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0),
           PropTreated = replace_na(PropTreated, 0)) # no non-herbicide method used that year
  
  # treatments by year
  non_herb_abu3 <- non_herb_abu2 %>%
    group_by(PermanentID, TaxonName, Species, GSYear) %>% # annual summary
    summarise(PropNonHerbTreated = sum(PropTreated), # add proportion lake treated for all treatments in a year (can exceed 1)
              AreaNonHerbTreated_ha = sum(AreaTreated_ha),
              NonHerbTreatmentDays = case_when(PropNonHerbTreated > 0 ~ as.numeric(n_distinct(TreatmentDate)),
                                               TRUE ~ 0),
              NonHerbTreated = as.numeric(PropNonHerbTreated > 0)) %>% 
    ungroup() %>%
    anti_join(unknown_ctrl) # remove year/lake combos in older dataset, where control type is unknown
  
  # combine non-herb and herb/invasion datasets
  abu_ctrl <- abu_herb_dat %>%
    left_join(non_herb_abu3) # only include data for lakes with surveys
  
  return(abu_ctrl)
  
}

