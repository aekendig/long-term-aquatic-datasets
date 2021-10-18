# function for cumulative herbicide
ctrl_lag_fun <- function(dat){

  outdat <- dat %>%
    group_by(PermanentID, TaxonName, Species, GSYear) %>% # annual summary
    summarise(PropTreated = sum(PropTreated),
              AreaTreated_ha = sum(AreaTreated_ha),
              TreatmentDays = n_distinct(TreatmentDate)) %>% # add proportion lake treated for all treatments in a year (can exceed 1)
    ungroup() %>%
    group_by(PermanentID, TaxonName, Species) %>% # summarize over lag years
    summarise(PropTreated = mean(PropTreated), # average proportion treated per year 
              AreaTreated_ha = mean(AreaTreated_ha),
              Treated = as.numeric(PropTreated > 0),
              TreatmentDays = mean(TreatmentDays)) %>% # average area treated per year
    ungroup()
  
  # return
  return(outdat)
}