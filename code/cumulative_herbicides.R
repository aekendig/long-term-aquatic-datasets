# function for cumulative herbicide
ctrl_lag_fun <- function(dat){

  outdat <- dat %>%
    group_by(PermanentID, TaxonName, Species, GSYear) %>% # annual summary
    summarise(PropTreated = sum(PropTreated)) %>% # add proportion lake treated for all treatments in a year (can exceed 1)
    ungroup() %>%
    group_by(PermanentID, TaxonName, Species) %>% # summarize over lag years
    summarise(PropTreated = mean(PropTreated), # average proportion treated per year 
              Treated = as.numeric(PropTreated > 0)) %>%
    ungroup()
  
  # return
  return(outdat)
}