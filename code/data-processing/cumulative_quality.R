# function for cumulative quality metrics
qual_lag_fun <- function(dat){
  
  outdat <- dat %>%
    group_by(PermanentID) %>% # annual summary
    summarise(P = mean(TP_ug_L, na.rm = T),
              N = mean(TN_ug_L, na.rm = T),
              chlor = mean(CHL_ug_L, na.rm = T),
              secchi = mean(Secchi, na.rm = T)) %>% # add proportion lake treated for all treatments in a year (can exceed 1)
    ungroup() %>%

  # return
  return(outdat)
}