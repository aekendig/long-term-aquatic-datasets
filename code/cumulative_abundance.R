# function for cumulative invasive abundance
# revise to remove specific dataset
inv_lag_fun <- function(GSYear, Lag){
  
  # parameters
  ModYear <- GSYear
  
  # filter dataset
  subdat <- inv_fwc2 %>% # use version with NA values
    filter(GSYear <= ModYear & GSYear >= (ModYear - Lag))
  
  # summarize
  outdat <- subdat %>%
    group_by(PermanentID, TaxonName) %>% # summarize over lag
    summarise(PropCovered = mean(PropCovered)) %>%
    ungroup() %>%
    pivot_wider(names_from = TaxonName,
                values_from = PropCovered) %>%
    rename(Hydrilla = "Hydrilla verticillata",
           WaterLettuce = "Pistia stratiotes",
           WaterHyacinth = "Eichhornia crassipes") %>%
    mutate(GSYear = ModYear,
           Lag = Lag)
  
  # return
  return(outdat)
}