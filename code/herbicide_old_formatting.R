herbicide_old_dataset <- function(ctrl_old, taxa){
  
  # old herbicide data
  ctrl_old2 <- ctrl_old %>%
    filter(Species %in% taxa$Species & TotalAcres > 0) %>%
    mutate(AreaTreated_ha = TotalAcres * 0.405,
           Area_ha = ShapeArea * 100,
           AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha, # make full area if it exceeds it
                                      TRUE ~ AreaTreated_ha),
           PropTreated = AreaTreated_ha / Area_ha,
           TreatmentMethod = "unknown",
           TreatmentID = paste("old", Year, substr(Species, 1, 1), TotalAcres, sep = "_"),
           TreatmentYear = Year,
           CtrlSet = "old",
           TreatmentDate = as.Date(paste0(TreatmentYear, "-07-01")),
           TreatmentMonth = 7,
           GSYear = case_when(TreatmentMonth >= 4 ~ TreatmentYear,
                              TreatmentMonth < 4 ~ TreatmentYear - 1)) %>%
    select(AreaOfInterestID, PermanentID, TreatmentYear, Species, Area_ha, AreaTreated_ha, PropTreated, TreatmentMethod, TreatmentMonth, TreatmentDate, TreatmentID, CtrlSet, GSYear)
  
  # combine herbicide data
  ctrl_old3 <- ctrl_old2 %>%
    full_join(ctrl_old2 %>% # all possible measurements
                select(PermanentID) %>%
                unique() %>%
                expand_grid(GSYear = min(ctrl_old2$TreatmentYear):max(ctrl_old2$TreatmentYear)) %>% 
                expand_grid(Species = unique(ctrl_old2$Species))) %>% # one row for each species
    mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0),
           PropTreated = replace_na(PropTreated, 0)) # no herbicide applied that year
  
  return(ctrl_old3)
}

