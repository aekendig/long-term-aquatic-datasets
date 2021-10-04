herbicide_new_dataset <- function(ctrl_new, taxa){
  
  # non-herbicide methods (from herbicide_initial_visualizations)
  non_herb <- c("Mechanical Harvester", 
                "Snagging (tree removal)", 
                "Aquatic Dye (for shading)", 
                "Grass Carp", "Hand Removal", 
                "Mechanical (Other)", 
                "Mechanical Shredder", 
                "Prescribed Fire")
  
  # new herbicide data
  # remove 38 cases in which control method is unknown
  ctrl_new2 <- ctrl_new %>%
    filter(Species %in% taxa$Species & 
             TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>% # herbicide control only
    left_join(herb_type %>%
                select(ControlMethod, ActiveIngredient)) %>% # can add mechanism and whether it's contact - need to be formated below
    group_by(AreaOfInterestID, PermanentID, Species, BeginDate, TreatmentID, TotalAcres, ShapeArea) %>%  # captures area treated for an event without duplication due to multiple herbicides
    mutate(AreaTreated_ha = TotalAcres * 0.405,
           Area_ha = ShapeArea * 100,
           AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha, # make full area if it exceeds it
                                      TRUE ~ AreaTreated_ha),
           PropTreated = AreaTreated_ha / Area_ha,
           TreatmentMethod = paste(sort(unique(ControlMethod)), collapse = " + "), # combines control methods
           ActiveIngredients = paste(sort(unique(ActiveIngredient)), collapse = ","),
           Endothall = if_else(str_detect("Endothall", ActiveIngredients) == T, 1, 0),
           TreatmentYear = year(BeginDate),
           TreatmentMonth = month(BeginDate),
           TreatmentID = as.character(TreatmentID),
           CtrlSet = "new",
           GSYear = case_when(TreatmentMonth >= 4 ~ TreatmentYear,
                              TreatmentMonth < 4 ~ TreatmentYear - 1)) %>%
    ungroup() %>%
    select(AreaOfInterestID, PermanentID, TreatmentYear, Species, Area_ha, AreaTreated_ha, PropTreated, TreatmentMethod, TreatmentMonth, BeginDate, TreatmentID, CtrlSet, GSYear) %>%
    rename(TreatmentDate = BeginDate) %>%
    unique() # some duplication due to different TotalHerbicideUsed (but nothing else) for the same TreatmentID
    
  
  return(ctrl_new2)
}

