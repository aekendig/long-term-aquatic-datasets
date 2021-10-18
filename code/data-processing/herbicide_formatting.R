herbicide_dataset <- function(ctrl_old, ctrl_new, taxa){
  
  # non-herbicide methods (from herbicide_initial_visualizations)
  non_herb <- c("Mechanical Harvester", 
                "Snagging (tree removal)", 
                "Aquatic Dye (for shading)", 
                "Grass Carp", "Hand Removal", 
                "Mechanical (Other)", 
                "Mechanical Shredder", 
                "Prescribed Fire")
  
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
  
  # new herbicide data
  ctrl_new2 <- ctrl_new %>%
    filter(Species %in% taxa$Species & 
             TotalAcres > 0 & !is.na(ControlMethod) & !(ControlMethod %in% non_herb)) %>% # herbicide control only
    group_by(AreaOfInterestID, PermanentID, Species, BeginDate, TreatmentID, TotalAcres, ShapeArea) %>%  # captures area treated for an event without duplication due to multiple herbicides
    mutate(AreaTreated_ha = TotalAcres * 0.405,
           Area_ha = ShapeArea * 100,
           AreaTreated_ha = case_when(AreaTreated_ha > Area_ha ~ Area_ha, # make full area if it exceeds it
                                      TRUE ~ AreaTreated_ha),
           PropTreated = AreaTreated_ha / Area_ha,
           TreatmentMethod = paste(sort(unique(ControlMethod)), collapse = " + "), # combines control methods
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
  
  # combine herbicide data
  ctrl <- ctrl_old2 %>%
    full_join(ctrl_new2) %>%
    full_join(ctrl_old2 %>% # all possible measurements
                select(PermanentID) %>%
                unique() %>%
                expand_grid(GSYear = min(ctrl_old2$TreatmentYear):max(ctrl_old2$TreatmentYear)) %>% 
                full_join(ctrl_new2 %>%
                            select(PermanentID) %>%
                            unique() %>%
                            expand_grid(GSYear = min(ctrl_new2$TreatmentYear):max(ctrl_new2$TreatmentYear))) %>%
                expand_grid(Species = unique(ctrl_new2$Species))) %>% # one row for each species
    mutate(AreaTreated_ha = replace_na(AreaTreated_ha, 0),
           PropTreated = replace_na(PropTreated, 0)) %>% # no herbicide applied that year
    full_join(ctrl_old2 %>% # all ID's across all years
                select(PermanentID) %>%
                unique() %>%
                full_join(ctrl_new2 %>%
                            select(PermanentID) %>%
                            unique()) %>%
                expand_grid(GSYear = min(ctrl_old2$TreatmentYear):max(ctrl_new2$TreatmentYear)) %>%
                expand_grid(Species = unique(ctrl_new2$Species))) 
  # leave AreaTreated_ha as NA for new rows from last join (lakes in old dataset not in new and vice versa)
  
  return(ctrl)
}

