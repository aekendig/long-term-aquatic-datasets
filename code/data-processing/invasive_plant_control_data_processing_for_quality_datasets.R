# format invasive plant and management datasets for water quality data processing scripts (sourced in those scripts)

#### set-up ####

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
qual_ctrl <- read_csv("intermediate-data/FWC_quality_control_formatted.csv")


#### edit data ####

# combine water hyacinth and lettuce percent covered
floating_cover <- inv_plant %>%
  filter(CommonName %in% c("Water hyacinth", "Water lettuce")) %>%
  group_by(PermanentID, GSYear) %>%
  summarize(across(.cols = ends_with ("PropCovered"), sum),
            InitPercCovered = sum(InitPercCovered),
            SpeciesAcres = sum(SpeciesAcres),
            EstAreaCoveredRaw_ha = sum(EstAreaCoveredRaw_ha)) %>%
  ungroup() %>%
  mutate(across(.cols = ends_with ("PropCovered"), ~ if_else(.x > 1, 1, .x)),
         InitPercCovered = if_else(InitPercCovered > 1, 1, InitPercCovered),
         CommonName = "floating plants")

# add floating cover
inv_plant2 <- inv_plant %>%
  full_join(floating_cover) %>%
  mutate(across(ends_with("PropCovered"), ~ .x * 100), # change proportion to percent
         Lag2APCsq = Lag2AvgPropCovered^2) %>% # square perc covered
  rename_with(str_replace, pattern = "PropCovered", replacement = "PercCovered") %>%
  select(CommonName, PermanentID, GSYear, ends_with("PercCovered"), Lag2APCsq, SpeciesAcres, EstAreaCoveredRaw_ha)

# Species names
inv_taxa <- tibble(Species = c("Hydrilla verticillata", "Floating Plants (Eichhornia and Pistia)", "Panicum repens", "Urochloa mutica", "Cyperus blepharoleptos"),
                   CommonName = c("Hydrilla", "floating plants", "Torpedograss", "Para grass", "Cuban bulrush"))

# make treatment year prior year (because we don't know actual dates pre-2010)
# use "Species" and "unique" to get floating plants combined
qual_ctrl2 <- qual_ctrl %>%
  mutate(TreatmentYear = Year - 1) %>%
  select(Species, PermanentID, TreatmentYear, LastTreatment, RecentTreatment, ends_with("Treated")) %>%
  unique() %>%
  left_join(inv_taxa)


#### identify waterbodies to use ####

# has the plant ever been established?
# by using EstAreaCoveredRaw_ha, there had to be more than one year per permanentID
perm_plant <- inv_plant2 %>%
  group_by(PermanentID, CommonName) %>%
  summarize(Established = as.numeric(sum(EstAreaCoveredRaw_ha) > 0)) %>%
  ungroup() %>%
  filter(Established > 0)

# has the plant ever been treated?
perm_ctrl <- qual_ctrl2 %>%
  group_by(PermanentID, Species) %>%
  summarize(Treatments = as.numeric(sum(Lag1Treated, na.rm = T) > 0)) %>%
  ungroup() %>%
  filter(Treatments > 0) %>%
  left_join(inv_taxa)

# waterbodies that have the species present and been managed at least once
perm_plant_ctrl <- inner_join(perm_plant, perm_ctrl) %>%
  select(PermanentID, CommonName)


