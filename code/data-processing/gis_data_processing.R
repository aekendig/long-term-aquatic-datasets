#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
gis_ed <- read_csv("gis/intermediate-data/Lakewatch_FWC_Waterbody_merge_join_edited.csv")
gis_man <- read_csv("gis/intermediate-data/Lakewatch_FWC_Waterbody_merge_join_manual.csv")
gis_fwri <- read_csv("gis/intermediate-data/FWRI_plants_edited.csv")
fwc_plant_new <- read_csv("original-data/FWCSurvey_2018_2019.csv")
lw_base <- read_csv("original-data/Lakewatch_Base_File_for_Amy_2020.csv")


#### edit data ####

# look at notes
unique(gis_ed$JoinNotes)
unique(gis_man$JoinNotes)

# waterbodies
n_distinct(gis_ed$PermanentID) # 1649
n_distinct(gis_man$PermanentID) # 59

# data structure for outer boundary
gis_ed %>%
  filter(JoinNotes == "outer boundary") %>%
  select(CoordSource, AreaOfInterest, AreaOfInterestID) %>%
  inner_join(gis_ed) %>%
  select(JoinNotes) %>%
  data.frame()

# convert areas to square km
# sum outer boundary and inner lake
outer_boundary_gis <- gis_ed %>%
  filter(JoinNotes == "outer boundary") %>%
  select(AreaOfInterest, AreaOfInterestID, County) %>%
  inner_join(gis_ed %>%
               filter(!is.na(PermanentID))) %>%
  mutate(ShapeArea = ifelse(ShapeSource == "FDEP", ShapeArea*10^-6, ShapeArea*10^4),
         CenterArea = ifelse(JoinNotes != "outer boundary" | is.na(JoinNotes), ShapeArea, NA)) %>%
  group_by(AreaOfInterest, AreaOfInterestID, County) %>%
  summarise(ShapeArea = sum(ShapeArea),
            CenterArea = sum(CenterArea, na.rm = T)) %>%
  ungroup()

# lakes in each dataset
filter(outer_boundary_gis, AreaOfInterestID > 0) # 4 FWC lakes
filter(outer_boundary_gis, AreaOfInterestID == 0) # 19 LW lakes

# match with FWC lakes
outer_boundary_gis %>%
  left_join(fwc_plant_new %>%
              mutate(FWCArea = WaterbodyAcres * 0.004,
                     AreaOfInterest = WaterbodyName) %>%
              select(AreaOfInterest, AreaOfInterestID, FWCArea) %>%
              unique()) %>%
  filter(!is.na(FWCArea))
# one is smaller than LW area, even with outer boundary
# two are more similar to center area

# match with LW lakes
lw_outer_boundary <- outer_boundary_gis %>%
  left_join(lw_base %>%
              mutate(LWArea = Surface_area_hectares * 0.01,
                     AreaOfInterest = Lake) %>%
              select(AreaOfInterest, County, LWArea) %>%
              unique()) %>%
  filter(!is.na(LWArea)) # some LW lakes missing area info
lw_outer_boundary
# LW area tends to be between center area and full area

cor.test(~ LWArea + ShapeArea, data = lw_outer_boundary) # 0.63
cor.test(~ LWArea + CenterArea, data = lw_outer_boundary) # 0.54

# remove lakes without shapes
# add manually matched lakes back in
# convert shape area to km squared
# remove shapes of outer boundary (complicated to combine them with center shapes)
gis <- gis_ed %>%
  filter(!is.na(PermanentID)) %>%
  full_join(gis_man) %>%
  mutate(ShapeArea = ifelse(ShapeSource == "FDEP", ShapeArea*10^-6, ShapeArea*10^4),
         County = toupper(County)) %>%
  filter(JoinNotes != "outer boundary" | is.na(JoinNotes))


#### evaluate data ####

# shape area vs. area_sqkm
ggplot(gis, aes(x = log(ShapeArea), y = log(Area_SqKm))) +
  geom_point()
# mostly well correlated

# compare to FWRI
gis_fwri2 <- gis_fwri %>%
  left_join(gis %>%
               rename_with(~ paste0("FWCLW_", .x)) %>%
               mutate(PermanentID = FWCLW_PermanentID) %>%
              select(PermanentID, FWCLW_ShapeSource, FWCLW_ShapeArea) %>%
              unique()) %>%
  mutate(ShapeSource = case_when(ShapeSource == "Lakewatch_FWC_control_plants_quality" ~ FWCLW_ShapeSource,
                                 ShapeSource == "FL_NHD_Waterbody_merge" ~ "NHD"),
         ShapeArea = ifelse(ShapeSource == "FDEP", ShapeArea*10^-6, ShapeArea*10^4))

gis_fwri2 %>%
  filter(!is.na(FWCLW_ShapeArea) & ShapeArea != FWCLW_ShapeArea) %>%
  select(AOI, ShapeArea, FWCLW_ShapeArea, ShapeSource, FWCLW_ShapeSource)

filter(gis_fwri2, is.na(FWCLW_ShapeSource)) %>%
  select(ShapeArea, FWCLW_ShapeArea)

# add FWRI to gis
gis2 <- gis %>%
  full_join(gis_fwri2 %>%
              select(-starts_with("FWCLW")))


#### outputs ####
write_csv(gis2, "intermediate-data/gis_fwc_lakewatch_fwri.csv")
