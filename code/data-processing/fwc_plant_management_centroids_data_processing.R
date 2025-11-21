#### description ####

# combine centroids from FDEP and FL NHD to use in latitudinal analysis

#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(sf)
library(janitor)
library(tigris)

# import
dat_nhd <- st_read("gis/intermediate-data/FL_NHD_FWC_plant_management_centroids.gpkg")
dat_fdep <- st_read("gis/intermediate-data/FDEP_FWC_plant_management_centroids.gpkg")
target_dat <- read_csv("intermediate-data/FWC_plant_management_target_analysis_formatted.csv")
methods_dat <- read_csv("intermediate-data/FWC_plant_management_methods_analysis_formatted.csv")


#### format data ####

# get location info from processed data
lakes <- target_dat %>% 
  distinct(PermanentID, AreaOfInterestID, AreaOfInterest, County, FType) %>% 
  rename(County_target = County) %>% 
  full_join(methods_dat %>% 
              distinct(PermanentID, AreaOfInterestID, AreaOfInterest, 
                       County, FType) %>% 
              rename(County_methods = County))

# duplicates
get_dupes(lakes, AreaOfInterestID)
get_dupes(lakes, PermanentID)

# select and rename columns
# combine data
# add attributes
dat <- dat_nhd %>% 
  select(PermanentI, GNIS_ID,GNIS_Name, FType, AreaSqKm, CoordSourc, ShapeSourc) %>% 
  rename(PermanentID = PermanentI,
         GNISID = GNIS_ID,
         GNISName = GNIS_Name,
         CoordSource = CoordSourc,
         ShapeSource = ShapeSourc) %>% 
  rbind(dat_fdep %>% 
          select(PermanentI, GNIS_ID,GNIS_NAME, FTYPE, AREASQKM, CoordSourc, 
                 ShapeSourc) %>% 
          rename(PermanentID = PermanentI,
                 GNISID = GNIS_ID,
                 GNISName = GNIS_NAME,
                 FType = FTYPE,
                 AreaSqKm = AREASQKM,
                 CoordSource = CoordSourc,
                 ShapeSource = ShapeSourc)) %>% 
  left_join(lakes)

# check for completeness
lakes %>% 
  distinct(AreaOfInterestID) %>% 
  anti_join(dat)

filter(dat, is.na(AreaOfInterestID))


### counties ####

# completely missing counties
filter(dat, is.na(County_target) & is.na(County_methods)) # none

# different counties
filter(dat, County_target != County_methods) # none

# one missing
filter(dat, is.na(County_target) | is.na(County_methods)) # missing from methods

# download counties
counties_fl <- counties("Florida", cb = TRUE)

# visualize
ggplot() +
  geom_sf(data = counties_fl, color="black",
        fill="white", size = 0.25) +
  geom_sf(data = dat, size = 0.25, color = "tomato1")

ggplot() +
  geom_sf(data = counties_fl, color="black",
          fill="white", size = 0.25) +
  geom_sf_text(data = counties_fl, aes(label = NAME), size = 2)

# select attributes needed
counties_fl2 <- counties_fl %>% 
  select(NAME) %>% 
  rename(County_map = NAME)

# get counties from intersection
dat2 <- st_join(dat, counties_fl2, join = st_within) %>% 
  mutate(County_target = str_to_title(County_target))

# check for disagreements
dat2 %>% 
  filter(County_target != County_map) %>% 
  select(AreaOfInterest, County_target, County_map)
# use the map counties -- thes are consistent with the centroids

# on the border with GA
filter(dat2, is.na(County_map)) # use County_target for these two


#### create table ####

# ftypes (https://www.usgs.gov/ngp-standards-and-specifications/national-hydrography-dataset-nhd-data-dictionary-feature-domains)
ftypes <- tibble(FType = c(390, 436, 466),
                 WaterbodyType = c("Lake/Pond", "Reservoir", "Swamp/Marsh"))
# Lafayette Lake in Leon county is the swamp/marsh
# it's shallow, but large (2.85 sq. miles), so okay to call a lake

# convert ftype
# add lat and long
# format county
dat3 <- dat2 %>% 
  left_join(ftypes) %>% 
  mutate(WaterbodyType = if_else(FType == 466 & 
                                   AreaOfInterest == "Lafayette Lake",
                                 "Lake/Pond", WaterbodyType),
         Longitude = st_coordinates(.)[, 1],
         Latitude = st_coordinates(.)[, 2],
         County = if_else(is.na(County_map), County_target, County_map))

# select columns for table
dat_tab <- dat3 %>% 
  select(AreaOfInterestID, AreaOfInterest, WaterbodyType, County, Longitude,
         Latitude) %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  arrange(AreaOfInterestID)

# save
write_csv(dat_tab, "intermediate-data/FWC_lakes_formatted.csv")
