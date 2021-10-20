#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)

# import data
qual <- read_csv("original-data/Lakewatch_1986-2019_TP_TN_Chl_Secchi_Color_Cond.csv")
gis <- read_csv("intermediate-data/gis_fwc_lakewatch_fwri.csv",
                col_types = list(wkt_geom = col_character(),
                                 AOI = col_character(),
                                 Lake = col_character()))


#### edit data ####

# duplication in GIS
gis %>%
  get_dupes() %>%
  data.frame()
# manually removed duplicate entry of Lake Emily (Lakewatch)

# format date
# make secchi disk data numeric
# add permanent ID based on Lake name and county
# remove missing ID (all have been checked, include rivers, creeks, etc.)
# remove duplicate rows
qual2 <- qual %>%
  mutate(Date = as.Date(Date, "%m/%d/%y"),
         SecchiCombined = ifelse(is.na(SECCHI_ft), SECCHI_2, SECCHI_ft) %>%
           tolower(),
         SecchiBottom = ifelse(str_detect(SecchiCombined, "bottom") == T, 1, 0),
         SecchiWeeds = ifelse(str_detect(SecchiCombined, "weeds") == T, 1, 0),
         Secchi = case_when(SecchiCombined == "." ~ NA_character_, 
                            SecchiCombined == "weeds" ~ NA_character_,
                            SecchiCombined == "weeds (surface)" ~ NA_character_,
                            SecchiCombined == "bottom" ~ NA_character_,
                            TRUE ~ SecchiCombined) %>%
           parse_number(),
         AreaOfInterest = Lake,
         County = toupper(County)) %>%
  rename("Secchi1_ft" = "SECCHI_ft",
         "Secchi2_ft" = "SECCHI_2") %>%
  left_join(gis %>%
              filter(CoordSource == "Lakewatch") %>%
              select(AreaOfInterest, County, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource)) %>%
  filter(!is.na(PermanentID)) %>% # non-lakes or lakes without clear geographic matches
  unique() %>%
  select(-AreaOfInterest)

# summarize quality data by growing season year
qual3 <- qual2 %>%
  mutate(GSYear = case_when(Month >= 4 ~ Year,
                            Month < 4 ~ Year - 1)) %>%
  group_by(PermanentID) %>%
  mutate(across(.cols = c(JoinNotes, ShapeSource, Lake),
                ~ paste(unique(.x), collapse = "/"))) %>% # combine all text info for Perm ID
  ungroup() %>%
  group_by(PermanentID, GSYear, Month) %>% # summarize by month (multiple stations per lake)
  summarise(across(.cols = c(TP_ug_L, TN_ug_L, CHL_ug_L, Secchi, Color_Pt_Co_Units, Cond_uS_cm), 
                   ~ mean(.x, na.rm = TRUE)), # average measurements
            across(.cols = c(Lake, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource), 
                   ~ unique(.x))) %>% # take unique values
  ungroup() %>%
  group_by(PermanentID, GSYear) %>% # summarize by year (multiple months of samples)
  summarise(across(.cols = c(TP_ug_L, TN_ug_L, CHL_ug_L, Secchi, Color_Pt_Co_Units, Cond_uS_cm), 
                   ~ mean(.x, na.rm = TRUE)),
            across(.cols = c(Lake, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource), 
                   ~ unique(.x)),
            n_qual = n()) %>%
  ungroup() 

# rows with all measurements = NA?
qual3 %>%
  filter(across(.cols = c(TP_ug_L, TN_ug_L, CHL_ug_L, Secchi, Color_Pt_Co_Units, Cond_uS_cm),
                ~ is.na(.x)))
# none

# distribution of months
ggplot(qual3, aes(x = n_qual)) +
  geom_histogram()


#### output ####
write_csv(qual3, "intermediate-data/LW_quality_formatted.csv")
