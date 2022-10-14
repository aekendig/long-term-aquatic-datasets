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

# dates per year
date_sum <- qual %>%
  group_by(County, Lake, Year) %>%
  summarize(Dates = n_distinct(Date)) %>%
  ungroup() 

ggplot(date_sum, aes(x = Dates)) +
  geom_histogram(binwidth = 1)

date_sum %>%
  summarize(min = min(Dates),
            max = max(Dates))

# waterbodies
qual %>%
  mutate(waterbody = paste(County, Lake)) %>%
  summarize(n = n_distinct(waterbody))

# format date
# make secchi disk data numeric
# remove duplicate rows (none)
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
  unique()

# multiple values per station
qual2_dupes <- get_dupes(qual2, County, Lake, Date, Station)
data.frame(qual2_dupes)
# 6 -- NA rows can be deleted

# make long by water quality metric
# remove NA values
qual3 <- qual2 %>%
  left_join(qual2_dupes %>%
              mutate(remove = if_else(is.na(TP_ug_L), 1, NA_real_))) %>%
  filter(is.na(remove)) %>%
  select(-c(remove, dupe_count, Secchi1_ft, Secchi2_ft, SecchiCombined)) %>%
  rename("Secchi_ft" = "Secchi") %>%
  pivot_longer(cols = c(TP_ug_L, TN_ug_L, CHL_ug_L, 
                        SecchiBottom, SecchiWeeds, Secchi_ft, 
                        Color_Pt_Co_Units, Cond_uS_cm),
               names_to = "QualityMetric",
               values_to = "QualityValue") %>%
  filter(!is.na(QualityValue))

# add permanent ID based on Lake name and county
# remove missing ID (all have been checked, include rivers, creeks, etc.)
qual4 <- qual3 %>%
  left_join(gis %>%
              filter(CoordSource == "Lakewatch") %>%
              select(AreaOfInterest, County, PermanentID, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource)) %>%
  filter(!is.na(PermanentID)) %>% # non-lakes or lakes without clear geographic matches
  select(-AreaOfInterest)

# summarize quality data by month
qual5 <- qual4 %>%
  group_by(PermanentID) %>%
  mutate(across(.cols = c(JoinNotes, ShapeSource, Lake),
                ~ paste(unique(.x), collapse = "/"))) %>% # combine all text info for Perm ID
  ungroup() %>%
  group_by(PermanentID, Year, Month, Date, QualityMetric) %>% # average across stations
  summarise(QualityValue = mean(QualityValue, na.rm = TRUE),
            StationsPerDate = n(),
            across(.cols = c(Lake, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource), 
                   ~ unique(.x))) %>% # take unique values
  ungroup() %>%
  group_by(PermanentID, Year, Month, QualityMetric) %>% # average across dates within a month (some months are highly sampled)
  summarise(QualityValue = mean(QualityValue, na.rm = TRUE),
            AvgStationsPerDate = mean(StationsPerDate),
            DatesSampled = n(),
            across(.cols = c(Lake, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource), 
                   ~ unique(.x))) %>%
  ungroup() 

# visualize station distribution (only part of above code was run)
max(qual$Station) # 6 is the max per Lakewatch-defined lake
ggplot(qual5, aes(x = StationsPerDate)) +
  geom_histogram(binwidth = 1)

# visulaize month distribution
ggplot(qual5, aes(x = Month)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ QualityMetric)
# grouping samples by quarters includes at least one highly sampled month

# year quarters
quarts <- tibble(Month = 1:12,
                 Quarter = rep(c(4, 1:3), each = 3))

# summarize quality data by quarter
qual6 <- qual5 %>%
  left_join(quarts) %>%
  group_by(PermanentID, Year, Quarter, QualityMetric) %>% # average across months in a quarter
  summarize(QualityValue = mean(QualityValue, na.rm = TRUE),
            AvgStationsPerDate = mean(AvgStationsPerDate),
            AvgDatesPerMonth = mean(DatesSampled),
            MonthsSampled = n_distinct(Month),
            across(.cols = c(Lake, GNISID, GNISName, Elevation, FType, FCode, ShapeArea, JoinNotes, ShapeSource), 
                   ~ unique(.x))) %>%
  ungroup()

# distribution of samples
qual6 %>%
  filter(QualityMetric %in% c("CHL_ug_L", "TN_ug_L", "TP_ug_L", "Secchi_ft", "Color_Pt_Co_Units", "Cond_uS_cm")) %>%
  ggplot(aes(x = MonthsSampled)) +
  geom_histogram(binwidth = 1) +
  facet_grid(Quarter ~ QualityMetric, scales = "free")

# some have a ton of dates (used to change summarizing above)
# (many_dates <- qual6 %>%
#   filter(DatesSampled > 24) %>%
#   select(PermanentID, Year, QualityMetric, DatesSampled) %>%
#   left_join(qual4)
# some months just have intense sampling


#### output ####
write_csv(qual6, "intermediate-data/LW_quality_formatted.csv")
