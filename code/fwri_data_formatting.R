#### info ####

# goal: import and edit FWRI data
# naming structure:
# camel case for column names
# suffix _FWC or _LW to indicate column source
# units after underscores


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)
library(readxl)


#### manual edits ####

# names of some Excel sheets weren't consistent with others, renamed to "Data" and "Plant Code List"
# started renaming columns (most are 'Site #', some are 'Site'),
# but edited code to deal with this
# Istokpoga 2017-2019 contained two surveys: one for lake and one for west shore, separated these
# 2019 folder contained "Copy of" Baldwin, KingsBay, and Panasoffkee files, removed these
# Orange 2017-2020 contained two surveys: one for open water zone, one for "all", separated these
# EastToho 2018, Hollingsworth 2018 had "1" character values in No_Access column, changed to numbers
# EastToho 2019, calculations in the first two rows, deleted rows
# formatted dates in Istokpoga 2015
# formatted dates in IstokpogaWestShore 2018


#### set-up ####

# years of data
fYears = 6

# initiate datasets
fwri_files <- list()
fwri <- list()
fwri_spp <- list()


#### loop through files ####

# loop through years
for(i in 1:fYears){
  # directory
  year_dir <- paste("original-data/fwri-point-intercept-data/", i + 2014, sep = "")
  
  # list of files for the year
  fwri_files[[i]] <- list.files(path = year_dir)
  
  # initiate dataset
  fwri_dat <- tibble(AOI = NA, Year = NA)
  fwri_spp_dat <- tibble(AOI = NA, Year = NA)
  
  # loop through files
  for(j in 1:length(fwri_files[[i]])){
    
    # import data
    temp_dat <- read_excel(paste(year_dir, fwri_files[[i]][j], sep = "/"), sheet = "Data", guess_max = 10000)
    temp_spp <- read_excel(paste(year_dir, fwri_files[[i]][j], sep = "/"), sheet = "Plant Code List", guess_max = 10000)

    # rename columns if separate characters
    if("Site #" %in% colnames(temp_dat)){
      temp_dat <- rename(temp_dat, "Site" = "Site #")
    }
    
    if("Site#" %in% colnames(temp_dat)){
      temp_dat <- rename(temp_dat, "Site" = "Site#")
    }
    
    if("Depth (ft)" %in% colnames(temp_dat)){
      temp_dat <- rename(temp_dat, "Depth_ft" = "Depth (ft)")
    }
    
    if("No Access" %in% colnames(temp_dat)){
      temp_dat <- rename(temp_dat, "NoAccess" = "No Access")
    }
    
    if("No_Access" %in% colnames(temp_dat)){
      temp_dat <- rename(temp_dat, "NoAccess" = "No_Access")
    }
    
    # move non-date characters to a different column
    if(is.POSIXct(temp_dat$Date) == F){
      temp_dat$OtherDate <- temp_dat$Date
      temp_dat <- temp_dat %>%
        select(-Date)
    }
    
    # change all dates to earliest date
    if("Date" %in% colnames(temp_dat)){
      temp_dat$Date <- min(temp_dat$Date, na.rm = T)
    }
    
    
    # add columns
    temp_dat2 <- temp_dat %>%
      mutate(AOI = str_split(fwri_files[[i]][j], "_")[[1]][1],
             Year = i + 2014,
             Samples = n(),
             Sites = length(unique(Site)))
    
    # add AOI and year to plant code list
    temp_spp2 <- temp_spp %>%
      mutate(AOI = str_split(fwri_files[[i]][j], "_")[[1]][1],
             Year = i + 2014)
    
    # join with previous files
    fwri_dat <- full_join(fwri_dat, temp_dat2)
    fwri_spp_dat <- full_join(fwri_spp_dat, temp_spp2)
  }
  
  # add data to list
  fwri[[i]] <- fwri_dat %>%
    filter(!is.na(AOI))
  
  fwri_spp[[i]] <- fwri_spp_dat %>%
    filter(!is.na(AOI))
}


#### collapse lists ####

# function to rename survey columns and make long
survey_revise_fun <- function(df){
  
  # put lat/long values in X/Y
  if("Latitude" %in% colnames(df)){
    df <- df %>%
      mutate(X = case_when(is.na(X) & !is.na(Longitude) ~ Longitude,
                           TRUE ~ X),
             Y = case_when(is.na(Y) & !is.na(Latitude) ~ Latitude,
                           TRUE ~ Y)) %>%
      select(-c(Latitude, Longitude))
  }
  
  # add OtherDate column
  if(!("OtherDate" %in% colnames(df))){
    df$OtherDate <- NA_character_
  }
  
  # add Depth column
  if(!("Depth_ft" %in% colnames(df))){
    df$Depth_ft <- NA
  }
  
  # add NoAccess column
  if(!("NoAccess" %in% colnames(df))){
    df$NoAccess <- NA
  }
  
  df2 <- df %>%
    select(AOI, Year, Lake, Date, OtherDate, Crew, Site, Y, X, Depth_ft, NoAccess, Samples, Sites, everything()) %>%
    mutate(Plants = select(., 14:ncol(.)) %>%
             rowSums(na.rm = T),
           NoData = as.numeric(Plants == 0)) %>%
    select(-Plants) %>%
    pivot_longer(cols = -c(AOI, Year, Lake, Date, OtherDate, Crew, Site, Y, X, Depth_ft, NoAccess, Samples, Sites),
                 names_to = "Code",
                 values_to = "Abundance") %>%
    filter(Abundance > 0)
  
  return(df2)
  
}

# rename columns and make long
fwri2 <- lapply(fwri, survey_revise_fun)

# collapse list into dataframe
surveys <- bind_rows(fwri2) %>%
  mutate(Lake = case_when(Lake == "Jackson" | is.na(Lake) ~ AOI,
                          TRUE ~ Lake))

# check for duplicate sites
surveys %>%
  filter(Samples != Sites) %>%
  select(AOI, Year, Samples, Sites) %>%
  unique()

# one dataset with duplication (confirmed by looking at Excel file)
ja18 <- surveys %>%
  filter(AOI == "Jackson(Leon)" & Year == 2018) %>%
  mutate(dups = duplicated(.))

nrow(filter(ja18, dups == T))

# remove duplicate rows
surveys2 <- surveys %>%
  unique()

nrow(surveys) - nrow(surveys2) # removed expected

# function to rename plant code list columns
list_rename_fun <- function(df){
  
  df2 <- df %>%
    rename_with(toupper) %>%
    rename("Year" = "YEAR",
           "Code" = "P_CODE",
           "CommonName" = "COMMON_NAME",
           "ScientificName" = "SCIENTIFIC_NAME",
           "FamilyScientific" = "FAMILY_SCIENTIFIC",
           "FamilyCommon" = "FAMILY_COMMON",
           "Type" = "TYPE")
  
  if("ECO_TYPE" %in% colnames(df2)){
    df2 <- df2 %>%
      rename("EcoType" = "ECO_TYPE")
  }
  
  if("FWC_ID" %in% colnames(df2)){
    df2 <- df2 %>%
      rename("FWCID" = "FWC_ID")
  }
  
  return(df2)
  
}

# rename plant code list columns
fwri_spp2 <- lapply(fwri_spp, list_rename_fun)

# collapse list into dataframe
plant_list <- bind_rows(fwri_spp2)

# remove lake/year info from plant list
plant_list2 <- plant_list %>%
  select(Code, ScientificName, CommonName) %>%
  mutate(across(where(is.character), ~na_if(., "NULL"))) %>% # get rid of NULL
  unique()


#### duplicates in plant list ####

# identify duplicates
dup_plants <- plant_list2 %>%
  mutate(dup_code = duplicated(Code),
         dup_sci = duplicated(ScientificName),
         dup_com = duplicated(CommonName)) %>%
  filter(dup_code == T | dup_sci == T | dup_com == T)

# filter for duplicate code
plant_list2 %>%
  filter(Code %in% filter(dup_plants, dup_code == T)$Code) %>%
  arrange(Code) %>%
  data.frame()

# filter for duplicate sci names
plant_list2 %>%
  filter(ScientificName %in% filter(dup_plants, dup_sci == T)$ScientificName) %>%
  arrange(ScientificName) %>%
  data.frame()

# filter for duplicate common names
plant_list2 %>%
  filter(CommonName %in% filter(dup_plants, dup_com == T)$CommonName) %>%
  arrange(CommonName) %>%
  data.frame()

# code change function
code_chg_fun <- function(dat) {
  
  dat_out <- dat %>%
    mutate(Code = case_when(Code == "EELG" ~ "AMEE",
                            Code == "BURE" ~ "AMBU",
                            Code == "BUTT" ~ "BUBU",
                            Code == "CDUC" ~ "CODU",
                            Code == "FWCR" ~ "FLWC",
                            Code %in% c("FYBL", "FLBT") ~ "FLBL",
                            Code == "FMUD" ~ "FLMU",
                            Code == "GIBR" ~ "GIBU",
                            Code == "GSAL" ~ "GISA",
                            Code == "HBLA" ~ "HUBL",
                            Code == "N/A" ~ "NOAC",
                            Code == "REFU" ~ "RUFU",
                            Code == "SSAG" ~ "SAGI",
                            Code == "UMSE" ~ "UMGR",
                            Code == "SALV" ~ "WAFE",
                            Code == "WMIL" ~ "WAMI",
                            TRUE ~ Code))
  
  return(dat_out)
  
}

# revise duplicates
plant_list3 <- plant_list2 %>%
  code_chg_fun %>%
  mutate(CommonName = case_when(Code == "BRPE" ~ "Brazilian peppertree",
                                Code == "ELEA" ~ "elephant-ear/wild taro",
                                Code == "FLBL" ~ "floating bladderwort or FL yellow bladderwort",
                                Code == "LFPR" ~ "large-flower primrose-willow",
                                Code == "MAST" ~ "man-made structure, such as a fishing pier or boat dock",
                                Code == "MELA" ~ "melaleuca",
                                Code == "NYMP" ~ "water-lily species",
                                Code == "NOAC" ~ "no access",
                                Code == "PASP" ~ "paspalum",
                                Code == "POND" ~ "large-leaf pondweed",
                                Code == "RUFU" ~ "rush fuirena",
                                Code == "SAGI" ~ "arrowhead",
                                Code == "SHPR" ~ "shrubby primrose willow",
                                Code == "TSBU" ~ "three-square bulrush",
                                Code == "TUSS" ~ "tussock; 1/more species growing on a buoyant organic mat",
                                Code == "TWMG" ~ "tropical white morning- glory or moonflower",
                                Code == "UMGR" ~ "umbrella-grass",
                                Code == "WAFE" ~ "water fern",
                                TRUE ~ CommonName),
         ScientificName = case_when(Code == "ELEA" ~ "Xanthosoma sagittifolium or Colocasia esculenta",
                                    Code == "FLBL" ~ "Utricularia inflata or Utricularia floridana",
                                    Code == "LFPR" ~ "Ludwigia grandiflora",
                                    Code == "MAST" ~ NA_character_,
                                    Code == "NOAC" ~ NA_character_,
                                    Code == "PANI" ~ "Panicum sp.",
                                    Code == "PASP" ~ "Paspalum fluitans",
                                    Code == "SJWO" ~ "Triadenum virginicum",
                                    Code == "SLSP" ~ "Eleocharis baldwinii",
                                    Code == "SPAT" ~ "Nuphar luteum",
                                    Code == "SWGR" ~ "Luziola fluitans (syn. Hydrochloa caroliniensis)",
                                    ScientificName == "Tussock" ~ NA_character_,
                                    Code == "UMGR" ~ "Fuirena squarrosa or Cyperus articulatus",
                                    Code == "WATE" ~ "Rorippa nasturtium-aquaticum",
                                    TRUE ~ ScientificName)) %>%
  unique()

# check that it worked
plant_list3 %>%
  mutate(dup_code = duplicated(Code),
         dup_sci = duplicated(ScientificName),
         dup_com = duplicated(CommonName)) %>%
  filter(dup_code == T | dup_sci == T | dup_com == T)


#### separate imports ####

# import data
kbfeb <- read_excel("original-data/fwri_KingsBay_PlantDataFebruary2015.xlsx")
kboct <- read_excel("original-data/fwri_KingsBay_PlantDataOctober2015.xlsx")
lf <- read_excel("original-data/fwri_LittleFish_PlantData2018.xlsx")
va <- read_excel("original-data/fwri_Virginia_PlantData2019.xlsx")

# rename columns
# add plant codes
# format no data rows
kbfeb2 <- kbfeb %>%
  rename(CN = "Species Common Name",
         ScientificName = "Species Scientific Name",
         Abundance = Density) %>%
  mutate(AOI = "KingsBay",
         Lake = "Kings Bay",
         Year = 2015,
         ScientificName = str_replace(ScientificName, "species", "sp.")) %>%
  left_join(plant_list3) %>%
  mutate(CommonName = case_when(is.na(CommonName) ~ CN,
                                TRUE ~ CommonName),
         Code = case_when(ScientificName == "No Plant" ~ "NoData",
                          TRUE ~ Code),
         ScientificName = case_when(Code == "NoData" ~ NA_character_,
                                    ScientificName == "Algae Unknown Brown - Phaeophyceae" ~ "Phaeophyceae sp.",
                                    ScientificName == "Algae Unknown Red - Rhodophyta" ~ "Rhodophyta sp.",
                                    ScientificName == "Algae Unknown Green - Chlorophyta" ~ "Chlorophyta sp.",
                                    ScientificName == "Vaucheria" ~ "Vaucheria sp.",
                                    TRUE ~ ScientificName),
         CommonName = case_when(Code == "NoData" ~ NA_character_,
                                TRUE ~ CommonName),
         Abundance = case_when(Code == "NoData" ~ 1,
                               TRUE ~ Abundance),
         Site = paste0(round(X, 6), round(Y, 6)) %>% as.factor() %>% as.numeric()) %>%
  group_by(Site) %>%
  mutate(Data = sum(!is.na(CommonName)) > 0) %>%
  ungroup() %>%
  filter(!(Code == "NoData" & Data == T) | is.na(Code)) %>%
  select(-c(CN, Data)) %>%
  unique() # remove duplicate rows

# check for missing codes
filter(kbfeb2, is.na(Code)) %>%
  select(ScientificName, CommonName) %>% 
  unique()

# rename columns
# add plant codes
# format no data rows
kboct2 <- kboct %>%
  rename(CN = "Species Common Name",
         ScientificName = "Species Scientific Name",
         Abundance = Density) %>%
  mutate(AOI = "KingsBay",
         Lake = "Kings Bay",
         Year = 2015,
         ScientificName = str_replace(ScientificName, "species", "sp."),
         ScientificName = case_when(ScientificName == "Potamogeton pectinatus" ~ "Stuckenia pectinatus",
                                    TRUE ~ ScientificName)) %>%
  left_join(plant_list3) %>%
  mutate(CommonName = case_when(is.na(CommonName) ~ CN,
                                TRUE ~ CommonName),
         Code = case_when(ScientificName == "No Plant" ~ "NoData",
                          TRUE ~ Code),
         ScientificName = case_when(Code == "NoData" ~ NA_character_,
                                    ScientificName == "Algae Unknown Brown - Phaeophyceae" ~ "Phaeophyceae sp.",
                                    ScientificName == "Algae Unknown Green - Chlorophyta" ~ "Chlorophyta sp.",
                                    TRUE ~ ScientificName),
         CommonName = case_when(Code == "NoData" ~ NA_character_,
                                TRUE ~ CommonName),
         Abundance = case_when(Code == "NoData" ~ 1,
                               TRUE ~ Abundance),
         Site = paste0(round(X, 6), round(Y, 6)) %>% as.factor() %>% as.numeric()) %>%
  group_by(Site) %>%
  mutate(Data = sum(!is.na(CommonName)) > 0) %>%
  ungroup() %>%
  filter(!(Code == "NoData" & Data == T) | is.na(Code)) %>%
  select(-c(CN, Data)) %>%
  unique() # remove duplicate rows

# check for missing codes
filter(kboct2, is.na(Code)) %>%
  select(ScientificName, CommonName) %>% 
  unique()

# add columns
# fill in missing dates
# make long
lf2 <- lf  %>%
  rename(NoAccess = No_Access) %>%
  mutate(Date = as.POSIXct("2018-06-29", tz = "America/New_York"),
         Year = 2018,
         AOI = "LittleFish") %>%
  select(AOI, Year, Lake, Date, Crew, Site, Y, X, NoAccess, everything()) %>%
  mutate(Plants = select(., 10:ncol(.)) %>%
           rowSums(na.rm = T),
         NoData = as.numeric(Plants == 0)) %>%
  select(-Plants) %>%
  pivot_longer(cols = -c(AOI, Year, Lake, Date, Crew, Site, Y, X, NoAccess),
               names_to = "Code",
               values_to = "Abundance") %>%
  filter(Abundance > 0) %>%
  code_chg_fun %>%
  left_join(plant_list3) %>%
  unique() # remove duplicate rows

# check for missing names
filter(lf2, is.na(CommonName) & Code != "NoData") %>%
  select(Code, ScientificName, CommonName) %>% 
  unique()

# add columns
# fill in missing dates
# make long
va2 <- va  %>%
  rename(NoAccess = No_Access) %>%
  mutate(Date = as.POSIXct("2019-06-07", tz = "America/New_York"),
         Year = 2019,
         AOI = "Virginia") %>%
  select(AOI, Year, Lake, Date, Crew, Site, Y, X, NoAccess, everything()) %>%
  mutate(Plants = select(., 10:ncol(.)) %>%
           rowSums(na.rm = T),
         NoData = as.numeric(Plants == 0)) %>%
  select(-Plants) %>%
  pivot_longer(cols = -c(AOI, Year, Lake, Date, Crew, Site, Y, X, NoAccess),
               names_to = "Code",
               values_to = "Abundance") %>%
  filter(Abundance > 0) %>%
  code_chg_fun %>%
  left_join(plant_list3) %>%
  unique() # remove duplicate rows

# check for missing names
filter(va2, is.na(CommonName) & Code != "NoData") %>%
  select(Code, ScientificName, CommonName) %>% 
  unique()


#### combine data ####

# add plant list to surveys based on lake and year
surveys3 <- surveys2 %>%
  code_chg_fun %>%
  left_join(plant_list3) %>%
  full_join(kboct2) %>%
  full_join(kbfeb2) %>%
  full_join(lf2) %>%
  full_join(va2)

# check rows
nrow(surveys3)
nrow(surveys2) + nrow(kboct2) + nrow(kbfeb2) + nrow(lf2) + nrow(va2)
unique(surveys3) # same number of rows


#### edit plant codes ####

# identify codes missing info
mis_codes <- surveys3 %>%
  filter(is.na(CommonName)) %>%
  select(Code) %>%
  unique()

# replacement codes
# remove duplicate entries
# manually add codes that were duplicates
rep_codes <- plant_list3 %>%
  filter(Code %in% c("CATA", "DUPO", "REMA", "UMGR", "WAPR")) %>%
  rename(CN = CommonName,
         SN = ScientificName)

# missing - use species name
mis_spp <- plant_list3 %>%
  filter(ScientificName %in% mis_codes$Code) %>%
  rename(Co = Code,
         CN = CommonName,
         SN = ScientificName) %>%
  mutate(Code = SN)

# are any common names?
mis_com <- plant_list3 %>%
  mutate(CN2 = toupper(CommonName)) %>%
  rename(Co = Code) %>%
  inner_join(mis_codes %>%
               mutate(CN2 = toupper(Code))) %>%
  full_join(plant_list3 %>%
              filter(Code %in% c("CAWE", "DOFE", "MOGL", "MOSS", "SAGI"))) %>%
  rename(CN = CommonName,
         SN = ScientificName) %>%
  select(-CN2)

# unknown species
unkn_spp <- surveys3 %>%
  filter(str_detect(Code, "UNK|Unknown") == T) %>%
  select(Code, CommonName, ScientificName) %>%
  unique() %>%
  mutate(num = 1:n(),
         Co = paste0("UNKN", num),
         CN = "unknown plant species") %>%
  select(-num)

# combine lists
code_fix <- rep_codes %>%
  full_join(mis_spp) %>%
  full_join(mis_com) %>%
  full_join(unkn_spp)
  
# remaining missing species
rem_sp <- mis_codes %>%
  anti_join(code_fix %>%
              select(Code)) %>%
  data.frame()


#### add missing info to dataset ####

surveys4 <- surveys3  %>%
  mutate(Code = case_when(str_detect(Code, "CATA") == T ~ "CATA", # duplicates in surveys, added numbers to name, includes an NA_
                          str_detect(Code, "DUPO") == T ~ "DUPO",  # duplicates in surveys, added numbers to name
                          str_detect(Code, "UMGR") == T ~ "UMGR",  # duplicates in surveys, added numbers to name
                          str_detect(Code, "WAPR") == T ~ "WAPR",  # duplicates in surveys, added numbers to name
                          str_detect(Code, "TUSS") == T ~ "TUSS", # tussock
                          Code == "Green Net" ~ "NOAC", # non-plant
                          Code == "Caesars Weed" ~ "CAWE", # added other info through mis_com
                          Code == "Triadica sebifera" ~ "CHTA",
                          Code == "CHAE_Chaetamorpha" ~ "CHAE",
                          ScientificName == "Chlorophyta sp." ~ "CHLO",
                          Code == "Sphagneticola trilobata" ~ "CROX",
                          Code %in% c("DogFennel", "DOGFENNEL", "Dog_Fennel") ~ "DOFE", # added other info through mis_com
                          Code == "EasternCottonwood" ~ "EACO",
                          Code == "ENTE_Enteramorpha" ~ "ENTE",
                          Code %in% c("fern spp.", "Fern spp.") == T ~ "FERN",
                          Code == "Cyperus haspan" ~ "HAFL",
                          Code == "HOPO_HornedPondweed" ~ "HOPO",
                          ScientificName == "Zannichellia palustris" ~ "HOPO",
                          Code == "Hypericum" ~ "HYPE",
                          Code == "Iris spp." ~ "IRIS",
                          Code == "Bay Tree sp." ~ "LAUR",
                          Code == "Lycopus sp." ~ "LYCO",
                          Code %in% c("Acer spp.", "Maple", "MAPLE", "Maple tree") ~ "MATR",
                          Code == "MorningGlory" ~ "MOGL",  # added other info through mis_com
                          Code %in% c("Sphagnum spp.", "Peat moss") ~ "MOSS",
                          str_detect(Code, "OAK|Oak") == T ~ "OATR",
                          Code %in% c("Lygodium microphullum", "Lygodium microphyllum", "OldWorldFern") ~ "OLWO",
                          ScientificName == "Phaeophyceae sp." ~ "PHAE",
                          Code == "RedMaple" ~ "REMA", # added other info through code_fix
                          ScientificName == "Rhodophyta sp." ~ "RHOD",
                          Code %in% c("Submersed Sagittaria", "SubmersedSagittaria") ~ "SAGI",  # added other info through mis_com
                          Code %in% c("Scleria spp.", "SCLERIA", "Scleria") ~ "SCLR",
                          Code == "ScarletSesbania" ~ "SCSE",
                          Code %in% c("SEDGE", "Sedge") ~ "SEDG",
                          Code == "Sea_Lettuc" ~ "SELE",
                          ScientificName == "Ulva lactuca" ~ "SELT",
                          Code == "SHRUB" ~ "SHRU",
                          Code == "Saccarum gigantium" ~ "SUPL",
                          Code == "VirginiaButtonweed" ~ "VIBU",
                          ScientificName == "Vaucheria sp." ~ "YEGA",
                          str_detect(Code, "N/A|NAIsland|NA Spoil|NaWeather|NATreeLineLand|NA_|No_Access") == T ~ "NOAC", # non-plant
                          TRUE ~ Code),
         CommonName = case_when(Code == "CHAE" ~ "Chaetomorpha",
                                Code == "CHTA" ~ "Chinese tallow tree",
                                Code == "CROX" ~ "creeping oxeye",
                                Code == "EACO" ~ "eastern cottonwood",
                                Code == "ENTE" ~ "Enteromorpha",
                                Code == "FERN" ~ "fern",
                                Code == "HAFL" ~ "haspan flatsedge",
                                Code == "HOPO" ~ "horned pondweed",
                                Code == "HYPE" ~ "St. John's wort",
                                Code == "IRIS" ~ "iris",
                                Code == "LAUR" ~ "Bay tree sp.",
                                Code == "LYCO" ~ "waterhorehound",
                                Code == "MATR" ~ "maple tree",
                                Code == "MOSS" ~ "sphagnum moss",
                                Code == "OATR" ~ "oak tree",
                                Code == "OLWO" ~ "Old world climbing fern",
                                Code == "SCLR" ~ "sedge",
                                Code == "SCSE" ~ "scarlet sesban",
                                Code == "SEDG" ~ "sedge",
                                Code == "SELE" ~ "sea lettuce",
                                Code == "SELT" ~ "sea lettuce",
                                Code == "SHRU" ~ "shrub",
                                Code == "SUPL" ~ "sugarcane plumegrass",
                                Code == "TUSS" ~ "tussock; 1/more species growing on a buoyant organic mat",
                                Code == "VIBU" ~ "Virginia buttonweed",
                                Code == "YEGA" ~ "yellow-green algae",
                                TRUE ~ CommonName),
         ScientificName = case_when(Code == "CHAE" ~ "Chaetomorpha",
                                    Code == "CHTA" ~ "Triadica sebifera",
                                    Code == "CROX" ~ "Sphagneticola trilobata",
                                    Code == "EACO" ~ "Populus deltoides",
                                    Code == "ENTE" ~ "Enteromorpha",
                                    Code == "FERN" ~ "fern sp.",
                                    Code == "HAFL" ~ "Cyperus haspan",
                                    Code == "HOPO" ~ "Zannichellia palustris",
                                    Code == "HYPE" ~ "Hypericum sp.",
                                    Code == "IRIS" ~ "Iris sp.",
                                    Code == "LAUR" ~ "Laurus sp.",
                                    Code == "LYCO" ~ "Lycopus sp.",
                                    Code == "MATR" ~ "Acer sp.",
                                    Code == "MOSS" ~ "Sphagnum sp.",
                                    Code == "OATR" ~ "Quercus sp.",
                                    Code == "OLWO" ~ "Lygodium microphyllum",
                                    Code == "SCSE" ~ "Sesbania punicea",
                                    Code == "SCLR" ~ "Scleria sp.",
                                    Code == "SELE" ~ "Ulva sp.",
                                    Code == "SUPL" ~ "Saccharum giganteum",
                                    Code == "VIBU" ~ "Diodia virginiana",
                                    TRUE ~ ScientificName)) %>%
  left_join(code_fix) %>%
  mutate(Code = ifelse(!is.na(Co), Co, Code),
         CommonName = ifelse(is.na(CommonName), CN, CommonName),
         ScientificName = ifelse(is.na(ScientificName), SN, ScientificName)) %>%
  select(-c(CN:Co))

# check updates
surveys4 %>%
  filter(is.na(CommonName)) %>%
  select(Code) %>%
  unique()
# a handful of codes have no species information: SCLE, PARA, JOSP, BALD, MIAD, Moss
# leave as codes to distinguish between other unknowns

# missing codes?
surveys4 %>%
  filter(is.na(Code)) %>%
  select(CommonName, ScientificName) %>%
  unique()

# duplicates?
surveys4 %>%
  select(Code, CommonName, ScientificName) %>%
  unique() %>%
  mutate(dup_code = duplicated(Code),
         dup_sci = duplicated(ScientificName),
         dup_sci = case_when(is.na(ScientificName) ~ F,
                             TRUE ~ dup_sci)) %>%
  filter(dup_code == T | dup_sci == T)


#### fix duplicate codes ####

# duplicates within surveys
site_date_dups <- surveys4 %>%
  group_by(AOI, Year, Date, OtherDate, Site) %>%
  mutate(dup_code = duplicated(Code)) %>%
  ungroup() %>%
  filter(dup_code == T) %>%
  select(AOI, Year, Date, OtherDate, Site, Code) 
# 48 lines, all King's Bay 2015

# choose highest abundance
surveys5 <- surveys4 %>%
  group_by(AOI, Year, Lake, Date, Crew, Site, Y, X, NoAccess, OtherDate, Depth_ft, Code, CommonName, ScientificName) %>%
  summarise(Abundance = max(Abundance)) %>%
  ungroup()

# check that it worked
surveys5 %>%
  group_by(AOI, Year, Date, OtherDate, Site) %>%
  mutate(dup_code = duplicated(Code)) %>%
  ungroup() %>%
  filter(dup_code == T) %>%
  select(AOI, Date, OtherDate, Site, Code) 


#### format abundance ####

# abundance
unique(surveys5$Abundance)

filter(surveys5, Abundance > 3)
# likely type-o's and should be 1

# remove 0 abundance and fix type-o's
surveys6 <- surveys5 %>%
  mutate(Abundance = case_when(Abundance > 3 ~ 1,
                               TRUE ~ Abundance))


#### format date ####

# look at dates
# add dates
surveys6 %>%
  filter(is.na(Date)) %>%
  select(AOI, OtherDate, Year) %>% 
  unique()

# add missing dates
# change timezone
surveys7 <- surveys6 %>%
  mutate(Date = case_when(AOI == "Haines" & Year == 2015 & is.na(Date) ~ as.POSIXct("2015-08-11", tz = "America/New_York"),
                          AOI == "Huntley" & Year == 2019 & is.na(Date) ~ as.POSIXct("2019-08-08", tz = "America/New_York"),
                          TRUE ~ as.POSIXct(as.character(Date), tz = "America/New_York"))) %>%
  select(-OtherDate)


#### format location ####

# duplicate AOI/Lake
dup_loc <- surveys7 %>%
  select(AOI, Lake) %>%
  unique() %>%
  mutate(dup_aoi = duplicated(AOI),
         dup_lake = duplicated(Lake))

dup_loc %>%
  filter(dup_aoi == T) %>%
  select(AOI) %>%
  unique() %>%
  inner_join(dup_loc)

dup_loc %>%
  filter(dup_lake == T) %>%
  select(Lake) %>%
  unique() %>%
  inner_join(dup_loc)
# these are fine

# standardize lake and AOI names
surveys8 <- surveys7 %>%
  mutate(Lake = case_when(AOI == "EastToho" ~ "East Tohopekaliga",
                          AOI == "Jackson(Highlands)" ~ "Jackson (Highlands)",
                          AOI == "Jackson(Leon)" ~ "Jackson (Leon)",
                          AOI == "Jackson(Osceola)" ~ "Jackson (Osceola)",
                          AOI == "JuneInWinter" ~ "June In Winter",
                          AOI == "KingsBay" ~ "Kings Bay",
                          AOI == "LittleConway" ~ "Little Lake Conway",
                          AOI == "NorthConway" ~ "North Conway",
                          AOI == "SouthConway" ~ "South Conway",
                          TRUE ~ Lake))

# AOI's
AOI_ID <- unique(surveys8$AOI)

# figure to check coordinates
pdf("output/FWRI_coordinate_check.pdf")

for(i in 1:length(AOI_ID)){
  
  dat <- surveys8 %>%
    filter(AOI == AOI_ID[i])
  
  print(ggplot(dat, aes(X, Y)) +
          geom_point() +
          facet_wrap(~ Year) +
          ggtitle(dat$AOI[1]))
}

dev.off()

# check lakes with incorrect coordinates
surveys8 %>%
  filter(AOI == "Cannon") %>%
  select(Year, X, Y) %>%
  unique() %>%
  ggplot(aes(X, Y)) +
  geom_point() +
  facet_wrap(~ Year, scales = "free")
# X and Y were probably flipped in 2016

surveys8 %>%
  filter(X > 0 | Y < 0) %>%
  group_by(Year, AOI) %>%
  summarise(X = mean(X),
            Y = mean(Y))
# Virginia is flipped too

surveys8 %>%
  filter(AOI == "Eustis") %>%
  select(Year, X, Y) %>%
  unique() %>%
  ggplot(aes(X, Y)) +
  geom_point() +
  facet_wrap(~ Year, scales = "free")
# 2019 is a different lake

# two lakes on the same shape
surveys8 %>%
  filter((AOI == "Eustis" & Year == 2019) | AOI == "EastToho") %>%
  select(AOI, Year, X, Y) %>%
  unique() %>%
  ggplot(aes(X, Y)) +
  geom_point() +
  facet_wrap(AOI ~ Year)

# compare data for two surveys
surveys8 %>%
  filter((AOI == "Eustis" & Year == 2019) | AOI == "EastToho") %>%
  group_by(AOI, Year) %>%
  summarise(species = sum(!is.na(ScientificName)),
            site = length(unique(Site)))

# update coordinates
surveys9 <- surveys8 %>%
  mutate(X_temp = X,
         X = case_when(X > 0 ~ Y,
                       TRUE ~ X),
         Y = case_when(Y < 0 ~ X_temp,
                       TRUE ~ Y),
         AOI = case_when(AOI == "Eustis" & Year == 2019 ~ "Eustis2",
                         TRUE ~ AOI)) %>%
  select(-X_temp)

# check fixes
surveys9 %>%
  filter(AOI %in% c("Virginia", "Cannon")) %>%
  select(AOI, Year, X, Y) %>%
  unique() %>%
  ggplot(aes(X, Y)) +
  geom_point() +
  facet_grid(AOI ~ Year, scales = "free")

# coordinates
surv_gis <- surveys9 %>%
  group_by(AOI, Lake) %>%
  summarise(longitude = mean(X),
            latitude = mean(Y)) %>%
  ungroup()

# output
write_csv(surv_gis, "gis/data/FWRI_Coordinates.csv")


#### add permanent ID ####

# import data
# see FWRI_lakes_map_methods for notes
surv_gis2 <- read_csv("gis/intermediate-data/FWRI_plants_edited.csv")

# check lakes with different spellings
surveys9 %>%
  filter(str_detect(Lake, "Conway") & str_detect(Lake, "North")) %>%
  select(AOI, Year, X, Y) %>%
  unique() %>%
  ggplot(aes(X, Y)) +
  geom_point() +
  facet_grid(AOI ~ Year, scales = "free")

surveys9 %>%
  filter(str_detect(Lake, "Conway") & str_detect(Lake, "South")) %>%
  select(AOI, Year, X, Y) %>%
  unique() %>%
  ggplot(aes(X, Y)) +
  geom_point() +
  facet_grid(AOI ~ Year, scales = "free")

# fix names
# add permanent ID's
surveys10 <- surveys9 %>%
  left_join(surv_gis2 %>%
              select(AOI, Lake, Permanent_, GNIS_ID, GNIS_Name, AreaSqKm, FType, FCode, ShapeSource)) %>%
  mutate(AOI = case_when(AOI == "Conway(NorthLobe)" ~ "NorthConway",
                         AOI == "Conway(SouthLobe)" ~ "SouthConway",
                         TRUE ~ AOI),
         Lake = case_when(Lake == "Conway (North Lobe)" ~ "North Conway",
                          Lake == "Conway (South Lobe)" ~ "South Conway",
                          TRUE ~ Lake)) %>%
  rename(PermanentID = Permanent_,
         GNISID = GNIS_ID,
         GNISName = GNIS_Name,
         Area_SqKm = AreaSqKm)


#### output ####

write_csv(surveys10, "intermediate-data/FWRI_plant_formatted.csv")

