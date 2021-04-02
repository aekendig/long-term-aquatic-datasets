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


#### set-up ####

# years of data
fYears = 5

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
    
    # add AOI and year to survey data
    temp_dat2 <- temp_dat %>%
      mutate(AOI = str_split(fwri_files[[i]][j], "_")[[1]][1],
             Year = i + 2014)
    
    # rename columns if separate characters
    if("Site #" %in% colnames(temp_dat2)){
      temp_dat2 <- rename(temp_dat2, "Site" = "Site #")
    }
    
    if("Site#" %in% colnames(temp_dat2)){
      temp_dat2 <- rename(temp_dat2, "Site" = "Site#")
    }
    
    if("Depth (ft)" %in% colnames(temp_dat2)){
      temp_dat2 <- rename(temp_dat2, "Depth_ft" = "Depth (ft)")
    }
    
    if("No Access" %in% colnames(temp_dat2)){
      temp_dat2 <- rename(temp_dat2, "NoAccess" = "No Access")
    }
    
    if("No_Access" %in% colnames(temp_dat2)){
      temp_dat2 <- rename(temp_dat2, "NoAccess" = "No_Access")
    }
    
    # move non-date characters to a different column
    if(is.POSIXct(temp_dat2$Date) == F){
      temp_dat2$OtherDate <- temp_dat2$Date
      temp_dat2 <- temp_dat2 %>%
        select(-Date)
    }
    
    # change date to median date
    if("Date" %in% colnames(temp_dat2)){
      temp_dat2$Date <- min(temp_dat2$Date, na.rm = T)
    }
    
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
  
  df2 <- df  %>%
    pivot_longer(cols = -c(AOI, Year, Lake, Date, OtherDate, Crew, Site, Y, X, Depth_ft, NoAccess),
                 names_to = "Code",
                 values_to = "Abundance") %>%
    filter(!is.na(Abundance))
  
  return(df2)
  
}

# rename columns and make long
fwri2 <- lapply(fwri, survey_revise_fun)

# collapse list into dataframe
surveys <- bind_rows(fwri2)

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

# add plant list to surveys based on lake and year
surveys2 <- surveys %>%
  left_join(plant_list)

# remove lake/year info from plant list
plant_list2 <- plant_list %>%
  select(-c(AOI, Year)) %>%
  unique()


#### edit plant codes ####

# identify codes missing info
mis_codes <- surveys2 %>%
  filter(is.na(CommonName)) %>%
  select(Code) %>%
  unique()

# find codes in plant list
# remove duplicate entries
# manually add codes that were duplicates
rep_codes <- plant_list2 %>%
  inner_join(mis_codes) %>%
  arrange(Code, FWCID) %>%
  mutate(dupCode = duplicated(Code)) %>%
  filter(dupCode == F) %>%
  select(-dupCode) %>%
  rename(CN = CommonName,
         SN = ScientificName,
         FS = FamilyScientific,
         FC = FamilyCommon,
         Ty = Type,
         FI = FWCID,
         ET = EcoType) %>%
  full_join(plant_list2 %>%
              filter(Code %in% c("CATA", "DUPO","UMGR", "WAPR") & !(is.na(FWCID))))

# are any species names?
mis_spp <- plant_list2 %>%
  filter(ScientificName %in% mis_codes$Code) %>%
  rename(Co = Code,
         CN = CommonName,
         SN = ScientificName,
         FS = FamilyScientific,
         FC = FamilyCommon,
         Ty = Type,
         FI = FWCID,
         ET = EcoType) %>%
  mutate(Code = SN)

# are any common names?
mis_com <- plant_list2 %>%
  mutate(CN2 = toupper(CommonName)) %>%
  rename(Co = Code) %>%
  inner_join(mis_codes %>%
               mutate(CN2 = toupper(Code))) %>%
  mutate(dups = duplicated(Co)) %>%
  filter(dups == F) %>%
  select(-c(dups, CN2)) %>%
  rename(CN = CommonName,
         SN = ScientificName,
         FS = FamilyScientific,
         FC = FamilyCommon,
         Ty = Type,
         FI = FWCID,
         ET = EcoType)

# remaining missing species
mis_codes %>%
  anti_join(rep_codes %>%
              select(Code)) %>%
  anti_join(mis_spp %>%
              select(Code)) %>%
  anti_join(mis_com %>%
              select(Code)) %>%
  data.frame()
#### start here: on 18 ####


#### add missing info to dataset ####
surveys3 <- surveys2 %>%
  mutate(Code = case_when(str_detect(Code, "UNK|Unknown") == T ~ "unknown",
                          str_detect(Code, "N/A|NAIsland|NA Spoil|NaWeather|NATreeLineLand|NA_") == T ~ "no access",
                          Code == "Triadica sebifera" ~ "CHTA",
                          str_detect(Code, "CATA") == T ~ "CATA",
                          Code %in% c("CHAE_Chaetamorpha", "Chaetomorpha") ~ "CHAE",
                          Code == "Sphagneticola trilobata" ~ "CROX",
                          str_detect(Code, "DUPO") == T ~ "DUPO",
                          Code %in% c("Enteromorpha", "ENTE_Enteramorpha") ~ "ENTE",
                          Code %in% c("fern spp.", "Fern spp.") == T ~ "FERN",
                          Code == "Cyperus haspan" ~ "HAFL",
                          str_detect(Code, "HOPO") == T ~ "HOPO",
                          Code == "Lycopus sp." ~ "LYCO",
                          Code %in% c("Acer spp.", "MAPLE", "Maple tree", "Red Maple") ~ "MATR",
                          str_detect(Code, "OAK|Oak") == T ~ "OATR",
                          Code %in% c("Lygodium microphullum", "Lygodium microphyllum", "OldWorldFern") ~ "OLWO",
                          Code == "Ludwigia peruviana" ~ "PEPR",
                          Code %in% c("Scleria spp.", "SCLERIA") ~ "SCLE",
                          Code %in% c("SEDGE", "Sedge") ~ "SEDG",
                          Code == "SHRUB" ~ "SHRU",
                          str_detect(Code, "TUSS") == T ~ "TUSS",
                          str_detect(Code, "UMGR") == T ~ "UMGR",
                          str_detect(Code, "WAPR") == T ~ "WAPR",
                          TRUE ~ Code),
         CommonName = case_when(Code == "CHAE" ~ "Chaetomorpha sp.",
                                Code == "CHTA" ~ "Chinese tallow tree",
                                Code == "CROX" ~ "creeping oxeye",
                                Code == "ENTE" ~ "Enteromorpha sp.",
                                Code == "FERN" ~ "fern",
                                Code == "HAFL" ~ "haspan flatsedge",
                                Code == "HOPO" ~ "horned pondweed",
                                Code == "LYCO" ~ "Lycopus sp.",
                                Code == "MATR" ~ "maple tree",
                                Code == "OATR" ~ "oak tree",
                                Code == "OLWO" ~ "Old world climbing fern",
                                Code == "PEPR" ~ "Peruvian primrose-willow",
                                Code == "SCLE" ~ "Scleria sp."
                                Code == "SEDG" ~ "sedge",
                                Code == "SHRU" ~ "shrub",
                                TRUE ~ CommonName),
         ScientificName = case_when(Code == "CHTA" ~ "Triadica sebifera",
                                    Code == "CROX" ~ "Sphagneticola trilobata",
                                    Code == "HAFL" ~ "Cyperus haspan",
                                    Code == "HOPO" ~ "Zannichellia palustris",
                                    Code == "OLWO" ~ "Lygodium microphyllum",
                                    Code == "PEPR" ~ "Ludwigia peruviana",
                                    TRUE ~ ScientificName)) %>%
  left_join(rep_codes) %>%
  mutate(Code = ifelse(!is.na(Co), Co, Code),
         CommonName = ifelse(is.na(CommonName), CN, CommonName),
         ScientificName = ifelse(is.na(ScientificName), SN, ScientificName),
         FamilyScientific = ifelse(is.na(FamilyScientific), FS, FamilyScientific),
         FamilyCommon = ifelse(is.na(FamilyCommon), FC, FamilyCommon),
         Type = ifelse(is.na(Type), Ty, Type),
         FWCID = ifelse(is.na(FWCID), FI, FWCID),
         EcoType = ifelse(is.na(EcoType), ET, EcoType))


# identify duplicate codes and names
dup_codes <- surveys2 %>%
  select(Code, CommonName, ScientificName) %>%
  unique() %>%
  mutate(dup_code = duplicated(Code),
         dup_common = duplicated(CommonName),
         dup_sci = duplicated(ScientificName)) %>%
  filter(dup_code == T | dup_common == T | dup_sci == T)

# filter for duplicates
dup_codes2 <- surveys2 %>%
  select(Code, CommonName, ScientificName) %>%
  unique() %>%
  filter(Code %in% dup_codes$Code | CommonName %in% dup_codes$CommonName | ScientificName %in% dup_codes$ScientificName)

# look at names
dup_codes2 %>%
  arrange(Code, CommonName, ScientificName) %>%
  data.frame()

dup_codes2 %>%
  arrange(ScientificName) %>%
  data.frame()

dup_codes2 %>%
  arrange(CommonName) %>%
  data.frame()

# find alternative names
plant_list %>%
  filter(str_detect(CommonName, "aple")) # maple trees

# revise duplicates
# add a column for if something is a plant
plant_list2 <- plant_list %>%
  mutate(Code = case_when(ScientificName == "Cephalanthus occidentalis" ~ "BUBU",
                          ScientificName == "Cyperus articulatus" ~ "UMSE",
                          ScientificName == "Fuirena scirpoidea" ~ "RUFU",
                          ScientificName == "Fuirena squarrosa" ~ "UMGR",
                          ScientificName == "Lemna minor" ~ "CODU",
                          ScientificName == "Myriophyllum sp." ~ "WAMI",
                          ScientificName == "Nasturtium floridanum" ~ "FLWC",
                          ScientificName == "Salvinia molesta" ~ "GISA",
                          ScientificName == "Salvinia sp." ~ "WAFE",
                          ScientificName == "Schoenoplectus californicus (syn. Scirpus californicus)" ~ "GIBU",
                          ScientificName == "Sparganium americanum" ~ "AMBU",
                          ScientificName == "Utricularia floridana" ~ "FYBL",
                          ScientificName == "Utricularia gibba" ~ "HUBL",
                          ScientificName == "Utricularia inflata" ~ "FLBT",
                          ScientificName == "Vallisneria americana" ~ "AMEE",
                          ScientificName == "Wolffiella gladiata" ~ "FLMU",
                          TRUE ~ Code),
         CommonName = case_when(Code == "BRPE" ~ "Brazilian peppertree",
                                Code == "LFPR" ~ "large-flower primrose-willow",
                                Code == "MAST" ~ "man-made structure, such as a fishing pier or boat dock",
                                Code == "MELA" ~ "Melaleuca",
                                Code == "N/A" ~ "No access",
                                Code == "NYMP" ~ "water-lily species",
                                Code == "PASP" ~ "Paspalum",
                                Code == "POND" ~ "large-leaf pondweed",
                                Code == "RUFU" ~ "Rush Fuirena",
                                Code == "SHPR" ~ "shrubby primrose willow",
                                Code == "TSBU" ~ "three-square bulrush",
                                Code == "TUSS" ~ "tussock; 1/more species growing on a buoyant organic mat",
                                Code == "TWMG" ~ "tropical white morning- glory or moonflower",
                                Code == "UMGR" ~ "umbrella-grass",
                                Code == "UMSE" ~ "umbrella sedge",
                                Code == "WAFE" ~ "water fern",
                                TRUE ~ CommonName),
         ScientificName = case_when(Code == "LFPR" ~ "Ludwigia grandiflora",
                                    Code == "PANI" ~ "Panicum sp.",
                                    Code == "PASP" ~ "Paspalum fluitans",
                                    Code == "SLSP" ~ "Eleocharis baldwinii",
                                    Code == "SPAT" ~ "Nuphar luteum",
                                    Code == "SWGR" ~ "Luziola fluitans (syn. Hydrochloa caroliniensis)",
                                    Code == "WATE" ~ "Rorippa nasturtium-aquaticum",
                                    ScientificName %in% c("NULL", "No Access", "Man-Made Structure", "Tussock") ~ NA_character_,
                                    TRUE ~ ScientificName),
         Plant = case_when(!is.na(ScientificName) | str_detect(CommonName, "asters|plant|algae|cover") == T ~ 1,
                           TRUE ~ 0)) %>%
  select(-c(AOI, Year)) %>%
  unique() %>%
  mutate(across(where(is.character), ~na_if(., "ULL")))

# check that revisions work
plant_list2 %>%
  filter(!is.na(ScientificName)) %>%
  select(Code, CommonName, ScientificName) %>%
  unique() %>%
  mutate(dup_code = duplicated(Code),
         dup_common = duplicated(CommonName),
         dup_sci = duplicated(ScientificName)) %>%
  filter(dup_code == T | dup_common == T | dup_sci == T)




# import separately (different format)
# fwri_KingsBay_PlantDataFebruary2015
# fwri_KingsBay_PlantDataOctober2015
# fwri_LittleFish_PlantData2018
# fwri_Virginia_PlantData2019


# rename "No Access", "Depth (ft)"