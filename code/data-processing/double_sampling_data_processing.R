#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)

# import data
fwri <- read_csv("intermediate-data/FWRI_plant_formatted.csv")
fwc <- read_csv("intermediate-data/FWC_plant_formatted_temporal_coverage.csv")
fwc_syn <- read_csv("intermediate-data/FWC_plant_survey_species_synonyms.csv")


#### find overlap between surveys ####

# rename columns
# format columns to match
# remove FWC outliers (incomplete surveys)
fwri2 <- fwri %>%
  rename_with(.cols = -c(PermanentID, ScientificName), 
              .fn = ~ paste0("fwri_", .x)) %>%
  rename(TaxonName = ScientificName) %>%
  mutate(fwri_Date = as.Date(fwri_Date),
         TaxonName = str_replace(TaxonName, "sp\\.", "spp\\."))
fwc2 <- fwc %>%
  filter(Outlier == 0) %>%
  rename(Year = SurveyYear,
         Date = SurveyDate,
         AOI = AreaOfInterest) %>%
  rename_with(.cols = -c(PermanentID, TaxonName), 
              .fn = ~ paste0("fwc_", .x)) %>%
  mutate(fwc_AOI = str_remove(fwc_AOI, ", Lake") %>%
           str_remove(" Lake"))

# isolate surveys 
fwri_surveys <- fwri2 %>%
  distinct(PermanentID, fwri_Year, fwri_AOI, fwri_Date)
fwc_surveys <- fwc2 %>%
  distinct(PermanentID, fwc_Year, fwc_AOI, fwc_Date)

# overlapping waterbodies
# calculate date difference
# require surveys to be within one-half year
survey_overlap <- inner_join(fwri_surveys, fwc_surveys, 
                             relationship = "many-to-many") %>%
  mutate(date_diff = abs(fwri_Date - fwc_Date)) %>%
  filter(date_diff <= 183)


#### troubleshoot spatial differences ####

# are any referring to different areas?
filter(survey_overlap, fwri_AOI != fwc_AOI) %>%
  distinct(PermanentID, fwri_AOI, fwc_AOI) %>%
  data.frame()

# North and South Conway split in FWRI
filter(fwc2, fwc_AOI == "Conway") %>% 
  distinct(fwc_AOI, fwc_WaterbodyAcres, fwc_ShapeArea*247.105)

# Boggy Creek is included as part of East Toho in FWC
filter(fwri2, fwri_AOI == "EastToho") %>% 
  distinct(fwri_X, fwri_Y) %>%
  ggplot(aes(x = fwri_X, y = fwri_Y)) +
  geom_point()
# omitted from FWRI

# Eustis 2 was a survey of East Toho saved under Eustis, but in the same year as another East Toho survey - omit

# Harris = Little Harris?
filter(fwc2, PermanentID == "120024301") %>% distinct(fwc_AOI, fwc_WaterbodyAcres)
# no
filter(fwri2, fwri_AOI == "Harris") %>% 
  distinct(fwri_X, fwri_Y) %>%
  ggplot(aes(x = fwri_X, y = fwri_Y)) +
  geom_point()
# Little Harris is included in the FWRI survey

# two Orange surveys?
filter(fwri2, str_detect(fwri_AOI, "Orange") == T) %>% distinct(fwri_Year, fwri_AOI)
# yes
filter(fwri2, str_detect(fwri_AOI, "Orange") == T) %>% 
  distinct(fwri_AOI, fwri_Year, fwri_X, fwri_Y) %>%
  ggplot(aes(x = fwri_X, y = fwri_Y)) +
  geom_point(size = 0.2) +
  facet_grid(fwri_Year ~ fwri_AOI)
# they're the same shape except in 2017 (Orange covers more than Open Water)
filter(fwri2, str_detect(fwri_AOI, "Orange") == T) %>%
  group_by(fwri_Year, fwri_AOI) %>%
  summarize(n_taxa = n_distinct(fwri_Code),
            .groups = "drop") %>%
  pivot_wider(names_from = fwri_AOI,
              values_from = n_taxa)
# likely same surveys except for 2017 - use Open water for more consistency over time

# are any of the permanent ID's affected by FWC temporal data processing included?
filter(survey_overlap, PermanentID %in% 
         c("107881197", "112029141", "112047993", "112049879", "167180956", "68792760")) %>%
  distinct(PermanentID, fwri_AOI, fwc_AOI)
# 5 out of 6 - manually checked if sub-sections should be included using FWRI points and GIS
# the shape of Butler in FWRI doesn't match GIS, I have a note about this in FWRI GIS methods,
# but I'm not sure if it was resolved or if I just left it in
# Silver Glen Spring Run excluded from George in FWRI survey
# Josephine Creek excluded from Josephine in FWRI survey
# Not sure what "Red Water Lake" with the same Permanent ID as Little Red Water Lake even means
# outfall canal excluded from Tarpon in FWRI survey

# remove locations that do not overlap
survey_overlap2 <- survey_overlap %>%
  filter(fwc_AOI != "Boggy Creek" & !(fwri_AOI %in% c("Eustis2", "Orange", "Butler")))

# indicate multiple AOI's within PermID
fwri_survey_overlap <- survey_overlap2 %>%
  distinct(PermanentID, fwri_Year, fwri_AOI, fwri_Date) %>%
  group_by(PermanentID) %>%
  mutate(fwri_AOIs = paste(sort(unique(fwri_AOI)), collapse = ", ")) %>%
  ungroup()

fwc_survey_overlap <- survey_overlap2 %>%
  distinct(PermanentID, fwc_Year, fwc_AOI, fwc_Date) %>%
  group_by(PermanentID) %>%
  mutate(fwc_AOIs = paste(sort(unique(fwc_AOI)), collapse = ", ")) %>%
  ungroup()

# update larger datasets based on overlap
fwri3 <- inner_join(fwri2, fwri_survey_overlap)
fwc3 <- inner_join(fwc2, fwc_survey_overlap)


#### troubleshoot temporal issues ####

# number of points in each survey
# number of days between surveys within one year
fwri_survey_dates <- fwri3 %>%
  distinct(PermanentID, fwri_Year, fwri_AOI, fwri_Date, fwri_X, fwri_Y) %>%
  count(PermanentID, fwri_AOI, fwri_Year, fwri_Date) %>%
  arrange(fwri_AOI, PermanentID, fwri_Year) %>%
  group_by(fwri_AOI, PermanentID, fwri_Year) %>%
  mutate(day_diff = fwri_Date - min(fwri_Date),
         tot_n = sum(n)) %>%
  ungroup()

# number of days between surveys within one year
fwc_survey_dates <- fwc3 %>%
  distinct(PermanentID, fwc_Year, fwc_AOI, fwc_Date) %>%
  arrange(fwc_AOI, PermanentID, fwc_Year) %>%
  group_by(fwc_AOI, PermanentID, fwc_Year) %>%
  mutate(day_diff = fwc_Date - min(fwc_Date)) %>%
  ungroup()

# do any later FWRI surveys add a ton of points?
filter(fwri_survey_dates, day_diff > 7 & (n/tot_n) > 0.1)
# manually checked all surveys within year and compared totals to other years
# all seem okay
# surveys that spanned December-January are all given the December year

# multiple surveys within one year for FWC?
filter(fwc_survey_dates, day_diff > 0)
# no


#### match taxa names ####

# taxa lists
fwri_tax <- distinct(fwri3, TaxonName)
fwc_tax <- fwc3 %>%
  filter(fwc_IsDetected == "Yes") %>%
  distinct(TaxonName)

# overlapping taxa
taxon_overlap <- inner_join(fwri_tax, fwc_tax)

# add synonyms to fwc taxa without matches
fwc_syn2 <- anti_join(fwc_tax, taxon_overlap) %>%
  left_join(fwc_syn %>%
              distinct(Taxon, syn_name) %>%
              rename(TaxonName = Taxon,
                     SynName = syn_name))

# are any taxa names in synonym list?
filter(fwc_syn2, TaxonName %in% fwc_syn$syn_name)
# yes, but it's also in the taxa list

# overlap with synonyms
(syn_overlap <- fwc_syn2 %>%
    filter(!is.na(SynName)) %>%
    inner_join(fwri_tax %>%
                 rename(SynName = TaxonName)))
# one

# update lists
taxon_overlap2 <- taxon_overlap %>%
  full_join(syn_overlap %>%
              select(TaxonName))
fwri4 <- fwri3  %>%
  left_join(syn_overlap %>%
              rename(NewTaxonName = TaxonName,
                     TaxonName = SynName)) %>%
  mutate(TaxonName = if_else(!is.na(NewTaxonName), NewTaxonName, TaxonName)) %>%
  select(-NewTaxonName)
fwri_tax2 <- distinct(fwri4, TaxonName)

# remaining unmatched taxa
tax_unmatched <- anti_join(fwc_tax, taxon_overlap2) %>%
  arrange(TaxonName) %>%
  rename(fwc_TaxonName = TaxonName) %>%
  expand_grid(anti_join(fwri_tax2, taxon_overlap2) %>%
                arrange(TaxonName) %>%
                rename(fwri_TaxonName = TaxonName))
view(tax_unmatched)
# in many cases, there is lumping in one survey (genus spp.) that doesn't have a consistent match in the other survey

# Are any species belonging to a group in the taxon overlap list still
# in the unmatched list? Their surveys would then be omitted from the 
# group by one of the monitoring programs.
genus_grps <- filter(taxon_overlap2, str_detect(TaxonName, "spp\\.") == T) %>%
  mutate(Genus = str_remove(TaxonName, " spp\\."))

tax_unmatched %>%
  distinct(fwc_TaxonName) %>%
  mutate(Genus = word(fwc_TaxonName, 1)) %>%
  inner_join(genus_grps)
# none in FWC

tax_unmatched %>%
  distinct(fwri_TaxonName) %>%
  mutate(Genus = word(fwri_TaxonName, 1)) %>%
  inner_join(genus_grps)
# none in FWRI

# species to combine in fwri survey
fwri_combine <- tax_unmatched %>%
  filter(fwri_TaxonName %in% c("Ludwigia octovalvis", "Ludwigia peruviana") &
           fwc_TaxonName == "Ludwigia octovalvis/peruviana") %>%
  rename(TaxonName = fwri_TaxonName)

# species to combine in fwc survey
fwc_combine <- tax_unmatched %>%
  filter(fwc_TaxonName %in% c("Utricularia floridana", "Utricularia inflata") &
           fwri_TaxonName == "Utricularia inflata or Utricularia floridana") %>%
  rename(TaxonName = fwc_TaxonName)

# check that fwc grouped taxa don't have areas
inner_join(fwc3, fwc_combine) %>% distinct(fwc_SpeciesAcres)
filter(fwc3, TaxonName %in% fwri_combine$fwc_TaxonName) %>% distinct(fwc_SpeciesAcres)
# no

# update lists
taxon_overlap3 <- taxon_overlap2 %>%
  full_join(fwri_combine %>%
              distinct(fwc_TaxonName) %>%
              rename(TaxonName = fwc_TaxonName)) %>%
  full_join(fwc_combine %>%
              distinct(fwri_TaxonName) %>%
              rename(TaxonName = fwri_TaxonName))
  
fwri5 <- fwri4  %>%
  left_join(fwri_combine) %>%
  mutate(TaxonName = if_else(!is.na(fwc_TaxonName), fwc_TaxonName, TaxonName)) %>%
  select(-fwc_TaxonName)

fwc4 <- fwc3 %>%
  left_join(fwc_combine) %>%
  mutate(TaxonName = if_else(!is.na(fwri_TaxonName), fwri_TaxonName, TaxonName)) %>%
  select(-fwri_TaxonName)

# species that might be the same (look up online)
# fwc: Cladium jamaicense; fwri: Cladium mariscus jamaicense - synonyms (Atlas of Florida Plants)
# fwc: Lemna/Spirodela spp.; fwri: Lemna spp. (no fwri Spirodela spp.) - both are considered "duckweed" (Wikipedia)
# fwc: Luziola fluitans; fwri: Luziola fluitans (syn. Hydrochloa caroliniensis) - synonyms (Atlas of Florida Plants)
# fwc: Nuphar advena; fwri: Nuphar luteum - Nuphar luteum lumps all the Nuphar species (https://plants.ifas.ufl.edu/plant-directory/nuphar-advena/)

# add synonyms (take out combine part)
fwc_syns <- tibble(TaxonName = "Cladium jamaicense", fwri_TaxonName = "Cladium mariscus jamaicense") %>%
  add_row(TaxonName = "Luziola fluitans", fwri_TaxonName = "Luziola fluitans (syn. Hydrochloa caroliniensis)")
fwri_syns <- tibble(fwc_TaxonName = "Cladium jamaicense", TaxonName = "Cladium mariscus jamaicense") %>%
  add_row(fwc_TaxonName = "Luziola fluitans", TaxonName = "Luziola fluitans (syn. Hydrochloa caroliniensis)")

# update lists
taxon_overlap4 <- taxon_overlap3 %>%
  full_join(fwc_syns %>%
              select(TaxonName))

fwri6 <- fwri5  %>%
  left_join(fwri_syns) %>%
  mutate(TaxonName = if_else(!is.na(fwc_TaxonName), fwc_TaxonName, TaxonName)) %>%
  select(-fwc_TaxonName)

# filter datasets based on overlapping taxa
fwri7 <- inner_join(fwri6, taxon_overlap4)
fwc5 <- inner_join(fwc4, taxon_overlap4)


#### summarize data ####

# get all fwri sample points (use data before subsetted for taxa) in waterbody
# sum sample points within same Permanent ID
# expand by taxon name
fwri_samples <- fwri6 %>%
  mutate(fwri_Point = paste(fwri_Y, fwri_X)) %>%
  group_by(fwri_Year, fwri_AOI, PermanentID, fwri_AOIs) %>%
  summarize(fwri_Points = n_distinct(fwri_Point),
            fwri_Date = max(fwri_Date),
            .groups = "drop") %>%
  group_by(fwri_Year, PermanentID, fwri_AOIs) %>%
  summarize(fwri_TotalPoints = sum(fwri_Points),
            fwri_Date = max(fwri_Date),
            .groups = "drop") %>%
  expand_grid(TaxonName = unique(fwri7$TaxonName))

# choose highest abundance at each sample point for species that are to be combined
# including date doesn't create more rows (each point surveyed one date)
# sum sample points within Permanent ID
# add all taxa in all waterbodies
fwri8 <- fwri7 %>%
  group_by(fwri_AOI, PermanentID, fwri_AOIs, fwri_Year, fwri_Site, fwri_Y, fwri_X, TaxonName) %>%
  summarise(fwri_Abundance = max(fwri_Abundance),
            .groups = "drop") %>%
  group_by(fwri_Year, PermanentID, fwri_AOIs, TaxonName) %>%
  summarize(fwri_IsDetected = if_else(sum(fwri_Abundance) > 0, "Yes", "No"),
            fwri_Abundance1 = sum(fwri_Abundance >= 1),
            fwri_Abundance2 = sum(fwri_Abundance >= 2),
            fwri_Abundance3 = sum(fwri_Abundance >= 3),
            .groups = "drop") %>%
  full_join(fwri_samples) %>%
  mutate(fwri_IsDetected = replace_na(fwri_IsDetected, "No"),
         fwri_Abundance1 = replace_na(fwri_Abundance1, 0),
         fwri_Abundance2 = replace_na(fwri_Abundance2, 0),
         fwri_Abundance3 = replace_na(fwri_Abundance3, 0),
         fwri_PAC1 = fwri_Abundance1/fwri_TotalPoints,
         fwri_PAC2 = fwri_Abundance2/fwri_TotalPoints,
         fwri_PAC3 = fwri_Abundance3/fwri_TotalPoints)

# presence in any AOI counts as detected
# sum acres across AOIs within permanent ID
fwc6 <- fwc5 %>%
  group_by(fwc_Year, PermanentID, fwc_ShapeArea, fwc_AOIs, TaxonName, fwc_Habitat) %>%
  summarize(fwc_IsDetected = if_else("Yes" %in% fwc_IsDetected, "Yes", "No"),
            fwc_SpeciesAcres = sum(fwc_SpeciesAcres, na.rm = T),
            fwc_WaterbodyAcres = sum(fwc_WaterbodyAcres),
            fwc_Date = max(fwc_Date),
            .groups = "drop") %>%
  mutate(fwc_PAC = fwc_SpeciesAcres/fwc_WaterbodyAcres)


#### screen for issues in detection or area ####

# permIDs
permids <- sort(unique(fwri8$PermanentID))

# richness
# initiate pdf
pdf("output/double_sampling_richness_raw.pdf")
# cycle through PermanentIDs
for(i in permids){
  
  # subset data for permanentID
  fwri_sub <- filter(fwri2, PermanentID == i & fwri_Abundance > 0) %>%
    group_by(fwri_Year, fwri_Lake) %>%
    summarize(Richness = n_distinct(fwri_CommonName),
              .groups = "drop")
  
  fwc_sub <- filter(fwc2, PermanentID == i & fwc_IsDetected == "Yes") %>%
    group_by(fwc_Year, fwc_AOI) %>%
    summarize(Richness = n_distinct(TaxonName),
              .groups = "drop")
  
  # figure of richness over time
  print(ggplot(fwri_sub, aes(x = fwri_Year, y = Richness)) +
          geom_point() +
          geom_line() +
          ggtitle(paste("FWRI", unique(fwri_sub$fwri_Lake))))
  
  print(ggplot(fwc_sub, aes(x = fwc_Year, y = Richness)) +
          geom_point() +
          geom_line() +
          ggtitle(paste("FWC", unique(fwc_sub$fwc_AOI))))
  
}
dev.off()

# export data
write_csv(fwri8, "intermediate-data/fwri_double_sampling_data.csv")
write_csv(fwc6, "intermediate-data/fwc_double_sampling_data.csv")
