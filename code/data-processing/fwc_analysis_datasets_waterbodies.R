#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
pho_dat <- read_csv("intermediate-data/FWC_phosphorus_analysis_formatted.csv")
nit_dat <- read_csv("intermediate-data/FWC_nitrogen_analysis_formatted.csv")
chl_dat <- read_csv("intermediate-data/FWC_chlorophyll_analysis_formatted.csv")
sec_dat <- read_csv("intermediate-data/FWC_secchi_analysis_formatted.csv")
nat_dat <- read_csv("intermediate-data/FWC_native_richness_analysis_formatted.csv")
inv_dat <- read_csv("intermediate-data/FWC_invasive_plant_analysis_waterbodies.csv")


#### edit data ####

# combine Permanent IDS
# add columns to indicate which quality metrics they include
qual_dat <- pho_dat %>%
  distinct(PermanentID) %>%
  mutate(Phosphorus = 1) %>%
  full_join(nit_dat %>%
              distinct(PermanentID) %>%
              mutate(Nitrogen = 1)) %>%
  full_join(chl_dat %>%
              distinct(PermanentID) %>%
              mutate(Chlorophyll = 1)) %>%
  full_join(sec_dat %>%
              distinct(PermanentID) %>%
              mutate(Secchi = 1)) %>%
  mutate(Phosphorus = replace_na(Phosphorus, 0),
         Nitrogen = replace_na(Nitrogen, 0),
         Chlorophyll = replace_na(Chlorophyll, 0),
         Secchi = replace_na(Secchi, 0))

# edit native richness data
nat_dat2 <- nat_dat %>%
  distinct(PermanentID) %>%
  mutate(NativeRichness = 1)

# are native and invasive data the same waterbodies
nat_dat2 %>%
  anti_join(inv_dat)

inv_dat %>%
  anti_join(inv_dat)
# yes


#### export ####

write_csv(qual_dat, "intermediate-data/FWC_water_quality_analyis_waterbodies.csv")
write_csv(nat_dat2, "intermediate-data/FWC_native_richness_analysis_waterbodies.csv")
