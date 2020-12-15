#### info ####

# goal: visualize FWC plant surveys


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(lubridate)

# import data
fwc_plant <- read_csv("original-data/FWC Plant Surveys.csv")
fwc_plant_new <- read_csv("original-data/FWCSurvey_2018_2019.csv")


#### edit data ####

# check for missing location data
sum(is.na(fwc_plant$WaterbodyName))
sum(is.na(fwc_plant$County))
sum(is.na(fwc_plant$SpeciesName))

sum(is.na(fwc_plant_new$WaterbodyName))
sum(is.na(fwc_plant_new$County))
sum(is.na(fwc_plant_new$SpeciesName))

# missing species
fwc_plant_new %>%
  filter(is.na(SpeciesName)) %>%
  data.frame()

# check acreage
fwc_plant %>%
  filter(WaterbodyAcres < SpeciesAcres) %>%
  select(WaterbodyName, Surveyor, SurveyDate, WaterbodyAcres, SpeciesAcres, SpeciesName) %>%
  arrange(SpeciesName, WaterbodyName, SurveyDate) %>%
  data.frame()

fwc_plant_new %>%
  filter(WaterbodyAcres < SpeciesAcres) %>%
  select(WaterbodyName, Surveyor, SurveyDate, WaterbodyAcres, SpeciesAcres, SpeciesName) %>%
  arrange(SpeciesName, WaterbodyName, SurveyDate) %>%
  data.frame()
# 30 cases where species acres exceeds waterbody acres -- check with Candice and Alex

# combine to save
fwc_acreage_question <- fwc_plant %>%
  filter(WaterbodyAcres < SpeciesAcres) %>%
  full_join(fwc_plant_new %>%
              filter(WaterbodyAcres < SpeciesAcres))

# export dataset with issue
write_csv(fwc_acreage_question, "output/FWC_SpeciesAcres_greater_than_WaterbodyAcres.csv")

# duplicate species data
fwc_plant %>%
  group_by(County, WaterbodyName, SurveyDate) %>%
  mutate(dup_species = duplicated(SpeciesName)) %>%
  ungroup() %>%
  filter(dup_species == T) %>%
  select(SpeciesName) %>%
  unique() %>%
  data.frame()
# some are not species names, can be assessed as separate species
# "Filamentous algae", "other", "spp.", "/"
# otherwise, use max value

fwc_plant_new %>%
  group_by(County, WaterbodyName, SurveyDate) %>%
  mutate(dup_species = duplicated(SpeciesName)) %>%
  ungroup() %>%
  filter(dup_species == T) %>%
  select(SpeciesName) %>%
  unique()
# none

# unique habitats for each species?
fwc_plant %>%
  filter(!is.na(SpeciesName)) %>%
  group_by(SpeciesName) %>%
  summarise(habitats = length(unique(Habitat))) %>%
  filter(habitats > 1)
# yes

fwc_plant_new %>%
  filter(!is.na(SpeciesName)) %>%
  group_by(SpeciesName) %>%
  summarise(habitats = length(unique(Habitat))) %>%
  filter(habitats > 1)
# yes

# combine dataframes (no overlapping rows)
# add columns
# sum cover for duplicate non-species
# use max cover for duplicate species
# remove data without plant identity
fwc_plant2 <- fwc_plant %>%
  mutate(date = as.Date(SurveyDate, "%m/%d/%Y"),
         WaterbodyName = ifelse(WaterbodyName == "Watermellon Pond", "Watermelon Pond", WaterbodyName)) %>%
  full_join(fwc_plant_new %>%
              mutate(date = as.Date(SurveyDate, "%m/%d/%y"))) %>%
  filter(!is.na(SpeciesName)) %>%
  group_by(County, WaterbodyName, WaterbodyAcres, SurveyDate, date, SpeciesName, Habitat) %>%
  summarise(SpeciesAcres = case_when(str_detect(SpeciesName, "Filamentous algae|other|spp.|/") == T ~ sum(SpeciesAcres, na.rm = T),
                                     TRUE ~ max(SpeciesAcres, na.rm = T))) %>%
  ungroup() %>%
  mutate(SpeciesAcres = ifelse(SpeciesAcres == -Inf, NA_real_, SpeciesAcres),
         lake_county = paste(WaterbodyName, County, sep = "_") %>%
           as.factor(),
         lake_group = cut(as.numeric(lake_county), breaks = 10),
         lake_county = fct_rev(lake_county),
         species_frequency = SpeciesAcres / WaterbodyAcres,
         year = year(date))

# summarize by lake
# redo when the over 100% species are understood
fwc_plant_lake <- fwc_plant2 %>% 
  group_by(County, WaterbodyName, lake_county, lake_group, date) %>%
  summarise(richness = length(SpeciesName)) %>%
  ungroup() %>%
  full_join(fwc_plant2 %>% 
              filter(!is.na(species_frequency) & species_frequency <= 1) %>%
              group_by(County, WaterbodyName, lake_county, lake_group, date) %>%
              summarise(shannon = -1 * sum(species_frequency * log(species_frequency)),
                        evenness = shannon / log(length(SpeciesName))) %>%
              ungroup() %>%
              mutate(evenness = ifelse(evenness == Inf, NA_real_, evenness))) %>%
  full_join(fwc_plant2 %>%
              filter(!is.na(species_frequency)) %>%
              group_by(County, WaterbodyName, lake_county, lake_group, date, Habitat) %>%
              summarise(frequency = sum(species_frequency)) %>%
              pivot_wider(names_from = Habitat,
                          values_from = frequency,
                          names_glue = "{Habitat}_frequency"))

# lake group
lake_grp = sort(unique(fwc_plant2$lake_group))


#### visualizations ####

# lake surveys over time
pdf("output/fwc_plant_survey_time_series.pdf")
for(i in 1:length(lake_grp)){
  
  fwc_plant_lake_sub <- fwc_plant_lake %>%
    filter(lake_group == lake_grp[i])
  
  print(ggplot(data = fwc_plant_lake_sub,
               aes(x = date,
                   y = lake_county,
                   color = lake_county)) +
          geom_point() +
          geom_line(size = 0.3) +
          xlab("Date") +
          ylab("Lake (Lake_County)") +
          theme_bw() +
          theme(legend.position = "none"))
  
}
dev.off()

# number of surveys
pdf("output/fwc_plant_surveys_histogram.pdf")
fwc_plant_lake %>%
  group_by(lake_county) %>%
  summarise(surveys = n()) %>%
  ungroup() %>%
  group_by(surveys) %>%
  summarise(lakes = length(lake_county)) %>%
  ggplot(aes(x = surveys, y = lakes)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = lakes), vjust = -0.25, size = 3) +
  ylab("Number of water bodies") +
  xlab("Number of surveys") +
  theme_bw()
dev.off()

# diversity histograms
pdf("output/fwc_plant_survey_diversity_histograms.pdf")
fwc_plant_lake %>%
  ggplot(aes(x = richness)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 5) +
  geom_vline(xintercept = mean(fwc_plant_lake$richness, na.rm = T),
             linetype = "dashed",
             color = "blue") +
  xlab("Species richness") +
  ylab("Number of surveys") +
  theme_bw()

fwc_plant_lake %>%
  ggplot(aes(x = evenness)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 1) +
  geom_vline(xintercept = mean(fwc_plant_lake$evenness, na.rm = T),
             linetype = "dashed",
             color = "blue") +
  xlab("Species evenness") +
  ylab("Number of surveys") +
  theme_bw()

fwc_plant_lake %>%
  ggplot(aes(x = shannon)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 1) +
  geom_vline(xintercept = mean(fwc_plant_lake$shannon, na.rm = T),
             linetype = "dashed",
             color = "blue") +
  xlab("Shannon diversity") +
  ylab("Number of surveys") +
  theme_bw()
dev.off()

# hydrilla density
pdf("output/fwc_plant_survey_hydrilla_time_series.pdf")
fwc_plant2 %>%
  filter(SpeciesName == "Hydrilla verticillata" & species_frequency <= 1) %>%
  group_by(lake_county, year) %>%
  summarise(frequency = mean(species_frequency)) %>%
  ggplot(aes(x = year, y = frequency)) +
  geom_line(aes(color = lake_county), alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean", color = "black") +
  ylab("Hydrilla cover") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

# water lettuce
pdf("output/fwc_plant_survey_pistia_time_series.pdf")
fwc_plant2 %>%
  filter(SpeciesName == "Pistia stratiotes" & species_frequency <= 1) %>%
  group_by(lake_county, year) %>%
  summarise(frequency = mean(species_frequency)) %>%
  ggplot(aes(x = year, y = frequency)) +
  geom_line(aes(color = lake_county), alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean", color = "black") +
  ylab("Water lettuce cover") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

# water hyacinth
# Eichhornia azurea (rooted water hyacinth) is also a species in the dataset, but there's only one record
pdf("output/fwc_plant_survey_eichhornia_time_series.pdf")
fwc_plant2 %>%
  filter(SpeciesName == "Eichhornia crassipes" & species_frequency <= 1) %>%
  group_by(lake_county, year) %>%
  summarise(frequency = mean(species_frequency)) %>%
  ggplot(aes(x = year, y = frequency)) +
  geom_line(aes(color = lake_county), alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean", color = "black") +
  ylab("Water hyacinth cover") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

# other species
pdf("output/fwc_plant_survey_other_plants_time_series.pdf")
fwc_plant2 %>%
  filter(!(SpeciesName %in% c("Eichhornia crassipes", "Pistia stratiotes", "Hydrilla verticillata")) &
           !is.na(species_frequency) & species_frequency <= 1) %>%
  group_by(lake_county, year) %>%
  summarise(frequency = mean(species_frequency)) %>%
  ggplot(aes(x = year, y = frequency)) +
  geom_line(aes(color = lake_county), alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean", color = "black") +
  ylab("Plant cover (excluding hydrilla, water lettuce, and water hyacinth)") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()