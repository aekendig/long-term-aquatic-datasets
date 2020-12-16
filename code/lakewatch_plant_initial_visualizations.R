#### info ####

# goal: visualize Lakewatch plant surveys


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
lw_plant <- read_csv("original-data/Lakewatch_Plant_Surveys.csv")


#### edit data ####

# check for missing location data
sum(is.na(lw_plant$Lake))
sum(is.na(lw_plant$County))
sum(is.na(lw_plant$Genus_species))

# check missing species
lw_plant %>%
  filter(is.na(Genus_species)) %>%
  data.frame()
# some have other data about the plants

# check frequency
lw_plant %>%
  mutate(rows_stations = N_rows / Stations * 100,
         rows_stations = case_when(rows_stations - round(rows_stations) == 0.5 ~ rows_stations + 0.5,
                                   TRUE ~ round(rows_stations))) %>%
  filter(rows_stations < round(Frequency_percent)) %>%
  data.frame()
# many cases where frequency doesn't match expected value
# when Frequency_percent is less than rows/stations * 100, other data were collected from the stations, but not species ID --> use the Frequency_percent value
# 20 cases where Frequency_percent is more than rows/stations * 100, off by a few points
# these are all from the same lake and sampling time, so there may be a mistake in the methods --> use rows_stations

# frequency data
lw_plant %>%
  filter(!is.na(Genus_species) & (is.na(N_rows) | is.na(Stations) | is.na(Frequency_percent)))
# none missing

# extract duplicate data
lw_plant %>%
  group_by(County, Lake, Year, Month, Day) %>%
  mutate(dup_species = duplicated(paste(Common_name, Genus_species, sep = "-"))) %>%
  ungroup() %>%
  filter(dup_species == T) %>%
  select(Genus_species) %>%
  unique() %>%
  data.frame()
# before I added Common_name, I got more because the species is unknown, but the common name is being used to distinguish them
# one species and it's unclear how to combine it - use larger value

# add columns
# don't need to correct duplicates if this dataset doesn't use those species
lw_plant2 <- lw_plant %>%
  mutate(date = as.Date(paste(Month, Day, Year, sep = "-"), "%m-%d-%Y"),
         lake_county = paste(Lake, County, sep = "_") %>%
           as.factor(),
         lake_group = cut(as.numeric(lake_county), breaks = 10),
         lake_county = fct_rev(lake_county),
         tot_biomass_kg_m2 = rowSums(.[c("Em_biomass_kg_m2", "Fl_biomass_kg_m2", "Sub_biomass_kg_m2")]),
         rows_stations = N_rows / Stations * 100,
         rows_stations = case_when(rows_stations - round(rows_stations) == 0.5 ~ rows_stations + 0.5,
                                   TRUE ~ round(rows_stations)),
         species_frequency = ifelse(rows_stations < round(Frequency_percent),
                                    rows_stations/100,
                                    Frequency_percent/100))

# summarize by lake
# sum duplicate unknown species in a genus
# use max cover for duplicate species
# add lake-level data in at end
lw_plant_lake <- lw_plant2 %>%
  filter(!is.na(Genus_species)) %>%
  group_by(County, Lake, Year, Month, Day, date, Common_name, Genus_species) %>%
  summarise(species_frequency = max(species_frequency)) %>%
  ungroup() %>%
  group_by(County, Lake, Year, Month, Day, date) %>%
  summarise(richness = n(),
            shannon = -1 * sum(species_frequency * log(species_frequency)),
            evenness = shannon / log(n())) %>%
  ungroup() %>%
  mutate(evenness = ifelse(evenness == Inf, NA_real_, evenness)) %>%
  full_join(lw_plant2 %>%
              select(County, Lake, Year, Month, Day, date, Stations, Em_fl_zone_width_ft, Em_biomass_kg_m2, Fl_biomass_kg_m2, Sub_biomass_kg_m2, Lake_depth_m, Percent_area_covered, Percent_volume_inhabited, lake_county, lake_group, tot_biomass_kg_m2) %>%
              unique())

# lake group
lake_grp = sort(unique(lw_plant2$lake_group))

# does the lake depth and em/fl zone width change within lakes across surveys?
lw_plant_lake %>%
  group_by(lake_county) %>%
  filter(n() > 1) %>%
  summarise(surveys = n(),
            Em_fl_zone_vals = length(unique(Em_fl_zone_width_ft)),
            Lake_depth_vals = length(unique(Lake_depth_m)))
# yes


#### Visualizations ####

# lake surveys over time
pdf("output/lakewatch_plant_survey_time_series.pdf")
for(i in 1:length(lake_grp)){
  
  lw_plant_lake_sub <- lw_plant_lake %>%
    filter(lake_group == lake_grp[i])
  
  print(ggplot(data = lw_plant_lake_sub,
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
pdf("output/lakewatch_plant_surveys_histogram.pdf")
lw_plant_lake %>%
  group_by(lake_county) %>%
  summarise(surveys = n()) %>%
  ungroup() %>%
  group_by(surveys) %>%
  summarise(lakes = length(lake_county)) %>%
  ggplot(aes(x = surveys, y = lakes)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = lakes), vjust = -0.25, size = 3) +
  ylab("Number of lakes") +
  xlab("Number of surveys") +
  theme_bw()
dev.off()

# lake depth histogram
pdf("output/lakewatch_plant_survey_depth_histogram.pdf")
lw_plant_lake %>%
  ggplot(aes(x = Lake_depth_m)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 1) +
  geom_vline(xintercept = mean(lw_plant_lake$Lake_depth_m, na.rm = T),
             linetype = "dashed",
             color = "blue") +
  xlab("Lake depth (m)") +
  ylab("Number of surveys") +
  theme_bw()
dev.off()

# emergent/floating zone histogram
pdf("output/lakewatch_plant_survey_em_fl_zone_histogram.pdf")
lw_plant_lake %>%
  mutate(Em_fl_zone_width_ft = Em_fl_zone_width_ft + 1) %>%
  ggplot(aes(x = Em_fl_zone_width_ft)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 0.1) +
  geom_vline(xintercept = mean(lw_plant_lake$Em_fl_zone_width_ft, na.rm = T),
             linetype = "dashed",
             color = "blue") +
  scale_x_log10() +
  xlab("Emergent/floating zone width (ft)") +
  ylab("Number of surveys") +
  theme_bw()
dev.off()

# biomass histogram
pdf("output/lakewatch_plant_survey_biomass_histogram.pdf")
lw_plant_lake %>%
  select(Em_biomass_kg_m2, Fl_biomass_kg_m2, Sub_biomass_kg_m2, tot_biomass_kg_m2) %>%
  pivot_longer(everything(),
               names_to = "zone",
               values_to = "biomass_kg_m2") %>%
  mutate(zone = substring(zone, 1, 3) %>%
           recode("Em_" = "emergent",
                  "Fl_" = "floating",
                  "Sub" = "submersed",
                  "tot" = "total"),
         log_biomass = log(biomass_kg_m2 + 1)) %>%
  ggplot(aes(x = log_biomass, fill = zone)) +
  geom_histogram(binwidth = 0.1) +
  facet_wrap(~ zone, scales = "free_y") +
  xlab("log(Biomass (kg/m2) + 1)") +
  ylab("Number of surveys") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

# area covered histogram
pdf("output/lakewatch_plant_survey_area_histogram.pdf")
lw_plant_lake %>%
  ggplot(aes(x = Percent_area_covered)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 5) +
  geom_vline(xintercept = mean(lw_plant_lake$Percent_area_covered, na.rm = T),
             linetype = "dashed",
             color = "blue") +
  xlab("Percent area covered") +
  ylab("Number of surveys") +
  theme_bw()
dev.off()

# volume histogram
pdf("output/lakewatch_plant_survey_volume_histogram.pdf")
lw_plant_lake %>%
  ggplot(aes(x = Percent_volume_inhabited)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 5) +
  geom_vline(xintercept = mean(lw_plant_lake$Percent_volume_inhabited, na.rm = T),
             linetype = "dashed",
             color = "blue") +
  xlab("Percent volume inhabited") +
  ylab("Number of surveys") +
  theme_bw()
dev.off()

# diversity histograms
pdf("output/lakewatch_plant_survey_diversity_histograms.pdf")
lw_plant_lake %>%
  ggplot(aes(x = richness)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 5) +
  geom_vline(xintercept = mean(lw_plant_lake$richness, na.rm = T),
             linetype = "dashed",
             color = "blue") +
  xlab("Species richness") +
  ylab("Number of surveys") +
  theme_bw()

lw_plant_lake %>%
  ggplot(aes(x = evenness)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 1) +
  geom_vline(xintercept = mean(lw_plant_lake$evenness, na.rm = T),
             linetype = "dashed",
             color = "blue") +
  xlab("Species evenness") +
  ylab("Number of surveys") +
  theme_bw()

lw_plant_lake %>%
  ggplot(aes(x = shannon)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 1) +
  geom_vline(xintercept = mean(lw_plant_lake$shannon, na.rm = T),
             linetype = "dashed",
             color = "blue") +
  xlab("Shannon diversity") +
  ylab("Number of surveys") +
  theme_bw()
dev.off()

# hydrilla density
pdf("output/lakewatch_plant_survey_hydrilla_time_series.pdf")
lw_plant2 %>%
  filter(Genus_species == "Hydrilla verticillata") %>%
  group_by(lake_county, Year) %>%
  summarise(frequency = mean(species_frequency)) %>%
  ggplot(aes(x = Year, y = frequency)) +
  geom_line(aes(color = lake_county), alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean", color = "black") +
  ylab("Hydrilla cover") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

# water lettuce
pdf("output/lakewatch_plant_survey_pistia_time_series.pdf")
lw_plant2 %>%
  filter(Genus_species == "Pistia stratiotes") %>%
  group_by(lake_county, Year) %>%
  summarise(frequency = mean(species_frequency)) %>%
  ggplot(aes(x = Year, y = frequency)) +
  geom_line(aes(color = lake_county), alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean", color = "black") +
  ylab("Water lettuce cover") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

# water hyacinth
# Eichhornia azurea (rooted water hyacinth) is also a species in the dataset, but there's only one record
pdf("output/lakewatch_plant_survey_eichhornia_time_series.pdf")
lw_plant2 %>%
  filter(Genus_species == "Eichhornia crassipes") %>%
  group_by(lake_county, Year) %>%
  summarise(frequency = mean(species_frequency)) %>%
  ggplot(aes(x = Year, y = frequency)) +
  geom_line(aes(color = lake_county), alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean", color = "black") +
  ylab("Water hyacinth cover") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

# other species
pdf("output/lakewatch_plant_survey_other_plants_time_series.pdf")
lw_plant2 %>%
  filter(!(Genus_species %in% c("Eichhornia crassipes", "Pistia stratiotes", "Hydrilla verticillata"))) %>%
  group_by(lake_county, Year) %>%
  summarise(frequency = mean(species_frequency)) %>%
  ggplot(aes(x = Year, y = frequency)) +
  geom_line(aes(color = lake_county), alpha = 0.5) +
  stat_summary(geom = "line", fun = "mean", color = "black") +
  ylab("Plant cover (excluding hydrilla, water lettuce, and water hyacinth)") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()