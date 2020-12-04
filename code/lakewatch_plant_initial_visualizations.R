#### info ####

# goal: visualize Lakewatch plant surveys


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
lw_plant <- read_csv("original-data/Lakewatch_Plant_Surveys.csv",
                     col_types = cols(.default = col_double(),
                                      County = col_character(),
                                      Lake = col_character(),
                                      Date = col_character(),
                                      Common_name = col_character(),
                                      Genus = col_character(),
                                      Species = col_character(),
                                      Genus_species = col_character(),
                                      Other_name = col_character()))
# specify column types because Other_name causes an error


#### edit data ####

# check for missing location data
sum(is.na(lw_plant$Lake))
sum(is.na(lw_plant$County))
sum(is.na(lw_plant$Genus_species))

# check missing species
lw_plant %>%
  filter(is.na(Genus_species)) %>%
  data.frame()

# check frequency
lw_plant %>%
  mutate(rows_stations = round(N_rows / Stations * 100)) %>%
  select(Stations, N_rows, Frequency_percent, rows_stations) %>%
  filter(rows_stations != Frequency_percent) %>%
  data.frame()
# many cases where frequency doesn't match expected value...asked Mark
# when this is figured out, add evenness to lake summary below and species-specific abundance to lw_plant2

# add columns
lw_plant2 <- lw_plant %>%
  mutate(date = as.Date(paste(Month, Day, Year, sep = "-"), "%m-%d-%Y"),
         lake_county = paste(Lake, County, sep = "_") %>%
           as.factor(),
         lake_group = cut(as.numeric(lake_county), breaks = 10),
         lake_county = fct_rev(lake_county),
         tot_biomass_kg_m2 = rowSums(.[c("Em_biomass_kg_m2", "Fl_biomass_kg_m2", "Sub_biomass_kg_m2")]))

# summarize by lake
lw_plant_lake <- lw_plant2 %>%
  group_by(County, Lake, Year, Month, Day, Date, Stations, Em_fl_zone_width_ft, Em_biomass_kg_m2, Fl_biomass_kg_m2, Sub_biomass_kg_m2, Lake_depth_m, Percent_area_covered, Percent_volume_inhabited, date, lake_county, lake_group, tot_biomass_kg_m2) %>%
  summarise(richness = sum(!is.na(unique(Genus_species)))) %>%
  ungroup()

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
