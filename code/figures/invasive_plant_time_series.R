#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(cowplot)

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")


#### time series by waterbody ####

# combine plant and control
inv_fwc <- inv_plant %>%
  inner_join(inv_ctrl)

# list of waterbodies
perm_ids <- inv_fwc %>%
  select(PermanentID, AreaName) %>%
  unique() %>%
  arrange(AreaName) %>%
  pull(PermanentID)

# loop through IDs and make figure
pdf("output/invasive_plant_time_series_by_waterbody.pdf")

for(i in 1:length(perm_ids)){
  
  # subset data
  subdat <- inv_fwc %>% filter(PermanentID == perm_ids[i])
  subdat_ctrl <- subdat %>% filter(Lag1Treated == 1)
  subdat_name <- subdat %>% select(AreaName) %>% pull() %>% unique()
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = PropCovered * 100, color = CommonName)) +
          geom_vline(data = subdat_ctrl, aes(xintercept = GSYear)) +
          geom_line() +
          geom_point() + 
          facet_wrap(~ CommonName, scales = "free_y", ncol = 1) +
          labs(x = "Year", y = "Percent area covered", title = subdat_name) +
          def_theme_paper +
          theme(strip.text = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 8)))

  }

dev.off()


#### edit data ####

# NEED TO UPDATE FOR NEW TIME_INT_FUN IN GENERIC-FUNCTIONS
# ADDED INPUTS OF TAXON AND DAT_IN
# see fwc_invasive_plant_analysis.R

# complete time intervals
inv_time_int <- inv_fwc %>%
  select(GSYear) %>%
  unique() %>%
  mutate(out = map(GSYear, time_int_fun)) %>%
  unnest(cols = out) %>%
  mutate(data_points = years_out * lakes) 

ggplot(inv_time_int, aes(x = data_points)) +
  geom_histogram()

# select largest number of datapoints
inv_time_int %>%
  filter(data_points == max(data_points))

# filter 
inv_fwc2 <- inv_fwc %>%
  filter(GSYear >= 1994 & GSYear < (1994 + 26)) %>%
  group_by(PermanentID) %>%
  mutate(NAVals = sum(is.na(EstAreaCovered_ha))) %>%
  ungroup() %>%
  filter(NAVals == 0)

# make sure it worked
length(unique(inv_fwc2$PermanentID))


#### visualize ####

# total area covered
inv_fwc2 %>%
  group_by(GSYear, CommonName) %>%
  summarise(AreaCovered = sum(EstAreaCovered_ha)) %>%
  ungroup() %>%
  mutate(LogArea = log10(AreaCovered),
         CommonName = tolower(CommonName)) %>%
  ggplot(aes(x = GSYear, y = LogArea, color = CommonName)) +
  geom_line(size = 1.5) +
  scale_color_viridis_d(end = 0.7, name = "Invasive species") +
  def_theme_paper +
  labs(x = "Year", y = expression(paste("Statewide area covered (", log[10], " ha)", sep = "")))

# number of lakes
lakes_fig <- inv_fwc2 %>%
  filter(SpeciesPresent == 1) %>%
  group_by(GSYear, CommonName) %>%
  summarise(Lakes = n_distinct(PermanentID)) %>%
  ungroup() %>%
  group_by(CommonName) %>%
  mutate(MeanLakes = mean(Lakes),
         TextLakes = if_else(CommonName == "Water hyacinth", 
                             MeanLakes + 12, MeanLakes - 7)) %>%
  ungroup() %>%
  mutate(CommonName = tolower(CommonName)) %>%
  ggplot(aes(x = GSYear, y = Lakes, color = CommonName)) +
  geom_hline(aes(yintercept = MeanLakes, color = CommonName), linetype = "dashed") +
  geom_line() +
  geom_text(x = 2012, aes(label = CommonName, y = TextLakes), check_overlap = T,
            size = paper_text_size) +
  scale_color_viridis_d(end = 0.7, guide = "none") +
  labs(x = "Year", y = "Waterbodies occupied (of 213)") +
  def_theme_paper +
  theme(axis.text.x = element_text(size = 8, color="black", hjust = 0.75),
        axis.title.y = element_text(size = 10, color="black", hjust = -0.4))

# area per lake
prop_fig <- inv_fwc2 %>%
  group_by(GSYear, CommonName) %>%
  summarise(PropPerLake = mean(PropCovered)) %>%
  ungroup() %>%
  group_by(CommonName) %>%
  mutate(MeanProp = mean(PropPerLake)) %>%
  ungroup() %>%
  ggplot(aes(x = GSYear, y = PropPerLake * 100, color = CommonName)) +
  geom_hline(aes(yintercept = MeanProp * 100, color = CommonName), linetype = "dashed") +
  geom_line() +
  scale_color_viridis_d(end = 0.7, guide = "none") +
  labs(x = "Year", y = "Avg. percent area covered") +
  def_theme_paper +
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 8, color="black", hjust = 0.75))

# combine
pdf("output/invasive_plant_time_series.pdf", width = 5, height = 2.5)
plot_grid(lakes_fig, prop_fig,
          nrow = 1,
          labels = LETTERS[1:2])
dev.off()