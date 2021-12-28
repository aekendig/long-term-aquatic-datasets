#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(cowplot)

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_fwc <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv",
                col_types = list(PrevPropCovered = col_double(),
                                 PrevPropCoveredAdj = col_double(),
                                 PrevAreaCoveredRaw_ha = col_double(),
                                 SurveyDays = col_double(),
                                 RatioCovered = col_double(),
                                 LogRatioCovered = col_double(),
                                 LogitPrevPropCovered = col_double(),
                                 LogRatioCovered = col_double(),
                                 LogitPrevPropCovered = col_double(),
                                 LogRatioCovered = col_double()))

# function to find longest time interval with the most lakes
time_int_fun <- function(year1){
  
  dat <- inv_fwc %>% # all possible surveys
    filter(TaxonName == "Hydrilla verticillata" & GSYear < 2020)
  
  dat2 <- dat %>%
    filter(GSYear >= year1 & is.na(EstAreaCovered_ha)) %>% # select missing years
    group_by(PermanentID) %>%
    summarise(year2 = min(GSYear)) %>% # identify first year missing data
    ungroup() %>%
    full_join(dat %>%
                select(PermanentID) %>%
                unique()) %>% # add all lakes (in case some had no missing data)
    mutate(year2 = replace_na(year2, max(dat$GSYear) + 1),   # assign last year (add one to count that year)
           years = year2 - year1) %>%
    group_by(years) %>%
    count() %>% # summarise number of lakes per timespan
    ungroup() %>%
    expand_grid(tibble(years_out = 0:(max(dat$GSYear) + 1 - year1))) %>% # expand for every years_out value
    filter(years >= years_out) %>%
    group_by(years_out) %>%
    summarise(lakes = sum(n))
  
  return(dat2)
}


#### edit data ####

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