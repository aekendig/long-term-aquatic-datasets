#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")


#### edit data ####

# select hydrilla
# percent difference
hydr_dat <- inv_plant %>% filter(CommonName == "Hydrilla") %>%
  mutate(PercDiffCovered = PropCovered * 100 - InitPercCovered)

# difference histogram of crashes
hydr_dat %>%
  filter(PercDiffCovered < 0) %>%
  ggplot(aes(x = PercDiffCovered)) +
  geom_histogram(binwidth = 1)

# waterbodies with crashes >= 50%
hydr_dat2 <- hydr_dat %>%
  filter(PercDiffCovered <= -50) %>%
  select(PermanentID) %>%
  unique() %>%
  inner_join(hydr_dat)

# waterbodies
n_distinct(hydr_dat2$PermanentID)

perm_ids <- hydr_dat2 %>%
  select(PermanentID, AreaName) %>%
  unique() %>%
  arrange(AreaName)

get_dupes(perm_ids, PermanentID)

# visualize
pdf("output/fwc_hydrilla_crash_exploratory_times_series.pdf")
for(i in 1:nrow(perm_ids)){
  
  subdat <- hydr_dat2 %>%
    filter(PermanentID == perm_ids[i])
  
  print(ggplot(subdat, aes(x = GSYear, y = PropCovered)) +
          geom_point() +
          geom_line() +
          labs(x = "Year", y = "Hydrilla PAC", title = unique(subdat$AreaName)) +
          def_theme_paper)
  
}
dev.off()

# waterbody names
sort(unique(inv_plant$AreaName))

# waterbodies of interest
area_names <- c("Istokpoga, Lake", "Weohyakapka, Lake")