#### set-up ####

# load packages
library(tidyverse)
library(plotly)
library(janitor)

# import data
dat <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")


#### examine waterbody area relationships ####

# isolate waterbodies
wbs <- dat %>%
  distinct(PermanentID, AreaName, WaterbodyList_ha, WaterbodySum_ha, Area_ha)

# figure
fig <- ggplot(wbs, aes(x = WaterbodySum_ha, y = Area_ha)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_text(aes(label = AreaName), size = 1.5)

ggplotly(fig)

# duplicated waterbodies
get_dupes(wbs, AreaName)
# some may be from different counties
# some have multiple WaterbodySum_ha values
# some have multiple Area_ha values (how were invasive plant values broken up?)
# maybe all FWC analyses should use AOI unless these issues were dealt with