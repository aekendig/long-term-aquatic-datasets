#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# figure settings
source("code/settings/figure_settings.R")

# import data
ctrl <- read_csv("intermediate-data/FWC_control_new_formatted.csv")
inv_dat <- read_csv("intermediate-data/FWC_invasive_plant_analysis_formatted.csv")


#### edit data ####

# split out floating plants
# format variables
ctrl2 <- ctrl  %>%
  filter(Species == "Floating Plants (Eichhornia and Pistia)") %>%
  mutate(Species = "Eichhornia crassipes") %>%
  full_join(ctrl %>%
              filter(Species == "Floating Plants (Eichhornia and Pistia)") %>%
              mutate(Species = "Pistia stratiotes")) %>%
  full_join(ctrl %>%
              filter(Species != "Floating Plants (Eichhornia and Pistia)")) %>%
  filter(Species %in% c("Eichhornia crassipes", "Pistia stratiotes", "Hydrilla verticillata")) %>%
  mutate(MethodHerbicide = fct_recode(MethodHerbicide,
                                      "chemical" = "herbicide",
                                      "mechanical/biocontrol/other" = "not herbicide"),
         TaxonName = Species,
         Species = fct_recode(Species, "Pontederia crassipes" = "Eichhornia crassipes"))

# select lakes/years in final analysis dataset
ctrl3 <- ctrl2 %>%
  inner_join(inv_dat %>%
             select(AreaOfInterest, AreaOfInterestID, PermanentID, SurveyDate, GSYear, TaxonName))

#### START HERE ####
# control method types by species (more detailed version of current supp figure)
# look for table I made for Candice to combine like methods
# application date by species

ctrl_sum <- ctrl2 %>%
  group_by(GSYear, MethodHerbicide) %>%
  summarize(Treatments = n()) %>%
  filter(MethodHerbicide != "unknown") %>%
  group_by(MethodHerbicide) %>%
  mutate(MeanTreat = mean(Treatments)) %>%
  ungroup()

# range of percentages for text
ctrl_sum %>%
  group_by(GSYear) %>%
  mutate(TotTreat = sum(Treatments)) %>%
  ungroup() %>%
  mutate(PercTreat = 100 * Treatments / TotTreat) %>%
  filter(MethodHerbicide == "chemical") %>%
  pull(PercTreat) %>%
  range()


#### visualize ####

pdf("output/control_time_series.pdf", width = 2.5, height = 2.5)
ggplot(ctrl_sum, aes(x = GSYear, y = Treatments, fill = MethodHerbicide)) +
  geom_area() +
  geom_text(x  = 2015, aes(label = MethodHerbicide, y = MeanTreat-150), 
            check_overlap = T,
            size = paper_text_size,
            color = "white") +
  scale_fill_brewer(type = "qual", palette = "Dark2", guide = "none") +
  scale_x_continuous(breaks = c(2011, 2015, 2019)) +
  labs(x = "Year", y = "Management activities") +
  def_theme_paper
dev.off()