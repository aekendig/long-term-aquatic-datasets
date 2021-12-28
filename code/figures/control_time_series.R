#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# figure settings
source("code/settings/figure_settings.R")

# import data
ctrl <- read_csv("intermediate-data/FWC_control_new_formatted.csv")


#### edit data ####

# change management type labels
ctrl2 <- ctrl %>%
  mutate(MethodHerbicide = fct_recode(MethodHerbicide,
                                      "mechanical/biocontrol/other" = "not herbicide"))


#### visualize ####

pdf("output/control_time_series.pdf", width = 2.5, height = 2.5)
ctrl2 %>%
  group_by(GSYear, MethodHerbicide) %>%
  summarize(Treatments = n()) %>%
  filter(MethodHerbicide != "unknown" & GSYear > 2010 & GSYear < 2020) %>%
  group_by(MethodHerbicide) %>%
  mutate(MeanTreat = mean(Treatments)) %>%
  ungroup() %>%
  ggplot(aes(x = GSYear, y = Treatments, fill = MethodHerbicide)) +
  geom_area() +
  geom_text(x  = 2015, aes(label = MethodHerbicide, y = MeanTreat-150), 
            check_overlap = T,
            size = paper_text_size,
            color = "white") +
  scale_fill_brewer(type = "qual", palette = "Dark2", guide = "none") +
  scale_x_continuous(breaks = c(2011, 2015, 2019)) +
  labs(x = "Year", y = "Treatments") +
  def_theme_paper
dev.off()