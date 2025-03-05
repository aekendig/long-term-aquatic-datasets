#### set-up ####

# load packages
library(tidyverse)

# import datasets
target_dat <- read_csv("intermediate-data/FWC_plant_management_target_analysis_formatted.csv")
methods_dat <- read_csv("intermediate-data/FWC_plant_management_methods_analysis_formatted.csv")
comp_dat <- read_csv("intermediate-data/fwc_double_sampling_data.csv")

# figure settings
source("code/code-for-pub/figure_settings.R")


#### timeline figure ####

# create dataset
timeline_dat <- target_dat %>%
  summarize(startYear = min(SurveyYear),
            endYear = max(SurveyYear),
            midYear = (endYear + startYear) / 2,
            waterbodies = n_distinct(AreaOfInterestID)) %>%
  mutate(dataset = "full") %>%
  full_join(methods_dat %>%
              summarize(startYear = min(SurveyYear),
                        endYear = max(SurveyYear),
                        midYear = (endYear + startYear) / 2,
                        waterbodies = n_distinct(AreaOfInterestID)) %>%
              mutate(dataset = "detailed\nmgmt.\nrecords")) %>%
  full_join(comp_dat %>%
              summarize(startYear = min(fwc_Year),
                        endYear = max(fwc_Year),
                        midYear = (endYear + startYear) / 2,
                        waterbodies = n_distinct(PermanentID)) %>%
              mutate(dataset = "methods\ncomparison")) %>%
  mutate(dataset = fct_relevel(dataset,
                               "methods\ncomparison", "detailed\nmgmt.\nrecords", "full"))

# figure
timeline_fig <- ggplot(timeline_dat, aes(y = dataset)) +
  geom_segment(aes(x = startYear, xend = endYear, color = dataset),
               linewidth = 15) +
  geom_text(aes(x = midYear, label = waterbodies), size = paper_text_size) +
  scale_color_manual(values = col_pal) +
  labs(x = "Year", y = "Dataset") +
  def_theme_paper +
  theme(legend.position = "none",
        axis.text.x = element_text(hjust = 0.8))

ggsave("output/dataset_timeline.png", timeline_fig, 
       width = 9, height = 5, units = "cm", dpi = 600)
