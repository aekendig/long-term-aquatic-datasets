#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# figure settings
source("code/settings/figure_settings.R")

# import data
comm_plant <- read_csv("intermediate-data/FWC_plant_community_formatted.csv")
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")


#### edit data ####

# summarize by permID
# select datasets with control data
comm_plant2 <- comm_plant %>%
  filter(Detected == 1 & PermanentID %in% inv_ctrl$PermanentID) %>%
  group_by(PermanentID) %>%
  summarize(Years = n_distinct(GSYear),
            Richness = n_distinct(TaxonName)) %>%
  ungroup() %>%
  full_join(comm_plant %>%
              filter(Detected == 1 & PermanentID %in% inv_ctrl$PermanentID) %>%
              select(PermanentID, GSYear, PreCtrl) %>%
              unique() %>% # remove duplication due to taxa
              group_by(PermanentID) %>%
              summarize(PropPostCtrl = sum(PreCtrl == "post ctrl data")/n()) %>%
              ungroup()) %>%
  full_join(comm_plant %>%
              filter(Detected == 1 & PermanentID %in% inv_ctrl$PermanentID) %>%
              select(PermanentID, TaxonName, Origin) %>%
              unique() %>% # remove duplication due to years
              group_by(PermanentID) %>%
              summarize(PropNative = sum(Origin == "Native")/n()) %>%
              ungroup())

# total species
tot_taxa <- n_distinct(comm_plant$TaxonName)


#### figure ####
pdf("output/fwc_plant_community_data.pdf", width = 6.5, height = 3)
ggplot(comm_plant2, aes(x = Years, y = Richness, size = PropPostCtrl, color = PropNative)) +
  geom_hline(yintercept = tot_taxa, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c(direction = -1, name = "Native taxa",) +
  scale_size(name = "Post control data") +
  labs(x = "Years of data", y = "Raw taxonomic richness") +
  def_theme_paper +
  theme(legend.position = c(0.15, 0.68),
        legend.box = "horizontal",
        legend.spacing.x = unit(0.01, "cm")) +
  guides(size = guide_legend(override.aes = list(color = "#20A387FF", alpha = 1), 
                             order = 2),
         color = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2)))
dev.off()