#### info ####

# goal: timeline figure


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# figure settings
def_theme <- theme_bw() +
  theme(axis.text = element_text(size = 8, color="black"),
        axis.title = element_text(size = 10, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.text = element_text(size = 10, color="black"),
        strip.background = element_blank())


#### create data ####

dat <- tibble(dataset = c("FWC plant surveys",
                          "FWC invasive plant management",
                          "FWC fish monitoring",
                          "FWC trophy catch program",
                          "FWC fish kill report",
                          "Lakewatch plant surveys",
                          "Lakewatch water quality") %>%
                rep(each = 2),
              year = c(1982, 2019,
                       1998, 2019,
                       2006, 2019,
                       2012, 2018,
                       1972, 2019,
                       1991, 2012,
                       1986, 2019)) %>%
  mutate(dataset = fct_relevel(dataset,
                               "Lakewatch water quality",
                               "Lakewatch plant surveys",
                               "FWC fish kill report",
                               "FWC trophy catch program",
                               "FWC fish monitoring",
                               "FWC invasive plant management", 
                               "FWC plant surveys")) %>%
  group_by(dataset) %>%
  mutate(years = paste(max(year) - min(year), " yr", sep = ""),
         minYear = ifelse(year == min(year), min(year) - 2, NA),
         maxYear = ifelse(year == max(year), max(year) + 2, NA),
         midYear = (min(year) + max(year)) / 2) %>%
  ungroup()


#### figure ####
pdf("output/dataset_timeline_figure.pdf", width = 6, height = 2)
ggplot(dat, aes(x = year, y = dataset, color = dataset)) +
  geom_line(size = 3) +
  geom_text(aes(x = minYear, label = year), color = "black", size = 2.5) +
  geom_text(aes(x = maxYear, label = year), color = "black", size = 2.5) +
  geom_label(aes(x = midYear, label = years), color = "black", fill = "white", size = 2.5, label.padding = unit(0.15, "lines"), label.size = 0) +
  xlab("Years of data") +
  scale_color_viridis_d() +
  def_theme +
  theme(legend.position = "none",
        axis.title.y = element_blank())
dev.off()