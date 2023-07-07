#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)
library(pals)
library(cowplot)
library(patchwork)

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_dat <- read_csv("intermediate-data/FWC_invasive_plant_analysis_formatted.csv")
nat_dat <- read_csv("intermediate-data/FWC_native_richness_analysis_formatted.csv")

# order of invasive species
inv_order <- c("hydrilla", "water hyacinth", "water lettuce", 
               "Cuban bulrush", "torpedograss", "para grass")

# figure aesthetics
shape_pal <- c(24, 22, 19)
col_pal <- kelly()[c(3:6, 10:11)]


#### invasive plant time series by waterbody ####

# list of waterbodies
aoi_ids <- inv_dat %>%
  select(AreaOfInterestID, AreaOfInterest) %>%
  unique() %>%
  arrange(AreaOfInterest) %>%
  pull(AreaOfInterestID)

# loop through IDs and make figure
pdf("output/invasive_plant_time_series_by_waterbody.pdf")

for(i in 1:length(aoi_ids)){
  
  # subset data
  subdat <- inv_dat %>% filter(AreaOfInterestID == aoi_ids[i])
  subdat_ctrl <- subdat %>% filter(Treated == 1)
  subdat_name <- subdat %>% select(AreaOfInterest) %>% pull() %>% unique()
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = PercCovered, color = CommonName)) +
          geom_vline(data = subdat_ctrl, aes(xintercept = GSYear), alpha = 0.5) +
          geom_line() +
          geom_point() + 
          facet_wrap(~ CommonName, scales = "free_y") +
          labs(x = "Year", y = "Percent area covered", title = subdat_name) +
          def_theme_paper +
          theme(strip.text = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 8)))

  }

dev.off()


#### invasive plant time series ####

inv_dat2 <- inv_dat %>%
  left_join(nat_dat %>%
              select(AreaOfInterestID, CommonName, GSYear) %>%
              mutate(native = 1)) %>%
  mutate(CommonName = if_else(CommonName == "Cuban bulrush", CommonName,
                             tolower(CommonName)) %>%
           fct_relevel(inv_order))

# summarize
inv_sum <- inv_dat2 %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(PercCovered),
            n = n_distinct(AreaOfInterestID)) %>%
  ungroup() %>%
  mutate(inv_group = if_else(CommonName == "Cuban bulrush" & GSYear < 2005, "Cuban bulrush 1",
                             CommonName))

inv_sum_nat <- inv_dat2 %>%
  filter(native == 1) %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(PercCovered),
            n = n_distinct(AreaOfInterestID)) %>%
  ungroup()

# are all lakes for a year included in subsets?
inv_sum_nat %>%
  select(GSYear, CommonName, n) %>%
  rename(n_nat = n) %>%
  full_join(inv_sum) %>%
  filter(n != n_nat)
# removed one lake each in 1999 and 2003 from native plant data

# figure
inv_fig <- ggplot(inv_sum, aes(x = GSYear, y = y, color = CommonName, group = inv_group)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                width = 0, size = 0.25) +
  geom_line(linewidth = 0.25) +
  geom_point(size = 1, stroke = 0.25) +
  scale_color_manual(values = col_pal, name = "Invasive") +
  labs(x = "Year", y = "Invasive plant cover (%)") +
  def_theme_paper +
  theme(axis.text.x = element_text(hjust = 0.8))


#### management time series ####

# summarize
mgmt_sum <- inv_dat2 %>%
  group_by(GSYear, CommonName) %>%
  summarize(s = sum(Treated),
            n = n_distinct(AreaOfInterestID)) %>%
  ungroup() %>%
  mutate(inv_group = if_else(CommonName == "Cuban bulrush" & GSYear < 2005, "Cuban bulrush 1",
                             CommonName),
         y = 100 * s/n)

# figure
mgmt_fig <- ggplot(mgmt_sum, aes(x = GSYear, y = y, color = CommonName, group = inv_group)) +
  geom_line(linewidth = 0.25) +
  geom_point(size = 1, stroke = 0.25) +
  scale_color_manual(values = col_pal, name = "Invasive") +
  labs(x = "Year", y = "Waterbodies managed (%)") +
  def_theme_paper +
  theme(axis.text.x = element_text(hjust = 0.8))


#### richness time series ####

# summarize
rich_sum <- nat_dat %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(Richness),
            n = n_distinct(AreaOfInterestID)) %>%
  ungroup() %>%
  mutate(CommonName = if_else(CommonName == "Cuban bulrush", CommonName,
                              tolower(CommonName)) %>%
           fct_relevel(inv_order),
         inv_group = case_when(CommonName == "Cuban bulrush" & GSYear < 2005 ~ "Cuban bulrush 1",
                               GSYear < 2000 ~ paste(CommonName, "early"),
                               TRUE ~ CommonName))

# figure
rich_fig <- inv_fig %+%
  rich_sum +
  labs(x = "Year", y = "Common native plant richness")


#### combine plots ####

combo_fig <- inv_fig + mgmt_fig + rich_fig +
  plot_layout(ncol = 1, guides = "collect")

ggsave("output/invasive_plant_time_series.pdf", combo_fig,
       width = 4, height = 6)


#### data check for map ####

# can I just use the invasive plant waterbodies?
nat_dat %>%
  distinct(AreaOfInterestID) %>%
  anti_join(inv_dat)

inv_dat %>%
  distinct(AreaOfInterestID) %>%
  anti_join(nat_dat)
# one lake is only used in invasive plant data