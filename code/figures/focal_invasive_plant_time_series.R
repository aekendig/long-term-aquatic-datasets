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
pho_dat <- read_csv("intermediate-data/FWC_phosphorus_analysis_formatted.csv")
nit_dat <- read_csv("intermediate-data/FWC_nitrogen_analysis_formatted.csv")
chl_dat <- read_csv("intermediate-data/FWC_chlorophyll_analysis_formatted.csv")
sec_dat <- read_csv("intermediate-data/FWC_secchi_analysis_formatted.csv", guess_max = 7000)

# order of invasive species
inv_order <- c("Hydrilla", "Water hyacinth", "Water lettuce")

# Quarter names
quart_name <- tibble(Quarter = 1:4,
                     QuarterF = c("Apr-Jun", "Jul-Sep", "Oct-Dec", "Jan-Mar") %>%
                       fct_relevel("Apr-Jun", "Jul-Sep", "Oct-Dec"))

# figure aesthetics
shape_pal <- c(24, 22, 21)
col_pal <- kelly()[3:5]


#### invasive plant time series by waterbody ####

# list of waterbodies
perm_ids <- inv_dat %>%
  select(PermanentID, AreaName) %>%
  unique() %>%
  arrange(AreaName) %>%
  pull(PermanentID)

# loop through IDs and make figure
pdf("output/invasive_plant_time_series_by_waterbody.pdf")

for(i in 1:length(perm_ids)){
  
  # subset data
  subdat <- inv_dat %>% filter(PermanentID == perm_ids[i])
  subdat_ctrl <- subdat %>% filter(Lag1Treated == 1)
  subdat_name <- subdat %>% select(AreaName) %>% pull() %>% unique()
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = PropCovered * 100, color = CommonName)) +
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
              select(PermanentID, CommonName, GSYear) %>%
              mutate(native = 1)) %>%
  left_join(chl_dat %>%
              select(PermanentID, CommonName, GSYear) %>%
              mutate(chlorophyll = 1)) %>%
  left_join(nit_dat %>%
              select(PermanentID, CommonName, GSYear) %>%
              mutate(nitrogen = 1)) %>%
  left_join(pho_dat %>%
              select(PermanentID, CommonName, GSYear) %>%
              mutate(phosphorus = 1)) %>%
  left_join(sec_dat %>%
              select(PermanentID, CommonName, GSYear) %>%
              mutate(secchi = 1)) %>%
  filter(CommonName %in% inv_order) %>%
  mutate(CommonName = fct_relevel(CommonName, inv_order))

# summarize
inv_sum <- inv_dat2 %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(PropCovered * 100),
            n = n_distinct(PermanentID)) %>%
  ungroup()

inv_sum_nat <- inv_dat2 %>%
  filter(native == 1) %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(PropCovered * 100),
            n = n_distinct(PermanentID)) %>%
  ungroup()

inv_sum_qual <- inv_dat2 %>%
  filter(chlorophyll == 1 | nitrogen == 1 | phosphorus == 1 | secchi == 1) %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(PropCovered * 100),
            n = n_distinct(PermanentID)) %>%
  ungroup()

# are all lakes for a year included in subsets?
inv_sum_nat %>%
  select(GSYear, CommonName, n) %>%
  anti_join(inv_sum) # yes

inv_sum_qual %>%
  select(GSYear, CommonName, n) %>%
  anti_join(inv_sum) # no

# indicate native data in inv_sum
# add mean from quality data
inv_sum2 <- inv_sum %>%
  left_join(inv_sum_nat %>%
              select(GSYear, CommonName) %>%
              mutate(Native = "invasve/native")) %>%
  mutate(Native = replace_na(Native, "invasive")) %>%
  left_join(inv_sum_qual %>%
              select(GSYear, CommonName, y) %>%
              rename(yQual = "y")) %>%
  mutate(Analyses = if_else(!is.na(yQual), paste0(Native, "/\nquality"),
                        Native))

# figure
inv_fig_col <- ggplot(inv_sum2, aes(x = GSYear, y = y, color = CommonName)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                width = 0, size = 0.25) +
  geom_line(size = 0.25) +
  geom_point(aes(shape = Analyses), 
             size = 1, stroke = 0.25, fill = "white") +
  geom_point(aes(y = yQual), size = 1, stroke = 0.25) + 
  scale_shape_manual(values = shape_pal, guide = "none") +
  scale_color_manual(values = col_pal, name = "Invasive") +
  labs(x = "Year", y = "Invasive PAC") +
  def_theme_paper +
  theme(axis.text.x = element_text(hjust = 0.8))

leg_col <- get_legend(inv_fig_col)

inv_fig_shape <- inv_fig_col +
  scale_shape_manual(values = shape_pal) +
  scale_color_manual(values = col_pal, guide = "none") +
  guides(shape = guide_legend(override.aes = list(fill = c("white", "white", "black"))))

leg_shape <- get_legend(inv_fig_shape)

inv_fig <- inv_fig_shape +
  theme(legend.position = "none")


#### management time series ####

# summarize
mgmt_sum <- inv_dat2 %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(Lag1Treated * 100),
            n = n_distinct(PermanentID)) %>%
  ungroup()

mgmt_sum_nat <- inv_dat2 %>%
  filter(native == 1) %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(Lag1Treated * 100),
            n = n_distinct(PermanentID)) %>%
  ungroup()

mgmt_sum_qual <- inv_dat2 %>%
  filter(chlorophyll == 1 | nitrogen == 1 | phosphorus == 1 | secchi == 1) %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(Lag1Treated * 100),
            n = n_distinct(PermanentID)) %>%
  ungroup()

# indicate native data in mgmt_sum
# add mean from quality data
mgmt_sum2 <- mgmt_sum %>%
  left_join(mgmt_sum_nat %>%
              select(GSYear, CommonName) %>%
              mutate(Native = "invasve/native")) %>%
  mutate(Native = replace_na(Native, "invasive")) %>%
  left_join(mgmt_sum_qual %>%
              select(GSYear, CommonName, y) %>%
              rename(yQual = "y")) %>%
  mutate(Analyses = if_else(!is.na(yQual), paste0(Native, "/\nquality"),
                            Native))

# figure
mgmt_fig <- inv_fig %+%
  mgmt_sum2 +
  labs(x = "Year", y = "Waterbodies managed (%)")


#### richness time series ####

# summarize
rich_sum <- nat_dat %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(Richness),
            n = n_distinct(PermanentID)) %>%
  ungroup() %>%
  filter(CommonName %in% inv_order) %>%
  mutate(CommonName = fct_relevel(CommonName, inv_order))

# figure
rich_fig <- ggplot(rich_sum, aes(x = GSYear, y = y, color = CommonName)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                width = 0, size = 0.25) +
  geom_line(size = 0.25) +
  geom_point(shape = shape_pal[2], size = 1, stroke = 0.25, fill = "white") +
  scale_color_manual(values = col_pal, guide = "none") +
  labs(x = "Year", y = "Native richness") +
  def_theme_paper +
  theme(axis.text.x = element_text(hjust = 0.8))


#### chlorophyll time series ####

# summarize
chl_sum <- chl_dat %>%
  group_by(GSYear, CommonName, Quarter) %>%
  summarize(mean_cl_boot(QualityValue),
            n = n_distinct(PermanentID)) %>%
  ungroup() %>%
  filter(CommonName %in% inv_order) %>%
  mutate(CommonName = fct_relevel(CommonName, inv_order[-6])) %>%
  left_join(quart_name) %>%
  filter(QuarterF == "Jul-Sep")

# figure
chl_fig <- ggplot(chl_sum, aes(x = GSYear, y = y, color = CommonName)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                width = 0, size = 0.25) +
  geom_line(size = 0.25) +
  geom_point(size = 1, stroke = 0.25) +
  scale_color_manual(values = col_pal, guide = "none") +
  labs(x = "Year",
       y = expression(paste("Chlorophyll ", italic(a), " (", mu, "g/L)", sep = ""))) +
  def_theme_paper +
  theme(axis.text.x = element_text(hjust = 0.8))


#### combine plots ####

combo_fig <- inv_fig + mgmt_fig + leg_shape + rich_fig + chl_fig + leg_col +
  plot_layout(ncol = 3, widths = c(1, 1, 0.5))

ggsave("output/focal_invasive_plant_time_series.png", combo_fig,
       width = 6, height = 4)
