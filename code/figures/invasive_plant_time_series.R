#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)
library(pals)
# library(cowplot)

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_dat <- read_csv("intermediate-data/FWC_invasive_plant_analysis_formatted.csv")
nat_dat <- read_csv("intermediate-data/FWC_native_richness_analysis_formatted.csv")
pho_dat <- read_csv("intermediate-data/FWC_phosphorus_analysis_formatted.csv")
nit_dat <- read_csv("intermediate-data/FWC_nitrogen_analysis_formatted.csv")
chl_dat <- read_csv("intermediate-data/FWC_chlorophyll_analysis_formatted.csv")
sec_dat <- read_csv("intermediate-data/FWC_secchi_analysis_formatted.csv")

# order of invasive species
inv_order <- c("Hydrilla", "Water hyacinth", "Water lettuce", 
               "Torpedograss", "Para grass", "Cuban bulrush")

# Quarter names
quart_name <- tibble(Quarter = 1:4,
                     QuarterF = c("Apr-Jun", "Jul-Sep", "Oct-Dec", "Jan-Mar") %>%
                       fct_relevel("Apr-Jun", "Jul-Sep", "Oct-Dec"))


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
  mutate(CommonName = fct_relevel(CommonName, inv_order))

# summarize
inv_sum <- inv_dat2 %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(PropCovered * 100),
            n = n_distinct(PermanentID)) %>%
  ungroup() %>%
  group_by(CommonName) %>%
  mutate(ylabel = max(ymax)) %>%
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

# figure
inv_fig <- ggplot(inv_sum, aes(x = GSYear, y = y)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, size = 0.25) +
  geom_point(size = 0.75, stroke = 0.25, color = "black") +
  geom_line(size = 0.25) +
  geom_point(data = inv_sum_nat, size = 0.75, stroke = 0.25, shape = 21, fill = "white") +
  geom_errorbar(data = inv_sum_qual, aes(ymin = ymin, ymax = ymax), 
                width = 0, size = 0.25, color = "#E69F00") +
  geom_point(data = inv_sum_qual, size = 0.5, shape = 18, color = "#E69F00") +
  geom_text(aes(label = CommonName, x = mean(GSYear), y = ylabel), 
            vjust = 1, size = tiny_text_size, check_overlap = T) +
  geom_text(aes(label = paste("N = ", n), x = max(GSYear), y = ylabel), 
            vjust = 1, hjust = 1,
            size = tiny_text_size - 0.5, check_overlap = T) +
  facet_wrap(~ CommonName, scales = "free_y", ncol = 1) +
  labs(x = "Year", y = "Invasive PAC") +
  scale_y_continuous(n.breaks = 3) +
  def_theme_tiny +
  theme(axis.line.y = element_line(color = "black", size = 0.25),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(1.75, "pt"))

ggsave("output/invasive_plant_time_series.png",
       plot = inv_fig, width = 1.75, height = 2)


#### management time series ####

# summarize
mgmt_sum <- inv_dat2 %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(Lag1Treated * 100),
            n = n_distinct(PermanentID)) %>%
  ungroup() %>%
  group_by(CommonName) %>%
  mutate(ylabel = if_else(str_detect(CommonName, "Water") == T, 
                          0, max(ymax))) %>%
  ungroup() %>%
  mutate(vjust = if_else(str_detect(CommonName, "Water") == T, 
                         -0.5, 1))

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

# figure
mgmt_fig <- ggplot(mgmt_sum, aes(x = GSYear, y = y)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, size = 0.25) +
  geom_point(size = 0.75, stroke = 0.25, color = "black") +
  geom_line(size = 0.25) +
  geom_point(data = mgmt_sum_nat, size = 0.75, stroke = 0.25, shape = 21, fill = "white") +
  geom_errorbar(data = inv_sum_qual, aes(ymin = ymin, ymax = ymax), 
                width = 0, size = 0.25, color = "#E69F00") +
  geom_point(data = mgmt_sum_qual, size = 0.5, shape = 18, color = "#E69F00") +
  geom_text(aes(label = CommonName, x = min(GSYear), y = ylabel, vjust = vjust), 
            hjust = 0, size = tiny_text_size, check_overlap = T) +
  geom_text(aes(label = paste("N = ", n), x = max(GSYear), y = 0), 
            vjust = -0.5, hjust = 1,
            size = tiny_text_size - 0.5, check_overlap = T) +
  facet_wrap(~ CommonName, scales = "free_y", ncol = 1) +
  labs(x = "Year", y = "Waterbodies managed (%)") +
  scale_y_continuous(n.breaks = 3) +
  def_theme_tiny +
  theme(axis.line.y = element_line(color = "black", size = 0.25),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(1.75, "pt"))

ggsave("output/management_time_series.png",
       plot = mgmt_fig, width = 1.75, height = 2)


#### richness time series ####

#### start here ####
# need to add size adjustments from above to figures
# too all figures below
# common name positions may need to be adjusted

# summarize
rich_sum <- nat_dat %>%
  group_by(GSYear, CommonName) %>%
  summarize(mean_cl_boot(Richness),
            n = n_distinct(PermanentID)) %>%
  ungroup() %>%
  group_by(CommonName) %>%
  mutate(ylabel = max(ymax)) %>%
  ungroup() %>%
  mutate(CommonName = fct_relevel(CommonName, inv_order))

# figure
ggplot(rich_sum, aes(x = GSYear, y = y)) +
  geom_hline(yintercept = min(rich_sum$ymin) - 0.5, size = 0.25) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
  geom_point(size = 2, color = "black") +
  geom_line() +
  geom_text(aes(label = CommonName, x = mean(GSYear), y = ylabel), 
            vjust = 1.5, size = paper_text_size, check_overlap = T) +
  geom_text(aes(label = paste(n, "waterbodies"), x = max(GSYear), y = ylabel), 
            vjust = 1.5, hjust = 1,
            size = paper_text_size - 1, check_overlap = T) +
  facet_wrap(~ CommonName, scales = "free_y", ncol = 1) +
  labs(x = "Year", y = "Native richness") +
  def_theme_paper +
  theme(axis.text = element_text(size = 6, color="black"),
        axis.line.y = element_line(color = "black", size = 0.25),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(1.5, "pt"))

ggsave("output/native_richness_time_series.png",
       plot = rich_fig, width = 1.75, height = 2)


#### chlorophyll time series ####

# summarize
chl_sum <- chl_dat %>%
  group_by(GSYear, CommonName, Quarter) %>%
  summarize(mean_cl_boot(QualityValue),
            n = n_distinct(PermanentID)) %>%
  ungroup() %>%
  group_by(CommonName) %>%
  mutate(ylabel = max(ymax)) %>%
  ungroup() %>%
  mutate(CommonName = fct_relevel(CommonName, inv_order[-5])) %>%
  left_join(quart_name)

# figure
ggplot(chl_sum, aes(x = GSYear, y = y)) +
  geom_hline(yintercept = 10, size = 0.25) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = QuarterF), width = 0) +
  geom_point(aes(color = QuarterF), size = 2) +
  geom_line(aes(color = QuarterF)) +
  geom_text(aes(label = CommonName, x = mean(GSYear), y = ylabel), 
            vjust = 1.5, size = paper_text_size, check_overlap = T) +
  geom_text(aes(label = paste(n, "waterbodies"), x = max(GSYear), y = ylabel), 
            vjust = 1.5, hjust = 1,
            size = paper_text_size - 1, check_overlap = T) +
  facet_wrap(~ CommonName, scales = "free_y", ncol = 1) +
  labs(x = "Year", 
       y = expression(paste("Chlorophyll ", italic(a), " (", mu, "g/L)", sep = ""))) +
  scale_color_manual(values = kelly()[3:6], name = "Quarter") +
  def_theme_paper +
  theme(axis.text = element_text(size = 6, color="black"),
        axis.line.y = element_line(color = "black", size = 0.25),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(1.5, "pt"),
        legend.position = c(0.23, 0.1)) +
  guides(color = guide_legend(nrow = 2))

ggsave("output/chlorophyll_time_series.png",
       plot = chl_fig, width = 1.75, height = 2)


#### nitrogen time series ####

# summarize
nit_sum <- nit_dat %>%
  group_by(GSYear, CommonName, Quarter) %>%
  summarize(mean_cl_boot(QualityValue),
            n = n_distinct(PermanentID)) %>%
  ungroup() %>%
  group_by(CommonName) %>%
  mutate(ylabel = max(ymax)) %>%
  ungroup() %>%
  mutate(CommonName = fct_relevel(CommonName, inv_order[-5])) %>%
  left_join(quart_name)

# figure
ggplot(nit_sum, aes(x = GSYear, y = y)) +
  geom_hline(yintercept = 700, size = 0.25) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = QuarterF), width = 0) +
  geom_point(aes(color = QuarterF), size = 2) +
  geom_line(aes(color = QuarterF)) +
  geom_text(aes(label = CommonName, x = mean(GSYear), y = ylabel), 
            vjust = 1.5, size = paper_text_size, check_overlap = T) +
  geom_text(aes(label = paste(n, "waterbodies"), x = max(GSYear), y = ylabel), 
            vjust = 1.5, hjust = 1,
            size = paper_text_size - 1, check_overlap = T) +
  facet_wrap(~ CommonName, scales = "free_y", ncol = 1) +
  labs(x = "Year", 
       y = expression(paste("Total nitrogen (", mu, "g/L)", sep = ""))) +
  scale_color_manual(values = kelly()[3:6], name = "Quarter") +
  def_theme_paper +
  theme(axis.text = element_text(size = 6, color="black"),
        axis.line.y = element_line(color = "black", size = 0.25),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(1.5, "pt"),
        legend.position = c(0.23, 0.1)) +
  guides(color = guide_legend(nrow = 2))

ggsave("output/nitrogen_time_series.png",
       plot = nit_fig, width = 1.75, height = 2)

#### phosphorus time series ####

# summarize
pho_sum <- pho_dat %>%
  group_by(GSYear, CommonName, Quarter) %>%
  summarize(mean_cl_boot(QualityValue),
            n = n_distinct(PermanentID)) %>%
  ungroup() %>%
  group_by(CommonName) %>%
  mutate(ylabel = max(ymax)) %>%
  ungroup() %>%
  mutate(CommonName = fct_relevel(CommonName, inv_order[-5])) %>%
  left_join(quart_name)

# figure
ggplot(pho_sum, aes(x = GSYear, y = y)) +
  geom_hline(yintercept = 30, size = 0.25) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = QuarterF), width = 0) +
  geom_point(aes(color = QuarterF), size = 2) +
  geom_line(aes(color = QuarterF)) +
  geom_text(aes(label = CommonName, x = mean(GSYear), y = ylabel), 
            vjust = 1.5, size = paper_text_size, check_overlap = T) +
  geom_text(aes(label = paste(n, "waterbodies"), x = max(GSYear), y = ylabel), 
            vjust = 1.5, hjust = 1,
            size = paper_text_size - 1, check_overlap = T) +
  facet_wrap(~ CommonName, scales = "free_y", ncol = 1) +
  labs(x = "Year", 
       y = expression(paste("Total phosphorus (", mu, "g/L)", sep = ""))) +
  scale_color_manual(values = kelly()[3:6], name = "Quarter") +
  def_theme_paper +
  theme(axis.text = element_text(size = 6, color="black"),
        axis.line.y = element_line(color = "black", size = 0.25),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(1.5, "pt"),
        legend.position = c(0.23, 0.1)) +
  guides(color = guide_legend(nrow = 2))

ggsave("output/phosphorus_time_series.png",
       plot = pho_fig, width = 1.75, height = 2)


#### Secchi time series ####

# summarize
sec_sum <- sec_dat %>%
  group_by(GSYear, CommonName, Quarter) %>%
  summarize(mean_cl_boot(QualityValue),
            n = n_distinct(PermanentID)) %>%
  ungroup() %>%
  group_by(CommonName) %>%
  mutate(ylabel = max(ymax)) %>%
  ungroup() %>%
  mutate(CommonName = fct_relevel(CommonName, inv_order[-5])) %>%
  left_join(quart_name)

# figure
ggplot(sec_sum, aes(x = GSYear, y = y)) +
  geom_hline(yintercept = 45, size = 0.25) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = QuarterF), width = 0) +
  geom_point(aes(color = QuarterF), size = 2) +
  geom_line(aes(color = QuarterF)) +
  geom_text(aes(label = CommonName, x = mean(GSYear), y = ylabel), 
            vjust = 1.5, size = paper_text_size, check_overlap = T) +
  geom_text(aes(label = paste(n, "waterbodies"), x = max(GSYear), y = ylabel), 
            vjust = 1.5, hjust = 1,
            size = paper_text_size - 1, check_overlap = T) +
  facet_wrap(~ CommonName, scales = "free_y", ncol = 1) +
  labs(x = "Year", 
       y = "Secchi disk depth (cm)") +
  scale_color_manual(values = kelly()[3:6], name = "Quarter") +
  def_theme_paper +
  theme(axis.text = element_text(size = 6, color="black"),
        axis.line.y = element_line(color = "black", size = 0.25),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(1.5, "pt"),
        legend.position = c(0.23, 0.1)) +
  guides(color = guide_legend(nrow = 2))

ggsave("output/secchi_time_series.png",
       plot = sec_fig, width = 1.75, height = 2)
