#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom)
library(janitor)
library(patchwork)

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")


#### edit data ####

# select hydrilla
# percent difference
# add management data
hydr_dat <- inv_plant %>% filter(CommonName == "Hydrilla") %>%
  mutate(PercCovered = PropCovered * 100,
         PercDiffCovered = PercCovered - InitPercCovered) %>% 
  inner_join(inv_ctrl)

# difference histogram of crashes
hydr_dat %>%
  filter(PercDiffCovered < 0) %>%
  ggplot(aes(x = PercDiffCovered)) +
  geom_histogram(binwidth = 1)

# waterbodies with crashes >= 30%
hydr_dat2 <- hydr_dat %>%
  filter(PercDiffCovered <= -30) %>%
  select(PermanentID) %>%
  unique() %>%
  inner_join(hydr_dat)

# waterbodies
n_distinct(hydr_dat2$PermanentID)

perm_ids <- hydr_dat2 %>%
  select(PermanentID, AreaName) %>%
  unique() %>%
  arrange(AreaName)

# visualize
pdf("output/fwc_hydrilla_crash_exploratory_times_series.pdf")
for(i in 1:nrow(perm_ids)){
  
  subdat <- hydr_dat2 %>%
    filter(PermanentID == perm_ids$PermanentID[i])
  
  print(ggplot(subdat, aes(x = GSYear, y = PercCovered)) +
          geom_point() +
          geom_line() +
          labs(x = "Year", y = "Hydrilla PAC", title = perm_ids$AreaName[i]) +
          def_theme_paper)
  
}
dev.off()

# remove years without previous data
# change prop treated to percent
# make lake name/county variable
hydr_dat3 <- hydr_dat2 %>%
  filter(!is.na(InitPercCovered)) %>%
  mutate(across(ends_with("PropTreated"), ~ .x * 100),
         AreaCounty = paste0(AreaName, " (", str_to_title(County), ")")) %>%
  rename_with(~ str_replace(.x, "Prop", "Perc"), ends_with("PropTreated"))


#### initial visualizations ####

ggplot(hydr_dat3, aes(x = Lag1PercTreated, y = PercDiffCovered,
                      color = PermanentID)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ AreaCounty, scales = "free") +
  def_theme_paper +
  theme(legend.position = "none")


#### fit models ####

# fit model to each lake
crash_mods <- hydr_dat3 %>%
  select(PermanentID, AreaCounty, GSYear, PercDiffCovered,
         Lag1PercTreated, Lag2PercTreated, Lag3PercTreated) %>%
  nest(data = c(GSYear, PercDiffCovered, Lag1PercTreated, Lag2PercTreated, Lag3PercTreated)) %>%
  mutate(fit1 = map(data, ~lm(PercDiffCovered ~ Lag1PercTreated, data = .)),
         fit2 = map(data, ~lm(PercDiffCovered ~ Lag2PercTreated, data = .)),
         fit3 = map(data, ~lm(PercDiffCovered ~ Lag3PercTreated, data = .)))

# all coefficients
crash_coef <- crash_mods %>%
  mutate(tidied1 = map(fit1, tidy),
         tidied2 = map(fit2, tidy),
         tidied3 = map(fit3, tidy)) %>%
  select(PermanentID, AreaCounty, tidied1, tidied2, tidied3) %>%
  unnest(c(tidied1, tidied2, tidied3),
         names_sep = "_") %>%
  filter(tidied1_term != "(Intercept)") %>% # remove intercepts
  pivot_longer(cols = tidied1_term:tidied3_p.value,
               names_to = c("lag", ".value"),
               names_pattern = "tidied(.)_(.*)")

# visualize coefficients
ggplot(crash_coef, aes(x = estimate, color = lag)) +
  geom_density()

# extreme values
filter(crash_coef, estimate < -5000 | estimate > 2000)
# unrealistic estimates from Rodman reservoir
# probably because it has very few treatments
# and they're very small

# non-finite values
filter(crash_coef, estimate == Inf | estimate == -Inf | is.na(estimate))
# Little Lake Sawgrass had no treatment

# remove unrealistic estimates
crash_coef2 <- crash_coef %>%
  filter(!(PermanentID %in% c("120053782", "164160423")))

# visualize coefficients
ggplot(crash_coef2, aes(x = estimate, color = lag)) +
  geom_density()

crash_coef2 %>%
  group_by(lag) %>%
  summarize(est = mean(estimate),
            est_sd = sd(estimate),
            p_sig = sum(p.value < 0.1 & estimate < 0))
# estimates become more negative with longer lag time
# fewer relationships are statistically significant

# remove lakes that couldn't be used to fit models properly
hydr_dat4 <- hydr_dat3 %>%
  filter(!(PermanentID %in% c("120053782", "164160423")))


#### visualization ####

#### start here ####
# I think it would help to have the slope labelled
# maybe in upper left
# and move P value to lower right?

# areas with high treatment percentage
high_area <- tibble(AreaCounty = c("Clarke, Lake (Palm Beach)",
                                  "Deer Lake (Polk)",
                                  "Eaton, Lake (Marion)",
                                  "Lawne, Lake (Orange)",
                                  "Martin Bayou (Bay)",
                                  "May, Lake (Polk)",
                                  "Sampson Lake (Bradford)",
                                  "Spring, Lake (Polk)",
                                  "Thomas, Lake (Polk)",
                                  "Tracy, Lake (Polk)"),
                    TreatAdj = c(50, 150, 50, 50, 50, 100, 50, 50, 100, 50))

# P values for figure
crash_coef3 <- crash_coef2 %>%
  filter(lag == 1) %>%
  left_join(hydr_dat4 %>%
              group_by(PermanentID) %>%
              summarize(Lag1PercTreated = max(Lag1PercTreated),
                     PercDiffCovered = max(PercDiffCovered)) %>%
              ungroup()) %>%
  left_join(high_area) %>%
  mutate(p.round = round_half_up(p.value, 3),
         p.text = paste0("italic('P')~`=`~", p.round),
         Lag1PercTreated = if_else(!is.na(TreatAdj), Lag1PercTreated - TreatAdj, Lag1PercTreated))

# subset for figure
area_names <- sort(crash_coef3$AreaCounty)
crash_coef3_a <- crash_coef3 %>% filter(AreaCounty %in% area_names[1:32])
crash_coef3_b <- crash_coef3 %>% filter(AreaCounty %in% area_names[33:63])
hydr_dat4_a <- hydr_dat4 %>% filter(AreaCounty %in% area_names[1:32])
hydr_dat4_b <- hydr_dat4 %>% filter(AreaCounty %in% area_names[33:63])

# time series fig
ts_fig_a <- ggplot(hydr_dat4_a, aes(x = GSYear, y = PercCovered,
                      color = AreaCounty)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.5) +
  facet_grid(rows = vars(AreaCounty)) +
  labs(x = "Year", title = "PAC") +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

ts_fig_b <- ggplot(hydr_dat4_b, aes(x = GSYear, y = PercCovered,
                                    color = AreaCounty)) +
  geom_point(size = 0.5) +
  geom_line(size = 0.5) +
  facet_grid(rows = vars(AreaCounty)) +
  labs(x = "Year", title = "PAC") +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

# model fit fig
mf_fig_a <- ggplot(hydr_dat4_a, aes(x = Lag1PercTreated, y = PercDiffCovered,
                      color = AreaCounty)) +
  geom_hline(yintercept = 0, size = 0.1, linetype = "dashed") +
  geom_smooth(method = "lm", formula = "y ~ x", size = 0.5) +
  geom_text(data = crash_coef3_a, 
            aes(label = p.text), parse = T,
            hjust = 0, vjust = 1, size = 2, color = "black") +
  facet_grid(rows = vars(AreaCounty)) +
  labs(x = "Percent area managed", title = "Annual difference in PAC") +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text.y = element_text(size = 6, color="black",
                                    angle = 0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

mf_fig_b <- ggplot(hydr_dat4_b, aes(x = Lag1PercTreated, y = PercDiffCovered,
                                    color = AreaCounty)) +
  geom_hline(yintercept = 0, size = 0.1, linetype = "dashed") +
  geom_smooth(method = "lm", formula = "y ~ x", size = 0.5) +
  geom_text(data = crash_coef3_b, 
            aes(label = p.text), parse = T,
            hjust = 0, vjust = 1, size = 2, color = "black") +
  facet_grid(rows = vars(AreaCounty)) +
  labs(x = "Percent area managed", title = "Annual difference in PAC") +
  def_theme_paper +
  theme(axis.line = element_line(color = "black", size = 0.25),
        legend.position = "none",
        strip.text.y = element_text(size = 6, color="black",
                                    angle = 0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 5, color="black"),
        panel.border = element_blank(),
        panel.spacing = unit(0, "pt"))

# combine and save
comb_fig <- ts_fig_a + mf_fig_a + ts_fig_b + mf_fig_b +
  plot_layout(nrow = 1)

ggsave("output/fwc_hydrilla_crash_time_series_prediction.png", comb_fig,
       device = "png", width = 7, height = 9, units = "in")
