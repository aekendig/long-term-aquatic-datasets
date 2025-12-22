#### set-up ####

# load packages
library(tidyverse)
library(GGally)
library(ggpmisc)
library(janitor)
library(patchwork)

# import datasets
target_dat <- read_csv("intermediate-data/FWC_plant_management_target_analysis_formatted.csv")
methods_dat <- read_csv("intermediate-data/FWC_plant_management_methods_analysis_formatted.csv")
loc_dat <- read_csv("intermediate-data/FWC_lakes_formatted.csv")

# load models
load("output/native_richness_target_model.rda")
load("output/native_richness_methods_model.rda")

# figure settings
source("code/code-for-pub/figure_settings.R")


#### target model variables ####

# variables in target mod
summary(target_mod2)

# waterbody-specific change in native richness
waterbodies <- sort(unique(target_dat$AreaOfInterestID))
target_changes <- NULL
for (i in waterbodies) {
  target_sub_dat <- filter(target_dat, AreaOfInterestID == i)
  model  <- lm(NativeRichness  ~  Time, data = target_sub_dat)
  bmodel <- broom::tidy(model) # convert to table
  bmodel$AreaOfInterestID <- i # add waterbody ID
  bmodel <- bmodel[2,] # get latitude row
  target_changes <- rbind(target_changes, bmodel) # combine
}
target_changes

# format data
target_dat2 <- target_dat %>% 
  distinct(AreaOfInterestID, HydrCov, FloatCov, HydrTrtFreq, FloatTrtFreq,
           OtherTrtFreq, HydrTrtArea, FloatTrtArea, OtherTrtArea) %>% 
  left_join(target_changes %>% 
              select(AreaOfInterestID, estimate) %>% 
              rename(NativeChange = estimate)) %>% 
  left_join(loc_dat %>% 
              select(AreaOfInterestID, Latitude))

# correlations
target_dat2 %>% 
  select(-AreaOfInterestID) %>% 
  ggpairs()

# regressions 
target_vars <- setdiff(names(target_dat2), c("Latitude", "AreaOfInterestID"))
target_lats <- NULL
for (j in target_vars) {
  model  <- lm(get(j)  ~  Latitude, data = target_dat2)
  bmodel <- broom::tidy(model) # convert to table
  bmodel$term[2] <- j # assign variable to the Latitude row
  bmodel <- bmodel[2,] # get latitude row
  target_lats <- rbind(target_lats, bmodel) # combine
}
target_lats # use to manually check plot stats

# figures
flt_pac_fig <- ggplot(target_dat2, aes(x = Latitude, y = FloatCov)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste0("Est.==~", round_half_up(after_stat(b_1), 2),
                                  "~~italic(P)==~", 
                                  round_half_up(after_stat(p.value), 3),
                                 "~~R^2==~", 
                                 round_half_up(after_stat(r.squared), 2))), 
               parse = TRUE,
               output.type = "numeric",
               label.x = "right",
               label.y = "top",
               color = "blue",
               size = 3.5,
               geom = "label_npc") +
  labs(y = "Floating plant cover",
       x = "Latitude") +
  def_theme_paper

hyd_frq_fig <- ggplot(target_dat2, aes(x = Latitude, y = HydrTrtFreq)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste0("Est.==~", round_half_up(after_stat(b_1), 2),
                                  "~~italic(P)==~", 
                                  round_half_up(after_stat(p.value), 10),
                                  "~~R^2==~", 
                                  round_half_up(after_stat(r.squared), 2))), 
               parse = TRUE,
               output.type = "numeric",
               label.x = "right",
               label.y = "top",
               color = "blue",
               size = 3.5,
               geom = "label_npc") +
  labs(y = "Hydrilla mgmt. frequency",
       x = "Latitude") +
  def_theme_paper

flt_frq_fig <- ggplot(target_dat2, aes(x = Latitude, y = FloatTrtFreq)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste0("Est.==~", round_half_up(after_stat(b_1), 2),
                                  "~~italic(P)==~", 
                                  round_half_up(after_stat(p.value), 5),
                                  "~~R^2==~", 
                                  round_half_up(after_stat(r.squared), 2))), 
               parse = TRUE,
               output.type = "numeric",
               label.x = "right",
               label.y = "top",
               color = "blue",
               size = 3.5,
               geom = "label_npc") +
  labs(y = "Floating plant mgmt. frequency",
       x = "Latitude") +
  def_theme_paper

oth_frq_fig <- ggplot(target_dat2, aes(x = Latitude, y = OtherTrtFreq)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste0("Est.==~", round_half_up(after_stat(b_1), 2),
                                  "~~italic(P)==~", 
                                  round_half_up(after_stat(p.value), 3),
                                  "~~R^2==~", 
                                  round_half_up(after_stat(r.squared), 2))), 
               parse = TRUE,
               output.type = "numeric",
               label.x = "right",
               label.y = "top",
               color = "blue",
               size = 3.5,
               geom = "label_npc") +
  labs(y = "Other plant mgmt. frequency",
       x = "Latitude") +
  def_theme_paper

hyd_ext_fig <- ggplot(target_dat2, aes(x = Latitude, y = HydrTrtArea)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste0("Est.==~", round_half_up(after_stat(b_1), 2),
                                  "~~italic(P)==~", 
                                  round_half_up(after_stat(p.value), 3),
                                  "~~R^2==~", 
                                  round_half_up(after_stat(r.squared), 2))), 
               parse = TRUE,
               output.type = "numeric",
               label.x = "right",
               label.y = "top",
               color = "blue",
               size = 3.5,
               geom = "label_npc") +
  labs(y = "Hydrilla mgmt. extent",
       x = "Latitude") +
  def_theme_paper

nat_targ_fig <- ggplot(target_dat2, aes(x = Latitude, y = NativeChange)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste0("Est.==~", round_half_up(after_stat(b_1), 2),
                                  "~~italic(P)==~", 
                                  round_half_up(after_stat(p.value), 3),
                                  "~~R^2==~", 
                                  round_half_up(after_stat(r.squared), 2))), 
               parse = TRUE,
               output.type = "numeric",
               label.x = "right",
               label.y = "top",
               color = "blue",
               size = 3.5,
               geom = "label_npc") +
  labs(y = "Change in native richness",
       x = "Latitude") +
  def_theme_paper

# combine figures
target_lat_fig <- hyd_frq_fig + hyd_ext_fig + flt_pac_fig + flt_frq_fig + 
  oth_frq_fig + nat_targ_fig +
  plot_layout(nrow = 3) &
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

ggsave("output/target_latitude_figure.png", target_lat_fig,
       width = 18, height = 20, units = "cm")


#### methods model variables ####

# variables in methods mod
summary(methods_mod2)

# waterbody-specific change in native richness
waterbodies2 <- sort(unique(methods_dat$AreaOfInterestID))
methods_changes <- NULL
for (i in waterbodies2) {
  methods_sub_dat <- filter(methods_dat, AreaOfInterestID == i)
  model  <- lm(NativeRichness  ~  Time, data = methods_sub_dat)
  bmodel <- broom::tidy(model) # convert to table
  bmodel$AreaOfInterestID <- i # add waterbody ID
  bmodel <- bmodel[2,] # get latitude row
  methods_changes <- rbind(methods_changes, bmodel) # combine
}
methods_changes

# format data
methods_dat2 <- methods_dat %>% 
  distinct(AreaOfInterestID, HydrPAC, FloatPAC, TrtAreaCon, TrtAreaSys,
           TrtAreaNon, TrtFreqCon, TrtFreqSys, TrtFreqNon) %>% 
  left_join(methods_changes %>% 
              select(AreaOfInterestID, estimate) %>% 
              rename(NativeChange = estimate)) %>% 
  left_join(loc_dat %>% 
              select(AreaOfInterestID, Latitude))

# correlations
methods_dat2 %>% 
  select(-AreaOfInterestID) %>% 
  ggpairs()

# regressions 
methods_vars <- setdiff(names(methods_dat2), c("Latitude", "AreaOfInterestID"))
methods_lats <- NULL
for (j in methods_vars) {
  model  <- lm(get(j)  ~  Latitude, data = methods_dat2)
  bmodel <- broom::tidy(model) # convert to table
  bmodel$term[2] <- j # assign variable to the Latitude row
  bmodel <- bmodel[2,] # get latitude row
  methods_lats <- rbind(methods_lats, bmodel) # combine
}
methods_lats # use to manually check plot stats

# figures
flt_pac_fig2 <- ggplot(methods_dat2, aes(x = Latitude, y = FloatPAC)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste0("Est.==~", round_half_up(after_stat(b_1), 2),
                                  "~~italic(P)==~", 
                                  round_half_up(after_stat(p.value), 3),
                                  "~~R^2==~", 
                                  round_half_up(after_stat(r.squared), 2))), 
               parse = TRUE,
               output.type = "numeric",
               label.x = "right",
               label.y = "top",
               color = "blue",
               size = 3.5,
               geom = "label_npc") +
  labs(y = "Floating plant cover",
       x = "Latitude") +
  def_theme_paper

cont_ext_fig <- ggplot(methods_dat2, aes(x = Latitude, y = TrtAreaCon)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste0("Est.==~", round_half_up(after_stat(b_1), 2),
                                  "~~italic(P)==~", 
                                  round_half_up(after_stat(p.value), 3),
                                  "~~R^2==~", 
                                  round_half_up(after_stat(r.squared), 2))), 
               parse = TRUE,
               output.type = "numeric",
               label.x = "right",
               label.y = "top",
               color = "blue",
               size = 3.5,
               geom = "label_npc") +
  labs(y = "Contact herbicide extent",
       x = "Latitude") +
  def_theme_paper

sys_ext_fig <- ggplot(methods_dat2, aes(x = Latitude, y = TrtAreaSys)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste0("Est.==~", round_half_up(after_stat(b_1), 2),
                                  "~~italic(P)==~", 
                                  round_half_up(after_stat(p.value), 3),
                                  "~~R^2==~", 
                                  round_half_up(after_stat(r.squared), 2))), 
               parse = TRUE,
               output.type = "numeric",
               label.x = "right",
               label.y = "top",
               color = "blue",
               size = 3.5,
               geom = "label_npc") +
  labs(y = "Systemic herbicide extent",
       x = "Latitude") +
  def_theme_paper

nat_meth_fig <- ggplot(methods_dat2, aes(x = Latitude, y = NativeChange)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = F) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste0("Est.==~", round_half_up(after_stat(b_1), 2),
                                  "~~italic(P)==~", 
                                  round_half_up(after_stat(p.value), 3),
                                  "~~R^2==~", 
                                  round_half_up(after_stat(r.squared), 2))), 
               parse = TRUE,
               output.type = "numeric",
               label.x = "right",
               label.y = "top",
               color = "blue",
               size = 3.5,
               geom = "label_npc") +
  labs(y = "Change in native richness",
       x = "Latitude") +
  def_theme_paper

# combine figures
methods_lat_fig <- flt_pac_fig2 + cont_ext_fig + sys_ext_fig + nat_meth_fig +
  plot_layout(nrow = 2) &
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

ggsave("output/methods_latitude_figure.png", methods_lat_fig,
       width = 18, height = 13.3, units = "cm")
