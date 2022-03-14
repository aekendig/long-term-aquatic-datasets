#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(plotly)
library(tidyverse)
library(fixest) # FE models
library(modelsummary)
library(inspectdf) # for inspect_cor

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")
nat_plant <- read_csv("intermediate-data/FWC_common_native_plants_formatted.csv")


#### richness-area relationship ####

# richness per waterbody for all years
nat_rich <- nat_plant %>%
  filter(Detected == 1) %>%
  group_by(PermanentID, AreaName, Area_ha) %>%
  summarize(Richness = n_distinct(TaxonName)) %>%
  ungroup() %>%
  mutate(LogRich = log(Richness),
         LogArea = log(Area_ha))

# are richness and area linearly related?
ggplot(nat_rich, aes(x = Area_ha, y = Richness)) +
  geom_point()

ggplot(nat_rich, aes(x = LogArea, y = LogRich)) +
  geom_point()

nat_rich %>%
  filter(Area_ha < 90000) %>%
  ggplot(aes(x = Area_ha, y = Richness)) +
  geom_point()
# no, richness saturates with area

# fit species richness-area relationship
rich_area_mod <- lm(LogRich ~ LogArea, data = nat_rich)
summary(rich_area_mod)
c <- exp(as.numeric(coef(rich_area_mod)[1]))
z <- as.numeric(coef(rich_area_mod)[2])

# simulate relationship
nat_rich_sim <- tibble(Area_ha = seq(min(nat_rich$Area_ha), 
                                     max(nat_rich$Area_ha), 
                                     length.out = 100)) %>%
  mutate(Richness = c*Area_ha^z)

ggplot(nat_rich, aes(x = Area_ha, y = Richness)) +
  geom_point() +
  geom_line(data = nat_rich_sim) +
  coord_cartesian(xlim = c(0, 20000))
# suggests that richness isn't saturated
# not a ton of data to evaluate it


#### edit native plant data ####

# summarize richness by waterbody and year
nat_plant2 <-  nat_plant %>%
  group_by(PermanentID, GSYear, Area_ha) %>%
  summarize(Richness = sum(Detected)) %>%
  ungroup() %>%
  group_by(PermanentID) %>%
  mutate(AvgRich = mean(Richness, na.rm = T)) %>% # average waterbody richness
  ungroup() %>%
  mutate(LogArea = log(Area_ha),
         LogRich = log(Richness))

# initial visualizations
plot_ly(nat_plant2, x = ~GSYear, y = ~Richness, color = ~PermanentID) %>%
  add_lines() %>% 
  layout(showlegend = FALSE)


#### edit other data ####

# check invasive plant data availability
inv_plant %>%
  inner_join(nat_plant2 %>%
               select(PermanentID, GSYear) %>%
               unique()) %>%
  ggplot(aes(x = GSYear, y = Lag6AvgPropCovered, color = PermanentID)) +
  geom_point() +
  facet_wrap(~ CommonName, scales = "free_y") +
  theme(legend.position = "none")
# Cuban bulrush can do smaller lags, but not larger

# Avg prop covered columns
inv_avg_cols <- tibble(cols = colnames(inv_plant)) %>%
  filter(str_detect(cols, "AvgProp") == T) %>%
  pull(cols)

# make long by lag
# make wide by inv. plant species
inv_plant2 <- inv_plant %>%
  select(PermanentID, GSYear, CommonName, all_of(inv_avg_cols)) %>%
  pivot_longer(cols = all_of(inv_avg_cols),
               names_to = "Lag",
               values_to = "AvgPropCovered") %>%
  filter(!is.na(AvgPropCovered)) %>% # remove missing data
  mutate(CommonName = fct_recode(CommonName,
                                 "WaterHyacinth" = "Water hyacinth",
                                 "WaterLettuce" = "Water lettuce",
                                 "ParaGrass" = "Para grass",
                                 "CubanBulrush" = "Cuban bulrush"),
         Lag = as.numeric(str_sub(Lag, 4, 4)),
         AvgPercCovered = AvgPropCovered * 100) %>%
  select(-AvgPropCovered) %>%
  pivot_wider(names_from = CommonName,
              values_from = AvgPercCovered,
              names_glue = "{CommonName}_AvgPercCovered") %>%
  mutate(Floating_AvgPercCovered = WaterHyacinth_AvgPercCovered + WaterLettuce_AvgPercCovered,
         Floating_AvgPercCovered = if_else(Floating_AvgPercCovered > 100, 100, Floating_AvgPercCovered))

# Control columns
inv_ctrl_cols <- tibble(cols = colnames(inv_ctrl)) %>%
  filter(str_detect(cols, "Treated") == T & 
           str_detect(cols, "Lag") == T & 
           str_detect(cols, "PropTreated") == F & 
           str_detect(cols, "All") == F) %>%
  pull(cols)

# make long by lag
# make wide by control target
inv_ctrl2 <- inv_ctrl %>%
  select(PermanentID, Species, GSYear, all_of(inv_ctrl_cols)) %>% 
  unique() %>% # remove duplication of floating plant treatment
  pivot_longer(cols = all_of(inv_ctrl_cols),
               names_to = "Lag",
               values_to = "TreatFreq") %>%
  filter(!is.na(TreatFreq)) %>% # remove missing data
  mutate(Species = fct_recode(Species,
                              "Floating" = "Floating Plants (Eichhornia and Pistia)",
                              "Hydrilla" = "Hydrilla verticillata",
                              "Torpedograss" = "Panicum repens",
                              "ParaGrass" = "Urochloa mutica",
                              "CubanBulrush" = "Cyperus blepharoleptos"),
         Lag = as.numeric(str_sub(Lag, 4, 4))) %>%
  pivot_wider(names_from = Species,
              values_from = TreatFreq,
              names_glue = "{Species}_TreatFreq")


#### combine data ####

# combine native, invasive, control
nat_dat <- nat_plant2 %>%
  inner_join(inv_plant2) %>%
  inner_join(inv_ctrl2)

# filter for focal invasive plants
nat_foc_dat <- nat_dat %>%
  filter(!is.na(Hydrilla_AvgPercCovered) &
           !is.na(Hydrilla_TreatFreq) &
           !is.na(Floating_AvgPercCovered) &
           !is.na(Floating_TreatFreq)) %>%
  select(-c(CubanBulrush_AvgPercCovered, CubanBulrush_TreatFreq,
            Torpedograss_AvgPercCovered, Torpedograss_TreatFreq,
            ParaGrass_AvgPercCovered, ParaGrass_TreatFreq))

nat_all_dat <- nat_dat %>%
  filter(!is.na(Hydrilla_AvgPercCovered) &
           !is.na(Hydrilla_TreatFreq) &
           !is.na(Floating_AvgPercCovered) &
           !is.na(Floating_TreatFreq) &
           !is.na(Torpedograss_AvgPercCovered) & 
           !is.na(Torpedograss_TreatFreq) &
           !is.na(ParaGrass_AvgPercCovered) & 
           !is.na(ParaGrass_TreatFreq))

nat_cb_dat <- nat_all_dat %>%
  filter(!is.na(CubanBulrush_AvgPercCovered) & 
           !is.na(CubanBulrush_TreatFreq))


#### to do ####
# subset of lakes and years where data are balanced


#### initial visualizations ####

# covariate correlations
nat_dat %>%
  select(ends_with("TreatFreq"), 
         ends_with("AvgPercCovered"),
         Lag) %>%
  group_by(Lag) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & corr >= 0.4) %>%
  data.frame()
# all floating plants, except
# Cuban bulrush with water hyacinth a few times
# para grass with water hyacinth

# invasive plant and richness
ggplot(nat_dat, aes(x = Hydrilla_AvgPercCovered, y = Richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(nat_dat, aes(x = Floating_AvgPercCovered, y = Richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(nat_dat, aes(x = ParaGrass_AvgPercCovered, y = Richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(nat_dat, aes(x = Torpedograss_AvgPercCovered, y = Richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)

ggplot(nat_dat, aes(x = CubanBulrush_AvgPercCovered, y = Richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Lag)


#### fit focal plant models ####

# subset data
nat_foc_dat1 <- filter(nat_foc_dat, Lag == 1)
nat_foc_dat2 <- filter(nat_foc_dat, Lag == 2)
nat_foc_dat3 <- filter(nat_foc_dat, Lag == 3)
nat_foc_dat4 <- filter(nat_foc_dat, Lag == 4)
nat_foc_dat5 <- filter(nat_foc_dat, Lag == 5)
nat_foc_dat6 <- filter(nat_foc_dat, Lag == 6)

# Poisson or negative binomial?
mean(nat_foc_dat1$Richness)
var(nat_foc_dat1$Richness)
# negative binomial
# fenegbin function suggested Poisson

# fit models for focal invasive plants
nat_foc_mod1 <- fepois(Richness ~ Floating_TreatFreq + Hydrilla_TreatFreq + Hydrilla_AvgPercCovered + Floating_AvgPercCovered | PermanentID + GSYear, data = nat_foc_dat1)
nat_foc_mod2 <- update(nat_foc_mod1, data = nat_foc_dat2)
nat_foc_mod3 <- update(nat_foc_mod1, data = nat_foc_dat3)
nat_foc_mod4 <- update(nat_foc_mod1, data = nat_foc_dat4)
nat_foc_mod5 <- update(nat_foc_mod1, data = nat_foc_dat5)
nat_foc_mod6 <- update(nat_foc_mod1, data = nat_foc_dat6)


#### fit all (-CB) plant models ####

# subset data
nat_all_dat1 <- filter(nat_all_dat, Lag == 1)
nat_all_dat2 <- filter(nat_all_dat, Lag == 2)
nat_all_dat3 <- filter(nat_all_dat, Lag == 3)
nat_all_dat4 <- filter(nat_all_dat, Lag == 4)
nat_all_dat5 <- filter(nat_all_dat, Lag == 5)
nat_all_dat6 <- filter(nat_all_dat, Lag == 6)

# fit models for focal invasive plants
nat_all_mod1 <- fepois(Richness ~ Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = nat_all_dat1)
nat_all_mod2 <- update(nat_all_mod1, data = nat_all_dat2)
nat_all_mod3 <- update(nat_all_mod1, data = nat_all_dat3)
nat_all_mod4 <- update(nat_all_mod1, data = nat_all_dat4)
nat_all_mod5 <- update(nat_all_mod1, data = nat_all_dat5)
nat_all_mod6 <- update(nat_all_mod1, data = nat_all_dat6)


#### fit all plant model ####

# subset data
nat_cb_dat1 <- filter(nat_cb_dat, Lag == 1)

# add Cuban bulrush
nat_cb_mod1 <- feols(Richness ~ CubanBulrush_TreatFreq + Floating_TreatFreq + Hydrilla_TreatFreq + ParaGrass_TreatFreq + Torpedograss_TreatFreq + CubanBulrush_AvgPercCovered + Floating_AvgPercCovered + Hydrilla_AvgPercCovered + ParaGrass_AvgPercCovered + Torpedograss_AvgPercCovered | PermanentID + GSYear, data = nat_cb_dat1)



#### model coefficient figures ####

# combine models
nat_foc_mods <- list(nat_foc_mod1, nat_foc_mod2, nat_foc_mod3, nat_foc_mod4, nat_foc_mod5, nat_foc_mod6)
nat_all_mods <- list(nat_all_mod1, nat_all_mod2, nat_all_mod3, nat_all_mod4, nat_all_mod5, nat_all_mod6)

# name models
names(nat_foc_mods) <- names(nat_all_mods) <- c("1", "2", "3", "4", "5", "6")

# rename coefficients
mgmt_coef_names <- c("CubanBulrush_TreatFreq" = "Cuban bulrush",
                     "Torpedograss_TreatFreq" = "Torpedograss",
                     "ParaGrass_TreatFreq" = "Para grass",
                     "Hydrilla_TreatFreq" = "Hydrilla",
                     "Floating_TreatFreq" = "Floating plant")

inv_coef_names <- c("CubanBulrush_AvgPercCovered" = "Cuban bulrush",
                    "Torpedograss_AvgPercCovered" = "Torpedograss",
                    "ParaGrass_AvgPercCovered" = "Para grass",
                    "Hydrilla_AvgPercCovered" = "Hydrilla",
                    "Floating_AvgPercCovered" = "Floating plant")

# figure
foc_mgmt_fig <- modelplot(nat_foc_mods,
                          coef_map = mgmt_coef_names,
                          background = list(geom_vline(xintercept = 0, color = "black",
                                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Management\nlag\n(years)") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Management") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

foc_inv_fig <- modelplot(nat_foc_mods,
                         coef_map = inv_coef_names,
                         background = list(geom_vline(xintercept = 0, color = "black",
                                                      size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Invasion\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

all_mgmt_fig <- modelplot(nat_all_mods,
                          coef_map = mgmt_coef_names,
                          background = list(geom_vline(xintercept = 0, color = "black",
                                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Management\nlag\n(years)") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Management") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

all_inv_fig <- modelplot(nat_all_mods,
                         coef_map = inv_coef_names,
                         background = list(geom_vline(xintercept = 0, color = "black",
                                                      size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1,
                        name = "Invasion\nyears\nincluded") +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered") +
  def_theme_paper +
  theme(legend.box.margin = margin(-10, 0, -10, -10),
        plot.margin = margin(3, 0, 3, 3)) +
  guides(color = guide_legend(reverse = TRUE))

cb_mgmt_fig <- modelplot(nat_cb_mod1,
          coef_map = mgmt_coef_names,
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed"))) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Management") +
  def_theme_paper

cb_inv_fig <- modelplot(nat_cb_mod1,
          coef_map = inv_coef_names,
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed"))) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       y = "Percent area covered") +
  def_theme_paper


#### export figures and models ####

# figures
ggsave("output/fwc_native_richness_focal_invasive_model.eps", foc_inv_fig,
       device = "eps", width = 3.5, height = 3, units = "in")
ggsave("output/fwc_native_richness_focal_control_model.eps", foc_mgmt_fig,
       device = "eps", width = 3.5, height = 3, units = "in")
ggsave("output/fwc_native_richness_all_invasive_model.eps", all_inv_fig,
       device = "eps", width = 3.5, height = 5, units = "in")
ggsave("output/fwc_native_richness_all_control_model.eps", all_mgmt_fig,
       device = "eps", width = 3.5, height = 5, units = "in")
ggsave("output/fwc_native_richness_recent_invasive_model.eps", cb_inv_fig,
       device = "eps", width = 3.5, height = 4, units = "in")
ggsave("output/fwc_native_richness_recent_control_model.eps", cb_mgmt_fig,
       device = "eps", width = 3.5, height = 4, units = "in")

# model objects
save(nat_foc_mods, file = "output/fwc_native_richness_focal_invasive_models.rda")
save(nat_all_mods, file = "output/fwc_native_richness_all_invasive_models.rda")
save(nat_cb_mod1, file = "output/fwc_native_richness_recent_invasive_model.rda")


#### model prediction figures ####

# function to create raw data
raw_dat_fun <- function(dat_in){
  
  dat_out <- dat_in %>%
    mutate(across(.cols = contains("TreatFreq"), ~ .x * 3))
  
  return(dat_out)
  
}

# function to create range of values
range_dat_fun <- function(dat_in, col_name){
  
  # set all values to zero
  zero_dat <- tibble(Floating_TreatFreq = 0,
                     Hydrilla_TreatFreq = 0,
                     ParaGrass_TreatFreq = 0,
                     Torpedograss_TreatFreq = 0,
                     Floating_AvgPercCovered = 0,
                     Hydrilla_AvgPercCovered = 0,
                     ParaGrass_AvgPercCovered = 0,
                     Torpedograss_AvgPercCovered = 0)
  
  # expand based on cover or treatment
  if(str_detect(col_name, "Covered") == T) {
    
    # min and max values
    min_value <- min(dat_in[, col_name])
    max_value <- max(dat_in[, col_name])
    
    # expand for column
    dat_out <- zero_dat %>%
      select(-{{col_name}}) %>%
      expand_grid(temp_col = seq(min_value, max_value, length.out = 10))
    
  } else {
    
    # expand for column
    dat_out <- zero_dat %>%
      select(-{{col_name}}) %>%
      expand_grid(temp_col = seq(0, 1, length.out = 10))
    
  }

  # rename column
  dat_out[, col_name] <- dat_out$temp_col
  
  # return
  return(dat_out %>% select(-temp_col))
  
}

# function to create predicted data
pred_dat_fun <- function(dat_in, mod){
  
  # column names
  col_names <- c("Floating_TreatFreq",
                 "Hydrilla_TreatFreq",
                 "ParaGrass_TreatFreq",
                 "Torpedograss_TreatFreq",
                 "Floating_AvgPercCovered",
                 "Hydrilla_AvgPercCovered",
                 "ParaGrass_AvgPercCovered",
                 "Torpedograss_AvgPercCovered")
  
  # set each variable to 0 and expand one variable at a time
  dat_out <- tibble(col_name = col_names) %>%
    mutate(out = pmap(., function(col_name) 
      range_dat_fun(dat_in = dat_in, col_name = col_name))) %>%
    unnest(cols = out) %>%
    select(-col_name) %>%
    unique() %>% # remove duplicate 0 rows
    expand_grid(dat_in %>%
                  select(PermanentID, GSYear) %>%
                  unique()) %>% # repeat for each waterbody and year
    mutate(Richness = predict(mod, newdata = ., type = "response"),
           across(.cols = contains("TreatFreq"), ~ .x * 3))
  
  return(dat_out)
  
}

# lag 1
raw_dat1 <- raw_dat_fun(nat_all_dat1)
pred_dat1 <- pred_dat_fun(nat_all_dat1, nat_all_mod1)

# make data long by cover
raw_cover_dat1 <- raw_dat1 %>%
  select(PermanentID, GSYear, contains("Covered"), Richness) %>%
  select(-c(CubanBulrush_AvgPercCovered, WaterHyacinth_AvgPercCovered, WaterLettuce_AvgPercCovered)) %>%
  pivot_longer(cols = contains("Covered"),
               names_to = "CommonName",
               values_to = "AvgPercCovered") %>%
  mutate(CommonName = str_remove(CommonName, "_AvgPercCovered") %>%
           fct_recode("Floating plants" = "Floating",
                      "Para grass" = "ParaGrass"),
         LogitPAC = logit(AvgPercCovered/100, adjust = 0.001))

pred_cover_dat1 <- pred_dat1 %>%
  select(PermanentID, GSYear, contains("Covered"), Richness) %>%
  pivot_longer(cols = contains("Covered"),
               names_to = "CommonName",
               values_to = "AvgPercCovered") %>%
  mutate(CommonName = str_remove(CommonName, "_AvgPercCovered") %>%
           fct_recode("Floating plants" = "Floating",
                      "Para grass" = "ParaGrass")) %>%
  filter(AvgPercCovered > 0) %>% # 0 values can have a range of other species treatment values
  full_join(pred_dat1 %>%
              filter(Hydrilla_AvgPercCovered == 0 & # use all zero for zero
                       ParaGrass_AvgPercCovered == 0 &
                       Torpedograss_AvgPercCovered == 0 &
                       Floating_AvgPercCovered == 0 & 
                       Hydrilla_TreatFreq == 0 &
                       ParaGrass_TreatFreq == 0 &
                       Torpedograss_TreatFreq == 0 &
                       Floating_TreatFreq == 0) %>%
              select(PermanentID, GSYear, contains("Covered"), Richness) %>%
              pivot_longer(cols = contains("Covered"),
                           names_to = "CommonName",
                           values_to = "AvgPercCovered") %>%
              mutate(CommonName = str_remove(CommonName, "_AvgPercCovered") %>%
                       fct_recode("Floating plants" = "Floating",
                                  "Para grass" = "ParaGrass"))) %>%
  mutate(LogitPAC = logit(AvgPercCovered/100, adjust = 0.001))

# make data long by treatment
raw_treat_dat1 <- raw_dat1 %>%
  select(PermanentID, GSYear, contains("Treat"), Richness) %>%
  select(-CubanBulrush_TreatFreq) %>%
  pivot_longer(cols = contains("Treat"),
               names_to = "CommonName",
               values_to = "TreatFreq") %>%
  mutate(CommonName = str_remove(CommonName, "_TreatFreq") %>%
           fct_recode("Floating plants" = "Floating",
                      "Para grass" = "ParaGrass"))

pred_treat_dat1 <- pred_dat1 %>%
  select(PermanentID, GSYear, contains("Treat"), Richness) %>%
  pivot_longer(cols = contains("Treat"),
               names_to = "CommonName",
               values_to = "TreatFreq") %>%
  mutate(CommonName = str_remove(CommonName, "_TreatFreq") %>%
           fct_recode("Floating plants" = "Floating",
                      "Para grass" = "ParaGrass")) %>%
  filter(TreatFreq > 0) %>% # 0 values can have a range of other species treatment values
  full_join(pred_dat1 %>%
              filter(Hydrilla_TreatFreq == 0 & # use all zero for zero
                       ParaGrass_TreatFreq == 0 &
                       Torpedograss_TreatFreq == 0 &
                       Floating_TreatFreq == 0 &
                       Hydrilla_AvgPercCovered == 0 &
                       ParaGrass_AvgPercCovered == 0 &
                       Torpedograss_AvgPercCovered == 0 &
                       Floating_AvgPercCovered == 0) %>%
              select(PermanentID, GSYear, contains("Treat"), Richness) %>%
              pivot_longer(cols = contains("Treat"),
                           names_to = "CommonName",
                           values_to = "TreatFreq") %>%
              mutate(CommonName = str_remove(CommonName, "_TreatFreq") %>%
                       fct_recode("Floating plants" = "Floating",
                                  "Para grass" = "ParaGrass")))

# figures
cover_fig1 <- ggplot(raw_cover_dat1, aes(x = AvgPercCovered, y = Richness)) +
  geom_point(aes(color = CommonName), alpha = 0.3) +
  stat_summary(data = pred_cover_dat1, aes(fill = CommonName),
               geom = "ribbon", fun.data = "mean_cl_boot", alpha = 0.5) +
  stat_summary(data = pred_cover_dat1, aes(color = CommonName),
               geom = "line", fun = "mean", size = 0.5) +
  facet_wrap(~ CommonName, scales = "free") +
  labs(x = "Average PAC", y = "Native taxonomic richness") +
  def_theme_paper +
  theme(legend.position = "none")

treat_fig1 <- ggplot(raw_treat_dat1, aes(x = TreatFreq, y = Richness)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0,
               aes(color = CommonName)) +
  stat_summary(geom = "point", fun = "mean", size = 2,
               aes(color = CommonName)) +
  stat_summary(data = pred_treat_dat1, aes(fill = CommonName),
               geom = "ribbon", fun.data = "mean_cl_boot", alpha = 0.5) +
  stat_summary(data = pred_treat_dat1, aes(color = CommonName),
               geom = "line", fun = "mean", size = 0.5) +
  facet_wrap(~ CommonName, scales = "free") +
  labs(x = "Years managed (out of 3)", y = "Native taxonomic richness") +
  def_theme_paper +
  theme(legend.position = "none")

# save
ggsave("output/fwc_native_richness_invasive_PAC_prediction.png", cover_fig1,
       device = "png", width = 5, height = 5, units = "in")
ggsave("output/fwc_native_richness_invasive_management_prediction.png", treat_fig1,
       device = "png", width = 5, height = 5, units = "in")

# data tables
nat_treat_sum <- pred_treat_dat1 %>%
  group_by(CommonName, TreatFreq) %>%
  summarize(Richness = mean(Richness)) %>%
  ungroup() %>%
  mutate(TreatFreq = fct_recode(as.factor(TreatFreq),
                              "None" = "0",
                              "One" = "1",
                              "Two" = "2",
                              "Three" = "3")) %>%
  pivot_wider(names_from = TreatFreq,
              values_from = Richness)

nat_cover_sum <- pred_cover_dat1 %>%
  group_by(CommonName) %>%
  summarize(AvgPercCovered = max(AvgPercCovered)) %>%
  ungroup() %>%
  full_join(tibble(CommonName = unique(pred_cover_dat1$CommonName),
                   AvgPercCovered = 0)) %>%
  left_join(pred_cover_dat1) %>%
  group_by(CommonName, AvgPercCovered) %>%
  summarize(Richness = mean(Richness)) %>%
  ungroup() %>%
  mutate(AvgPercCovered = if_else(AvgPercCovered == 0, "Min", "Max")) %>%
  pivot_wider(names_from = AvgPercCovered,
              values_from = Richness)

# save data table
write_csv(nat_cover_sum, "output/fwc_native_richness_invasive_PAC_prediction.csv")
write_csv(nat_treat_sum, "output/fwc_native_richness_invasive_management_prediction.csv")