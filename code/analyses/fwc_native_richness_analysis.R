#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(plotly)
library(tidyverse)
library(fixest) # FE models
library(modelsummary)
library(inspectdf) # for inspect_cor
library(patchwork)

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_analysis_formatted.csv") # plant and control data, continuous data
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
  summarize(Richness = sum(Detected),
            PrevRichness = sum(PrevDetected)) %>%
  ungroup() %>%
  mutate(RichnessDiff = Richness - PrevRichness) %>%
  filter(!is.na(RichnessDiff))

# initial visualizations
plot_ly(nat_plant2, x = ~GSYear, y = ~Richness, color = ~PermanentID) %>%
  add_lines() %>% 
  layout(showlegend = FALSE)


#### combine data ####

# combine native, invasive, control
nat_dat <- inv_plant %>%
  left_join(nat_plant2)

# identify missing data
nat_dat %>%
  filter(is.na(Richness)) %>%
  group_by(GSYear, CommonName) %>%
  summarize(Lakes = n_distinct(PermanentID))
# same number of lakes are missing each year
# dataset is just cut early because of no native
# plant sampling 2000-2001

# remove missing data
nat_dat2 <- nat_dat %>%
  filter(!is.na(Richness)) %>%
  mutate(across(ends_with("AvgPropCovered"), ~ .x * 100)) %>%
  rename_with(str_replace, pattern = "AvgPropCovered", replacement = "AvgPercCovered")

# split by species
hydr_dat <- filter(nat_dat2, CommonName == "Hydrilla") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))
wale_dat <- filter(nat_dat2, CommonName == "Water lettuce") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))
wahy_dat <- filter(nat_dat2, CommonName == "Water hyacinth") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))
torp_dat <- filter(nat_dat2, CommonName == "Torpedograss") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))
cubu_dat <- filter(nat_dat2, CommonName == "Cuban bulrush") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))
pagr_dat <- filter(nat_dat2, CommonName == "Para grass") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))


#### initial visualizations ####

# covariate correlations
nat_dat2 %>%
  select(CommonName, Lag1Treated, Lag1AvgPercCovered, SurveyorExperience, PrevRichness) %>%
  group_by(CommonName) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05) %>%
  data.frame()

# richness diff distribution
ggplot(nat_dat2, aes(x = RichnessDiff)) +
  geom_histogram() +
  facet_wrap(~ CommonName)

# coefficients and richness
ggplot(nat_dat2, aes(x = Lag1AvgPercCovered, y = RichnessDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(nat_dat2, aes(x = Lag1Treated, y = RichnessDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(nat_dat2, aes(x = PrevRichness, y = RichnessDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")


#### model-fitting functions ####

# data filter function
dat_mod_filt <- function(treat_col, inv_col, dat_in){
  
  dat_mod <- dat_in %>%
    filter(!is.na(Lag1Treated) & !is.na(Lag2Treated) & !is.na(Lag3Treated) & !is.na(Lag4Treated) & !is.na(Lag5Treated) & !is.na(Lag6Treated) & !is.na(!!sym(inv_col))) %>%
    mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
           Treated = !!sym(treat_col),
           AvgPercCovered = !!sym(inv_col),
           PrevRichness_c = PrevRichness - mean(PrevRichness),
           AvgPercCovered_c = AvgPercCovered - mean(AvgPercCovered))
  
  return(dat_mod)
  
}

# function to fit models
mod_fit <- function(dat_in, inv_col){
  
  # subset data
  dat_mod1 <- dat_mod_filt("Lag1Treated", inv_col, dat_in)
  dat_mod2 <- dat_mod_filt("Lag2Treated", inv_col, dat_in)
  dat_mod3 <- dat_mod_filt("Lag3Treated", inv_col, dat_in)
  dat_mod4 <- dat_mod_filt("Lag4Treated", inv_col, dat_in)
  dat_mod5 <- dat_mod_filt("Lag5Treated", inv_col, dat_in)
  dat_mod6 <- dat_mod_filt("Lag6Treated", inv_col, dat_in)
  
  # with initial richness
  mod1 <- feols(RichnessDiff ~ PrevRichness_c + AvgPercCovered_c + Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod1)
  mod2 <- update(mod1, data = dat_mod2)
  mod3 <- update(mod1, data = dat_mod3)
  mod4 <- update(mod1, data = dat_mod4)
  mod5 <- update(mod1, data = dat_mod5)
  mod6 <- update(mod1, data = dat_mod6)
  
  # without initial richness
  mod7 <- feols(RichnessDiff ~ AvgPercCovered_c + Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod1)
  mod8 <- update(mod7, data = dat_mod2)
  mod9 <- update(mod7, data = dat_mod3)
  mod10 <- update(mod7, data = dat_mod4)
  mod11 <- update(mod7, data = dat_mod5)
  mod12 <- update(mod7, data = dat_mod6)
  
  # with initial richness without management
  mod1b <- feols(RichnessDiff ~ PrevRichness_c + AvgPercCovered_c + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod1)
  mod2b <- update(mod1b, data = dat_mod2)
  mod3b <- update(mod1b, data = dat_mod3)
  mod4b <- update(mod1b, data = dat_mod4)
  mod5b <- update(mod1b, data = dat_mod5)
  mod6b <- update(mod1b, data = dat_mod6)
  
  # without initial richness without management
  mod7b <- feols(RichnessDiff ~ AvgPercCovered_c + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod1)
  mod8b <- update(mod7b, data = dat_mod2)
  mod9b <- update(mod7b, data = dat_mod3)
  mod10b <- update(mod7b, data = dat_mod4)
  mod11b <- update(mod7b, data = dat_mod5)
  mod12b <- update(mod7b, data = dat_mod6)
  
  
  # output
  return(list(mod1, mod2, mod3, mod4, mod5, mod6,
              mod7, mod8, mod9, mod10, mod11, mod12,
              mod1b, mod2b, mod3b, mod4b, mod5b, mod6b,
              mod7b, mod8b, mod9b, mod10b, mod11b, mod12b))
  
}


#### fit models ####

# fit models with all lags
hydr_lag1_mods <- mod_fit(hydr_dat, "Lag1AvgPercCovered")
wahy_lag1_mods <- mod_fit(wahy_dat, "Lag1AvgPercCovered")
wale_lag1_mods <- mod_fit(wale_dat, "Lag1AvgPercCovered")
torp_lag1_mods <- mod_fit(torp_dat, "Lag1AvgPercCovered")
cubu_lag1_mods <- mod_fit(cubu_dat, "Lag1AvgPercCovered")
pagr_lag1_mods <- mod_fit(pagr_dat, "Lag1AvgPercCovered")

hydr_lag6_mods <- mod_fit(hydr_dat, "Lag6AvgPercCovered")
wahy_lag6_mods <- mod_fit(wahy_dat, "Lag6AvgPercCovered")
wale_lag6_mods <- mod_fit(wale_dat, "Lag6AvgPercCovered")
torp_lag6_mods <- mod_fit(torp_dat, "Lag6AvgPercCovered")
pagr_lag6_mods <- mod_fit(pagr_dat, "Lag6AvgPercCovered")


# name models
names(hydr_lag1_mods) <- names(wahy_lag1_mods) <- names(wale_lag1_mods) <- names(torp_lag1_mods) <- names(cubu_lag1_mods) <- names(pagr_lag1_mods) <- names(hydr_lag6_mods) <- names(wahy_lag6_mods) <- names(wale_lag6_mods) <- names(torp_lag6_mods) <- names(pagr_lag6_mods) <- rep(c("1", "2", "3", "4", "5", "6"), 4)


#### coefficient figures and tables ####

# rename coefficients
coef_names <- c("SurveyorExperience_s" = "Surveyor experience",
                "PrevRichness_c" = "Initial richness",
                "Treated" = "Management", 
                "AvgPercCovered_c" = "Invasive PAC")

# ggplot function
plot_fun <- function(models){
  
  plot_out <- modelplot(models,
                        coef_map = coef_names,
                        background = list(geom_vline(xintercept = 0, color = "black",
                                                     size = 0.5, linetype = "dashed"))) +
    scale_color_viridis_d(direction = -1) +
    scale_x_continuous(labels = scale_fun_1) +
    def_theme_paper +
    theme(plot.title = element_text(size = 9))
  
  return(plot_out)
  
}

# panel plot function
panel_plot_fun <- function(mods1, mods2, mods3, 
                           spp1, spp2, spp3,
                           filename){
  
  # focal panels
  fig1 <- plot_fun(mods1) +
    labs(x = "",
         title = paste("(A)", spp1)) +
    def_theme_paper +
    theme(legend.position = "none")
  
  fig2 <- plot_fun(mods2) +
    labs(x = expression(paste("Estimate"%+-%" 95% CI", sep = "")),
         title = paste("(B)", spp2)) +
    theme(legend.position = "none",
          axis.text.y = element_blank())
  
  fig3 <- plot_fun(mods3) +
    labs(x = "",
         title = paste("(C)", spp3)) +
    theme(axis.text.y = element_blank(),
          legend.box.margin = margin(-10, 0, -10, -10)) +
    scale_color_viridis_d(direction = -1, name = "Management\nlag\n(years)") +
    guides(color = guide_legend(reverse = TRUE))
  
  comb_fig <- fig1 + fig2 + fig3 + plot_annotation(
    theme = theme(plot.margin = margin(5, -5, 0, -10),
                  plot.title = element_text(size = 10, hjust = 0.5)),
    title = "Effects on annual difference in native richness")
  
  ggsave(filename, comb_fig,
         device = "eps", width = 6.5, height = 3.5, units = "in")
  
}

# with initial richness and treatment
panel_plot_fun(hydr_lag1_mods[1:6], wahy_lag1_mods[1:6], wale_lag1_mods[1:6],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_native_richness_init_treat_lag1PAC_model.eps")
panel_plot_fun(cubu_lag1_mods[1:6], pagr_lag1_mods[1:6], torp_lag1_mods[1:6],
               "Cuban bulrush", "Para grass", "Torpedograss",
               "output/fwc_non_focal_native_richness_init_treat_lag1PAC_model.eps")
panel_plot_fun(hydr_lag6_mods[1:6], wahy_lag6_mods[1:6], wale_lag6_mods[1:6],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_native_richness_init_treat_lag6PAC_model.eps")
panel_plot_fun(hydr_lag6_mods[1:6], pagr_lag6_mods[1:6], torp_lag6_mods[1:6],
               "ignore", "Para grass", "Torpedograss",
               "output/fwc_non_focal_native_richness_init_treat_lag6PAC_model.eps")

# without initial richness, with treatment
panel_plot_fun(hydr_lag1_mods[7:12], wahy_lag1_mods[7:12], wale_lag1_mods[7:12],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_native_richness_treat_lag1PAC_model.eps")
panel_plot_fun(cubu_lag1_mods[7:12], pagr_lag1_mods[7:12], torp_lag1_mods[7:12],
               "Cuban bulrush", "Para grass", "Torpedograss",
               "output/fwc_non_focal_native_richness_treat_lag1PAC_model.eps")
panel_plot_fun(hydr_lag6_mods[7:12], wahy_lag6_mods[7:12], wale_lag6_mods[7:12],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_native_richness_treat_lag6PAC_model.eps")
panel_plot_fun(hydr_lag6_mods[7:12], pagr_lag6_mods[7:12], torp_lag6_mods[7:12],
               "ignore", "Para grass", "Torpedograss",
               "output/fwc_non_focal_native_richness_treat_lag6PAC_model.eps")

# with initial richness, without treatment
panel_plot_fun(hydr_lag1_mods[13:18], wahy_lag1_mods[13:18], wale_lag1_mods[13:18],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_native_richness_init_lag1PAC_model.eps")
panel_plot_fun(cubu_lag1_mods[13:18], pagr_lag1_mods[13:18], torp_lag1_mods[13:18],
               "Cuban bulrush", "Para grass", "Torpedograss",
               "output/fwc_non_focal_native_richness_init_lag1PAC_model.eps")
panel_plot_fun(hydr_lag6_mods[13:18], wahy_lag6_mods[13:18], wale_lag6_mods[13:18],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_native_richness_init_lag6PAC_model.eps")
panel_plot_fun(hydr_lag6_mods[13:18], pagr_lag6_mods[13:18], torp_lag6_mods[13:18],
               "ignore", "Para grass", "Torpedograss",
               "output/fwc_non_focal_native_richness_init_lag6PAC_model.eps")

# without initial richness and treatment
panel_plot_fun(hydr_lag1_mods[19:24], wahy_lag1_mods[19:24], wale_lag1_mods[19:24],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_native_richness_lag1PAC_model.eps")
panel_plot_fun(cubu_lag1_mods[19:24], pagr_lag1_mods[19:24], torp_lag1_mods[19:24],
               "Cuban bulrush", "Para grass", "Torpedograss",
               "output/fwc_non_focal_native_richness_lag1PAC_model.eps")
panel_plot_fun(hydr_lag6_mods[19:24], wahy_lag6_mods[19:24], wale_lag6_mods[19:24],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_native_richness_lag6PAC_model.eps")
panel_plot_fun(hydr_lag6_mods[19:24], pagr_lag6_mods[19:24], torp_lag6_mods[19:24],
               "ignore", "Para grass", "Torpedograss",
               "output/fwc_non_focal_native_richness_lag6PAC_model.eps")



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