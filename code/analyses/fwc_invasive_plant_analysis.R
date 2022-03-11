#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(GGally)
library(fixest) # FE models
library(modelsummary)
library(patchwork)
library(car)
library(glmmTMB)

# figure settings
source("code/settings/figure_settings.R")

# functions
source("code/generic-functions/proportion_transformations.R")
source("code/generic-functions/continuous_time_interval.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")


#### edit data ####

# combine datasets
inv_dat <- inv_plant %>%
  filter(!is.na(InitPercCovered)) %>% # need two consecutive years
  inner_join(inv_ctrl) %>%
  mutate(InitPACBin = case_when(InitPercCovered < 1 ~ "< 1%",
                                InitPercCovered >= 1 & InitPercCovered < 10 ~ "1%-10%",
                                TRUE ~ "≥ 10%") %>%
           fct_relevel("< 1%", "1%-10%"),
         PercChangeCovered = PropCovered * 100 - InitPercCovered,
         PropCoveredLogit = logit(PropCovered, adjust = 0.001))

# lakes that have the species present and been managed at least once
# by using EstAreaCoveredRaw_ha, there had to be more than one year per permanentID
perm_tax <- inv_dat %>%
  group_by(PermanentID, TaxonName) %>%
  summarize(Treatments = as.numeric(sum(Lag1Treated, na.rm = T) > 0),
            Established = as.numeric(sum(EstAreaCoveredRaw_ha) > 0)) %>%
  ungroup() %>%
  filter(Treatments > 0 & Established > 0)

# filter for lakes
inv_dat2 <- inv_dat %>%
  inner_join(perm_tax %>%
               select(PermanentID, TaxonName))

# will need surveyor experience
filter(inv_dat2, is.na(MinSurveyorExperience))
# none missing

# check data availability
inv_dat2 %>%
  filter(PropCovered > 0) %>%
  ggplot(aes(x = GSYear, y = PropCovered, color = PermanentID)) +
  geom_vline(xintercept = 2014) +
  geom_point() +
  facet_wrap(~ CommonName, scales = "free_y") +
  theme(legend.position = "none")
# Cuban bulrush is missing a lot of data before 2014
# Para grass and torpedograss start in 2000

# complete time intervals
inv_time_int <- inv_dat2 %>%
  select(GSYear, CommonName) %>% 
  unique() %>%
  mutate(out = pmap(., function(GSYear, CommonName) 
                    time_int_fun(year1 = GSYear, taxon = CommonName, dat_in = inv_dat2))) %>%
  unnest(cols = out)

inv_time_int %>%
  select(CommonName, years_out, data_points) %>%
  unique() %>%
  ggplot(aes(x = data_points)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

# select largest number of datapoints
inv_time_int2 <- inv_time_int %>%
  group_by(CommonName) %>%
  mutate(max_data_points = max(data_points)) %>%
  ungroup() %>%
  filter(data_points == max_data_points)

# filter dataset for lakes with complete surveys
inv_dat3 <- inv_dat2 %>%
  inner_join(inv_time_int2 %>%
               mutate(MinGSYear = GSYear,
                      MaxGSYear = GSYear + years_out - 1) %>% # 1 year added to years_out to count initial year
               select(CommonName, PermanentID, MinGSYear, MaxGSYear)) %>%
  filter(GSYear >= MinGSYear & GSYear <= MaxGSYear)

# taxa
inv_taxa <- sort(unique(inv_dat3$CommonName))

# loop through taxa
pdf("output/invasive_plant_continuous_time_series_by_taxon.pdf")

for(i in 1:length(inv_taxa)){
  
  # subset data
  subdat <- inv_dat3 %>% filter(CommonName == inv_taxa[i])
  subdat_ctrl <- subdat %>% filter(Lag1Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = PropCovered * 100, color = PermanentID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          labs(x = "Year", y = "Percent area covered", title = inv_taxa[i]) +
          def_theme_paper +
          theme(strip.text = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()

# combine water hyacinth and lettuce initial percent covered
floating_cover <- inv_dat %>%
  filter(CommonName %in% c("Water hyacinth", "Water lettuce")) %>%
  group_by(PermanentID, GSYear) %>%
  summarize(InitPercCovered = sum(InitPercCovered)) %>%
  ungroup() %>%
  mutate(InitPercCoveredFloat = if_else(InitPercCovered > 100, 100, InitPercCovered)) %>%
  expand_grid(CommonName = c("Water hyacinth", "Water lettuce")) %>%
  select(PermanentID, GSYear, InitPercCoveredFloat, CommonName)

# add floating cover
inv_dat4 <- inv_dat3 %>%
  left_join(floating_cover)

# split by species
hydr_dat <- filter(inv_dat4, CommonName == "Hydrilla") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))
wale_dat <- filter(inv_dat4, CommonName == "Water lettuce") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))
wahy_dat <- filter(inv_dat4, CommonName == "Water hyacinth") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))
torp_dat <- filter(inv_dat4, CommonName == "Torpedograss") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))
cubu_dat <- filter(inv_dat4, CommonName == "Cuban bulrush") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))
pagr_dat <- filter(inv_dat4, CommonName == "Para grass") %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience))


#### initial visualizations ####

# all and species-specific control correlated?
cor.test(~ Lag1Treated + Lag1AllTreated, data = hydr_dat) # 0.6
cor.test(~ Lag1Treated + Lag1AllTreated, data = wale_dat) # 0.8
cor.test(~ Lag1Treated + Lag1AllTreated, data = torp_dat) # 0.8 
cor.test(~ Lag1Treated + Lag1AllTreated, data = cubu_dat) # 0.3
cor.test(~ Lag1Treated + Lag1AllTreated, data = pagr_dat) # 0.2
# yes

# covariate correlations
hydr_dat %>%
  select(Lag1Treated, InitPercCovered, MinSurveyorExperience) %>%
  ggpairs()

wahy_dat %>%
  select(Lag1Treated, InitPercCovered, InitPercCoveredFloat, MinSurveyorExperience) %>%
  ggpairs() # one high cover value

wale_dat %>%
  select(Lag1Treated, InitPercCovered,InitPercCoveredFloat,  MinSurveyorExperience) %>%
  ggpairs() # one high cover value

torp_dat %>%
  select(Lag1Treated, InitPercCovered, MinSurveyorExperience) %>%
  ggpairs()

cubu_dat %>%
  select(Lag1Treated, InitPercCovered, MinSurveyorExperience) %>%
  ggpairs()

pagr_dat %>%
  select(Lag1Treated, InitPercCovered, MinSurveyorExperience) %>%
  ggpairs()

# response distributions
ggplot(inv_dat4, aes(x = PercChangeCovered)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_dat4, aes(x = PropCoveredLogit)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

# initial and responses
ggplot(inv_dat4, aes(x = InitPercCovered, y = PropCoveredLogit,
                    color = as.factor(Lag1Treated))) +
  geom_point() +
  facet_wrap(~ CommonName, scales = "free") +
  geom_smooth()
# positive, saturating
# paragrass has very few treatments

ggplot(inv_dat4, aes(x = InitPercCovered, y = PercChangeCovered,
                     color = as.factor(Lag1Treated))) +
  geom_point() +
  facet_wrap(~ CommonName, scales = "free") +
  geom_smooth()
# negative, linear

# treated and change in prop
ggplot(inv_dat4, aes(x = Lag6Treated, y = PropCoveredLogit)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ CommonName, scales = "free") 
# hydrilla slope is very positive
# weakens with increasing lag

ggplot(inv_dat4, aes(x = InitPercCovered, y = Lag1Treated)) +
  geom_point() +
  stat_smooth(method = "glm") +
  facet_wrap(~ CommonName, scales = "free") 
# yes for hydrilla and cuban bulrush (weaker with higher lag)
# no for water hyacinth and lettuce -- need to be combined?

ggplot(inv_dat4, aes(x = InitPercCoveredFloat, y = Lag1Treated)) +
  geom_point() +
  stat_smooth(method = "glm") +
  facet_wrap(~ CommonName, scales = "free") 
# helps water hyacinth, not water lettuce

# data availability for lags
inv_dat4 %>%
  select(PermanentID, GSYear, CommonName,
         Lag1Treated, Lag2Treated, Lag3Treated, Lag4Treated, Lag5Treated, Lag6Treated) %>%
  pivot_longer(cols = starts_with("Lag"),
               names_to = "Lag",
               values_to = "Treated") %>%
  filter(!is.na(Treated)) %>%
  ggplot(aes(x = Lag)) +
  geom_bar() +
  facet_wrap(~ CommonName, scales = "free") 

# initial PAC distributions
inv_dat4 %>%
  ggplot(aes(x = InitPercCovered)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free") 

# initial PAC distributions
# tried different thresholds before choosing these
ggplot(inv_dat4, aes(x = InitPACBin)) +
  geom_bar() +
  facet_wrap(~ CommonName, scales = "free") 
# para grass doesn't get above 10%

# before/after comparison
ggplot(inv_dat4, aes(x = Lag1Treated, y = PercChangeCovered)) +
  geom_hline(yintercept = 0, size = 0.25) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_dat4, aes(x = Lag6Treated, y = PercChangeCovered)) +
  geom_hline(yintercept = 0, size = 0.25) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ CommonName, scales = "free")


#### model-fitting functions ####

# data filter function
dat_mod_filt <- function(treat_col, dat_in){
  
  dat_mod <- dat_in %>%
    filter(!is.na(Lag1Treated) & !is.na(Lag2Treated) & !is.na(Lag3Treated) & !is.na(Lag4Treated) & !is.na(Lag5Treated) & !is.na(Lag6Treated)) %>%
    mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
           Treated = !!sym(treat_col))
  
  return(dat_mod)
  
}

# function to fit models
mod_fit <- function(dat_in){
  
  # subset data
  dat_mod1 <- dat_mod_filt("Lag1Treated", dat_in)
  dat_mod2 <- dat_mod_filt("Lag2Treated", dat_in)
  dat_mod3 <- dat_mod_filt("Lag3Treated", dat_in)
  dat_mod4 <- dat_mod_filt("Lag4Treated", dat_in)
  dat_mod5 <- dat_mod_filt("Lag5Treated", dat_in)
  dat_mod6 <- dat_mod_filt("Lag6Treated", dat_in)
  
  # fit models
  mod1 <- feols(PropCoveredLogit ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod1)
  mod2 <- feols(PropCoveredLogit ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod2)
  mod3 <- feols(PropCoveredLogit ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod3)
  mod4 <- feols(PropCoveredLogit ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod4)
  mod5 <- feols(PropCoveredLogit ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod5)
  mod6 <- feols(PropCoveredLogit ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod6)

  mod7 <- feols(PercChangeCovered ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod1)
  mod8 <- feols(PercChangeCovered ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod2)
  mod9 <- feols(PercChangeCovered ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod3)
  mod10 <- feols(PercChangeCovered ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod4)
  mod11 <- feols(PercChangeCovered ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod5)
  mod12 <- feols(PercChangeCovered ~ InitPercCovered * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod6)

  # output
  return(list(mod1, mod2, mod3, mod4, mod4, mod6,
              mod7, mod8, mod9, mod10, mod11, mod12))

}


#### fit models ####

# fit models with all lags
hydr_mods <- mod_fit(hydr_dat)
wahy_mods <- mod_fit(wahy_dat)
wale_mods <- mod_fit(wale_dat)
torp_mods <- mod_fit(torp_dat)
cubu_mods <- mod_fit(cubu_dat)
pagr_mods <- mod_fit(pagr_dat)

# separate by response
hydr_mods_prop <- hydr_mods[1:6]
wahy_mods_prop <- wahy_mods[1:6]
wale_mods_prop <- wale_mods[1:6]
torp_mods_prop <- torp_mods[1:6]
cubu_mods_prop <- cubu_mods[1:6]
pagr_mods_prop <- pagr_mods[1:6]

hydr_mods_change <- hydr_mods[7:12]
wahy_mods_change <- wahy_mods[7:12]
wale_mods_change <- wale_mods[7:12]
torp_mods_change <- torp_mods[7:12]
cubu_mods_change <- cubu_mods[7:12]
pagr_mods_change <- pagr_mods[7:12]

# name models
names(hydr_mods_prop) <- names(wahy_mods_prop) <- names(wale_mods_prop) <- names(torp_mods_prop) <- names(cubu_mods_prop) <- names(pagr_mods_prop) <- names(hydr_mods_change) <- names(wahy_mods_change) <- names(wale_mods_change) <- names(torp_mods_change) <- names(cubu_mods_change) <- names(pagr_mods_change) <- c("1", "2", "3", "4", "5", "6")


#### coefficient figures and tables ####

# rename coefficients
coef_names <- c("SurveyorExperience_s" = "Surveyor experience",
                "InitPercCovered:Treated" = "Management:PAC",
                "InitPercCovered" = "Initial PAC",
                "Treated" = "Management")

# ggplot function
plot_fun <- function(models){
  
  plot_out <- modelplot(models,
                        coef_map = coef_names,
                        background = list(geom_vline(xintercept = 0, color = "black",
                                                     size = 0.5, linetype = "dashed"))) +
    scale_color_viridis_d(direction = -1) +
    scale_x_continuous(labels = scale_fun_1) +
    def_theme_paper
  
  return(plot_out)
  
}

# focal panels
hydr_fig_prop <- plot_fun(hydr_mods_prop) +
  labs(x = "",
       title = "(A) Hydrilla") +
  def_theme_paper +
  theme(legend.position = "none")

wahy_fig_prop <- plot_fun(wahy_mods_prop) +
  labs(x = expression(paste("Estimate"%+-%" 95% CI", sep = "")),
       title = "(B) Water hyacinth") +
  theme(legend.position = "none",
        axis.text.y = element_blank())

wale_fig_prop <- plot_fun(wale_mods_prop) +
  labs(x = "",
       title = "(C) Water lettuce") +
  theme(axis.text.y = element_blank(),
        legend.box.margin = margin(-10, 0, -10, -10)) +
  scale_color_viridis_d(direction = -1, name = "Management\nlag\n(years)") +
  guides(color = guide_legend(reverse = TRUE))

hydr_fig_change <- plot_fun(hydr_mods_change) +
  labs(x = "",
       title = "(A) Hydrilla") +
  def_theme_paper +
  theme(legend.position = "none")

wahy_fig_change <- plot_fun(wahy_mods_change) +
  labs(x = expression(paste("Change in PAC"%+-%" 95% CI", sep = "")),
       title = "(B) Water hyacinth") +
  theme(legend.position = "none",
        axis.text.y = element_blank())

wale_fig_change <- plot_fun(wale_mods_change) +
  labs(x = "",
       title = "(C) Water lettuce") +
  theme(axis.text.y = element_blank(),
        legend.box.margin = margin(-10, 0, -10, -10)) +
  scale_color_viridis_d(direction = -1, name = "Management\nlag\n(years)") +
  guides(color = guide_legend(reverse = TRUE))

# combine focal panels
foc_fig_prop <- hydr_fig_prop + wahy_fig_prop + wale_fig_prop + plot_annotation(theme = theme(plot.margin = margin(0, -5, 0, -10)))
foc_fig_change <- hydr_fig_change + wahy_fig_change + wale_fig_change + plot_annotation(theme = theme(plot.margin = margin(0, -5, 0, -10)))

# save focal panels
ggsave("output/fwc_focal_invasive_prop_covered_treatment_model.eps", foc_fig_prop,
       device = "eps", width = 6.5, height = 3.5, units = "in")

ggsave("output/fwc_focal_invasive_PAC_change_treatment_model.eps", foc_fig_change,
       device = "eps", width = 6.5, height = 3.5, units = "in")

# non-focal panels
cubu_fig_prop <- plot_fun(cubu_mods_prop) +
  labs(x = "", title = "(A) Cuban bulrush") +
  theme(legend.position = "none")

pagr_fig_prop <- plot_fun(pagr_mods_prop) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       title = "(B) para grass") +
  theme(legend.position = "none",
        axis.text.y = element_blank())

torp_fig_prop <- plot_fun(torp_mods_prop) +
  scale_color_viridis_d(direction = -1, name = "Management\nlag\n(years)") +
  labs(x = "",
       title = "(C) torpedograss") +
  theme(axis.text.y = element_blank(),
        legend.box.margin = margin(-10, 0, -10, -10)) +
  guides(color = guide_legend(reverse = TRUE))

cubu_fig_change <- plot_fun(cubu_mods_change) +
  labs(x = "", title = "(A) Cuban bulrush") +
  theme(legend.position = "none")

pagr_fig_change <- plot_fun(pagr_mods_change) +
  labs(x = expression(paste("Change in PAC "%+-%" 95% CI", sep = "")),
       title = "(B) para grass") +
  theme(legend.position = "none",
        axis.text.y = element_blank())

torp_fig_change <- plot_fun(torp_mods_change) +
  scale_color_viridis_d(direction = -1, name = "Management\nlag\n(years)") +
  labs(x = "",
       title = "(C) torpedograss") +
  theme(axis.text.y = element_blank(),
        legend.box.margin = margin(-10, 0, -10, -10)) +
  guides(color = guide_legend(reverse = TRUE))

# combine figures
nonfoc_fig_prop <- cubu_fig_prop + pagr_fig_prop + torp_fig_prop + plot_annotation(theme = theme(plot.margin = margin(0, -5, 0, -10)))
nonfoc_fig_change <- cubu_fig_change + pagr_fig_change + torp_fig_change + plot_annotation(theme = theme(plot.margin = margin(0, -5, 0, -10)))

# save non-focal panels
ggsave("output/fwc_non_focal_invasive_prop_covered_treatment_model.eps", nonfoc_fig_prop,
       device = "eps", width = 6.5, height = 3.5, units = "in")
ggsave("output/fwc_non_focal_invasive_PAC_change_treatment_model.eps", nonfoc_fig_change,
       device = "eps", width = 6.5, height = 3.5, units = "in")


#### finalize models ####

# don't need surveyor experience
# only using year1 lag
hydr_dat1 <- hydr_dat %>% filter(!is.na(Lag1Treated))
wahy_dat1 <- wahy_dat %>% filter(!is.na(Lag1Treated))
wale_dat1 <- wale_dat %>% filter(!is.na(Lag1Treated))
cubu_dat1 <- cubu_dat %>% filter(!is.na(Lag1Treated))
torp_dat1 <- torp_dat %>% filter(!is.na(Lag1Treated))
pagr_dat1 <- pagr_dat %>% filter(!is.na(Lag1Treated))

# fit models
hydr_mod <- feols(PercChangeCovered ~ Lag1Treated | PermanentID + GSYear, data = hydr_dat1)
wahy_mod <- feols(PercChangeCovered ~ Lag1Treated | PermanentID + GSYear, data = wahy_dat1)
wale_mod <- feols(PercChangeCovered ~ Lag1Treated | PermanentID + GSYear, data = wale_dat1)
cubu_mod <- feols(PercChangeCovered ~ Lag1Treated | PermanentID + GSYear, data = cubu_dat1)
torp_mod <- feols(PercChangeCovered ~ Lag1Treated | PermanentID + GSYear, data = torp_dat1)
pagr_mod <- feols(PercChangeCovered ~ Lag1Treated | PermanentID + GSYear, data = pagr_dat1)

# export models
save(hydr_mod, file = "output/fwc_hydrilla_treatment_model.rda")
save(wahy_mod, file = "output/fwc_water_hyacinth_treatment_model.rda")
save(wale_mod, file = "output/fwc_water_lettuce_treatment_model.rda")
save(cubu_mod, file = "output/fwc_cuban_bulrush_treatment_model.rda")
save(pagr_mod, file = "output/fwc_para_grass_treatment_model.rda")
save(torp_mod, file = "output/fwc_torpedograss_lettuce_treatment_model.rda")


#### model prediction figures ####

# function to create predicted data
pred_dat <- function(dat_in, mod){
  
  # prediction dataset
  dat_out <- dat_in %>%
    select(PermanentID, GSYear) %>%
    expand_grid(Lag1Treated = seq(0, 1, length.out = 4)) %>%
    mutate(Pred = predict(mod, newdata = .),
           CommonName = unique(dat_in$CommonName))

  return(dat_out)
  
}

# combine datasets
foc_raw_dat <- hydr_dat1 %>%
  full_join(wahy_dat1) %>%
  full_join(wale_dat1)

foc_pred_dat <- pred_dat(hydr_dat1, hydr_mod) %>%
  full_join(pred_dat(wahy_dat1, wahy_mod)) %>%
  full_join(pred_dat(wale_dat1, wale_mod))

# figures
foc_pred_fig <- ggplot(foc_raw_dat, aes(x = Lag1Treated, y = PercChangeCovered)) +
  geom_hline(yintercept = 0, size = 0.25) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0,
               aes(color = CommonName)) +
  stat_summary(geom = "point", fun = "mean", size = 2,
               aes(color = CommonName)) +
  stat_summary(data = foc_pred_dat, aes(y = Pred, fill = CommonName),
               geom = "ribbon", fun.data = "mean_cl_boot", alpha = 0.5) +
  stat_summary(data = foc_pred_dat, aes(y = Pred, color = CommonName),
               geom = "line", fun = "mean", size = 1) +
  scale_x_continuous(breaks = c(0, 0.33, 0.67, 1)) +
  labs(x = "Management frequency", y = "Change in PAC") +
  def_theme_paper +
  facet_wrap(~ CommonName, scales = "free")

# doesn't allow ribbon to be saved for eps
ggsave("output/fwc_focal_invasive_PAC_change_treatment_prediction.eps", foc_pred_fig,
       device = "eps", width = 6.5, height = 2.5, units = "in")


#### older code ####

#### treatment, surveyor, clarity model figure ####

# # fit models
# hydr_qual_mods <- mod_qual_fit(hydr_qual_dat)
# wahy_qual_mods <- mod_qual_fit(wahy_qual_dat)
# wale_qual_mods <- mod_qual_fit(wale_qual_dat)
# torp_qual_mods <- mod_qual_fit(torp_qual_dat)
# alwe_qual_mods <- mod_qual_fit(alwe_qual_dat) # 'PrevPercCovered_c:Lag0Treated' removed because of collinearity
# # wita_qual_mods <- mod_qual_fit(wita_qual_dat)
# cubu_qual_mods <- mod_qual_fit(cubu_qual_dat)
# pagr_qual_mods <- mod_qual_fit(pagr_qual_dat)
# wafe_qual_mods <- mod_qual_fit(wafe_qual_dat) # all 'PrevPercCovered_c:LagXTreated' removed because of collinearity
# 
# # name models
# names(hydr_qual_mods) <- names(wahy_qual_mods) <- names(wale_qual_mods) <- names(torp_qual_mods) <- names(cubu_qual_mods) <- names(pagr_qual_mods) <- c("1", "2", "3", "4", "5")
# 
# # rename coefficients
# coef_qual_names <- c("Clarity_s" = "Water clarity", 
#                      "SurveyorExperience_s" = "Surveyor experience",
#                      "PrevPercCovered_c:Lag0Treated" = "Treatment x abundance",
#                      "PrevPercCovered_c:Lag1Treated" = "Treatment x abundance",
#                      "PrevPercCovered_c:Lag2Treated" = "Treatment x abundance",
#                      "PrevPercCovered_c:Lag3Treated" = "Treatment x abundance",
#                      "PrevPercCovered_c:Lag4Treated" = "Treatment x abundance",
#                      "PrevPercCovered_c:Lag5Treated" = "Treatment x abundance",
#                      "PrevPercCovered_c" = "Initial abundance (%)",
#                      "Lag0Treated" = "Treatment frequency",
#                      "Lag1Treated" = "Treatment frequency",
#                      "Lag2Treated" = "Treatment frequency",
#                      "Lag3Treated" = "Treatment frequency",
#                      "Lag4Treated" = "Treatment frequency",
#                      "Lag5Treated" = "Treatment frequency")
# 
# # panels
# hydr_qual_fig <- modelplot(hydr_qual_mods,
#                       coef_map = coef_qual_names,
#                       background = list(geom_vline(xintercept = 0, color = "black",
#                                                    size = 0.5, linetype = "dashed"))) +
#   scale_color_viridis_d(direction = -1) +
#   labs(x = "",
#        title = "(A) hydrilla") +
#   def_theme_paper +
#   theme(legend.position = "none")
# 
# wahy_qual_fig <- modelplot(wahy_qual_mods,
#                       coef_map = coef_qual_names,
#                       background = list(geom_vline(xintercept = 0, color = "black",
#                                                    size = 0.5, linetype = "dashed"))) +
#   scale_color_viridis_d(direction = -1) +
#   labs(x = "",
#        title = "(B) water hyacinth") +
#   def_theme_paper +
#   theme(legend.position = "none",
#         axis.text.y = element_blank())
# 
# wale_qual_fig <- modelplot(wale_qual_mods,
#                       coef_map = coef_qual_names,
#                       background = list(geom_vline(xintercept = 0, color = "black",
#                                                    size = 0.5, linetype = "dashed"))) +
#   scale_color_viridis_d(direction = -1) +
#   labs(x = "",
#        title = "(C) water lettuce") +
#   def_theme_paper +
#   theme(axis.text.y = element_blank(),
#         legend.position = "none")
# 
# cubu_qual_fig <- modelplot(cubu_qual_mods,
#                       coef_map = coef_qual_names,
#                       background = list(geom_vline(xintercept = 0, color = "black",
#                                                    size = 0.5, linetype = "dashed"))) +
#   scale_color_viridis_d(direction = -1) +
#   labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
#        title = "(D) Cuban bulrush") +
#   def_theme_paper +
#   theme(legend.position = "none",
#         axis.text.y = element_blank())
# 
# pagr_qual_fig <- modelplot(pagr_qual_mods,
#                       coef_map = coef_qual_names,
#                       background = list(geom_vline(xintercept = 0, color = "black",
#                                                    size = 0.5, linetype = "dashed"))) +
#   scale_color_viridis_d(direction = -1) +
#   labs(x = "",
#        title = "(E) para grass") +
#   def_theme_paper +
#   theme(legend.position = "none",
#         axis.text.y = element_blank())
# 
# torp_qual_fig <- modelplot(torp_qual_mods,
#                       coef_map = coef_qual_names,
#                       background = list(geom_vline(xintercept = 0, color = "black",
#                                                    size = 0.5, linetype = "dashed"))) +
#   scale_color_viridis_d(direction = -1) +
#   labs(x = "",
#        title = "(F) torpedograss") +
#   def_theme_paper +
#   theme(legend.position = "none",
#         axis.text.y = element_blank())
# 
# # wita_qual_fig <- modelplot(wita_qual_mods,
# #                       coef_map = coef_qual_names,
# #                       background = list(geom_vline(xintercept = 0, color = "black",
# #                                                    size = 0.5, linetype = "dashed"))) +
# #   scale_color_viridis_d(name = "Years\nof\ndata", direction = -1) +
# #   labs(x = "",
# #        title = "(G) wild taro") +
# #   def_theme_paper +
# #   theme(axis.text.y = element_blank(),
# #         legend.box.margin = margin(-10, 0, -10, -10)) +
# #   guides(color = guide_legend(reverse = TRUE))
# 
# # combine figures
# pdf("output/fwc_invasive_plant_treatment_clarity_model.pdf", width = 13, height = 5)
# plot_grid(hydr_qual_fig, wahy_qual_fig, wale_qual_fig, cubu_qual_fig, pagr_qual_fig, torp_qual_fig, 
#           nrow = 1,
#           rel_widths = c(1, rep(0.57, 5), 0.75))
# dev.off()


#### model prediction figure ####

# function to format data
pred_fig <- function(dat_in, mod_list){
  
  # subset data
  raw_dat <- dat_in %>%
    filter(!is.na(Lag0Treated) & !is.na(Lag1Treated) & !is.na(Lag2Treated) & !is.na(Lag3Treated) & !is.na(Lag4Treated) & !is.na(Lag5Treated) & !is.na(SurveyorExperience)) %>%
    mutate(PrevPercCovered_c = PrevPercCovered - mean(PrevPercCovered))
  
  # prediction dataset
  pred_dat <- tibble(PrevPercCovered = c(0, 15, 30)) %>%
    mutate(PrevPercCovered_c = PrevPercCovered - mean(raw_dat$PrevPercCovered),
           Clarity_s = 0) %>%
    expand_grid(tibble(Lag5Treated = c(0, 1/6, 2/6, 1/2, 4/6, 5/6, 1))) %>%
    expand_grid(raw_dat %>%
                  select(PermanentID, GSYear) %>%
                  unique()) %>%
    mutate(SurveyorExperience_s = 0,
           PrevPercCovered_f = as.factor(PrevPercCovered) %>%
             fct_rev()) %>%
    mutate(Pred = predict(mod_list[[6]], newdata = .))
  
  fig_out <- ggplot(pred_dat, aes(x = Lag5Treated, y = Pred)) +
    geom_point(data = filter(raw_dat, PrevPercCovered <= 1), 
               alpha = 0.2, position = position_jitter(width = 0.01),
               shape = 21, color = "transparent",
               aes(y = LogRatioCovered, fill = PrevPercCovered)) +
    geom_point(data = filter(raw_dat, PrevPercCovered > 1), 
               alpha = 0.5, position = position_jitter(width = 0.01),
               shape = 21, color = "transparent",
               aes(y = LogRatioCovered, fill = PrevPercCovered)) +
    stat_summary(geom = "line", size = 1.3, fun = "median", aes(color = PrevPercCovered_f)) +
    scale_color_viridis_d(name = "Modeled\ninitial\nabundance (%)", begin = 0.7) +
    scale_fill_viridis_c(name = "Observed\ninitial\nabundance (%)", direction = -1) +
    labs(x = "", y = "Change in abundance") +
    def_theme_paper
  
  return(fig_out)
  
}

# figures
hydr_pred_fig <- pred_fig(hydr_dat, hydr_mods) +
  labs(title = "(A) hydrilla") +
  theme(legend.position = "none")
wahy_pred_fig <- pred_fig(wahy_dat, wahy_mods) +
  labs(title = "(B) water hyacinth", x = "Treatment frequency") +
  theme(legend.position = "none",
        axis.title.y = element_blank())
wale_pred_fig <- pred_fig(wale_dat, wale_mods) +
  labs(title = "(C) water lettuce") +
  theme(axis.title.y = element_blank())

hydr_qual_pred_fig <- pred_fig(hydr_qual_dat, hydr_qual_mods) +
  labs(title = "(A) hydrilla") +
  theme(legend.position = "none")
wahy_qual_pred_fig <- pred_fig(wahy_qual_dat, wahy_qual_mods) +
  labs(title = "(B) water hyacinth", x = "Treatment frequency") +
  theme(legend.position = "none",
        axis.title.y = element_blank())
wale_qual_pred_fig <- pred_fig(wale_qual_dat, wale_qual_mods) +
  labs(title = "(C) water lettuce") +
  theme(axis.title.y = element_blank())

# combine figures
pdf("output/fwc_invasive_plant_treatment_predictions.pdf", width = 6.5, height = 3)
plot_grid(hydr_pred_fig, wahy_pred_fig, wale_pred_fig,
          nrow = 1,
          rel_widths = c(1, 0.9, 1.25))
dev.off()

pdf("output/fwc_invasive_plant_treatment_clarity_predictions.pdf", width = 6.5, height = 3)
plot_grid(hydr_qual_pred_fig, wahy_qual_pred_fig, wale_qual_pred_fig,
          nrow = 1,
          rel_widths = c(1, 0.9, 1.25))
dev.off()



### model comparison ####

# compare models (use delta 4 as cut-off)
mod_lag_comp(hydr_mods) # 4, 5 similar
mod_lag_comp(wahy_mods) # all except 1 similar
mod_lag_comp(wale_mods) # 5 best
mod_lag_comp(cubu_mods) # 2 best
mod_lag_comp(pagr_mods) # all similar
mod_lag_comp(torp_mods) # 0 and 2 best
mod_lag_comp(wita_mods) # 1-5 similar

mod_lag_comp(hydr_qual_mods) # all except 1 similar
mod_lag_comp(wahy_qual_mods) # 3, 4, 5 similar
mod_lag_comp(wale_qual_mods) # 4, 5 similar
mod_lag_comp(cubu_qual_mods) # 0, 1 best
mod_lag_comp(pagr_qual_mods) # all similar
mod_lag_comp(torp_qual_mods) # all similar
mod_lag_comp(wita_qual_mods) # all similar

# output AIC values
aic_out <- tibble(Lag = mod_lag_comp(hydr_mods)$Lag,
                  Hydrilla = mod_lag_comp(hydr_mods)$deltaAIC,
                  Waterhyacinth = mod_lag_comp(wahy_mods)$deltaAIC,
                  Waterlettuce = mod_lag_comp(wale_mods)$deltaAIC,
                  Cubanbulrush = mod_lag_comp(cubu_mods)$deltaAIC,
                  Paragrass = mod_lag_comp(pagr_mods)$deltaAIC,
                  Torpedograss = mod_lag_comp(torp_mods)$deltaAIC,
                  Wildtaro = mod_lag_comp(wita_mods)$deltaAIC)

write_csv(aic_out, "output/fwc_invasive_plant_delta_aic.csv")

aic_qual_out <- tibble(Lag = mod_lag_comp(hydr_qual_mods)$Lag,
                       Hydrilla = mod_lag_comp(hydr_qual_mods)$deltaAIC,
                       Waterhyacinth = mod_lag_comp(wahy_qual_mods)$deltaAIC,
                       Waterlettuce = mod_lag_comp(wale_qual_mods)$deltaAIC,
                       Cubanbulrush = mod_lag_comp(cubu_qual_mods)$deltaAIC,
                       Paragrass = mod_lag_comp(pagr_qual_mods)$deltaAIC,
                       Torpedograss = mod_lag_comp(torp_qual_mods)$deltaAIC,
                       Wildtaro = mod_lag_comp(wita_qual_mods)$deltaAIC)

write_csv(aic_qual_out, "output/fwc_invasive_plant_clarity_delta_aic.csv")
