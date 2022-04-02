#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(inspectdf) # inspect_cor
library(modelsummary) # modelplot
library(patchwork) # combining figures
library(plm) # panel data models
library(pglm) # panel data models
library(sandwich) # vcovHC
library(lmtest) # coeftest
library(glmmTMB) # random effects
library(pals) # color palettes               

# figure settings
source("code/settings/figure_settings.R")

# functions
source("code/generic-functions/model_structure_comparison.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_analysis_formatted.csv") # plant and control data, continuous data
lw_chl <- read_csv("intermediate-data/LW_chlorophyll_formatted.csv")
lwwa_chl <- read_csv("intermediate-data/LW_water_atlas_chlorophyll_formatted.csv")


#### edit data ####

# require longest reasonable lag
# longer lags mean fewer years are needed from chlorophyll dataset
# add max years column
inv_plant2 <- inv_plant %>%
  filter((CommonName != "Cuban bulrush" & !is.na(Lag6AvgPropCovered) & !is.na(Lag6Treated)) |
           (CommonName == "Cuban bulrush" & !is.na(Lag3AvgPropCovered) & !is.na(Lag3Treated))) %>%
  group_by(CommonName) %>%
  mutate(maxYears = n_distinct(GSYear)) %>%
  ungroup()

# check
inv_plant2 %>%
  group_by(CommonName) %>%
  summarize(maxYears = unique(maxYears),
            lastYear = max(GSYear),
            firstYear = min(GSYear),
            locations = n_distinct(PermanentID))

# year distribution
ggplot(lwwa_chl, aes(x = GSYear)) +
  geom_bar() +
  facet_wrap(~ Quarter)
# some are low for 2019

# combine chlorophyll, invasive, control
# select waterbodies sampled throughout
chl_dat <- lwwa_chl %>%
  filter(!is.na(PrevValue)) %>%
  inner_join(inv_plant2 %>%
               filter(GSYear >= 2005 & GSYear < 2019) %>% # year cut-offs from data exploration below
               group_by(CommonName) %>%
               mutate(maxYears = n_distinct(GSYear)) %>% # recalculate max years
               ungroup()) %>%
  group_by(CommonName, PermanentID, Quarter) %>%
  mutate(nYears = n_distinct(GSYear)) %>% # years per waterbody
  ungroup() %>%
  filter(nYears == maxYears) %>%
  mutate(ValueDiff = QualityValue - PrevValue,  # change over time
         across(ends_with("AvgPropCovered"), ~ .x * 100)) %>%
  rename_with(str_replace, pattern = "AvgPropCovered", replacement = "AvgPercCovered")

# sample sizes
chl_dat %>%
  group_by(CommonName, Quarter) %>%
  summarize(Years = n_distinct(GSYear),
            Waterbodies = n_distinct(PermanentID)) %>%
  data.frame()
# torpedograss was completely missing before date cut-offs
# its distribution is highly concentrated in center of state
# https://nas.er.usgs.gov/viewer/omap.aspx?SpeciesID=1124
# para grass only has one or two waterbodies, but it only has 10 total

# why are there no data matching with torpedograss?
inv_plant2 %>%
  filter(CommonName == "Torpedograss") %>%
  select(PermanentID, GSYear) %>%
  unique() %>%
  inner_join(lwwa_chl %>%
               select(PermanentID, GSYear) %>%
               unique()) %>%
  group_by(GSYear) %>%
  count() # a lot of waterbodies each year

inv_plant2 %>%
  filter(CommonName == "Torpedograss") %>%
  select(PermanentID, GSYear) %>%
  unique() %>%
  inner_join(lwwa_chl %>%
               select(PermanentID, GSYear) %>%
               unique()) %>%
  group_by(PermanentID) %>%
  summarize(Years = n_distinct(GSYear),
            firstYear = min(GSYear),
            lastYear = max(GSYear),
            yearDiff = lastYear - firstYear + 1) %>%
  data.frame() 
# most are missing sampling in 2004

# taxa
inv_taxa <- sort(unique(chl_dat$CommonName))

# loop through taxa
pdf("output/chlorophyll_continuous_time_series_by_taxon.pdf")

for(i in 1:length(inv_taxa)){
  
  # subset data
  subdat <- chl_dat %>% filter(CommonName == inv_taxa[i])
  subdat_ctrl <- subdat %>% filter(Lag1Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = QualityValue, color = PermanentID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          facet_wrap(~ Quarter) + 
          labs(x = "Year", y = "Chlorophyll a (ug/L)", title = inv_taxa[i]) +
          def_theme_paper +
          theme(plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()

# remove paragrass (only 1-2 waterbodies)
# use same waterbodies in all 4 quarters
chl_dat2 <- chl_dat %>%
  filter(CommonName != "Para grass") %>%
  group_by(CommonName, PermanentID, GSYear) %>%
  mutate(nQuart = n_distinct(Quarter)) %>%
  ungroup() %>%
  filter(nQuart == 4)

# sample sizes
chl_dat2 %>%
  group_by(CommonName, Quarter) %>%
  summarize(Years = n_distinct(GSYear),
            Waterbodies = n_distinct(PermanentID),
            N = Years * Waterbodies) %>%
  data.frame()

# split by species
hydr_dat <- filter(chl_dat2, CommonName == "Hydrilla")
wale_dat <- filter(chl_dat2, CommonName == "Water lettuce")
wahy_dat <- filter(chl_dat2, CommonName == "Water hyacinth")
torp_dat <- filter(chl_dat2, CommonName == "Torpedograss")
cubu_dat <- filter(chl_dat2, CommonName == "Cuban bulrush")


#### initial visualizations ####

# covariate correlations
chl_dat2 %>%
  select(Quarter, CommonName, 
         Lag3Treated, Lag3AvgPercCovered, MinSurveyorExperience, PrevValue) %>%
  group_by(Quarter, CommonName) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & abs(corr) >= 0.4) %>%
  data.frame()
# water lettuce: prev value is correlated with PAC
# water lettuce: surveyor experience is correlated with treat

# response distributions
ggplot(chl_dat2, aes(x = ValueDiff)) +
  geom_histogram() +
  facet_grid(CommonName ~ Quarter, scales = "free")
# normal and relatively wide around zero

ggplot(chl_dat2, aes(x = QualityValue)) +
  geom_histogram() +
  facet_grid(CommonName ~ Quarter, scales = "free")
# skewed

# coefficients and QualityValue
ggplot(chl_dat2, aes(x = Lag3AvgPercCovered, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# generally close to zero

ggplot(chl_dat2, aes(x = Lag3AvgPercCovered, y = QualityValue)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# strong positive for water hyacinth and lettuce
# negative for hydrilla and torpedograss
# switches over time for Cuban bulrush

ggplot(chl_dat2, aes(x = Lag3Treated, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# negligible

ggplot(chl_dat2, aes(x = Lag1Treated, y = QualityValue)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# positive for floating plants in quarter 3
# negative for torpedograss in quarters 2 and 3

ggplot(chl_dat2, aes(x = PrevValue, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(CommonName ~ Quarter, scales = "free")
# consistently negative

ggplot(chl_dat2, aes(x = PrevValue, y = QualityValue)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(CommonName ~ Quarter, scales = "free")
# consistently positive


#### evaluate model structure ####

# function to fit models for each species
mod_structure_fits <- function(dat_in){
  
  # create fixed effects data frame
  # choose a quarter
  # cannot include quarter as an interaction in fixed-effect models because then there are duplicate
  # individual-time rows
  dat_fix <- dat_in %>%
    filter(Quarter == 3) %>%
    mutate(PrevValue_s = (PrevValue - mean(PrevValue)) / sd(PrevValue))  %>%
    ungroup() %>%
    pdata.frame(index = c("PermanentID", "GSYear"))
  # each waterbody is an individual
  
  # simple lm
  mod_lm <- lm(QualityValue ~ Lag1AvgPercCovered + Lag1Treated, data = dat_in)
  
  # random effects
  mod_ran_loc <- glmmTMB(QualityValue ~ Lag1AvgPercCovered + Lag1Treated + (1|PermanentID), data = dat_in)
  mod_ran_yr <- glmmTMB(QualityValue ~ Lag1AvgPercCovered + Lag1Treated + (1|GSYear), data = dat_in)
  mod_ran_loc_yr <- glmmTMB(QualityValue ~ Lag1AvgPercCovered + Lag1Treated + (1|PermanentID) + (1|GSYear), data = dat_in)
  
  # fixed effects
  mod_fix_loc <- plm(QualityValue ~ Lag1AvgPercCovered + Lag1Treated, data = dat_fix,
                      model = "within")
  mod_fix_loc_yr <- plm(QualityValue ~ Lag1AvgPercCovered + Lag1Treated, data = dat_fix,
                         model = "within", effect = "twoways")
  
  # use initial richness to account for reverse causality
  mod_init_fix_loc_yr <- plm(QualityValue ~ PrevValue_s + Lag1AvgPercCovered + Lag1Treated, data = dat_fix,
                              model = "within", effect = "twoways")
  
  # use richness difference to account for reverse causality
  mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag1AvgPercCovered + Lag1Treated, data = dat_fix,
                             index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
  
  # richness difference and initial richness
  mod_init_diff_fix_loc_yr <- plm(ValueDiff ~ PrevValue_s + Lag1AvgPercCovered + Lag1Treated, data = dat_fix,
                                  index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
  
  # return list of models
  return(list(lm = mod_lm,
              ran_loc = mod_ran_loc,
              ran_yr = mod_ran_yr,
              ran_loc_yr = mod_ran_loc_yr,
              fix_loc = mod_fix_loc,
              fix_loc_yr = mod_fix_loc_yr,
              init_fix_loc_yr = mod_init_fix_loc_yr,
              diff_fix_loc_yr = mod_diff_fix_loc_yr,
              init_diff_fix_loc_yr = mod_init_diff_fix_loc_yr))
  
}

# fit models for each species
hydr_mod_struc <- mod_structure_fits(hydr_dat)
wahy_mod_struc <- mod_structure_fits(wahy_dat)
wale_mod_struc <- mod_structure_fits(wale_dat)

# compare model estimates
hydr_mod_comp <- mod_structure_comp(simp_mods = hydr_mod_struc[1], 
                                    ran_mods = hydr_mod_struc[2:4],
                                    fix_mods = hydr_mod_struc[5:9])
wahy_mod_comp <- mod_structure_comp(simp_mods = wahy_mod_struc[1], 
                                    ran_mods = wahy_mod_struc[2:4],
                                    fix_mods = wahy_mod_struc[5:9]) 
wale_mod_comp <- mod_structure_comp(simp_mods = wale_mod_struc[1], 
                                    ran_mods = wale_mod_struc[2:4],
                                    fix_mods = wale_mod_struc[5:9]) 

# combine species
mod_comp <- hydr_mod_comp %>%
  mutate(Species = "hydrilla") %>%
  full_join(wahy_mod_comp %>%
              mutate(Species = "water hyacinth")) %>%
  full_join(wale_mod_comp %>%
              mutate(Species = "water lettuce")) %>%
  mutate(coefficients = str_replace(coefficients, "Lag1Treated", "management"),
         coefficients = str_replace(coefficients, "Lag1AvgPercCovered", "PAC"),
         across(!c(coefficients, Species), ~ round(.x, digits = 3))) %>%
  relocate(Species)

write_csv(mod_comp, "output/fwc_chlorophyll_model_structure_comparison.csv")

# model comparison notes:
# simple model: PAC reduces chlorophyll (hydrilla) or increases it (floating),
# management increases chlorophyll (hydrilla and hyacinth) or reduces it (lettuce)
# location random effect: hydrilla PAC now increases chlorophyll and management reduces
# floating PAC and management increase it
# year random effect: similar to simple model
# location + year random: similar to location random
# location fixed effect: similar to location random
# location + year fixed: similar to location random
# initial quality: same direction as location random, smaller estimates
# difference response: hydrilla PAC and management increase chlorophyll change
# floating PAC increases it, management decreases it
# difference response + initial quality: estimates almost identical to initial quality

# test fixed effects (seems like year isn't necessary)
# have to refit because data need to be accessible (not "dat_fix")
hydr_mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag1AvgPercCovered + Lag1Treated, 
                                data = filter(hydr_dat, Quarter == 3), 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(hydr_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(hydr_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

wahy_mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag1AvgPercCovered + Lag1Treated, 
                                data = filter(wahy_dat, Quarter == 3), 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(wahy_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(wahy_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

wale_mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag1AvgPercCovered + Lag1Treated, 
                                data = filter(wale_dat, Quarter == 3), 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(wale_mod_diff_fix_loc_yr, effect = "time", type = "bp") # marginal
plmtest(wale_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

# use annual difference without initial value and with waterbody and year fixed effects


#### model-fitting functions ####

# data filter function
dat_mod_filt <- function(treat_col, inv_col, dat_in){
  
  dat_mod <- dat_in %>%
    mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
           Treated = !!sym(treat_col),
           AvgPercCovered = !!sym(inv_col),
           AvgPercCovered_c = AvgPercCovered - mean(AvgPercCovered))
  
  return(dat_mod)
  
}

# function to fit models
mod_fit <- function(dat_in){
  
  # focal species
  foc_sp <- unique(dat_in$CommonName)
  
  # subset data
  dat_mod1 <- dat_mod_filt("Lag1Treated", "Lag1AvgPercCovered", dat_in)
  dat_mod2 <- dat_mod_filt("Lag2Treated", "Lag2AvgPercCovered", dat_in)
  dat_mod3 <- dat_mod_filt("Lag3Treated", "Lag3AvgPercCovered", dat_in)
  dat_mod4 <- dat_mod_filt("Lag4Treated", "Lag4AvgPercCovered", dat_in)
  dat_mod5 <- dat_mod_filt("Lag5Treated", "Lag5AvgPercCovered", dat_in)
  dat_mod6 <- dat_mod_filt("Lag6Treated", "Lag6AvgPercCovered", dat_in)
  
  # fit models
  if(foc_sp == "Cuban bulrush") {
    
    mod1 <- plm(ValueDiff ~ AvgPercCovered_c + Treated + SurveyorExperience_s, data = dat_mod1,
                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
    mod2 <- update(mod1, data = dat_mod2)
    mod3 <- update(mod1, data = dat_mod3)
    
    mods_out <- list(mod1, mod2, mod3)
    
  } else {
    
    mod1 <- plm(ValueDiff ~ AvgPercCovered_c + Treated + SurveyorExperience_s, data = dat_mod1,
                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
    mod2 <- update(mod1, data = dat_mod2)
    mod3 <- update(mod1, data = dat_mod3)
    mod4 <- update(mod1, data = dat_mod4)
    mod5 <- update(mod1, data = dat_mod5)
    mod6 <- update(mod1, data = dat_mod6)
    
    mods_out <- list(mod1, mod2, mod3, mod4, mod5, mod6)
    
  }
  
  # output
  return(mods_out)
  
}


#### fit models ####

# fit models with all lags
hydr_mods_q1 <- mod_fit(filter(hydr_dat, Quarter == 1))
wahy_mods_q1 <- mod_fit(filter(wahy_dat, Quarter == 1))
wale_mods_q1 <- mod_fit(filter(wale_dat, Quarter == 1))
torp_mods_q1 <- mod_fit(filter(torp_dat, Quarter == 1))
cubu_mods_q1 <- mod_fit(filter(cubu_dat, Quarter == 1))

hydr_mods_q2 <- mod_fit(filter(hydr_dat, Quarter == 2))
wahy_mods_q2 <- mod_fit(filter(wahy_dat, Quarter == 2))
wale_mods_q2 <- mod_fit(filter(wale_dat, Quarter == 2))
torp_mods_q2 <- mod_fit(filter(torp_dat, Quarter == 2))
cubu_mods_q2 <- mod_fit(filter(cubu_dat, Quarter == 2))

hydr_mods_q3 <- mod_fit(filter(hydr_dat, Quarter == 3))
wahy_mods_q3 <- mod_fit(filter(wahy_dat, Quarter == 3))
wale_mods_q3 <- mod_fit(filter(wale_dat, Quarter == 3))
torp_mods_q3 <- mod_fit(filter(torp_dat, Quarter == 3))
cubu_mods_q3 <- mod_fit(filter(cubu_dat, Quarter == 3))

hydr_mods_q4 <- mod_fit(filter(hydr_dat, Quarter == 4))
wahy_mods_q4 <- mod_fit(filter(wahy_dat, Quarter == 4))
wale_mods_q4 <- mod_fit(filter(wale_dat, Quarter == 4))
torp_mods_q4 <- mod_fit(filter(torp_dat, Quarter == 4))
cubu_mods_q4 <- mod_fit(filter(cubu_dat, Quarter == 4))

# name models
names(hydr_mods_q1) <- names(wahy_mods_q1) <- names(wale_mods_q1) <- names(torp_mods_q1) <- names(hydr_mods_q2) <- names(wahy_mods_q2) <- names(wale_mods_q2) <- names(torp_mods_q2) <- names(hydr_mods_q3) <- names(wahy_mods_q3) <- names(wale_mods_q3) <- names(torp_mods_q3) <- names(hydr_mods_q4) <- names(wahy_mods_q4) <- names(wale_mods_q4) <- names(torp_mods_q4) <- c("1", "2", "3", "4", "5", "6")
  
names(cubu_mods_q1) <- names(cubu_mods_q2) <- names(cubu_mods_q3) <- names(cubu_mods_q4) <- c("1", "2", "3")


#### coefficient figures and tables ####

# rename coefficients
coef_names <- c("SurveyorExperience_s" = "Surveyor experience",
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
    scale_color_viridis_d(direction = -1, name = "Lag\n(years)") +
    guides(color = guide_legend(reverse = TRUE))
  
  comb_fig <- fig1 + fig2 + fig3 + plot_annotation(
    theme = theme(plot.margin = margin(5, -5, 0, -10),
                  plot.title = element_text(size = 10, hjust = 0.5)),
    title = "Effects on annual difference in chlorophyll a")
  
  ggsave(filename, comb_fig,
         device = "eps", width = 6.5, height = 3.5, units = "in")
  
}

panel_plot_non_foc_fun <- function(mods1, mods2,
                                   spp1, spp2,
                                   filename){
  
  # focal panels
  fig1 <- plot_fun(mods1) +
    labs(x = "",
         title = paste("(A)", spp1)) +
    theme(legend.position = "none")
  
  fig2 <- plot_fun(mods2) +
    labs(x = "",
         title = paste("(B)", spp2)) +
    theme(axis.text.y = element_blank(),
          legend.box.margin = margin(-10, 0, -10, -10)) +
    scale_color_viridis_d(direction = -1, name = "Lag\n(years)") +
    guides(color = guide_legend(reverse = TRUE))
  
  comb_fig <- (fig1 + fig2) + 
    # labs(x = expression(paste("Estimate"%+-%" 95% CI", sep = ""))) +
    # xlab(label = "Estimate2") +
    plot_annotation(
      caption = expression(paste("Estimate"%+-%" 95% CI", sep = "")),
      theme = theme(plot.caption = element_text(size = 9, color="black", hjust = 0.6, vjust = 10),
                    plot.margin = margin(5, -5, -5, -10),
                    plot.title = element_text(size = 10, hjust = 0.5)),
      title = "Effects on annual difference in chlorophyll a")
  
  ggsave(filename, comb_fig,
         device = "eps", width = 4.7, height = 3.5, units = "in")
  
}

# figures
panel_plot_fun(hydr_mods_q1, wahy_mods_q1, wale_mods_q1,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_chlorophyll_quarter1_diff_model.eps")
panel_plot_fun(hydr_mods_q2, wahy_mods_q2, wale_mods_q2,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_chlorophyll_quarter2_diff_model.eps")
panel_plot_fun(hydr_mods_q3, wahy_mods_q3, wale_mods_q3,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_chlorophyll_quarter3_diff_model.eps")
panel_plot_fun(hydr_mods_q4, wahy_mods_q4, wale_mods_q4,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_chlorophyll_quarter4_diff_model.eps")

panel_plot_non_foc_fun(cubu_mods_q1, torp_mods_q1,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_chlorophyll_quarter1_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q2, torp_mods_q2,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_chlorophyll_quarter2_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q3, torp_mods_q3,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_chlorophyll_quarter3_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q4, torp_mods_q4,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_chlorophyll_quarter4_diff_model.eps")

# lag/quarter notes
# most effects are probably non significant except some year 1
# chlorophyll change measured in Q1 (Apr-Jun) or Q4 (Jan-Mar) 
# increased with floating plant abundance with 1-year lag
# primary production correlation?
# hydrilla management with 1-year lag also increased
# chlorophyll change measured in Q4


#### start here ####

#### finalize models ####

# SE with heteroskedasticity and autocorrelation
coeftest(hydr_mods[[3]], vcov = vcovHC, type = "HC3") # none sig
coeftest(wahy_mods[[3]], vcov = vcovHC, type = "HC3") # treatment sig
coeftest(wale_mods[[3]], vcov = vcovHC, type = "HC3") # treatment sig
coeftest(cubu_mods[[3]], vcov = vcovHC, type = "HC3") # none sig
coeftest(pagr_mods[[3]], vcov = vcovHC, type = "HC3") # none sig
coeftest(torp_mods[[3]], vcov = vcovHC, type = "HC3") # treatment and PAC sig

# don't need surveyor experience
# only using lag 3
hydr_dat3 <- hydr_dat %>% filter(!is.na(Lag3Treated) & !is.na(Lag3AvgPercCovered)) %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
         AvgPercCovered_c = Lag3AvgPercCovered - mean(Lag3AvgPercCovered)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
wahy_dat3 <- wahy_dat %>% filter(!is.na(Lag3Treated) & !is.na(Lag3AvgPercCovered)) %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
         AvgPercCovered_c = Lag3AvgPercCovered - mean(Lag3AvgPercCovered)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
wale_dat3 <- wale_dat %>% filter(!is.na(Lag3Treated) & !is.na(Lag3AvgPercCovered)) %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
         AvgPercCovered_c = Lag3AvgPercCovered - mean(Lag3AvgPercCovered)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
cubu_dat3 <- cubu_dat %>% filter(!is.na(Lag3Treated) & !is.na(Lag3AvgPercCovered)) %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
         AvgPercCovered_c = Lag3AvgPercCovered - mean(Lag3AvgPercCovered)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
torp_dat3 <- torp_dat %>% filter(!is.na(Lag3Treated) & !is.na(Lag3AvgPercCovered)) %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
         AvgPercCovered_c = Lag3AvgPercCovered - mean(Lag3AvgPercCovered)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
pagr_dat3 <- pagr_dat %>% filter(!is.na(Lag3Treated) & !is.na(Lag3AvgPercCovered)) %>%
  mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
         AvgPercCovered_c = Lag3AvgPercCovered - mean(Lag3AvgPercCovered)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))

# fit models
hydr_nat_rich_mod <- plm(RichnessDiff ~ AvgPercCovered_c + Lag3Treated, data = hydr_dat3,
                         index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
wahy_nat_rich_mod <- update(hydr_nat_rich_mod, data = wahy_dat3)
wale_nat_rich_mod <- update(hydr_nat_rich_mod, data = wale_dat3)
cubu_nat_rich_mod <- update(hydr_nat_rich_mod, data = cubu_dat3)
torp_nat_rich_mod <- update(hydr_nat_rich_mod, data = torp_dat3)
pagr_nat_rich_mod <- update(hydr_nat_rich_mod, data = pagr_dat3)

# check fits
# overall p-values matches management estimate - use adjusted below
summary(hydr_nat_rich_mod)
summary(wahy_nat_rich_mod)
summary(wale_nat_rich_mod)
summary(cubu_nat_rich_mod)
summary(torp_nat_rich_mod)
summary(pagr_nat_rich_mod)
# all are balanced

# SE with heteroscedasticity and autocorrelation
(hydr_nat_rich_mod_se <- coeftest(hydr_nat_rich_mod, vcov = vcovHC, type = "HC3"))
(wahy_nat_rich_mod_se <- coeftest(wahy_nat_rich_mod, vcov = vcovHC, type = "HC3")) # treated sig
(wale_nat_rich_mod_se <- coeftest(wale_nat_rich_mod, vcov = vcovHC, type = "HC3")) # treated sig
(cubu_nat_rich_mod_se <- coeftest(cubu_nat_rich_mod, vcov = vcovHC, type = "HC3")) # none
(pagr_nat_rich_mod_se <- coeftest(pagr_nat_rich_mod, vcov = vcovHC, type = "HC3")) # cover sig
(torp_nat_rich_mod_se <- coeftest(torp_nat_rich_mod, vcov = vcovHC, type = "HC3")) # all sig

# add fitted values to pdata.frame (important to match rows)
# convert to regular dataframe
hydr_fit <- mutate(hydr_dat3, Fitted = as.numeric(hydr_dat3$RichnessDiff - hydr_nat_rich_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit <- mutate(wahy_dat3, Fitted = as.numeric(wahy_dat3$RichnessDiff - wahy_nat_rich_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit <- mutate(wale_dat3, Fitted = as.numeric(wale_dat3$RichnessDiff - wale_nat_rich_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit <- mutate(cubu_dat3, Fitted = as.numeric(cubu_dat3$RichnessDiff - cubu_nat_rich_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
pagr_fit <- mutate(pagr_dat3, Fitted = as.numeric(pagr_dat3$RichnessDiff - pagr_nat_rich_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit <- mutate(torp_dat3, Fitted = as.numeric(torp_dat3$RichnessDiff - torp_nat_rich_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)

# fitted vs. observed
ggplot(hydr_fit, aes(x = Fitted, y = PercDiffCovered)) + geom_point()
ggplot(wahy_fit, aes(x = Fitted, y = PercDiffCovered)) + geom_point()
ggplot(wale_fit, aes(x = Fitted, y = PercDiffCovered)) + geom_point()
ggplot(cubu_fit, aes(x = Fitted, y = PercDiffCovered)) + geom_point()
ggplot(pagr_fit, aes(x = Fitted, y = PercDiffCovered)) + geom_point()
ggplot(torp_fit, aes(x = Fitted, y = PercDiffCovered)) + geom_point()
# generally poor fits

# export models
save(hydr_nat_rich_mod, file = "output/fwc_hydrilla_native_richness_model.rda")
save(wahy_nat_rich_mod, file = "output/fwc_water_hyacinth_native_richness_model.rda")
save(wale_nat_rich_mod, file = "output/fwc_water_lettuce_native_richness_model.rda")
save(cubu_nat_rich_mod, file = "output/fwc_cuban_bulrush_native_richness_model.rda")
save(pagr_nat_rich_mod, file = "output/fwc_para_grass_native_richness_model.rda")
save(torp_nat_rich_mod, file = "output/fwc_torpedograss_native_richness_model.rda")

# combine SE tables
foc_mod_se <- tidy(hydr_nat_rich_mod_se) %>%
  mutate(Invasive = "hydrilla",
         R2 = r.squared(hydr_nat_rich_mod),
         Waterbodies = n_distinct(hydr_dat3$PermanentID),
         Years = n_distinct(hydr_dat3$GSYear),
         N = nrow(hydr_dat3)) %>%
  full_join(tidy(wahy_nat_rich_mod_se) %>%
              mutate(Invasive = "water hyacinth",
                     R2 = r.squared(wahy_nat_rich_mod),
                     Waterbodies = n_distinct(wahy_dat3$PermanentID),
                     Years = n_distinct(wahy_dat3$GSYear),
                     N = nrow(wahy_dat3))) %>%
  full_join(tidy(wale_nat_rich_mod_se) %>%
              mutate(Invasive = "water lettuce",
                     R2 = r.squared(wale_nat_rich_mod),
                     Waterbodies = n_distinct(wale_dat3$PermanentID),
                     Years = n_distinct(wale_dat3$GSYear),
                     N = nrow(wale_dat3))) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "AvgPercCovered_c",
                           "management" = "Lag3Treated")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive)

non_foc_mod_se <- tidy(cubu_nat_rich_mod_se) %>%
  mutate(Invasive = "Cuban bulrush",
         R2 = r.squared(cubu_nat_rich_mod),
         Waterbodies = n_distinct(cubu_dat3$PermanentID),
         Years = n_distinct(cubu_dat3$GSYear),
         N = nrow(cubu_dat3)) %>%
  full_join(tidy(pagr_nat_rich_mod_se) %>%
              mutate(Invasive = "para grass",
                     R2 = r.squared(pagr_nat_rich_mod),
                     Waterbodies = n_distinct(pagr_dat3$PermanentID),
                     Years = n_distinct(pagr_dat3$GSYear),
                     N = nrow(pagr_dat3))) %>%
  full_join(tidy(torp_nat_rich_mod_se) %>%
              mutate(Invasive = "torpedograss",
                     R2 = r.squared(torp_nat_rich_mod),
                     Waterbodies = n_distinct(torp_dat3$PermanentID),
                     Years = n_distinct(torp_dat3$GSYear),
                     N = nrow(torp_dat3))) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "AvgPercCovered_c",
                           "management" = "Lag3Treated")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive)

# export
write_csv(foc_mod_se, "output/fwc_focal_native_richness_model_summary.csv")
write_csv(non_foc_mod_se, "output/fwc_non_focal_native_richness_model_summary.csv")