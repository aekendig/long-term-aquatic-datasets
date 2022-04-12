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
lw_pho <- read_csv("intermediate-data/LW_phosphorus_formatted.csv")
lwwa_pho <- read_csv("intermediate-data/LW_water_atlas_phosphorus_formatted.csv")
uninv <- read_csv("output/fwc_uninvaded_permID.csv") # lakes with no recorded invasion


#### edit data ####

# require longest reasonable lag
# longer lags mean fewer years are needed from phosphorus dataset
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
ggplot(lwwa_pho, aes(x = GSYear)) +
  geom_bar() +
  facet_wrap(~ Quarter)
# some are low for 2019

# combine phosphorus, invasive, control
# select waterbodies sampled throughout
pho_dat <- lwwa_pho %>%
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
(pho_samp_sum <- pho_dat %>%
  group_by(CommonName, Quarter) %>%
  summarize(Years = n_distinct(GSYear),
            Waterbodies = n_distinct(PermanentID),
            minYear = min(GSYear),
            maxYear = max(GSYear)) %>%
  data.frame())

# why are there no data matching with torpedograss?
inv_plant2 %>%
  filter(CommonName == "Torpedograss") %>%
  select(PermanentID, GSYear) %>%
  unique() %>%
  inner_join(lwwa_pho %>%
               select(PermanentID, GSYear) %>%
               unique()) %>%
  group_by(GSYear) %>%
  count() # a lot of waterbodies each year

inv_plant2 %>%
  filter(CommonName == "Torpedograss") %>%
  select(PermanentID, GSYear) %>%
  unique() %>%
  inner_join(lwwa_pho %>%
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
inv_taxa <- sort(unique(pho_dat$CommonName))

# loop through taxa
pdf("output/phosphorus_continuous_time_series_by_taxon.pdf")

for(i in 1:length(inv_taxa)){
  
  # subset data
  subdat <- pho_dat %>% filter(CommonName == inv_taxa[i])
  subdat_ctrl <- subdat %>% filter(Lag1Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = QualityValue, color = PermanentID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          facet_wrap(~ Quarter) + 
          labs(x = "Year", y = "Total phosphorus (ug/L)", title = inv_taxa[i]) +
          def_theme_paper +
          theme(plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()

# look at high values
pho_dat %>%
  filter(QualityValue > 1000) %>%
  select(PermanentID, GSYear, QualityValue) %>%
  unique() %>%
  inner_join(lwwa_pho)
# checked data in phosphorus_data_processing.R and they seem fine

# remove paragrass (only 1-2 waterbodies)
# use same waterbodies in all 4 quarters
pho_dat2 <- pho_dat %>%
  filter(CommonName != "Para grass") %>%
  group_by(CommonName, PermanentID, GSYear) %>%
  mutate(nQuart = n_distinct(Quarter)) %>%
  ungroup() %>%
  filter(nQuart == 4)

# sample sizes
pho_dat2 %>%
  group_by(CommonName, Quarter) %>%
  summarize(Years = n_distinct(GSYear),
            Waterbodies = n_distinct(PermanentID),
            N = Years * Waterbodies) %>%
  data.frame()

# split by species
hydr_dat <- filter(pho_dat2, CommonName == "Hydrilla")
wale_dat <- filter(pho_dat2, CommonName == "Water lettuce")
wahy_dat <- filter(pho_dat2, CommonName == "Water hyacinth")
torp_dat <- filter(pho_dat2, CommonName == "Torpedograss")
cubu_dat <- filter(pho_dat2, CommonName == "Cuban bulrush")

# add water quality to uninvaded dataset
# select years to match invasion dataset
uninv2 <- lwwa_pho %>%
  filter(!is.na(PrevValue)) %>%
  inner_join(uninv) %>%
  left_join(pho_samp_sum %>%
              select(CommonName, minYear, maxYear) %>%
              unique()) %>%
  filter(GSYear >= minYear & GSYear <= maxYear)


#### initial visualizations ####

# covariate correlations
pho_dat2 %>%
  select(Quarter, CommonName, 
         Lag3Treated, Lag3AvgPercCovered, PrevValue) %>%
  group_by(Quarter, CommonName) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & abs(corr) >= 0.4) %>%
  data.frame()
# floating plants: prev value and PAC

# response distributions
ggplot(pho_dat2, aes(x = ValueDiff)) +
  geom_histogram() +
  facet_grid(CommonName ~ Quarter, scales = "free")
# normal and highly clustered around zero

ggplot(pho_dat2, aes(x = QualityValue)) +
  geom_histogram() +
  facet_grid(CommonName ~ Quarter, scales = "free")
# skewed

# coefficients and QualityValue
ggplot(pho_dat2, aes(x = Lag3AvgPercCovered, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# generally close to zero

ggplot(pho_dat2, aes(x = Lag3AvgPercCovered, y = QualityValue)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# slight negative for hydrilla
# strong positive for floating plants

ggplot(pho_dat2, aes(x = Lag3Treated, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# negligible

ggplot(pho_dat2, aes(x = Lag1Treated, y = QualityValue)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# positive for hydrilla (maybe correlated with hydrilla abundance)
# others seem close to zero

ggplot(pho_dat2, aes(x = PrevValue, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(CommonName ~ Quarter, scales = "free")
# consistently negative

ggplot(pho_dat2, aes(x = PrevValue, y = QualityValue)) +
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
wahy_mod_struc <- mod_structure_fits(wahy_dat) # model convergence error for glmmTMB
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

write_csv(mod_comp, "output/fwc_phosphorus_model_structure_comparison.csv")

# model comparison notes:
# simple model: hydrilla PAC and management decreased
# floating PAC and management increased
# random location: all PAC increased, hydrilla management decreased, floating management increased
# random year: same as simple model
# random location + year: similar to random location
# fixed location: hydrilla PAC no effect, floating PAC decreased
# hydrilla management decreased, floating management increased
# fixed lcoation + year: similar to fixed location
# initial value: similar to random location
# difference response: hydrilla PAC no effect, floating PAC increased
# hydrilla management decreased, floating management increased
# difference response + initial: similar to random location

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
plmtest(wahy_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # not sig

wale_mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag1AvgPercCovered + Lag1Treated, 
                                data = filter(wale_dat, Quarter == 3), 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(wale_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(wale_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # not sig

# use annual difference without initial value and with waterbody and year fixed effects


#### model-fitting functions ####

# data filter function
dat_mod_filt <- function(treat_col, inv_col, dat_in){
  
  dat_mod <- dat_in %>%
    mutate(Treated = !!sym(treat_col),
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
    
    mod1 <- plm(ValueDiff ~ AvgPercCovered_c + Treated, data = dat_mod1,
                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
    mod2 <- update(mod1, data = dat_mod2)
    mod3 <- update(mod1, data = dat_mod3)
    
    mods_out <- list(mod1, mod2, mod3)
    
  } else {
    
    mod1 <- plm(ValueDiff ~ AvgPercCovered_c + Treated, data = dat_mod1,
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
coef_names <- c("Treated" = "Management", 
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
    title = "Effects on annual difference in total phosphorus")
  
  ggsave(filename, comb_fig,
         device = "eps", width = 6.5, height = 3, units = "in")
  
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
      title = "Effects on annual difference in total phosphorus")
  
  ggsave(filename, comb_fig,
         device = "eps", width = 4.7, height = 3, units = "in")
  
}

# figures
panel_plot_fun(hydr_mods_q1, wahy_mods_q1, wale_mods_q1,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_phosphorus_quarter1_diff_model.eps")
panel_plot_fun(hydr_mods_q2, wahy_mods_q2, wale_mods_q2,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_phosphorus_quarter2_diff_model.eps")
panel_plot_fun(hydr_mods_q3, wahy_mods_q3, wale_mods_q3,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_phosphorus_quarter3_diff_model.eps")
panel_plot_fun(hydr_mods_q4, wahy_mods_q4, wale_mods_q4,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_phosphorus_quarter4_diff_model.eps")

panel_plot_non_foc_fun(cubu_mods_q1, torp_mods_q1,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_phosphorus_quarter1_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q2, torp_mods_q2,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_phosphorus_quarter2_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q3, torp_mods_q3,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_phosphorus_quarter3_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q4, torp_mods_q4,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_phosphorus_quarter4_diff_model.eps")

# lag/quarter notes
# generally positive PAC effects and management effects around zero
# lag 3 generally represents conservative trends


#### finalize models ####

# data format function
dat_mod_fin <- function(dat_in, quarter){
  
  # center invasive PAC
  # make p data frame
  dat_out <- dat_in %>%
    filter(Quarter == quarter) %>%
    mutate(AvgPercCovered_c = Lag3AvgPercCovered - mean(Lag3AvgPercCovered),
           Treated = Lag3Treated) %>%
    pdata.frame(index = c("PermanentID", "GSYear"))
  
  return(dat_out)
    
}

# format data
hydr_dat3_q1 <- dat_mod_fin(hydr_dat, 1)
wahy_dat3_q1 <- dat_mod_fin(wahy_dat, 1)
wale_dat3_q1 <- dat_mod_fin(wale_dat, 1)
cubu_dat3_q1 <- dat_mod_fin(cubu_dat, 1)
torp_dat3_q1 <- dat_mod_fin(torp_dat, 1)

hydr_dat3_q2 <- dat_mod_fin(hydr_dat, 2)
wahy_dat3_q2 <- dat_mod_fin(wahy_dat, 2)
wale_dat3_q2 <- dat_mod_fin(wale_dat, 2)
cubu_dat3_q2 <- dat_mod_fin(cubu_dat, 2)
torp_dat3_q2 <- dat_mod_fin(torp_dat, 2)

hydr_dat3_q3 <- dat_mod_fin(hydr_dat, 3)
wahy_dat3_q3 <- dat_mod_fin(wahy_dat, 3)
wale_dat3_q3 <- dat_mod_fin(wale_dat, 3)
cubu_dat3_q3 <- dat_mod_fin(cubu_dat, 3)
torp_dat3_q3 <- dat_mod_fin(torp_dat, 3)

hydr_dat3_q4 <- dat_mod_fin(hydr_dat, 4)
wahy_dat3_q4 <- dat_mod_fin(wahy_dat, 4)
wale_dat3_q4 <- dat_mod_fin(wale_dat, 4)
cubu_dat3_q4 <- dat_mod_fin(cubu_dat, 4)
torp_dat3_q4 <- dat_mod_fin(torp_dat, 4)

# fit models
hydr_pho_mod_q1 <- plm(ValueDiff ~ AvgPercCovered_c + Treated, data = hydr_dat3_q1,
                       index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
wahy_pho_mod_q1 <- update(hydr_pho_mod_q1, data = wahy_dat3_q1)
wale_pho_mod_q1 <- update(hydr_pho_mod_q1, data = wale_dat3_q1)
cubu_pho_mod_q1 <- update(hydr_pho_mod_q1, data = cubu_dat3_q1)
torp_pho_mod_q1 <- update(hydr_pho_mod_q1, data = torp_dat3_q1)

hydr_pho_mod_q2 <- update(hydr_pho_mod_q1, data = hydr_dat3_q2)
wahy_pho_mod_q2 <- update(hydr_pho_mod_q2, data = wahy_dat3_q2)
wale_pho_mod_q2 <- update(hydr_pho_mod_q2, data = wale_dat3_q2)
cubu_pho_mod_q2 <- update(hydr_pho_mod_q2, data = cubu_dat3_q2)
torp_pho_mod_q2 <- update(hydr_pho_mod_q2, data = torp_dat3_q2)

hydr_pho_mod_q3 <- update(hydr_pho_mod_q1, data = hydr_dat3_q3)
wahy_pho_mod_q3 <- update(hydr_pho_mod_q3, data = wahy_dat3_q3)
wale_pho_mod_q3 <- update(hydr_pho_mod_q3, data = wale_dat3_q3)
cubu_pho_mod_q3 <- update(hydr_pho_mod_q3, data = cubu_dat3_q3)
torp_pho_mod_q3 <- update(hydr_pho_mod_q3, data = torp_dat3_q3)

hydr_pho_mod_q4 <- update(hydr_pho_mod_q1, data = hydr_dat3_q4)
wahy_pho_mod_q4 <- update(hydr_pho_mod_q4, data = wahy_dat3_q4)
wale_pho_mod_q4 <- update(hydr_pho_mod_q4, data = wale_dat3_q4)
cubu_pho_mod_q4 <- update(hydr_pho_mod_q4, data = cubu_dat3_q4)
torp_pho_mod_q4 <- update(hydr_pho_mod_q4, data = torp_dat3_q4)

# SE with heteroscedasticity and autocorrelation
coeftest(hydr_pho_mod_q1, vcov = vcovHC, type = "HC3")
coeftest(wahy_pho_mod_q1, vcov = vcovHC, type = "HC3") # +PAC
coeftest(wale_pho_mod_q1, vcov = vcovHC, type = "HC3") # +PAC 
coeftest(cubu_pho_mod_q1, vcov = vcovHC, type = "HC3") 
coeftest(torp_pho_mod_q1, vcov = vcovHC, type = "HC3") # +PAC marginal 

coeftest(hydr_pho_mod_q2, vcov = vcovHC, type = "HC3") # +treat marginal
coeftest(wahy_pho_mod_q2, vcov = vcovHC, type = "HC3") # +PAC 
coeftest(wale_pho_mod_q2, vcov = vcovHC, type = "HC3") 
coeftest(cubu_pho_mod_q2, vcov = vcovHC, type = "HC3")
coeftest(torp_pho_mod_q2, vcov = vcovHC, type = "HC3") # +treat 

coeftest(hydr_pho_mod_q3, vcov = vcovHC, type = "HC3")
coeftest(wahy_pho_mod_q3, vcov = vcovHC, type = "HC3") # +PAC 
coeftest(wale_pho_mod_q3, vcov = vcovHC, type = "HC3") # +PAC +treat, both marginal 
coeftest(cubu_pho_mod_q3, vcov = vcovHC, type = "HC3") # +treat 
coeftest(torp_pho_mod_q3, vcov = vcovHC, type = "HC3") 

coeftest(hydr_pho_mod_q4, vcov = vcovHC, type = "HC3")
coeftest(wahy_pho_mod_q4, vcov = vcovHC, type = "HC3") # +PAC 
coeftest(wale_pho_mod_q4, vcov = vcovHC, type = "HC3") # +PAC 
coeftest(cubu_pho_mod_q4, vcov = vcovHC, type = "HC3") 
coeftest(torp_pho_mod_q4, vcov = vcovHC, type = "HC3") # +treat marginal 

# add fitted values to pdata.frame (important to match rows)
# convert to regular dataframe
hydr_fit_q1 <- mutate(hydr_dat3_q1, Fitted = as.numeric(hydr_dat3_q1$ValueDiff - hydr_pho_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q1 <- mutate(wahy_dat3_q1, Fitted = as.numeric(wahy_dat3_q1$ValueDiff - wahy_pho_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q1 <- mutate(wale_dat3_q1, Fitted = as.numeric(wale_dat3_q1$ValueDiff - wale_pho_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q1 <- mutate(cubu_dat3_q1, Fitted = as.numeric(cubu_dat3_q1$ValueDiff - cubu_pho_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q1 <- mutate(torp_dat3_q1, Fitted = as.numeric(torp_dat3_q1$ValueDiff - torp_pho_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_fit_q2 <- mutate(hydr_dat3_q2, Fitted = as.numeric(hydr_dat3_q2$ValueDiff - hydr_pho_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q2 <- mutate(wahy_dat3_q2, Fitted = as.numeric(wahy_dat3_q2$ValueDiff - wahy_pho_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q2 <- mutate(wale_dat3_q2, Fitted = as.numeric(wale_dat3_q2$ValueDiff - wale_pho_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q2 <- mutate(cubu_dat3_q2, Fitted = as.numeric(cubu_dat3_q2$ValueDiff - cubu_pho_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q2 <- mutate(torp_dat3_q2, Fitted = as.numeric(torp_dat3_q2$ValueDiff - torp_pho_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_fit_q3 <- mutate(hydr_dat3_q3, Fitted = as.numeric(hydr_dat3_q3$ValueDiff - hydr_pho_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q3 <- mutate(wahy_dat3_q3, Fitted = as.numeric(wahy_dat3_q3$ValueDiff - wahy_pho_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q3 <- mutate(wale_dat3_q3, Fitted = as.numeric(wale_dat3_q3$ValueDiff - wale_pho_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q3 <- mutate(cubu_dat3_q3, Fitted = as.numeric(cubu_dat3_q3$ValueDiff - cubu_pho_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q3 <- mutate(torp_dat3_q3, Fitted = as.numeric(torp_dat3_q3$ValueDiff - torp_pho_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_fit_q4 <- mutate(hydr_dat3_q4, Fitted = as.numeric(hydr_dat3_q4$ValueDiff - hydr_pho_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q4 <- mutate(wahy_dat3_q4, Fitted = as.numeric(wahy_dat3_q4$ValueDiff - wahy_pho_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q4 <- mutate(wale_dat3_q4, Fitted = as.numeric(wale_dat3_q4$ValueDiff - wale_pho_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q4 <- mutate(cubu_dat3_q4, Fitted = as.numeric(cubu_dat3_q4$ValueDiff - cubu_pho_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q4 <- mutate(torp_dat3_q4, Fitted = as.numeric(torp_dat3_q4$ValueDiff - torp_pho_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)

# fitted vs. observed
ggplot(hydr_fit_q1, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(wahy_fit_q1, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(wale_fit_q1, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(cubu_fit_q1, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(torp_fit_q1, aes(x = Fitted, y = ValueDiff)) + geom_point()

ggplot(hydr_fit_q2, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(wahy_fit_q2, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(wale_fit_q2, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(cubu_fit_q2, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(torp_fit_q2, aes(x = Fitted, y = ValueDiff)) + geom_point()

ggplot(hydr_fit_q3, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(wahy_fit_q3, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(wale_fit_q3, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(cubu_fit_q3, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(torp_fit_q3, aes(x = Fitted, y = ValueDiff)) + geom_point()

ggplot(hydr_fit_q4, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(wahy_fit_q4, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(wale_fit_q4, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(cubu_fit_q4, aes(x = Fitted, y = ValueDiff)) + geom_point()
ggplot(torp_fit_q4, aes(x = Fitted, y = ValueDiff)) + geom_point()
# look okay

# combine models
hydr_pho_mods <- list(hydr_pho_mod_q1, hydr_pho_mod_q2, hydr_pho_mod_q3, hydr_pho_mod_q4)
wahy_pho_mods <- list(wahy_pho_mod_q1, wahy_pho_mod_q2, wahy_pho_mod_q3, wahy_pho_mod_q4)
wale_pho_mods <- list(wale_pho_mod_q1, wale_pho_mod_q2, wale_pho_mod_q3, wale_pho_mod_q4)
cubu_pho_mods <- list(cubu_pho_mod_q1, cubu_pho_mod_q2, cubu_pho_mod_q3, cubu_pho_mod_q4)
torp_pho_mods <- list(torp_pho_mod_q1, torp_pho_mod_q2, torp_pho_mod_q3, torp_pho_mod_q4)

# export models
save(hydr_pho_mods, file = "output/fwc_hydrilla_phosphorus_models.rda")
save(wahy_pho_mods, file = "output/fwc_water_hyacinth_phosphorus_models.rda")
save(wale_pho_mods, file = "output/fwc_water_lettuce_phosphorus_models.rda")
save(cubu_pho_mods, file = "output/fwc_cuban_bulrush_phosphorus_models.rda")
save(torp_pho_mods, file = "output/fwc_torpedograss_phosphorus_models.rda")

# load models
load("output/fwc_hydrilla_phosphorus_models.rda")
load("output/fwc_water_hyacinth_phosphorus_models.rda")
load("output/fwc_water_lettuce_phosphorus_models.rda")
load("output/fwc_cuban_bulrush_phosphorus_models.rda")
load("output/fwc_torpedograss_phosphorus_models.rda")

# process model SE
mod_se_fun <- function(models, dat, spp){
  
  dat_out <- tidy(coeftest(models[[1]], vcov = vcovHC, type = "HC3")) %>%
    mutate(Quarter = "Apr-Jun",
           R2 = r.squared(models[[1]])) %>%
    full_join(tidy(coeftest(models[[2]], vcov = vcovHC, type = "HC3")) %>%
                mutate(Quarter = "Jul-Sep",
                       R2 = r.squared(models[[2]]))) %>%
    full_join(tidy(coeftest(models[[3]], vcov = vcovHC, type = "HC3")) %>%
                mutate(Quarter = "Oct-Dec",
                       R2 = r.squared(models[[3]]))) %>%
    full_join(tidy(coeftest(models[[4]], vcov = vcovHC, type = "HC3")) %>%
                mutate(Quarter = "Jan-Mar",
                       R2 = r.squared(models[[4]]))) %>%
    mutate(Invasive = spp,
           Waterbodies = n_distinct(dat$PermanentID),
           Years = n_distinct(dat$GSYear),
           N = Waterbodies * Years)
  
  return(dat_out)
  
}

# combine SE tables
foc_mod_se <- mod_se_fun(hydr_pho_mods, hydr_dat, "hydrilla") %>%
  full_join(mod_se_fun(wahy_pho_mods, wahy_dat, "water hyacinth")) %>%
  full_join(mod_se_fun(wale_pho_mods, wale_dat, "water lettuce")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "AvgPercCovered_c",
                           "management" = "Treated")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Quarter)

non_foc_mod_se <- mod_se_fun(cubu_pho_mods, cubu_dat, "Cuban bulrush") %>%
  full_join(mod_se_fun(torp_pho_mods, torp_dat, "torpedograss")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "AvgPercCovered_c",
                           "management" = "Treated")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Quarter)

# export
write_csv(foc_mod_se, "output/fwc_focal_phosphorus_model_summary.csv")
write_csv(non_foc_mod_se, "output/fwc_non_focal_phosphorus_model_summary.csv")


#### values for text ####

# summarize uninvaded
uninv_sum <- uninv2 %>%
  group_by(CommonName, Quarter) %>%
  summarize(UninvAvg = mean(PrevValue),
            UninvN = n()) %>%
  ungroup()

# translate model coefficients
mod_coef_fun <- function(models, spp){
  
  dat_out <- tibble(Invasive = spp,
                    Quarter = c("Apr-Jun", "Jul-Sep", "Oct-Dec", "Jan-Mar"),
                    DiffAvg = sapply(models, function(x) mean(fixef(x))),
                    PACEffect = sapply(models, function(x) coef(x)[1]),
                    TreatEffect = DiffAvg + sapply(models, function(x) coef(x)[2]))
  
  return(dat_out)
  
}

# identify significant effects
foc_sig <- foc_mod_se %>%
  filter(P < 0.1) %>%
  select(Invasive, Quarter, Term) %>%
  left_join(mod_coef_fun(hydr_pho_mods, "hydrilla") %>%
              full_join(mod_coef_fun(wahy_pho_mods, "water hyacinth")) %>%
              full_join(mod_coef_fun(wale_pho_mods, "water lettuce"))) %>%
  mutate(PACEffect = if_else(Term == "management", NA_real_, PACEffect),
         TreatEffect = if_else(Term == "invasive PAC", NA_real_, TreatEffect),
         Metric = "total phosphorus") %>%
  left_join(hydr_dat %>%
              group_by(CommonName, Quarter) %>%
              summarize(Average = mean(PrevValue)) %>%
              ungroup() %>%
              full_join(wahy_dat %>%
                          group_by(CommonName, Quarter) %>%
                          summarize(Average = mean(PrevValue)) %>%
                          ungroup()) %>%
              full_join(wale_dat %>%
                          group_by(CommonName, Quarter) %>%
                          summarize(Average = mean(PrevValue)) %>%
                          ungroup()) %>%
              left_join(uninv_sum) %>%
              mutate(Quarter = case_when(Quarter == 1 ~ "Apr-Jun", 
                                         Quarter == 2 ~ "Jul-Sep", 
                                         Quarter == 3 ~ "Oct-Dec", 
                                         Quarter == 4 ~ "Jan-Mar"),
                     CommonName = tolower(CommonName)) %>%
              rename(Invasive = CommonName))

write_csv(foc_sig, "output/fwc_focal_invasive_phosphorus_significant.csv")

non_foc_sig <- non_foc_mod_se %>%
  filter(P < 0.1) %>%
  select(Invasive, Quarter, Term) %>%
  left_join(mod_coef_fun(cubu_pho_mods, "Cuban bulrush") %>%
              full_join(mod_coef_fun(torp_pho_mods, "torpedograss"))) %>%
  mutate(PACEffect = if_else(Term == "management", NA_real_, PACEffect),
         TreatEffect = if_else(Term == "invasive PAC", NA_real_, TreatEffect),
         Metric = "total phosphorus") %>%
  left_join(cubu_dat %>%
              group_by(CommonName, Quarter) %>%
              summarize(Average = mean(PrevValue)) %>%
              ungroup() %>%
              full_join(torp_dat %>%
                          group_by(CommonName, Quarter) %>%
                          summarize(Average = mean(PrevValue)) %>%
                          ungroup()) %>%
              left_join(uninv_sum) %>%
              mutate(Quarter = case_when(Quarter == 1 ~ "Apr-Jun", 
                                         Quarter == 2 ~ "Jul-Sep", 
                                         Quarter == 3 ~ "Oct-Dec", 
                                         Quarter == 4 ~ "Jan-Mar"),
                     CommonName = fct_recode(CommonName, 
                                             "torpedograss" = "Torpedograss")) %>%
              rename(Invasive = CommonName))
  
write_csv(non_foc_sig, "output/fwc_non_focal_invasive_phosphorus_significant.csv")
