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
lw_sec <- read_csv("intermediate-data/LW_secchi_formatted.csv")
lwwa_sec <- read_csv("intermediate-data/LW_water_atlas_secchi_formatted.csv")
uninv <- read_csv("output/fwc_uninvaded_permID.csv") # lakes with no recorded invasion


#### edit data ####

# require longest reasonable lag
# longer lags mean fewer years are needed from secchi dataset
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
ggplot(lwwa_sec, aes(x = GSYear)) +
  geom_bar() +
  facet_wrap(~ Quarter)
# some are low for 2019

# combine secchi, invasive, control
# select waterbodies sampled throughout
sec_dat <- lwwa_sec %>%
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
(sec_samp_sum <- sec_dat %>%
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
  inner_join(lwwa_sec %>%
               select(PermanentID, GSYear) %>%
               unique()) %>%
  group_by(GSYear) %>%
  count() # a lot of waterbodies each year

inv_plant2 %>%
  filter(CommonName == "Torpedograss") %>%
  select(PermanentID, GSYear) %>%
  unique() %>%
  inner_join(lwwa_sec %>%
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
inv_taxa <- sort(unique(sec_dat$CommonName))

# loop through taxa
pdf("output/secchi_continuous_time_series_by_taxon.pdf")

for(i in 1:length(inv_taxa)){
  
  # subset data
  subdat <- sec_dat %>% filter(CommonName == inv_taxa[i])
  subdat_ctrl <- subdat %>% filter(Lag1Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = QualityValue, color = PermanentID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          facet_wrap(~ Quarter) + 
          labs(x = "Year", y = "secchi a (ug/L)", title = inv_taxa[i]) +
          def_theme_paper +
          theme(plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()

# remove paragrass (only 1-2 waterbodies)
# use same waterbodies in all 4 quarters
sec_dat2 <- sec_dat %>%
  filter(CommonName != "Para grass") %>%
  group_by(CommonName, PermanentID, GSYear) %>%
  mutate(nQuart = n_distinct(Quarter)) %>%
  ungroup() %>%
  filter(nQuart == 4)

# sample sizes
sec_dat2 %>%
  group_by(CommonName, Quarter) %>%
  summarize(Years = n_distinct(GSYear),
            Waterbodies = n_distinct(PermanentID),
            N = Years * Waterbodies) %>%
  data.frame()

# split by species
hydr_dat <- filter(sec_dat2, CommonName == "Hydrilla")
wale_dat <- filter(sec_dat2, CommonName == "Water lettuce")
wahy_dat <- filter(sec_dat2, CommonName == "Water hyacinth")
torp_dat <- filter(sec_dat2, CommonName == "Torpedograss")
cubu_dat <- filter(sec_dat2, CommonName == "Cuban bulrush")

# add water quality to uninvaded dataset
# select years to match invasion dataset
uninv2 <- lwwa_sec %>%
  filter(!is.na(PrevValue)) %>%
  inner_join(uninv) %>%
  left_join(sec_samp_sum %>%
              select(CommonName, minYear, maxYear) %>%
              unique()) %>%
  filter(GSYear >= minYear & GSYear <= maxYear)


#### initial visualizations ####

# covariate correlations
sec_dat2 %>%
  select(Quarter, CommonName, 
         Lag3Treated, Lag3AvgPercCovered, PrevValue) %>%
  group_by(Quarter, CommonName) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & abs(corr) >= 0.4) %>%
  data.frame()
# torpedograss: previous value and PAC are highly correlated

# response distributions
ggplot(sec_dat2, aes(x = ValueDiff)) +
  geom_histogram() +
  facet_grid(CommonName ~ Quarter, scales = "free")
# normal and relatively wide around zero

ggplot(sec_dat2, aes(x = QualityValue)) +
  geom_histogram() +
  facet_grid(CommonName ~ Quarter, scales = "free")
# skewed

# coefficients and QualityValue
ggplot(sec_dat2, aes(x = Lag3AvgPercCovered, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# generally close to zero

ggplot(sec_dat2, aes(x = Lag3AvgPercCovered, y = QualityValue)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# strong positive for torpedograss
# negative for water hyacinth and lettuce

ggplot(sec_dat2, aes(x = Lag3Treated, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# negligible

ggplot(sec_dat2, aes(x = Lag1Treated, y = QualityValue)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Quarter ~ CommonName, scales = "free")
# minimal changes

ggplot(sec_dat2, aes(x = PrevValue, y = ValueDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(CommonName ~ Quarter, scales = "free")
# consistently negative

ggplot(sec_dat2, aes(x = PrevValue, y = QualityValue)) +
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

write_csv(mod_comp, "output/fwc_secchi_model_structure_comparison.csv")

# model comparison notes:
# simple model: PAC increases secchi (hydrilla) or reduces it (floating),
# management reduces secchi (all three)
# location random effect: hydrilla PAC now reduces secchi, water hyacinth increases, 
# and management increases for those two
# floating PAC and management increase it
# year random effect: similar to simple model
# location + year random: similar to location random
# location fixed effect: similar to location random
# location + year fixed: similar to location random
# initial quality: hydrilla and water lettuce PAC and management decrease secchi
# water hyacinth PAC and management increase
# difference response: similar to simple model
# difference response + initial quality: estimates almost identical to initial quality

# test fixed effects (seems like year isn't necessary)
# have to refit because data need to be accessible (not "dat_fix")
hydr_mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag1AvgPercCovered + Lag1Treated, 
                                data = filter(hydr_dat, Quarter == 3), 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(hydr_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(hydr_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

wahy_mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag1AvgPercCovered + Lag1Treated, 
                                data = filter(wahy_dat, Quarter == 3), 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(wahy_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(wahy_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

wale_mod_diff_fix_loc_yr <- plm(ValueDiff ~ Lag1AvgPercCovered + Lag1Treated, 
                                data = filter(wale_dat, Quarter == 3), 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(wale_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(wale_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

# use annual difference without initial PAC and year fixed effects


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
                index = c("PermanentID", "GSYear"), model = "within")
    mod2 <- update(mod1, data = dat_mod2)
    mod3 <- update(mod1, data = dat_mod3)
    
    mods_out <- list(mod1, mod2, mod3)
    
  } else {
    
    mod1 <- plm(ValueDiff ~ AvgPercCovered_c + Treated, data = dat_mod1,
                index = c("PermanentID", "GSYear"), model = "within")
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
    title = "Effects on annual difference in secchi depth (ft)")
  
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
      title = "Effects on annual difference in secchi depth (ft)")
  
  ggsave(filename, comb_fig,
         device = "eps", width = 4.7, height = 3, units = "in")
  
}

# figures
panel_plot_fun(hydr_mods_q1, wahy_mods_q1, wale_mods_q1,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_secchi_quarter1_diff_model.eps")
panel_plot_fun(hydr_mods_q2, wahy_mods_q2, wale_mods_q2,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_secchi_quarter2_diff_model.eps")
panel_plot_fun(hydr_mods_q3, wahy_mods_q3, wale_mods_q3,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_secchi_quarter3_diff_model.eps")
panel_plot_fun(hydr_mods_q4, wahy_mods_q4, wale_mods_q4,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_secchi_quarter4_diff_model.eps")

panel_plot_non_foc_fun(cubu_mods_q1, torp_mods_q1,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_secchi_quarter1_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q2, torp_mods_q2,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_secchi_quarter2_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q3, torp_mods_q3,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_secchi_quarter3_diff_model.eps")
panel_plot_non_foc_fun(cubu_mods_q4, torp_mods_q4,
                       "Cuban bulrush", "Torpedograss",
                       "output/fwc_non_focal_secchi_quarter4_diff_model.eps")

# lag/quarter notes
# very small PAC effects
# management becomes more negative with time, but not sig
# negative torpedograss effect (something weird about torpedograss)
# negative management for CB, positive for torpedograss


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
hydr_sec_mod_q1 <- plm(ValueDiff ~ AvgPercCovered_c + Treated, data = hydr_dat3_q1,
                       index = c("PermanentID", "GSYear"), model = "within")
wahy_sec_mod_q1 <- update(hydr_sec_mod_q1, data = wahy_dat3_q1)
wale_sec_mod_q1 <- update(hydr_sec_mod_q1, data = wale_dat3_q1)
cubu_sec_mod_q1 <- update(hydr_sec_mod_q1, data = cubu_dat3_q1)
torp_sec_mod_q1 <- update(hydr_sec_mod_q1, data = torp_dat3_q1)

hydr_sec_mod_q2 <- update(hydr_sec_mod_q1, data = hydr_dat3_q2)
wahy_sec_mod_q2 <- update(hydr_sec_mod_q2, data = wahy_dat3_q2)
wale_sec_mod_q2 <- update(hydr_sec_mod_q2, data = wale_dat3_q2)
cubu_sec_mod_q2 <- update(hydr_sec_mod_q2, data = cubu_dat3_q2)
torp_sec_mod_q2 <- update(hydr_sec_mod_q2, data = torp_dat3_q2)

hydr_sec_mod_q3 <- update(hydr_sec_mod_q1, data = hydr_dat3_q3)
wahy_sec_mod_q3 <- update(hydr_sec_mod_q3, data = wahy_dat3_q3)
wale_sec_mod_q3 <- update(hydr_sec_mod_q3, data = wale_dat3_q3)
cubu_sec_mod_q3 <- update(hydr_sec_mod_q3, data = cubu_dat3_q3)
torp_sec_mod_q3 <- update(hydr_sec_mod_q3, data = torp_dat3_q3)

hydr_sec_mod_q4 <- update(hydr_sec_mod_q1, data = hydr_dat3_q4)
wahy_sec_mod_q4 <- update(hydr_sec_mod_q4, data = wahy_dat3_q4)
wale_sec_mod_q4 <- update(hydr_sec_mod_q4, data = wale_dat3_q4)
cubu_sec_mod_q4 <- update(hydr_sec_mod_q4, data = cubu_dat3_q4)
torp_sec_mod_q4 <- update(hydr_sec_mod_q4, data = torp_dat3_q4)

# SE with heteroscedasticity and autocorrelation
coeftest(hydr_sec_mod_q1, vcov = vcovHC, type = "HC3") # PAC
coeftest(wahy_sec_mod_q1, vcov = vcovHC, type = "HC3") # PAC marg
coeftest(wale_sec_mod_q1, vcov = vcovHC, type = "HC3") # PAC
coeftest(cubu_sec_mod_q1, vcov = vcovHC, type = "HC3") 
coeftest(torp_sec_mod_q1, vcov = vcovHC, type = "HC3") # PAC

coeftest(hydr_sec_mod_q2, vcov = vcovHC, type = "HC3")
coeftest(wahy_sec_mod_q2, vcov = vcovHC, type = "HC3") 
coeftest(wale_sec_mod_q2, vcov = vcovHC, type = "HC3") # PAC
coeftest(cubu_sec_mod_q2, vcov = vcovHC, type = "HC3")
coeftest(torp_sec_mod_q2, vcov = vcovHC, type = "HC3") # PAC

coeftest(hydr_sec_mod_q3, vcov = vcovHC, type = "HC3")
coeftest(wahy_sec_mod_q3, vcov = vcovHC, type = "HC3") 
coeftest(wale_sec_mod_q3, vcov = vcovHC, type = "HC3") 
coeftest(cubu_sec_mod_q3, vcov = vcovHC, type = "HC3") 
coeftest(torp_sec_mod_q3, vcov = vcovHC, type = "HC3") # PAC + mgmt

coeftest(hydr_sec_mod_q4, vcov = vcovHC, type = "HC3")
coeftest(wahy_sec_mod_q4, vcov = vcovHC, type = "HC3") # mgmt marg
coeftest(wale_sec_mod_q4, vcov = vcovHC, type = "HC3") 
coeftest(cubu_sec_mod_q4, vcov = vcovHC, type = "HC3") 
coeftest(torp_sec_mod_q4, vcov = vcovHC, type = "HC3") # PAC

# add fitted values to pdata.frame (important to match rows)
# convert to regular dataframe
hydr_fit_q1 <- mutate(hydr_dat3_q1, Fitted = as.numeric(hydr_dat3_q1$ValueDiff - hydr_sec_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q1 <- mutate(wahy_dat3_q1, Fitted = as.numeric(wahy_dat3_q1$ValueDiff - wahy_sec_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q1 <- mutate(wale_dat3_q1, Fitted = as.numeric(wale_dat3_q1$ValueDiff - wale_sec_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q1 <- mutate(cubu_dat3_q1, Fitted = as.numeric(cubu_dat3_q1$ValueDiff - cubu_sec_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q1 <- mutate(torp_dat3_q1, Fitted = as.numeric(torp_dat3_q1$ValueDiff - torp_sec_mod_q1$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_fit_q2 <- mutate(hydr_dat3_q2, Fitted = as.numeric(hydr_dat3_q2$ValueDiff - hydr_sec_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q2 <- mutate(wahy_dat3_q2, Fitted = as.numeric(wahy_dat3_q2$ValueDiff - wahy_sec_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q2 <- mutate(wale_dat3_q2, Fitted = as.numeric(wale_dat3_q2$ValueDiff - wale_sec_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q2 <- mutate(cubu_dat3_q2, Fitted = as.numeric(cubu_dat3_q2$ValueDiff - cubu_sec_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q2 <- mutate(torp_dat3_q2, Fitted = as.numeric(torp_dat3_q2$ValueDiff - torp_sec_mod_q2$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_fit_q3 <- mutate(hydr_dat3_q3, Fitted = as.numeric(hydr_dat3_q3$ValueDiff - hydr_sec_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q3 <- mutate(wahy_dat3_q3, Fitted = as.numeric(wahy_dat3_q3$ValueDiff - wahy_sec_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q3 <- mutate(wale_dat3_q3, Fitted = as.numeric(wale_dat3_q3$ValueDiff - wale_sec_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q3 <- mutate(cubu_dat3_q3, Fitted = as.numeric(cubu_dat3_q3$ValueDiff - cubu_sec_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q3 <- mutate(torp_dat3_q3, Fitted = as.numeric(torp_dat3_q3$ValueDiff - torp_sec_mod_q3$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_fit_q4 <- mutate(hydr_dat3_q4, Fitted = as.numeric(hydr_dat3_q4$ValueDiff - hydr_sec_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
wahy_fit_q4 <- mutate(wahy_dat3_q4, Fitted = as.numeric(wahy_dat3_q4$ValueDiff - wahy_sec_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
wale_fit_q4 <- mutate(wale_dat3_q4, Fitted = as.numeric(wale_dat3_q4$ValueDiff - wale_sec_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_fit_q4 <- mutate(cubu_dat3_q4, Fitted = as.numeric(cubu_dat3_q4$ValueDiff - cubu_sec_mod_q4$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_fit_q4 <- mutate(torp_dat3_q4, Fitted = as.numeric(torp_dat3_q4$ValueDiff - torp_sec_mod_q4$residuals)) %>%
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

# combine models
hydr_sec_mods <- list(hydr_sec_mod_q1, hydr_sec_mod_q2, hydr_sec_mod_q3, hydr_sec_mod_q4)
wahy_sec_mods <- list(wahy_sec_mod_q1, wahy_sec_mod_q2, wahy_sec_mod_q3, wahy_sec_mod_q4)
wale_sec_mods <- list(wale_sec_mod_q1, wale_sec_mod_q2, wale_sec_mod_q3, wale_sec_mod_q4)
cubu_sec_mods <- list(cubu_sec_mod_q1, cubu_sec_mod_q2, cubu_sec_mod_q3, cubu_sec_mod_q4)
torp_sec_mods <- list(torp_sec_mod_q1, torp_sec_mod_q2, torp_sec_mod_q3, torp_sec_mod_q4)

# export models
save(hydr_sec_mods, file = "output/fwc_hydrilla_secchi_models.rda")
save(wahy_sec_mods, file = "output/fwc_water_hyacinth_secchi_models.rda")
save(wale_sec_mods, file = "output/fwc_water_lettuce_secchi_models.rda")
save(cubu_sec_mods, file = "output/fwc_cuban_bulrush_secchi_models.rda")
save(torp_sec_mods, file = "output/fwc_torpedograss_secchi_models.rda")

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
foc_mod_se <- mod_se_fun(hydr_sec_mods, hydr_dat, "hydrilla") %>%
  full_join(mod_se_fun(wahy_sec_mods, wahy_dat, "water hyacinth")) %>%
  full_join(mod_se_fun(wale_sec_mods, wale_dat, "water lettuce")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "AvgPercCovered_c",
                           "management" = "Treated")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Quarter)

non_foc_mod_se <- mod_se_fun(cubu_sec_mods, cubu_dat, "Cuban bulrush") %>%
  full_join(mod_se_fun(torp_sec_mods, torp_dat, "torpedograss")) %>%
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
write_csv(foc_mod_se, "output/fwc_focal_secchi_model_summary.csv")
write_csv(non_foc_mod_se, "output/fwc_non_focal_secchi_model_summary.csv")


#### values for text ####

# summarize uninvaded
uninv_sum <- uninv2 %>%
  group_by(CommonName, Quarter) %>%
  summarize(PrevValueUninv = mean(PrevValue),
            UninvN = n()) %>%
  ungroup()

# focal summaries
foc_sum <- tibble(CommonName = c("Hydrilla", "Water hyacinth", "Water hyacinth", "Water lettuce", "Water lettuce"),
                      Quarter = c(1, 1, 4, 1, 2),
                      DiffNone = c(mean(fixef(hydr_sec_mod_q1)), mean(fixef(wahy_sec_mod_q1)), mean(fixef(wahy_sec_mod_q4)), mean(fixef(wale_sec_mod_q1)), mean(fixef(wale_sec_mod_q2))),
                      PAC = as.numeric(c(coef(hydr_sec_mod_q1)[1], coef(wahy_sec_mod_q1)[1], coef(wahy_sec_mod_q4)[1], coef(wale_sec_mod_q1)[1], coef(wale_sec_mod_q2)[1])),
                  Treat= as.numeric(c(coef(hydr_sec_mod_q1)[2], coef(wahy_sec_mod_q1)[2], coef(wahy_sec_mod_q4)[2], coef(wale_sec_mod_q1)[2], coef(wale_sec_mod_q2)[2]))) %>%
  mutate(DiffPAC = DiffNone + PAC,
         DiffTreat = DiffNone + Treat) %>%
  left_join(hydr_dat %>%
              group_by(CommonName, Quarter) %>%
              summarize(PrevValue = mean(PrevValue)) %>%
              ungroup() %>%
              full_join(wahy_dat %>%
                          group_by(CommonName, Quarter) %>%
                          summarize(PrevValue = mean(PrevValue)) %>%
                          ungroup()) %>%
              full_join(wale_dat %>%
                          group_by(CommonName, Quarter) %>%
                          summarize(PrevValue = mean(PrevValue)) %>%
                          ungroup())) %>%
  left_join(uninv_sum) %>%
  mutate(across(.cols = c(DiffNone, PAC, Treat, DiffPAC, DiffTreat, PrevValue, PrevValueUninv), ~ .x * 30.48)) # convert from ft to cm

write_csv(foc_sum, "output/fwc_focal_invasive_secchi_prediction.csv")

# non-focal summaries
non_foc_sum <- tibble(CommonName = "Torpedograss",
                      Quarter = 1:4,
                      DiffNone = c(mean(fixef(torp_sec_mod_q1)), mean(fixef(torp_sec_mod_q2)), mean(fixef(torp_sec_mod_q3)), mean(fixef(torp_sec_mod_q4))),
                      PAC = as.numeric(c(coef(torp_sec_mod_q1)[1], coef(torp_sec_mod_q2)[1], coef(torp_sec_mod_q3)[1], coef(torp_sec_mod_q4)[1]))) %>%
  mutate(DiffPAC = DiffNone + PAC) %>%
  left_join(torp_dat %>%
              group_by(CommonName, Quarter) %>%
              summarize(PrevValue = mean(PrevValue)) %>%
              ungroup()) %>%
  left_join(uninv_sum) %>%
  mutate(across(.cols = c(DiffNone, PAC, DiffPAC, PrevValue, PrevValueUninv), ~ .x * 30.48)) # convert from ft to cm

write_csv(non_foc_sum, "output/fwc_non_focal_invasive_secchi_prediction.csv")
