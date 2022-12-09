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
library(car) # logit

# figure settings
source("code/settings/figure_settings.R")

# functions
source("code/generic-functions/model_structure_comparison.R")

# import data
qual <- read_csv("intermediate-data/water_quality_data_formatted.csv")
uninv <- read_csv("intermediate-data/water_quality_fwc_uninvaded_permID.csv")


#### edit data ####

# focal water quality
dat <- filter(qual, QualityMetric == "TP_ug_L")

# taxa
inv_taxa <- sort(unique(dat$CommonName))

# loop through taxa
pdf("output/phosphorus_time_series_by_taxon.pdf")

for(i in inv_taxa){
  
  # subset data
  subdat <- dat %>% filter(CommonName == i)
  subdat_ctrl <- subdat %>% filter(Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = QualityMean, color = PermanentID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          labs(x = "Year", y = "total phosphorus (ug/L)", title = i) +
          def_theme_paper +
          theme(plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()

# save data
write_csv(dat, "intermediate-data/FWC_phosphorus_analysis_formatted.csv")

# split by species (para grass doesn't have enough data)
hydr_dat <- filter(dat, CommonName == "Hydrilla")
torp_dat <- filter(dat, CommonName == "Torpedograss")
cubu_dat <- filter(dat, CommonName == "Cuban bulrush")
flpl_dat <- filter(dat, CommonName == "floating plants")


#### initial visualizations ####

# covariate correlations
dat %>%
  select(CommonName, PercCovered, Treated,
         PropTreated, RecentTreatment) %>%
  group_by(CommonName) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & abs(corr) >= 0.4) %>%
  data.frame()
# hydrilla: PropTreated and PercCovered 0.45

# response distributions
ggplot(dat, aes(x = QualityMean)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# skewed

ggplot(dat, aes(x = log(QualityMean))) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = QualityMax)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# skewed

ggplot(dat, aes(x = log(QualityMax))) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = QualityMin)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# skewed

ggplot(dat, aes(x = log(QualityMin))) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = QualityCV)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# skewed

ggplot(dat, aes(x = log(QualityCV))) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

# coefficients and response values
ggplot(dat, aes(x = logit(PercCovered), y = LogQualityMean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = logit(PercCovered), y = LogQualityMean, color = PermanentID)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")

dat %>%
  filter(logit(PercCovered) < -3) %>%
  ggplot(aes(x = logit(PercCovered), y = LogQualityMean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = Treated, y = LogQualityMean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat, aes(x = Lag3Treated, y = LogQualityMean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# all missing values are from beginning of dataset

ggplot(dat, aes(x = RecentTreatment, y = LogQualityMean)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  # geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# P increases the more recently lake has been treated for floating plants

ggplot(dat, aes(x = RecentTreatment, y = LogQualityMean, color = PermanentID)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")
# floating plant trend may be driven by lake/treatment confounding


#### evaluate model structure ####

# function to fit models for each species
mod_structure_fits <- function(dat_in){
  
  # create fixed effects data frame
  # individual-time rows
  # each waterbody is an individual
  dat_fix <- dat_in %>%
    pdata.frame(index = c("PermanentID", "GSYear"))
  
  # simple lm
  mod_lm <- lm(LogQualityMean ~ PercCovered + Treated, data = dat_fix)
  
  # random effects
  mod_ran_loc <- glmmTMB(LogQualityMean ~ PercCovered + Treated + (1|PermanentID), data = dat_fix)
  mod_ran_yr <- glmmTMB(LogQualityMean ~ PercCovered + Treated + (1|GSYear), data = dat_fix)
  mod_ran_loc_yr <- glmmTMB(LogQualityMean ~ PercCovered + Treated + (1|PermanentID) + (1|GSYear), data = dat_fix)
  
  # fixed effects
  mod_fix_loc <- plm(LogQualityMean ~ PercCovered + Treated, data = dat_fix,
                      model = "within")
  mod_fix_loc_yr <- plm(LogQualityMean ~ PercCovered + Treated, data = dat_fix,
                         model = "within", effect = "twoways")
  
  # add quadratic term to account for multiple processes
  # mod_fix_quad <- plm(LogQualityMean ~ PercCovered + Lag3APCsq + Lag3Treated, data = dat_fix,
  #                     model = "within", effect = "twoways")
  # tried this, need more years to do a squared term
  
  # use last treatment
  mod_fix_rec <- plm(LogQualityMean ~ PercCovered + RecentTreatment, data = dat_fix,
                      model = "within", effect = "twoways")
  
  # treatment-plant interactions
  mod_fix_int <- plm(LogQualityMean ~ PercCovered * Treated, data = dat_fix,
                      model = "within", effect = "twoways")
  mod_fix_rec_int <- plm(LogQualityMean ~ PercCovered * RecentTreatment, data = dat_fix,
                     model = "within", effect = "twoways")
  
  # return list of models
  return(list(lm = mod_lm,
              ran_loc = mod_ran_loc,
              ran_yr = mod_ran_yr,
              ran_loc_yr = mod_ran_loc_yr,
              fix_loc = mod_fix_loc,
              fix_loc_yr = mod_fix_loc_yr,
              rec = mod_fix_rec,
              int = mod_fix_int,
              rec_int = mod_fix_rec_int))
  
}

# fit models for each species
hydr_mod_struc <- mod_structure_fits(hydr_dat)
flpl_mod_struc <- mod_structure_fits(flpl_dat)
cubu_mod_struc <- mod_structure_fits(cubu_dat)
torp_mod_struc <- mod_structure_fits(torp_dat)

# compare model estimates
hydr_mod_comp <- mod_structure_comp(simp_mods = hydr_mod_struc[1], 
                                    ran_mods = hydr_mod_struc[2:4],
                                    fix_mods = hydr_mod_struc[5:9])
flpl_mod_comp <- mod_structure_comp(simp_mods = flpl_mod_struc[1], 
                                    ran_mods = flpl_mod_struc[2:4],
                                    fix_mods = flpl_mod_struc[5:9])
cubu_mod_comp <- mod_structure_comp(simp_mods = cubu_mod_struc[1], 
                                    ran_mods = cubu_mod_struc[2:4],
                                    fix_mods = cubu_mod_struc[5:9])
torp_mod_comp <- mod_structure_comp(simp_mods = torp_mod_struc[1], 
                                    ran_mods = torp_mod_struc[2:4],
                                    fix_mods = torp_mod_struc[5:9])

# combine species
mod_comp <- hydr_mod_comp %>%
  mutate(Species = "hydrilla") %>%
  full_join(flpl_mod_comp %>%
              mutate(Species = "floating plants")) %>%
  full_join(cubu_mod_comp %>%
              mutate(Species = "Cuban bulrush")) %>%
  full_join(torp_mod_comp %>%
              mutate(Species = "torpedograss")) %>%
  mutate(across(!c(coefficients, Species), ~ round(.x, digits = 3))) %>%
  relocate(Species)

write_csv(mod_comp, "output/fwc_phosphorus_model_structure_comparison.csv")

# model comparison notes:
# hydrilla and treatment increase P until location is accounted for, then they decrease
# recent treatment has a weaker effect than treated
# interaction is minimal

# accounting for location reduces the effects of floating plants and treatment on P
# both have positive effect
# recent treatment has a stronger effect
# treatment eliminates the positive effect of the plants

# accounting for location and time don't seem to affect estimates much
# recent treatment is also similar to treated
# treatment and plants both reduce P
# synergistic effect of the two (even more reduction)

# accounting for location reduces the effects
# both reduce P
# rectn similar to treated
# interaction accounts for most of the treated effects

# test fixed effects (seems like year isn't necessary)
# have to refit because data need to be accessible (not "dat_fix")
hydr_mod_diff_fix_loc_yr <- plm(LogQualityMean ~ PercCovered * RecentTreatment, 
                                data = hydr_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(hydr_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(hydr_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

flpl_mod_diff_fix_loc_yr <- plm(LogQualityMean ~ PercCovered * RecentTreatment, 
                                data = flpl_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(flpl_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(flpl_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

cubu_mod_diff_fix_loc_yr <- plm(LogQualityMean ~ PercCovered * RecentTreatment, 
                                data = cubu_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(cubu_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(cubu_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

torp_mod_diff_fix_loc_yr <- plm(LogQualityMean ~ PercCovered * RecentTreatment, 
                                data = torp_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(torp_mod_diff_fix_loc_yr, effect = "time", type = "bp") # not sig
plmtest(torp_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

# use waterbody fixed effects
# decided to leave in year effects because there should in theory be a year effect


# #### model-fitting functions ####
# 
# # data filter function
# dat_mod_filt <- function(treat_col, inv_col, dat_in){
# 
#   dat_mod <- dat_in %>%
#     mutate(Treated = !!sym(treat_col),
#            AvgPercCovered = !!sym(inv_col))
# 
#   return(dat_mod)
# 
# }
# 
# # function to fit models
# mod_fit <- function(dat_in){
# 
#   # focal species
#   foc_sp <- unique(dat_in$CommonName)
# 
#   # subset data
#   dat_mod1 <- dat_mod_filt("LastTreatment", "Lag1AvgPercCovered", dat_in)
#   dat_mod2 <- dat_mod_filt("LastTreatment", "Lag2AvgPercCovered", dat_in)
#   dat_mod3 <- dat_mod_filt("LastTreatment", "Lag3AvgPercCovered", dat_in)
#   dat_mod4 <- dat_mod_filt("LastTreatment", "Lag4AvgPercCovered", dat_in)
#   dat_mod5 <- dat_mod_filt("LastTreatment", "Lag5AvgPercCovered", dat_in)
# 
#   # fit models
#   if(foc_sp == "Cuban bulrush") {
# 
#     mod1 <- plm(logQual ~ AvgPercCovered * Treated, data = dat_mod1,
#                 index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
#     mod2 <- update(mod1, data = dat_mod2)
#     mod3 <- update(mod1, data = dat_mod3)
# 
#     mods_out <- list(mod1, mod2, mod3)
# 
#   } else {
# 
#     mod1 <- plm(logQual ~ AvgPercCovered * Treated, data = dat_mod1,
#                 index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
#     mod2 <- update(mod1, data = dat_mod2)
#     mod3 <- update(mod1, data = dat_mod3)
#     mod4 <- update(mod1, data = dat_mod4)
#     mod5 <- update(mod1, data = dat_mod5)
#     mod6 <- update(mod1, data = dat_mod6)
# 
#     mods_out <- list(mod1, mod2, mod3, mod4, mod5, mod6)
# 
#   }
# 
#   # output
#   return(mods_out)
# 
# }
# 
# 
# #### fit models ####
# 
# # fit models with all lags
# hydr_mods_q1 <- mod_fit(filter(hydr_dat, Quarter == 1))
# flpl_mods_q1 <- mod_fit(filter(flpl_dat, Quarter == 1))
# torp_mods_q1 <- mod_fit(filter(torp_dat, Quarter == 1))
# cubu_mods_q1 <- mod_fit(filter(cubu_dat, Quarter == 1))
# 
# hydr_mods_q2 <- mod_fit(filter(hydr_dat, Quarter == 2))
# flpl_mods_q2 <- mod_fit(filter(flpl_dat, Quarter == 2))
# torp_mods_q2 <- mod_fit(filter(torp_dat, Quarter == 2))
# cubu_mods_q2 <- mod_fit(filter(cubu_dat, Quarter == 2))
# 
# hydr_mods_q3 <- mod_fit(filter(hydr_dat, Quarter == 3))
# flpl_mods_q3 <- mod_fit(filter(flpl_dat, Quarter == 3))
# torp_mods_q3 <- mod_fit(filter(torp_dat, Quarter == 3))
# cubu_mods_q3 <- mod_fit(filter(cubu_dat, Quarter == 3))
# 
# hydr_mods_q4 <- mod_fit(filter(hydr_dat, Quarter == 4))
# flpl_mods_q4 <- mod_fit(filter(flpl_dat, Quarter == 4))
# torp_mods_q4 <- mod_fit(filter(torp_dat, Quarter == 4))
# cubu_mods_q4 <- mod_fit(filter(cubu_dat, Quarter == 4))
# 
# # name models
# names(hydr_mods_q1) <- names(flpl_mods_q1) <- names(torp_mods_q1) <- names(hydr_mods_q2) <- names(flpl_mods_q2) <- names(torp_mods_q2) <- names(hydr_mods_q3) <- names(flpl_mods_q3) <- names(torp_mods_q3) <- names(hydr_mods_q4) <- names(flpl_mods_q4) <- names(torp_mods_q4) <- c("1", "2", "3", "4", "5", "6")
# 
# names(cubu_mods_q1) <- names(cubu_mods_q2) <- names(cubu_mods_q3) <- names(cubu_mods_q4) <- c("1", "2", "3")
# 
# 
# #### coefficient figures and tables ####
# 
# # rename coefficients
# coef_names <- c("Treated" = "Management",
#                 "AvgPercCovered" = "Invasive PAC",
#                 "AvgPercCovered:Treated" = "interaction")
# 
# # ggplot function
# plot_fun <- function(models){
# 
#   plot_out <- modelplot(models,
#                         coef_map = coef_names,
#                         background = list(geom_vline(xintercept = 0, color = "black",
#                                                      size = 0.5, linetype = "dashed"))) +
#     scale_color_viridis_d(direction = -1) +
#     scale_x_continuous(labels = scale_fun_1) +
#     def_theme_paper +
#     theme(plot.title = element_text(size = 9))
# 
#   return(plot_out)
# 
# }
# 
# # panel plot function
# panel_plot_fun <- function(mods1, mods2,
#                                    spp1, spp2,
#                                    filename){
# 
#   # focal panels
#   fig1 <- plot_fun(mods1) +
#     labs(x = "",
#          title = paste("(A)", spp1)) +
#     theme(legend.position = "none")
# 
#   fig2 <- plot_fun(mods2) +
#     labs(x = "",
#          title = paste("(B)", spp2)) +
#     theme(axis.text.y = element_blank(),
#           legend.box.margin = margin(-10, 0, -10, -10)) +
#     scale_color_viridis_d(direction = -1, name = "Lag\n(years)") +
#     guides(color = guide_legend(reverse = TRUE))
# 
#   comb_fig <- (fig1 + fig2) +
#     plot_annotation(
#       caption = expression(paste("Estimate"%+-%" 95% CI", sep = "")),
#       theme = theme(plot.caption = element_text(size = 9, color="black", hjust = 0.6, vjust = 10),
#                     plot.margin = margin(5, -5, -5, -10),
#                     plot.title = element_text(size = 10, hjust = 0.5)),
#       title = "Effects on total phosphorus")
# 
#   ggsave(filename, comb_fig,
#          device = "eps", width = 4.7, height = 3, units = "in")
# 
# }
# 
# # figures
# panel_plot_fun(hydr_mods_q1, flpl_mods_q1,
#                "Hydrilla", "floating plants",
#                "output/fwc_focal_phosphorus_quarter1_diff_model.eps")
# panel_plot_fun(hydr_mods_q2, flpl_mods_q2,
#                "Hydrilla", "floating plants",
#                "output/fwc_focal_phosphorus_quarter2_diff_model.eps")
# panel_plot_fun(hydr_mods_q3, flpl_mods_q3,
#                "Hydrilla", "floating plants",
#                "output/fwc_focal_phosphorus_quarter3_diff_model.eps")
# panel_plot_fun(hydr_mods_q4, flpl_mods_q4,
#                "Hydrilla", "floating plants",
#                "output/fwc_focal_phosphorus_quarter4_diff_model.eps")
# 
# panel_plot_fun(cubu_mods_q1, torp_mods_q1,
#                        "Cuban bulrush", "Torpedograss",
#                        "output/fwc_non_focal_phosphorus_quarter1_diff_model.eps")
# 
# panel_plot_fun(cubu_mods_q2, torp_mods_q2,
#                        "Cuban bulrush", "Torpedograss",
#                        "output/fwc_non_focal_phosphorus_quarter2_diff_model.eps")
# panel_plot_fun(cubu_mods_q3, torp_mods_q3,
#                        "Cuban bulrush", "Torpedograss",
#                        "output/fwc_non_focal_phosphorus_quarter3_diff_model.eps")
# panel_plot_fun(cubu_mods_q4, torp_mods_q4,
#                        "Cuban bulrush", "Torpedograss",
#                        "output/fwc_non_focal_phosphorus_quarter4_diff_model.eps")


#### finalize models ####

# format data
hydr_dat_fix <- hydr_dat %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
flpl_dat_fix <- flpl_dat %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
cubu_dat_fix <- cubu_dat %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
torp_dat_fix <- torp_dat %>%
  pdata.frame(index = c("PermanentID", "GSYear"))

# fit models
hydr_mean_mod <- plm(LogQualityMean ~ PercCovered * RecentTreatment, data = hydr_dat_fix,
                       index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
flpl_mean_mod <- update(hydr_mean_mod, data = flpl_dat_fix)
cubu_mean_mod <- update(hydr_mean_mod, data = cubu_dat_fix)
torp_mean_mod <- update(hydr_mean_mod, data = torp_dat_fix)

hydr_max_mod <- plm(LogQualityMax ~ PercCovered * RecentTreatment, data = hydr_dat_fix,
                     index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
flpl_max_mod <- update(hydr_max_mod, data = flpl_dat_fix)
cubu_max_mod <- update(hydr_max_mod, data = cubu_dat_fix)
torp_max_mod <- update(hydr_max_mod, data = torp_dat_fix)

hydr_min_mod <- plm(LogQualityMin ~ PercCovered * RecentTreatment, data = hydr_dat_fix,
                     index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
flpl_min_mod <- update(hydr_min_mod, data = flpl_dat_fix)
cubu_min_mod <- update(hydr_min_mod, data = cubu_dat_fix)
torp_min_mod <- update(hydr_min_mod, data = torp_dat_fix)

hydr_cv_mod <- plm(LogQualityCV ~ PercCovered * RecentTreatment, data = hydr_dat_fix,
                     index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
flpl_cv_mod <- update(hydr_cv_mod, data = flpl_dat_fix)
cubu_cv_mod <- update(hydr_cv_mod, data = cubu_dat_fix)
torp_cv_mod <- update(hydr_cv_mod, data = torp_dat_fix)

# add fitted values to pdata.frame (important to match rows)
# convert to regular dataframe
hydr_mean_fit <- mutate(hydr_dat_fix, Fitted = as.numeric(hydr_dat_fix$LogQualityMean - hydr_mean_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_mean_fit <- mutate(flpl_dat_fix, Fitted = as.numeric(flpl_dat_fix$LogQualityMean - flpl_mean_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_mean_fit <- mutate(cubu_dat_fix, Fitted = as.numeric(cubu_dat_fix$LogQualityMean - cubu_mean_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_mean_fit <- mutate(torp_dat_fix, Fitted = as.numeric(torp_dat_fix$LogQualityMean - torp_mean_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_max_fit <- mutate(hydr_dat_fix, Fitted = as.numeric(hydr_dat_fix$LogQualityMax - hydr_max_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_max_fit <- mutate(flpl_dat_fix, Fitted = as.numeric(flpl_dat_fix$LogQualityMax - flpl_max_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_max_fit <- mutate(cubu_dat_fix, Fitted = as.numeric(cubu_dat_fix$LogQualityMax - cubu_max_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_max_fit <- mutate(torp_dat_fix, Fitted = as.numeric(torp_dat_fix$LogQualityMax - torp_max_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_min_fit <- mutate(hydr_dat_fix, Fitted = as.numeric(hydr_dat_fix$LogQualityMin - hydr_min_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_min_fit <- mutate(flpl_dat_fix, Fitted = as.numeric(flpl_dat_fix$LogQualityMin - flpl_min_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_min_fit <- mutate(cubu_dat_fix, Fitted = as.numeric(cubu_dat_fix$LogQualityMin - cubu_min_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_min_fit <- mutate(torp_dat_fix, Fitted = as.numeric(torp_dat_fix$LogQualityMin - torp_min_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)

hydr_cv_fit <- mutate(hydr_dat_fix, Fitted = as.numeric(hydr_dat_fix$LogQualityCV - hydr_cv_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
flpl_cv_fit <- mutate(flpl_dat_fix, Fitted = as.numeric(flpl_dat_fix$LogQualityCV - flpl_cv_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
cubu_cv_fit <- mutate(cubu_dat_fix, Fitted = as.numeric(cubu_dat_fix$LogQualityCV - cubu_cv_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)
torp_cv_fit <- mutate(torp_dat_fix, Fitted = as.numeric(torp_dat_fix$LogQualityCV - torp_cv_mod$residuals)) %>%
  as.data.frame(keep.attributes = F)

# fitted vs. observed
ggplot(hydr_mean_fit, aes(x = Fitted, y = LogQualityMean)) + geom_point()
ggplot(flpl_mean_fit, aes(x = Fitted, y = LogQualityMean)) + geom_point()
ggplot(cubu_mean_fit, aes(x = Fitted, y = LogQualityMean)) + geom_point()
ggplot(torp_mean_fit, aes(x = Fitted, y = LogQualityMean)) + geom_point()

ggplot(hydr_max_fit, aes(x = Fitted, y = LogQualityMax)) + geom_point()
ggplot(flpl_max_fit, aes(x = Fitted, y = LogQualityMax)) + geom_point()
ggplot(cubu_max_fit, aes(x = Fitted, y = LogQualityMax)) + geom_point()
ggplot(torp_max_fit, aes(x = Fitted, y = LogQualityMax)) + geom_point()

ggplot(hydr_min_fit, aes(x = Fitted, y = LogQualityMin)) + geom_point()
ggplot(flpl_min_fit, aes(x = Fitted, y = LogQualityMin)) + geom_point()
ggplot(cubu_min_fit, aes(x = Fitted, y = LogQualityMin)) + geom_point()
ggplot(torp_min_fit, aes(x = Fitted, y = LogQualityMin)) + geom_point()

ggplot(hydr_cv_fit, aes(x = Fitted, y = LogQualityCV)) + geom_point()
ggplot(flpl_cv_fit, aes(x = Fitted, y = LogQualityCV)) + geom_point()
ggplot(cubu_cv_fit, aes(x = Fitted, y = LogQualityCV)) + geom_point()
ggplot(torp_cv_fit, aes(x = Fitted, y = LogQualityCV)) + geom_point()
# much weaker relationships with CV

# export models
save(hydr_mean_mod, file = "output/fwc_hydrilla_mean_phosphorus_model.rda")
save(flpl_mean_mod, file = "output/fwc_floating_plant_mean_phosphorus_model.rda")
save(cubu_mean_mod, file = "output/fwc_cuban_bulrush_mean_phosphorus_model.rda")
save(torp_mean_mod, file = "output/fwc_torpedograss_mean_phosphorus_model.rda")

save(hydr_max_mod, file = "output/fwc_hydrilla_max_phosphorus_model.rda")
save(flpl_max_mod, file = "output/fwc_floating_plant_max_phosphorus_model.rda")
save(cubu_max_mod, file = "output/fwc_cuban_bulrush_max_phosphorus_model.rda")
save(torp_max_mod, file = "output/fwc_torpedograss_max_phosphorus_model.rda")

save(hydr_min_mod, file = "output/fwc_hydrilla_min_phosphorus_model.rda")
save(flpl_min_mod, file = "output/fwc_floating_plant_min_phosphorus_model.rda")
save(cubu_min_mod, file = "output/fwc_cuban_bulrush_min_phosphorus_model.rda")
save(torp_min_mod, file = "output/fwc_torpedograss_min_phosphorus_model.rda")

save(hydr_cv_mod, file = "output/fwc_hydrilla_cv_phosphorus_model.rda")
save(flpl_cv_mod, file = "output/fwc_floating_plant_cv_phosphorus_model.rda")
save(cubu_cv_mod, file = "output/fwc_cuban_bulrush_cv_phosphorus_model.rda")
save(torp_cv_mod, file = "output/fwc_torpedograss_cv_phosphorus_model.rda")

# process model SE
mod_se_fun <- function(model, dat, spp){
  
  dat_out <- tidy(coeftest(model, vcov = vcovHC, type = "HC3")) %>%
    mutate(R2 = r.squared(model),
           Invasive = spp,
           Waterbodies = n_distinct(dat$PermanentID),
           Years = n_distinct(dat$GSYear),
           N = Waterbodies * Years)
  
  return(dat_out)
  
}

# combine SE tables
foc_mod_se <- mod_se_fun(hydr_mean_mod, hydr_dat, "hydrilla") %>%
  mutate(Response = "mean") %>%
  full_join(mod_se_fun(hydr_max_mod, hydr_dat, "hydrilla") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(hydr_min_mod, hydr_dat, "hydrilla") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(hydr_cv_mod, hydr_dat, "hydrilla") %>%
              mutate(Response = "cv")) %>%
  full_join(mod_se_fun(flpl_mean_mod, flpl_dat, "floating plant") %>%
              mutate(Response = "mean")) %>%
  full_join(mod_se_fun(flpl_max_mod, flpl_dat, "floating plant") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(flpl_min_mod, flpl_dat, "floating plant") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(flpl_cv_mod, flpl_dat, "floating plant") %>%
              mutate(Response = "cv")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "PercCovered",
                           "management" = "RecentTreatment",
                           "invasive PAC:management" = "PercCovered:RecentTreatment")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Response)

non_foc_mod_se <- mod_se_fun(cubu_mean_mod, cubu_dat, "Cuban bulrush") %>%
  mutate(Response = "mean") %>%
  full_join(mod_se_fun(cubu_max_mod, cubu_dat, "Cuban bulrush") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(cubu_min_mod, cubu_dat, "Cuban bulrush") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(cubu_cv_mod, cubu_dat, "Cuban bulrush") %>%
              mutate(Response = "cv")) %>%
  full_join(mod_se_fun(torp_mean_mod, torp_dat, "torpedograss") %>%
              mutate(Response = "mean")) %>%
  full_join(mod_se_fun(torp_max_mod, torp_dat, "torpedograss") %>%
              mutate(Response = "max")) %>%
  full_join(mod_se_fun(torp_min_mod, torp_dat, "torpedograss") %>%
              mutate(Response = "min")) %>%
  full_join(mod_se_fun(torp_cv_mod, torp_dat, "torpedograss") %>%
              mutate(Response = "cv")) %>%
  mutate(term = fct_recode(term,
                           "invasive PAC" = "PercCovered",
                           "management" = "RecentTreatment",
                           "invasive PAC:management" = "PercCovered:RecentTreatment")) %>%
  rename(Term = term,
         Estimate = estimate,
         SE = std.error,
         t = statistic,
         P = p.value) %>%
  relocate(Invasive, Response)

# export
write_csv(foc_mod_se, "output/fwc_focal_phosphorus_model_summary.csv")
write_csv(non_foc_mod_se, "output/fwc_non_focal_phosphorus_model_summary.csv")



#### values for text ####

# summarize uninvaded
uninv_sum <- uninv2 %>%
  group_by(CommonName, Quarter) %>%
  summarize(UninvAvg = mean(QualityValue),
            UninvN = n()) %>%
  ungroup()

# translate model coefficients
mod_coef_fun <- function(models, spp){
  
  dat_out <- tibble(Invasive = spp,
                    Quarter = c("Apr-Jun", "Jul-Sep", "Oct-Dec", "Jan-Mar"),
                    Intercept = sapply(models, function(x) mean(fixef(x))),
                    BetaPAC = sapply(models, function(x) coef(x)[1]),
                    BetaTreat = sapply(models, function(x) coef(x)[2]))
  
  return(dat_out)
  
}

# identify significant effects
foc_sig <- foc_mod_se %>%
  filter(P < 0.1) %>%
  select(Invasive, Quarter, Term) %>%
  left_join(mod_coef_fun(hydr_pho_mods, "hydrilla") %>%
              full_join(mod_coef_fun(wahy_pho_mods, "water hyacinth")) %>%
              full_join(mod_coef_fun(wale_pho_mods, "water lettuce"))) %>%
  mutate(BetaPAC = if_else(Term == "management", NA_real_, BetaPAC),
         BetaTreat = if_else(Term == "invasive PAC", NA_real_, BetaTreat),
         Metric = "total phosphorus") %>%
  left_join(uninv_sum %>%
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
  mutate(BetaPAC = if_else(Term == "management", NA_real_, BetaPAC),
         BetaTreat = if_else(Term == "invasive PAC", NA_real_, BetaTreat),
         Metric = "total phosphorus") %>%
  left_join(uninv_sum %>%
              mutate(Quarter = case_when(Quarter == 1 ~ "Apr-Jun", 
                                         Quarter == 2 ~ "Jul-Sep", 
                                         Quarter == 3 ~ "Oct-Dec", 
                                         Quarter == 4 ~ "Jan-Mar"),
                     CommonName = fct_recode(CommonName, 
                                             "torpedograss" = "Torpedograss")) %>%
              rename(Invasive = CommonName))

write_csv(non_foc_sig, "output/fwc_non_focal_invasive_phosphorus_significant.csv")
