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

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv")
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")


#### edit data ####

# combine datasets
inv_dat <- inv_plant %>%
  filter(!is.na(LogRatioCovered)) %>% # need two consecutive years
  inner_join(inv_ctrl) %>%
  mutate(InitPACBin = case_when(InitPercCoveredAdj < 1 ~ "< 1%",
                                InitPercCoveredAdj >= 1 & InitPercCoveredAdj < 10 ~ "1%-10%",
                                TRUE ~ "≥ 10%") %>%
           fct_relevel("< 1%", "1%-10%"),
         PercChangePAC = InitPercCoveredAdj * (exp(LogRatioCovered) - 1)) # final % - initial %

# will need surveyor experience
filter(inv_dat, is.na(MinSurveyorExperience))
# none missing

# check data availability
inv_dat %>%
  filter(PropCovered > 0) %>%
  ggplot(aes(x = GSYear, y = PropCovered, color = PermanentID)) +
  geom_vline(xintercept = 2014) +
  geom_point() +
  facet_wrap(~ CommonName, scales = "free_y") +
  theme(legend.position = "none")
# Cuban bulrush is missing a lot of data before 2014
# Para grass and torpedograss start in 2000

# split by species
hydr_dat <- filter(inv_dat, CommonName == "Hydrilla")
wale_dat <- filter(inv_dat, CommonName == "Water lettuce")
wahy_dat <- filter(inv_dat, CommonName == "Water hyacinth")
torp_dat <- filter(inv_dat, CommonName == "Torpedograss")
cubu_dat <- filter(inv_dat, CommonName == "Cuban bulrush" & GSYear >= 2014)
pagr_dat <- filter(inv_dat, CommonName == "Para grass")

# check Cuban bulrush
cubu_dat %>%
  filter(PropCovered > 0) %>%
  ggplot(aes(x = GSYear, y = PropCovered, color = PermanentID)) +
  geom_line() +
  facet_wrap(~ CommonName, scales = "free_y") +
  theme(legend.position = "none")


#### initial visualizations ####

# all and species-specific control correlated?
cor.test(~ Lag1Treated + Lag1AllTreated, data = hydr_dat) # 0.5
cor.test(~ Lag1Treated + Lag1AllTreated, data = wale_dat) # 0.6
cor.test(~ Lag1Treated + Lag1AllTreated, data = torp_dat) # 0.2 
cor.test(~ Lag1Treated + Lag1AllTreated, data = cubu_dat) # 0.3
cor.test(~ Lag1Treated + Lag1AllTreated, data = pagr_dat) # 0.04
# yes

# covariate correlations
hydr_dat %>%
  select(Lag1Treated, InitPercCoveredAdj, MinSurveyorExperience) %>%
  ggpairs()

wahy_dat %>%
  select(Lag1Treated, InitPercCoveredAdj, MinSurveyorExperience) %>%
  ggpairs() # one high cover value

wahy_dat %>%
  select(Lag1Treated, InitPercCoveredAdj, MinSurveyorExperience) %>%
  ggpairs() # one high cover value

torp_dat %>%
  select(Lag1Treated, InitPercCoveredAdj, MinSurveyorExperience) %>%
  ggpairs()

cubu_dat %>%
  select(Lag1Treated, InitPercCoveredAdj, MinSurveyorExperience) %>%
  ggpairs()

pagr_dat %>%
  select(Lag1Treated, InitPercCoveredAdj, MinSurveyorExperience) %>%
  ggpairs()

# change in cover
ggplot(inv_dat, aes(x = InitPercCovered, y = PropCovered)) +
  geom_point() +
  facet_wrap(~ CommonName, scales = "free") +
  geom_smooth(method = "lm")

# log ratio prop covered
ggplot(hydr_dat, aes(x = InitPercCoveredAdj, y = LogRatioCovered, 
                     color = as.factor(Lag1Treated))) +
  geom_point()

ggplot(wale_dat, aes(x = InitPercCoveredAdj, y = LogRatioCovered, 
                     color = as.factor(Lag1Treated))) +
  geom_point()

ggplot(wahy_dat, aes(x = InitPercCoveredAdj, y = LogRatioCovered, 
                     color = as.factor(Lag1Treated))) +
  geom_point()

ggplot(torp_dat, aes(x = InitPercCoveredAdj, y = LogRatioCovered, 
                     color = as.factor(Lag1Treated))) +
  geom_point()
# a lot untreated

ggplot(cubu_dat, aes(x = InitPercCoveredAdj, y = LogRatioCovered, 
                     color = as.factor(Lag1Treated))) +
  geom_point()

ggplot(pagr_dat, aes(x = InitPercCoveredAdj, y = LogRatioCovered, 
                     color = as.factor(Lag1Treated))) +
  geom_point()
# a lot untreated

# treated and change in prop
ggplot(hydr_dat, aes(x = Lag6Treated, y = LogRatioCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

ggplot(wale_dat, aes(x = Lag6Treated, y = LogRatioCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

ggplot(wahy_dat, aes(x = Lag6Treated, y = LogRatioCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

ggplot(torp_dat, aes(x = Lag6Treated, y = LogRatioCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

ggplot(cubu_dat, aes(x = Lag6Treated, y = LogRatioCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

ggplot(pagr_dat, aes(x = Lag6Treated, y = LogRatioCovered)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean")

# data availability for lags
hydr_dat %>%
  select(PermanentID, GSYear,
         Lag1Treated, Lag2Treated, Lag3Treated, Lag4Treated, Lag5Treated, Lag6Treated) %>%
  pivot_longer(cols = starts_with("Lag"),
               names_to = "Lag",
               values_to = "Treated") %>%
  filter(!is.na(Treated)) %>%
  ggplot(aes(x = Lag)) +
  geom_bar()

# initial PAC distributions
hydr_dat %>%
  mutate(InitPACBin = cut_number(InitPercCoveredAdj, n = 3)) %>%
  ggplot(aes(x = InitPACBin)) +
  geom_bar()
# intervals are too small to be meaningful

# initial PAC distributions
# tried different thresholds before choosing these
ggplot(hydr_dat, aes(x = InitPACBin)) +
  geom_bar()
ggplot(wale_dat, aes(x = InitPACBin)) +
  geom_bar()
ggplot(wahy_dat, aes(x = InitPACBin)) +
  geom_bar()
ggplot(torp_dat, aes(x = InitPACBin)) +
  geom_bar()
ggplot(cubu_dat, aes(x = InitPACBin)) +
  geom_bar()
ggplot(pagr_dat, aes(x = InitPACBin)) +
  geom_bar()

# % change in PAC
ggplot(hydr_dat, aes(x = PercChangePAC)) +
  geom_histogram()
ggplot(wale_dat, aes(x = PercChangePAC)) +
  geom_histogram()
ggplot(wahy_dat, aes(x = PercChangePAC)) +
  geom_histogram()
ggplot(torp_dat, aes(x = PercChangePAC)) +
  geom_histogram()
ggplot(cubu_dat, aes(x = PercChangePAC)) +
  geom_histogram()
ggplot(pagr_dat, aes(x = PercChangePAC)) +
  geom_histogram()


#### model-fitting functions ####

# data filter function
dat_mod_filt <- function(treat_col, dat_in){
  
  dat_mod <- dat_in %>%
    filter(!is.na(LogRatioCovered) & !is.na(InitPercCoveredAdj) & !is.na(MinSurveyorExperience) & !is.na(Lag1Treated) & !is.na(Lag2Treated) & !is.na(Lag3Treated) & !is.na(Lag4Treated) & !is.na(Lag5Treated) & !is.na(Lag6Treated)) %>%
    mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
           InitPercCoveredAdj_c = InitPercCoveredAdj - mean(InitPercCoveredAdj),
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
  mod1 <- feols(LogRatioCovered ~ InitPercCoveredAdj_c * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod1)
  mod2 <- feols(LogRatioCovered ~ InitPercCoveredAdj_c * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod2)
  mod3 <- feols(LogRatioCovered ~ InitPercCoveredAdj_c * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod3)
  mod4 <- feols(LogRatioCovered ~ InitPercCoveredAdj_c * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod4)
  mod5 <- feols(LogRatioCovered ~ InitPercCoveredAdj_c * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod5)
  mod6 <- feols(LogRatioCovered ~ InitPercCoveredAdj_c * Treated + SurveyorExperience_s | PermanentID + GSYear, data = dat_mod6)

  
  # output
  return(list(mod1, mod2, mod3, mod4, mod4, mod6))

}


#### fit models ####

# fit models
hydr_mods <- mod_fit(hydr_dat)
wahy_mods <- mod_fit(wahy_dat)
wale_mods <- mod_fit(wale_dat)
torp_mods <- mod_fit(torp_dat)
cubu_mods <- mod_fit(cubu_dat)
pagr_mods <- mod_fit(pagr_dat)

# name models
names(hydr_mods) <- names(wahy_mods) <- names(wale_mods) <- names(torp_mods) <- names(cubu_mods) <- names(pagr_mods) <- c("1", "2", "3", "4", "5", "6")


#### coefficient figures and tables ####

# rename coefficients
coef_names <- c("SurveyorExperience_s" = "Surveyor experience",
                "InitPercCoveredAdj_c:Treated" = "Management:PAC",
                "InitPercCoveredAdj_c" = "Initial PAC",
                "Treated" = "Management frequency")

# focal panels
hydr_fig <- modelplot(hydr_mods,
          coef_map = coef_names,
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1) +
  scale_x_continuous(labels = scale_fun_1) +
  labs(x = "",
       title = "(A) Hydrilla") +
  def_theme_paper +
  theme(legend.position = "none")

wahy_fig <- modelplot(wahy_mods,
          coef_map = coef_names,
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       title = "(B) Water hyacinth") +
  def_theme_paper +
  theme(legend.position = "none",
        axis.text.y = element_blank())

wale_fig <- modelplot(wale_mods,
          coef_map = coef_names,
          background = list(geom_vline(xintercept = 0, color = "black",
                                       size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1, name = "Management\nyears\nincluded") +
  labs(x = "",
       title = "(C) Water lettuce") +
  def_theme_paper +
  theme(axis.text.y = element_blank(),
        legend.box.margin = margin(-10, 0, -10, -10)) +
  guides(color = guide_legend(reverse = TRUE))

# combine focal panels
foc_fig <- hydr_fig + wahy_fig + wale_fig + plot_annotation(theme = theme(plot.margin = margin(0, -5, 0, -10)))
ggsave("output/fwc_focal_invasive_plant_treatment_model.eps", foc_fig,
       device = "eps", width = 6.5, height = 3.5, units = "in")

# export models
save(hydr_mods, file = "output/fwc_hydrilla_treatment_models.rda")
save(wahy_mods, file = "output/fwc_water_hyacinth_treatment_models.rda")
save(wale_mods, file = "output/fwc_water_lettuce_treatment_models.rda")

# non-focal panels
cubu_fig <- modelplot(cubu_mods,
                      coef_map = coef_names,
                      background = list(geom_vline(xintercept = 0, color = "black",
                                                   size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1) +
  labs(x = "", title = "(A) Cuban bulrush") +
  def_theme_paper +
  theme(legend.position = "none")

pagr_fig <- modelplot(pagr_mods,
                      coef_map = coef_names,
                      background = list(geom_vline(xintercept = 0, color = "black",
                                                   size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1) +
  labs(x = expression(paste("Estimate "%+-%" 95% CI", sep = "")),
       title = "(B) para grass") +
  def_theme_paper +
  theme(legend.position = "none",
        axis.text.y = element_blank())

torp_fig <- modelplot(torp_mods,
                      coef_map = coef_names,
                      background = list(geom_vline(xintercept = 0, color = "black",
                                                   size = 0.5, linetype = "dashed"))) +
  scale_color_viridis_d(direction = -1, name = "Management\nyears\nincluded") +
  labs(x = "",
       title = "(C) torpedograss") +
  def_theme_paper +
  theme(axis.text.y = element_blank(),
        legend.box.margin = margin(-10, 0, -10, -10)) +
  guides(color = guide_legend(reverse = TRUE))

# combine figures
nonfoc_fig <- cubu_fig + pagr_fig + torp_fig + plot_annotation(theme = theme(plot.margin = margin(0, -5, 0, -10)))
ggsave("output/fwc_non_focal_invasive_plant_treatment_model.eps", nonfoc_fig,
       device = "eps", width = 6.5, height = 3.5, units = "in")

# export models
save(cubu_mods, file = "output/fwc_cuban_bulrush_treatment_models.rda")
save(pagr_mods, file = "output/fwc_para_grass_treatment_models.rda")
save(torp_mods, file = "output/fwc_torpedograss_lettuce_treatment_models.rda")


#### model prediction figures ####

# function to format data
pred_fig <- function(treat_col, dat_in, mod){
  
  # edit data
  raw_dat <- dat_mod_filt(treat_col, dat_in) %>%
    mutate(PropChangePAC = exp(LogRatioCovered) - 1)
  
  # treatment frequency
  treat_n <- as.numeric(str_sub(treat_col, 4, 4))
  
  # prediction dataset
  pred_dat <- raw_dat %>%
    group_by(InitPACBin) %>%
    mutate(InitPercCoveredAdjMean = mean(InitPercCoveredAdj)) %>% # set to mean for group
    ungroup() %>%
    mutate(InitPercCoveredAdj_c = InitPercCoveredAdjMean - mean(raw_dat$InitPercCoveredAdj)) %>% # center on overall mean
    select(-Treated) %>%
    expand_grid(tibble(Treated = seq(0, 1, length.out = treat_n + 1))) %>% # repeat all observations across all treatment frequencies
    mutate(Pred = predict(mod, newdata = .),
           PercChangePAC = InitPercCoveredAdj * (exp(Pred) - 1),
           PropChangePAC = exp(Pred) - 1)
  
  # fig_out <- ggplot(pred_dat, aes(x = as.factor(round(Treated, 2)), y = PercChangePAC)) +
  #   geom_hline(yintercept = 0, size = 0.25) +
  #   # stat_summary(aes(color = InitPACBin),
  #   #              geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  #   # stat_summary(aes(color = InitPACBin),
  #   #              geom = "point", fun = "mean", size = 2) +
  #   # stat_summary(data = pred_dat, aes(color = InitPACBin),
  #   #              geom = "line", fun = "mean", size = 1) +
  #   ggdist::stat_halfeye(aes(color = InitPACBin),
  #                        ## custom bandwidth
  #                        adjust = .5, 
  #                        ## adjust height
  #                        width = .6, 
  #                        ## move geom to the right
  #                        justification = -.2, 
  #                        ## remove slab interval
  #                        .width = 0, 
  #                        point_colour = NA) +
  #   scale_color_viridis_d(name = "Initial\nabundance",
  #                         direction = -1) +
  #   labs(x = "Management frequency", y = "Change in PAC") +
  #   def_theme_paper
  # 
  return(pred_dat)
  
}

#### start here ####

# why are estimates so different from raw data?

# figures
# hydr_pred_fig <- 
  pred_fig("Lag6Treated", hydr_dat, hydr_mods[[6]])
  
  pred_fig("Lag6Treated", wahy_dat, hydr_mods[[6]])

test_dat <- pred_fig("Lag6Treated", hydr_dat, hydr_mods[[6]])

ggplot(test_dat, aes(x = PercChangePAC)) +
  geom_histogram() +
  facet_grid(Treated ~ InitPACBin, scales = "free") +
  def_theme

test_dat2 <- dat_mod_filt("Lag6Treated", hydr_dat) %>%
  group_by(InitPACBin) %>%
  mutate(InitPercCoveredAdjMean = mean(InitPercCoveredAdj)) %>% # set to mean for group
  ungroup() %>%
  mutate(InitPercCoveredAdj_c = InitPercCoveredAdjMean - mean(dat_mod_filt("Lag6Treated", hydr_dat)$InitPercCoveredAdj)) %>% # center on overall mean
  select(-Treated) %>%
  expand_grid(tibble(Treated = seq(0, 1, length.out = 6 + 1))) %>% # repeat all observations across all treatment frequencies
  mutate(Pred = predict(hydr_mods[[6]], newdata = .),
         PercChangePAC = InitPercCoveredAdj * (exp(Pred) - 1),
         PropChangePAC = exp(Pred) - 1)

ggplot(test_dat2, aes(x = PercChangePAC)) +
  geom_histogram() +
  facet_grid(Treated ~ InitPACBin, scales = "free") +
  def_theme


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
