#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(plotly)
library(tidyverse)
library(inspectdf) # inspect_cor
library(modelsummary) # modelplot
library(patchwork) # combining figures
library(car) # logit
library(plm) # panel data models
library(glmmTMB) # random effects
library(sandwich) # vcovHC
library(lmtest) # coeftest


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
  filter(!is.na(InitPercCovered)) %>% 
  inner_join(inv_ctrl) %>%
  mutate(PercDiffCovered = PropCovered * 100 - InitPercCovered,
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
  geom_vline(xintercept = 2013) +
  geom_point() +
  facet_wrap(~ CommonName, scales = "free_y") +
  theme(legend.position = "none")
# Cuban bulrush is missing a lot of data before 2013
# Para grass and torpedograss start in 1999

# complete time intervals
# surveys were not conducted every year on every lake for every species
# control data is implicitly complete -- missing interpretted as no control
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
  filter(GSYear >= MinGSYear & GSYear <= MaxGSYear & !is.na(Lag1Treated))

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

# export
write_csv(inv_dat4, "intermediate-data/FWC_invasive_plant_analysis_formatted.csv")

# split by species
hydr_dat <- filter(inv_dat4, CommonName == "Hydrilla")
wale_dat <- filter(inv_dat4, CommonName == "Water lettuce")
wahy_dat <- filter(inv_dat4, CommonName == "Water hyacinth")
torp_dat <- filter(inv_dat4, CommonName == "Torpedograss")
cubu_dat <- filter(inv_dat4, CommonName == "Cuban bulrush")
pagr_dat <- filter(inv_dat4, CommonName == "Para grass")


#### initial visualizations ####

# all and species-specific control correlated?
inv_dat4 %>%
  select(CommonName, Lag1Treated, Lag1AllTreated, MinSurveyorExperience, InitPercCovered) %>%
  group_by(CommonName) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & corr > 0.4) %>%
  data.frame()
# treated and all treated correlated for focal invasive species (expected)

# response distributions
ggplot(inv_dat4, aes(x = PercDiffCovered)) +
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

ggplot(inv_dat4, aes(x = InitPercCovered, y = PercDiffCovered,
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

# before/after comparison
ggplot(inv_dat4, aes(x = Lag1Treated, y = PercDiffCovered)) +
  geom_hline(yintercept = 0, size = 0.25) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_dat4, aes(x = Lag6Treated, y = PercDiffCovered)) +
  geom_hline(yintercept = 0, size = 0.25) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ CommonName, scales = "free")


#### evaluate model structure ####

# work from the premise that treatments, when done correctly, harm their target plants

# filter for data with Lag1Treated
hydr_lag1_dat <- hydr_dat %>% filter(!is.na(Lag1Treated))

# simple lm
mod_lm <- lm(PropCoveredLogit ~ Lag1Treated, data = hydr_lag1_dat)
summary(mod_lm)
# est. = 2.1

# random effects
mod_ran_loc <- glmmTMB(PropCoveredLogit ~ Lag1Treated + (1|PermanentID), data = hydr_lag1_dat)
summary(mod_ran_loc)
# reduces est. to 1.5
mod_ran_yr <- glmmTMB(PropCoveredLogit ~ Lag1Treated + (1|GSYear), data = hydr_lag1_dat)
summary(mod_ran_yr)
# keeps est. at 2.1
mod_ran_loc_yr <- glmmTMB(PropCoveredLogit ~ Lag1Treated + (1|PermanentID) + (1|GSYear), data = hydr_lag1_dat)
summary(mod_ran_loc_yr)
# reduces est. to 1.5
# location variance helps explain some of the positive relationship between treatment and cover

# fixed effects
mod_fix_loc <- plm(PropCoveredLogit ~ Lag1Treated, data = hydr_lag1_dat, 
                   index = c("PermanentID", "GSYear"), model = "within")
summary(mod_fix_loc)
# reduces est. to 1.4
mod_fix_loc_yr <- plm(PropCoveredLogit ~ Lag1Treated, data = hydr_lag1_dat, 
                      index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
summary(mod_fix_loc_yr)
# est. still 1.4

# use treatment that is more temporally removed from cover
hydr_lag3_dat <- hydr_dat %>% filter(!is.na(Lag3Treated)) # note that N changes
hydr_lag6_dat <- hydr_dat %>% filter(!is.na(Lag6Treated))
mod_fix_loc_yr3 <- plm(PropCoveredLogit ~ Lag3Treated, data = hydr_lag3_dat, 
                      index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
summary(mod_fix_loc_yr3)
# reduces est. to 0.5
mod_fix_loc_yr6 <- plm(PropCoveredLogit ~ Lag6Treated, data = hydr_lag6_dat, 
                       index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
summary(mod_fix_loc_yr6)
# reduces est. to -0.4
# is this significant when standard errors are corrected for heteroskedasticity and autocorrelation?
coeftest(mod_fix_loc_yr6, vcov = vcovHC, type = "HC1")
# no longer significant for < 0.05

# use initial cover to account for reverse causality
hydr_init_lag1_dat <- hydr_lag1_dat %>%
  mutate(InitPercCovered_c = InitPercCovered - mean(InitPercCovered))
mod_init_fix_loc_yr <- plm(PropCoveredLogit ~ InitPercCovered_c * Lag1Treated, data = hydr_init_lag1_dat,
                           index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
summary(mod_init_fix_loc_yr)
# reduces treatment est. to 1
# for each percentage increase in cover, treatment effect decreases by 0.03
# should get a negative treatment effect for ~40% increase in PAC above mean

# use cover difference to account for reverse causality
mod_diff_fix_loc_yr <- plm(PercDiffCovered ~ Lag1Treated, data = hydr_init_lag1_dat,
                           index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
summary(mod_diff_fix_loc_yr)
# reduces treatment est. to -2.1
# is this significant when standard errors are corrected for heteroskedasticity and autocorrelation?
coeftest(mod_diff_fix_loc_yr, vcov = vcovHC, type = "HC1")
# yes, highly significant

# initial cover increases treatment frequency and decreases difference
# how does including change estimate?
mod_init_diff_fix_loc_yr <- plm(PercDiffCovered ~ InitPercCovered_c * Lag1Treated, data = hydr_init_lag1_dat,
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
summary(mod_init_diff_fix_loc_yr)
# initial cover explains declines in difference
# treatment est. increases to 3
# treatment effect declines with higher initial cover
# would need to be ~25% above average PAC to see negative effect of treatment


#### model-fitting functions ####

# data filter function
dat_mod_filt <- function(treat_col, dat_in){
  
  dat_mod <- dat_in %>% 
    filter(!is.na(InitPercCovered)) %>%
    filter(!is.na(Lag1Treated) & !is.na(Lag2Treated) & !is.na(Lag3Treated) & !is.na(Lag4Treated) & !is.na(Lag5Treated) & !is.na(Lag6Treated)) %>%
    mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
           Treated = !!sym(treat_col),
           InitPercCovered_c = InitPercCovered - mean(InitPercCovered))
  
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
  
  # prop covered with initial PAC
  mod1 <- plm(PropCoveredLogit ~ InitPercCovered_c * Treated + SurveyorExperience_s, data = dat_mod1,
              index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
  mod2 <- update(mod1, data = dat_mod2)
  mod3 <- update(mod1, data = dat_mod3)
  mod4 <- update(mod1, data = dat_mod4)
  mod5 <- update(mod1, data = dat_mod5)
  mod6 <- update(mod1, data = dat_mod6)

  # perc diff with initial PAC
  mod7 <- plm(PercDiffCovered ~ InitPercCovered_c * Treated + SurveyorExperience_s, data = dat_mod1,
              index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
  mod8 <- update(mod7, data = dat_mod2)
  mod9 <- update(mod7, data = dat_mod3)
  mod10 <- update(mod7, data = dat_mod4)
  mod11 <- update(mod7, data = dat_mod5)
  mod12 <- update(mod7, data = dat_mod6)
  
  # prop covered without initial PAC
  mod1b <- plm(PropCoveredLogit ~ Treated + SurveyorExperience_s, data = dat_mod1,
               index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
  mod2b <- update(mod1b, data = dat_mod2)
  mod3b <- update(mod1b, data = dat_mod3)
  mod4b <- update(mod1b, data = dat_mod4)
  mod5b <- update(mod1b, data = dat_mod5)
  mod6b <- update(mod1b, data = dat_mod6)
  
  # perc diff without initial PAC
  mod7b <- plm(PercDiffCovered ~ Treated + SurveyorExperience_s, data = dat_mod1,
               index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
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
hydr_mods <- mod_fit(hydr_dat)
wahy_mods <- mod_fit(wahy_dat)
wale_mods <- mod_fit(wale_dat)
torp_mods <- mod_fit(torp_dat)
cubu_mods <- mod_fit(cubu_dat)
pagr_mods <- mod_fit(pagr_dat)

# name models
names(hydr_mods) <- names(wahy_mods) <- names(wale_mods) <- names(torp_mods) <- names(cubu_mods) <- names(pagr_mods) <- rep(c("1", "2", "3", "4", "5", "6"), 4)


#### coefficient figures and tables ####

# rename coefficients
coef_names <- c("SurveyorExperience_s" = "Surveyor experience",
                "InitPercCovered_c:Treated" = "Management:PAC",
                "InitPercCovered_c" = "Initial PAC",
                "Treated" = "Management")

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
                           filename, plotname){
  
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
    title = plotname
  )
  
  ggsave(filename, comb_fig,
         device = "eps", width = 6.5, height = 3.5, units = "in")
  
}

# prop covered with initial PAC
panel_plot_fun(hydr_mods[1:6], wahy_mods[1:6], wale_mods[1:6],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_invasive_prop_covered_init_model.eps",
               "Effects on invasive plant PAC (log-odds)")
panel_plot_fun(cubu_mods[1:6], pagr_mods[1:6], torp_mods[1:6],
               "Cuban bulrush", "Para grass", "Torpedograss",
               "output/fwc_non_focal_invasive_prop_covered_init_model.eps",
               "Effects on invasive plant PAC (log-odds)")

# perc diff with initial PAC
panel_plot_fun(hydr_mods[7:12], wahy_mods[7:12], wale_mods[7:12],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_invasive_PAC_diff_init_model.eps",
               "Effects on annual difference in invasive plant PAC")
panel_plot_fun(cubu_mods[7:12], pagr_mods[7:12], torp_mods[7:12],
               "Cuban bulrush", "Para grass", "Torpedograss",
               "output/fwc_non_focal_invasive_PAC_diff_init_model.eps",
               "Effects on annual difference in invasive plant PAC")

# prop covered without initial PAC
panel_plot_fun(hydr_mods[13:18], wahy_mods[13:18], wale_mods[13:18],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_invasive_prop_covered_model.eps",
               "Effects on invasive plant PAC (log-odds)")
panel_plot_fun(cubu_mods[13:18], pagr_mods[13:18], torp_mods[13:18],
               "Cuban bulrush", "Para grass", "Torpedograss",
               "output/fwc_non_focal_invasive_prop_covered_model.eps",
               "Effects on invasive plant PAC (log-odds)")

# perc diff without initial PAC
panel_plot_fun(hydr_mods[19:24], wahy_mods[19:24], wale_mods[19:24],
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_invasive_PAC_diff_model.eps",
               "Effects on annual difference in invasive plant PAC")
panel_plot_fun(cubu_mods[19:24], pagr_mods[19:24], torp_mods[19:24],
               "Cuban bulrush", "Para grass", "Torpedograss",
               "output/fwc_non_focal_invasive_PAC_diff_model.eps",
               "Effects on annual difference in invasive plant PAC")


#### finalize models ####

# SE with heteroskedasticity and autocorrelation
coeftest(hydr_mods[[19]], vcov = vcovHC, type = "HC1") # treatment sig, surveyor not
coeftest(wahy_mods[[19]], vcov = vcovHC, type = "HC1") # treatment sig, surveyor not
coeftest(wale_mods[[19]], vcov = vcovHC, type = "HC1") # neither sig
coeftest(cubu_mods[[19]], vcov = vcovHC, type = "HC1") # neither sig
coeftest(pagr_mods[[19]], vcov = vcovHC, type = "HC1") # neither sig
coeftest(torp_mods[[19]], vcov = vcovHC, type = "HC1") # neither sig

# don't need surveyor experience
# only using year1 lag
hydr_dat1 <- hydr_dat %>% filter(!is.na(Lag1Treated)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
wahy_dat1 <- wahy_dat %>% filter(!is.na(Lag1Treated)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
wale_dat1 <- wale_dat %>% filter(!is.na(Lag1Treated)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
cubu_dat1 <- cubu_dat %>% filter(!is.na(Lag1Treated)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
torp_dat1 <- torp_dat %>% filter(!is.na(Lag1Treated)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))
pagr_dat1 <- pagr_dat %>% filter(!is.na(Lag1Treated)) %>%
  pdata.frame(index = c("PermanentID", "GSYear"))

# fit models
hydr_mod <- plm(PercDiffCovered ~ Lag1Treated, data = hydr_dat1,
                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
wahy_mod <- update(hydr_mod, data = wahy_dat1)
wale_mod <- update(hydr_mod, data = wale_dat1)
cubu_mod <- update(hydr_mod, data = cubu_dat1)
torp_mod <- update(hydr_mod, data = torp_dat1)
pagr_mod <- update(hydr_mod, data = pagr_dat1)

# SE with heteroskedasticity and autocorrelation
coeftest(hydr_mod, vcov = vcovHC, type = "HC1") # sig
coeftest(wahy_mod, vcov = vcovHC, type = "HC1") # sig
coeftest(wale_mod, vcov = vcovHC, type = "HC1") # marginal
coeftest(cubu_mod, vcov = vcovHC, type = "HC1") # not
coeftest(pagr_mod, vcov = vcovHC, type = "HC1") # not
coeftest(torp_mod, vcov = vcovHC, type = "HC1") # not

# adjusted R2
summary(hydr_mod) # -0.06
summary(wahy_mod) # -0.06
summary(wale_mod) # -0.06
summary(cubu_mod) # -0.2
summary(pagr_mod) # -0.2
summary(torp_mod) # -0.07

# add fitted values
hydr_dat2 <- hydr_dat1 %>%
  mutate(Fitted = as.numeric(hydr_dat1$PercDiffCovered - residuals(hydr_mod)))

# start here
# repeate above with other datasets
# make figure like this:
# https://github.com/ycroissant/plm/issues/16


# fitted vs. observed
ggplot(hydr_dat2, aes(x = Fitted, y = PercDiffCovered)) + geom_point()
plot(as.numeric(hydr_dat1$PercDiffCovered - residuals(hydr_mod)), hydr_dat1$PercDiffCovered)
plot(as.numeric(wahy_dat1$PercDiffCovered - residuals(wahy_mod)), wahy_dat1$PercDiffCovered)
plot(as.numeric(wale_dat1$PercDiffCovered - residuals(wale_mod)), wale_dat1$PercDiffCovered)
plot(as.numeric(cubu_dat1$PercDiffCovered - residuals(cubu_mod)), cubu_dat1$PercDiffCovered)
plot(as.numeric(pagr_dat1$PercDiffCovered - residuals(pagr_mod)), pagr_dat1$PercDiffCovered)
plot(as.numeric(torp_dat1$PercDiffCovered - residuals(torp_mod)), torp_dat1$PercDiffCovered)
# poor fits


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
    tibble(PermanentID = names(fixef(mod))) %>% 
    mutate(FixedEffects = fixef(mod)) %>%
    expand_grid(Lag1Treated = seq(0, 1, length.out = 4)) %>%
    mutate(PercDiffCovered = FixedEffects + Lag1Treated * as.numeric(coef(mod)),
           CommonName = unique(dat_in$CommonName),
           PercChangeCovered = 100 * PercDiffCovered / mean(dat_in$InitPercCovered),
           Treated = Lag1Treated * 3)

  return(dat_out)
  
}

# combine datasets
foc_raw_dat <- hydr_dat1 %>%
  full_join(wahy_dat1) %>%
  full_join(wale_dat1) %>%
  group_by(CommonName) %>%
  mutate(PercChangeCovered = 100 * PercDiffCovered / mean(InitPercCovered),
         Treated = Lag1Treated * 3) %>%
  ungroup()

foc_pred_dat <- pred_dat(hydr_dat1, hydr_mod) %>%
  full_join(pred_dat(wahy_dat1, wahy_mod)) %>%
  full_join(pred_dat(wale_dat1, wale_mod))

foc_fit_dat <- hydr_dat1 %>%
  mutate(Fitted = as.numeric(PercDiffCovered - residuals(hydr_mod))) %>%
  full_join(wahy_dat1 %>%
              mutate(Fitted = as.numeric(PercDiffCovered - residuals(wahy_mod)))) %>%
  full_join(wale_dat1 %>%
              mutate(Fitted = as.numeric(PercDiffCovered - residuals(wale_mod)))) %>%
  mutate(Treated = Lag1Treated * 3) 

ggplot(foc_fit_dat, aes(x = Treated, y = Fitted)) +
  geom_hline(yintercept = 0, size = 0.25) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0,
               aes(color = CommonName)) +
  stat_summary(geom = "point", fun = "mean", size = 2,
               aes(color = CommonName)) +
  facet_wrap(~ CommonName, scales = "free") +
  labs(x = "Years managed (out of 3)", y = "Percent change in PAC\n(relative to avg. invasion)") +
  def_theme_paper +
  theme(legend.position = "none")

# figures
foc_pred_fig <- ggplot(foc_raw_dat, aes(x = Treated, y = PercChangeCovered)) +
  geom_hline(yintercept = 0, size = 0.25) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0,
               aes(color = CommonName)) +
  stat_summary(geom = "point", fun = "mean", size = 2,
               aes(color = CommonName)) +
  stat_summary(data = foc_pred_dat, aes(fill = CommonName),
               geom = "ribbon", fun.data = "mean_cl_boot", alpha = 0.5) +
  stat_summary(data = foc_pred_dat, aes(color = CommonName),
               geom = "line", fun = "mean", size = 0.5) +
  facet_wrap(~ CommonName, scales = "free") +
  labs(x = "Years managed (out of 3)", y = "Percent change in PAC\n(relative to avg. invasion)") +
  def_theme_paper +
  theme(legend.position = "none")

# save
ggsave("output/fwc_focal_invasive_PAC_change_treatment_prediction.png", foc_pred_fig,
       device = "png", width = 6.5, height = 2.5, units = "in")

# data tables
foc_sum <- foc_raw_dat %>%
  group_by(CommonName) %>%
  summarize(InitPercCovered = mean(InitPercCovered)) %>%
  ungroup() %>%
  full_join(foc_pred_dat %>%
              group_by(CommonName, Treated) %>%
              summarize(PercChangeCovered = mean(PercChangeCovered)) %>%
              ungroup() %>%
              mutate(Treated = fct_recode(as.factor(Treated),
                                          "None" = "0",
                                          "One" = "1",
                                          "Two" = "2",
                                          "Three" = "3")) %>%
              pivot_wider(names_from = Treated,
                          values_from = PercChangeCovered))

# save data table
write_csv(foc_sum, "output/fwc_focal_invasive_PAC_change_treatment_prediction.csv")


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
