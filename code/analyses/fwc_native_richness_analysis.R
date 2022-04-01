#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(plotly)
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
# require previous year's richness
nat_plant2 <-  nat_plant %>%
  group_by(PermanentID, GSYear, Area_ha) %>%
  summarize(Richness = sum(Detected),
            PrevRichness = sum(PrevDetected)) %>%
  ungroup() %>%
  mutate(RichnessDiff = Richness - PrevRichness) %>%
  filter(!is.na(PrevRichness))

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

# taxa
inv_taxa <- sort(unique(nat_dat2$CommonName))

# loop through taxa
pdf("output/native_richness_continuous_time_series_by_taxon.pdf")

for(i in 1:length(inv_taxa)){
  
  # subset data
  subdat <- nat_dat2 %>% filter(CommonName == inv_taxa[i])
  subdat_ctrl <- subdat %>% filter(Lag1Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = Richness, color = PermanentID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          labs(x = "Year", y = "Native taxonomic richness", title = inv_taxa[i]) +
          def_theme_paper +
          theme(strip.text = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()

# split by species
hydr_dat <- filter(nat_dat2, CommonName == "Hydrilla")
wale_dat <- filter(nat_dat2, CommonName == "Water lettuce")
wahy_dat <- filter(nat_dat2, CommonName == "Water hyacinth")
torp_dat <- filter(nat_dat2, CommonName == "Torpedograss")
cubu_dat <- filter(nat_dat2, CommonName == "Cuban bulrush")
pagr_dat <- filter(nat_dat2, CommonName == "Para grass")


#### initial visualizations ####

# covariate correlations
nat_dat2 %>%
  select(CommonName, Lag1Treated, Lag1AvgPercCovered, MinSurveyorExperience, PrevRichness) %>%
  group_by(CommonName) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & corr >= 0.4) %>%
  data.frame()
# previous richness and surveyor experience for para grass

# richness diff distribution
ggplot(nat_dat2, aes(x = RichnessDiff)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free_y")

ggplot(nat_dat2, aes(x = Richness)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free_y")

# coefficients and richness
ggplot(nat_dat2, aes(x = Lag1AvgPercCovered, y = RichnessDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# generally close to zero

ggplot(nat_dat2, aes(x = Lag1AvgPercCovered, y = Richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# strong positive for hydrilla and torpedograss, negative for others

ggplot(nat_dat2, aes(x = Lag6Treated, y = RichnessDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# negligible

ggplot(nat_dat2, aes(x = Lag1Treated, y = Richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# positive for all except para grass

ggplot(nat_dat2, aes(x = PrevRichness, y = RichnessDiff)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# consistently negative, but not strong

ggplot(nat_dat2, aes(x = PrevRichness, y = Richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free")
# very strong positive


#### evaluate model structure ####

# Poisson or negative binomial
mean(wahy_dat$Richness)
var(wahy_dat$Richness)

# function to fit models for each species
mod_structure_fits <- function(dat_in){
  
  # create fixed effects data frame
  # make a global variable so that pglm can access it (probably a pglm bug)
  dat_fix <<- dat_in %>%
    mutate(PrevRichness_c = PrevRichness - mean(PrevRichness))  %>%
    ungroup() %>%
    pdata.frame(index = c("PermanentID", "GSYear"))
  # each waterbody is an individual
  
  # simple glm
  mod_glm <- glm(Richness ~ Lag1AvgPercCovered + Lag1Treated, family = poisson, data = dat_in)

  # random effects
  mod_ran_loc <- glmmTMB(Richness ~ Lag1AvgPercCovered + Lag1Treated + (1|PermanentID), family = poisson, data = dat_in)
  mod_ran_yr <- glmmTMB(Richness ~ Lag1AvgPercCovered + Lag1Treated + (1|GSYear), family = poisson, data = dat_in)
  mod_ran_loc_yr <- glmmTMB(Richness ~ Lag1AvgPercCovered + Lag1Treated + (1|PermanentID) + (1|GSYear), family = poisson, data = dat_in)

  # fixed effects
  mod_fix_loc <- pglm(Richness ~ Lag1AvgPercCovered + Lag1Treated, family = poisson, data = dat_fix,
                     model = "within")
  mod_fix_loc_yr <- pglm(Richness ~ Lag1AvgPercCovered + Lag1Treated, family = poisson, data = dat_fix,
                        model = "within", effect = "twoways")

  # use initial richness to account for reverse causality
  mod_init_fix_loc_yr <- pglm(Richness ~ PrevRichness_c + Lag1AvgPercCovered + Lag1Treated, family = poisson, data = dat_fix,
                              model = "within", effect = "twoways")

  # use richness difference to account for reverse causality
  mod_diff_fix_loc_yr <- plm(RichnessDiff ~ Lag1AvgPercCovered + Lag1Treated, data = dat_fix,
                             index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")

  # richness difference and initial richness
  mod_init_diff_fix_loc_yr <- plm(RichnessDiff ~ PrevRichness_c + Lag1AvgPercCovered + Lag1Treated, data = dat_fix,
                                  index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")

  # return list of models
  return(list(glm = mod_glm,
              ran_loc = mod_ran_loc,
              ran_yr = mod_ran_yr,
              ran_loc_yr = mod_ran_loc_yr,
              fix_loc = mod_fix_loc,
              fix_loc_yr = mod_fix_loc_yr,
              init_fix_loc_yr = mod_init_fix_loc_yr,
              diff_fix_loc_yr = mod_diff_fix_loc_yr,
              init_diff_fix_loc_yr = mod_init_diff_fix_loc_yr))

}

# filter for data with Lag1Treated
nat_dat2 %>% filter(is.na(Lag1Treated) | is.na(Lag1AvgPercCovered)) # none are missing

# fit models for each species
hydr_mod_struc <- mod_structure_fits(hydr_dat)
wahy_mod_struc <- mod_structure_fits(wahy_dat)
wale_mod_struc <- mod_structure_fits(wale_dat)

# remove dat_fix
rm("dat_fix")

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

write_csv(mod_comp, "output/fwc_native_richness_model_structure_comparison.csv")

# model comparison notes:
# global intercept -> hydrilla PAC increased and floating plant PAC decreased
  # all management increased
# random effect waterbody -> reduced all estimates except water hyacinth PAC
# random effect year -> only slight change
# fixed effect waterbody -> similar to random effect waterbody
# fixed effect year -> no change
# previous richness -> same direction for all compared to fixed waterbody, smaller magnitudes
# richness difference -> all estimates became negative
  # floating management effect is similar for water hyacinth and water lettuce models (expected)
# previous richness + richness difference -> positive effects of hydrilla and management
  # negative effects of floating plants and management
  # larger magnitudes than previous models
  # floating management effects differ between water hyacinth and lettuce models
  # previous richness estimates are consistent

# test fixed effects (seems like year isn't necessary)
# have to refit because data need to be accessible (not "dat_fix")
hydr_mod_diff_fix_loc_yr <- plm(RichnessDiff ~ Lag1AvgPercCovered + Lag1Treated, 
                                data = hydr_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(hydr_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(hydr_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

wahy_mod_diff_fix_loc_yr <- plm(RichnessDiff ~ Lag1AvgPercCovered + Lag1Treated, 
                                data = wahy_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(wahy_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(wahy_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

wale_mod_diff_fix_loc_yr <- plm(RichnessDiff ~ Lag1AvgPercCovered + Lag1Treated, 
                                data = wale_dat, 
                                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
plmtest(wale_mod_diff_fix_loc_yr, effect = "time", type = "bp") # sig
plmtest(wale_mod_diff_fix_loc_yr, effect = "individual", type = "bp") # sig

# use annual difference without initial PAC and with waterbody and year fixed effects



#### model-fitting functions ####

# data filter function
dat_mod_filt <- function(treat_col, inv_col, dat_in){
  
  dat_mod <- dat_in %>%
    filter(!is.na(Lag1Treated) & !is.na(Lag2Treated) & !is.na(Lag3Treated) & !is.na(Lag4Treated) & !is.na(Lag5Treated) & !is.na(Lag6Treated) & !is.na(Lag1AvgPercCovered) & !is.na(Lag2AvgPercCovered) & !is.na(Lag3AvgPercCovered) & !is.na(Lag4AvgPercCovered) & !is.na(Lag5AvgPercCovered) & !is.na(Lag6AvgPercCovered)) %>%
    mutate(SurveyorExperience_s = (MinSurveyorExperience - mean(MinSurveyorExperience)) / sd(MinSurveyorExperience),
           Treated = !!sym(treat_col),
           AvgPercCovered = !!sym(inv_col),
           PrevRichness_c = PrevRichness - mean(PrevRichness),
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
  # only one year available for Cuban bulrush
  if(foc_sp == "Cuban bulrush") {
    
    mod1 <- lm(RichnessDiff ~ AvgPercCovered_c + Treated + SurveyorExperience_s, data = dat_mod1)
    mod2 <- update(mod1, data = dat_mod2)
    mod3 <- update(mod1, data = dat_mod3)
    mod4 <- update(mod1, data = dat_mod4)
    mod5 <- update(mod1, data = dat_mod5)
    mod6 <- update(mod1, data = dat_mod6)
    
  } else {
    
    mod1 <- plm(RichnessDiff ~ AvgPercCovered_c + Treated + SurveyorExperience_s, data = dat_mod1,
                index = c("PermanentID", "GSYear"), model = "within", effect = "twoways")
    mod2 <- update(mod1, data = dat_mod2)
    mod3 <- update(mod1, data = dat_mod3)
    mod4 <- update(mod1, data = dat_mod4)
    mod5 <- update(mod1, data = dat_mod5)
    mod6 <- update(mod1, data = dat_mod6)
    
  }

  
  # output
  return(list(mod1, mod2, mod3, mod4, mod5, mod6))
  
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
names(hydr_mods) <- names(wahy_mods) <- names(wale_mods) <- names(torp_mods) <- names(cubu_mods) <- names(pagr_mods) <- c("1", "2", "3", "4", "5", "6")


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
    title = "Effects on annual difference in native richness")
  
  ggsave(filename, comb_fig,
         device = "eps", width = 6.5, height = 3.5, units = "in")
  
}

# figures
panel_plot_fun(hydr_mods, wahy_mods, wale_mods,
               "Hydrilla", "Water hyacinth", "Water lettuce",
               "output/fwc_focal_native_richness_diff_model.eps")
panel_plot_fun(cubu_mods, pagr_mods, torp_mods,
               "Cuban bulrush", "Para grass", "Torpedograss",
               "output/fwc_non_focal_native_richness_diff_model.eps")


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


#### model prediction figures ####

# extract fixed effects and treatment coefficients
# fitted values (as used in fwc_invasive_plant_analysis) include surveyor and PAC effects
# combine data
foc_fit_dat <- hydr_fit %>%
  full_join(wahy_fit) %>%
  full_join(wale_fit) %>%
  as_tibble() %>%
  select(CommonName, PermanentID, GSYear, Lag3Treated, RichnessDiff) %>%
  full_join(tibble(PermanentID = names(fixef(hydr_nat_rich_mod)),
                   fixef = as.numeric(fixef(hydr_nat_rich_mod)),
                   coef = as.numeric(coef(hydr_nat_rich_mod)[2]),
                   CommonName = "Hydrilla",
                   PanelName = "(A) hydrilla management") %>%
              full_join(tibble(PermanentID = names(fixef(wahy_nat_rich_mod)),
                               fixef = as.numeric(fixef(wahy_nat_rich_mod)),
                               coef = as.numeric(coef(wahy_nat_rich_mod)[2]),
                               CommonName = "Water hyacinth",
                               PanelName = "(B) water hyacinth management")) %>%
              full_join(tibble(PermanentID = names(fixef(wale_nat_rich_mod)),
                               fixef = as.numeric(fixef(wale_nat_rich_mod)),
                               coef = as.numeric(coef(wale_nat_rich_mod)[2]),
                               CommonName = "Water lettuce",
                               PanelName = "(C) water lettuce management"))) %>%
  mutate(Treated = Lag3Treated * 3,
         Fitted = fixef + coef * Treated) # replace Fitted with treatment-only effect

non_foc_fit_dat <- cubu_fit %>%
  full_join(pagr_fit) %>%
  full_join(torp_fit) %>%
  as_tibble() %>%
  select(CommonName, PermanentID, GSYear, Lag3Treated, RichnessDiff) %>%
  full_join(tibble(PermanentID = names(fixef(cubu_nat_rich_mod)),
                   fixef = as.numeric(fixef(cubu_nat_rich_mod)),
                   coef = as.numeric(coef(cubu_nat_rich_mod)[2]),
                   CommonName = "Cuban bulrush",
                   PanelName = "(A) Cuban bulrush management") %>%
              full_join(tibble(PermanentID = names(fixef(pagr_nat_rich_mod)),
                               fixef = as.numeric(fixef(pagr_nat_rich_mod)),
                               coef = as.numeric(coef(pagr_nat_rich_mod)[2]),
                               CommonName = "Para grass",
                               PanelName = "(B) para grass management")) %>%
              full_join(tibble(PermanentID = names(fixef(torp_nat_rich_mod)),
                               fixef = as.numeric(fixef(torp_nat_rich_mod)),
                               coef = as.numeric(coef(torp_nat_rich_mod)[2]),
                               CommonName = "Torpedograss",
                               PanelName = "(C) torpedograss management"))) %>%
  mutate(Treated = Lag3Treated * 3,
         Fitted = fixef + coef * Treated) # replace Fitted with treatment-only effect

# raw data (each waterbody represented)
ggplot(foc_fit_dat, aes(x = Treated, color = PermanentID)) +
  geom_point(aes(y = RichnessDiff), alpha = 0.3) +
  geom_line(aes(y = Fitted)) +
  facet_wrap(~ PanelName, scales = "free") +
  labs(x = "Years managed (out of 3)", y = "Annual difference in native richness") +
  def_theme_paper +
  theme(legend.position = "none")

ggplot(non_foc_fit_dat, aes(x = Treated, color = PermanentID)) +
  geom_point(aes(y = RichnessDiff), alpha = 0.3) +
  geom_line(aes(y = Fitted)) +
  facet_wrap(~ PanelName, scales = "free") +
  labs(x = "Years managed (out of 3)", y = "Annual difference in native richness") +
  def_theme_paper +
  theme(legend.position = "none")

# summarize raw data
foc_pred_fig <- ggplot(foc_fit_dat, aes(x = Treated)) +
  geom_point(aes(y = Fitted, color = PermanentID), size = 0.5, alpha = 0.05) +
  geom_line(aes(y = Fitted, color = PermanentID), alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0,
               aes(y = RichnessDiff)) +
  stat_summary(geom = "point", fun = "mean", size = 2,
               aes(y = RichnessDiff)) +
  facet_wrap(~ PanelName, scales = "free") +
  labs(x = "Years managed (out of 3)", y = "Annual difference in native richness") +
  scale_color_manual(values = rep(kelly(), 9)) +
  def_theme_paper +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, color = "black", hjust = 0))

non_foc_pred_fig <- ggplot(non_foc_fit_dat, aes(x = Treated)) +
  geom_point(aes(y = Fitted, color = PermanentID), size = 0.5, alpha = 0.05) +
  geom_line(aes(y = Fitted, color = PermanentID), alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0,
               aes(y = RichnessDiff)) +
  stat_summary(geom = "point", fun = "mean", size = 2,
               aes(y = RichnessDiff)) +
  facet_wrap(~ PanelName, scales = "free") +
  labs(x = "Years managed (out of 3)", y = "Annual difference in native richness") +
  scale_color_manual(values = rep(kelly(), 9)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3)) +
  def_theme_paper +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, color = "black", hjust = 0))

# save
ggsave("output/fwc_focal_invasive_native_richness_treatment_prediction.png", foc_pred_fig,
       device = "png", width = 6.5, height = 2.5, units = "in")

ggsave("output/fwc_non_focal_invasive_native_richness_treatment_prediction.png", non_foc_pred_fig,
       device = "png", width = 6.5, height = 2.5, units = "in")


#### values for text ####

# data tables
foc_sum <- tibble(CommonName = c("Hydrilla", "Water hyacinth", "Water lettuce"),
                   DiffNone = c(mean(fixef(hydr_nat_rich_mod)), mean(fixef(wahy_nat_rich_mod)), mean(fixef(wale_nat_rich_mod))),
                   Treat = as.numeric(c(coef(hydr_nat_rich_mod)[2], coef(wahy_nat_rich_mod)[2], coef(wale_nat_rich_mod)[2]))) %>%
  mutate(DiffThree = DiffNone + Treat,
         YearsNone = 1/DiffNone,
         YearsThree = 1/DiffThree)

non_foc_sum <- tibble(CommonName = c("Cuban bulrush", "Para grass", "Torpedograss"),
                      DiffNone = c(mean(fixef(cubu_nat_rich_mod)), mean(fixef(pagr_nat_rich_mod)), mean(fixef(torp_nat_rich_mod))),
                      Treat = as.numeric(c(coef(cubu_nat_rich_mod)[2], coef(pagr_nat_rich_mod)[2], coef(torp_nat_rich_mod)[2])),
                      PAC = as.numeric(c(coef(cubu_nat_rich_mod)[1], coef(pagr_nat_rich_mod)[1], coef(torp_nat_rich_mod)[1]))) %>%
  mutate(DiffThree = DiffNone + Treat,
         YearsNone = 1/DiffNone,
         YearsThree = 1/DiffThree,
         YearsPAC = 1/(DiffNone + PAC)) %>%
  left_join(cubu_dat3 %>%
              full_join(pagr_dat3) %>%
              full_join(torp_dat3) %>%
              group_by(CommonName) %>%
              summarize(InitPercCovered = mean(InitPercCovered)) %>%
              ungroup())

# save data table
write_csv(foc_sum, "output/fwc_focal_invasive_native_richness_treatment_prediction.csv")
write_csv(non_foc_sum, "output/fwc_non_focal_invasive_native_richness_treatment_prediction.csv")
