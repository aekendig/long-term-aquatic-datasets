#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(modelsummary) # modelplot
library(patchwork) # combining figures
library(plm) # panel data models
library(sandwich) # vcovHC
library(lmtest) # coeftest
library(pals) # color palettes
library(janitor)

# figure settings
source("code/settings/figure_settings.R")

# import data
inv_plant <- read_csv("intermediate-data/FWC_only_invasive_plant_formatted.csv")
inv_ctrl <- read_csv("intermediate-data/FWC_invasive_control_formatted.csv")

# color palette
pal <- c("#000000", "#56B4E9")


#### edit data ####

# change column name of ctrl
inv_ctrl2 <- inv_ctrl %>%
  rename(TaxonName = Species)

# ctrl data range
range(inv_ctrl2$GSYear) # 1998-2020
range(inv_plant$GSYear) # 1982-2019

# combine datasets
inv_dat <- inner_join(inv_plant, inv_ctrl2)

# species presence
pres_tax <- inv_dat %>%
  group_by(AreaOfInterestID, TaxonName) %>%
  summarize(Present = sum(EstAreaCoveredRaw_ha > 0)) %>%
  ungroup() 

# lakes that have the species detected
pres_tax_inv <- pres_tax %>%
  filter(Present > 0)

# uninvaded lakes
pres_tax_uninv <- pres_tax %>%
  filter(Present == 0)

write_csv(pres_tax_uninv, "output/fwc_uninvaded_AOI.csv")

# require non-NA previous year
# require presence in at least one year
inv_dat2 <- inv_dat %>%
  filter(!is.na(PrevPercCovered)) %>%
  inner_join(pres_tax_inv %>%
               select(AreaOfInterestID, TaxonName)) %>%
  mutate(AreaOfInterestID = as.factor(AreaOfInterestID),
         Treatment = if_else(Treated == 0, "Not managed", "Managed") %>%
           fct_relevel("Not managed"))

# check data availability
inv_dat2 %>%
  filter(PercCovered > 0) %>%
  ggplot(aes(x = GSYear, y = PercCovered, color = AreaOfInterestID)) +
  geom_vline(xintercept = 2013) +
  geom_point() +
  facet_wrap(~ CommonName, scales = "free_y") +
  theme(legend.position = "none")
# Cuban bulrush is missing a lot of data before 2013

# species
tax_spp <- sort(unique(inv_dat2$TaxonName))

# loop through taxa
pdf("output/invasive_plant_time_series_by_taxon.pdf")

for(i in 1:length(tax_spp)){
  
  # subset data
  subdat <- inv_dat2 %>% filter(TaxonName == tax_spp[i])
  subdat_ctrl <- subdat %>% filter(Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = PercCovered, 
                           color = AreaOfInterestID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          labs(x = "Year", y = "Percent area covered", title = tax_spp[i]) +
          def_theme_paper +
          theme(strip.text = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()

# export
write_csv(inv_dat2, "intermediate-data/FWC_invasive_plant_analysis_formatted.csv")

inv_dat2 %>%
  distinct(AreaOfInterestID, PermanentID) %>%
  mutate(plant_management = 1) %>%
  write_csv("intermediate-data/FWC_invasive_plant_analysis_waterbodies.csv")

# split by species
hydr_dat <- filter(inv_dat2, CommonName == "Hydrilla")
wahy_dat <- filter(inv_dat2, CommonName == "Water hyacinth")
wale_dat <- filter(inv_dat2, CommonName == "Water lettuce")
torp_dat <- filter(inv_dat2, CommonName == "Torpedograss")
cubu_dat <- filter(inv_dat2, CommonName == "Cuban bulrush")
pagr_dat <- filter(inv_dat2, CommonName == "Para grass")


#### initial visualizations ####

# distributions
ggplot(inv_dat2, aes(x = PercCovLogRatio)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_dat2, aes(x = RecentTreatment)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# many points have never been treated

# treated and change in prop
ggplot(inv_dat2, aes(x = Treated, y = PercCovLogRatio)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ CommonName, scales = "free") 
# all decline with treatment

ggplot(inv_dat2, aes(x = Lag6Treated, y = PercCovLogRatio)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ CommonName, scales = "free") 
# less clear

ggplot(inv_dat2, aes(x = RecentTreatment, y = PercCovLogRatio)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ CommonName, scales = "free") 
# barely perceptible

# higher treatment with higher initial cover?
ggplot(inv_dat2, aes(x = PrevPercCovered, y = Treated)) +
  geom_point() +
  stat_smooth(method = "glm") +
  facet_wrap(~ CommonName, scales = "free") 
# yes 

# data availability for lags
inv_dat2 %>%
  rename(Lag1Treated = Treated) %>%
  select(AreaOfInterestID, GSYear, CommonName,
         Lag1Treated, Lag2Treated, Lag3Treated, Lag4Treated, Lag5Treated, Lag6Treated) %>%
  pivot_longer(cols = starts_with("Lag"),
               names_to = "Lag",
               values_to = "Treated",
               names_pattern = "Lag(.)Treated") %>%
  filter(!is.na(Treated)) %>%
  ggplot(aes(x = Lag)) +
  geom_bar() +
  facet_wrap(~ CommonName, scales = "free") 

# variation in growth rate
ggplot(inv_dat2, aes(x = AreaOfInterestID, y = PercCovLogRatio)) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ CommonName, scales = "free") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# years per waterbody
inv_dat2 %>%
  group_by(CommonName, AreaOfInterestID) %>%
  summarize(years = n_distinct(GSYear)) %>%
  ungroup() %>%
  ggplot(aes(x = years)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")
# include in supplement

# sample sizes
inv_dat2 %>%
  group_by(CommonName, Treated) %>%
  summarize(mean = mean(PercCovLogRatio),
            var = var(PercCovLogRatio),
            n = length(PercCovLogRatio))
# variance of hydrilla in treated years is very high


#### evaluate model structure ####

# sources: vignette("A_plmPackage"), https://www.princeton.edu/~otorres/Panel101R.pdf

# pooling model: same intercept and coefficients for all units and time
# within model: different fixed intercepts for each unit, same coefficients
# random effects model: different random intercepts for each unit, same coefficients
# variable coefficients model: different fixed or random intercept and coefficients for each

# same coefficients for each unit?
pooltest(PercCovLogRatio ~ Treated, data = hydr_dat, 
         index = c("AreaOfInterestID", "GSYear"), model = "within") # error
hydrw <- plm(PercCovLogRatio ~ Treated, data = hydr_dat, 
             index = c("AreaOfInterestID", "GSYear"), model = "within")
hydrv <- pvcm(PercCovLogRatio ~ Treated, data = hydr_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "within")
# fails here, some units don't have enough data to estimate coefficients

# variable intercept?
hydrp <- plm(PercCovLogRatio ~ Treated, data = hydr_dat, 
             index = c("AreaOfInterestID", "GSYear"), model = "pooling")
hydrl <- lm(PercCovLogRatio ~ Treated, data = hydr_dat)
summary(hydrw)
summary(hydrp)
summary(hydrl) # same as above
# need to account for variation among units in growth rate (see viz fig)

# random or fixed effects
# null hyp is that unit-level errors are uncorrelated with regressors
# not sig: use random effects, sig: use fixed effects
phtest(PercCovLogRatio ~ Treated, data = hydr_dat, 
       index = c("AreaOfInterestID", "GSYear")) # fixed
phtest(PercCovLogRatio ~ Treated, data = wahy_dat, 
       index = c("AreaOfInterestID", "GSYear")) # fixed
phtest(PercCovLogRatio ~ Treated, data = wale_dat, 
       index = c("AreaOfInterestID", "GSYear")) # random
phtest(PercCovLogRatio ~ Treated, data = cubu_dat, 
       index = c("AreaOfInterestID", "GSYear")) # fixed
phtest(PercCovLogRatio ~ Treated, data = pagr_dat, 
       index = c("AreaOfInterestID", "GSYear")) # random
phtest(PercCovLogRatio ~ Treated, data = torp_dat, 
       index = c("AreaOfInterestID", "GSYear")) # random

# individual and time effects for within?
plmtest(PercCovLogRatio ~ Treated, data = hydr_dat, 
        index = c("AreaOfInterestID", "GSYear"), model = "within",
        effect = "twoways", type = "bp") # sig
plmtest(PercCovLogRatio ~ Treated, data = wahy_dat, 
        index = c("AreaOfInterestID", "GSYear"), model = "within",
        effect = "twoways", type = "bp") # sig
plmtest(PercCovLogRatio ~ Treated, data = cubu_dat, 
        index = c("AreaOfInterestID", "GSYear"), model = "within",
        effect = "twoways", type = "bp") # sig

# individual and time effects for random?
waler1 <- plm(PercCovLogRatio ~ Treated, data = wale_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "twoways")
waler2 <- plm(PercCovLogRatio ~ Treated, data = wale_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "individual")
waler3 <- plm(PercCovLogRatio ~ Treated, data = wale_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "time")
waler4 <- lme4::lmer(PercCovLogRatio ~ Treated + (1|AreaOfInterestID), data = wale_dat)
summary(waler1)
summary(waler2)
summary(waler3)
summary(waler4) # same as 2, singular boundary error
# no variance among individuals, but there is some variance among years

pagrr1 <- plm(PercCovLogRatio ~ Treated, data = pagr_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "twoways")
pagrr2 <- plm(PercCovLogRatio ~ Treated, data = pagr_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "individual")
pagrr3 <- plm(PercCovLogRatio ~ Treated, data = pagr_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "time")
summary(pagrr1)
summary(pagrr2)
summary(pagrr3) # no individual variance, small time variance

torpr1 <- plm(PercCovLogRatio ~ Treated, data = torp_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "twoways")
torpr2 <- plm(PercCovLogRatio ~ Treated, data = torp_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "individual")
torpr3 <- plm(PercCovLogRatio ~ Treated, data = torp_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "time")
summary(torpr1)
summary(torpr2)
summary(torpr3) # no individual variance, small time variance


#### fit models ####

hydr_mod <- plm(PercCovLogRatio ~ Treated, data = hydr_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "within",
                effect = "twoways")
summary(hydr_mod)

wahy_mod <- plm(PercCovLogRatio ~ Treated, data = wahy_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "within",
                effect = "twoways")
summary(wahy_mod)

wale_mod <- plm(PercCovLogRatio ~ Treated, data = wale_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "random",
                effect = "time")
summary(wale_mod)

# try wale with fixed to see result
wale_mod2 <- plm(PercCovLogRatio ~ Treated, data = wale_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "within",
                effect = "twoways")
summary(wale_mod2) # not sig

cubu_mod <- plm(PercCovLogRatio ~ Treated, data = cubu_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "within",
                effect = "twoways")
summary(cubu_mod)

pagr_mod <- plm(PercCovLogRatio ~ Treated, data = pagr_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "random",
                effect = "time")
summary(pagr_mod)

torp_mod <- plm(PercCovLogRatio ~ Treated, data = torp_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "random",
                effect = "time")
summary(torp_mod)

# SE with heteroskedasticity and autocorrelation
(hydr_mod_p <- coeftest(hydr_mod, vcov = vcovHC, type = "HC3")) # sig
(wahy_mod_p <- coeftest(wahy_mod, vcov = vcovHC, type = "HC3")) # sig
(wale_mod_p <- coeftest(wale_mod, vcov = vcovHC, type = "HC3")) # not
(cubu_mod_p <- coeftest(cubu_mod, vcov = vcovHC, type = "HC3")) # sig
(pagr_mod_p <- coeftest(pagr_mod, vcov = vcovHC, type = "HC3")) # not
(torp_mod_p <- coeftest(torp_mod, vcov = vcovHC, type = "HC3")) # not

# result: negative effect of treatment for all species, sig for 3


#### figures ####

# trial figures
ggplot(inv_dat2, aes(x = Treatment, y = PercCovLogRatio)) +
  geom_point(size = 0.1, alpha = 0.5, 
             position = position_jitter(width = 0.1, height = 0)) +
  geom_boxplot() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_dat2, aes(x = Treatment, y = PercCovLogRatio)) +
  geom_point(size = 0.1, alpha = 0.5, 
             position = position_jitter(width = 0.1, height = 0)) +
  geom_violin(fill = NA) +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_dat2, aes(x = Treatment, y = PercCovLogRatio)) +
  geom_point(size = 0.1, alpha = 0.5, color = "yellow", 
             position = position_jitter(width = 0.1, height = 0)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ CommonName, scales = "free")

ggplot(inv_dat2, aes(x = Treatment, y = exp(PercCovLogRatio))) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ CommonName, scales = "free")

# figure function
treat_fig_fun <- function(dat_in, p_val, panel_title, file_name) {
  
  if(p_val < 0.001) {
    
    p_val <- formatC(p_val, format = "e", digits = 1)
    
  } else {
    
    p_val <- formatC(p_val, format = "g", digits = 1)
    
  }
  
  fig_p_val <- paste("p =", p_val)
  
  dat_sum <- dat_in %>%
    group_by(Treatment) %>%
    summarize(mean = mean(PercCovLogRatio),
              n = length(PercCovLogRatio),
              ymin = as.numeric(mean_cl_boot(PercCovLogRatio)[2]),
              ymax = as.numeric(mean_cl_boot(PercCovLogRatio)[3])) %>%
    ungroup() %>%
    mutate(samps = paste("n =", n),
           samps_y = min(ymin))
  
  raw_dat_fig <- ggplot(dat_in, aes(x = Treatment, y = PercCovLogRatio)) +
    geom_point(size = 0.1, alpha = 0.5, 
               position = position_jitter(width = 0.1, height = 0)) +
    labs(y = "Growth rate (log ratio)",
         title = panel_title) +
    def_theme_paper +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 9, color="black"))
  
  if(panel_title == "(A) hydrilla") {

    raw_dat_fig <- raw_dat_fig +
      scale_y_continuous(breaks = c(-5, 0, 5))

  }
  
  sum_dat_fig <- ggplot(dat_sum, aes(x = Treatment, y = mean)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
    geom_point(size = 2, shape = 21,
               aes(fill = Treatment)) +
    geom_text(aes(label = samps, y = samps_y),
              size = paper_text_size, vjust = 1.2) +
    annotate(geom = "text", label = fig_p_val, size = paper_text_size, 
             x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5) +
    labs(y = "Growth rate (log ratio)",
         title = panel_title) +
    scale_y_continuous(expand = expansion(mult = 0.1)) + 
    scale_fill_manual(values = pal, guide = "none") +
    def_theme_paper +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 9, color="black"))
  
  sum_dat_fig_pres <- ggplot(dat_sum, aes(x = Treatment, y = mean)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
    geom_point(size = 2) +
    geom_text(aes(label = samps, y = samps_y),
              size = 4, vjust = 1.2) +
    annotate(geom = "text", label = fig_p_val, size = 4, 
             x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5) +
    labs(y = "Growth rate (log ratio)") +
    def_theme +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 14, color="black"))
  
  ggsave(filename = file_name, plot = sum_dat_fig_pres, device = "jpeg",
         width = 5, height = 4)
  
  return(list(sum_dat_fig, raw_dat_fig))
  
}

# figures
hydr_fig <- treat_fig_fun(hydr_dat, hydr_mod_p[[4]], "(A) hydrilla",
                          "output/hydrilla_treatment_fig_presentation.jpg")
wahy_fig <- treat_fig_fun(wahy_dat, wahy_mod_p[[4]], "(B) water hyacinth",
                          "output/water_hyacinth_treatment_fig_presentation.jpg")
wale_fig <- treat_fig_fun(wale_dat, wale_mod_p[[8]], "(C) water lettuce",
                          "output/water_lettuce_treatment_fig_presentation.jpg")
cubu_fig <- treat_fig_fun(cubu_dat, cubu_mod_p[[4]], "(A) Cuban bulrush",
                          "output/cuban_bulrush_treatment_fig_presentation.jpg")
pagr_fig <- treat_fig_fun(pagr_dat, pagr_mod_p[[8]], "(B) para grass",
                          "output/paragrass_treatment_fig_presentation.jpg")
torp_fig <- treat_fig_fun(torp_dat, torp_mod_p[[8]], "(C) torpedograss",
                          "output/torpedograss_treatment_fig_presentation.jpg")

# combine
foc_figs <- hydr_fig[[1]] + theme(axis.title.y = element_blank(),
                                  axis.text.x = element_blank()) + 
  wahy_fig[[1]] + theme(axis.text.x = element_blank()) +
  wale_fig[[1]] + theme(axis.title.y = element_blank()) + 
  plot_layout(ncol = 1)
ggsave("output/fwc_focal_invasive_growth_treatment.png", foc_figs,
       device = "png", width = 3, height = 8, units = "in")

non_foc_figs <- cubu_fig[[1]] + theme(axis.title.y = element_blank(),
                                      axis.text.x = element_blank()) +  
  pagr_fig[[1]] + theme(axis.text.x = element_blank()) +
  torp_fig[[1]] + theme(axis.title.y = element_blank()) + 
  plot_layout(ncol = 1)
ggsave("output/fwc_non_focal_invasive_growth_treatment.png", non_foc_figs,
       device = "png", width = 3, height = 8, units = "in")

foc_raw_figs <- hydr_fig[[2]] + theme(axis.title.y = element_blank(),
                                      axis.text.x = element_blank()) + 
  wahy_fig[[2]] +  theme(axis.text.x = element_blank()) +
  wale_fig[[2]] + theme(axis.title.y = element_blank()) + 
  plot_layout(ncol = 1)
ggsave("output/fwc_focal_raw_invasive_growth_treatment.png", foc_raw_figs,
       device = "png", width = 3, height = 8, units = "in")

non_foc_raw_figs <- cubu_fig[[2]] + theme(axis.title.y = element_blank(),
                                          axis.text.x = element_blank()) + 
  pagr_fig[[2]] + theme(axis.text.x = element_blank()) +
  torp_fig[[2]] + theme(axis.title.y = element_blank()) + 
  plot_layout(ncol = 1)
ggsave("output/fwc_non_focal_raw_invasive_growth_treatment.png", non_foc_raw_figs,
       device = "png", width = 3, height = 8, units = "in")

#### tables ####

# combine summary
foc_mod_sum <- tibble(Species = "hydrilla",
                      Estimate = hydr_mod_p[1, 1],
                      SE = hydr_mod_p[1, 2],
                      t = hydr_mod_p[1, 3],
                      P = hydr_mod_p[1, 4],
                      R2 = r.squared(hydr_mod),
                      Waterbodies = n_distinct(hydr_dat$AreaOfInterestID),
                      Years = paste(range(count(hydr_dat, AreaOfInterestID)$n), collapse = "-"),
                      N = nrow(hydr_dat),
                      Model = "fixed effects") %>%
  full_join(tibble(Species = "water hyacinth",
                   Estimate = wahy_mod_p[1, 1],
                   SE = wahy_mod_p[1, 2],
                   t = wahy_mod_p[1, 3],
                   P = wahy_mod_p[1, 4],
                   R2 = r.squared(wahy_mod),
                   Waterbodies = n_distinct(wahy_dat$AreaOfInterestID),
                   Years = paste(range(count(wahy_dat, AreaOfInterestID)$n), collapse = "-"),
                   N = nrow(wahy_dat),
                   Model = "fixed effects")) %>%
  full_join(tibble(Species = "water lettuce",
                   Estimate = wale_mod_p[2, 1],
                   SE = wale_mod_p[2, 2],
                   t = wale_mod_p[2, 3],
                   P = wale_mod_p[2, 4],
                   R2 = r.squared(wale_mod),
                   Waterbodies = n_distinct(wale_dat$AreaOfInterestID),
                   Years = paste(range(count(wahy_dat, AreaOfInterestID)$n), collapse = "-"),
                   N = nrow(wale_dat),
                   Model = "random effects"))

non_foc_mod_sum <- tibble(Species = "Cuban bulrush",
                          Estimate = cubu_mod_p[1, 1],
                          SE = cubu_mod_p[1, 2],
                          t = cubu_mod_p[1, 3],
                          P = cubu_mod_p[1, 4],
                          R2 = r.squared(cubu_mod),
                          Waterbodies = n_distinct(cubu_dat$AreaOfInterestID),
                          Years = paste(range(count(cubu_dat, AreaOfInterestID)$n), collapse = "-"),
                          N = nrow(cubu_dat),
                          Model = "fixed effects") %>%
  full_join(tibble(Species = "para grass",
                   Estimate = pagr_mod_p[2, 1],
                   SE = pagr_mod_p[2, 2],
                   t = pagr_mod_p[2, 3],
                   P = pagr_mod_p[2, 4],
                   R2 = r.squared(pagr_mod),
                   Waterbodies = n_distinct(pagr_dat$AreaOfInterestID),
                   Years = paste(range(count(wahy_dat, AreaOfInterestID)$n), collapse = "-"),
                   N = nrow(pagr_dat),
                   Model = "random effects")) %>%
  full_join(tibble(Species = "torpedograss",
                   Estimate = torp_mod_p[2, 1],
                   SE = torp_mod_p[2, 2],
                   t = torp_mod_p[2, 3],
                   P = torp_mod_p[2, 4],
                   R2 = r.squared(torp_mod),
                   Waterbodies = n_distinct(torp_dat$AreaOfInterestID),
                   Years = paste(range(count(wahy_dat, AreaOfInterestID)$n), collapse = "-"),
                   N = nrow(torp_dat),
                   Model = "random effects"))

# export
write_csv(foc_mod_sum, "output/fwc_focal_treatment_model_summary.csv")
write_csv(non_foc_mod_sum, "output/fwc_non_focal_treatment_model_summary.csv")


#### values for text ####

# data tables
foc_sum <- tibble(CommonName = c("Hydrilla", "Water hyacinth", "Water lettuce"),
                  Incpt = c(mean(fixef(hydr_mod)), mean(fixef(wahy_mod)), 
                            as.numeric(wale_mod$coefficients["(Intercept)"])),
                  Beta = as.numeric(c(coef(hydr_mod), coef(wahy_mod), coef(wale_mod)[2]))) %>%
  mutate(IncptBeta = Incpt + Beta,
         NoTreat = 100 * (exp(Incpt) - 1),
         Treat = 100 * (exp(IncptBeta) - 1))

non_foc_sum <- tibble(CommonName = c("Cuban bulrush", "Para grass", "Torpedograss"),
                      Incpt = c(mean(fixef(cubu_mod)), as.numeric(pagr_mod$coefficients["(Intercept)"]), 
                                as.numeric(torp_mod$coefficients["(Intercept)"])),
                      Beta = as.numeric(c(coef(cubu_mod), coef(pagr_mod)[2], coef(torp_mod)[2]))) %>%
  mutate(IncptBeta = Incpt + Beta,
         NoTreat = 100 * (exp(Incpt) - 1),
         Treat = 100 * (exp(IncptBeta) - 1))

# save data table
write_csv(foc_sum, "output/fwc_focal_invasive_growth_treatment_prediction.csv")
write_csv(non_foc_sum, "output/fwc_non_focal_invasive_growth_treatment_prediction.csv")

