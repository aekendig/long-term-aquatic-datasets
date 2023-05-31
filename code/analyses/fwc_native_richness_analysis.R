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
library(janitor)
library(inspectdf) # inspect_cor

# figure settings
source("code/settings/figure_settings.R")

# import data
dat <- read_csv("intermediate-data/FWC_common_native_richness_invaded_data_formatted.csv")

# color palette
pal <- c("#000000", "#56B4E9")


#### edit data ####

# remove missing values
dat2 <- filter(dat, !is.na(LogRatioRichness)) %>%
  mutate(AreaOfInterestID = as.factor(AreaOfInterestID),
         Treatment = fct_relevel(Treatment, "Not managed"))

# taxa
inv_taxa <- sort(unique(dat2$CommonName))

# loop through taxa
pdf("output/native_richness_time_series_by_taxon.pdf")

for(i in inv_taxa){
  
  # subset data
  subdat <- dat2 %>% filter(CommonName == i)
  subdat_ctrl <- subdat %>% filter(Treated > 0)
  
  # make figure
  print(ggplot(subdat, aes(x = GSYear, y = Richness, color = AreaOfInterestID)) +
          geom_line() +
          geom_point(data = subdat_ctrl) + 
          labs(x = "Year", y = "Native richness", title = i) +
          def_theme_paper +
          theme(plot.title = element_text(hjust = 0.5, size = 8),
                legend.position = "none"))
  
}

dev.off()

# split data by invasive species
hydr_dat <- filter(dat2, CommonName == "Hydrilla")
wale_dat <- filter(dat2, CommonName == "Water lettuce")
wahy_dat <- filter(dat2, CommonName == "Water hyacinth") 
cubu_dat <- filter(dat2, CommonName == "Cuban bulrush")
pagr_dat <- filter(dat2, CommonName == "Para grass")
torp_dat <- filter(dat2, CommonName == "Torpedograss")


#### initial visualizations ####

# covariate correlations
dat2 %>%
  select(CommonName, PrevPercCovered, Treated) %>%
  group_by(CommonName) %>%
  inspect_cor() %>% 
  ungroup() %>%
  filter(p_value < 0.05 & abs(corr) >= 0.4) %>%
  data.frame()

# response distributions
ggplot(dat2, aes(x = LogRatioRichness)) +
  geom_histogram() +
  facet_wrap(~ CommonName, scales = "free")

ggplot(dat2, aes(x = PrevPercCovered, y = LogRatioRichness, 
                color = AreaOfInterestID)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, linewidth = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")

ggplot(dat2, aes(x = Treated, y = LogRatioRichness)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ CommonName, scales = "free") 

ggplot(dat2, aes(x = Lag6Treated, y = LogRatioRichness)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean") +
  facet_wrap(~ CommonName, scales = "free") 

# higher treatment with higher initial richness?
ggplot(dat2, aes(x = PrevRichness, y = Treated)) +
  geom_point() +
  stat_smooth(method = "glm") +
  facet_wrap(~ CommonName, scales = "free") 
# yes


#### evaluate model structure ####

# sources: vignette("A_plmPackage"), https://www.princeton.edu/~otorres/Panel101R.pdf

# pooling model: same intercept and coefficients for all units and time
# within model: different fixed intercepts for each unit, same coefficients
# random effects model: different random intercepts for each unit, same coefficients
# variable coefficients model: different fixed or random intercept and coefficients for each

# same coefficients for each unit?
pooltest(LogRatioRichness ~ Treated * PrevPercCovered, data = hydr_dat, 
         index = c("AreaOfInterestID", "GSYear"), model = "within") # error
hydrw <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = hydr_dat, 
             index = c("AreaOfInterestID", "GSYear"), model = "within")
hydrv <- pvcm(LogRatioRichness ~ Treated * PrevPercCovered, data = hydr_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "within")
# fails here, some units don't have enough data to estimate coefficients

# variable intercept?
hydrp <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = hydr_dat, 
             index = c("AreaOfInterestID", "GSYear"), model = "pooling")
summary(hydrw)
summary(hydrp)
# need to account for variation among units in growth rate (see viz fig)

# random or fixed effects
# null hyp is that unit-level errors are uncorrelated with regressors
# not sig: use random effects, sig: use fixed effects
phtest(LogRatioRichness ~ Treated * PrevPercCovered, data = hydr_dat, 
       index = c("AreaOfInterestID", "GSYear")) # random
phtest(LogRatioRichness ~ Treated * PrevPercCovered, data = wahy_dat, 
       index = c("AreaOfInterestID", "GSYear")) # random
phtest(LogRatioRichness ~ Treated * PrevPercCovered, data = wale_dat, 
       index = c("AreaOfInterestID", "GSYear")) # random
phtest(LogRatioRichness ~ Treated * PrevPercCovered, data = cubu_dat, 
       index = c("AreaOfInterestID", "GSYear")) # random
phtest(LogRatioRichness ~ Treated * PrevPercCovered, data = pagr_dat, 
       index = c("AreaOfInterestID", "GSYear")) # fixed
phtest(LogRatioRichness ~ Treated * PrevPercCovered, data = torp_dat, 
       index = c("AreaOfInterestID", "GSYear")) # fixed

# individual and time effects for within?
plmtest(LogRatioRichness ~ Treated * PrevPercCovered, data = pagr_dat, 
        index = c("AreaOfInterestID", "GSYear"), model = "within",
        effect = "twoways", type = "bp") # sig
plmtest(LogRatioRichness ~ Treated * PrevPercCovered, data = torp_dat, 
        index = c("AreaOfInterestID", "GSYear"), model = "within",
        effect = "twoways", type = "bp") # sig

# individual and time effects for random?
hydrr1 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = hydr_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "twoways")
hydrr2 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = hydr_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "individual")
hydrr3 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = hydr_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "time")
summary(hydrr1)
summary(hydrr2)
summary(hydrr3)
# no variance among individuals, but there is some variance among years

wahyr1 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = wahy_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "twoways")
wahyr2 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = wahy_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "individual")
wahyr3 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = wahy_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "time")
summary(wahyr1)
summary(wahyr2)
summary(wahyr3)
# no variance among individuals, but there is some variance among years

waler1 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = wale_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "twoways")
waler2 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = wale_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "individual")
waler3 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = wale_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "time")
summary(waler1)
summary(waler2)
summary(waler3)
# no variance among individuals, but there is some variance among years

cubur1 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = cubu_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "twoways")
cubur2 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = cubu_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "individual")
cubur3 <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = cubu_dat, 
              index = c("AreaOfInterestID", "GSYear"), model = "random",
              effect = "time")
summary(cubur1)
summary(cubur2)
summary(cubur3)
# no variance among individuals, but there is some variance among years


#### fit models ####

hydr_mod <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = hydr_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "random",
                effect = "time")
summary(hydr_mod)

wahy_mod <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = wahy_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "random",
                effect = "time")
summary(wahy_mod)

wale_mod <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = wale_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "random",
                effect = "time")
summary(wale_mod)

cubu_mod <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = cubu_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "random",
                effect = "time")
summary(cubu_mod)

pagr_mod <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = pagr_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "within",
                effect = "twoways")
summary(pagr_mod)

torp_mod <- plm(LogRatioRichness ~ Treated * PrevPercCovered, data = torp_dat, 
                index = c("AreaOfInterestID", "GSYear"), model = "within",
                effect = "twoways")
summary(torp_mod)

# SE with heteroskedasticity and autocorrelation
(hydr_mod_p <- coeftest(hydr_mod, vcov = vcovHC, type = "HC3")) 
(wahy_mod_p <- coeftest(wahy_mod, vcov = vcovHC, type = "HC3")) 
(wale_mod_p <- coeftest(wale_mod, vcov = vcovHC, type = "HC3")) 
(cubu_mod_p <- coeftest(cubu_mod, vcov = vcovHC, type = "HC3")) 
(pagr_mod_p <- coeftest(pagr_mod, vcov = vcovHC, type = "HC3")) 
(torp_mod_p <- coeftest(torp_mod, vcov = vcovHC, type = "HC3")) 


#### treatment figures ####

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
    summarize(mean = mean(LogRatioRichness),
              n = length(LogRatioRichness),
              ymin = as.numeric(mean_cl_boot(LogRatioRichness)[2]),
              ymax = as.numeric(mean_cl_boot(LogRatioRichness)[3])) %>%
    ungroup() %>%
    mutate(samps = paste("n =", n),
           samps_y = min(ymin))
  
  raw_dat_fig <- ggplot(dat_in, aes(x = Treatment, y = LogRatioRichness)) +
    geom_point(size = 0.1, alpha = 0.5, 
               position = position_jitter(width = 0.1, height = 0)) +
    labs(y = "Change in native richness (log ratio)",
         title = panel_title) +
    def_theme_paper +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 9, color="black"))
  
  sum_dat_fig <- ggplot(dat_sum, aes(x = Treatment, y = mean)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                  width = 0) +
    geom_point(size = 2, shape = 21,
               aes(fill = Treatment)) +
    geom_text(aes(label = samps, y = samps_y),
              size = paper_text_size, vjust = 1.2) +
    annotate(geom = "text", label = fig_p_val, size = paper_text_size, 
             x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5) +
    labs(y = "Change in native richness (log ratio)",
         title = panel_title) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult = 0.1)) + 
    def_theme_paper +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 9, color="black"),
          legend.position = "none")
  
  sum_dat_fig_pres <- ggplot(dat_sum, aes(x = Treatment, y = mean)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
    geom_point(size = 2) +
    geom_text(aes(label = samps, y = samps_y),
              size = 4, vjust = 1.2) +
    annotate(geom = "text", label = fig_p_val, size = 4, 
             x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5) +
    labs(y = "Change in native richness (log ratio)") +
    def_theme +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 14, color="black"))
  
  ggsave(filename = file_name, plot = sum_dat_fig_pres, device = "jpeg",
         width = 5, height = 4)
  
  return(list(sum_dat_fig, raw_dat_fig))
  
}

# figures
hydr_fig <- treat_fig_fun(hydr_dat, hydr_mod_p[2, 4], "(A) hydrilla",
                          "output/richness_hydrilla_treatment_fig_presentation.jpg")
wahy_fig <- treat_fig_fun(wahy_dat, wahy_mod_p[2, 4], "(B) water hyacinth",
                          "output/richness_water_hyacinth_treatment_fig_presentation.jpg")
wale_fig <- treat_fig_fun(wale_dat, wale_mod_p[2, 4], "(C) water lettuce",
                          "output/richness_water_lettuce_treatment_fig_presentation.jpg")
cubu_fig <- treat_fig_fun(cubu_dat, cubu_mod_p[2, 4], "(A) Cuban bulrush",
                          "output/richness_cuban_bulrush_treatment_fig_presentation.jpg")
pagr_fig <- treat_fig_fun(pagr_dat, pagr_mod_p[1, 4], "(B) paragrass",
                          "output/richness_paragrass_treatment_fig_presentation.jpg")
torp_fig <- treat_fig_fun(torp_dat, torp_mod_p[1, 4], "(C) torpedograss",
                          "output/richness_torpedograss_treatment_fig_presentation.jpg")

# combine
foc_figs <- hydr_fig[[1]] + wahy_fig[[1]] + wale_fig[[1]] + plot_layout(ncol = 1)
ggsave("output/fwc_focal_richness_treatment.png", foc_figs,
       device = "png", width = 3, height = 7, units = "in")

non_foc_figs <- cubu_fig[[1]] + pagr_fig[[1]] + torp_fig[[1]] + plot_layout(ncol = 1)
ggsave("output/fwc_non_focal_richness_treatment.png", non_foc_figs,
       device = "png", width = 3, height = 8, units = "in")

foc_raw_figs <- hydr_fig[[2]] + wahy_fig[[2]] + wale_fig[[2]] + plot_layout(ncol = 1)
ggsave("output/fwc_focal_raw_richness_treatment.png", foc_raw_figs,
       device = "png", width = 3, height = 7, units = "in")

non_foc_raw_figs <- cubu_fig[[2]] + pagr_fig[[2]] + torp_fig[[2]] + plot_layout(ncol = 1)
ggsave("output/fwc_non_focal_raw_richness_treatment.png", non_foc_raw_figs,
       device = "png", width = 3, height = 8, units = "in")


#### PAC figures ####

# figure function
PAC_fig_fun <- function(dat_in, p_vals, mod, panel_title, file_name) {
  
  if(p_vals[1] < 0.001) {
    
    p_vals1 <- formatC(p_vals[1], format = "e", digits = 1)
    
  } else {
    
    p_vals1 <- formatC(p_vals[1], format = "g", digits = 1)
    
  }
  
  if(p_vals[2] < 0.001) {
    
    p_vals2 <- formatC(p_vals[2], format = "e", digits = 1)
    
  } else {
    
    p_vals2 <- formatC(p_vals[2], format = "g", digits = 1)
    
  } 
  
  fig_p_vals <- paste0("PAC p = ", p_vals1, ", interaction p = ", p_vals2)

  dat_in2 <- dat_in %>%
    mutate(Pred = predict(mod, newdata = .),
           LogitPrevPercCov = car::logit(PrevPercCovered/100, adjust = 0.001))
  
  fig <- ggplot(dat_in2, aes(x = PrevPercCovered, y = LogRatioRichness)) +
    geom_point(alpha = 0.2, size = 0.5, aes(color = Treatment)) +
    geom_line(aes(y = Pred, color = Treatment)) +
    annotate(geom = "text", label = fig_p_vals, size = paper_text_size, 
             x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5) +
    scale_y_continuous(expand = expansion(mult = 0.1)) + 
    scale_color_manual(values = pal) +
    labs(x = "Invasive plant PAC") +
    def_theme_paper +
    theme(axis.title.y = element_blank(),
          legend.title = element_blank())
  
  fig_pres <- ggplot(dat_in2, aes(x = PrevPercCovered, y = LogRatioRichness)) +
    geom_point(alpha = 0.2, aes(color = Treatment)) +
    geom_line(aes(y = Pred, color = Treatment)) +
    annotate(geom = "text", label = fig_p_vals, size = paper_text_size, 
             x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5) +
    scale_color_manual(values = pal) +
    labs(x = "Invasive plant PAC",
         y = "Change in native richness (log ratio)") +
    def_theme +
    theme(legend.title = element_blank())
  
  ggsave(filename = file_name, plot = fig_pres, device = "jpeg",
         width = 5, height = 4)
  
  return(fig)
  
}

hydr_fig2 <- PAC_fig_fun(hydr_dat, hydr_mod_p[3:4, 4], hydr_mod, "(A) hydrilla",
                          "output/richness_hydrilla_PAC_fig_presentation.jpg")
wahy_fig2 <- PAC_fig_fun(hydr_dat, wahy_mod_p[3:4, 4], wahy_mod, "(B) water hyacinth",
                         "output/richness_water_hyacinth_PAC_fig_presentation.jpg")
wale_fig2 <- PAC_fig_fun(hydr_dat, wale_mod_p[3:4, 4], hydr_mod, "(C) water lettuce",
                         "output/richness_water_lettuce_PAC_fig_presentation.jpg")

cubu_fig2 <- PAC_fig_fun(cubu_dat, cubu_mod_p[3:4, 4], cubu_mod, "(A) Cuban bulrush",
                         "output/richness_cuban_bulrush_PAC_fig_presentation.jpg")
pagr_fig2 <- PAC_fig_fun(pagr_dat, pagr_mod_p[2:3, 4], pagr_mod, "(B) paragrass",
                         "output/richness_paragrass_PAC_fig_presentation.jpg")
torp_fig2 <- PAC_fig_fun(torp_dat, as.numeric(torp_mod_p[2:3, 4]), torp_mod, "(B) torpedograss",
                         "output/richness_torpedograss_PAC_fig_presentation.jpg")

# combine and save
foc_figs_PAC <- hydr_fig[[1]] + theme(axis.title.y = element_blank(),
                                      axis.text.x = element_blank()) + 
  hydr_fig2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  wahy_fig[[1]] + theme(axis.title.y = element_text(hjust = -0.01),
                        axis.text.x = element_blank()) +
  wahy_fig2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  wale_fig[[1]] + theme(axis.title.y = element_blank()) + 
  wale_fig2 + 
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, 0, 0, 0)) &
  guides(fill = "none")
ggsave("output/fwc_focal_PAC_richness.png", foc_figs_PAC,
       device = "png", width = 6, height = 7, units = "in")

non_foc_figs_PAC <- cubu_fig[[1]] + theme(axis.title.y = element_blank(),
                                      axis.text.x = element_blank()) + 
  cubu_fig2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  pagr_fig[[1]] + theme(axis.title.y = element_text(hjust = -0.01),
                        axis.text.x = element_blank()) +
  pagr_fig2 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  torp_fig[[1]] + theme(axis.title.y = element_blank()) + 
  torp_fig2 + 
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, 0, 0, 0)) &
  guides(fill = "none")
ggsave("output/fwc_non_focal_PAC_richness.png", non_foc_figs_PAC,
       device = "png", width = 6, height = 7, units = "in")


#### cover figures ####

cov_fig_fun <- function(dat_in, p_vals, mod, panel_title, file_name) {
  
  if(p_vals[1] < 0.001) {
    
    p_vals1 <- formatC(p_vals[1], format = "e", digits = 1)
    
  } else {
    
    p_vals1 <- formatC(p_vals[1], format = "g", digits = 1)
    
  }
  
  if(p_vals[2] < 0.001) {
    
    p_vals2 <- formatC(p_vals[2], format = "e", digits = 1)
    
  } else {
    
    p_vals2 <- formatC(p_vals[2], format = "g", digits = 1)
    
  } 
  
  fig_p_vals <- paste0("% cover p = ", p_vals1, ", interaction p = ", p_vals2)
  
  dat_in2 <- dat_in %>%
    mutate(Pred = predict(mod, newdata = .),
           LogitPrevPercCov = car::logit(PrevPercCovered/100, adjust = 0.001))
  
  fig <- ggplot(dat_in2, aes(x = PrevPercCovered, y = LogRatioRichness)) +
    geom_point(alpha = 0.2, size = 0.5, aes(color = Treatment)) +
    geom_line(aes(y = Pred, color = Treatment)) +
    annotate(geom = "text", label = fig_p_vals, size = paper_text_size, 
             x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5) +
    scale_y_continuous(expand = expansion(mult = 0.1)) + 
    scale_color_manual(values = pal) +
    labs(x = "Invasive plant cover (%)") +
    def_theme_paper +
    theme(axis.title.y = element_blank(),
          legend.title = element_blank())
  
  fig_pres <- ggplot(dat_in2, aes(x = PrevPercCovered, y = LogRatioRichness)) +
    geom_point(alpha = 0.2, aes(color = Treatment)) +
    geom_line(aes(y = Pred, color = Treatment)) +
    annotate(geom = "text", label = fig_p_vals, size = paper_text_size, 
             x = -Inf, y = Inf, hjust = -0.15, vjust = 1.5) +
    scale_color_manual(values = pal) +
    labs(x = "Invasive plant cover (%)",
         y = "Change in native richness (log ratio)") +
    def_theme +
    theme(legend.title = element_blank())
  
  ggsave(filename = file_name, plot = fig_pres, device = "jpeg",
         width = 5, height = 4)
  
  return(fig)
  
}

hydr_fig3 <- cov_fig_fun(hydr_dat, hydr_mod_p[3:4, 4], hydr_mod, "(A) hydrilla",
                         "output/richness_hydrilla_cov_fig_presentation.jpg")
wahy_fig3 <- cov_fig_fun(hydr_dat, wahy_mod_p[3:4, 4], wahy_mod, "(B) water hyacinth",
                         "output/richness_water_hyacinth_cov_fig_presentation.jpg")
wale_fig3 <- cov_fig_fun(hydr_dat, wale_mod_p[3:4, 4], hydr_mod, "(C) water lettuce",
                         "output/richness_water_lettuce_cov_fig_presentation.jpg")

cubu_fig3 <- cov_fig_fun(cubu_dat, cubu_mod_p[3:4, 4], cubu_mod, "(A) Cuban bulrush",
                         "output/richness_cuban_bulrush_cov_fig_presentation.jpg")
pagr_fig3 <- cov_fig_fun(pagr_dat, pagr_mod_p[2:3, 4], pagr_mod, "(B) paragrass",
                         "output/richness_paragrass_cov_fig_presentation.jpg")
torp_fig3 <- cov_fig_fun(torp_dat, as.numeric(torp_mod_p[2:3, 4]), torp_mod, "(B) torpedograss",
                         "output/richness_torpedograss_cov_fig_presentation.jpg")

# combine and save
foc_figs_cov <- hydr_fig[[1]] + theme(axis.title.y = element_blank(),
                                      axis.text.x = element_blank()) + 
  hydr_fig3 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  wahy_fig[[1]] + theme(axis.title.y = element_text(hjust = -0.01),
                        axis.text.x = element_blank()) +
  wahy_fig3 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  wale_fig[[1]] + theme(axis.title.y = element_blank()) + 
  wale_fig3 + 
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, 0, 0, 0)) &
  guides(fill = "none")
ggsave("output/fwc_focal_cov_richness.png", foc_figs_cov,
       device = "png", width = 6, height = 7, units = "in")

non_foc_figs_cov <- cubu_fig[[1]] + theme(axis.title.y = element_blank(),
                                          axis.text.x = element_blank()) + 
  cubu_fig3 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  pagr_fig[[1]] + theme(axis.title.y = element_text(hjust = -0.01),
                        axis.text.x = element_blank()) +
  pagr_fig3 + theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  torp_fig[[1]] + theme(axis.title.y = element_blank()) + 
  torp_fig3 + 
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.margin = margin(-10, 0, 0, 0)) &
  guides(fill = "none")
ggsave("output/fwc_non_focal_cov_richness.png", non_foc_figs_cov,
       device = "png", width = 6, height = 7, units = "in")


#### tables ####

hydr_mod_sum <- tibble(Coefficient = c("intercept", "hydrilla management", "hydrilla PAC", 
                                      "interaction", "error", "random: year"),
                      Estimate = c(hydr_mod_p[, 1], summary(hydr_mod)$ercomp$sigma2),
                      SE = c(hydr_mod_p[, 2], NA, NA),
                      t = c(hydr_mod_p[, 3], NA, NA),
                      P = c(hydr_mod_p[, 4], NA, NA),
                      R2 = summary(hydr_mod)$r.squared[2],
                      Waterbodies = n_distinct(hydr_dat$AreaOfInterestID),
                      Years = paste(range(count(hydr_dat, AreaOfInterestID)$n), collapse = "-"),
                      N = nrow(hydr_dat))

wale_mod_sum <- tibble(Coefficient = c("intercept", "water lettuce management", "water lettuce PAC", 
                                       "interaction", "error", "random: year"),
                       Estimate = c(wale_mod_p[, 1], summary(wale_mod)$ercomp$sigma2),
                       SE = c(wale_mod_p[, 2], NA, NA),
                       t = c(wale_mod_p[, 3], NA, NA),
                       P = c(wale_mod_p[, 4], NA, NA),
                       R2 = summary(wale_mod)$r.squared[2],
                       Waterbodies = n_distinct(wale_dat$AreaOfInterestID),
                       Years = paste(range(count(wale_dat, AreaOfInterestID)$n), collapse = "-"),
                       N = nrow(wale_dat))

wahy_mod_sum <- tibble(Coefficient = c("intercept", "water hyacinth management", "water hyacinth PAC", 
                                       "interaction", "error", "random: year"),
                       Estimate = c(wahy_mod_p[, 1], summary(wahy_mod)$ercomp$sigma2),
                       SE = c(wahy_mod_p[, 2], NA, NA),
                       t = c(wahy_mod_p[, 3], NA, NA),
                       P = c(wahy_mod_p[, 4], NA, NA),
                       R2 = summary(wahy_mod)$r.squared[2],
                       Waterbodies = n_distinct(wahy_dat$AreaOfInterestID),
                       Years = paste(range(count(wahy_dat, AreaOfInterestID)$n), collapse = "-"),
                       N = nrow(wahy_dat))

#### start here ####
# export tables and add to BAA
# extract values below for text


#### values for text ####

# data tables
foc_sum <- tibble(CommonName = c("Hydrilla", "Water hyacinth", "Water lettuce"),
                  Incpt = c(mean(fixef(hydr_mod)), mean(fixef(wahy_mod)), 
                            as.numeric(wale_mod$coefficients["(Intercept)"])),
                  Beta = as.numeric(c(coef(hydr_mod), coef(wahy_mod), coef(wale_mod)[2]))) %>%
  mutate(IncptBeta = Incpt + Beta,
         NoTreat = 100 * (exp(Incpt) - 1),
         Treat = 100 * (exp(IncptBeta) - 1))
