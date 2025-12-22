#### set-up ####

# load packages
library(tidyverse)
library(broom.mixed)
library(glmmTMB)
library(DHARMa)
library(GGally)
library(betareg)
library(janitor)
library(patchwork)
library(performance)
library(ggtext)

# figure settings
source("code/code-for-pub/figure_settings.R")

# import data
target_dat <- read_csv("intermediate-data/FWC_plant_management_target_analysis_formatted.csv")
target_taxa_dat <- read_csv("intermediate-data/FWC_plant_management_target_taxa_analysis_formatted.csv")
methods_dat <- read_csv("intermediate-data/FWC_plant_management_methods_analysis_formatted.csv")
methods_taxa_dat <- read_csv("intermediate-data/FWC_plant_management_methods_taxa_analysis_formatted.csv")
target_dat_new <- read_csv("intermediate-data/FWC_plant_management_new_target_analysis_formatted.csv")


#### examine target dataset ####

# select variables by waterbodies
# create variables for management efficacy analysis
target_dat_var <- target_dat %>%
  mutate(HydrTrt = if_else(HydrTrtFreq == 0, 0, 1),
         FloatTrt = if_else(FloatTrtFreq == 0, 0, 1),
         HydrProp = HydrCov / 100,
         FloatProp = FloatCov / 100) %>%
  distinct(AreaOfInterestID, HydrCovC, FloatCovC, 
           HydrTrtFreqC, HydrTrtAreaC, FloatTrtFreqC, 
           FloatTrtAreaC, OtherTrtFreqC, OtherTrtAreaC,
           HydrTrt, FloatTrt, HydrProp, FloatProp, LatitudeC)

# correlations among explanatory variables
target_dat_var %>%
  select(-c(AreaOfInterestID, HydrTrt, FloatTrt, HydrProp, FloatProp)) %>%
  ggpairs()
# correlations >= 0.4 and sig
# OtherTrtFreq and FloatTrtFreq: 0.6
# OtherTrtArea and FloatTrtArea: 0.4
# HydrTrtArea and HydrTrtFreq: 0.4


#### examine methods data ####

# select variables by waterbodies
methods_dat_var <- methods_dat %>%
  distinct(AreaOfInterestID, HydrPACc, FloatPACc, 
           TrtAreaConC, TrtAreaSysC, TrtAreaNonC, TrtAreaHerbC, 
           TrtFreqConC, TrtFreqSysC, TrtFreqNonC, TrtFreqHerbC, 
           TrtMonth, TrtMonthStd, LatitudeC)

# correlations among explanatory variables
methods_dat_var %>%
  select(-c(AreaOfInterestID, TrtMonth)) %>%
  ggpairs()
# correlations >= 0.4 and sig
# TrtAreaSys & TrtAreaCon: 0.4
# TrtFreqSys & TrtFreqCon: 0.6
# Herb with Con and Sys (expected)


#### survey date ####

# change in survey dates over time
day_mod1 <- glmmTMB(SurveyDay ~ Time + (1|AreaOfInterestID), data = target_dat,
                    family = poisson)
day_res1 <- simulateResiduals(day_mod1, n = 1000)
plot(day_res1) # need to refit

day_mod2 <- update(day_mod1, family = "nbinom2")
day_res2 <- simulateResiduals(day_mod2, n = 1000)
plot(day_res2) # better

summary(day_mod2)

write_csv(tidy(day_mod2), "output/survey_date_model_summary.csv")


#### management efficacy models ####

# hydrilla binary model
hydr_mod1 <- betareg(HydrProp ~ HydrTrt + LatitudeC, data = target_dat_var)
plot(hydr_mod1)
summary(hydr_mod1)

# save
write_csv(tidy(hydr_mod1), 
          "output/hydrilla_mgmt_efficacy_binary_model_summary.csv")

# hydrilla continuous model
hydr_mod2 <- betareg(HydrProp ~ HydrTrtFreqC + HydrTrtAreaC + LatitudeC, 
                     data = target_dat_var)
plot(hydr_mod2)
summary(hydr_mod2)

# save
write_csv(tidy(hydr_mod2), 
          "output/hydrilla_mgmt_efficacy_continuous_model_summary.csv")

# floating plant binary model
float_mod1 <- betareg(FloatProp ~ FloatTrt + LatitudeC, data = target_dat_var)
plot(float_mod1)
summary(float_mod1)

# save
write_csv(tidy(float_mod1), 
          "output/floating_mgmt_efficacy_binary_model_summary.csv")

# floating plant continuous model
float_mod2 <- betareg(FloatProp ~ FloatTrtFreqC + FloatTrtAreaC + LatitudeC, 
                      data = target_dat_var)
plot(float_mod2)
summary(float_mod2)

# save
write_csv(tidy(float_mod2), 
          "output/floating_mgmt_efficacy_continuous_model_summary.csv")


#### native richness model function ####

# combine model variables without time component
mod_dat_var <- target_dat_var %>%
  select(-LatitudeC) %>% # only use methods dat to show latitude results
  full_join(methods_dat_var)

# function for native plant model predictions
pred_fun <- function(variable, xaxis, model, dat){
  
  zero_dat <- tibble(Time = min(dat$Time),
                     HydrCovC = 0, 
                     FloatCovC = 0,
                     HydrTrtFreqC = 0, 
                     FloatTrtFreqC = 0, 
                     OtherTrtFreqC = 0,
                     HydrTrtAreaC = 0, 
                     FloatTrtAreaC = 0, 
                     OtherTrtAreaC = 0,
                     HydrPACc = 0, 
                     FloatPACc = 0,
                     TrtAreaConC = 0, 
                     TrtAreaSysC = 0, 
                     TrtAreaNonC = 0,
                     TrtFreqConC = 0, 
                     TrtFreqSysC = 0, 
                     TrtFreqNonC = 0, 
                     TrtMonthStd = 0,
                     LatitudeC = 0) 
  
  zero_dat_pred <- zero_dat %>%
    mutate(Pred = predict(model, newdata = ., type = "response", 
                          re.form = NA)) %>%
    pull(Pred)
  
  fin_dat = zero_dat %>%
    mutate(Time = max(dat$Time)) %>%
    select(-{{variable}}) %>%
    expand_grid(mod_dat_var %>%
                  filter(!is.na({{variable}}))  %>%
                  distinct({{variable}})) %>%
    mutate(Pred = predict(model, newdata = ., type = "response", re.form = NA),
           PredSE = predict(model, newdata = ., type = "response", re.form = NA,
                            se.fit = T)$se.fit,
           ZeroPred = zero_dat_pred,
           PredDiff = Pred - ZeroPred) %>%
    left_join(dat %>%
                distinct({{variable}}, {{xaxis}}))
  
}


#### native richness target model ####

# response variable
ggplot(target_dat, aes(x = NativeRichness)) +
  geom_density()
mean(target_dat$NativeRichness)
var(target_dat$NativeRichness)

# fit model
# negative binomial is very slow to fit -- use Poisson to evaluate
target_mod1 <- glmmTMB(NativeRichness ~ Time + Time:(HydrCovC + 
                                                       FloatCovC + 
                                                       HydrTrtFreqC +
                                                       FloatTrtFreqC + 
                                                       OtherTrtFreqC + 
                                                       HydrTrtAreaC +
                                                       FloatTrtAreaC + 
                                                       OtherTrtAreaC +
                                                       LatitudeC) +
                         (1|AreaOfInterestID),
                       data = target_dat, family = "poisson",
                       control = glmmTMBControl(optimizer = optim,
                                                optArgs = list(method = "BFGS")))

# assess model
target_res1 <- simulateResiduals(target_mod1, n = 1000)
plot(target_res1)
summary(target_mod1)

# evaluate correlated variables
check_collinearity(target_mod1)
# all low

# correlated: floating and other
target_mod1a <- update(target_mod1, .~. - Time:OtherTrtFreqC - 
                         Time:OtherTrtAreaC)
summary(target_mod1)
summary(target_mod1a)

target_mod1b <- update(target_mod1, .~. - Time:FloatTrtFreqC - 
                         Time:FloatTrtAreaC)
summary(target_mod1)
summary(target_mod1b)

# check hydrilla correlations
target_mod1c <- update(target_mod1, .~. -Time:HydrTrtFreqC)
summary(target_mod1)
summary(target_mod1c)

target_mod1d <- update(target_mod1, .~. -Time:HydrTrtAreaC)
summary(target_mod1)
summary(target_mod1d)

# fit with negative binomial
target_mod2 <- update(target_mod1, family = "nbinom2")

# save model
save(target_mod2, file = "output/native_richness_target_model.rda")
load("output/native_richness_target_model.rda")

# assess model
target_res2 <- simulateResiduals(target_mod2, n = 1000)
plot(target_res2)
summary(target_mod2)

# table
write_csv(tidy(target_mod2), "output/native_richness_target_model_summary.csv")

# data for figures
hyd_cov_fig_dat <- pred_fun(HydrCovC, HydrCov, target_mod2, target_dat)
hyd_frq_fig_dat <- pred_fun(HydrTrtFreqC, HydrTrtFreq, target_mod2, target_dat)
hyd_ext_fig_dat <- pred_fun(HydrTrtAreaC, HydrTrtArea, target_mod2, target_dat)
flt_cov_fig_dat <- pred_fun(FloatCovC, FloatCov, target_mod2, target_dat)
flt_frq_fig_dat <- pred_fun(FloatTrtFreqC, FloatTrtFreq, target_mod2, target_dat)
flt_ext_fig_dat <- pred_fun(FloatTrtAreaC, FloatTrtArea, target_mod2, target_dat)
oth_frq_fig_dat <- pred_fun(OtherTrtFreqC, OtherTrtFreq, target_mod2, target_dat)
oth_ext_fig_dat <- pred_fun(OtherTrtAreaC, OtherTrtArea, target_mod2, target_dat)

# figures
hyd_cov_fig <- ggplot(hyd_cov_fig_dat, aes(x = HydrCov, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line(linetype = "dashed") +
  labs(x = "Hydrilla cover",
       y = "Change in native richness") +
  def_theme_paper

hyd_frq_fig <- ggplot(hyd_frq_fig_dat, aes(x = HydrTrtFreq, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Hydrilla mgmt.\nfrequency",
       y = "Change in native richness") +
  def_theme_paper

hyd_ext_fig <- ggplot(hyd_ext_fig_dat, aes(x = HydrTrtArea, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line(linetype = "dashed") +
  labs(x = "Hydrilla mgmt.\nextent",
       y = "Change in native richness") +
  def_theme_paper

flt_cov_fig <- ggplot(flt_cov_fig_dat, aes(x = FloatCov, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Floating plant cover",
       y = "Change in native richness") +
  def_theme_paper

flt_frq_fig <- ggplot(flt_frq_fig_dat, aes(x = FloatTrtFreq, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Floating plant mgmt.\nfrequency",
       y = "Change in native richness") +
  def_theme_paper

flt_ext_fig <- ggplot(flt_ext_fig_dat, aes(x = FloatTrtArea, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line(linetype = "dashed") +
  labs(x = "Floating plant mgmt.\nextent",
       y = "Change in native richness") +
  def_theme_paper

oth_frq_fig <- ggplot(oth_frq_fig_dat, aes(x = OtherTrtFreq, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Other plant mgmt.\nfrequency",
       y = "Change in native richness") +
  def_theme_paper

oth_ext_fig <- ggplot(oth_ext_fig_dat, aes(x = OtherTrtArea, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line(linetype = "dashed") +
  labs(x = "Other plant mgmt.\nextent",
       y = "Change in native richness") +
  def_theme_paper


#### taxon-specific target models ####

# select native taxa
# add columns
target_taxa_dat2 <- target_taxa_dat %>%
  filter(Origin == "Native") %>%
  mutate(Habitat = tolower(Habitat) %>%
           fct_recode("emergent" = "emersed"))

# check taxa for all 1's or 0's
taxa_target_sum <- target_taxa_dat2 %>%
  group_by(TaxonName) %>%
  summarize(Detections = sum(Detected) / n(),
            .groups = "drop")

arrange(taxa_target_sum, Detections)
arrange(taxa_target_sum, desc(Detections))
# check the two extremes

# list of taxa
target_taxa <- sort(unique(target_taxa_dat2$TaxonName))

# set i
i <- 1

# subset data
target_taxa_sub <- filter(target_taxa_dat2, TaxonName == target_taxa[i])

# fit model
target_taxa_mod <- glmmTMB(Detected ~ Time + Time:(HydrCovC + FloatCovC +
                                                     HydrTrtFreqC + 
                                                     FloatTrtFreqC + 
                                                     OtherTrtFreqC + 
                                                     HydrTrtAreaC + 
                                                     FloatTrtAreaC + 
                                                     OtherTrtAreaC +
                                                     LatitudeC) + 
                             (1|AreaOfInterestID),
                           data = target_taxa_sub, family = "binomial")

# use this to manually cycle through taxa (some models had warnings below)
# i <- i + 1
# taxa with model convergence issues:

# # look at models
summary(target_taxa_mod)
target_taxa_res <- simulateResiduals(target_taxa_mod, n = 1000)
plot(target_taxa_res, title = target_taxa[i])

# save
save(target_taxa_mod, file = paste0("output/target_model_",
                                     str_to_lower(target_taxa[i]) %>%
                                       str_replace_all(" ", "_"),
                                     ".rda"))

# save results when i <- 1
target_taxa_coefs <- tidy(target_taxa_mod) %>%
  mutate(TaxonName = target_taxa[1]) %>%
  relocate(TaxonName)

# loop through taxa
pdf("output/taxa_specific_target_models.pdf")

for(i in target_taxa[-1]) {
  
  # subset data
  target_taxa_sub <- filter(target_taxa_dat2, TaxonName == i)
  
  # fit model
  target_taxa_mod <- glmmTMB(Detected ~ Time + Time:(HydrCovC + FloatCovC +
                                                       HydrTrtFreqC + 
                                                       FloatTrtFreqC + 
                                                       OtherTrtFreqC + 
                                                       HydrTrtAreaC + 
                                                       FloatTrtAreaC + 
                                                       OtherTrtAreaC +
                                                       LatitudeC) + 
                               (1|AreaOfInterestID),
                             data = target_taxa_sub, family = "binomial")
  
  # extract and plot residuals
  target_taxa_res <- simulateResiduals(target_taxa_mod, n = 1000)
  plot(target_taxa_res, title = i)
  
  # save model
  save(target_taxa_mod, file = paste0("output/target_model_",
                                       str_to_lower(i) %>%
                                         str_remove_all("\\.|\\(|\\)")%>%
                                         str_replace_all("/|, | ", "_"),
                                       ".rda"))
  
  target_taxa_coefs <- target_taxa_coefs %>%
    full_join(tidy(target_taxa_mod) %>%
                mutate(TaxonName = i) %>%
                relocate(TaxonName))
}

dev.off()

# save
write_csv(target_taxa_coefs, "output/target_model_taxa_coefficients.csv")
# import if needed
# target_taxa_coefs <- read_csv("output/target_model_taxa_coefficients.csv")

# correct p-values
# remove taxon that the estimates didn't work
target_taxa_coefs2 <- target_taxa_coefs %>%
  filter(str_detect(term, "Time:")) %>%
  mutate(q.value = p.adjust(p.value, method = "fdr"),
         lower = estimate - 1.96 * std.error,
         upper = estimate + 1.96 * std.error,
         odds_perc = 100 * (exp(estimate) - 1),
         lower_odds_perc = 100 * (exp(lower) - 1),
         upper_odds_perc = 100 * (exp(upper) - 1))

# save
write_csv(target_taxa_coefs2, 
          "output/target_model_taxa_interaction_coefficients.csv")
# import if needed
# target_taxa_coefs2 <- read_csv("output/target_model_taxa_interaction_coefficients.csv")

# significant interactions
target_taxa_coefs3 <- target_taxa_coefs2 %>%
  filter(q.value < 0.05 & term != "Time:LatitudeC") %>% 
  select(TaxonName, estimate, std.error, q.value, term) %>%
  left_join(target_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat) %>%
  mutate(Habitat = fct_relevel(Habitat, "floating"))

# significant coefficients table
target_taxa_tab <- target_taxa_coefs3 %>%
  mutate(across(.cols = c(estimate, std.error), 
                .fns = ~round_half_up(.x, 4) %>%
                  format(scientific = F)),
         q.value = round_half_up(q.value, 3),
         q.value = if_else(q.value == 0, "< 0.001", as.character(q.value)),
         response = paste0(estimate, " (", std.error, ")\n", q.value),
         term = str_remove(term, "Time:")) %>%
  select(Habitat, TaxonName, response, term) %>%
  pivot_wider(names_from = "term",
              values_from = "response") %>%
  arrange(Habitat, TaxonName) %>%
  relocate(c(HydrCovC, HydrTrtFreqC, HydrTrtAreaC, 
             FloatCovC, FloatTrtFreqC, FloatTrtAreaC, 
             OtherTrtFreqC, OtherTrtAreaC), 
           .after = TaxonName) %>%
  mutate(across(.cols = c(HydrCovC, HydrTrtFreqC, HydrTrtAreaC, 
                          FloatCovC, FloatTrtFreqC, FloatTrtAreaC, 
                          OtherTrtFreqC, OtherTrtAreaC),
                .fns = ~replace_na(.x, "")),
         TaxonName = str_replace_all(TaxonName, " ", "\n") %>%
           str_replace_all("\\/", "\\/\n"))

write_csv(target_taxa_tab, "output/target_model_taxa_interactions_table.csv")

# format estimates for figure
target_taxa_est <- target_taxa_coefs3 %>%
  mutate(term = fct_recode(term,
                           "hydrilla\ncover" = "Time:HydrCovC",
                           "hydrilla\nmgmt.\nfrequency" = "Time:HydrTrtFreqC",
                           "hydrilla\nmgmt.\nextent" = "Time:HydrTrtAreaC",
                           "floating\nplant\ncover" = "Time:FloatCovC",
                           "floating\nplant\nmgmt.\nfrequency" = "Time:FloatTrtFreqC",
                           "floating\nplant\nmgmt.\nextent" = "Time:FloatTrtAreaC",
                           "other\nplant\nmgmt.\nfrequency" = "Time:OtherTrtFreqC",
                           "other\nplant\nmgmt.\nextent" = "Time:OtherTrtAreaC"),
         sign = if_else(estimate > 0, "pos", "neg")) %>%
  count(term, sign, Habitat) %>%
  complete(term, sign, Habitat) %>%
  mutate(n = replace_na(n, 0),
         n = if_else(sign == "neg", -1 * n, n),
         term = fct_relevel(term, "hydrilla\ncover",
                            "hydrilla\nmgmt.\nfrequency",
                            "hydrilla\nmgmt.\nextent",
                            "floating\nplant\ncover",
                            "floating\nplant\nmgmt.\nfrequency",
                            "floating\nplant\nmgmt.\nextent",
                            "other\nplant\nmgmt.\nfrequency",
                            "other\nplant\nmgmt.\nextent"))

# counts for text
target_taxa_est %>%
  group_by(term, sign) %>% 
  summarize(sum = sum(n))

# figure
taxa_target_fig <- ggplot(target_taxa_est,
                          aes(x = term, y = n, fill = Habitat)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  scale_fill_manual(values = rev(col_pal), name = "Growth form") +
  labs(x = "Covariate interaction with time",
       y = "Number of native taxa (in direction of response)") +
  def_theme_paper +
  theme(legend.position = "inside",
        legend.position.inside = c(0.5, 0.78))

# combine figures
target_fig1 <- hyd_cov_fig + hyd_frq_fig + hyd_ext_fig + flt_cov_fig +
  plot_layout(nrow = 1, axis_titles = "collect_y")
target_fig2 <- flt_frq_fig + flt_ext_fig + oth_frq_fig + oth_ext_fig +
  plot_layout(nrow = 1, axis_titles = "collect_y")
target_fig <- target_fig1 / target_fig2 / taxa_target_fig +
  plot_layout(nrow = 3, heights = c(0.6, 0.6, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
                  title = "Figure 3") + 
  theme(plot.tag = element_text(size = 12, hjust = 0, vjust = 0),
        plot.title = element_text(size = 12))

ggsave("output/target_figure.png", target_fig, dpi = 600,
       width = 18, height = 20, units = "cm")


#### native richness methods model ####

# response variable
ggplot(methods_dat, aes(x = NativeRichness)) +
  geom_density()
mean(methods_dat$NativeRichness)
var(methods_dat$NativeRichness)

# fit model
methods_mod1 <- glmmTMB(NativeRichness ~ Time + Time:(HydrPACc + FloatPACc +
                                                        TrtAreaConC + 
                                                        TrtAreaSysC + 
                                                        TrtAreaNonC +
                                                        TrtFreqConC + 
                                                        TrtFreqSysC + 
                                                        TrtFreqNonC +
                                                        TrtMonthStd +
                                                        LatitudeC) +
                         (1|AreaOfInterestID),
                       data = methods_dat,
                       family = poisson,
                       control=glmmTMBControl(optimizer=optim,
                                              optArgs=list(method="BFGS")))

# assess model
methods_res1 <- simulateResiduals(methods_mod1, n = 1000)
plot(methods_res1)
summary(methods_mod1)

# evaluate correlated variables
check_collinearity(methods_mod1)
# all are low

# look at correlated variables
methods_mod1a <- update(methods_mod1, .~. - Time:TrtFreqSysC)
summary(methods_mod1)
summary(methods_mod1a)

methods_mod1b <- update(methods_mod1, .~. - Time:TrtFreqConC)
summary(methods_mod1)
summary(methods_mod1b)
# don't need to use herbicide variable

# fit with negative binomial
methods_mod2 <- update(methods_mod1, family = "nbinom1")
# nbinom2 couldn't converge

# save model
save(methods_mod2, file = "output/native_richness_methods_model.rda")
load("output/native_richness_methods_model.rda")

# assess model
methods_res2 <- simulateResiduals(methods_mod2, n = 1000)
plot(methods_res2)
summary(methods_mod2)

# table
write_csv(tidy(methods_mod2), 
          "output/native_richness_methods_model_summary.csv")

# data for figures
hyd_pac_fig_dat <- pred_fun(HydrPACc, HydrPAC, methods_mod2, methods_dat)
flt_pac_fig_dat <- pred_fun(FloatPACc, FloatPAC, methods_mod2, methods_dat)
cont_freq_fig_dat <- pred_fun(TrtFreqConC, TrtFreqCon, methods_mod2, methods_dat)
cont_ext_fig_dat <- pred_fun(TrtAreaConC, TrtAreaCon, methods_mod2, methods_dat)
sys_freq_fig_dat <- pred_fun(TrtFreqSysC, TrtFreqSys, methods_mod2, methods_dat)
sys_ext_fig_dat <- pred_fun(TrtAreaSysC, TrtAreaSys, methods_mod2, methods_dat)
month_fig_dat <- pred_fun(TrtMonthStd, TrtMonth, methods_mod2, methods_dat)
meth_lat_fig_dat <- pred_fun(LatitudeC, Latitude, methods_mod2, methods_dat)

# figures
hyd_pac_fig <- ggplot(hyd_pac_fig_dat, aes(x = HydrPAC, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line(linetype = "dashed") +
  labs(x = "Hydrilla cover",
       y = "Change in native richness") +
  def_theme_paper

flt_pac_fig <- ggplot(flt_pac_fig_dat, aes(x = FloatPAC, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Floating plant cover",
       y = "Change in native richness") +
  def_theme_paper

cont_freq_fig <- ggplot(cont_freq_fig_dat, aes(x = TrtFreqCon, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line(linetype = "dashed") +
  labs(x = "Contact herbicide\nfrequency",
       y = "Change in native richness") +
  def_theme_paper

cont_ext_fig <- ggplot(cont_ext_fig_dat, aes(x = TrtAreaCon, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Contact herbicide\nextent",
       y = "Change in native richness") +
  def_theme_paper

sys_freq_fig <- ggplot(sys_freq_fig_dat, aes(x = TrtFreqSys, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line(linetype = "dashed") +
  labs(x = "Systemic herbicide\nfrequency",
       y = "Change in native richness") +
  def_theme_paper

sys_ext_fig <- ggplot(sys_ext_fig_dat, aes(x = TrtAreaSys, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Systemic herbicide\nextent",
       y = "Change in native richness") +
  def_theme_paper

month_labs <- c(3, 6, 9, 12)

month_fig <- ggplot(month_fig_dat, aes(x = TrtMonth, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line(linetype = "dashed") +
  scale_x_continuous(breaks = month_labs, labels = month.abb[month_labs]) +
  labs(x = "Average mgmt.\nmonth",
       y = "Change in native richness") +
  def_theme_paper

meth_lat_fig <- ggplot(meth_lat_fig_dat, aes(x = Latitude, y = PredDiff)) +
  geom_ribbon(aes(ymin = PredDiff - PredSE, ymax = PredDiff + PredSE), 
              alpha = 0.3, color = NA) +
  geom_line() +
  labs(x = "Latitude",
       y = "Change in native richness") +
  def_theme_paper


#### taxon-specific methods models ####

# select native taxa
methods_taxa_dat2 <- methods_taxa_dat %>%
  filter(Origin == "Native") %>%
  mutate(Habitat = tolower(Habitat) %>%
           fct_recode("emergent" = "emersed"))

# check taxa for all 1's or 0's
taxa_methods_sum <- methods_taxa_dat2 %>%
  group_by(TaxonName) %>%
  summarize(Detections = sum(Detected) / n(),
            .groups = "drop")

arrange(taxa_methods_sum, Detections)
arrange(taxa_methods_sum, desc(Detections))
# check the two extremes

# list of taxa
methods_taxa <- sort(unique(methods_taxa_dat2$TaxonName))

# set i
i <- 1

# subset data
methods_taxa_sub <- filter(methods_taxa_dat2, TaxonName == methods_taxa[i])

# fit model
methods_taxa_mod <- glmmTMB(Detected ~ Time + Time:(HydrPACc + FloatPACc +
                                                      TrtAreaConC + 
                                                      TrtAreaSysC + 
                                                      TrtAreaNonC +
                                                      TrtFreqConC + 
                                                      TrtFreqSysC + 
                                                      TrtFreqNonC +
                                                      TrtMonthStd +
                                                      LatitudeC) +
                              (1|AreaOfInterestID),
                            data = methods_taxa_sub, family = "binomial")

# use this to manually cycle through taxa (some models had warnings below)
# i <- i + 1
# taxa with model convergence issues:
# 22: Fontinalis spp.
# 43: Najas marina
# 70: Sparganium americanum

# # look at models
summary(methods_taxa_mod)
methods_taxa_res <- simulateResiduals(methods_taxa_mod, n = 1000)
plot(methods_taxa_res, title = methods_taxa[i])

# save
save(methods_taxa_mod, file = paste0("output/methods_model_",
                                     str_to_lower(methods_taxa[i]) %>%
                                       str_replace_all(" ", "_"),
                                     ".rda"))

# save results when i <- 1
methods_taxa_coefs <- tidy(methods_taxa_mod) %>%
  mutate(TaxonName = methods_taxa[1]) %>%
  relocate(TaxonName)

# loop through taxa
pdf("output/taxa_specific_methods_models.pdf")

for(i in methods_taxa[-c(1, 22, 43, 70)]) {
  
  # subset data
  methods_taxa_sub <- filter(methods_taxa_dat2, TaxonName == i)
  
  # fit model
  methods_taxa_mod <- glmmTMB(Detected ~ Time + Time:(HydrPACc + FloatPACc +
                                                        TrtAreaConC + 
                                                        TrtAreaSysC + 
                                                        TrtAreaNonC +
                                                        TrtFreqConC + 
                                                        TrtFreqSysC + 
                                                        TrtFreqNonC +
                                                        TrtMonthStd +
                                                        LatitudeC) +
                                (1|AreaOfInterestID),
                              data = methods_taxa_sub, family = "binomial")
  
  # extract and plot residuals
  methods_taxa_res <- simulateResiduals(methods_taxa_mod, n = 1000)
  plot(methods_taxa_res, title = i)
  
  # save model
  save(methods_taxa_mod, file = paste0("output/methods_model_",
                                       str_to_lower(i) %>%
                                         str_remove_all("\\.|\\(|\\)")%>%
                                         str_replace_all("/|, | ", "_"),
                                       ".rda"))
  
  # save model coefficients
  methods_taxa_coefs <- methods_taxa_coefs %>%
    full_join(tidy(methods_taxa_mod) %>%
                mutate(TaxonName = i) %>%
                relocate(TaxonName))
}

dev.off()

# save
write_csv(methods_taxa_coefs, "output/methods_model_taxa_coefficients.csv")
# import if needed
# methods_taxa_coefs <- read_csv("output/methods_model_taxa_coefficients.csv")

# correct p-values for interaction terms
methods_taxa_coefs2 <- methods_taxa_coefs %>%
  filter(str_detect(term, "Time:")) %>%
  mutate(q.value = p.adjust(p.value, method = "fdr"),
         lower = estimate - 1.96 * std.error,
         upper = estimate + 1.96 * std.error,
         odds_perc = 100 * (exp(estimate) - 1),
         lower_odds_perc = 100 * (exp(lower) - 1),
         upper_odds_perc = 100 * (exp(upper) - 1))

# save
write_csv(methods_taxa_coefs2, "output/methods_model_taxa_interaction_coefficients.csv")
# import if needed
# methods_taxa_coefs2 <- read_csv("output/methods_model_taxa_interaction_coefficients.csv")

# significant interactions
methods_taxa_coefs3 <- methods_taxa_coefs2 %>%
  filter(q.value < 0.05 & str_detect(term, "Non") == F) %>% 
  select(TaxonName, estimate, std.error, q.value, term) %>%
  left_join(methods_taxa_dat2 %>%
              distinct(TaxonName, Habitat)) %>%
  relocate(Habitat) %>%
  mutate(Habitat = fct_relevel(Habitat, "floating"))

# significant coefficients table
methods_taxa_tab <- methods_taxa_coefs3 %>%
  mutate(across(.cols = c(estimate, std.error, q.value), .fns = ~round_half_up(.x, 3)),
         q.value = if_else(q.value == 0, "< 0.001", as.character(q.value)),
         response = paste0(estimate, " (", std.error, ")\n", q.value),
         term = str_remove(term, "Time:") %>% str_remove("Trt")) %>%
  select(Habitat, TaxonName, response, term) %>%
  pivot_wider(names_from = "term",
              values_from = "response") %>%
  arrange(Habitat, TaxonName) %>%
  relocate(c(HydrPACc, FloatPACc, FreqConC, AreaConC, FreqSysC, AreaSysC, 
             MonthStd, LatitudeC), .after = TaxonName) %>%
  mutate(across(.cols = c(HydrPACc, FloatPACc, FreqConC, AreaConC, FreqSysC, 
                          AreaSysC, MonthStd, LatitudeC),
                .fns = ~replace_na(.x, "")),
         TaxonName = str_replace_all(TaxonName, " ", "\n") %>%
           str_replace_all("\\/", "\\/\n"))

write_csv(methods_taxa_tab, "output/methods_model_taxa_interactions_table.csv")

# format estimates for figure
method_taxa_est <- methods_taxa_coefs3 %>%
  mutate(term = fct_recode(term,
                           "hydrilla\ncover" = "Time:HydrPACc",
                           "floating\nplant\ncover" = "Time:FloatPACc",
                           "contact\nherbicide\nextent" = "Time:TrtAreaConC",
                           "systemic\nherbicide\nextent" = "Time:TrtAreaSysC",
                           "contact\nherbicide\nfrequency" = "Time:TrtFreqConC",
                           "systemic\nherbicide\nfrequency" = "Time:TrtFreqSysC",
                           "average\nmanagement\nmonth" = "Time:TrtMonthStd",
                           "latitude" = "Time:LatitudeC"),
         sign = if_else(estimate > 0, "pos", "neg")) %>%
  count(term, sign, Habitat) %>%
  complete(term, sign, Habitat) %>%
  mutate(n = replace_na(n, 0),
         n = if_else(sign == "neg", -1 * n, n),
         term = fct_relevel(term, 
                            "hydrilla\ncover",
                            "floating\nplant\ncover",
                            "contact\nherbicide\nfrequency",
                            "contact\nherbicide\nextent",
                            "systemic\nherbicide\nfrequency",
                            "systemic\nherbicide\nextent",
                            "average\nmanagement\nmonth",
                            "latitude"))

# counts for text
method_taxa_est %>%
  group_by(term, sign) %>% 
  summarize(sum = sum(n))

# figure
taxa_method_fig <- ggplot(method_taxa_est,
                          aes(x = term, y = n, fill = Habitat)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  scale_fill_manual(values = rev(col_pal), name = "Growth form") +
  labs(x = "Covariate interaction with time",
       y = "Number of native taxa (in direction of response)") +
  def_theme_paper +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.78))

# combine figures
method_fig1 <- hyd_pac_fig + flt_pac_fig + cont_freq_fig + cont_ext_fig +
  plot_layout(nrow = 1, axis_titles = "collect_y")
method_fig2 <- sys_freq_fig + sys_ext_fig + month_fig + meth_lat_fig +
  plot_layout(nrow = 1, axis_titles = "collect_y")
method_fig <- method_fig1 / method_fig2 / taxa_method_fig +
  plot_layout(heights = c(0.6, 0.6, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
                  title = "Figure 5") + 
  theme(plot.tag = element_text(size = 12, hjust = 0, vjust = 0),
        plot.title = element_text(size = 12))

ggsave("output/method_figure.png", method_fig,
       width = 18, height = 20, units = "cm", dpi = 600)


#### variations of target model ####

# remove waterbodies with no management
# re-center variables
target_dat_mgmt <- target_dat %>%
  filter(!(HydrTrtFreq == 0 & FloatTrtFreq == 0 & OtherTrtFreq == 0)) %>%
  mutate(HydrCovC = HydrCov - mean(HydrCov),
         FloatCovC = FloatCov - mean(FloatCov),
         HydrTrtFreqC = HydrTrtFreq - mean(HydrTrtFreq),
         FloatTrtFreqC = FloatTrtFreq - mean(FloatTrtFreq),
         OtherTrtFreqC = OtherTrtFreq - mean(OtherTrtFreq),
         HydrTrtAreaC = HydrTrtArea - mean(HydrTrtArea),
         FloatTrtAreaC = FloatTrtArea - mean(FloatTrtArea),
         OtherTrtAreaC = OtherTrtArea - mean(OtherTrtArea),
         LatitudeC = Latitude - mean(Latitude))

# compare number of waterbodies
n_distinct(target_dat_mgmt$AreaOfInterestID)
n_distinct(target_dat$AreaOfInterestID)

# fit model
target_mod_mgmt <- update(target_mod2, data = target_dat_mgmt)
summary(target_mod_mgmt)
summary(target_mod2)

# save model
save(target_mod_mgmt, file = "output/native_richness_target_model_mgmt_only.rda")
write_csv(tidy(target_mod_mgmt), "output/target_model_all_mgmt_summary.csv")

# compare number of waterbodies target dataset from "new" management period 
# (used for methods model)
n_distinct(target_dat_new$AreaOfInterestID)
n_distinct(methods_dat$AreaOfInterestID)
# same

# compare sample size
nrow(target_dat_new)
nrow(methods_dat)
# same

# check for all waterbodies
anti_join(methods_dat %>% distinct(AreaOfInterestID, SurveyYear),
          target_dat_new %>% distinct(AreaOfInterestID, SurveyYear))

# fit model
target_mod_new <- glmmTMB(NativeRichness ~ Time + Time:(HydrPACc + 
                                                          FloatPACc + 
                                                          HydrTrtFreqC +
                                                          FloatTrtFreqC + 
                                                          OtherTrtFreqC + 
                                                          HydrTrtAreaC +
                                                          FloatTrtAreaC + 
                                                          OtherTrtAreaC +
                                                          LatitudeC) +
                            (1|AreaOfInterestID),
                          data = target_dat_new, family = "nbinom2",
                          control = glmmTMBControl(optimizer = optim,
                                                   optArgs = list(method = "BFGS"))) 
summary(target_mod_new) 
summary(target_mod2)

# save model
save(target_mod_new, file = "output/native_richness_target_model_newer_data.rda")
write_csv(tidy(target_mod_new), "output/target_model_new_mgmt_summary.csv")
