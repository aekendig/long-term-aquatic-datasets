#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(janitor)
library(cowplot)
library(grid)
library(gridExtra)
library(glmmTMB)
library(DHARMa)
library(lubridate)
library(broom.mixed) 
library(dotwhisker)

# figure settings
source("code/settings/figure_settings.R")

# import data
fwri <- read_csv("intermediate-data/FWRI_invasive_plant_formatted.csv")
fwc <- read_csv("intermediate-data/FWC_invasive_plant_formatted.csv",
                col_types = list(PrevPropCovered = col_double(),
                                 PrevPropCoveredAdj = col_double(),
                                 PrevAreaCoveredRaw_ha = col_double(),
                                 SurveyDays = col_double(),
                                 RatioCovered = col_double(),
                                 LogRatioCovered = col_double(),
                                 LogitPrevPropCovered = col_double(),
                                 LogRatioCovered = col_double(),
                                 LogitPrevPropCovered = col_double(),
                                 LogRatioCovered = col_double()))
qual <- read_csv("intermediate-data/LW_quality_formatted.csv")


#### combine data ####

# summarize qual (historical data)
qual_sum <- qual %>%
  group_by(PermanentID) %>%
  summarize(Secchi_hist = mean(Secchi, na.rm = T)) %>%
  ungroup()

# combine fwc and fwri
comb <- fwc %>%
  select(PermanentID, AreaName, Area_ha, GSYear, SurveyDate, CommonName, PropCovered, SurveyorExperience, SurveyorExperienceB, Edensa_Present, Nguad_Present) %>%
  filter(!is.na(PropCovered)) %>%
  rename_with(.cols = c(AreaName, Area_ha, SurveyDate, PropCovered, SurveyorExperience, SurveyorExperienceB, Edensa_Present, Nguad_Present), 
              ~ paste0("FWC_", .x)) %>%
  inner_join(fwri %>%
               select(PermanentID, AreaName, Area_ha, GSYear, SurveyDate, CommonName, 
                      PropCovered1, PropCovered2, PropCovered3, Nguad_Present) %>%
               filter(!is.na(PropCovered1)) %>%
               rename_with(.cols = c(AreaName, Area_ha, SurveyDate, Nguad_Present, 
                                     PropCovered1, PropCovered2, PropCovered3), 
                           ~ paste0("FWRI_", .x))) %>%
  left_join(qual %>%
              select(PermanentID, GSYear, Secchi, n_qual)) %>% # current year turbidity
  left_join(qual_sum) %>% # historical turbidity
  mutate(DateDiff = as.numeric(difftime(FWC_SurveyDate,
                                        FWRI_SurveyDate, 
                                        units = "days"))) %>%
  pivot_longer(cols = starts_with("FWRI_PropCovered"),
               names_to = "FWRI_PropType",
               names_prefix = "FWRI_PropCovered",
               values_to = "FWRI_PropCovered") %>%
  mutate(FWRI_PropType = fct_recode(as.character(FWRI_PropType),
                                    "present" = "1",
                                    "moderate-to-dense" = "2",
                                    "dense" = "3"),
         FWC_SurveyMonth = 11 - month(FWC_SurveyDate), # max month = 11
         CommonName = tolower(CommonName),
         CoverDiff = FWC_PropCovered - FWRI_PropCovered,
         AbsCoverDiff = abs(CoverDiff))

# data points
comb %>% 
  mutate(lake_year = paste0(PermanentID, GSYear)) %>% 
  summarise(points = n_distinct(lake_year),
            lakes = n_distinct(PermanentID))

# check that areas match
comb %>%
  select(PermanentID, FWRI_Area_ha, FWC_Area_ha) %>%
  unique() %>%
  filter(FWRI_Area_ha != FWC_Area_ha)
# yes

# quality samples
comb %>%
  select(PermanentID, GSYear, n_qual) %>%
  unique() %>%
  mutate(n_qual = replace_na(n_qual, 0)) %>%
  ggplot(aes(x = n_qual)) +
  geom_histogram(binwidth = 1)
# most do not have water quality data
# fewer may have Secchi data

# missing current Secchi data
comb %>%
  filter(is.na(Secchi)) %>%
  select(PermanentID, FWC_AreaName, FWRI_AreaName) %>%
  unique() # 69 lakes

# missing historical Secchi data
comb %>%
  filter(is.na(Secchi_hist)) %>%
  select(PermanentID, FWC_AreaName, FWRI_AreaName) %>%
  unique() # 7 lakes


#### correlations ####

comb_cor <- comb %>%
  group_by(CommonName, FWRI_PropType) %>%
  summarise(Cor = cor.test(FWC_PropCovered, FWRI_PropCovered)$estimate,
            P_value = cor.test(FWC_PropCovered, FWRI_PropCovered)$p.value) %>%
  ungroup() %>%
  mutate(Cor_text = as.character(round_half_up(Cor, digits = 2)),
         Star = if_else(P_value < 0.05, "*", ""),
         Cor_P = paste0(Cor_text, Star))


#### correlation figure ####

pdf("output/fwri_fwc_correlations.pdf", width = 4, height = 4)
ggplot(comb_cor, aes(x = CommonName, y = Cor, 
                     fill = as.factor(FWRI_PropType))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Cor_P), 
            position = position_dodge(width = 0.9),
            vjust = 0,
            size = paper_text_size) +
  scale_fill_viridis_d(name = "FWRI population type") +
  labs(x = "Plant species", y = "Pearson correlation coefficient") +
  def_theme_paper +
  theme(legend.position = c(0.77, 0.83))
dev.off()


#### cover difference figure ####

pdf("output/fwri_fwc_proportion_difference_histogram.pdf", width = 6.5, height = 6.5)
ggplot(comb, aes(x = CoverDiff)) +
  geom_histogram(binwidth = 0.01) +
  facet_grid(CommonName ~ FWRI_PropType, scales = "free") +
  labs(x = "Annual survey - FWRI survey",
       y = "Lake-year combinations") +
  def_theme_paper
dev.off()


#### cover difference regressions ####

# function to process data
mod_dat_filt <- function(Species, PropType){
  
  dat_out <- comb %>%
    filter(CommonName == Species & 
             FWRI_PropType == PropType) %>%
    mutate(LakeArea_s = (FWRI_Area_ha - min(FWRI_Area_ha)) / sd(FWRI_Area_ha),
           DateDiff_s = (DateDiff - min(DateDiff)) / sd(DateDiff),
           SurveyorExperience_s = (FWC_SurveyorExperience - mean(FWC_SurveyorExperience)) / sd(FWC_SurveyorExperience),
           SurveyorExperienceB = fct_relevel(FWC_SurveyorExperienceB, "medium", "low"))
  
  dat_out_qual <- dat_out %>%
    # filter(!is.na(Secchi)) %>% # limits dataset to use yearly Secchi
    filter(!is.na(Secchi_hist)) %>%
    # mutate(Turbidity_s = (max(Secchi) - Secchi) / sd(Secchi))
    mutate(Turbidity_s = (max(Secchi_hist) - Secchi_hist) / sd(Secchi_hist))
  
  
  return(list(dat_out, dat_out_qual))
  
}

# filter for hydrilla and each pop type
hydr1 <- mod_dat_filt("hydrilla", "present")[[1]]
hydr2 <- mod_dat_filt("hydrilla", "moderate-to-dense")[[1]]
hydr3 <- mod_dat_filt("hydrilla", "dense")[[1]]
wahy1 <- mod_dat_filt("water hyacinth", "present")[[1]]
wahy2 <- mod_dat_filt("water hyacinth", "moderate-to-dense")[[1]]
wahy3 <- mod_dat_filt("water hyacinth", "dense")[[1]]
wale1 <- mod_dat_filt("water lettuce", "present")[[1]]
wale2 <- mod_dat_filt("water lettuce", "moderate-to-dense")[[1]]
wale3 <- mod_dat_filt("water lettuce", "dense")[[1]]

# filter for water quality data
hydr1_qual <- mod_dat_filt("hydrilla", "present")[[2]]
hydr2_qual <- mod_dat_filt("hydrilla", "moderate-to-dense")[[2]]
hydr3_qual <- mod_dat_filt("hydrilla", "dense")[[2]]
wahy1_qual <- mod_dat_filt("water hyacinth", "present")[[2]]
wahy2_qual <- mod_dat_filt("water hyacinth", "moderate-to-dense")[[2]]
wahy3_qual <- mod_dat_filt("water hyacinth", "dense")[[2]]
wale1_qual <- mod_dat_filt("water lettuce", "present")[[2]]
wale2_qual <- mod_dat_filt("water lettuce", "moderate-to-dense")[[2]]
wale3_qual <- mod_dat_filt("water lettuce", "dense")[[2]]

# visualize scaled variables
ggplot(hydr1, aes(x = LakeArea_s)) +
  geom_histogram()

ggplot(hydr1, aes(x = DateDiff_s)) +
  geom_histogram()

ggplot(hydr1, aes(x = SurveyorExperience_s)) +
  geom_histogram()

ggplot(hydr1, aes(x = FWC_SurveyMonth)) +
  geom_histogram()

ggplot(hydr1, aes(x = FWRI_Nguad_Present)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = LakeArea_s)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = DateDiff_s)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = SurveyorExperience_s)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = Turbidity_s)) +
  geom_histogram()

# check correlation between variables
cor.test(~ LakeArea_s + DateDiff_s, data = hydr1) # not sig
cor.test(~ LakeArea_s + SurveyorExperience_s, data = hydr1) # 0.2
cor.test(~ SurveyorExperience_s + DateDiff_s, data = hydr1) # 0.2
cor.test(~ LakeArea_s + FWC_SurveyMonth, data = hydr1) # -0.3
cor.test(~ FWC_SurveyMonth + SurveyorExperience_s, data = hydr1) # -0.2
cor.test(~ FWC_SurveyMonth + DateDiff_s, data = hydr1) # -0.7

cor.test(~ LakeArea_s + DateDiff_s, data = hydr1_qual) # not sig
cor.test(~ LakeArea_s + SurveyorExperience_s, data = hydr1_qual) # 0.2
cor.test(~ SurveyorExperience_s + DateDiff_s, data = hydr1_qual) # 0.2
cor.test(~ LakeArea_s + FWC_SurveyMonth, data = hydr1_qual) # -0.3
cor.test(~ FWC_SurveyMonth + SurveyorExperience_s, data = hydr1_qual) # -0.2
cor.test(~ FWC_SurveyMonth + DateDiff_s, data = hydr1_qual) # -0.7
cor.test(~ LakeArea_s + Turbidity_s, data = hydr1_qual) # 0.3
cor.test(~ Turbidity_s + SurveyorExperience_s, data = hydr1_qual) # not sig
cor.test(~ Turbidity_s + DateDiff_s, data = hydr1_qual) # not sig

# models
# opted for glmmTMB because we can't include time-invariant variables like LakeArea_s and Turbidity_s
# with a fixed effect model (collinear with lake fixed effect)
# glmmTMB also has flexible families if needed

# hyrilla models full dataset
hydr_mod1 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + FWRI_Nguad_Present + (1|GSYear) + (1|PermanentID), data = hydr1)
summary(hydr_mod1)
plot(simulateResiduals(hydr_mod1))

hydr_mod2 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + FWRI_Nguad_Present + (1|GSYear) + (1|PermanentID), data = hydr2)
summary(hydr_mod2)
plot(simulateResiduals(hydr_mod2))

hydr_mod3 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + FWRI_Nguad_Present + (1|GSYear) + (1|PermanentID), data = hydr3)
summary(hydr_mod3)
plot(simulateResiduals(hydr_mod3))

# hydrilla models water Turbidity dataset
hydr_qual_mod1 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + FWRI_Nguad_Present + Turbidity_s + (1|GSYear) + (1|PermanentID), data = hydr1_qual)
summary(hydr_qual_mod1)
plot(simulateResiduals(hydr_qual_mod1))

hydr_qual_mod2 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + FWRI_Nguad_Present + Turbidity_s + (1|GSYear) + (1|PermanentID), data = hydr2_qual)
summary(hydr_qual_mod2)
plot(simulateResiduals(hydr_qual_mod2))

hydr_qual_mod3 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + FWRI_Nguad_Present + Turbidity_s + (1|GSYear) + (1|PermanentID), data = hydr3_qual)
summary(hydr_qual_mod3)
plot(simulateResiduals(hydr_qual_mod3))

# water hyacinth models full dataset
wahy_mod1 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wahy1)
summary(wahy_mod1)
plot(simulateResiduals(wahy_mod1))

wahy_mod2 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wahy2)
summary(wahy_mod2)
plot(simulateResiduals(wahy_mod2))

wahy_mod3 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wahy3)
summary(wahy_mod3)
plot(simulateResiduals(wahy_mod3))

# water hyacinth models water Turbidity dataset
wahy_qual_mod1 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + Turbidity_s + (1|GSYear) + (1|PermanentID), data = wahy1_qual)
summary(wahy_qual_mod1)
plot(simulateResiduals(wahy_qual_mod1))

wahy_qual_mod2 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + Turbidity_s + (1|GSYear) + (1|PermanentID), data = wahy2_qual)
summary(wahy_qual_mod2)
plot(simulateResiduals(wahy_qual_mod2))

wahy_qual_mod3 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + Turbidity_s + (1|GSYear) + (1|PermanentID), data = wahy3_qual)
summary(wahy_qual_mod3)
plot(simulateResiduals(wahy_qual_mod3))

# water lettuce models full dataset
wale_mod1 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wale1)
summary(wale_mod1)
plot(simulateResiduals(wale_mod1))

wale_mod2 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wale2)
summary(wale_mod2)
plot(simulateResiduals(wale_mod2))

wale_mod3 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wale3)
summary(wale_mod3)
plot(simulateResiduals(wale_mod3))

# water lettuce models water Turbidity dataset
wale_qual_mod1 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + Turbidity_s + (1|GSYear) + (1|PermanentID), data = wale1_qual)
summary(wale_qual_mod1)
plot(simulateResiduals(wale_qual_mod1))

wale_qual_mod2 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + Turbidity_s + (1|GSYear) + (1|PermanentID), data = wale2_qual)
summary(wale_qual_mod2)
plot(simulateResiduals(wale_qual_mod2))

wale_qual_mod3 <- glmmTMB(CoverDiff ~ FWC_SurveyMonth + LakeArea_s + SurveyorExperience_s + Turbidity_s + (1|GSYear) + (1|PermanentID), data = wale3_qual)
summary(wale_qual_mod3)
plot(simulateResiduals(wale_qual_mod3))


#### values for text ####

# survyeor experience
mean(hydr2_qual$FWC_SurveyorExperience)
sd(hydr2_qual$FWC_SurveyorExperience)
sd(hydr2_qual$FWC_Area_ha)
sd(wahy2_qual$FWC_Area_ha)
sd(hydr2_qual$Secchi_hist)

# model summaries
summary(hydr_qual_mod2)
summary(wahy_qual_mod2)
summary(wale_qual_mod2)


#### regression coefficient plot ####

# terms
mod_terms <- tibble(term = c("(Intercept)", "FWC_SurveyMonth", "LakeArea_s",
                             "SurveyorExperience_s", "FWRI_Nguad_Present",
                             "Turbidity_s"),
                    variable = c("intercept", "month", "area", "surveyor",
                                 "so. naiad", "turbidity"))
# extract model summary
hydr_qual_tid2 <- tidy(hydr_qual_mod2) %>%
  mutate(model = "hydrilla")

wahy_qual_tid2 <- tidy(wahy_qual_mod2) %>%
  mutate(model = "water hyacinth")

wale_qual_tid2 <- tidy(wale_qual_mod2) %>%
  mutate(model = "water lettuce")

mod_qual_tid <- hydr_qual_tid2 %>%
  full_join(wahy_qual_tid2) %>%
  full_join(wale_qual_tid2) %>%
  left_join(mod_terms)%>%
  filter(effect == "fixed") %>%
  mutate(term = fct_relevel(variable, "intercept", "area",
                            "month", "surveyor",
                            "turbidity", "so. naiad"))

# figure
pdf("output/fwri_fwc_comparison_coefficients.pdf", width = 4, height = 5)
dwplot(mod_qual_tid,
       vline = geom_vline(
         xintercept = 0,
         colour = "grey60",
         linetype = 2)) +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.77, 0.8))
dev.off()


#### correlation with lake area figure ####

# FWRI present
area_fig1 <- ggplot(filter(comb, FWRI_PropType == "present"), 
       aes(FWRI_PropCovered, FWC_PropCovered)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.6, aes(color = log(FWRI_Area_ha))) +
  facet_wrap(~ CommonName, scales = "free") +
  scale_color_distiller(name = "Lake area\n(log-ha)", 
                        direction = 1,
                        palette = "Blues") +
  labs(x = "FWRI survey: proportion occupied (present)") +
  def_theme_paper +
  theme(axis.title.y = element_blank())

# FWRI moderate
area_fig2 <- ggplot(filter(comb, FWRI_PropType == "moderate-to-dense"), 
       aes(FWRI_PropCovered, FWC_PropCovered)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.6, aes(color = log(FWRI_Area_ha))) +
  facet_wrap(~ CommonName, scales = "free") +
  scale_color_distiller(name = "Lake area\n(log-ha)", 
                        direction = 1,
                        palette = "Blues") +
  labs(x = "FWRI survey: proportion occupied (moderate-to-dense)") +
  def_theme_paper +
  theme(legend.position = "none",
        axis.title.y = element_blank())

# FWRI dense
area_fig3 <- ggplot(filter(comb, FWRI_PropType == "dense"), 
       aes(FWRI_PropCovered, FWC_PropCovered)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.6, aes(color = log(FWRI_Area_ha))) +
  facet_wrap(~ CommonName, scales = "free") +
  scale_color_distiller(name = "Lake area\n(log-ha)", 
                        direction = 1,
                        palette = "Blues") +
  labs(x = "FWRI survey: proportion occupied (dense)") +
  def_theme_paper +
  theme(legend.position = "none",
        axis.title.y = element_blank())

# axes titles
area_figy <- textGrob("Annual survey: proportion occupied",
                      gp = gpar(fontsize = 10), rot = 90)

# combine panels
area_pan <- plot_grid(area_fig1 + theme(legend.position = "none"), 
                      area_fig2, area_fig3,
                      nrow = 3)

# add legend
area_leg <- get_legend(area_fig1)
area_fig <- plot_grid(area_pan, area_leg,
                      nrow = 1,
                      rel_widths = c(1, 0.13))

# add axes
pdf("output/fwri_fwc_correlations_lake_area.pdf", 
    width = 6.5, height = 6)
grid.arrange(arrangeGrob(area_fig, 
                         left = area_figy))
dev.off()


#### explore data ####

filter(comb, CommonName == "Water hyacinth" &
         (FWRI_PropCovered1 > 0.15 | FWC_PropCovered > 0.15))

