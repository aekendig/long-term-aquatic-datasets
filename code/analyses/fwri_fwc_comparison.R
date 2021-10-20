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

comb <- fwc %>%
  select(PermanentID, AreaName, Area_ha, GSYear, SurveyDate, CommonName, PropCovered, SurveyorExperience, SurveyorExperienceB) %>%
  filter(!is.na(PropCovered)) %>%
  rename_with(.cols = c(AreaName, Area_ha, SurveyDate, PropCovered, SurveyorExperience, SurveyorExperienceB), 
              ~ paste0("FWC_", .x)) %>%
  inner_join(fwri %>%
               select(PermanentID, AreaName, Area_ha, GSYear, SurveyDate, CommonName, 
                      PropCovered1, PropCovered2, PropCovered3) %>%
               filter(!is.na(PropCovered1)) %>%
               rename_with(.cols = c(AreaName, Area_ha, SurveyDate, 
                                     PropCovered1, PropCovered2, PropCovered3), 
                           ~ paste0("FWRI_", .x))) %>%
  left_join(qual %>%
              select(PermanentID, GSYear, Secchi, n_qual)) %>%
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
    filter(!is.na(Secchi)) %>%
    mutate(Clarity_s = (max(Secchi) - Secchi) / sd(Secchi))
  
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

# visualize scales variables
ggplot(hydr1, aes(x = LakeArea_s)) +
  geom_histogram()

ggplot(hydr1, aes(x = DateDiff_s)) +
  geom_histogram()

ggplot(hydr1, aes(x = SurveyorExperience_s)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = LakeArea_s)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = DateDiff_s)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = SurveyorExperience_s)) +
  geom_histogram()

ggplot(hydr1_qual, aes(x = Clarity_s)) +
  geom_histogram()

# check correlation between variables
cor.test(~ LakeArea_s + DateDiff_s, data = hydr1)
cor.test(~ LakeArea_s + SurveyorExperience_s, data = hydr1)
cor.test(~ SurveyorExperience_s + DateDiff_s, data = hydr1)
# uncorrelated or small correlations (< 0.2)

cor.test(~ LakeArea_s + DateDiff_s, data = hydr1_qual)
cor.test(~ LakeArea_s + SurveyorExperience_s, data = hydr1_qual)
cor.test(~ SurveyorExperience_s + DateDiff_s, data = hydr1_qual)
cor.test(~ LakeArea_s + Clarity_s, data = hydr1_qual) # largest cor: 0.4
cor.test(~ Clarity_s + SurveyorExperience_s, data = hydr1_qual)
cor.test(~ Clarity_s + DateDiff_s, data = hydr1_qual)

# hyrilla models full dataset
hydr_mod1 <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s * SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = hydr1)
summary(hydr_mod1)
hydr_mod1b <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = hydr1)
summary(hydr_mod1b)
plot(simulateResiduals(hydr_mod1b))

hydr_mod2 <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s * SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = hydr2)
summary(hydr_mod2)
hydr_mod2b <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = hydr2)
summary(hydr_mod2b)
plot(simulateResiduals(hydr_mod2b))

hydr_mod3 <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s * SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = hydr3)
summary(hydr_mod3)
hydr_mod3b <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = hydr3)
summary(hydr_mod3b)
plot(simulateResiduals(hydr_mod3b))

# hydrilla models water clarity dataset
hydr_qual_mod1 <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s * ( LakeArea_s + Clarity_s) + (1|GSYear) + (1|PermanentID), data = hydr1_qual)
summary(hydr_qual_mod1)
hydr_qual_mod1b <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s + LakeArea_s + Clarity_s, data = hydr1_qual)
summary(hydr_qual_mod1b)
plot(simulateResiduals(hydr_qual_mod1b))

hydr_qual_mod2 <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s * ( LakeArea_s + Clarity_s) + (1|GSYear) + (1|PermanentID), data = hydr2_qual)
summary(hydr_qual_mod2)
hydr_qual_mod2b <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s + LakeArea_s + Clarity_s, data = hydr2_qual)
summary(hydr_qual_mod2b)
plot(simulateResiduals(hydr_qual_mod2b))

hydr_qual_mod3 <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s * ( LakeArea_s + Clarity_s) + (1|GSYear) + (1|PermanentID), data = hydr3_qual)
summary(hydr_qual_mod3)
hydr_qual_mod3b <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s + LakeArea_s + Clarity_s, data = hydr3_qual)
summary(hydr_qual_mod3b)
plot(simulateResiduals(hydr_qual_mod3b))

# water hyacinth models full dataset
wahy_mod1 <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s * SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wahy1)
summary(wahy_mod1)
wahy_mod1b <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wahy1)
summary(wahy_mod1b)
plot(simulateResiduals(wahy_mod1b))

wahy_mod2 <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s * SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wahy2)
summary(wahy_mod2)
# sig interaction
plot(simulateResiduals(wahy_mod2))

wahy_mod3 <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s * SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wahy3)
summary(wahy_mod3)
# sig interaction
plot(simulateResiduals(wahy_mod3))

# water hyacinth models water clarity dataset
wahy_qual_mod1 <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s * ( LakeArea_s + Clarity_s) + (1|GSYear) + (1|PermanentID), data = wahy1_qual)
# can't converge
wahy_qual_mod1b <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s + LakeArea_s + Clarity_s, data = wahy1_qual)
summary(wahy_qual_mod1b)
plot(simulateResiduals(wahy_qual_mod1b))

wahy_qual_mod2 <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s * ( LakeArea_s + Clarity_s) + (1|GSYear) + (1|PermanentID), data = wahy2_qual)
# can't converge
wahy_qual_mod2b <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s + LakeArea_s + Clarity_s, data = wahy2_qual)
summary(wahy_qual_mod2b)
plot(simulateResiduals(wahy_qual_mod2b))

wahy_qual_mod3 <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s * ( LakeArea_s + Clarity_s) + (1|GSYear) + (1|PermanentID), data = wahy3_qual)
summary(wahy_qual_mod3)
wahy_qual_mod3b <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s + LakeArea_s + Clarity_s, data = wahy3_qual)
summary(wahy_qual_mod3b)
plot(simulateResiduals(wahy_qual_mod3b))

# water lettuce models full dataset
wale_mod1 <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s * SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wale1)
summary(wale_mod1)
wale_mod1b <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wale1)
summary(wale_mod1b)
plot(simulateResiduals(wale_mod1b))

wale_mod2 <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s * SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wale2)
summary(wale_mod2)
wale_mod2b <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wale2)
summary(wale_mod2b)
plot(simulateResiduals(wale_mod2b))

wale_mod3 <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s * SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wale3)
summary(wale_mod3)
wale_mod3b <- glmmTMB(CoverDiff ~ DateDiff_s + LakeArea_s + SurveyorExperience_s + (1|GSYear) + (1|PermanentID) , data = wale3)
summary(wale_mod3b)
plot(simulateResiduals(wale_mod3b))

# water lettuce models water clarity dataset
wale_qual_mod1 <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s * ( LakeArea_s + Clarity_s) + (1|GSYear) + (1|PermanentID), data = wale1_qual)
summary(wale_qual_mod1)
wale_qual_mod1b <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s + LakeArea_s + Clarity_s, data = wale1_qual)
summary(wale_qual_mod1b)
plot(simulateResiduals(wale_qual_mod1b))

wale_qual_mod2 <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s * ( LakeArea_s + Clarity_s) + (1|GSYear) + (1|PermanentID), data = wale2_qual)
summary(wale_qual_mod2)
wale_qual_mod2b <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s + LakeArea_s + Clarity_s, data = wale2_qual)
summary(wale_qual_mod2b)
plot(simulateResiduals(wale_qual_mod2b))

wale_qual_mod3 <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s * ( LakeArea_s + Clarity_s) + (1|GSYear) + (1|PermanentID), data = wale3_qual)
summary(wale_qual_mod3)
wale_qual_mod3b <- glmmTMB(CoverDiff ~ DateDiff_s + SurveyorExperience_s + LakeArea_s + Clarity_s, data = wale3_qual)
summary(wale_qual_mod3b)
plot(simulateResiduals(wale_qual_mod3b))


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

