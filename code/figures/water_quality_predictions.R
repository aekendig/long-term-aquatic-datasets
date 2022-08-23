#### setup ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(pals)

# figure settings
source("code/settings/figure_settings.R")

# functions
source("code/generic-functions/cut_mean.R")

# import significant water quality results
foc_chl_sig <- read_csv("output/fwc_focal_invasive_chlorophyll_significant.csv")
foc_sec_sig <- read_csv("output/fwc_focal_invasive_secchi_significant.csv")
foc_nit_sig <- read_csv("output/fwc_focal_invasive_nitrogen_significant.csv")
foc_pho_sig <- read_csv("output/fwc_focal_invasive_phosphorus_significant.csv")

# import model summaries
foc_chl <- read_csv("output/fwc_focal_chlorophyll_model_summary.csv")
foc_sec <- read_csv("output/fwc_focal_secchi_model_summary.csv")
foc_nit <- read_csv("output/fwc_focal_nitrogen_model_summary.csv")
foc_pho <- read_csv("output/fwc_focal_phosphorus_model_summary.csv")

# import raw data
chl_dat <- read_csv("intermediate-data/FWC_chlorophyll_analysis_formatted.csv")
sec_dat <- read_csv("intermediate-data/FWC_secchi_analysis_formatted.csv", guess_max = 7000)
nit_dat <- read_csv("intermediate-data/FWC_nitrogen_analysis_formatted.csv")
pho_dat <- read_csv("intermediate-data/FWC_phosphorus_analysis_formatted.csv")

# color palette
col_pal <- kelly()[c(3:6, 10:11)]


#### edit data ####

# quarter names
quart_name <- tibble(Quarter = 1:4,
                     QuarterF = c("Apr-Jun", "Jul-Sep", "Oct-Dec", "Jan-Mar") %>%
                       fct_relevel("Apr-Jun", "Jul-Sep", "Oct-Dec"))

# quarter color palette
quart_col_pal <- col_pal[1:4]
names(quart_col_pal) <- levels(quart_name$QuarterF)

# invasive species names
foc_inv_names <- tibble(Invasive = c("hydrilla", "water hyacinth", "water lettuce"),
                        CommonName = c("Hydrilla", "Water hyacinth", "Water lettuce")) %>%
  mutate(PanelNamePAC = paste(c("(A)", "(B)", "(C)"), Invasive),
         PanelNameTreat = paste(PanelNamePAC, "management"))

# sig tables assume P<0.1
# change the threshold if needed
p_thresh <- 0.1

# raw data
chl_dat2 <- chl_dat %>%
  left_join(quart_name) %>%
  left_join(foc_inv_names) %>%
  mutate(Treated = Lag3Treated * 3)

sec_dat2 <- sec_dat %>%
  left_join(quart_name) %>%
  left_join(foc_inv_names) %>%
  mutate(Treated = Lag3Treated * 3)

nit_dat2 <- nit_dat %>%
  left_join(quart_name) %>%
  left_join(foc_inv_names) %>%
  mutate(Treated = Lag3Treated * 3)

pho_dat2 <- pho_dat %>%
  left_join(quart_name) %>%
  left_join(foc_inv_names) %>%
  mutate(Treated = Lag3Treated * 3)

# select values based on p_thresh
# add raw data
foc_chl_sig2 <- foc_chl %>%
  filter(P < p_thresh) %>%
  select(Invasive, Quarter, Term, Estimate) %>%
  mutate(Term = fct_recode(Term, "PAC" = "invasive PAC",
                           "Treat" = "management")) %>%
  pivot_wider(names_from = Term, # do it this way because foc_chl_sig has them on separate rows
              values_from = Estimate,
              names_glue = "Beta{Term}") %>%
  left_join(foc_chl_sig %>%
              select(Invasive, Quarter, Intercept, Metric) %>%
              unique()) %>%
  rename(QuarterF = Quarter) %>%
  left_join(chl_dat2 %>%
              select(QuarterF, Invasive, PanelNamePAC, PanelNameTreat, Lag3AvgPercCovered) %>%
              unique())
  
foc_sec_sig2 <- foc_sec %>%
  filter(P < p_thresh) %>%
  select(Invasive, Quarter, Term, Estimate) %>%
  mutate(Term = fct_recode(Term, "PAC" = "invasive PAC",
                           "Treat" = "management")) %>%
  pivot_wider(names_from = Term,
              values_from = Estimate,
              names_glue = "Beta{Term}") %>%
  left_join(foc_sec_sig %>%
              select(Invasive, Quarter, Intercept, Metric) %>%
              unique()) %>%
  rename(QuarterF = Quarter) %>%
  left_join(sec_dat2 %>%
              select(QuarterF, Invasive, PanelNamePAC, PanelNameTreat, Lag3AvgPercCovered) %>%
              unique())

foc_nit_sig2 <- foc_nit %>%
  filter(P < p_thresh) %>%
  select(Invasive, Quarter, Term, Estimate) %>%
  mutate(Term = fct_recode(Term, "PAC" = "invasive PAC",
                           "Treat" = "management")) %>%
  pivot_wider(names_from = Term,
              values_from = Estimate,
              names_glue = "Beta{Term}") %>%
  left_join(foc_nit_sig %>%
              select(Invasive, Quarter, Intercept, Metric) %>%
              unique()) %>%
  rename(QuarterF = Quarter) %>%
  left_join(nit_dat2 %>%
              select(QuarterF, Invasive, PanelNamePAC, PanelNameTreat, Lag3AvgPercCovered) %>%
              unique())

foc_pho_sig2 <- foc_pho %>%
  filter(P < p_thresh) %>%
  select(Invasive, Quarter, Term, Estimate) %>%
  mutate(Term = fct_recode(Term, "PAC" = "invasive PAC")) %>%
  pivot_wider(names_from = Term,
              values_from = Estimate,
              names_glue = "Beta{Term}") %>%
  left_join(foc_pho_sig %>%
              select(Invasive, Quarter, Intercept, Metric) %>%
              unique()) %>%
  rename(QuarterF = Quarter) %>%
  left_join(pho_dat2 %>%
              select(QuarterF, Invasive, PanelNamePAC, PanelNameTreat, Lag3AvgPercCovered) %>%
              unique())

# combine datasets
# estimate fitted values
# standardize fitted values
qual_sig <- foc_chl_sig2 %>%
  full_join(foc_sec_sig2) %>%
  full_join(foc_nit_sig2) %>%
  full_join(foc_pho_sig2) %>%
  mutate(QuarterF = fct_relevel(QuarterF, "Apr-Jun", "Jul-Sep", "Oct-Dec"))

qual_PAC_fit <- qual_sig %>%
  filter(!is.na(BetaPAC)) %>%
  mutate(FittedPAC = Intercept + BetaPAC * Lag3AvgPercCovered,  # PAC-only effect
         expFittedPAC = exp(FittedPAC)) %>%
  group_by(Invasive, Metric, QuarterF) %>%
  mutate(FittedPACStd = (FittedPAC - mean(FittedPAC)) / sd(FittedPAC),
         AvgPAC = mean(Lag3AvgPercCovered),
         AvgPAC1SD = sd(Lag3AvgPercCovered)) %>%
  ungroup() %>%
  mutate(Lag3PACMin = AvgPAC - AvgPAC1SD,
         Lag3PACMax = AvgPAC + AvgPAC1SD)

qual_treat_fit <- qual_sig %>%
  filter(!is.na(BetaTreat)) %>%
  select(-Lag3AvgPercCovered) %>%
  unique() %>%
  expand_grid(Lag3Treated = c(0, 1/3, 2/3, 1)) %>%
  mutate(Treated = Lag3Treated * 3,
         FittedTreat = Intercept + BetaTreat * Lag3Treated, # treatment-only effect
         expFittedTreat = exp(FittedTreat)) %>%
  group_by(Invasive, Metric, QuarterF) %>%
  mutate(FittedTreatStd = (FittedTreat - mean(FittedTreat)) / sd(FittedTreat)) %>%
  ungroup()

# summarized chlorophyll data
foc_chl_dat <- chl_dat2 %>%
  inner_join(foc_chl_sig2 %>%
               select(Invasive, QuarterF, BetaPAC, BetaTreat) %>%
               unique()) %>%
  mutate(QuarterF = fct_relevel(QuarterF, "Apr-Jun", "Jul-Sep", "Oct-Dec"))

foc_chl_dat_PAC <- foc_chl_dat %>%
  filter(!is.na(BetaPAC)) %>%
  group_by(Invasive, PanelNamePAC, QuarterF) %>%
  mutate(BinPAC = cut_interval(Lag3AvgPercCovered, n = 3)) %>%
  group_by(Invasive, PanelNamePAC, QuarterF, BinPAC) %>%
  summarize(BinPACMid = cut_mean(BinPAC),
            BinPACMean = mean(Lag3AvgPercCovered),
            BinN = length(Lag3AvgPercCovered),
            BinPACSE = sd(Lag3AvgPercCovered) / sqrt(BinN),
            QualityMean = mean(QualityValue),
            QualitySE = sd(QualityValue) / sqrt(length(QualityValue))) %>%
  ungroup()

foc_chl_dat_Treat <- foc_chl_dat %>%
  filter(!is.na(BetaTreat))

# summarized Secchi data
foc_sec_dat <- sec_dat2 %>%
  inner_join(foc_sec_sig2 %>%
               select(Invasive, QuarterF, BetaPAC, BetaTreat) %>%
               unique()) %>%
  mutate(QuarterF = fct_relevel(QuarterF, "Apr-Jun", "Jul-Sep", "Oct-Dec"))

foc_sec_dat_PAC <- foc_sec_dat %>%
  filter(!is.na(BetaPAC)) %>%
  group_by(Invasive, PanelNamePAC, QuarterF) %>%
  mutate(BinPAC = cut_interval(Lag3AvgPercCovered, n = 3)) %>%
  group_by(Invasive, PanelNamePAC, QuarterF, BinPAC) %>%
  summarize(BinPACMid = cut_mean(BinPAC),
            BinPACMean = mean(Lag3AvgPercCovered),
            BinN = length(Lag3AvgPercCovered),
            BinPACSE = sd(Lag3AvgPercCovered) / sqrt(BinN),
            QualityMean = mean(QualityValue),
            QualitySE = sd(QualityValue) / sqrt(length(QualityValue))) %>%
  ungroup()

foc_sec_dat_Treat <- foc_sec_dat %>%
  filter(!is.na(BetaTreat))

# summarized nitrogen data
foc_nit_dat <- nit_dat2 %>%
  inner_join(foc_nit_sig2 %>%
               select(Invasive, QuarterF, BetaPAC, BetaTreat) %>%
               unique()) %>%
  mutate(QuarterF = fct_relevel(QuarterF, "Apr-Jun", "Jul-Sep", "Oct-Dec"))

foc_nit_dat_PAC <- foc_nit_dat %>%
  filter(!is.na(BetaPAC)) %>%
  group_by(Invasive, PanelNamePAC, QuarterF) %>%
  mutate(BinPAC = cut_interval(Lag3AvgPercCovered, n = 3)) %>%
  group_by(Invasive, PanelNamePAC, QuarterF, BinPAC) %>%
  summarize(BinPACMid = cut_mean(BinPAC),
            BinPACMean = mean(Lag3AvgPercCovered),
            BinN = length(Lag3AvgPercCovered),
            BinPACSE = sd(Lag3AvgPercCovered) / sqrt(BinN),
            QualityMean = mean(QualityValue),
            QualitySE = sd(QualityValue) / sqrt(length(QualityValue))) %>%
  ungroup()

foc_nit_dat_Treat <- foc_nit_dat %>%
  filter(!is.na(BetaTreat))

# summarized phosphorus data
foc_pho_dat <- pho_dat2 %>%
  inner_join(foc_pho_sig2 %>%
               select(Invasive, QuarterF, BetaPAC) %>%
               unique()) %>%
  mutate(QuarterF = fct_relevel(QuarterF, "Apr-Jun", "Jul-Sep", "Oct-Dec"))

foc_pho_dat_PAC <- foc_pho_dat %>%
  filter(!is.na(BetaPAC)) %>%
  group_by(Invasive, PanelNamePAC, QuarterF) %>%
  mutate(BinPAC = cut_interval(Lag3AvgPercCovered, n = 3)) %>%
  group_by(Invasive, PanelNamePAC, QuarterF, BinPAC) %>%
  summarize(BinPACMid = cut_mean(BinPAC),
            BinPACMean = mean(Lag3AvgPercCovered),
            BinN = length(Lag3AvgPercCovered),
            BinPACSE = sd(Lag3AvgPercCovered) / sqrt(BinN),
            QualityMean = mean(QualityValue),
            QualitySE = sd(QualityValue) / sqrt(length(QualityValue))) %>%
  ungroup()


#### all-metric figures ####

qual_PAC_fit %>%
  #filter(Lag3AvgPercCovered >= Lag3PACMin & Lag3AvgPercCovered <= Lag3PACMax) %>%
  ggplot(aes(x = Lag3AvgPercCovered, y = FittedPACStd)) +
  geom_line(aes(color = Metric, linetype = QuarterF), alpha = 0.7) +
  facet_wrap(~ PanelNamePAC, scales = "free") +
  labs(x = "3-year average PAC", y = "Standardized water quality") +
  scale_color_manual(values = col_pal) +
  scale_linetype(name = "Quarter") +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))

qual_treat_fit %>%
  ggplot(aes(x = Treated, y = FittedTreatStd)) +
  geom_line(aes(color = Metric, linetype = QuarterF), alpha = 0.7) +
  facet_wrap(~ PanelNameTreat, scales = "free") +
  labs(x = "Years managed (out of 3)", y = "Standardized water quality") +
  scale_color_manual(values = col_pal) +
  scale_linetype(name = "Quarter") +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))


#### PAC figures by metric ####

foc_chl_PAC_fig <- qual_PAC_fit %>%
  filter(Metric == "chlorophyll a") %>%
  ggplot(aes(color = QuarterF)) +
  geom_line(aes(x = Lag3AvgPercCovered, y = expFittedPAC), alpha = 0.5) +
  geom_errorbar(data = foc_chl_dat_PAC,
                aes(x = BinPACMean, y = QualityMean,
                    ymin = QualityMean - QualitySE, ymax = QualityMean + QualitySE),
                width = 0) +
  geom_errorbarh(data = foc_chl_dat_PAC,
                 aes(y = QualityMean,
                     xmin = BinPACMean - BinPACSE, xmax = BinPACMean + BinPACSE),
                 height = 0) +
  geom_point(data = foc_chl_dat_PAC,
             aes(x = BinPACMean, y = QualityMean),
             size = 1) +
  facet_wrap(~ PanelNamePAC, scales = "free",
             labeller =  labeller(PanelNamePAC = c("(B) water hyacinth" = "(A) water hyacinth",
                                                   "(C) water lettuce" = "(B) water lettuce"))) +
  labs(x = "3-year average PAC",
       y = expression(paste("Chlorophyll ", italic(a), " (", mu, "g/L)", sep = ""))) +
  scale_color_manual(values = quart_col_pal, name = "Quarter", limits = force) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))

foc_sec_PAC_fig <- qual_PAC_fit %>%
  filter(Metric == "Secchi disk depth") %>%
  ggplot(aes(color = QuarterF)) +
  geom_line(aes(x = Lag3AvgPercCovered, y = expFittedPAC), alpha = 0.5) +
  geom_errorbar(data = foc_sec_dat_PAC,
                aes(x = BinPACMean, y = QualityMean,
                    ymin = QualityMean - QualitySE, ymax = QualityMean + QualitySE),
                width = 0) +
  geom_errorbarh(data = foc_sec_dat_PAC,
                 aes(y = QualityMean,
                     xmin = BinPACMean - BinPACSE, xmax = BinPACMean + BinPACSE),
                 height = 0) +
  geom_point(data = foc_sec_dat_PAC,
             aes(x = BinPACMean, y = QualityMean),
             size = 1) +
  facet_wrap(~ PanelNamePAC, scales = "free",
             labeller =  labeller(PanelNamePAC = c("(A) hydrilla" = "(A) hydrilla",
                                                   "(C) water lettuce" = "(B) water lettuce"))) +
  labs(x = "3-year average PAC",
       y = "Secchi disk depth (ft)") +
  scale_color_manual(values = quart_col_pal, name = "Quarter", limits = force) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))
  

foc_nit_PAC_fig <- qual_PAC_fit %>%
  filter(Metric == "total nitrogen") %>%
  ggplot(aes(color = QuarterF)) +
  geom_line(aes(x = Lag3AvgPercCovered, y = expFittedPAC), alpha = 0.5) +
  geom_errorbar(data = foc_nit_dat_PAC,
                aes(x = BinPACMean, y = QualityMean,
                    ymin = QualityMean - QualitySE, ymax = QualityMean + QualitySE),
                width = 0) +
  geom_errorbarh(data = foc_nit_dat_PAC,
                 aes(y = QualityMean,
                     xmin = BinPACMean - BinPACSE, xmax = BinPACMean + BinPACSE),
                 height = 0) +
  geom_point(data = foc_nit_dat_PAC,
             aes(x = BinPACMean, y = QualityMean),
             size = 1) +
  facet_wrap(~ PanelNamePAC, scales = "free",
             labeller =  labeller(PanelNamePAC = c("(B) water hyacinth" = "(A) water hyacinth",
                                                   "(C) water lettuce" = "(B) water lettuce"))) +
  labs(x = "3-year average PAC",
       y = expression(paste("Total nitrogen (", mu, "g/L)", sep = ""))) +
  scale_color_manual(values = quart_col_pal, name = "Quarter", limits = force) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))

foc_pho_PAC_fig <- qual_PAC_fit %>%
  filter(Metric == "total phosphorus") %>%
  ggplot(aes(color = QuarterF)) +
  geom_line(aes(x = Lag3AvgPercCovered, y = expFittedPAC), alpha = 0.5) +
  geom_errorbar(data = foc_pho_dat_PAC,
                aes(x = BinPACMean, y = QualityMean,
                    ymin = QualityMean - QualitySE, ymax = QualityMean + QualitySE),
                width = 0) +
  geom_errorbarh(data = foc_pho_dat_PAC,
                 aes(y = QualityMean,
                     xmin = BinPACMean - BinPACSE, xmax = BinPACMean + BinPACSE),
                 height = 0) +
  geom_point(data = foc_pho_dat_PAC,
             aes(x = BinPACMean, y = QualityMean),
             size = 1) +
  facet_wrap(~ PanelNamePAC, scales = "free") +
  labs(x = "3-year average PAC",
       y = expression(paste("Total phosphorus (", mu, "g/L)", sep = ""))) +
  scale_color_manual(values = quart_col_pal, name = "Quarter", limits = force) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))


#### treatment figures by metric ####

foc_chl_treat_fig <- qual_treat_fit %>%
  filter(Metric == "chlorophyll a") %>%
  ggplot(aes(color = QuarterF)) +
  geom_line(aes(x = Treated, y = expFittedTreat), alpha = 0.5) +
  stat_summary(data = foc_chl_dat_Treat, aes(x = Treated, y = QualityValue),
               geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(data = foc_chl_dat_Treat, aes(x = Treated, y = QualityValue),
               geom = "point", fun = "mean", size = 1) +
  facet_wrap(~ PanelNameTreat, scales = "free",
             labeller =  labeller(PanelNameTreat = c("(A) hydrilla management" = "(A) hydrilla management",
                                                   "(C) water lettuce management" = "(B) water lettuce management"))) +
  labs(x = "Years managed (out of 3)",
       y = expression(paste("Chlorophyll ", italic(a), " (", mu, "g/L)", sep = ""))) +
  scale_color_manual(values = quart_col_pal, name = "Quarter", limits = force) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))

foc_sec_treat_fig <- qual_treat_fit %>%
  filter(Metric == "Secchi disk depth") %>%
  ggplot(aes(color = QuarterF)) +
  geom_line(aes(x = Treated, y = expFittedTreat), alpha = 0.5) +
  stat_summary(data = foc_sec_dat_Treat, aes(x = Treated, y = QualityValue),
               geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(data = foc_sec_dat_Treat, aes(x = Treated, y = QualityValue),
               geom = "point", fun = "mean", size = 1) +
  facet_wrap(~ PanelNameTreat, scales = "free",
             labeller =  labeller(PanelNameTreat = c("(A) hydrilla management" = "hydrilla management"))) +
  labs(x = "Years managed (out of 3)",
       y = "Secchi disk depth (ft)") +
  scale_color_manual(values = quart_col_pal, name = "Quarter", limits = force) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))

foc_nit_treat_fig <- qual_treat_fit %>%
  filter(Metric == "total nitrogen") %>%
  ggplot(aes(color = QuarterF)) +
  geom_line(aes(x = Treated, y = expFittedTreat), alpha = 0.5) +
  stat_summary(data = foc_nit_dat_Treat, aes(x = Treated, y = QualityValue),
               geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(data = foc_nit_dat_Treat, aes(x = Treated, y = QualityValue),
               geom = "point", fun = "mean", size = 1) +
  facet_wrap(~ PanelNameTreat, scales = "free",
             labeller =  labeller(PanelNameTreat = c("(A) hydrilla management" = "hydrilla management"))) +
  labs(x = "Years managed (out of 3)",
       y = expression(paste("Total nitrogen (", mu, "g/L)", sep = ""))) +
  scale_color_manual(values = quart_col_pal, name = "Quarter", limits = force) +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))


#### save figures ####

ggsave("output/fwc_focal_invasive_chlorophyll_PAC_prediction.png", foc_chl_PAC_fig,
       device = "png", width = 5, height = 2.5, units = "in")
ggsave("output/fwc_focal_invasive_secchi_PAC_prediction.png", foc_sec_PAC_fig,
       device = "png", width = 5, height = 2.5, units = "in")
ggsave("output/fwc_focal_invasive_nitrogen_PAC_prediction.png", foc_nit_PAC_fig,
       device = "png", width = 5, height = 2.5, units = "in")
ggsave("output/fwc_focal_invasive_phosphorus_PAC_prediction.png", foc_pho_PAC_fig,
       device = "png", width = 7, height = 2.5, units = "in")

ggsave("output/fwc_focal_invasive_chlorophyll_treatment_prediction.png", foc_chl_treat_fig,
       device = "png", width = 5, height = 2.5, units = "in")
ggsave("output/fwc_focal_invasive_secchi_treatment_prediction.png", foc_sec_treat_fig,
       device = "png", width = 3, height = 2.5, units = "in")
ggsave("output/fwc_focal_invasive_nitrogen_treatment_prediction.png", foc_nit_treat_fig,
       device = "png", width = 3, height = 2.5, units = "in")
