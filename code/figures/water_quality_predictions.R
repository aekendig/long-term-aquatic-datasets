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

# invasive species names
foc_inv_names <- tibble(Invasive = c("hydrilla", "water hyacinth", "water lettuce"),
                        CommonName = c("Hydrilla", "Water hyacinth", "Water lettuce")) %>%
  mutate(PanelNamePAC = paste(c("(A)", "(B)", "(C)"), Invasive),
         PanelNameTreat = paste(PanelNamePAC, "management"))

# sig tables assume P<0.1
# change the threshold if needed
p_thresh <- 0.1

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
  left_join(chl_dat %>%
              select(Quarter, CommonName, Lag3AvgPercCovered) %>%
              unique() %>%
              left_join(quart_name) %>%
              left_join(foc_inv_names))
  
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
  left_join(sec_dat %>%
              select(Quarter, CommonName, Lag3AvgPercCovered) %>%
              left_join(quart_name) %>%
              left_join(foc_inv_names))

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
  left_join(nit_dat %>%
              select(Quarter, CommonName, Lag3AvgPercCovered) %>%
              left_join(quart_name) %>%
              left_join(foc_inv_names))

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
  left_join(pho_dat %>%
              select(Quarter, CommonName, Lag3AvgPercCovered) %>%
              left_join(quart_name) %>%
              left_join(foc_inv_names))

# combine datasets
# estimate fitted values
# standardize fitted values
qual_sig <- foc_chl_sig2 %>%
  full_join(foc_sec_sig2) %>%
  full_join(foc_nit_sig2) %>%
  full_join(foc_pho_sig2) 

qual_PAC_fit <- qual_sig %>%
  filter(!is.na(BetaPAC)) %>%
  mutate(FittedPAC = Intercept + BetaPAC * Lag3AvgPercCovered,  # treatment-only effect
         expFittedPAC = exp(FittedPAC)) %>%
  group_by(CommonName, Metric, QuarterF) %>%
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
  group_by(CommonName, Metric, QuarterF) %>%
  mutate(FittedTreatStd = (FittedTreat - mean(FittedTreat)) / sd(FittedTreat)) %>%
  ungroup()

# raw data
foc_chl_dat <- chl_dat %>%
  inner_join(foc_chl_sig2 %>%
               select(CommonName, Quarter, BetaPAC, BetaTreat) %>%
               unique()) %>%
  group_by(CommonName, Quarter) %>%
  mutate(BinPAC = cut_number(log(Lag3AvgPercCovered + 1), n = 3)) %>%
  group_by(CommonName, Quarter, BinPAC) %>%
  mutate(BinPACMean = cut_mean(BinPAC)) %>%
  ungroup() %>%
  left_join(quart_name) %>%
  left_join(foc_inv_names) %>%
  mutate(Treated = Lag3Treated * 3)

foc_sec_dat <- sec_dat %>%
  inner_join(foc_sec_sig2 %>%
               select(CommonName, Quarter, BetaPAC, BetaTreat) %>%
               unique()) %>%
  group_by(CommonName, Quarter) %>%
  mutate(BinPAC = cut_number(log(Lag3AvgPercCovered + 1), n = 3)) %>%
  group_by(CommonName, Quarter, BinPAC) %>%
  mutate(BinPACMean = cut_mean(BinPAC)) %>%
  ungroup() %>%
  left_join(quart_name) %>%
  left_join(foc_inv_names) %>%
  mutate(Treated = Lag3Treated * 3)

foc_nit_dat <- nit_dat %>%
  inner_join(foc_nit_sig2 %>%
               select(CommonName, Quarter, BetaPAC, BetaTreat) %>%
               unique()) %>%
  group_by(CommonName, Quarter) %>%
  mutate(BinPAC = cut_number(log(Lag3AvgPercCovered + 1), n = 3)) %>%
  group_by(CommonName, Quarter, BinPAC) %>%
  mutate(BinPACMean = cut_mean(BinPAC)) %>%
  ungroup() %>%
  left_join(quart_name) %>%
  left_join(foc_inv_names) %>%
  mutate(Treated = Lag3Treated * 3)

foc_pho_dat <- pho_dat %>%
  inner_join(foc_pho_sig2 %>%
               select(CommonName, Quarter, BetaPAC) %>%
               unique()) %>%
  group_by(CommonName, Quarter) %>%
  mutate(BinPAC = cut_number(log(Lag3AvgPercCovered + 1), n = 3)) %>%
  group_by(CommonName, Quarter, BinPAC) %>%
  mutate(BinPACMean = cut_mean(BinPAC)) %>%
  ungroup() %>%
  left_join(quart_name) %>%
  left_join(foc_inv_names) %>%
  mutate(Treated = Lag3Treated * 3)


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


#### figures by metric ####

#### start here ####
# apply below to all metrics and make treatment figures

# PAC figures
qual_PAC_fit %>%
  filter(Metric == "chlorophyll a") %>%
  ggplot(aes(color = QuarterF)) +
  geom_line(aes(x = log(Lag3AvgPercCovered + 1), y = expFittedPAC), alpha = 0.5) +
  geom_point(data = filter(foc_chl_dat, !is.na(BetaPAC)),
             aes(x = log(Lag3AvgPercCovered + 1), y = QualityValue),
             alpha = 0.5) +
  # stat_summary(data = filter(foc_chl_dat, !is.na(BetaPAC)),
  #              aes(x = log(BinPACMean + 1), y = QualityValue),
  #              geom = "errorbar", fun.data = "mean_se", width = 0) +
  # stat_summary(data = filter(foc_chl_dat, !is.na(BetaPAC)),
  #              aes(x = log(BinPACMean + 1), y = QualityValue),
  #              geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ PanelNamePAC, scales = "free",
             labeller =  labeller(PanelNamePAC = c("(B) water hyacinth" = "(A) water hyacinth",
                                                   "(C) water lettuce" = "(B) water lettuce"))) +
  labs(x = "3-year average PAC)",
       y = expression(paste("Chlorophyll ", italic(a), " (", mu, "g/L)", sep = ""))) +
  scale_color_manual(values = col_pal, name = "Quarter") +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))

qual_PAC_fit %>%
  filter(Metric == "Secchi disk depth") %>%
  ggplot(aes(color = QuarterF)) +
  geom_line(aes(x = Lag3AvgPercCovered, y = expFittedPAC), alpha = 0.5) +
  stat_summary(data = filter(foc_sec_dat, !is.na(BetaPAC)),
               aes(x = BinPACMean, y = QualityValue),
               geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(data = filter(foc_sec_dat, !is.na(BetaPAC)),
               aes(x = BinPACMean, y = QualityValue),
               geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ PanelNamePAC, scales = "free") +
  labs(x = "3-year average PAC)",
       y = "Secchi disk depth (ft)") +
  scale_color_manual(values = col_pal, name = "Quarter") +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0))
  

qual_PAC_fit %>%
  filter(Metric == "total nitrogen")

qual_PAC_fit %>%
  filter(Metric == "total phosphorus")