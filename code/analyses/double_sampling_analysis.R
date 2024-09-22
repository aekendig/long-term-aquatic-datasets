#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom.mixed)
library(glmmTMB)
library(DHARMa)
library(coin)

# figure settings
source("code/settings/figure_settings.R")

# import data
fwri <- read_csv("intermediate-data/fwri_double_sampling_data.csv")
fwc <- read_csv("intermediate-data/fwc_double_sampling_data.csv")


#### format data ####

# combine data
dat <- fwri %>%
  rename(Year = fwri_Year) %>%
  full_join(fwc %>%
              rename(Year = fwc_Year)) %>%
  mutate(DateDiff = fwc_Date - fwri_Date,
         YearF = as.factor(Year))

# sample sizes
n_distinct(dat$TaxonName) # 62 taxa
n_distinct(dat$PermanentID) # 86 lakes
n_distinct(dat$Year) # 5 years
distinct(dat, PermanentID, Year) %>% nrow() # 210 lake-years
211*62 - 62 # each row is each species in each lake-year combo


#### richness ####

# number of taxa observed by each survey
dat_rich <- dat %>%
  group_by(PermanentID, Year, YearF, DateDiff) %>%
  summarize(NTaxa = n_distinct(TaxonName),
            DetectedBoth = sum(fwri_IsDetected == "Yes" & fwc_IsDetected == "Yes"),
            FwriOnly = sum(fwri_IsDetected == "Yes" & fwc_IsDetected == "No"),
            FwcOnly = sum(fwri_IsDetected == "No" & fwc_IsDetected == "Yes"),
            FWRI = sum(fwri_IsDetected == "Yes"),
            FWC = sum(fwc_IsDetected == "Yes"),
            MethodDiff = FWC - FWRI,
            DetectedTotal = sum(fwri_IsDetected == "Yes" | fwc_IsDetected == "Yes"),
            .groups = "drop")

# waterbodies
dat_rich_permID <- sort(unique(dat_rich$PermanentID))

# double check for issues
pdf("output/double_sampling_richness_over_time.pdf")
# cycle through waterbodies
for(i in dat_rich_permID){
  
  # select waterbody
  dat_sub <- filter(dat_rich, PermanentID == i)
  
  # make graph
  print(ggplot(dat_sub, aes(x = Year, y = FWC)) +
          geom_point(color = "blue") +
          geom_line(color = "blue") +
          geom_point(aes(y = FWRI), color = "red") +
          geom_line(aes(y = FWRI), color = "red") +
          ggtitle(i))
}
dev.off()
# these all look okay

# unique years
sort(unique(dat_rich$Year))
# checked these years in double_sampling_richness_raw, and these all look okay

# visualize
ggplot(dat_rich, aes(x = FWRI, y = FWC)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  def_theme
# FWC consistently detected more taxa than FWRI

ggplot(dat_rich) +
  geom_density(aes(x = FwriOnly), color = "blue") +
  geom_density(aes(x = FwcOnly), color = "red") +
  def_theme
# FWC often detects several taxa that FWRI doesn't
# FWRI often detects none or a couple that FWC doesn't

ggplot(dat_rich, aes(x = DetectedTotal, y = DetectedBoth)) +
  geom_point() +
  def_theme
# for lower richness waterbody-lake surveys, the number of 
# taxa detected by both surveys is ~constant, but when 
# richness is high, the number of taxa detected by both
# surveys increases with higher richness

ggplot(dat_rich, aes(x = DateDiff, y = MethodDiff)) +
  geom_point() +
  def_theme
# date doesn't seem related

ggplot(dat_rich, aes(x = MethodDiff)) +
  geom_histogram() +
  def_theme

# analysis
rich_mod <- glmmTMB(MethodDiff ~ DateDiff + (1|PermanentID) + (1|YearF),
                    data = dat_rich)
summary(rich_mod)
rich_mod_res <- simulateResiduals(rich_mod, n = 1000)
plot(rich_mod_res)
confint(rich_mod)


#### richness over time ####

# change in richness
dat_rich_change <- dat_rich %>%
  arrange(PermanentID, Year) %>%
  group_by(PermanentID) %>%
  mutate(YearDiff = Year - lag(Year),
         FwriChange = (FWRI - lag(FWRI)) / YearDiff,
         FwcChange = (FWC - lag(FWC)) / YearDiff) %>%
  ungroup() %>%
  filter(!is.na(FwriChange)) %>%
  mutate(MethodDiff = FwcChange - FwriChange)

# sample sizes
n_distinct(dat_rich_change$PermanentID) # 45
nrow(dat_rich_change) # 124

# ranges of values
unique(dat_rich_change$YearDiff)
range(dat_rich_change$FwriChange)
range(dat_rich_change$FwcChange)

# averages
mean(dat_rich_change$FwriChange)
mean(dat_rich_change$FwcChange)

# visualize
ggplot(dat_rich_change, aes(x = FwriChange, y = FwcChange)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  def_theme

ggplot(dat_rich_change) +
  geom_density(aes(x = FwriChange), color = "blue") +
  geom_density(aes(x = FwcChange), color = "red") +
  def_theme
# more FWC changes ~ 0
# FWRI changes may be detecting species that were not detected
# in prior/later survey

ggplot(dat_rich_change, aes(x = MethodDiff)) +
  geom_histogram() +
  def_theme

# analysis
rich_change_mod <- glmmTMB(MethodDiff ~ 1 + (1|PermanentID) + (1|YearF),
                           data = dat_rich_change)
summary(rich_change_mod)
rich_change_mod_res <- simulateResiduals(rich_change_mod, n = 1000)
plot(rich_change_mod_res)
confint(rich_change_mod, parm = "(Intercept)")


#### taxon-specific detections ####

# make data long
dat_long <- dat %>%
  select(TaxonName, PermanentID, Year, YearF, fwri_IsDetected, fwc_IsDetected) %>%
  pivot_longer(cols = c(fwri_IsDetected, fwc_IsDetected),
               values_to = "Detection",
               names_to = "Method") %>%
  mutate(Method = str_remove(Method, "_IsDetected"),
         Detection = if_else(Detection == "Yes", 1, 0))

# check taxa for all 1's or 0's
# got convergence warning when all taxa models were fit
taxon_sum <- dat_long %>%
  group_by(TaxonName) %>%
  summarize(Detections = sum(Detection) / n(),
            .groups = "drop")

arrange(taxon_sum, Detections)
arrange(taxon_sum, desc(Detections))
# check the two extremes

# list of taxa
taxa <- sort(unique(dat$TaxonName))

# set i
i <- 1

# subset data
dat_sub <- filter(dat_long, TaxonName == taxa[i])

# analysis
mod_sub <- glmmTMB(Detection ~ Method + (1|PermanentID) + (1|YearF),
                   data = dat_sub, family = "binomial")

# use this to manually cycle through taxa (some models had warnings below)
# i <- i + 1

# taxa with model convergence issues
# 13: Crinum americanum - only 2 presences in FWRI
# 32: Nelumbo lutea - not sure why this model isn't converging

# look at models
summary(mod_sub)
mod_sub_res <- simulateResiduals(mod_sub, n = 1000)
plot(mod_sub_res, title = taxa[i])

# save results when i <- 1
taxa_coefs <- tidy(mod_sub) %>%
  filter(term %in% c("(Intercept)", "Methodfwri")) %>%
  mutate(TaxonName = taxa[1]) %>%
  relocate(TaxonName)

# loop through taxa
pdf("output/double_sampling_taxon_detection.pdf")

for(i in taxa[-c(13, 32)]) {
  
  # subset data
  dat_sub <- filter(dat_long, TaxonName == i)
  
  # fit model
  mod_sub <- glmmTMB(Detection ~ Method + (1|PermanentID) + (1|YearF),
                     data = dat_sub, family = "binomial")
  
  # extract and plot residuals
  mod_sub_res <- simulateResiduals(mod_sub, n = 1000)
  print(plot(mod_sub_res, title = i))
  
  # save model coefficients
  taxa_coefs <- taxa_coefs %>%
    full_join(tidy(mod_sub) %>%
                filter(term %in% c("(Intercept)", "Methodfwri")) %>%
                mutate(TaxonName = i))
}

dev.off()

# remove redundancy from first taxon
# correct method p-values
taxa_coefs_method <- taxa_coefs[-c(1:2),] %>%
  filter(term == "Methodfwri") %>%
  mutate(q.value = p.adjust(p.value, method = "fdr"),
         lower = estimate - 1.96 * std.error,
         upper = estimate + 1.96 * std.error,
         odds_perc = 100 * (exp(estimate) - 1),
         lower_odds_perc = 100 * (exp(lower) - 1),
         upper_odds_perc = 100 * (exp(upper) - 1))

# taxa with significant differences by q-value
taxa_coefs_method %>%
  filter(q.value < 0.05) %>%
  mutate(method_higher = if_else(estimate > 0, "FWRI", "FWC")) %>%
  group_by(method_higher) %>%
  summarize(n = n(),
            min_odds_perc = min(odds_perc),
            max_odds_perc = max(odds_perc))

# check with p-values
taxa_coefs_method %>%
  filter(p.value < 0.05) %>%
  mutate(method_higher = if_else(estimate > 0, "FWRI", "FWC")) %>%
  group_by(method_higher) %>%
  summarize(n = n(),
            min_odds_perc = min(odds_perc),
            max_odds_perc = max(odds_perc))
# same

# look at three taxa
taxa_coefs_method %>%
  filter(q.value < 0.05 & estimate > 0)

taxa_coefs_method %>%
  filter(p.value < 0.05 & estimate > 0)


#### invasive species abundance ####

# filter for focal invasive species
dat_inv <- dat %>%
  filter(TaxonName %in% c("Hydrilla verticillata", "Pistia stratiotes", "Eichhornia crassipes")) %>%
  mutate(fwc_PAC = 100 * fwc_PAC,
         fwri_PAC1 = 100 * fwri_PAC1,
         fwri_PAC2 = 100 * fwri_PAC2,
         fwri_PAC3 = 100 * fwri_PAC3,
         MethodDiff1 = fwc_PAC - fwri_PAC1,
         MethodDiff2 = fwc_PAC - fwri_PAC2,
         MethodDiff3 = fwc_PAC - fwri_PAC3)

# visualize
ggplot(dat_inv, aes(x = fwri_PAC1, y = fwc_PAC)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~ TaxonName, scales = "free") +
  def_theme

ggplot(dat_inv, aes(x = fwri_PAC2, y = fwc_PAC)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~ TaxonName, scales = "free") +
  def_theme

ggplot(dat_inv, aes(x = fwri_PAC3, y = fwc_PAC)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~ TaxonName, scales = "free") +
  def_theme
# correlation looks like it depends on the FWRI abundance
# measurement, which has different optimization for different
# taxa

ggplot(dat_inv) +
  geom_density(aes(x = MethodDiff1), color = "blue") +
  geom_density(aes(x = MethodDiff2), color = "red") +
  geom_density(aes(x = MethodDiff3), color = "purple") +
  facet_wrap(~ TaxonName, scales = "free") +
  def_theme
# most differences are very small

# split data by taxa
dat_hydr <- filter(dat_inv, TaxonName == "Hydrilla verticillata")
dat_pist <- filter(dat_inv, TaxonName == "Pistia stratiotes")
dat_eich <- filter(dat_inv, TaxonName == "Eichhornia crassipes")

# fit models
hydr_mod1 <- glmmTMB(MethodDiff1 ~ DateDiff + (1|PermanentID) + (1|YearF),
                     data = dat_hydr)
summary(hydr_mod1) # model wasn't fit
hydr_mod1b <- update(hydr_mod1, .~.-DateDiff)
summary(hydr_mod1b)
hydr_mod_res1 <- simulateResiduals(hydr_mod1b, n = 1000)
plot(hydr_mod_res1)
ggplot(tibble(x = hydr_mod_res1$scaledResiduals), aes(x = x)) +
  geom_density() # residuals are highly concentrated
# may be overfit, remove random effect with very low estimate
hydr_mod1c <- update(hydr_mod1b, .~.-(1|YearF))
summary(hydr_mod1c)
hydr_mod_res1 <- simulateResiduals(hydr_mod1c, n = 1000)
plot(hydr_mod_res1) # no change
hydr_mod1d <- glmmTMB(MethodDiff1 ~ DateDiff + (1|PermanentID) + (1|YearF),
                      data = dat_hydr, family = "t_family")
summary(hydr_mod1d)
hydr_mod_res1 <- simulateResiduals(hydr_mod1d, n = 1000)
plot(hydr_mod_res1) 
mean(dat_hydr$MethodDiff1) # very different from model estimate

# paired tests and correlations
shapiro.test(dat_hydr$MethodDiff1) # not normal, know this from above
median(dat_hydr$MethodDiff1)
range(dat_hydr$MethodDiff1)
wilcoxsign_test(dat_hydr$fwc_PAC ~ dat_hydr$fwri_PAC1,
                distribution = "exact")
cor.test(dat_hydr$fwc_PAC, dat_hydr$fwri_PAC1)

median(dat_hydr$MethodDiff2)
range(dat_hydr$MethodDiff2)
wilcoxsign_test(dat_hydr$fwc_PAC ~ dat_hydr$fwri_PAC2,
                distribution = "exact")
cor.test(dat_hydr$fwc_PAC, dat_hydr$fwri_PAC2)

median(dat_hydr$MethodDiff3)
range(dat_hydr$MethodDiff3)
wilcoxsign_test(dat_hydr$fwc_PAC ~ dat_hydr$fwri_PAC3,
                distribution = "exact")
cor.test(dat_hydr$fwc_PAC, dat_hydr$fwri_PAC3)

shapiro.test(dat_eich$MethodDiff1)
median(dat_eich$MethodDiff1)
range(dat_eich$MethodDiff1)
wilcoxsign_test(dat_eich$fwc_PAC ~ dat_eich$fwri_PAC1,
                distribution = "exact")
cor.test(dat_eich$fwc_PAC, dat_eich$fwri_PAC1)

median(dat_eich$MethodDiff2)
range(dat_eich$MethodDiff2)
wilcoxsign_test(dat_eich$fwc_PAC ~ dat_eich$fwri_PAC2,
                distribution = "exact")
cor.test(dat_eich$fwc_PAC, dat_eich$fwri_PAC2)

median(dat_eich$MethodDiff3)
range(dat_eich$MethodDiff3)
wilcoxsign_test(dat_eich$fwc_PAC ~ dat_eich$fwri_PAC3,
                distribution = "exact")
cor.test(dat_eich$fwc_PAC, dat_eich$fwri_PAC3)

shapiro.test(dat_pist$MethodDiff1)
median(dat_pist$MethodDiff1)
range(dat_pist$MethodDiff1)
wilcoxsign_test(dat_pist$fwc_PAC ~ dat_pist$fwri_PAC1,
                distribution = "exact")
cor.test(dat_pist$fwc_PAC, dat_pist$fwri_PAC1)

median(dat_pist$MethodDiff2)
range(dat_pist$MethodDiff2)
wilcoxsign_test(dat_pist$fwc_PAC ~ dat_pist$fwri_PAC2,
                distribution = "exact")
cor.test(dat_pist$fwc_PAC, dat_pist$fwri_PAC2)

median(dat_pist$MethodDiff3)
range(dat_pist$MethodDiff3)
wilcoxsign_test(dat_pist$fwc_PAC ~ dat_pist$fwri_PAC3,
                distribution = "exact")
cor.test(dat_pist$fwc_PAC, dat_pist$fwri_PAC3)
