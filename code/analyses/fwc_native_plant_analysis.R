#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(inspectdf) # inspect_cor
library(broom) # glance, tidy
library(pals) # color palettes
library(patchwork)
library(gllvm)
library(Hmsc)
library(corrplot)

# figure settings
source("code/settings/figure_settings.R")

# functions
source("code/generic-functions/proportion_transformations.R")

# import data
dat <- read_csv("intermediate-data/FWC_common_native_plants_invaded_data_formatted.csv")
dat_full <- read_csv("intermediate-data/FWC_common_native_plants_invasive_species_data_formatted.csv")


#### edit data ####

# # summarize
# dat2 <- dat %>%
#   group_by(CommonName, PermanentID, TaxonName, Habitat) %>%
#   summarize(AvgPAC = mean(PercCovered * 100),
#             TreatFreq = mean(Treated),
#             YearsDetected = sum(Detected),
#             YearsSurveyed = n_distinct(GSYear),
#             LogYearsSurveyed = log(YearsSurveyed)) %>%
#   ungroup() %>%
#   group_by(CommonName, TaxonName) %>%
#   mutate(AvgPAC_c = AvgPAC - mean(AvgPAC),
#          TreatFreq_c = TreatFreq - mean(TreatFreq),
#          AvgPAC_s = AvgPAC_c / sd(AvgPAC),
#          TreatFreq_s = TreatFreq_c / sd(TreatFreq)) %>%
#   ungroup() %>%
#   mutate(Habitat = tolower(Habitat))
#
# # split data by invasive species
# hydr_dat <- filter(dat2, CommonName == "Hydrilla") %>%
#   arrange(PermanentID, TaxonName)
# flpl_dat <- filter(dat2, CommonName == "floating plants") %>%
#   arrange(PermanentID, TaxonName)
# cubu_dat <- filter(dat2, CommonName == "Cuban bulrush") %>%
#   arrange(PermanentID, TaxonName)
# pagr_dat <- filter(dat2, CommonName == "Para grass") %>%
#   arrange(PermanentID, TaxonName)
# torp_dat <- filter(dat2, CommonName == "Torpedograss") %>%
#   arrange(PermanentID, TaxonName)
# 
# # site-by-species matrices
# # remove taxa with no presences
# hydr_mat <- hydr_dat %>%
#   select(PermanentID, TaxonName, YearsDetected) %>%
#   pivot_wider(names_from = TaxonName,
#               values_from = YearsDetected) %>%
#   select(-PermanentID) %>%
#   select(where(~ sum(.) != 0)) %>% 
#   as.matrix()
# flpl_mat <- flpl_dat %>%
#   select(PermanentID, TaxonName, YearsDetected) %>%
#   pivot_wider(names_from = TaxonName,
#               values_from = YearsDetected) %>%
#   select(-PermanentID) %>%
#   select(where(~ sum(.) != 0)) %>% 
#   as.matrix()
# cubu_mat <- cubu_dat %>%
#   select(PermanentID, TaxonName, YearsDetected) %>%
#   pivot_wider(names_from = TaxonName,
#               values_from = YearsDetected) %>%
#   select(-PermanentID) %>%
#   select(where(~ sum(.) != 0)) %>% 
#   as.matrix()
# pagr_mat <- pagr_dat %>%
#   select(PermanentID, TaxonName, YearsDetected) %>%
#   pivot_wider(names_from = TaxonName,
#               values_from = YearsDetected) %>%
#   select(-PermanentID) %>%
#   select(where(~ sum(.) != 0)) %>% 
#   as.matrix()
# torp_mat <- torp_dat %>%
#   select(PermanentID, TaxonName, YearsDetected) %>%
#   pivot_wider(names_from = TaxonName,
#               values_from = YearsDetected) %>%
#   select(-PermanentID) %>%
#   select(where(~ sum(.) != 0)) %>% 
#   as.matrix()
# 
# # covariates
# hydr_cov <- hydr_dat %>%
#   distinct(PermanentID, AvgPAC_s, TreatFreq_s, LogYearsSurveyed) %>%
#   arrange(PermanentID) %>%
#   select(-PermanentID) %>%
#   data.frame()
# flpl_cov <- flpl_dat %>%
#   distinct(PermanentID, AvgPAC_s, TreatFreq_s, LogYearsSurveyed) %>%
#   arrange(PermanentID) %>%
#   select(-PermanentID) %>%
#   data.frame()
# cubu_cov <- cubu_dat %>%
#   distinct(PermanentID, AvgPAC_s, TreatFreq_s, LogYearsSurveyed) %>%
#   arrange(PermanentID) %>%
#   select(-PermanentID) %>%
#   data.frame()
# pagr_cov <- pagr_dat %>%
#   distinct(PermanentID, AvgPAC_s, TreatFreq_s, LogYearsSurveyed) %>%
#   arrange(PermanentID) %>%
#   select(-PermanentID) %>%
#   data.frame()
# torp_cov <- torp_dat %>%
#   distinct(PermanentID, AvgPAC_s, TreatFreq_s, LogYearsSurveyed) %>%
#   arrange(PermanentID) %>%
#   select(-PermanentID) %>%
#   data.frame()
# 
# # correlations
# cor.test(~ AvgPAC_s + TreatFreq_s, data = hydr_cov) # sig, 0.3
# cor.test(~ AvgPAC_s + TreatFreq_s, data = flpl_cov) # not sig
# cor.test(~ AvgPAC_s + TreatFreq_s, data = cubu_cov) # sig, 0.2
# cor.test(~ AvgPAC_s + TreatFreq_s, data = pagr_cov) # not sig
# cor.test(~ AvgPAC_s + TreatFreq_s, data = torp_cov) # not sig

# split data by invasive species
hydr_dat <- filter(dat, CommonName == "Hydrilla") %>%
  arrange(PermanentID, GSYear, TaxonName)
flpl_dat <- filter(dat, CommonName == "floating plants") %>%
  arrange(PermanentID, GSYear, TaxonName)
cubu_dat <- filter(dat, CommonName == "Cuban bulrush") %>%
  arrange(PermanentID, GSYear, TaxonName)
pagr_dat <- filter(dat, CommonName == "Para grass") %>%
  arrange(PermanentID, GSYear, TaxonName)
torp_dat <- filter(dat, CommonName == "Torpedograss") %>%
  arrange(PermanentID, GSYear, TaxonName)

# site-by-species matrices
# remove taxa with no presences
hydr_mat <- hydr_dat %>%
  select(PermanentID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(PermanentID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
flpl_mat <- flpl_dat %>%
  select(PermanentID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(PermanentID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
cubu_mat <- cubu_dat %>%
  select(PermanentID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(PermanentID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
pagr_mat <- pagr_dat %>%
  select(PermanentID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(PermanentID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()
torp_mat <- torp_dat %>%
  select(PermanentID, GSYear, TaxonName, Detected) %>%
  pivot_wider(names_from = TaxonName,
              values_from = Detected) %>%
  select(-c(PermanentID, GSYear)) %>%
  select(where(~ sum(.) != 0)) %>% 
  as.matrix()

# covariates
hydr_cov <- hydr_dat %>%
  distinct(PermanentID, GSYear, PercCovered, RecentTreatment) %>%
  arrange(PermanentID, GSYear) %>%
  select(-c(PermanentID, GSYear)) %>%
  data.frame()
flpl_cov <- flpl_dat %>%
  distinct(PermanentID, GSYear, PercCovered, RecentTreatment) %>%
  arrange(PermanentID, GSYear) %>%
  select(-c(PermanentID, GSYear)) %>%
  data.frame()
cubu_cov <- cubu_dat %>%
  distinct(PermanentID, GSYear, PercCovered, RecentTreatment) %>%
  arrange(PermanentID, GSYear) %>%
  select(-c(PermanentID, GSYear)) %>%
  data.frame()
pagr_cov <- pagr_dat %>%
  distinct(PermanentID, GSYear, PercCovered, RecentTreatment) %>%
  arrange(PermanentID, GSYear) %>%
  select(-c(PermanentID, GSYear)) %>%
  data.frame()
torp_cov <- torp_dat %>%
  distinct(PermanentID, GSYear, PercCovered, RecentTreatment) %>%
  arrange(PermanentID, GSYear) %>%
  select(-c(PermanentID, GSYear)) %>%
  data.frame()

# study design
hydr_stud <- hydr_dat %>%
  distinct(PermanentID, GSYear) %>%
  arrange(PermanentID, GSYear) %>%
  mutate(PermanentID = as.factor(PermanentID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
flpl_stud <- flpl_dat %>%
  distinct(PermanentID, GSYear) %>%
  arrange(PermanentID, GSYear) %>%
  mutate(PermanentID = as.factor(PermanentID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
cubu_stud <- cubu_dat %>%
  distinct(PermanentID, GSYear) %>%
  arrange(PermanentID, GSYear) %>%
  mutate(PermanentID = as.factor(PermanentID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
pagr_stud <- pagr_dat %>%
  distinct(PermanentID, GSYear) %>%
  arrange(PermanentID, GSYear) %>%
  mutate(PermanentID = as.factor(PermanentID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()
torp_stud <- torp_dat %>%
  distinct(PermanentID, GSYear) %>%
  arrange(PermanentID, GSYear) %>%
  mutate(PermanentID = as.factor(PermanentID),
         GSYear = as.factor(GSYear)) %>%
  data.frame()

# traits
hydr_hab <- hydr_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(hydr_hab) <- hydr_hab$TaxonName
hydr_trait <- select(hydr_hab, -TaxonName)
flpl_hab <- flpl_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(flpl_hab) <- flpl_hab$TaxonName
flpl_trait <- select(flpl_hab, -TaxonName)
cubu_hab <- cubu_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(cubu_hab) <- cubu_hab$TaxonName
cubu_trait <- select(cubu_hab, -TaxonName)
pagr_hab <- pagr_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(pagr_hab) <- pagr_hab$TaxonName
pagr_trait <- select(pagr_hab, -TaxonName)
torp_hab <- torp_dat %>%
  distinct(TaxonName, Habitat) %>%
  arrange(TaxonName) %>%
  data.frame()
rownames(torp_hab) <- torp_hab$TaxonName
torp_trait <- select(torp_hab, -TaxonName)


#### fit models with gllvm ####

# # fit Poisson with 2 LVs
# hydr_pois <- gllvm(y = hydr_mat, 
#                    # X = hydr_cov,
#                    # formula = ~ AvgPAC_c * TreatFreq,
#                    family = poisson(),
#                    num.lv = 2,
#                    row_eff = "random",
#                    offset = hydr_cov$LogYearsSurveyed)
# plot(hydr_pois, var.colors = 1)
# # got error with formula included: In nlminb(objr$par, objr$fn, objr$gr, control = list(rel.tol = reltol,  :
# # NA/NaN function evaluation
# # unusual shape in residuals plot, Normal Q-Q has high deviations at ends
# 
# hydr_nb <- gllvm(y = hydr_mat, 
#                    family = "negative.binomial",
#                    num.lv = 2,
#                    row_eff = "random",
#                    offset = hydr_cov$LogYearsSurveyed)
# plot(hydr_nb, var.colors = 1)
# # Normal Q-Q has high deviations at ends still
# # try normal distribution with log-transformed values
# 
# hydr_zip <- gllvm(y = hydr_mat, 
#                    family = "ZIP",
#                    num.lv = 2,
#                    row_eff = "random",
#                    offset = hydr_cov$LogYearsSurveyed)
# plot(hydr_zip, var.colors = 1)
# # same normality issues as poisson
# 
# # compare AIC
# AIC(hydr_pois)
# AIC(hydr_nb) # lowest
# AIC(hydr_zip)
# 
# # try more LVs
# hydr_nb2 <- gllvm(y = hydr_mat, 
#                  family = "negative.binomial",
#                  num.lv = 3,
#                  row_eff = "random",
#                  offset = hydr_cov$LogYearsSurveyed)
# par(mfrow = c(3, 2))
# plot(hydr_nb, var.colors = 1)
# par(mfrow = c(3, 2))
# plot(hydr_nb2, var.colors = 1)
# # better normal QQ
# 
# AIC(hydr_nb)
# AIC(hydr_nb2) # lower
# 
# hydr_nb3 <- gllvm(y = hydr_mat, 
#                   family = "negative.binomial",
#                   num.lv = 4,
#                   row_eff = "random",
#                   offset = hydr_cov$LogYearsSurveyed)
# par(mfrow = c(3, 2))
# plot(hydr_nb2, var.colors = 1)
# par(mfrow = c(3, 2))
# plot(hydr_nb3, var.colors = 1)
# # similar plots
# 
# AIC(hydr_nb2) # lower
# AIC(hydr_nb3)
# 
# # add predictors
# # tried with Avg_c/TreatFreq, Avg_s/Treat_s
# # all had unusually small SE for at least one var
# # and a ton of significant effects (doesn't seem realistic)
# hydr_nb4 <- gllvm(y = hydr_mat, 
#                   X = hydr_cov,
#                   formula = ~ AvgPAC_c * TreatFreq_c,
#                   family = "negative.binomial",
#                   num.lv = 3,
#                   row_eff = "random",
#                   offset = hydr_cov$LogYearsSurveyed)
# par(mfrow = c(3, 2))
# plot(hydr_nb4)
# summary(hydr_nb4)
# par(mfrow = c(1, 1))
# coefplot(hydr_nb4, which.Xcoef = "AvgPAC_s")
# coefplot(hydr_nb4, which.Xcoef = "TreatFreq_s")
# coefplot(hydr_nb4, which.Xcoef = "AvgPAC_s:TreatFreq_s")


#### fit site-summarized models with Hmsc ####
# 
# # hydrilla data, lognormal poisson
# hydr_lpois = Hmsc(Y = hydr_mat,
#                   XData = hydr_cov, 
#                   XFormula = ~LogYearsSurveyed + AvgPAC_s * TreatFreq_s, # no offset option, include years surveyed
#                   XScale = F, # already scaled
#                   distr = "lognormal poisson")
# 
# # MCMC settings
# nChains = 2
# thin = 5
# samples = 100
# transient = 50*thin
# verbose = 50*thin
# 
# # fit model
# hydr_lpois_fit = sampleMcmc(hydr_lpois,
#                             thin = thin,
#                             samples = samples,
#                             transient = transient,
#                             nChains = nChains,
#                             nParallel = nChains,
#                             verbose = verbose)
# 
# # posterior samples
# hydr_lpois_post = convertToCodaObject(hydr_lpois_fit)
# 
# # effective sample sizes
# # total samples = samples * chains
# # and Gelman diagnostics (potential scale reduction factors)
# # one indicates better convergence among chains
# hist(effectiveSize(hydr_lpois_post$Beta), main="ess(beta)")
# hist(gelman.diag(hydr_lpois_post$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")
# 
# # explanatory power, compares posterior distribution to observed values
# hydr_lpois_preds = computePredictedValues(hydr_lpois_fit)
# hydr_lpois_eval = evaluateModelFit(hM = hydr_lpois_fit, predY = hydr_lpois_preds)
# hist(hydr_lpois_eval$RMSE,
#      main = paste0("Mean = ", round(mean(hydr_lpois_eval$RMSE),2)),
#      xlab = "root-mean-square error")
# hist(hydr_lpois_eval$SR2,
#      main = paste0("Mean = ", round(mean(hydr_lpois_eval$SR2),2)),
#      xlab = "pseudo R2")
# hist(hydr_lpois_eval$O.AUC,
#      main = paste0("Mean = ", round(mean(hydr_lpois_eval$O.AUC, na.rm = T),2)),
#      xlab = "prese/abs area under curve")
# 
# # design matrix
# head(hydr_lpois_fit$X)
# 
# # variance partitioning
# hydr_lpois_vp = computeVariancePartitioning(hydr_lpois_fit)
# plotVariancePartitioning(hydr_lpois_fit, VP = hydr_lpois_vp) # manually construct plot, legend blocks image
# 
# # estimates - back-transform these to original scale?
# hydr_lpois_post_beta = getPostEstimate(hydr_lpois_fit, parName="Beta")
# hist(hydr_lpois_post_beta$mean[3,],
#      main = paste0("Mean = ", round(mean(hydr_lpois_post_beta$mean[3,]),2)),
#      xlab = "PAC")
# hist(hydr_lpois_post_beta$mean[4,],
#      main = paste0("Mean = ", round(mean(hydr_lpois_post_beta$mean[4,]),2)),
#      xlab = "Treatment")
# hist(hydr_lpois_post_beta$mean[5,],
#      main = paste0("Mean = ", round(mean(hydr_lpois_post_beta$mean[5,]),2)),
#      xlab = "Interaction")
# plotBeta(hydr_lpois_fit, post = hydr_lpois_post_beta, param = "Support", supportLevel = 0.95)
# plotBeta(hydr_lpois_fit, post = hydr_lpois_post_beta, param = "Mean", supportLevel = 0.95)
# 
# # residual associations among species - need random effects?
# hydr_lpois_omega_cor = computeAssociations(hydr_lpois_fit)
# support_level = 0.95
# hydr_lpois_omega_cor2 = ((hydr_lpois_omega_cor[[1]]$support > support_level) +
#                 (hydr_lpois_omega_cor[[1]]$support < (1-support_level)) > 0) *
#   hydr_lpois_omega_cor[[1]]$mean
# corrplot(hydr_lpois_omega_cor2, method = "color",
#          col = colorRampPalette(c("blue","white","red"))(200),
#          tl.cex = 0.6, tl.col = "black",
#          title = paste("random effect level:", m_fit$rLNames[1]), # not sure what this does
#          mar = c(0,0,1,0))
# # correlated responses to missing covariates or
# # species interactions
# 
# # species richness change over time
# hydr_lpois_pac_gradient = constructGradient(hydr_lpois_fit, focalVariable = "AvgPAC_s")
# hydr_lpois_pac_pred_grad = predict(hydr_lpois_fit, 
#                                    XData = hydr_lpois_pac_gradient$XDataNew,
#                     studyDesign = hydr_lpois_pac_gradient$studyDesignNew,
#                     ranLevels = hydr_lpois_pac_gradient$rLNew)
# plotGradient(hydr_lpois_fit, hydr_lpois_pac_gradient, pred = hydr_lpois_pac_pred_grad, measure = "S",
#              showData = TRUE, jigger = 0.2)
# 
# # estimates and 95% CI
# summary(hydr_lpois_post$Beta)
#
# something seems off about not having omega values
# could include growth form as a trait
# is there a way to include years surveyed as an offset?
# prediction gradient seems off, should all points be on one line
# very high error at end

#### fit site-by-year models with Hmsc ####

# model structure
hydr_mod = Hmsc(Y = hydr_mat,
                XData = hydr_cov, 
                XFormula = ~PercCovered * RecentTreatment,
                TrData = hydr_trait,
                TrFormula = ~Habitat,
                distr = "probit",
                studyDesign = hydr_stud,
                ranLevels = list(GSYear = HmscRandomLevel(units = hydr_stud$GSYear), 
                                 PermanentID = HmscRandomLevel(units = hydr_stud$PermanentID)))

# MCMC settings
nChains = 2
thin = 10
samples = 100
transient = 10*thin
verbose = 10*thin

# fit model
hydr_fit = sampleMcmc(hydr_mod,
                      thin = thin,
                      samples = samples,
                      transient = transient,
                      nChains = nChains,
                      nParallel = nChains,
                      verbose = verbose)

# posterior samples
hydr_post = convertToCodaObject(hydr_fit)

# effective sample sizes
# total samples = samples * chains
# and Gelman diagnostics (potential scale reduction factors)
# one indicates better convergence among chains
hist(effectiveSize(hydr_post$Beta), main="ess(beta)")
hist(gelman.diag(hydr_post$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")
hist(effectiveSize(hydr_post$Omega[[1]]), main="ess(omega)")
hist(gelman.diag(hydr_post$Omega[[1]], multivariate=FALSE)$psrf, main="psrf(omega)") # slow

# explanatory power, compares posterior distribution to observed values
hydr_preds = computePredictedValues(hydr_fit)
hydr_eval = evaluateModelFit(hM = hydr_fit, predY = hydr_preds)
hist(hydr_eval$RMSE,
     main = paste0("Mean = ", round(mean(hydr_eval$RMSE),2)),
     xlab = "root-mean-square error")
hist(hydr_eval$TjurR2,
     main = paste0("Mean = ", round(mean(hydr_eval$TjurR2),2)),
     xlab = "R2")
hist(hydr_eval$AUC,
     main = paste0("Mean = ", round(mean(hydr_eval$AUC),2)),
     xlab = "prese/abs area under curve")

# design matrix
head(hydr_fit$X)

# variance partitioning
hydr_vp = computeVariancePartitioning(hydr_fit)
plotVariancePartitioning(hydr_fit, VP = hydr_vp)

# estimates - back-transform these to original scale?
post_beta = getPostEstimate(hydr_fit, parName="Beta")
hist(post_beta$mean[2,],
     main = paste0("Mean = ", round(mean(post_beta$mean[2,]),2)),
     xlab = "PAC")
hist(post_beta$mean[3,],
     main = paste0("Mean = ", round(mean(post_beta$mean[3,]),2)),
     xlab = "Treatment")
hist(post_beta$mean[4,],
     main = paste0("Mean = ", round(mean(post_beta$mean[4,]),2)),
     xlab = "Interaction")
plotBeta(hydr_fit, post = post_beta, param = "Support", supportLevel = 0.95)
plotBeta(hydr_fit, post = post_beta, param = "Mean", supportLevel = 0.95) # magnitudes are very different, show estimates separately

# residual associations among species
hydr_omega_cor = computeAssociations(hydr_fit)
support_level = 0.95
hydr_omega_cor2 = ((hydr_omega_cor[[1]]$support > support_level) +
                     (hydr_omega_cor[[1]]$support < (1-support_level)) > 0) *
  hydr_omega_cor[[1]]$mean
corrplot(hydr_omega_cor2, method = "color",
         col = colorRampPalette(c("blue","white","red"))(200),
         tl.cex = 0.6, tl.col = "black",
         title = paste("random effect level:", hydr_fit$rLNames[1]), # not sure what this does
         mar = c(0,0,1,0))
# correlated responses to missing covariates or
# species interactions

# species richness change over PAC (can do treatment too)
hydr_pac_gradient = constructGradient(hydr_fit, focalVariable = "PercCovered")
hydr_pac_pred = predict(hydr_fit, XData = hydr_pac_gradient$XDataNew,
                        studyDesign = hydr_pac_gradient$studyDesignNew,
                        ranLevels = hydr_pac_gradient$rLNew)
plotGradient(hydr_fit, hydr_pac_gradient, pred = hydr_pac_pred, measure = "S",
             showData = TRUE, jigger = 0.2)

# trait association
hydr_postGamma = getPostEstimate(hydr_fit, parName = "Gamma")
plotGamma(hydr_fit, post = hydr_postGamma, 
          param = "Support", 
          supportLevel = 0.95)

# estimates and 95% CI
summary(hydr_post$Beta)


#### initial visualizations ####

# covariate correlations
comb_dat %>%
  select(CommonName, PermanentID, TreatFreq, AvgPAC_c) %>%
  unique() %>%
  select(-PermanentID) %>%
  group_by(CommonName) %>%
  inspect_cor() %>% 
  ungroup()
# treatment not significantly correlated with PAC 
# except for hydrilla, but corr = 0.2

ggplot(comb_dat, aes(AvgPAC, YearsDetected, color = TaxonName)) +
  geom_point(size = 0.75, alpha = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")

ggplot(comb_dat, aes(TreatFreq, YearsDetected, color = TaxonName)) +
  geom_point(size = 0.75, alpha = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  facet_wrap(~ CommonName, scales = "free") +
  theme(legend.position = "none")


#### fit models ####

# species with too many 1's or 0's
comb_dat %>%
  group_by(TaxonName, CommonName) %>%
  summarize(YearsDetected = sum(YearsDetected),
            YearsSurveyed = sum(YearsSurveyed)) %>%
  ungroup() %>%
  filter(YearsDetected >= (YearsSurveyed - 20) | YearsDetected <= 20)
# para grass models might not work, some have a few or a lot of detections

# apply model to each invasive plant and taxon
# plant_mods <- comb_dat %>%
#   select(CommonName, TaxonName, Habitat, YearsDetected, YearsUndetected,
#          AvgPAC_c, TreatFreq) %>%
#   nest(data = c(YearsDetected, YearsUndetected, AvgPAC_c, TreatFreq)) %>%
#   mutate(fit = map(data, ~glm(cbind(YearsDetected, YearsUndetected) ~ AvgPAC_c + TreatFreq,
#                               data = ., family = binomial)))
# model fit error

# convert warnings to errors to break loop
# options(warn = 2)

# common/taxa list
# update as pairs break loop
comm_tax <- comb_dat %>%
  select(CommonName, TaxonName) %>%
  unique() %>%
  mutate(SppPair = paste(CommonName, TaxonName, sep = ", ")) %>%
  filter(!(SppPair %in% c("Para grass, Cicuta maculata",
                          "Para grass, Crinum americanum",
                          "Para grass, Salix spp.",
                          "Water hyacinth, Mayaca fluviatilis",
                          "Water hyacinth, Myriophyllum laxum/pinnatum",
                          "Water lettuce, Fuirena spp.",
                          "Water lettuce, Mayaca fluviatilis",
                          "Water lettuce, Micranthemum umbrosum",
                          "Water lettuce, Myriophyllum laxum/pinnatum",
                          "Water lettuce, Nitella spp.",
                          "Water lettuce, Nymphoides aquatica")))

# find issue model
# commented out because it takes a while to run
# for(i in 1:nrow(comm_tax)) {
#     
#     print(comm_tax$SppPair[i])
#     
#     print(glm(cbind(YearsDetected, YearsUndetected) ~ AvgPAC_c + TreatFreq,
#               data = filter(comb_dat, 
#                             CommonName == comm_tax$CommonName[i] & 
#                               TaxonName == comm_tax$TaxonName[i]), 
#               family = binomial))
#   
# }

# remove problematic species pairs
comb_dat2 <- comb_dat %>%
  inner_join(comm_tax)

# refit models
plant_mods2 <- comb_dat2 %>%
  select(CommonName, TaxonName, Habitat, YearsDetected, YearsUndetected, 
         AvgPAC_c, TreatFreq) %>%
  nest(data = c(YearsDetected, YearsUndetected, AvgPAC_c, TreatFreq)) %>%
  mutate(fit = map(data, ~glm(cbind(YearsDetected, YearsUndetected) ~ AvgPAC_c + TreatFreq, 
                              data = ., family = binomial)))


#### coefficients ####

# all coefficients
plant_coef <- plant_mods2 %>%
  mutate(tidied = map(fit, tidy)) %>%
  select(CommonName, TaxonName, Habitat, tidied) %>%
  unnest(tidied)

# models with significant PAC effects
plant_coef_PAC <- plant_coef %>%
  mutate(Sig = if_else(term == "AvgPAC_c" & p.value < 0.1, "yes", "no")) %>%
  filter(term == "AvgPAC_c") %>%
  rename(Coef = estimate) %>%
  select(CommonName, TaxonName, Habitat, Coef, Sig) %>%
  full_join(plant_coef %>%
              filter(term == "(Intercept)") %>%
              select(CommonName, TaxonName, estimate) %>%
              rename(Intercept = estimate))

# models with significant treatment effects
plant_coef_treat <- plant_coef %>%
  mutate(Sig = if_else(term == "TreatFreq" & p.value < 0.1, "yes", "no")) %>%
  filter(term == "TreatFreq") %>%
  rename(Coef = estimate) %>%
  select(CommonName, TaxonName, Habitat, Coef, Sig) %>%
  full_join(plant_coef %>%
              filter(term == "(Intercept)") %>%
              select(CommonName, TaxonName, estimate) %>%
              rename(Intercept = estimate))


#### summary table ####

# count taxa per category
plant_coef_summ <- plant_coef %>%
  filter(term %in% c("AvgPAC_c", "TreatFreq")) %>%
  mutate(Sig = if_else(p.value < 0.1, "yes", "no"),
         Dir = if_else(estimate > 0, "pos", "neg"),
         term = fct_recode(term, "PAC" = "AvgPAC_c",
                           "Treat" = "TreatFreq"),
         CountGroup = paste(term, Sig, Dir, sep = "_"),
         Habitat = fct_relevel(Habitat, "floating")) %>%
  group_by(CommonName, Habitat, CountGroup) %>%
  summarize(Taxa = n_distinct(TaxonName)) %>%
  ungroup() %>%
  pivot_wider(names_from = CountGroup,
              values_from = Taxa) %>%
  mutate(across(.cols = where(is.integer), ~ replace_na(.x, 0)),
         PAC_no = PAC_no_neg + PAC_no_pos,
         Treat_no = Treat_no_neg + Treat_no_pos,
         CommonName = tolower(CommonName),
         CommonName = str_replace(CommonName, "cuban", "Cuban")) %>%
  select(CommonName, Habitat, PAC_yes_neg, PAC_yes_pos, PAC_no, Treat_yes_neg, Treat_yes_pos, Treat_no)

# split by invasive group
foc_plant_coef_summ <- plant_coef_summ %>%
  filter(CommonName %in% c("hydrilla", "water hyacinth", "water lettuce"))

non_foc_plant_coef_summ <- plant_coef_summ %>%
  filter(CommonName %in% c("Cuban bulrush", "para grass", "torpedograss"))

write_csv(foc_plant_coef_summ, "output/fwc_focal_invasive_native_detected_PAC_treatment_summary.csv")
write_csv(non_foc_plant_coef_summ, "output/fwc_non_focal_invasive_native_detected_PAC_treatment_summary.csv")


#### predicted values ####

# add coefficients to full data
comb_dat_PAC <- comb_dat2 %>%
  left_join(plant_coef_PAC) %>%
  mutate(Pred = logit2prob(Intercept + Coef * AvgPAC_c) * YearsSurveyed,
         Sig = fct_relevel(Sig, "yes"),
         PanelName = case_when(CommonName == "Hydrilla" ~ "(A) hydrilla",
                               CommonName == "Water hyacinth" ~ "(B) water hyacinth",
                               CommonName == "Water lettuce" ~ "(C) water lettuce",
                               CommonName == "Cuban bulrush" ~ "(A) Cuban bulrush",
                               CommonName == "Para grass" ~ "(B) para grass",
                               CommonName == "Torpedograss" ~ "(C) torpedograss"),
         Habitat = fct_relevel(Habitat, "floating"),
         Xval = AvgPAC) 

comb_dat_treat <- comb_dat2 %>%
  left_join(plant_coef_treat) %>%
  mutate(Pred = logit2prob(Intercept + Coef * TreatFreq) * YearsSurveyed,
         Sig = fct_relevel(Sig, "yes"),
         Treat = TreatFreq * YearsSurveyed,
         PanelName = case_when(CommonName == "Hydrilla" ~ "(A) hydrilla management",
                               CommonName == "Water hyacinth" ~ "(B) water hyacinth management",
                               CommonName == "Water lettuce" ~ "(C) water lettuce management",
                               CommonName == "Cuban bulrush" ~ "(A) Cuban bulrush management",
                               CommonName == "Para grass" ~ "(B) para grass management",
                               CommonName == "Torpedograss" ~ "(C) torpedograss management"),
         Habitat = fct_relevel(Habitat, "floating"),
         Xval = Treat)

# split by focal/non-focal
foc_dat_PAC <- comb_dat_PAC %>%
  filter(CommonName %in% c("Hydrilla", "Water hyacinth", "Water lettuce"))

foc_dat_treat <- comb_dat_treat %>%
  filter(CommonName %in% c("Hydrilla", "Water hyacinth", "Water lettuce"))

non_foc_dat_PAC <- comb_dat_PAC %>%
  filter(CommonName %in% c("Cuban bulrush", "Para grass", "Torpedograss"))

non_foc_dat_treat <- comb_dat_treat %>%
  filter(CommonName %in% c("Cuban bulrush", "Para grass", "Torpedograss"))



#### figures ####

# PAC figures
foc_fig_PAC <- ggplot(foc_dat_PAC, aes(x = AvgPAC, y = Pred, group = TaxonName, color = Habitat)) +
  geom_line(aes(alpha = Sig, size = Sig)) +
  facet_wrap(~ PanelName, scales = "free") +
  scale_alpha_manual(values = c(1, 0.5), name = "P < 0.1") +
  scale_size_manual(values = c(0.4, 0.25), name = "P < 0.1") +
  labs(x = "Invasive plant PAC", y = "Native abundance\n(years detected)") +
  scale_color_manual(values = brewer.set2(n = 3), name = "Native\nhabitat") +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0)) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 0.4, alpha = 1)))

non_foc_fig_PAC <- ggplot(non_foc_dat_PAC, aes(x = AvgPAC, y = Pred, group = TaxonName, color = Habitat)) +
  geom_line(aes(alpha = Sig, size = Sig)) +
  facet_wrap(~ PanelName, scales = "free") +
  scale_alpha_manual(values = c(1, 0.5), name = "P < 0.1") +
  scale_size_manual(values = c(0.4, 0.25), name = "P < 0.1") +
  labs(x = "Invasive plant PAC", y = "Native abundance\n(years detected)") +
  scale_color_manual(values = brewer.set2(n = 3), name = "Native\nhabitat") +
  def_theme_paper +
  theme(strip.text = element_text(size = 9, color = "black", hjust = 0)) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 0.4, alpha = 1)))

# can't adjust strip text titles to fit management in title
# make separate panels
# figure template
panel_template <- function(dat) {
  
  dat2 <- dat %>%
    mutate(label_y = max(Pred),
           label_x = (max(Xval) + min(Xval)) / 2)
  
  ggplot(dat2, aes(x = Xval, y = Pred, group = TaxonName, color = Habitat)) +
    geom_line(aes(alpha = Sig, size = Sig)) +
    facet_wrap(~ Habitat) +
    geom_text(aes(x = label_x, y = label_y, label = Habitat), 
              check_overlap = T, size = paper_text_size, 
              hjust = 0.5, vjust = 1, color = "black") +
    scale_alpha_manual(values = c(1, 0.5)) +
    scale_size_manual(values = c(0.4, 0.25)) +
    scale_color_manual(values = brewer.set2(n = 3)) +
    def_theme_paper +
    theme(legend.position = "none",
          strip.text = element_blank(),
          plot.title = element_text(size = 9, hjust = 0))
  
}

# focal PAC
hydr_fig_PAC <- panel_template(filter(foc_dat_PAC, CommonName == "Hydrilla")) +
  labs(y = "", title = "(A) hydrilla") +
  theme(axis.title.x = element_blank())
wahy_fig_PAC <- panel_template(filter(foc_dat_PAC, CommonName == "Water hyacinth")) +
  labs(y = "", title = "(B) water hyacinth") +
  theme(axis.title.x = element_blank())
wale_fig_PAC <- panel_template(filter(foc_dat_PAC, CommonName == "Water lettuce")) +
  labs(y = "", x = "Invasive plant PAC", title = "(C) water lettuce")

foc_fig_PAC <- hydr_fig_PAC + wahy_fig_PAC + wale_fig_PAC + 
  plot_layout(ncol = 1) + ylab("Native abundance (years detected)") +
  theme(axis.title.y = element_text(size = 9, hjust = -1.2))

# non-focal PAC
cubu_fig_PAC <- panel_template(filter(non_foc_dat_PAC, CommonName == "Cuban bulrush")) +
  labs(y = "", title = "(A) Cuban bulrush") +
  theme(axis.title.x = element_blank())
pagr_fig_PAC <- panel_template(filter(non_foc_dat_PAC, CommonName == "Para grass")) +
  labs(y = "", title = "(B) para grass") +
  theme(axis.title.x = element_blank())
torp_fig_PAC <- panel_template(filter(non_foc_dat_PAC, CommonName == "Torpedograss")) +
  labs(x = "Invasive plant PAC", title = "(C) torpedograss")

non_foc_fig_PAC <- cubu_fig_PAC + pagr_fig_PAC + torp_fig_PAC + 
  plot_layout(ncol = 1) + ylab("Native abundance (years detected)") +
  theme(axis.title.y = element_text(size = 9, hjust = -1.2))

# focal treatment
hydr_fig_treat <- panel_template(filter(foc_dat_treat, CommonName == "Hydrilla")) +
  labs(y = "", title = "(A) hydrilla management") +
  theme(axis.title.x = element_blank())
wahy_fig_treat <- panel_template(filter(foc_dat_treat, CommonName == "Water hyacinth")) +
  labs(y = "", title = "(B) water hyacinth management") +
  theme(axis.title.x = element_blank())
wale_fig_treat <- panel_template(filter(foc_dat_treat, CommonName == "Water lettuce")) +
  labs(y = "", x = "Years managed", title = "(C) water lettuce management")

foc_fig_treat <- hydr_fig_treat + wahy_fig_treat + wale_fig_treat + 
  plot_layout(ncol = 1) + ylab("Native abundance (years detected)") +
  theme(axis.title.y = element_text(size = 9, hjust = -1.2))

# non-focal treatment
cubu_fig_treat <- panel_template(filter(non_foc_dat_treat, CommonName == "Cuban bulrush")) +
  labs(y = "", title = "(A) Cuban bulrush management") +
  theme(axis.title.x = element_blank())
pagr_fig_treat <- panel_template(filter(non_foc_dat_treat, CommonName == "Para grass")) +
  labs(y = "", title = "(B) para grass management") +
  theme(axis.title.x = element_blank())
torp_fig_treat <- panel_template(filter(non_foc_dat_treat, CommonName == "Torpedograss")) +
  labs(x = "Years managed", title = "(C) torpedograss management")

non_foc_fig_treat <- cubu_fig_treat + pagr_fig_treat + torp_fig_treat + 
  plot_layout(ncol = 1) + ylab("Native abundance (years detected)") +
  theme(axis.title.y = element_text(size = 9, hjust = -1.2))

# save figures
ggsave("output/fwc_focal_invasive_native_detected_PAC_prediction.png", foc_fig_PAC,
       device = "png", width = 4, height = 5, units = "in")
ggsave("output/fwc_focal_invasive_native_detected_treatment_prediction.png", foc_fig_treat,
       device = "png", width = 4, height = 5, units = "in")
ggsave("output/fwc_non_focal_invasive_native_detected_PAC_prediction.png", non_foc_fig_PAC,
       device = "png", width = 4, height = 5, units = "in")
ggsave("output/fwc_non_focal_invasive_native_detected_treatment_prediction.png", non_foc_fig_treat,
       device = "png", width = 4, height = 5, units = "in")


#### older code ####

plant_coef <- plant_coef1 %>%
  filter(!(TaxonName %in% c("Mayaca fluviatilis", "Myriophyllum laxum/pinnatum"))) %>%
  full_join(plant_coef2 %>%
              filter(TaxonName %in% c("Mayaca fluviatilis", "Myriophyllum laxum/pinnatum"))) %>%
  mutate(sig = case_when(p.value < 0.1 & p.value >= 0.05 ~ "< 0.1",
                         p.value < 0.05 & p.value >= 0.01 ~ "< 0.05",
                         p.value < 0.01 & p.value >= 0.001 ~ "< 0.01",
                         p.value < 0.001 ~ "< 0.001",
                         TRUE ~ "\u2265 0.1"),
         term = fct_recode(term, "intercept" = "(Intercept)",
                           "floating\ntreat." = "FloatingTrt",
                           "hydrilla" = "Hydrilla",
                           "hydrilla\ntreat." = "HydrillaTrt",
                           "water\nhyacinth" = "WaterHyacinth",
                           "water\nlettuce" = "WaterLettuce") %>%
           fct_relevel("intercept", "hydrilla", "water\nhyacinth",
                       "water\nlettuce", "hydrilla\ntreat.", "floating\ntreat.")) %>%
  left_join(nat_plant3 %>%
              select(TaxonName, Habitat) %>%
              unique()) %>%
  mutate(HabitatLabel = fct_recode(Habitat,
                                   "(B) emersed taxa" = "Emersed",
                                   "(C) floating taxa" = "Floating",
                                   "(D) submersed taxa" = "Submersed"))

# divide coefficients by group
plant_coef_fig <- plant_coef %>%
  mutate(HabitatLabel = "(A) all taxa") %>%
  full_join(plant_coef)


#### one-sample t-test ####

# use all coefficients
plant_coef_test <- plant_coef %>%
  mutate(HabitatLabel = "(A) all taxa")

# t-tests
plant_coef_test2 <- plant_coef_test %>%
  select(term, estimate) %>%
  nest(data = estimate) %>%
  mutate(fit = map(data, ~ t.test(.)))

# summary
plant_coef_test3 <- plant_coef_test2 %>%
  mutate(tidied = map(fit, tidy)) %>%
  select(term, tidied) %>%
  unnest(tidied) %>%
  mutate(sig = case_when(p.value < 0.1 & p.value >= 0.05 ~ ".",
                           p.value < 0.05 & p.value >= 0.01 ~ "*",
                           p.value < 0.01 & p.value >= 0.001 ~ "**",
                           p.value < 0.001 ~ "***",
                           TRUE ~ NA_character_),
         HabitatLabel = "(A) all taxa")

# save
write_csv(plant_coef_test3, "output/fwc_native_plant_ttests.csv")


#### figure ####
cairo_pdf("output/fwc_native_plant_coefficients.pdf", width = 6.5, height = 6.5)
ggplot(plant_coef_fig , aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = NA, outlier.shape = 4, size = 0.25) +
  geom_point(aes(fill = sig), shape = 21, 
             position = position_jitter(width = 0.1), size = 1, stroke = 0.3) +
  geom_text(data = plant_coef_test3, aes(y = 5.1, label = sig)) +
  scale_fill_viridis_d(name = expression(paste(italic(P), " value", sep = ""))) +
  labs(y = "Estimate (log-odds)") +
  facet_wrap(~ HabitatLabel, scales = "free_y") +
  def_theme_paper +
  theme(strip.text = element_text(size = 11, color="black", hjust = 0),
        axis.title.x = element_blank(),
        legend.position = c(0.05, 0.4),
        legend.box = "horizontal",
        legend.spacing.x = unit(-0.1, "cm"),
        legend.box.background = element_rect(fill = NA, color = NA),
        legend.key.height = unit(2, "mm"))
dev.off()


#### table ####

plant_table <- plant_coef %>%
  relocate(TaxonName, Habitat) %>%
  mutate(term = str_replace(term, "\n", " ") %>%
           fct_relevel("intercept", "hydrilla", "water hyacinth",
                       "water lettuce", "hydrilla treat.", "floating treat."),
         Habitat = tolower(Habitat),
         across(where(is.double) & !p.value, ~ round_half_up(., digits = 3))) %>%
  rename("Taxon" = "TaxonName",
         "Parameter" = "term",
         "Estimate" = "estimate",
         "Std.error" = "std.error",
         "Z" = "statistic",
         "P" = "p.value",
         "Significance" = "sig") %>%
  select(-c(HabitatLabel)) %>%
  arrange(Parameter, Taxon)

write_csv(plant_table, "output/fwc_native_plant_coefficients.csv")
