#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(rstan)

# avoid recompiling unchanged Stan program
rstan_options(auto_write = TRUE)

# do in parallel
options(mc.cores = parallel::detectCores())

# import data
dat <- read_csv("intermediate-data/simulated_data_pre_herbicide_model.csv")


#### prepare data ####

# select first year value for each lake
y0 <- apply(dat, 1, function(x) first(na.omit(x)))

# make matrix
y_mat <- as.matrix(dat)

# make vector
y_vec <- y_mat[!is.na(y_mat)]

# column and row indices
y_indx <- which(!is.na(y_mat), arr.ind = TRUE)  # index of the non-NAs
y_col <- as.vector(y_indx[, "col"])
y_row <- as.vector(y_indx[, "row"])

# combine input data
input_dat <- list(TT = ncol(dat),
                  N = nrow(dat),
                  n_pos = length(y_vec),
                  col_indx_pos = y_col,
                  row_indx_pos = y_row,
                  y = y_vec,
                  y0 = y0)


#### fit model 1 ####

# model file
pre_herb_mod <- 'models/pre_herbicide_model_1.stan'

# time model fit
time1 <- Sys.time()

# fit model
pre_herb_fit <- rstan::stan(file = pre_herb_mod, 
                            data = input_dat, 
                            pars = c("sd_q", "x", "sd_r","u", "x0"), 
                            iter = 4000, chains = 3, thin = 1, warmup = 1000)

time2 <- Sys.time()

# time
print(time2 - time1)


#### diagnostics model 1 ####

# summary
print(pre_herb_fit)

# pairs plots
pairs(pre_herb_fit, pars = c("lp__", "sd_q"))
# slight funnel shape, but not where divergences are concentrated
pairs(pre_herb_fit, pars = c("lp__", "sd_r"))
pairs(pre_herb_fit, pars = c("lp__", "x0"))
# transitions are below diagonal: can't fix by increasing adapt_delta

# exctract draws
params_fit <- as.tibble(extract(pre_herb_fit, permuted=FALSE))

# rename columns
params_fit2 <- params_fit %>% 
  rename_with(~ gsub(":", "_", .x, fixed = T)) %>%
  mutate(iter = 1:nrow(params_fit))

# sd q (easiest to work with)
ggplot(params_fit2, aes(iter, log(chain_1.sd_q))) +
  geom_point(color = "blue") +
  geom_point(aes(y = log(chain_2.sd_q)), color = "green") +
  geom_point(aes(y = log(chain_3.sd_q)), color = "yellow")

# true value
set.seed(20)
true_q <- abs(rt(n = 1, df = 3)) 

# running mean
running_means <- sapply(params_fit2$iter, function(n) mean(log(params_fit2$chain_1.sd_q)[1:n]))
params_fit3 <- params_fit2 %>%
  mutate(sd_q_1_run = running_means)

ggplot(params_fit3, aes(iter, sd_q_1_run)) +
  geom_point() +
  geom_hline(yintercept = log(true_q))

# divergent transitions
divergent <- get_sampler_params(pre_herb_fit, inc_warmup=FALSE)[[1]][,'divergent__']

params_fit4 <- params_fit3 %>%
  mutate(divergent = divergent)

ggplot(params_fit4, aes(chain_1.lp__, log(chain_1.sd_q))) +
  geom_point(aes(color = as.factor(divergent)), alpha = 0.5)
# not at low values of sd_q, but at high values of lp__

ggplot(params_fit4, aes(chain_2.lp__, log(chain_2.sd_q))) +
  geom_point(aes(color = as.factor(divergent)), alpha = 0.5)
# distributed throughout

ggplot(params_fit4, aes(chain_3.lp__, log(chain_3.sd_q))) +
  geom_point(aes(color = as.factor(divergent)), alpha = 0.5)
# distributed throughout

stan_diag(pre_herb_fit, info = "sample") 
# don't know how to interpret this, second plot is accept_stat, which doesn't have a clear definition


#### fit model 2 ####

# model file
pre_herb_mod_2 <- 'models/pre_herbicide_model_2.stan'

# time model fit
time1 <- Sys.time()

# fit model
pre_herb_fit_2 <- rstan::stan(file = pre_herb_mod_2, 
                            data = input_dat, 
                            pars = c("sd_q", "pro_dev", "sd_r","u", "x0"), 
                            iter = 4000, chains = 3, thin = 1, warmup = 1000)

time2 <- Sys.time()

# time
print(time2 - time1)
