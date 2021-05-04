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
y0 <- dat[ , which()]