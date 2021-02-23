# load packages
library(rstan)

# do in parallel
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())