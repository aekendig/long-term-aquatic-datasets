# load packages
library(rstan)

# avoid recompiling unchanged Stan program
rstan_options(auto_write = TRUE)

# do in parallel
options(mc.cores = parallel::detectCores())